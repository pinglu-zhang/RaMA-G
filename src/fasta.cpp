#include "ramag/fasta.hpp"

#include "ramag/runtime.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <system_error>
#include <unordered_set>

#include <zlib.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace ramag {
namespace {

[[nodiscard]] std::string PathForMessage(const std::filesystem::path& path) {
    return path.empty() ? std::string{"<empty path>"} : path.string();
}

[[nodiscard]] std::string Lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

[[nodiscard]] bool HasGzipSuffix(const std::filesystem::path& path) {
    const std::string name = Lowercase(path.filename().string());
    constexpr std::array<std::string_view, 3> suffixes{".gz", ".bgz", ".bgzf"};
    return std::any_of(suffixes.begin(), suffixes.end(), [&](std::string_view suffix) {
        return name.size() >= suffix.size() &&
               std::string_view{name}.substr(name.size() - suffix.size()) == suffix;
    });
}

[[nodiscard]] bool HasUnsupportedCompressedSuffix(
    const std::filesystem::path& path) {
    const std::string name = Lowercase(path.filename().string());
    constexpr std::array<std::string_view, 4> suffixes{
        ".bz2", ".xz", ".zst", ".zip"};
    return std::any_of(suffixes.begin(), suffixes.end(), [&](std::string_view suffix) {
        return name.size() >= suffix.size() &&
               std::string_view{name}.substr(name.size() - suffix.size()) == suffix;
    });
}

[[nodiscard]] bool HasGzipMagic(const std::array<unsigned char, 6>& bytes,
                                std::size_t size) {
    return size >= 2 && bytes[0] == 0x1fU && bytes[1] == 0x8bU;
}

[[nodiscard]] bool HasUnsupportedCompressedMagic(
    const std::array<unsigned char, 6>& bytes, std::size_t size) {
    const bool bzip2 = size >= 3 && bytes[0] == 'B' && bytes[1] == 'Z' && bytes[2] == 'h';
    const bool xz = size >= 6 && bytes[0] == 0xfdU && bytes[1] == 0x37U &&
                    bytes[2] == 0x7aU && bytes[3] == 0x58U && bytes[4] == 0x5aU &&
                    bytes[5] == 0x00U;
    const bool zstd = size >= 4 && bytes[0] == 0x28U && bytes[1] == 0xb5U &&
                      bytes[2] == 0x2fU && bytes[3] == 0xfdU;
    const bool zip = size >= 4 && bytes[0] == 'P' && bytes[1] == 'K' &&
                     (bytes[2] == 0x03U || bytes[2] == 0x05U || bytes[2] == 0x07U) &&
                     (bytes[3] == 0x04U || bytes[3] == 0x06U || bytes[3] == 0x08U);
    return bzip2 || xz || zstd || zip;
}

[[noreturn]] void ThrowLineError(const std::filesystem::path& path,
                                 std::uint64_t line_number,
                                 const std::string& message) {
    std::ostringstream stream;
    stream << "invalid FASTA '" << PathForMessage(path) << "' at line " << line_number
           << ": " << message;
    throw FastaError(stream.str());
}

[[nodiscard]] char NormalizeBase(unsigned char base) {
    const char upper = static_cast<char>(std::toupper(base));
    switch (upper) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return upper;
        default:
            // The public contract intentionally treats every other alphabetic
            // symbol (IUPAC ambiguity, U, X, etc.) as one seed hard break.
            return 'N';
    }
}

[[nodiscard]] char ComplementBase(char base) noexcept {
    switch (base) {
        case 'A':
        case 'a':
            return 'T';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';
        default:
            return 'N';
    }
}

[[nodiscard]] std::string TrimmedHeaderText(std::string_view line) {
    std::size_t begin = 1;
    while (begin < line.size() &&
           std::isspace(static_cast<unsigned char>(line[begin])) != 0) {
        ++begin;
    }
    std::size_t end = line.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(line[end - 1])) != 0) {
        --end;
    }
    return std::string{line.substr(begin, end - begin)};
}

[[nodiscard]] std::string HeaderId(std::string_view header) {
    const std::size_t separator = header.find_first_of(" \t\v\f\r\n");
    return std::string{header.substr(0, separator)};
}

void CheckedIncrement(Length& value, Length increment, const std::string& what) {
    if (increment > std::numeric_limits<Length>::max() - value) {
        throw FastaError(what + " exceeds the 64-bit RaMA-G coordinate range");
    }
    value += increment;
}

void ForEachPlainLine(const std::filesystem::path& path,
                      const std::function<void(std::string)>& callback) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw FastaError("cannot open FASTA '" + PathForMessage(path) + "'");
    std::string line;
    while (std::getline(input, line)) {
        CheckInterruption("input");
        callback(std::move(line));
    }
    if (input.bad()) {
        throw FastaError("I/O error while reading FASTA '" + PathForMessage(path) + "'");
    }
}

void ForEachGzipLine(const std::filesystem::path& path,
                     const std::function<void(std::string)>& callback) {
    gzFile input = gzopen(path.c_str(), "rb");
    if (input == nullptr) {
        throw FastaError("cannot open gzip FASTA '" + PathForMessage(path) + "'");
    }
    std::array<char, 64U * 1024U> buffer{};
    std::string line;
    try {
        while (true) {
            CheckInterruption("input");
            const int available = gzread(input, buffer.data(),
                                         static_cast<unsigned int>(buffer.size()));
            if (available < 0) {
                int code = Z_OK;
                const char* message = gzerror(input, &code);
                throw FastaError("invalid or truncated gzip FASTA '" +
                                 PathForMessage(path) + "': " +
                                 (message == nullptr ? "zlib read error" : message));
            }
            if (available == 0) break;
            for (int index = 0; index < available; ++index) {
                const char value = buffer[static_cast<std::size_t>(index)];
                if (value == '\n') {
                    callback(std::move(line));
                    line.clear();
                } else {
                    if (line.size() == std::numeric_limits<std::size_t>::max()) {
                        throw FastaError("gzip FASTA line is too large");
                    }
                    line.push_back(value);
                }
            }
        }
        if (!line.empty()) callback(std::move(line));
    } catch (...) {
        gzclose(input);
        throw;
    }
    const int close_result = gzclose(input);
    if (close_result != Z_OK) {
        throw FastaError("gzip integrity validation failed for '" +
                         PathForMessage(path) + "'");
    }
}

}  // namespace

bool IsCanonicalBase(char base) noexcept {
    switch (base) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'a':
        case 'c':
        case 'g':
        case 't':
            return true;
        default:
            return false;
    }
}

std::string ReverseComplement(std::string_view bases) {
    std::string result;
    result.resize(bases.size());
    for (std::size_t index = 0; index < bases.size(); ++index) {
        result[index] = ComplementBase(bases[bases.size() - index - 1]);
    }
    return result;
}

std::vector<std::string> ReverseComplements(
    const std::vector<SequenceRecord>& sequences,
    std::uint32_t worker_threads) {
    if (worker_threads == 0) {
        throw FastaError("reverse-complement worker_threads must be greater than zero");
    }

    std::vector<std::string> results(sequences.size());
    for (std::size_t sequence_index = 0; sequence_index < sequences.size();
         ++sequence_index) {
        results[sequence_index].resize(sequences[sequence_index].bases.size());
    }

    struct WorkChunk {
        std::size_t sequence_index{};
        std::size_t output_begin{};
        std::size_t output_end{};
    };
    constexpr std::size_t kChunkBases = 1U * 1024U * 1024U;
    std::vector<WorkChunk> chunks;
    std::uint64_t total_bases = 0;
    for (std::size_t sequence_index = 0; sequence_index < sequences.size();
         ++sequence_index) {
        const std::size_t length = sequences[sequence_index].bases.size();
        if (length > std::numeric_limits<std::uint64_t>::max() - total_bases) {
            throw FastaError(
                "reverse-complement collection exceeds the 64-bit size range");
        }
        total_bases += static_cast<std::uint64_t>(length);
        for (std::size_t begin = 0; begin < length;) {
            const std::size_t width = std::min(kChunkBases, length - begin);
            chunks.push_back(WorkChunk{sequence_index, begin, begin + width});
            begin += width;
        }
    }
    if (chunks.size() >
        static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())) {
        throw FastaError("reverse-complement task count exceeds the OpenMP loop range");
    }

    const auto process_chunk = [&](std::size_t chunk_index) noexcept {
        const WorkChunk chunk = chunks[chunk_index];
        const std::string& input = sequences[chunk.sequence_index].bases;
        std::string& output = results[chunk.sequence_index];
        for (std::size_t output_index = chunk.output_begin;
             output_index < chunk.output_end; ++output_index) {
            output[output_index] =
                ComplementBase(input[input.size() - output_index - 1]);
        }
    };

    constexpr std::uint64_t kParallelThresholdBases = 1U * 1024U * 1024U;
#if defined(_OPENMP)
    const bool may_enter_parallel =
        worker_threads > 1 && chunks.size() > 1 &&
        total_bases >= kParallelThresholdBases && omp_in_parallel() == 0;
    if (may_enter_parallel) {
        const std::uint64_t planned = std::min<std::uint64_t>(
            {worker_threads, static_cast<std::uint64_t>(chunks.size()),
             static_cast<std::uint64_t>(std::numeric_limits<int>::max())});
        const int omp_workers = static_cast<int>(planned);
#pragma omp parallel for schedule(static) num_threads(omp_workers) default(none) \
    shared(chunks, process_chunk)
        for (std::ptrdiff_t chunk_index = 0;
             chunk_index < static_cast<std::ptrdiff_t>(chunks.size());
             ++chunk_index) {
            process_chunk(static_cast<std::size_t>(chunk_index));
        }
        return results;
    }
#else
    static_cast<void>(total_bases);
#endif
    for (std::size_t chunk_index = 0; chunk_index < chunks.size(); ++chunk_index) {
        process_chunk(chunk_index);
    }
    return results;
}

FastaCompression DetectFastaCompression(const std::filesystem::path& path) {
    if (path.empty()) throw FastaError("FASTA path is empty");
    std::ifstream input(path, std::ios::binary);
    if (!input) throw FastaError("cannot open FASTA '" + PathForMessage(path) + "'");
    std::array<unsigned char, 6> magic{};
    input.read(reinterpret_cast<char*>(magic.data()),
               static_cast<std::streamsize>(magic.size()));
    const std::size_t size = static_cast<std::size_t>(input.gcount());
    if (HasUnsupportedCompressedSuffix(path) ||
        HasUnsupportedCompressedMagic(magic, size)) {
        throw UnsupportedFastaFormat(
            "only plain and gzip FASTA are supported: '" + PathForMessage(path) + "'");
    }
    if (HasGzipMagic(magic, size)) return FastaCompression::Gzip;
    if (HasGzipSuffix(path)) {
        throw FastaError("FASTA has a gzip suffix but no gzip magic: '" +
                         PathForMessage(path) + "'");
    }
    return FastaCompression::Plain;
}

FastaData ReadFasta(const std::filesystem::path& path, SequenceId first_id) {
    const auto compression = DetectFastaCompression(path);
    FastaData result;
    std::error_code absolute_error;
    result.source = std::filesystem::absolute(path, absolute_error);
    if (absolute_error) {
        result.source = path;
    }

    std::unordered_set<std::string> seen_names;
    std::uint64_t line_number = 0;
    SequenceRecord current;
    bool have_record = false;
    bool saw_any_content = false;

    auto finish_record = [&]() {
        if (!have_record) {
            return;
        }
        if (current.bases.empty()) {
            throw FastaError("FASTA record '" + current.name + "' in '" +
                             PathForMessage(path) + "' has an empty sequence");
        }
        result.sequences.push_back(std::move(current));
        current = SequenceRecord{};
        have_record = false;
    };

    const auto consume_line = [&](std::string line) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            if (!have_record) {
                ThrowLineError(path, line_number, "blank content before the first header");
            }
            // Blank lines inside a record do not contribute bases.  A record
            // made only of blank lines is still rejected by finish_record().
            return;
        }
        saw_any_content = true;

        if (line.front() == '>') {
            finish_record();
            const std::string header = TrimmedHeaderText(line);
            if (header.empty()) {
                ThrowLineError(path, line_number, "header has no identifier");
            }
            const std::string name = HeaderId(header);
            if (name.empty()) {
                ThrowLineError(path, line_number, "header has no identifier");
            }
            if (!seen_names.insert(name).second) {
                ThrowLineError(path, line_number, "duplicate record identifier '" + name + "'");
            }
            const std::uint64_t offset = static_cast<std::uint64_t>(result.sequences.size());
            if (offset > static_cast<std::uint64_t>(
                             std::numeric_limits<SequenceId>::max() - first_id)) {
                throw FastaError("FASTA contains too many records for 32-bit sequence ids");
            }
            current.numeric_id = static_cast<SequenceId>(first_id + offset);
            current.name = name;
            current.header = header;
            have_record = true;
            return;
        }

        if (!have_record) {
            ThrowLineError(path, line_number, "sequence data appears before a header");
        }

        if (line.size() > std::numeric_limits<std::size_t>::max() - current.bases.size()) {
            throw FastaError("FASTA record '" + current.name + "' is too large");
        }
        current.bases.reserve(current.bases.size() + line.size());
        for (std::size_t column = 0; column < line.size(); ++column) {
            const unsigned char base = static_cast<unsigned char>(line[column]);
            if (std::isalpha(base) == 0) {
                std::ostringstream message;
                message << "invalid sequence character at column " << (column + 1);
                ThrowLineError(path, line_number, message.str());
            }
            const char normalized = NormalizeBase(base);
            current.bases.push_back(normalized);
            CheckedIncrement(result.total_bases, 1, "FASTA base count");
            if (normalized == 'N') {
                CheckedIncrement(result.ambiguous_bases, 1, "FASTA ambiguous-base count");
            }
        }
    };
    if (compression == FastaCompression::Gzip) {
        ForEachGzipLine(path, consume_line);
    } else {
        ForEachPlainLine(path, consume_line);
    }
    if (!saw_any_content) {
        throw FastaError("FASTA is empty: '" + PathForMessage(path) + "'");
    }
    finish_record();
    if (result.sequences.empty()) {
        throw FastaError("FASTA contains no records: '" + PathForMessage(path) + "'");
    }
    return result;
}

void ValidateUncompressedFasta(const std::filesystem::path& path) {
    if (DetectFastaCompression(path) != FastaCompression::Plain) {
        throw UnsupportedFastaFormat(
            "SeqPro mmap input requires uncompressed FASTA: '" +
            PathForMessage(path) + "'");
    }
}

}  // namespace ramag
