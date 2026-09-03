#include "ramag/cli.hpp"
#include "ramag/fasta.hpp"
#include "ramag/pipeline.hpp"
#include "ramag/reference_index.hpp"
#include "ramag/runtime.hpp"
#include "ramag/writers.hpp"

#include <filesystem>
#include <iostream>
#include <new>
#include <sstream>
#include <string>

namespace {

std::string QuoteArgument(std::string_view argument) {
  if (argument.find_first_of(" \t\n\"\\") == std::string_view::npos) {
    return std::string(argument);
  }
  std::string quoted{"\""};
  for (const char value : argument) {
    if (value == '\"' || value == '\\') {
      quoted.push_back('\\');
    }
    quoted.push_back(value);
  }
  quoted.push_back('\"');
  return quoted;
}

std::string Invocation(int argc, const char* const* argv) {
  std::ostringstream result;
  for (int index = 0; index < argc; ++index) {
    if (index != 0) {
      result << ' ';
    }
    result << QuoteArgument(argv[index]);
  }
  return result.str();
}

std::filesystem::path RunningBinaryPath(const char* argv_zero) {
#if defined(__linux__)
  std::error_code link_error;
  const auto proc_path = std::filesystem::read_symlink("/proc/self/exe", link_error);
  if (!link_error && !proc_path.empty()) {
    return proc_path;
  }
#endif
  std::error_code absolute_error;
  auto fallback = std::filesystem::absolute(argv_zero, absolute_error);
  return absolute_error ? std::filesystem::path{argv_zero} : fallback;
}

}  // namespace

int main(int argc, char** argv) {
  ramag::InstallSignalHandlers();
  try {
    const auto parsed = ramag::ParseCommandLine(argc, argv);
    if (parsed.show_help) {
      std::cout << ramag::HelpText();
      return 0;
    }
    if (parsed.show_version) {
      std::cout << ramag::VersionText() << '\n';
      return 0;
    }
    auto binary_path = RunningBinaryPath(argv[0]);
    if (parsed.command == ramag::CommandKind::Index) {
      ramag::ConfigureOpenMpRuntime(parsed.index_spec);
      if (parsed.print_effective_config) {
        std::cout << ramag::EffectiveIndexConfigText(parsed.index_spec);
        return 0;
      }
      const auto paths = ramag::RunReferenceIndexPipeline(
          parsed.index_spec, Invocation(argc, argv), binary_path);
      std::cout << "RaMA-G reference index completed: index=" << paths.index
                << "; marker=" << paths.complete << '\n';
      return 0;
    }
    ramag::ConfigureOpenMpRuntime(parsed.run_spec);
    if (parsed.print_effective_config) {
      std::cout << ramag::EffectiveConfigText(parsed.run_spec);
      return 0;
    }
    const auto outcome = ramag::RunAlignmentPipeline(
        parsed.run_spec, Invocation(argc, argv), binary_path);
    std::cout << "RaMA-G completed: " << outcome.statistics.alignment_count
              << " alignment(s); marker=" << outcome.paths.complete << '\n';
    return 0;
  } catch (const ramag::InterruptedError& error) {
    std::cerr << "ramag: interrupted: " << error.what() << '\n';
    return error.ExitCode();
  } catch (const ramag::UnsupportedSeedMode& error) {
    std::cerr << "ramag: unsupported: " << error.what() << '\n';
    return 4;
  } catch (const ramag::UnsupportedFastaFormat& error) {
    std::cerr << "ramag: unsupported: " << error.what() << '\n';
    return 4;
  } catch (const ramag::CliError& error) {
    std::cerr << "ramag: cli/config: " << error.what()
              << "\nRun 'ramag --help' for usage.\n";
    return 2;
  } catch (const ramag::FastaError& error) {
    std::cerr << "ramag: input: " << error.what() << '\n';
    return 3;
  } catch (const ramag::DependencyError& error) {
    std::cerr << "ramag: dependency/index: " << error.what() << '\n';
    return 5;
  } catch (const ramag::WriterError& error) {
    std::cerr << "ramag: output/validation: " << error.what() << '\n';
    return 7;
  } catch (const ramag::AlignmentError& error) {
    std::cerr << "ramag: compute/resource: " << error.what() << '\n';
    return 6;
  } catch (const std::bad_alloc& error) {
    std::cerr << "ramag: compute/resource: memory allocation failed: "
              << error.what() << '\n';
    return 6;
  } catch (const std::exception& error) {
    std::cerr << "ramag: internal: " << error.what() << '\n';
    return 8;
  }
}
