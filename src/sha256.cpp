#include "ramag/sha256.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ramag {
namespace {

class Sha256 {
 public:
  void Update(const unsigned char* data, std::size_t size) {
    if (size > std::numeric_limits<std::uint64_t>::max() - total_) {
      throw std::overflow_error("SHA-256 input exceeds uint64 size");
    }
    total_ += static_cast<std::uint64_t>(size);
    while (size != 0) {
      const auto copied = std::min(size, block_.size() - buffered_);
      std::copy_n(data, copied, block_.begin() + static_cast<std::ptrdiff_t>(buffered_));
      data += copied;
      size -= copied;
      buffered_ += copied;
      if (buffered_ == block_.size()) {
        Transform();
        buffered_ = 0;
      }
    }
  }

  void Update(std::string_view text) {
    Update(reinterpret_cast<const unsigned char*>(text.data()), text.size());
  }

  void UpdateLength(std::uint64_t value) {
    std::array<unsigned char, 8> bytes{};
    for (std::size_t index = 0; index < bytes.size(); ++index) {
      bytes[index] = static_cast<unsigned char>((value >> (index * 8U)) & 0xffU);
    }
    Update(bytes.data(), bytes.size());
  }

  std::string Finalize() {
    if (total_ > std::numeric_limits<std::uint64_t>::max() / 8U) {
      throw std::overflow_error("SHA-256 bit length overflow");
    }
    const auto bit_length = total_ * 8U;
    block_[buffered_++] = 0x80U;
    if (buffered_ > 56U) {
      std::fill(block_.begin() + static_cast<std::ptrdiff_t>(buffered_), block_.end(), 0U);
      Transform();
      buffered_ = 0;
    }
    std::fill(block_.begin() + static_cast<std::ptrdiff_t>(buffered_),
              block_.begin() + 56, 0U);
    for (std::size_t index = 0; index < 8U; ++index) {
      block_[56U + index] = static_cast<unsigned char>(
          (bit_length >> ((7U - index) * 8U)) & 0xffU);
    }
    Transform();
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (const auto word : state_) output << std::setw(8) << word;
    return output.str();
  }

 private:
  static constexpr std::array<std::uint32_t, 64> k{
      0x428a2f98U,0x71374491U,0xb5c0fbcfU,0xe9b5dba5U,0x3956c25bU,0x59f111f1U,0x923f82a4U,0xab1c5ed5U,
      0xd807aa98U,0x12835b01U,0x243185beU,0x550c7dc3U,0x72be5d74U,0x80deb1feU,0x9bdc06a7U,0xc19bf174U,
      0xe49b69c1U,0xefbe4786U,0x0fc19dc6U,0x240ca1ccU,0x2de92c6fU,0x4a7484aaU,0x5cb0a9dcU,0x76f988daU,
      0x983e5152U,0xa831c66dU,0xb00327c8U,0xbf597fc7U,0xc6e00bf3U,0xd5a79147U,0x06ca6351U,0x14292967U,
      0x27b70a85U,0x2e1b2138U,0x4d2c6dfcU,0x53380d13U,0x650a7354U,0x766a0abbU,0x81c2c92eU,0x92722c85U,
      0xa2bfe8a1U,0xa81a664bU,0xc24b8b70U,0xc76c51a3U,0xd192e819U,0xd6990624U,0xf40e3585U,0x106aa070U,
      0x19a4c116U,0x1e376c08U,0x2748774cU,0x34b0bcb5U,0x391c0cb3U,0x4ed8aa4aU,0x5b9cca4fU,0x682e6ff3U,
      0x748f82eeU,0x78a5636fU,0x84c87814U,0x8cc70208U,0x90befffaU,0xa4506cebU,0xbef9a3f7U,0xc67178f2U};

  void Transform() {
    std::array<std::uint32_t, 64> w{};
    for (std::size_t i = 0; i < 16; ++i) {
      const auto o = i * 4U;
      w[i] = (static_cast<std::uint32_t>(block_[o]) << 24U) |
             (static_cast<std::uint32_t>(block_[o + 1U]) << 16U) |
             (static_cast<std::uint32_t>(block_[o + 2U]) << 8U) |
             static_cast<std::uint32_t>(block_[o + 3U]);
    }
    for (std::size_t i = 16; i < w.size(); ++i) {
      const auto s0 = std::rotr(w[i - 15], 7) ^ std::rotr(w[i - 15], 18) ^ (w[i - 15] >> 3U);
      const auto s1 = std::rotr(w[i - 2], 17) ^ std::rotr(w[i - 2], 19) ^ (w[i - 2] >> 10U);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    auto a=state_[0], b=state_[1], c=state_[2], d=state_[3];
    auto e=state_[4], f=state_[5], g=state_[6], h=state_[7];
    for (std::size_t i = 0; i < w.size(); ++i) {
      const auto s1 = std::rotr(e,6) ^ std::rotr(e,11) ^ std::rotr(e,25);
      const auto ch = (e & f) ^ (~e & g);
      const auto t1 = h + s1 + ch + k[i] + w[i];
      const auto s0 = std::rotr(a,2) ^ std::rotr(a,13) ^ std::rotr(a,22);
      const auto maj = (a & b) ^ (a & c) ^ (b & c);
      const auto t2 = s0 + maj;
      h=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    state_[0]+=a; state_[1]+=b; state_[2]+=c; state_[3]+=d;
    state_[4]+=e; state_[5]+=f; state_[6]+=g; state_[7]+=h;
  }

  std::array<std::uint32_t, 8> state_{
      0x6a09e667U,0xbb67ae85U,0x3c6ef372U,0xa54ff53aU,
      0x510e527fU,0x9b05688cU,0x1f83d9abU,0x5be0cd19U};
  std::array<unsigned char, 64> block_{};
  std::size_t buffered_{};
  std::uint64_t total_{};
};

}  // namespace

std::string NormalizedReferenceSha256(
    std::span<const SequenceRecord> records) {
  Sha256 hash;
  hash.Update("RaMA-G-normalized-reference-v1");
  hash.UpdateLength(static_cast<std::uint64_t>(records.size()));
  for (const auto& record : records) {
    hash.UpdateLength(static_cast<std::uint64_t>(record.name.size()));
    hash.Update(record.name);
    hash.UpdateLength(static_cast<std::uint64_t>(record.header.size()));
    hash.Update(record.header);
    hash.UpdateLength(static_cast<std::uint64_t>(record.bases.size()));
    hash.Update(record.bases);
  }
  return hash.Finalize();
}

std::string Sha256Hex(std::string_view content) {
  Sha256 hash;
  hash.Update(content);
  return hash.Finalize();
}

std::string FileSha256(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open file for SHA-256: " + path.string());
  }
  Sha256 hash;
  std::array<unsigned char, 1024U * 1024U> buffer{};
  while (input) {
    input.read(reinterpret_cast<char*>(buffer.data()),
               static_cast<std::streamsize>(buffer.size()));
    const auto count = input.gcount();
    if (count > 0) hash.Update(buffer.data(), static_cast<std::size_t>(count));
  }
  if (input.bad()) {
    throw std::runtime_error("I/O failure while hashing file: " + path.string());
  }
  return hash.Finalize();
}

}  // namespace ramag
