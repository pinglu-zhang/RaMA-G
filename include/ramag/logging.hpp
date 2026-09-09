#pragma once
#include <exception>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
namespace ramag {
// CLI-owned logger. Library callers opt in explicitly; no global registration.
class RunLogger {
 public:
  explicit RunLogger(const std::filesystem::path& work);
  ~RunLogger();
  RunLogger(const RunLogger&) = delete;
  RunLogger& operator=(const RunLogger&) = delete;
  void Info(std::string_view message, bool terminal = true);
  void Error(std::string_view message);
  void Flush();
  const std::string& Id() const;
  std::filesystem::path Path() const;
 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
int FailureExitCode(const std::exception& error) noexcept;
}  // namespace ramag
