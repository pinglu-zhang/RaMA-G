#pragma once
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <exception>
#include <filesystem>
#include <memory>



namespace ramag {
struct ProgressOptions;

class InterruptedError : public std::runtime_error {
 public:
  InterruptedError(int signal_number, std::string stage);
  [[nodiscard]] int SignalNumber() const noexcept { return signal_number_; }
  [[nodiscard]] int ExitCode() const noexcept;

 private:
  int signal_number_{};
};

// Current RSS, or zero when the platform cannot provide it.
[[nodiscard]] std::uint64_t CurrentResidentBytes();
void InstallSignalHandlers();
void CheckInterruption(std::string_view stage);

class RunLogger;

class ProgressSession {
 public:
  ProgressSession(const ProgressOptions& options, std::string run_id,
                  std::uint32_t threads, RunLogger* logger = nullptr);
  ~ProgressSession();
  ProgressSession(const ProgressSession&) = delete;
  ProgressSession& operator=(const ProgressSession&) = delete;

  void Stage(std::string stage, std::uint64_t completed = 0,
             std::uint64_t total = 0, std::string detail = {});
  void Update(std::uint64_t completed, std::uint64_t total = 0,
              std::string detail = {});
  void Finish(std::string detail = {});

 private:
  struct Impl;
  Impl* implementation_{};
};

}  // namespace ramag

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
