#pragma once

#include "ramag/cli.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

namespace ramag {

class InterruptedError : public std::runtime_error {
 public:
  InterruptedError(int signal_number, std::string stage);
  [[nodiscard]] int SignalNumber() const noexcept { return signal_number_; }
  [[nodiscard]] int ExitCode() const noexcept;

 private:
  int signal_number_{};
};

void InstallSignalHandlers();
void CheckInterruption(std::string_view stage);

class ProgressSession {
 public:
  ProgressSession(const ProgressOptions& options, std::string run_id,
                  std::uint32_t threads);
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
