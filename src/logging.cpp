#include "ramag/logging.hpp"
#include "ramag/runtime.hpp"
#include "ramag/fasta.hpp"
#include "ramag/sufkit_adapter.hpp"
#include "ramag/writers.hpp"
#include <chrono>
#include <random>
#include <sstream>
#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
namespace ramag {
struct RunLogger::Impl {
  std::string id;
  std::filesystem::path path;
  std::shared_ptr<spdlog::logger> file, console;
};
RunLogger::RunLogger(const std::filesystem::path& work) : impl_(std::make_unique<Impl>()) {
  std::ostringstream id;
  id << std::hex << std::chrono::steady_clock::now().time_since_epoch().count()
     << '-' << std::random_device{}();
  impl_->id = id.str();
  const auto directory = work / "runs" / impl_->id;
  try {
    std::filesystem::create_directories(work / "runs");
    if (!std::filesystem::create_directory(directory)) throw WriterError("run directory collision");
    impl_->path = directory / "run.log";
    auto file = std::make_shared<spdlog::sinks::basic_file_sink_mt>(impl_->path.string(), false);
    auto console = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    file->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
    console->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [%l] %v%$");
    impl_->file = std::make_shared<spdlog::logger>(impl_->id+"-file", file);
    impl_->console = std::make_shared<spdlog::logger>(impl_->id+"-console", console);
    impl_->file->set_error_handler([](const std::string& message) { throw WriterError("log write failed: " + message); });
    impl_->file->flush_on(spdlog::level::info);
    Info("run_id=" + impl_->id + " log=" + impl_->path.string());
  } catch (const std::exception& error) {
    throw WriterError("cannot initialize run log: " + std::string(error.what()));
  }
}
RunLogger::~RunLogger() { try { Flush(); } catch (...) {} }
void RunLogger::Info(std::string_view message, bool terminal) {
  impl_->file->info("{}", message);
  if (terminal) impl_->console->info("{}", message);
}
void RunLogger::Error(std::string_view message) {
  impl_->file->error("{}", message); impl_->console->error("{}", message); Flush();
}
void RunLogger::Flush() { if(impl_->file) impl_->file->flush(); if(impl_->console) impl_->console->flush(); }
const std::string& RunLogger::Id() const { return impl_->id; }
std::filesystem::path RunLogger::Path() const { return impl_->path; }
int FailureExitCode(const std::exception& error) noexcept {
  if (const auto* e = dynamic_cast<const InterruptedError*>(&error)) return e->ExitCode();
  if (dynamic_cast<const UnsupportedSeedMode*>(&error) || dynamic_cast<const UnsupportedFastaFormat*>(&error)) return 4;
  if (dynamic_cast<const CliError*>(&error)) return 2;
  if (dynamic_cast<const FastaError*>(&error)) return 3;
  if (dynamic_cast<const DependencyError*>(&error)) return 5;
  if (dynamic_cast<const WriterError*>(&error)) return 7;
  if (dynamic_cast<const AlignmentError*>(&error) || dynamic_cast<const std::bad_alloc*>(&error)) return 6;
  return 8;
}
}  // namespace ramag
