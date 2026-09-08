#include "ramag/runtime.hpp"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <csignal>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <thread>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace ramag {
namespace {

volatile std::sig_atomic_t requested_signal = 0;
volatile std::sig_atomic_t snapshot_requested = 0;

extern "C" void HandleSignal(int signal_number) {
#if defined(SIGUSR1)
  if (signal_number == SIGUSR1) {
    snapshot_requested = 1;
    return;
  }
#endif
  if (requested_signal == 0) requested_signal = signal_number;
}

std::uint64_t CurrentRssBytes() {
#if defined(__linux__)
  std::ifstream input("/proc/self/status");
  std::string key;
  while (input >> key) {
    if (key == "VmRSS:") {
      std::uint64_t kib = 0;
      std::string unit;
      input >> kib >> unit;
      return kib * 1024U;
    }
    std::string rest;
    std::getline(input, rest);
  }
#endif
  return 0;
}

bool ProgressEnabled(ProgressMode mode) {
  if (mode == ProgressMode::On) return true;
  if (mode == ProgressMode::Off) return false;
#if defined(__unix__) || defined(__APPLE__)
  return isatty(fileno(stderr)) != 0;
#else
  return false;
#endif
}

}  // namespace

InterruptedError::InterruptedError(int signal_number, std::string stage)
    : std::runtime_error("interrupted by signal " + std::to_string(signal_number) +
                         " during " + stage),
      signal_number_(signal_number) {}

int InterruptedError::ExitCode() const noexcept {
  return signal_number_ == SIGINT ? 130 : signal_number_ == SIGTERM ? 143 : 128 + signal_number_;
}

void InstallSignalHandlers() {
  requested_signal = 0;
  snapshot_requested = 0;
  std::signal(SIGINT, HandleSignal);
  std::signal(SIGTERM, HandleSignal);
#if defined(SIGUSR1)
  std::signal(SIGUSR1, HandleSignal);
#endif
}

std::uint64_t CurrentResidentBytes() { return CurrentRssBytes(); }

void CheckInterruption(std::string_view stage) {
  const int value = static_cast<int>(requested_signal);
  if (value != 0) throw InterruptedError(value, std::string(stage));
}

struct ProgressSession::Impl {
  using Clock = std::chrono::steady_clock;
  ProgressOptions options;
  std::string run_id;
  std::uint32_t threads{};
  bool enabled{};
  std::mutex mutex;
  std::condition_variable condition;
  std::thread worker;
  bool stopping{};
  std::string stage{"starting"};
  std::string detail;
  std::uint64_t completed{};
  std::uint64_t total{};
  Clock::time_point started{Clock::now()};
  Clock::time_point stage_started{started};

  Impl(ProgressOptions value, std::string id, std::uint32_t worker_count)
      : options(value), run_id(std::move(id)), threads(worker_count),
        enabled(ProgressEnabled(options.mode)) {
    if (enabled) Print("start");
    // The watcher also services SIGUSR1 and reports pending interruption while
    // a third-party routine owns the main thread.  It therefore remains active
    // when periodic progress is disabled.
    worker = std::thread([this] { Run(); });
  }

  ~Impl() {
    {
      std::lock_guard lock(mutex);
      stopping = true;
    }
    condition.notify_all();
    if (worker.joinable()) worker.join();
  }

  void Run() {
    std::unique_lock lock(mutex);
    int reported_signal = 0;
    while (!stopping) {
      const auto interval = std::chrono::seconds(options.interval_seconds);
      condition.wait_for(lock, std::chrono::milliseconds(200));
      if (stopping) break;
      const bool snapshot = snapshot_requested != 0;
      if (snapshot) snapshot_requested = 0;
      const auto now = Clock::now();
      static thread_local auto last = started;
      const bool periodic = enabled && now - last >= interval;
      const int current_signal = static_cast<int>(requested_signal);
      const bool new_signal = current_signal != 0 && current_signal != reported_signal;
      if (snapshot || periodic || new_signal) {
        last = now;
        if (new_signal) reported_signal = current_signal;
        lock.unlock();
        Print(snapshot ? "snapshot" : new_signal ? "interrupt-pending" : "heartbeat");
        lock.lock();
      }
    }
  }

  void Print(std::string_view event) {
    std::string stage_copy;
    std::string detail_copy;
    std::uint64_t completed_copy = 0;
    std::uint64_t total_copy = 0;
    Clock::time_point stage_time;
    {
      std::lock_guard lock(mutex);
      stage_copy = stage;
      detail_copy = detail;
      completed_copy = completed;
      total_copy = total;
      stage_time = stage_started;
    }
    const auto now = Clock::now();
    std::ostringstream output;
    output << "[ramag progress] run_id=" << run_id << " event=" << event
           << " stage=" << stage_copy << " total_seconds=" << std::fixed
           << std::setprecision(1) << std::chrono::duration<double>(now - started).count()
           << " stage_seconds=" << std::chrono::duration<double>(now - stage_time).count();
    if (total_copy != 0) output << " completed=" << completed_copy << '/' << total_copy;
    else if (completed_copy != 0) output << " completed=" << completed_copy;
    output << " threads=" << threads << " rss_bytes=" << CurrentRssBytes();
    if (!detail_copy.empty()) output << " detail=" << detail_copy;
    std::cerr << output.str() << '\n';
  }
};

ProgressSession::ProgressSession(const ProgressOptions& options,
                                 std::string run_id, std::uint32_t threads)
    : implementation_(new Impl(options, std::move(run_id), threads)) {}

ProgressSession::~ProgressSession() { delete implementation_; }

void ProgressSession::Stage(std::string stage, std::uint64_t completed,
                            std::uint64_t total, std::string detail) {
  if (implementation_ == nullptr) return;
  {
    std::lock_guard lock(implementation_->mutex);
    implementation_->stage = std::move(stage);
    implementation_->completed = completed;
    implementation_->total = total;
    implementation_->detail = std::move(detail);
    implementation_->stage_started = Impl::Clock::now();
  }
  if (implementation_->enabled) implementation_->Print("stage");
  CheckInterruption(implementation_->stage);
  const char* paused_stage = std::getenv("RAMAG_TEST_PAUSE_STAGE");
  if (paused_stage != nullptr && implementation_->stage == paused_stage) {
    for (int attempt = 0; attempt < 1000; ++attempt) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
      CheckInterruption(implementation_->stage);
    }
  }
}

void ProgressSession::Update(std::uint64_t completed, std::uint64_t total,
                             std::string detail) {
  if (implementation_ == nullptr) return;
  {
    std::lock_guard lock(implementation_->mutex);
    implementation_->completed = completed;
    implementation_->total = total;
    implementation_->detail = std::move(detail);
  }
  CheckInterruption(implementation_->stage);
}

void ProgressSession::Finish(std::string detail) {
  Stage("complete", 1, 1, std::move(detail));
  if (implementation_->enabled) implementation_->Print("finish");
}

}  // namespace ramag
