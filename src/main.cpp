#include "ramag/logging.hpp"
#include <memory>
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

#if defined(__linux__) && defined(__ELF__)
#include <cerrno>
#include <cstring>
#include <sched.h>
namespace {
// ELF preinit runs before shared-library constructors, including libgomp.
// Only POD state and the affinity system interface are used here.
cpu_set_t launch_cpu_mask;
int launch_cpu_status = -1;
int launch_cpu_error = 0;
__attribute__((no_sanitize("address", "undefined")))
void CaptureLaunchCpuMask() {
  launch_cpu_status = sched_getaffinity(0, sizeof(launch_cpu_mask),
                                       &launch_cpu_mask);
  if (launch_cpu_status != 0) launch_cpu_error = errno;
}
using PreinitFunction = void (*)();
__attribute__((section(".preinit_array"), used))
PreinitFunction capture_launch_cpu_mask = CaptureLaunchCpuMask;
}
#endif

namespace {
ramag::CpuAffinityInfo LaunchCpuAffinity() {
#if defined(__linux__) && defined(__ELF__)
  if (launch_cpu_status != 0) {
    throw ramag::CliError("cannot capture pre-runtime CPU affinity: " +
                         std::string(std::strerror(launch_cpu_error)));
  }
  ramag::CpuAffinityInfo info;
  info.supported = true;
  info.source = "elf-preinit-sched-affinity";
  for (std::uint32_t cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (CPU_ISSET(static_cast<int>(cpu), &launch_cpu_mask)) {
      info.logical_cpus.push_back(cpu);
    }
  }
  return info;
#else
  return ramag::CurrentCpuAffinity();
#endif
}
}

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
#if defined(__linux__) && defined(__ELF__)
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
  std::unique_ptr<ramag::RunLogger> logger;
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
      const auto launch_affinity = LaunchCpuAffinity();
      ramag::ConfigureIndexRuntime(parsed.index_spec, launch_affinity);
      if (parsed.print_effective_config) {
        std::cout << ramag::EffectiveIndexConfigText(parsed.index_spec);
        return 0;
      }
      logger = std::make_unique<ramag::RunLogger>(parsed.index_spec.work_dir);
      logger->Info(ramag::VersionText(), false);
      logger->Info(ramag::EffectiveIndexConfigText(parsed.index_spec), false);
      logger->Info("invocation=" + Invocation(argc, argv), false);
      const auto paths = ramag::RunReferenceIndexPipeline(
          parsed.index_spec, Invocation(argc, argv), binary_path,
          launch_affinity, logger.get());
      logger->Info("status=success exit_code=0 index=" + paths.index.string());
      logger->Flush();
      return 0;
    }
    ramag::ConfigureOpenMpRuntime(parsed.run_spec);
#if defined(__linux__)
    // Restore only the executable's pre-libgomp authorized CPU mask. The
    // reusable library never changes its caller's affinity or signal handlers.
    (void)LaunchCpuAffinity();
    if (static_cast<unsigned>(CPU_COUNT(&launch_cpu_mask)) < parsed.run_spec.threads) {
      throw ramag::CliError("pairwise core: requested threads exceed launch allowed CPUs");
    }
    if (sched_setaffinity(0, sizeof(launch_cpu_mask), &launch_cpu_mask) != 0) {
      throw ramag::CliError("pairwise core: could not restore launch CPU affinity");
    }
#endif
    if (parsed.print_effective_config) {
      std::cout << ramag::EffectiveConfigText(parsed.run_spec);
      return 0;
    }
    logger = std::make_unique<ramag::RunLogger>(parsed.run_spec.work_dir);
    logger->Info(ramag::VersionText(), false);
    logger->Info(ramag::EffectiveConfigText(parsed.run_spec), false);
    logger->Info("invocation=" + Invocation(argc, argv), false);
    const auto outcome = ramag::RunAlignmentPipeline(
        parsed.run_spec, Invocation(argc, argv), binary_path, logger.get());
    (void)outcome;
    logger->Flush();
    return 0;
  } catch (const std::exception& error) {
    const int code = ramag::FailureExitCode(error);
    const std::string message = "status=" + std::string(dynamic_cast<const ramag::InterruptedError*>(&error) ? "interrupted" : "failed") + " exit_code=" + std::to_string(code) + " message=" + error.what();
    try {
      if (logger) logger->Error(message);
      else std::cerr << "ramag: " << message << '\n';
    } catch (...) { std::cerr << "ramag: " << message << '\n'; }
    return code;
  }
}
