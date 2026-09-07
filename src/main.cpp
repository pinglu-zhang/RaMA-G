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
      const auto launch_affinity = LaunchCpuAffinity();
      ramag::ConfigureIndexRuntime(parsed.index_spec, launch_affinity);
      if (parsed.print_effective_config) {
        std::cout << ramag::EffectiveIndexConfigText(parsed.index_spec);
        return 0;
      }
      const auto paths = ramag::RunReferenceIndexPipeline(
          parsed.index_spec, Invocation(argc, argv), binary_path,
          launch_affinity);
      std::cout << "RaMA-G reference index completed: index=" << paths.index
                << "; marker=" << paths.complete << '\n';
      return 0;
    }
    ramag::ConfigureOpenMpRuntime(parsed.run_spec);
#if RAMAG_USE_PAIRWISE_CORE && defined(__linux__)
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
