cmake_minimum_required(VERSION 3.22)

project(
  RaMA-G
  VERSION 0.1.0
  DESCRIPTION "Deterministic pairwise whole-genome alignment"
  LANGUAGES C CXX)

string(TOUPPER "${CMAKE_BUILD_TYPE}" RAMAG_BUILD_TYPE_UPPER)
set(RAMAG_EFFECTIVE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${RAMAG_BUILD_TYPE_UPPER}}")

include(CTest)
include(GNUInstallDirs)
find_package(ZLIB REQUIRED)

# Linux x86-64 is the v1 production platform.  OpenMP is therefore a required
# RaMA-G build dependency there: accepting --threads while silently compiling
# pragmas away would make both the product contract and benchmark provenance
# misleading.  Other platforms may still build the scalar portability path.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  find_package(OpenMP REQUIRED COMPONENTS C CXX)
else()
  find_package(OpenMP QUIET COMPONENTS C CXX)
endif()
if(OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
  set(RAMAG_OPENMP_ENABLED_VALUE 1)
  string(REPLACE ";" "," RAMAG_OPENMP_RUNTIME_NAME
         "${OpenMP_CXX_LIB_NAMES}")
  if(RAMAG_OPENMP_RUNTIME_NAME STREQUAL "")
    set(RAMAG_OPENMP_RUNTIME_NAME "${CMAKE_CXX_COMPILER_ID}-openmp")
  endif()
  if(NOT OpenMP_CXX_SPEC_DATE)
    set(OpenMP_CXX_SPEC_DATE 0)
  endif()
  string(APPEND RAMAG_EFFECTIVE_CXX_FLAGS
         " openmp=${OpenMP_CXX_FLAGS}")
else()
  set(RAMAG_OPENMP_ENABLED_VALUE 0)
  set(RAMAG_OPENMP_RUNTIME_NAME "unavailable")
  set(OpenMP_CXX_SPEC_DATE 0)
endif()

option(RAMAG_BUILD_TESTS "Build RaMA-G tests" ON)
option(RAMAG_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(RAMAG_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(RAMAG_ENABLE_PINNED_DEPENDENCIES
       "Fetch/validate the pinned sufkit and SeqPro snapshots for production adapters" ON)
set(RAMAG_MUMMER4_BIN_DIR "" CACHE PATH
    "Directory containing original delta-filter/show-coords/show-snps/show-diff/dnadiff")

# Compile-time-only experimental gap-extension selection.  The default scalar
# build deliberately resolves neither KSW2 nor Rust/Block Aligner, so ordinary
# RaMA-G users retain the existing dependency and runtime contract.
set(RAMAG_INTERNAL_EXTENSION_BACKEND "scalar" CACHE STRING
    "Internal extension backend: scalar, ksw2-exact, ksw2-band-auto, block-exact, or block-adaptive")
set_property(CACHE RAMAG_INTERNAL_EXTENSION_BACKEND PROPERTY STRINGS
             scalar ksw2-exact ksw2-band-auto block-exact block-adaptive)
set(RAMAG_KSW2_COMMIT "289609bd9e5381a13b16239d0a7703f1ff03f9ca"
    CACHE STRING "Required KSW2 commit for experimental extension backends")
set(RAMAG_KSW2_SOURCE_DIR "" CACHE PATH "Clean KSW2 source override")
set(RAMAG_BLOCK_ALIGNER_COMMIT "4fcf630cf775de5b578fe63971f210e1dc958791"
    CACHE STRING "Required Block Aligner commit for experimental extension backends")
set(RAMAG_BLOCK_ALIGNER_SOURCE_DIR "" CACHE PATH "Clean Block Aligner source override")
set(RAMAG_BLOCK_ALIGNER_SIMD "sse2" CACHE STRING
    "Block Aligner SIMD implementation for internal experiments: sse2 or avx2")
set_property(CACHE RAMAG_BLOCK_ALIGNER_SIMD PROPERTY STRINGS sse2 avx2)
set(RAMAG_BLOCK_ALIGNER_RUST_TOOLCHAIN "nightly-2024-06-22" CACHE STRING
    "Pinned Rust toolchain used for the experimental Block Aligner C ABI")
set(RAMAG_BLOCK_ALIGNER_RUSTUP_HOME "" CACHE PATH
    "Isolated RUSTUP_HOME containing the pinned Block Aligner toolchain")
set(RAMAG_BLOCK_ALIGNER_CARGO_HOME "" CACHE PATH
    "Optional isolated CARGO_HOME for Block Aligner (defaults inside the build tree)")
option(RAMAG_INTERNAL_EXTENSION_METRICS
       "Collect per-call experimental extension timing and work estimates" OFF)
option(RAMAG_BUILD_EXTENSION_BENCHMARK
       "Build the internal deterministic gap-extension benchmark" OFF)
mark_as_advanced(
  RAMAG_INTERNAL_EXTENSION_BACKEND
  RAMAG_KSW2_COMMIT RAMAG_KSW2_SOURCE_DIR
  RAMAG_BLOCK_ALIGNER_COMMIT RAMAG_BLOCK_ALIGNER_SOURCE_DIR
  RAMAG_BLOCK_ALIGNER_SIMD RAMAG_BLOCK_ALIGNER_RUST_TOOLCHAIN
  RAMAG_BLOCK_ALIGNER_RUSTUP_HOME RAMAG_BLOCK_ALIGNER_CARGO_HOME
  RAMAG_INTERNAL_EXTENSION_METRICS RAMAG_BUILD_EXTENSION_BENCHMARK)

set(RAMAG_SUFKIT_COMMIT "bdb67c6de5daddd8a005640de73d96549d2575f4"
    CACHE STRING "Required sufkit commit")
set(RAMAG_SEQPRO_COMMIT "6781cadcf81a0da53d7573444594c1484947017c"
    CACHE STRING "Required SeqPro commit")
set(RAMAG_SUFKIT_SOURCE_DIR "" CACHE PATH "Clean sufkit source override")
set(RAMAG_SEQPRO_SOURCE_DIR "" CACHE PATH "Clean SeqPro source override")

# Dependency identities are part of RaMA-G's reproducibility contract.  Do not
# accept branch names, moving tags, abbreviated object names, or mixed-case
# spellings that could make manifest/string equality ambiguous.
foreach(_ramag_dependency_commit IN ITEMS
    RAMAG_SUFKIT_COMMIT RAMAG_SEQPRO_COMMIT
    RAMAG_KSW2_COMMIT RAMAG_BLOCK_ALIGNER_COMMIT)
  string(LENGTH "${${_ramag_dependency_commit}}"
         _ramag_dependency_commit_length)
  if(NOT _ramag_dependency_commit_length EQUAL 40 OR
     NOT "${${_ramag_dependency_commit}}" MATCHES "^[0-9a-f]+$")
    message(FATAL_ERROR
      "${_ramag_dependency_commit} must be a full 40-character lowercase "
      "Git commit, got '${${_ramag_dependency_commit}}'")
  endif()
endforeach()

set(_ramag_extension_backends
    scalar ksw2-exact ksw2-band-auto block-exact block-adaptive)
if(NOT RAMAG_INTERNAL_EXTENSION_BACKEND IN_LIST _ramag_extension_backends)
  message(FATAL_ERROR
    "RAMAG_INTERNAL_EXTENSION_BACKEND must be one of: ${_ramag_extension_backends}; "
    "got '${RAMAG_INTERNAL_EXTENSION_BACKEND}'")
endif()
if(NOT RAMAG_BLOCK_ALIGNER_SIMD STREQUAL "sse2" AND
   NOT RAMAG_BLOCK_ALIGNER_SIMD STREQUAL "avx2")
  message(FATAL_ERROR
    "RAMAG_BLOCK_ALIGNER_SIMD must be sse2 or avx2, got "
    "'${RAMAG_BLOCK_ALIGNER_SIMD}'")
endif()

if(MSVC)
  string(APPEND RAMAG_EFFECTIVE_CXX_FLAGS " cxx_std_20 /W4 /permissive-")
  if(RAMAG_WARNINGS_AS_ERRORS)
    string(APPEND RAMAG_EFFECTIVE_CXX_FLAGS " /WX")
  endif()
else()
  string(APPEND RAMAG_EFFECTIVE_CXX_FLAGS
         " cxx_std_20 -Wall -Wextra -Wpedantic -Wconversion -Wshadow")
  if(RAMAG_WARNINGS_AS_ERRORS)
    string(APPEND RAMAG_EFFECTIVE_CXX_FLAGS " -Werror")
  endif()
  if(RAMAG_ENABLE_SANITIZERS)
    string(APPEND RAMAG_EFFECTIVE_CXX_FLAGS
           " -fsanitize=address,undefined -fno-omit-frame-pointer")
  endif()
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(RaMAGitSnapshot)

ramag_resolve_git_revision(RAMAG_GIT_COMMIT "${CMAKE_CURRENT_SOURCE_DIR}")

if(RAMAG_ENABLE_PINNED_DEPENDENCIES)
  include(RaMADependencies)
  ramag_configure_pinned_dependencies()
  FetchContent_GetProperties(sufkit SOURCE_DIR RAMAG_RESOLVED_SUFKIT_SOURCE_DIR)
  FetchContent_GetProperties(seqpro SOURCE_DIR RAMAG_RESOLVED_SEQPRO_SOURCE_DIR)
endif()

include(RaMAExtensionBackends)
ramag_configure_extension_backend()

function(ramag_apply_project_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive-)
    if(RAMAG_WARNINGS_AS_ERRORS)
      target_compile_options(${target} PRIVATE /WX)
    endif()
  else()
    target_compile_options(${target} PRIVATE
      -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    if(RAMAG_WARNINGS_AS_ERRORS)
      target_compile_options(${target} PRIVATE -Werror)
    endif()
  endif()
endfunction()

set(RAMAG_SUFKIT_DIVSUFSORT_OPENMP_VALUE 0)
if(RAMAG_ENABLE_PINNED_DEPENDENCIES AND OpenMP_C_FOUND)
  foreach(_ramag_divsufsort_target IN ITEMS
      sufkit_divsufsort32 sufkit_divsufsort64)
    if(TARGET ${_ramag_divsufsort_target})
      # These object targets are part of the pinned sufkit build.  Applying the
      # C OpenMP imported target activates bundled libdivsufsort's guarded
      # _OPENMP implementation without changing the sufkit source checkout.
      target_link_libraries(${_ramag_divsufsort_target}
                            PRIVATE OpenMP::OpenMP_C)
      set(RAMAG_SUFKIT_DIVSUFSORT_OPENMP_VALUE 1)
    endif()
  endforeach()
endif()

add_library(ramag_core
  src/alignment.cpp
  src/cli.cpp
  src/extension_backend.cpp
  src/fasta.cpp
  src/manifest.cpp
  src/pipeline.cpp
  src/reference_index.cpp
  src/runtime.cpp
  src/seqpro_adapter.cpp
  src/sha256.cpp
  src/sufkit_adapter.cpp
  src/writers.cpp)
add_library(RaMA-G::core ALIAS ramag_core)

target_compile_features(ramag_core PUBLIC cxx_std_20)
target_link_libraries(ramag_core PRIVATE ZLIB::ZLIB)
target_include_directories(
  ramag_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_definitions(
  ramag_core
  PUBLIC
    RAMAG_VERSION="${PROJECT_VERSION}"
    RAMAG_GIT_COMMIT="${RAMAG_GIT_COMMIT}"
    RAMAG_SUFKIT_COMMIT="${RAMAG_SUFKIT_COMMIT}"
    RAMAG_SEQPRO_COMMIT="${RAMAG_SEQPRO_COMMIT}"
    RAMAG_EXTENSION_BACKEND="${RAMAG_INTERNAL_EXTENSION_BACKEND}"
    RAMAG_KSW2_COMMIT="${RAMAG_KSW2_COMMIT}"
    RAMAG_BLOCK_ALIGNER_COMMIT="${RAMAG_BLOCK_ALIGNER_COMMIT}"
  PRIVATE
    RAMAG_CXX_COMPILER="${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}"
    RAMAG_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
    RAMAG_COMPILE_FLAGS="${RAMAG_EFFECTIVE_CXX_FLAGS}"
    RAMAG_OPENMP_ENABLED=${RAMAG_OPENMP_ENABLED_VALUE}
    RAMAG_OPENMP_RUNTIME="${RAMAG_OPENMP_RUNTIME_NAME}"
    RAMAG_OPENMP_SPEC_DATE=${OpenMP_CXX_SPEC_DATE}
    RAMAG_SUFKIT_DIVSUFSORT_OPENMP=${RAMAG_SUFKIT_DIVSUFSORT_OPENMP_VALUE})

if(RAMAG_INTERNAL_EXTENSION_METRICS OR
   NOT RAMAG_INTERNAL_EXTENSION_BACKEND STREQUAL "scalar")
  target_compile_definitions(ramag_core PRIVATE RAMAG_EXTENSION_METRICS=1)
else()
  target_compile_definitions(ramag_core PRIVATE RAMAG_EXTENSION_METRICS=0)
endif()

ramag_link_extension_backend(ramag_core)

# Keep OpenMP target-scoped to RaMA-G.  In particular, do not mutate global
# compiler flags or silently change MUMmer4/other comparison binaries.
if(OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
  target_link_libraries(ramag_core PRIVATE OpenMP::OpenMP_CXX)
endif()

if(RAMAG_ENABLE_PINNED_DEPENDENCIES)
  target_compile_definitions(ramag_core PRIVATE
    RAMAG_HAVE_SUFKIT=1
    RAMAG_HAVE_SEQPRO=1
    RAMAG_DEPENDENCY_MODE="pinned-production-adapters"
    RAMAG_SUFKIT_SOURCE_PATH="${RAMAG_RESOLVED_SUFKIT_SOURCE_DIR}"
    RAMAG_SEQPRO_SOURCE_PATH="${RAMAG_RESOLVED_SEQPRO_SOURCE_DIR}")
  target_link_libraries(ramag_core PRIVATE sufkit::sufkit SeqPro::seqpro)
else()
  target_compile_definitions(ramag_core PRIVATE
    RAMAG_HAVE_SUFKIT=0
    RAMAG_HAVE_SEQPRO=0
    RAMAG_DEPENDENCY_MODE="bundled-oracle-baseline"
    RAMAG_SUFKIT_SOURCE_PATH="not-linked"
    RAMAG_SEQPRO_SOURCE_PATH="not-linked")
endif()

ramag_apply_project_warnings(ramag_core)

if(RAMAG_ENABLE_SANITIZERS AND NOT MSVC)
  target_compile_options(ramag_core PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(ramag_core PUBLIC -fsanitize=address,undefined)
endif()

add_executable(ramag src/main.cpp)
target_link_libraries(ramag PRIVATE RaMA-G::core)
ramag_apply_project_warnings(ramag)

if(RAMAG_BUILD_EXTENSION_BENCHMARK)
  add_executable(ramag_extension_benchmark
                 benchmarks/extension_benchmark.cpp)
  target_link_libraries(ramag_extension_benchmark PRIVATE RaMA-G::core)
  ramag_apply_project_warnings(ramag_extension_benchmark)
endif()

install(TARGETS ramag RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES LICENSE THIRD_PARTY_NOTICES.md DESTINATION "${CMAKE_INSTALL_DATADIR}/ramag")

if(BUILD_TESTING AND RAMAG_BUILD_TESTS AND
   NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
  message(STATUS
    "RaMA-G test sources are not included; skipping test targets")
endif()

if(BUILD_TESTING AND RAMAG_BUILD_TESTS AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
  add_executable(ramag_unit_tests tests/unit_tests.cpp)
  target_link_libraries(ramag_unit_tests PRIVATE RaMA-G::core)
  ramag_apply_project_warnings(ramag_unit_tests)
  add_test(NAME ramag.unit COMMAND ramag_unit_tests)

  add_test(
    NAME ramag.openmp.config
    COMMAND "${CMAKE_COMMAND}"
      -DRAMAG_EXE=$<TARGET_FILE:ramag>
      -DOPENMP_C_FOUND=${OpenMP_C_FOUND}
      -DOPENMP_CXX_FOUND=${OpenMP_CXX_FOUND}
      -DEXPECTED_OPENMP_RUNTIME=${RAMAG_OPENMP_RUNTIME_NAME}
      -DEXPECTED_OPENMP_SPEC_DATE=${OpenMP_CXX_SPEC_DATE}
      -DEXPECTED_DIVSUFSORT_OPENMP=${RAMAG_SUFKIT_DIVSUFSORT_OPENMP_VALUE}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/openmp_config.cmake)
  add_test(
    NAME ramag.openmp.required
    COMMAND "${CMAKE_COMMAND}"
      -DSOURCE_ROOT=${CMAKE_CURRENT_SOURCE_DIR}
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/openmp-required-configure
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/openmp_required.cmake)

  add_executable(ramag_seed_oracle_tests tests/seed_oracle_tests.cpp)
  target_link_libraries(ramag_seed_oracle_tests PRIVATE RaMA-G::core)
  ramag_apply_project_warnings(ramag_seed_oracle_tests)
  add_test(NAME ramag.seed_oracle COMMAND ramag_seed_oracle_tests)

  add_executable(ramag_adapter_differential_tests
                 tests/adapter_differential_tests.cpp)
  target_link_libraries(ramag_adapter_differential_tests PRIVATE RaMA-G::core)
  ramag_apply_project_warnings(ramag_adapter_differential_tests)
  add_test(NAME ramag.adapter_differential COMMAND ramag_adapter_differential_tests)

  add_executable(ramag_chaining_tests tests/chaining_tests.cpp)
  target_link_libraries(ramag_chaining_tests PRIVATE RaMA-G::core)
  ramag_apply_project_warnings(ramag_chaining_tests)
  if(RAMAG_ENABLE_SANITIZERS)
    add_test(NAME ramag.chaining COMMAND ramag_chaining_tests --skip-scale)
  else()
    add_test(NAME ramag.chaining COMMAND ramag_chaining_tests)
  endif()
  set_tests_properties(ramag.chaining PROPERTIES
    TIMEOUT 180
    LABELS "chaining;scale")

  foreach(gate_mode IN ITEMS clean dirty wrong-commit abbreviated-commit missing)
    add_test(
      NAME ramag.dependency_gate.${gate_mode}
      COMMAND "${CMAKE_COMMAND}"
        -DGATE_MODE=${gate_mode}
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/dependency-gate-${gate_mode}
        -DMODULE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/cmake/RaMAGitSnapshot.cmake
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/dependency_override_gate.cmake)
  endforeach()
  set_tests_properties(
    ramag.dependency_gate.dirty
    ramag.dependency_gate.wrong-commit
    ramag.dependency_gate.abbreviated-commit
    ramag.dependency_gate.missing
    PROPERTIES WILL_FAIL TRUE)

  add_test(
    NAME ramag.integration.smoke
    COMMAND "${CMAKE_COMMAND}"
      -DRAMAG_EXE=$<TARGET_FILE:ramag>
      -DEXPECTED_SUFKIT_COMMIT=${RAMAG_SUFKIT_COMMIT}
      -DEXPECTED_DIVSUFSORT_OPENMP=${RAMAG_SUFKIT_DIVSUFSORT_OPENMP_VALUE}
      -DEXPECT_PINNED_ADAPTERS=${RAMAG_ENABLE_PINNED_DEPENDENCIES}
      -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/integration-smoke
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration_smoke.cmake)

  if(UNIX)
    add_test(
      NAME ramag.integration.signals
      COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/signal_smoke.sh
        $<TARGET_FILE:ramag>
        ${CMAKE_CURRENT_BINARY_DIR}/signal-smoke)
  endif()

  find_program(RAMAG_SAMTOOLS_EXECUTABLE samtools)
  if(RAMAG_SAMTOOLS_EXECUTABLE)
    add_test(
      NAME ramag.integration.samtools
      COMMAND "${RAMAG_SAMTOOLS_EXECUTABLE}" view -h -o /dev/null
              ${CMAKE_CURRENT_BINARY_DIR}/integration-smoke/result.sam)
    set_tests_properties(ramag.integration.samtools PROPERTIES
      DEPENDS ramag.integration.smoke)
  endif()

  if(RAMAG_MUMMER4_BIN_DIR)
    foreach(tool IN ITEMS delta-filter show-coords show-snps show-diff dnadiff mummer)
      if(NOT EXISTS "${RAMAG_MUMMER4_BIN_DIR}/${tool}")
        message(FATAL_ERROR "Missing MUMmer4 compatibility tool: ${RAMAG_MUMMER4_BIN_DIR}/${tool}")
      endif()
    endforeach()
    add_test(
      NAME ramag.integration.mummer4
      COMMAND "${CMAKE_COMMAND}"
        -DRAMAG_EXE=$<TARGET_FILE:ramag>
        -DMUMMER_BIN_DIR=${RAMAG_MUMMER4_BIN_DIR}
        -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/mummer4-compat
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/mummer4_compat.cmake)
    # Local libtool `.libs` utilities depend on the sibling uninstalled
    # libumdmummer shared object.  Make that runtime dependency explicit for
    # the compatibility test instead of relying on the caller's shell.
    set_tests_properties(ramag.integration.mummer4 PROPERTIES
      ENVIRONMENT_MODIFICATION
        "LD_LIBRARY_PATH=path_list_prepend:${RAMAG_MUMMER4_BIN_DIR}")

    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(
      NAME ramag.integration.mummer4.mum_oracle
      COMMAND "${Python3_EXECUTABLE}"
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/mummer4_mum_oracle.py
        --mummer ${RAMAG_MUMMER4_BIN_DIR}/mummer
        --test-root ${CMAKE_CURRENT_BINARY_DIR}/mummer4-mum-oracle)
    set_tests_properties(ramag.integration.mummer4.mum_oracle PROPERTIES
      ENVIRONMENT_MODIFICATION
        "LD_LIBRARY_PATH=path_list_prepend:${RAMAG_MUMMER4_BIN_DIR}")
  endif()
endif()
