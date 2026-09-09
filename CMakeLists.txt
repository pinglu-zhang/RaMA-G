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
set(RAMAG_INTERNAL_LCP_ENCODING "byte-coded" CACHE STRING
    "Internal index LCP storage: profile-default, raw, byte-coded")
set_property(CACHE RAMAG_INTERNAL_LCP_ENCODING PROPERTY STRINGS profile-default raw byte-coded)
mark_as_advanced(RAMAG_INTERNAL_LCP_ENCODING)
if(NOT RAMAG_INTERNAL_LCP_ENCODING MATCHES "^(profile-default|raw|byte-coded)$")
  message(FATAL_ERROR "Invalid RAMAG_INTERNAL_LCP_ENCODING: ${RAMAG_INTERNAL_LCP_ENCODING}")
endif()
foreach(obsolete IN ITEMS RAMAG_INTERNAL_ALIGNMENT_CORE RAMAG_INTERNAL_EXTENSION_BACKEND RAMAG_BUILD_PAIRWISE_LIBRARY)
  if(DEFINED ${obsolete})
    message(FATAL_ERROR "${obsolete} has been removed. Use a fresh build without this option; pairwise + KSW2 is now the only production core and is always installed.")
  endif()
endforeach()
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
       "Fetch/validate the pinned sufkit snapshot for production adapters" ON)
set(RAMAG_MUMMER4_BIN_DIR "" CACHE PATH
    "Directory containing original delta-filter/show-coords/show-snps/show-diff/dnadiff")

set(RAMAG_SUFKIT_COMMIT "f8c4c386ee883e45ad0f973efc4c8e1148b0068a"
    CACHE STRING "Required sufkit commit")
set(RAMAG_SUFKIT_SOURCE_DIR "" CACHE PATH "Clean sufkit source override")

# Dependency identities are part of RaMA-G's reproducibility contract.  Do not
# accept branch names, moving tags, abbreviated object names, or mixed-case
# spellings that could make dependency identity ambiguous.
foreach(_ramag_dependency_commit IN ITEMS
    RAMAG_SUFKIT_COMMIT)
  string(LENGTH "${${_ramag_dependency_commit}}"
         _ramag_dependency_commit_length)
  if(NOT _ramag_dependency_commit_length EQUAL 40 OR
     NOT "${${_ramag_dependency_commit}}" MATCHES "^[0-9a-f]+$")
    message(FATAL_ERROR
      "${_ramag_dependency_commit} must be a full 40-character lowercase "
      "Git commit, got '${${_ramag_dependency_commit}}'")
  endif()
endforeach()

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
endif()

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
  src/fasta.cpp
  src/pipeline.cpp
  src/reference_index.cpp
  src/runtime.cpp
  src/logging.cpp
  src/sufkit_adapter.cpp
  src/writers.cpp)
add_library(RaMA-G::core ALIAS ramag_core)

include(cmake/RaMAGPairwiseCore.cmake)
target_link_libraries(ramag_core PRIVATE ramag_pairwise)

target_compile_features(ramag_core PUBLIC cxx_std_20)
find_package(Threads REQUIRED)
target_link_libraries(ramag_core PRIVATE ZLIB::ZLIB Threads::Threads)
target_include_directories(ramag_core SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/third_party/spdlog/include ${PROJECT_SOURCE_DIR}/third_party/kseq)
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
  PRIVATE
    RAMAG_CXX_COMPILER="${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}"
    RAMAG_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
    RAMAG_INTERNAL_LCP_ENCODING="${RAMAG_INTERNAL_LCP_ENCODING}"
    RAMAG_COMPILE_FLAGS="${RAMAG_EFFECTIVE_CXX_FLAGS}"
    RAMAG_OPENMP_ENABLED=${RAMAG_OPENMP_ENABLED_VALUE}
    RAMAG_OPENMP_RUNTIME="${RAMAG_OPENMP_RUNTIME_NAME}"
    RAMAG_OPENMP_SPEC_DATE=${OpenMP_CXX_SPEC_DATE}
    RAMAG_SUFKIT_DIVSUFSORT_OPENMP=${RAMAG_SUFKIT_DIVSUFSORT_OPENMP_VALUE})

# Keep OpenMP target-scoped to RaMA-G.  In particular, do not mutate global
# compiler flags or silently change MUMmer4/other comparison binaries.
if(OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
  target_link_libraries(ramag_core PRIVATE OpenMP::OpenMP_CXX)
endif()

if(RAMAG_ENABLE_PINNED_DEPENDENCIES)
  target_compile_definitions(ramag_core PRIVATE
    RAMAG_HAVE_SUFKIT=1
    RAMAG_DEPENDENCY_MODE="pinned-production-adapters"
    RAMAG_SUFKIT_SOURCE_PATH="${RAMAG_RESOLVED_SUFKIT_SOURCE_DIR}")
  target_link_libraries(ramag_core PRIVATE sufkit::sufkit)
else()
  target_compile_definitions(ramag_core PRIVATE
    RAMAG_HAVE_SUFKIT=0
    RAMAG_DEPENDENCY_MODE="bundled-oracle-baseline"
    RAMAG_SUFKIT_SOURCE_PATH="not-linked")
endif()

ramag_apply_project_warnings(ramag_core)

if(RAMAG_ENABLE_SANITIZERS AND NOT MSVC)
  target_compile_options(ramag_core PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(ramag_core PUBLIC -fsanitize=address,undefined)
endif()

add_executable(ramag src/main.cpp)
target_link_libraries(ramag PRIVATE RaMA-G::core)
ramag_apply_project_warnings(ramag)

install(TARGETS ramag RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES LICENSE THIRD_PARTY_NOTICES.md DESTINATION "${CMAKE_INSTALL_DATADIR}/ramag")

if(BUILD_TESTING AND RAMAG_BUILD_TESTS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit_tests.cpp")
  if(TARGET sufkit::sufkit)
    add_executable(ramag_index_storage_tests tests/index_storage_tests.cpp)
    target_link_libraries(ramag_index_storage_tests PRIVATE RaMA-G::core sufkit::sufkit)
    ramag_apply_project_warnings(ramag_index_storage_tests)
    add_test(NAME ramag.index_storage COMMAND ramag_index_storage_tests)
  endif()
  add_executable(ramag_logging_input_tests tests/logging_input_tests.cpp)
  target_include_directories(ramag_logging_input_tests SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/third_party/spdlog/include)
  target_link_libraries(ramag_logging_input_tests PRIVATE RaMA-G::core Threads::Threads)
  ramag_apply_project_warnings(ramag_logging_input_tests)
  add_test(NAME ramag.logging_input COMMAND ramag_logging_input_tests)
  add_executable(ramag_unified_core_tests tests/unified_core_tests.cpp)
  target_link_libraries(ramag_unified_core_tests PRIVATE RaMA-G::core RaMAG::pairwise)
  ramag_apply_project_warnings(ramag_unified_core_tests)
  add_test(NAME ramag.unified.core COMMAND ramag_unified_core_tests)
  if(TARGET ramag_pairwise AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/pairwise_core_tests.cpp")
    add_library(ramag_pairwise_reference OBJECT tests/pairwise_reference_gap.cpp)
    target_compile_features(ramag_pairwise_reference PRIVATE cxx_std_20)
    target_include_directories(ramag_pairwise_reference PRIVATE third_party/ksw2)
    target_link_libraries(ramag_pairwise_reference PRIVATE RaMAG::pairwise)
    add_executable(ramag_pairwise_core_tests tests/pairwise_core_tests.cpp $<TARGET_OBJECTS:ramag_pairwise_reference>)
    target_link_libraries(ramag_pairwise_core_tests PRIVATE RaMAG::pairwise)
    ramag_apply_project_warnings(ramag_pairwise_core_tests)
    add_test(NAME ramag.pairwise.core COMMAND ramag_pairwise_core_tests)
  endif()
  add_library(ramag_legacy_oracle STATIC tests/legacy_alignment_oracle.cpp)
  target_link_libraries(ramag_legacy_oracle PRIVATE RaMA-G::core)
  target_compile_features(ramag_legacy_oracle PRIVATE cxx_std_20)
  ramag_apply_project_warnings(ramag_legacy_oracle)
  add_executable(ramag_unit_tests tests/unit_tests.cpp)
  target_link_libraries(ramag_unit_tests PRIVATE ramag_legacy_oracle RaMA-G::core)
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

  if(TARGET ramag_pairwise)
    add_executable(ramag_coverage_replay tests/coverage_replay.cpp)
    target_link_libraries(ramag_coverage_replay PRIVATE ramag_pairwise RaMA-G::core)
    target_include_directories(ramag_coverage_replay PRIVATE ${PROJECT_SOURCE_DIR}/third_party/ksw2)
    # This private harness includes the source-attributed kernel, so use the
    # same warning contract as that kernel (without weakening production).
    if(NOT MSVC)
      target_compile_options(ramag_coverage_replay PRIVATE -Wall -Wextra -Wno-unused-function -Wno-unused-parameter)
      if(RAMAG_WARNINGS_AS_ERRORS)
        target_compile_options(ramag_coverage_replay PRIVATE -Werror)
      endif()
    endif()
    add_executable(ramag_pairwise_memory_tests tests/pairwise_memory_tests.cpp)
    target_compile_definitions(ramag_pairwise_memory_tests PRIVATE RAMAG_INTERNAL_COVERAGE_RECOVERY=$<BOOL:${RAMAG_INTERNAL_COVERAGE_RECOVERY}>)
    target_link_libraries(ramag_pairwise_memory_tests PRIVATE ramag_pairwise)
    target_include_directories(ramag_pairwise_memory_tests PRIVATE ${PROJECT_SOURCE_DIR}/third_party/ksw2)
    if(OpenMP_CXX_FOUND)
      target_link_libraries(ramag_pairwise_memory_tests PRIVATE OpenMP::OpenMP_CXX)
    endif()
    add_test(NAME ramag.pairwise.memory COMMAND ramag_pairwise_memory_tests)
  endif()

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
  target_link_libraries(ramag_chaining_tests PRIVATE ramag_legacy_oracle RaMA-G::core)
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
      -DRAMAG_ALIGNMENT_CORE=pairwise
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
