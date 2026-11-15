cmake_minimum_required(VERSION 3.22)

project(
  RaMA-G
  VERSION 0.1.1
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

set(RAMAG_SUFKIT_COMMIT "028075e6f2d622fcbdcf76b153bbde0f069ca64f"
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
  src/genome/alignment.cpp
  src/app/cli.cpp
  src/io/fasta.cpp
  src/app/pipeline.cpp
  src/index/reference_index.cpp
  src/app/runtime.cpp

  src/index/sufkit_adapter.cpp
  src/io/writers.cpp)
add_library(RaMA-G::core ALIAS ramag_core)

include(cmake/RaMAGPairwiseCore.cmake)
target_link_libraries(ramag_core PRIVATE ramag_pairwise)

target_compile_features(ramag_core PUBLIC cxx_std_20)
find_package(Threads REQUIRED)
target_link_libraries(ramag_core PRIVATE ZLIB::ZLIB Threads::Threads)
include(cmake/RaMACentromere.cmake)
target_include_directories(ramag_core PRIVATE "${PROJECT_SOURCE_DIR}/include")
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

add_executable(ramag src/app/main.cpp)
target_link_libraries(ramag PRIVATE RaMA-G::core)
target_include_directories(ramag PRIVATE "${PROJECT_SOURCE_DIR}/include")
ramag_apply_project_warnings(ramag)

install(TARGETS ramag RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES LICENSE THIRD_PARTY_NOTICES.md DESTINATION "${CMAKE_INSTALL_DATADIR}/ramag")

if(BUILD_TESTING AND RAMAG_BUILD_TESTS)
  # Public synthetic checks must not depend on internal fixture availability.
  add_subdirectory(tests/public)
  if(EXISTS "${PROJECT_SOURCE_DIR}/cmake/RaMAInternalTests.cmake")
    include(cmake/RaMAInternalTests.cmake)
  endif()
endif()
