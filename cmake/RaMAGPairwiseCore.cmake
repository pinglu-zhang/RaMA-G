# Licensed dependencies live in third_party/; no network or Rust is
# required by this library. No Sufkit/SeqPro/graph dependency is exposed.
add_library(ramag_pairwise STATIC src/pairwise_core.cpp third_party/ksw2/ksw2_extz2_sse.c)
option(RAMAG_INTERNAL_COVERAGE_RECOVERY "Enable experimental residual recovery in one-to-one CLI" OFF)
set(RAMAG_INTERNAL_GAP_FILL "off" CACHE STRING "Experimental guarded gap fill: off|exact-gap|ksw2-gap")
set_property(CACHE RAMAG_INTERNAL_GAP_FILL PROPERTY STRINGS off exact-gap ksw2-gap)
set(_ramag_gap_modes off exact-gap ksw2-gap)
list(FIND _ramag_gap_modes "${RAMAG_INTERNAL_GAP_FILL}" _ramag_gap_mode)
if(_ramag_gap_mode EQUAL -1)
  message(FATAL_ERROR "RAMAG_INTERNAL_GAP_FILL must be off|exact-gap|ksw2-gap")
endif()
if(_ramag_gap_mode GREATER 0 AND NOT RAMAG_INTERNAL_COVERAGE_RECOVERY)
  message(FATAL_ERROR "Gap fill requires RAMAG_INTERNAL_COVERAGE_RECOVERY=ON")
endif()
mark_as_advanced(RAMAG_INTERNAL_GAP_FILL)
target_compile_definitions(ramag_pairwise PRIVATE RAMAG_INTERNAL_GAP_FILL=${_ramag_gap_mode})
mark_as_advanced(RAMAG_INTERNAL_COVERAGE_RECOVERY)
target_compile_definitions(ramag_pairwise PRIVATE
  RAMAG_INTERNAL_COVERAGE_RECOVERY=$<BOOL:${RAMAG_INTERNAL_COVERAGE_RECOVERY}>)
option(RAMAG_INTERNAL_COVERAGE_DIAGNOSTICS "Enable private coverage capture in CLI bridge" OFF)
mark_as_advanced(RAMAG_INTERNAL_COVERAGE_DIAGNOSTICS)
target_compile_definitions(ramag_pairwise PRIVATE
  RAMAG_INTERNAL_COVERAGE_DIAGNOSTICS=$<BOOL:${RAMAG_INTERNAL_COVERAGE_DIAGNOSTICS}>)
add_library(RaMAG::pairwise ALIAS ramag_pairwise)
set_target_properties(ramag_pairwise PROPERTIES EXPORT_NAME pairwise POSITION_INDEPENDENT_CODE ON)
target_compile_features(ramag_pairwise PUBLIC cxx_std_20)
target_include_directories(ramag_pairwise PRIVATE ${PROJECT_SOURCE_DIR}/third_party/ksw2)
target_include_directories(ramag_pairwise PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(OpenMP_CXX_FOUND)
  target_link_libraries(ramag_pairwise PRIVATE OpenMP::OpenMP_CXX)
endif()
if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i.86)$")
  message(FATAL_ERROR "The frozen pairwise KSW2 kernel requires x86 SSE2")
endif()
if(NOT MSVC)
  # The source-attributed upstream algorithms intentionally retain their local
  # names and integer interfaces. Unused oracle helpers are retained for tests.
  target_compile_options(ramag_pairwise PRIVATE -Wall -Wextra -Wno-unused-function -Wno-unused-parameter -msse2)
  if(RAMAG_WARNINGS_AS_ERRORS)
    target_compile_options(ramag_pairwise PRIVATE -Werror)
  endif()
  if(RAMAG_ENABLE_SANITIZERS)
    target_compile_options(ramag_pairwise PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(ramag_pairwise PUBLIC -fsanitize=address,undefined)
  endif()
endif()
include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/RaMAGPairwiseConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/RaMAGPairwiseConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RaMAGPairwise)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/RaMAGPairwiseConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(TARGETS ramag_pairwise EXPORT RaMAGPairwiseTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(EXPORT RaMAGPairwiseTargets NAMESPACE RaMAG:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RaMAGPairwise)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/RaMAGPairwiseConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/RaMAGPairwiseConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RaMAGPairwise)
install(FILES include/ramag/pairwise_core.hpp include/ramag/alignment.hpp include/ramag/model.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ramag)
install(FILES third_party/ksw2/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/ramag/third_party/ksw2)
install(FILES third_party/kseq/kseq.h DESTINATION ${CMAKE_INSTALL_DATADIR}/ramag/third_party/kseq)
install(FILES third_party/spdlog/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/ramag/third_party/spdlog)
