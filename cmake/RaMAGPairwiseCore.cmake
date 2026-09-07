# Licensed dependencies live in third_party/; no network or Rust is
# required by this library. No Sufkit/SeqPro/graph dependency is exposed.
add_library(ramag_pairwise STATIC src/pairwise_core.cpp third_party/ksw2/ksw2_extz2_sse.c)
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
install(FILES third_party/attribution/LICENSE.pairwise third_party/ksw2/LICENSE third_party/attribution/pairwise-source.json DESTINATION ${CMAKE_INSTALL_DATADIR}/ramag/third_party)
