# Historical whole-genome experiment variants are no longer product options.
foreach(_retired IN ITEMS
    RAMAG_INTERNAL_COVERAGE_RECOVERY
    RAMAG_INTERNAL_GAP_FILL
    RAMAG_INTERNAL_PRESERVE_LINK_CANDIDATES
    RAMAG_INTERNAL_LINK_PRECHECK
    RAMAG_INTERNAL_GUARDED_RECOVERY
    RAMAG_INTERNAL_POST_SELECTION
    RAMAG_INTERNAL_COVERAGE_DIAGNOSTICS
    RAMAG_INTERNAL_FRAGMENT_LOCAL_EXCHANGE
    RAMAG_INTERNAL_FRAGMENT_SHORT_EXACT
    RAMAG_INTERNAL_SUPPLEMENT_PRUNED_CLUSTERS
    RAMAG_INTERNAL_SUPPLEMENT_ALTERNATIVE_CHAIN
    RAMAG_INTERNAL_STRAND_AWARE_MERGE
    RAMAG_INTERNAL_STRAND_AWARE_DIAGONAL
    RAMAG_INTERNAL_HALF_OPEN_CHAINING
    RAMAG_INTERNAL_SUPPLEMENT_SHORT_CHAIN
    RAMAG_INTERNAL_EXTEND_SELECTED_ENDPOINTS
)
  if(DEFINED ${_retired})
    message(FATAL_ERROR "${_retired} has been retired. Use a fresh build directory without historical experiment options.")
  endif()
endforeach()

# Licensed fixed-end KSW2 and the selected main-path core are always built.
add_library(ramag_pairwise STATIC src/genome/core.cpp src/genome/extension.cpp src/genome/linking.cpp third_party/ksw2/ksw2_extz2_sse.c)
add_library(RaMAG::pairwise ALIAS ramag_pairwise)
set_target_properties(ramag_pairwise PROPERTIES EXPORT_NAME pairwise POSITION_INDEPENDENT_CODE ON)
target_compile_features(ramag_pairwise PUBLIC cxx_std_20)
target_include_directories(ramag_pairwise PRIVATE ${PROJECT_SOURCE_DIR}/third_party/ksw2 ${PROJECT_SOURCE_DIR}/include)
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
install(FILES include/ramag/pairwise_core.hpp include/ramag/alignment.hpp include/ramag/model.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ramag)
install(FILES third_party/ksw2/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/ramag/third_party/ksw2)
install(FILES third_party/kseq/kseq.h DESTINATION ${CMAKE_INSTALL_DATADIR}/ramag/third_party/kseq)
install(FILES third_party/spdlog/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/ramag/third_party/spdlog)
