include(FetchContent)

function(_ramag_populate_clean_dependency name repository commit override label)
  if(override)
    ramag_require_clean_commit("${override}" "${commit}" "${label}")
    string(TOUPPER "${name}" upper_name)
    set("FETCHCONTENT_SOURCE_DIR_${upper_name}" "${override}")
  endif()

  FetchContent_Declare(
    ${name}
    GIT_REPOSITORY "${repository}"
    GIT_TAG "${commit}"
    GIT_SHALLOW FALSE
    GIT_PROGRESS TRUE)
  FetchContent_GetProperties(${name})
  if(NOT ${name}_POPULATED)
    FetchContent_Populate(${name})
  endif()
  FetchContent_GetProperties(${name} SOURCE_DIR resolved_source)
  ramag_require_clean_commit(
    "${resolved_source}" "${commit}" "resolved ${label}")
  set("RAMAG_RESOLVED_${name}_SOURCE_DIR" "${resolved_source}" PARENT_SCOPE)
endfunction()

function(ramag_configure_extension_backend)
  set(RAMAG_EXTENSION_TARGET "" PARENT_SCOPE)
  set(RAMAG_EXTENSION_SOURCE "not-linked" PARENT_SCOPE)

  if(RAMAG_INTERNAL_EXTENSION_BACKEND MATCHES "^ksw2-")
    if(RAMAG_INTERNAL_EXTENSION_BACKEND STREQUAL "ksw2-exact" AND
       NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64)$")
      message(FATAL_ERROR "The pinned KSW2 SSE backend requires x86-64")
    endif()
    _ramag_populate_clean_dependency(
      ksw2 https://github.com/lh3/ksw2.git "${RAMAG_KSW2_COMMIT}"
      "${RAMAG_KSW2_SOURCE_DIR}" "KSW2")
    if(RAMAG_INTERNAL_EXTENSION_BACKEND STREQUAL "ksw2-exact")
      set(ksw2_implementation
          "${RAMAG_RESOLVED_ksw2_SOURCE_DIR}/ksw2_extz2_sse.c")
    else()
      # extz2_sse's score and traceback can disagree when a global optimum
      # leaves a narrow fixed band.  KSW2's standard global implementation
      # enforces the band in both DP and traceback and therefore supplies the
      # correctness-first automatic-band route.
      set(ksw2_implementation
          "${RAMAG_RESOLVED_ksw2_SOURCE_DIR}/ksw2_gg.c")
    endif()
    add_library(ramag_ksw2 STATIC "${ksw2_implementation}")
    target_include_directories(ramag_ksw2 SYSTEM PUBLIC
      "${RAMAG_RESOLVED_ksw2_SOURCE_DIR}")
    if(RAMAG_INTERNAL_EXTENSION_BACKEND STREQUAL "ksw2-exact")
      target_compile_definitions(ramag_ksw2 PRIVATE KSW_SSE2_ONLY=1)
    endif()
    if(RAMAG_ENABLE_SANITIZERS AND NOT MSVC)
      target_compile_options(ramag_ksw2 PRIVATE
        -fsanitize=address,undefined -fno-omit-frame-pointer)
    endif()
    set_target_properties(ramag_ksw2 PROPERTIES POSITION_INDEPENDENT_CODE ON)
    set(RAMAG_EXTENSION_TARGET ramag_ksw2 PARENT_SCOPE)
    set(RAMAG_EXTENSION_SOURCE
        "${RAMAG_RESOLVED_ksw2_SOURCE_DIR}" PARENT_SCOPE)
    return()
  endif()

  if(RAMAG_INTERNAL_EXTENSION_BACKEND MATCHES "^block-")
    if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64)$")
      message(FATAL_ERROR "The configured Block Aligner SIMD backend requires x86-64")
    endif()
    find_program(RAMAG_CARGO_EXECUTABLE NAMES cargo REQUIRED)
    get_filename_component(block_cargo_bin_dir
      "${RAMAG_CARGO_EXECUTABLE}" DIRECTORY)
    find_program(RAMAG_RUSTC_EXECUTABLE NAMES rustc
      HINTS "${block_cargo_bin_dir}" REQUIRED)
    if(RAMAG_BLOCK_ALIGNER_RUSTUP_HOME)
      set(block_rustup_home "${RAMAG_BLOCK_ALIGNER_RUSTUP_HOME}")
    elseif(DEFINED ENV{RUSTUP_HOME} AND NOT "$ENV{RUSTUP_HOME}" STREQUAL "")
      set(block_rustup_home "$ENV{RUSTUP_HOME}")
    else()
      message(FATAL_ERROR
        "A Block Aligner build requires an isolated RUSTUP_HOME containing "
        "${RAMAG_BLOCK_ALIGNER_RUST_TOOLCHAIN}; set "
        "RAMAG_BLOCK_ALIGNER_RUSTUP_HOME or the RUSTUP_HOME environment variable")
    endif()
    get_filename_component(block_rustup_home "${block_rustup_home}" ABSOLUTE)
    if(RAMAG_BLOCK_ALIGNER_CARGO_HOME)
      get_filename_component(block_cargo_home
        "${RAMAG_BLOCK_ALIGNER_CARGO_HOME}" ABSOLUTE)
    else()
      set(block_cargo_home
        "${CMAKE_CURRENT_BINARY_DIR}/block-aligner-cargo-home")
    endif()
    execute_process(
      COMMAND "${CMAKE_COMMAND}" -E env
        "CARGO_HOME=${block_cargo_home}"
        "RUSTUP_HOME=${block_rustup_home}"
        "RUSTUP_TOOLCHAIN=${RAMAG_BLOCK_ALIGNER_RUST_TOOLCHAIN}"
        "${RAMAG_CARGO_EXECUTABLE}" --version
      RESULT_VARIABLE block_cargo_status
      OUTPUT_VARIABLE block_cargo_version
      ERROR_VARIABLE block_cargo_error
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_STRIP_TRAILING_WHITESPACE)
    if(NOT block_cargo_status EQUAL 0)
      message(FATAL_ERROR
        "The pinned Block Aligner Rust toolchain is unavailable: "
        "${block_cargo_error}")
    endif()
    execute_process(
      COMMAND "${CMAKE_COMMAND}" -E env
        "CARGO_HOME=${block_cargo_home}"
        "RUSTUP_HOME=${block_rustup_home}"
        "RUSTUP_TOOLCHAIN=${RAMAG_BLOCK_ALIGNER_RUST_TOOLCHAIN}"
        "${RAMAG_RUSTC_EXECUTABLE}" --version
      RESULT_VARIABLE block_rustc_status
      OUTPUT_VARIABLE block_rustc_version
      ERROR_VARIABLE block_rustc_error
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_STRIP_TRAILING_WHITESPACE)
    if(NOT block_rustc_status EQUAL 0)
      message(FATAL_ERROR
        "The pinned Block Aligner rustc is unavailable: ${block_rustc_error}")
    endif()
    _ramag_populate_clean_dependency(
      block_aligner https://github.com/Daniel-Liu-c0deb0t/block-aligner.git
      "${RAMAG_BLOCK_ALIGNER_COMMIT}"
      "${RAMAG_BLOCK_ALIGNER_SOURCE_DIR}" "Block Aligner")

    # Cargo writes Cargo.lock beside the selected manifest.  Upstream's
    # dependency-free C ABI crate intentionally has no checked-in lockfile, so
    # build a byte-for-byte source staging copy and freeze a lockfile there.
    # The validated Git source remains clean and is the provenance authority.
    set(block_stage "${CMAKE_CURRENT_BINARY_DIR}/block-aligner-source")
    file(REMOVE_RECURSE "${block_stage}")
    file(COPY "${RAMAG_RESOLVED_block_aligner_SOURCE_DIR}/"
         DESTINATION "${block_stage}"
         PATTERN ".git" EXCLUDE
         PATTERN "target" EXCLUDE)
    set(block_target "${CMAKE_CURRENT_BINARY_DIR}/block-aligner-target")
    set(block_library "${block_target}/release/libblock_aligner_c.a")
    set(block_feature "simd_${RAMAG_BLOCK_ALIGNER_SIMD}")
    add_custom_command(
      OUTPUT "${block_library}"
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${block_cargo_home}"
      COMMAND "${CMAKE_COMMAND}" -E env
        "CARGO_HOME=${block_cargo_home}"
        "RUSTUP_HOME=${block_rustup_home}"
        "CARGO_TARGET_DIR=${block_target}"
        "RUSTUP_TOOLCHAIN=${RAMAG_BLOCK_ALIGNER_RUST_TOOLCHAIN}"
        "${RAMAG_CARGO_EXECUTABLE}" generate-lockfile
        --manifest-path "${block_stage}/c/Cargo.toml" --offline
      COMMAND "${CMAKE_COMMAND}" -E env
        "CARGO_HOME=${block_cargo_home}"
        "RUSTUP_HOME=${block_rustup_home}"
        "CARGO_TARGET_DIR=${block_target}"
        "RUSTUP_TOOLCHAIN=${RAMAG_BLOCK_ALIGNER_RUST_TOOLCHAIN}"
        "${RAMAG_CARGO_EXECUTABLE}" build
        --manifest-path "${block_stage}/c/Cargo.toml"
        --release --locked --offline --no-default-features
        --features "${block_feature}"
      DEPENDS
        "${block_stage}/c/Cargo.toml"
        "${block_stage}/src/lib.rs"
      VERBATIM)
    add_custom_target(ramag_block_aligner_build DEPENDS "${block_library}")
    add_library(ramag_block_aligner STATIC IMPORTED GLOBAL)
    set_target_properties(ramag_block_aligner PROPERTIES
      IMPORTED_LOCATION "${block_library}"
      INTERFACE_INCLUDE_DIRECTORIES "${block_stage}/c")
    add_dependencies(ramag_block_aligner ramag_block_aligner_build)
    set(RAMAG_EXTENSION_TARGET ramag_block_aligner PARENT_SCOPE)
    set(RAMAG_EXTENSION_SOURCE
        "${RAMAG_RESOLVED_block_aligner_SOURCE_DIR}" PARENT_SCOPE)
    set(RAMAG_BLOCK_ALIGNER_CARGO_VERSION
        "${block_cargo_version}" PARENT_SCOPE)
    set(RAMAG_BLOCK_ALIGNER_RUSTC_VERSION
        "${block_rustc_version}" PARENT_SCOPE)
  endif()
endfunction()

function(ramag_link_extension_backend target)
  if(RAMAG_INTERNAL_EXTENSION_BACKEND STREQUAL "scalar")
    target_compile_definitions(${target} PRIVATE
      RAMAG_EXTENSION_SCALAR=1
      RAMAG_EXTENSION_DEPENDENCY_SOURCE="not-linked")
    return()
  endif()

  if(RAMAG_INTERNAL_EXTENSION_BACKEND MATCHES "^ksw2-")
    target_sources(${target} PRIVATE src/extension_ksw2.cpp)
    target_link_libraries(${target} PRIVATE "${RAMAG_EXTENSION_TARGET}")
    if(RAMAG_INTERNAL_EXTENSION_BACKEND STREQUAL "ksw2-exact")
      target_compile_definitions(${target} PRIVATE RAMAG_EXTENSION_KSW2_EXACT=1)
    else()
      target_compile_definitions(${target} PRIVATE RAMAG_EXTENSION_KSW2_BAND_AUTO=1)
    endif()
  elseif(RAMAG_INTERNAL_EXTENSION_BACKEND MATCHES "^block-")
    find_package(Threads REQUIRED)
    target_sources(${target} PRIVATE src/extension_block.cpp)
    target_link_libraries(${target} PRIVATE
      "${RAMAG_EXTENSION_TARGET}" Threads::Threads ${CMAKE_DL_LIBS})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      target_link_libraries(${target} PRIVATE m)
    endif()
    if(RAMAG_INTERNAL_EXTENSION_BACKEND STREQUAL "block-exact")
      target_compile_definitions(${target} PRIVATE RAMAG_EXTENSION_BLOCK_EXACT=1)
    else()
      target_compile_definitions(${target} PRIVATE RAMAG_EXTENSION_BLOCK_ADAPTIVE=1)
    endif()
  endif()
  target_compile_definitions(${target} PRIVATE
    RAMAG_EXTENSION_DEPENDENCY_SOURCE="${RAMAG_EXTENSION_SOURCE}"
    RAMAG_BLOCK_ALIGNER_SIMD="${RAMAG_BLOCK_ALIGNER_SIMD}"
    RAMAG_BLOCK_ALIGNER_RUST_TOOLCHAIN="${RAMAG_BLOCK_ALIGNER_RUST_TOOLCHAIN}"
    RAMAG_BLOCK_ALIGNER_CARGO_VERSION="${RAMAG_BLOCK_ALIGNER_CARGO_VERSION}"
    RAMAG_BLOCK_ALIGNER_RUSTC_VERSION="${RAMAG_BLOCK_ALIGNER_RUSTC_VERSION}")
endfunction()
