include(FetchContent)

function(ramag_configure_pinned_dependencies)
  if(RAMAG_SUFKIT_SOURCE_DIR)
    ramag_require_clean_commit(
      "${RAMAG_SUFKIT_SOURCE_DIR}" "${RAMAG_SUFKIT_COMMIT}" "sufkit")
    set(FETCHCONTENT_SOURCE_DIR_SUFKIT "${RAMAG_SUFKIT_SOURCE_DIR}")
  endif()
  if(RAMAG_SEQPRO_SOURCE_DIR)
    ramag_require_clean_commit(
      "${RAMAG_SEQPRO_SOURCE_DIR}" "${RAMAG_SEQPRO_COMMIT}" "SeqPro")
    set(FETCHCONTENT_SOURCE_DIR_SEQPRO "${RAMAG_SEQPRO_SOURCE_DIR}")
  endif()

  set(SUFKIT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(SUFKIT_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
  set(SUFKIT_ENABLE_SEQPRO OFF CACHE BOOL "" FORCE)
  set(SEQPRO_BUILD_TOOLS OFF CACHE BOOL "" FORCE)

  # SeqPro deliberately follows the standard BUILD_TESTING switch instead of
  # defining a project-specific test option.  Keep dependency tests out of the
  # RaMA-G product test registry; each pinned library is validated separately
  # from its clean checkout, then restore the parent value before RaMA-G adds
  # its own tests.
  set(_ramag_parent_build_testing "${BUILD_TESTING}")
  set(BUILD_TESTING OFF CACHE BOOL "Build project tests" FORCE)

  FetchContent_Declare(
    sufkit
    GIT_REPOSITORY https://github.com/malabz/sufkit.git
    GIT_TAG "${RAMAG_SUFKIT_COMMIT}"
    GIT_SHALLOW FALSE
    GIT_PROGRESS TRUE)
  FetchContent_Declare(
    seqpro
    GIT_REPOSITORY https://github.com/malabz/seqpro.git
    GIT_TAG "${RAMAG_SEQPRO_COMMIT}"
    GIT_SHALLOW FALSE
    GIT_PROGRESS TRUE)
  FetchContent_MakeAvailable(sufkit seqpro)
  set(BUILD_TESTING "${_ramag_parent_build_testing}"
      CACHE BOOL "Build project tests" FORCE)

  # Validate the sources FetchContent actually resolved, not only the two
  # RaMA-G-specific source override variables.  This also rejects a dirty
  # reused _deps checkout and a caller-provided FETCHCONTENT_SOURCE_DIR_* that
  # bypasses the normal override path.
  FetchContent_GetProperties(sufkit SOURCE_DIR _ramag_sufkit_source_dir)
  FetchContent_GetProperties(seqpro SOURCE_DIR _ramag_seqpro_source_dir)
  ramag_require_clean_commit(
    "${_ramag_sufkit_source_dir}" "${RAMAG_SUFKIT_COMMIT}" "resolved sufkit")
  ramag_require_clean_commit(
    "${_ramag_seqpro_source_dir}" "${RAMAG_SEQPRO_COMMIT}" "resolved SeqPro")

  if(NOT TARGET sufkit::sufkit)
    message(FATAL_ERROR "Pinned sufkit did not define sufkit::sufkit")
  endif()
  if(NOT TARGET SeqPro::seqpro)
    message(FATAL_ERROR "Pinned SeqPro did not define SeqPro::seqpro")
  endif()
endfunction()
