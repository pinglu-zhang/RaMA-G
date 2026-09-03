function(ramag_resolve_git_revision output directory)
  execute_process(
    COMMAND git -C "${directory}" rev-parse HEAD
    RESULT_VARIABLE git_result
    OUTPUT_VARIABLE git_output
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT git_result EQUAL 0 OR NOT git_output MATCHES "^[0-9a-fA-F]+$")
    set(${output} "unknown" PARENT_SCOPE)
    return()
  endif()

  execute_process(
    COMMAND git -C "${directory}" status --porcelain --untracked-files=all
    RESULT_VARIABLE status_result
    OUTPUT_VARIABLE dirty_state
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT status_result EQUAL 0)
    set(${output} "${git_output}-status-unknown" PARENT_SCOPE)
  elseif(dirty_state STREQUAL "")
    set(${output} "${git_output}" PARENT_SCOPE)
  else()
    set(${output} "${git_output}-dirty" PARENT_SCOPE)
  endif()
endfunction()

function(ramag_require_clean_commit directory expected label)
  if(NOT IS_DIRECTORY "${directory}")
    message(FATAL_ERROR "${label} source override does not exist: ${directory}")
  endif()
  execute_process(
    COMMAND git -C "${directory}" rev-parse HEAD
    RESULT_VARIABLE head_result
    OUTPUT_VARIABLE actual_head
    ERROR_VARIABLE head_error
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT head_result EQUAL 0)
    message(FATAL_ERROR "${label} source override is not a Git checkout: ${head_error}")
  endif()
  if(NOT actual_head STREQUAL expected)
    message(FATAL_ERROR
      "${label} source override commit mismatch: expected ${expected}, got ${actual_head}")
  endif()
  execute_process(
    COMMAND git -C "${directory}" status --porcelain --untracked-files=all
    RESULT_VARIABLE status_result
    OUTPUT_VARIABLE dirty_state
    ERROR_VARIABLE status_error
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT status_result EQUAL 0)
    message(FATAL_ERROR "Cannot inspect ${label} source override: ${status_error}")
  endif()
  if(NOT dirty_state STREQUAL "")
    message(FATAL_ERROR
      "${label} source override must be clean; refusing an uncommitted dependency tree")
  endif()
endfunction()
