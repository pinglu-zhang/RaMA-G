# Linux is the production platform of the pinned WFA2 C implementation.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(_rama_default ON)
else()
  set(_rama_default OFF)
endif()
option(RAMAG_ENABLE_CENTROMERE "Build the source-derived RaMA/WFA2 backend" ${_rama_default})
option(RAMAG_CENTROMERE_ANCHORS_ONLY "Temporarily stop align -c after source RaMA anchors" ON)
option(RAMAG_CENTROMERE_DIRECTIONS "Experimental length-driven bidirectional evidence for align -c" OFF)
option(RAMAG_CENTROMERE_BLOCK_ALIGNMENT "Experimental initial-direction blocks followed by recursive RaMA/WFA" OFF)
option(RAMAG_CENTROMERE_ORIENTATION_CONNECTED_GAPS "Experimental continuous-window and known-block connections for v9" OFF)
option(RAMAG_CENTROMERE_ORIENTATION_BRIDGE_UNKNOWN_GAPS "Experimental short unknown gaps with same-direction endpoints and no opposing evidence" OFF)
option(RAMAG_CENTROMERE_ORIENTATION_AUDIT_KNOWN_BLOCKS "Experimental opposing-window audit of old block internal gaps" OFF)
option(RAMAG_CENTROMERE_RECOVER_LONG_GAP_ENDS "Experimental bounded same-direction recovery at long unknown gap ends" OFF)
option(RAMAG_CENTROMERE_RECOVER_SPARSE_PIVOTS "Experimental frozen sparse root pivots inside reliable same-direction gaps" OFF)
option(RAMAG_CENTROMERE_RECOVER_SHORT_GAP_ENDS "Experimental continuous-window partial recovery at short gap ends" OFF)
option(RAMAG_CENTROMERE_RECOVER_MIXED_GAP_ENDS "Experimental independent pure-window recovery at opposite-direction gap ends" OFF)
option(RAMAG_CENTROMERE_RECOVER_SPARSE_GAP_ENDS "Experimental frozen reachable sparse-root recovery at rejected long edge ends" OFF)
option(RAMAG_CENTROMERE_PREFER_STRONGER_SPARSE_GAP_END "Experimental stronger accepted real-root tip when bounded sparse gap ends meet" OFF)
option(RAMAG_CENTROMERE_UNIFY_REAL_ANCHOR_ENDS "Experimental shared original-real-anchor endpoint recovery and same-direction union" OFF)
option(RAMAG_CENTROMERE_REVIEW_LOCAL_COMPETITIVE_GAPS "Experimental local corridor review of competitive same-direction gaps" OFF)
option(RAMAG_CENTROMERE_REVIEW_ISOLATED_OPPOSITE_ANCHORS "Experimental contextual veto review without dropping opposite roots" OFF)
option(RAMAG_CENTROMERE_REVIEW_COLLINEAR_FRAGMENT_CONNECTIONS "Experimental two-flank collinear witness for unknown direction gaps" OFF)
option(RAMAG_CENTROMERE_REVIEW_COLLINEAR_WINDOW_GAPS "Audit short window breaks with existing two-flank collinear witnesses" OFF)
option(RAMAG_CENTROMERE_REVIEW_ALTERNATE_ORIGINAL_ENDPOINTS "Review bounded original endpoint alternatives after the nearest pair fails" OFF)
option(RAMAG_CENTROMERE_REVIEW_NEUTRAL_WINDOW_GAPS "Fully audit uncontested gaps without own-direction window votes" OFF)
option(RAMAG_CENTROMERE_REVIEW_REGION_CHAINS "Experimental complete region-level collinear review" OFF)
option(RAMAG_CENTROMERE_SUPPLEMENT_REGION_ANCHORS "One frozen bounded original-RaMA local search wave" OFF)
option(RAMAG_CENTROMERE_REVIEW_UNIFIED_GAPS "Experimental frozen multi-hypothesis gap review" OFF)
option(RAMAG_CENTROMERE_REVIEW_GAP_INTERIOR_OPPOSITION "Experimental unknown-gap scope for opposite evidence eligibility" OFF)
option(RAMAG_CENTROMERE_REVIEW_KNOWN_BLOCK_COUNTERCHAINS "Experimental final local counterchain review of inherited known blocks" OFF)
option(RAMAG_CENTROMERE_REPAIR_INHERITED_DIRECTIONS "Experimental one-pass inherited direction repair from existing roots" OFF)
option(RAMAG_CENTROMERE_REPAIR_DRIFT_CONSISTENT_GAPS "Experimental unique-root drift-consistent long direction edges" OFF)
set(RAMAG_CENTROMERE_UNIFIED_QUERY_LIMIT 100000 CACHE STRING "Maximum query span in one unified gap search")
set(RAMAG_CENTROMERE_UNIFIED_REFERENCE_LIMIT 512000 CACHE STRING "Maximum reference span in one unified gap search")
set(RAMAG_CENTROMERE_UNIFIED_INPUT_BUDGET 4000000 CACHE STRING "Frozen total two-direction input budget for unified gap search")
foreach(_parameter QUERY_LIMIT REFERENCE_LIMIT INPUT_BUDGET)
  if(NOT RAMAG_CENTROMERE_UNIFIED_${_parameter} MATCHES "^[0-9]+$")
    message(FATAL_ERROR "Unified gap bounds must be nonnegative integers")
  endif()
endforeach()
if(RAMAG_CENTROMERE_REPAIR_DRIFT_CONSISTENT_GAPS AND NOT RAMAG_CENTROMERE_REPAIR_INHERITED_DIRECTIONS)
  message(FATAL_ERROR "Drift-consistent repair requires inherited-direction repair")
endif()
if(RAMAG_CENTROMERE_REPAIR_INHERITED_DIRECTIONS AND NOT RAMAG_CENTROMERE_REVIEW_KNOWN_BLOCK_COUNTERCHAINS)
  message(FATAL_ERROR "Inherited-direction repair requires the frozen known-counterchain baseline")
endif()
if(RAMAG_CENTROMERE_REVIEW_KNOWN_BLOCK_COUNTERCHAINS AND NOT RAMAG_CENTROMERE_REVIEW_GAP_INTERIOR_OPPOSITION)
  message(FATAL_ERROR "Known-block counterchain review requires the frozen gap-interior baseline")
endif()
if(RAMAG_CENTROMERE_REVIEW_GAP_INTERIOR_OPPOSITION AND NOT RAMAG_CENTROMERE_REVIEW_UNIFIED_GAPS)
  message(FATAL_ERROR "Gap-interior opposition review requires unified gap review")
endif()
if(RAMAG_CENTROMERE_REVIEW_UNIFIED_GAPS AND NOT RAMAG_CENTROMERE_SUPPLEMENT_REGION_ANCHORS)
  message(FATAL_ERROR "Unified gap review requires the frozen complete-region baseline")
endif()
if(RAMAG_CENTROMERE_SUPPLEMENT_REGION_ANCHORS AND NOT RAMAG_CENTROMERE_REVIEW_REGION_CHAINS)
  message(FATAL_ERROR "Local region anchor supplementation requires region-chain review")
endif()
if(RAMAG_CENTROMERE_REVIEW_REGION_CHAINS AND NOT RAMAG_CENTROMERE_REVIEW_NEUTRAL_WINDOW_GAPS)
  message(FATAL_ERROR "Region chain review requires the frozen neutral-window layers")
endif()
if(RAMAG_CENTROMERE_REVIEW_NEUTRAL_WINDOW_GAPS AND NOT RAMAG_CENTROMERE_REVIEW_ALTERNATE_ORIGINAL_ENDPOINTS)
  message(FATAL_ERROR "Neutral window review requires bounded alternate original endpoint review")
endif()
if(RAMAG_CENTROMERE_REVIEW_ALTERNATE_ORIGINAL_ENDPOINTS AND NOT RAMAG_CENTROMERE_REVIEW_COLLINEAR_WINDOW_GAPS)
  message(FATAL_ERROR "Alternate endpoint review requires frozen collinear window-gap review")
endif()
if(RAMAG_CENTROMERE_REVIEW_COLLINEAR_WINDOW_GAPS AND NOT RAMAG_CENTROMERE_REVIEW_COLLINEAR_FRAGMENT_CONNECTIONS)
  message(FATAL_ERROR "Collinear window-gap review requires collinear fragment review")
endif()
if(RAMAG_CENTROMERE_REVIEW_COLLINEAR_FRAGMENT_CONNECTIONS AND NOT RAMAG_CENTROMERE_REVIEW_ISOLATED_OPPOSITE_ANCHORS)
  message(FATAL_ERROR "Collinear fragment review requires frozen isolated review layers")
endif()
if(RAMAG_CENTROMERE_REVIEW_ISOLATED_OPPOSITE_ANCHORS AND NOT RAMAG_CENTROMERE_REVIEW_LOCAL_COMPETITIVE_GAPS)
  message(FATAL_ERROR "Isolated opposite review requires frozen local review layers")
endif()
if(RAMAG_CENTROMERE_REVIEW_LOCAL_COMPETITIVE_GAPS AND NOT RAMAG_CENTROMERE_UNIFY_REAL_ANCHOR_ENDS)
  message(FATAL_ERROR "Local competitive gap review requires frozen unified endpoints")
endif()
if(RAMAG_CENTROMERE_UNIFY_REAL_ANCHOR_ENDS AND NOT RAMAG_CENTROMERE_PREFER_STRONGER_SPARSE_GAP_END)
  message(FATAL_ERROR "Unified real endpoints require the frozen v9q layers")
endif()
if(RAMAG_CENTROMERE_PREFER_STRONGER_SPARSE_GAP_END AND NOT RAMAG_CENTROMERE_RECOVER_SPARSE_GAP_ENDS)
  message(FATAL_ERROR "Stronger sparse-gap endpoint preference requires sparse-gap endpoint recovery")
endif()
if(RAMAG_CENTROMERE_RECOVER_SPARSE_GAP_ENDS AND NOT RAMAG_CENTROMERE_RECOVER_MIXED_GAP_ENDS)
  message(FATAL_ERROR "Sparse-gap endpoint recovery requires mixed-gap orientation evidence")
endif()
if(RAMAG_CENTROMERE_RECOVER_MIXED_GAP_ENDS AND NOT RAMAG_CENTROMERE_RECOVER_SHORT_GAP_ENDS)
  message(FATAL_ERROR "Mixed-gap endpoint recovery requires short-gap orientation evidence")
endif()
if(RAMAG_CENTROMERE_RECOVER_SHORT_GAP_ENDS AND NOT RAMAG_CENTROMERE_RECOVER_SPARSE_PIVOTS)
  message(FATAL_ERROR "Short-gap endpoint recovery requires sparse-pivot orientation evidence")
endif()
set(RAMAG_CENTROMERE_BLOCK_PLANNER "initial-direction-blocks-v1" CACHE STRING "Experimental centromere block planner")
set_property(CACHE RAMAG_CENTROMERE_BLOCK_PLANNER PROPERTY STRINGS initial-direction-blocks-v1 joint-local-blocks-v2 local-coherent-frontier-v3 hierarchical-direction-regions-v4 direction-regions-only-v5 fast-local-direction-v6 count-cluster-directions-v8 orientation-evidence-v9)
if(NOT RAMAG_CENTROMERE_BLOCK_PLANNER MATCHES "^(initial-direction-blocks-v1|joint-local-blocks-v2|local-coherent-frontier-v3|hierarchical-direction-regions-v4|direction-regions-only-v5|fast-local-direction-v6|count-cluster-directions-v8|orientation-evidence-v9)$")
  message(FATAL_ERROR "Unknown centromere block planner")
endif()
if(RAMAG_CENTROMERE_ORIENTATION_CONNECTED_GAPS AND NOT RAMAG_CENTROMERE_BLOCK_PLANNER STREQUAL "orientation-evidence-v9")
  message(FATAL_ERROR "Connected orientation gaps require orientation-evidence-v9")
endif()
if(RAMAG_CENTROMERE_ORIENTATION_AUDIT_KNOWN_BLOCKS AND NOT RAMAG_CENTROMERE_BLOCK_PLANNER STREQUAL "orientation-evidence-v9")
  message(FATAL_ERROR "Known-block orientation audit requires orientation-evidence-v9")
endif()
if(RAMAG_CENTROMERE_ORIENTATION_BRIDGE_UNKNOWN_GAPS AND NOT RAMAG_CENTROMERE_ORIENTATION_CONNECTED_GAPS)
  message(FATAL_ERROR "Unknown-gap bridging requires connected orientation gaps")
endif()
if(RAMAG_CENTROMERE_RECOVER_LONG_GAP_ENDS AND
   (NOT RAMAG_CENTROMERE_BLOCK_PLANNER STREQUAL "orientation-evidence-v9" OR
    NOT RAMAG_CENTROMERE_ORIENTATION_CONNECTED_GAPS))
  message(FATAL_ERROR "Long-gap endpoint recovery requires connected orientation-evidence-v9")
endif()
if(RAMAG_CENTROMERE_RECOVER_SPARSE_PIVOTS AND
   (NOT RAMAG_CENTROMERE_BLOCK_PLANNER STREQUAL "orientation-evidence-v9" OR
    NOT RAMAG_CENTROMERE_ORIENTATION_CONNECTED_GAPS OR
    NOT RAMAG_CENTROMERE_ORIENTATION_BRIDGE_UNKNOWN_GAPS OR
    NOT RAMAG_CENTROMERE_ORIENTATION_AUDIT_KNOWN_BLOCKS OR
    NOT RAMAG_CENTROMERE_RECOVER_LONG_GAP_ENDS))
  message(FATAL_ERROR "Sparse root pivots require connected, unknown-gap, audited and long-gap orientation-evidence-v9")
endif()
if(RAMAG_CENTROMERE_BLOCK_PLANNER MATCHES "^(count-cluster-directions-v8|orientation-evidence-v9)$" AND
   (NOT RAMAG_ENABLE_CENTROMERE OR NOT RAMAG_CENTROMERE_ANCHORS_ONLY OR RAMAG_CENTROMERE_BLOCK_ALIGNMENT OR RAMAG_CENTROMERE_DIRECTIONS))
  message(FATAL_ERROR "v8/v9 are direction-only: backend and anchors-only ON, block alignment and legacy directions OFF")
endif()
foreach(_pair WINDOW STEP MIN_ANCHORS DOMINANCE_RATIO CLUSTER_GAP INFERRED_GAP)
  if(_pair STREQUAL "WINDOW" OR _pair STREQUAL "CLUSTER_GAP")
    set(_default 10000)
  elseif(_pair STREQUAL "STEP")
    set(_default 5000)
  elseif(_pair STREQUAL "MIN_ANCHORS")
    set(_default 3)
  elseif(_pair STREQUAL "DOMINANCE_RATIO")
    set(_default 2)
  else()
    set(_default 20000)
  endif()
  set(RAMAG_CENTROMERE_COUNT_${_pair} "${_default}" CACHE STRING "Internal count-cluster ${_pair}")
  if(NOT RAMAG_CENTROMERE_COUNT_${_pair} MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "Count-cluster settings must be positive integers")
  endif()
endforeach()
if(RAMAG_CENTROMERE_BLOCK_PLANNER MATCHES "^(hierarchical-direction-regions-v4|direction-regions-only-v5|fast-local-direction-v6)$")
  set(_rama_wfa_budget "536870912")
else()
  set(_rama_wfa_budget "8589934592")
endif()
set(RAMAG_CENTROMERE_BLOCK_WFA_BYTES "${_rama_wfa_budget}" CACHE STRING "Maximum additional WFA allocation bytes per gap")
set(RAMAG_CENTROMERE_BLOCK_WFA_SECONDS "120" CACHE STRING "Cooperative time limit per WFA gap")
set(RAMAG_CENTROMERE_BLOCK_STAGE_SECONDS "1800" CACHE STRING "Cooperative time limit for block alignment")
set(RAMAG_CENTROMERE_BLOCK_PLANNER_BYTES "10000000000" CACHE STRING "Maximum live planner capacity bytes")
foreach(_limit WFA_BYTES WFA_SECONDS STAGE_SECONDS PLANNER_BYTES)
  if(NOT RAMAG_CENTROMERE_BLOCK_${_limit} MATCHES "^[0-9]+$")
    message(FATAL_ERROR "Block alignment budgets must be unsigned integers")
  endif()
endforeach()
if(RAMAG_CENTROMERE_BLOCK_ALIGNMENT AND (NOT RAMAG_ENABLE_CENTROMERE OR RAMAG_CENTROMERE_ANCHORS_ONLY OR RAMAG_CENTROMERE_DIRECTIONS))
  message(FATAL_ERROR "Block alignment requires the backend, anchors-only OFF and directions OFF")
endif()
set(RAMAG_CENTROMERE_DIRECTION_SEARCH "rama-first-hit-focused-v1" CACHE STRING "Experimental direction search strategy")
set_property(CACHE RAMAG_CENTROMERE_DIRECTION_SEARCH PROPERTY STRINGS rama-first-hit-focused-v1 legacy-multiscale-v1)
if(NOT RAMAG_CENTROMERE_DIRECTION_SEARCH MATCHES "^(rama-first-hit-focused-v1|legacy-multiscale-v1)$")
  message(FATAL_ERROR "Unknown RAMAG_CENTROMERE_DIRECTION_SEARCH")
endif()
if(RAMAG_CENTROMERE_DIRECTIONS AND NOT RAMAG_ENABLE_CENTROMERE)
  message(FATAL_ERROR "Direction evidence requires RAMAG_ENABLE_CENTROMERE")
endif()
if(RAMAG_CENTROMERE_DIRECTIONS AND NOT RAMAG_CENTROMERE_ANCHORS_ONLY)
  message(FATAL_ERROR "Direction evidence is anchors-only; enable RAMAG_CENTROMERE_ANCHORS_ONLY")
endif()
set(RAMAG_WFA2_SOURCE_DIR "" CACHE PATH "Clean pinned WFA2 source override")
add_library(ramag_centromere STATIC src/centromere_alignment.cpp src/centromere.cpp
  src/centromere_blocks.cpp
  src/centromere_joint_planner.cpp
  src/centromere_region_refinement.cpp
  src/centromere_direction.cpp src/centromere_direction_geometry.cpp src/centromere_direction_output.cpp)
add_library(RaMAG::centromere ALIAS ramag_centromere)
target_compile_features(ramag_centromere PUBLIC cxx_std_20)
target_compile_definitions(ramag_centromere PRIVATE
  RAMAG_CENTROMERE_BLOCK_ALIGNMENT=$<BOOL:${RAMAG_CENTROMERE_BLOCK_ALIGNMENT}>
  RAMAG_CENTROMERE_JOINT_BLOCKS=$<STREQUAL:${RAMAG_CENTROMERE_BLOCK_PLANNER},joint-local-blocks-v2>
  RAMAG_CENTROMERE_LOCAL_FRONTIER=$<STREQUAL:${RAMAG_CENTROMERE_BLOCK_PLANNER},local-coherent-frontier-v3>
  RAMAG_CENTROMERE_HIERARCHICAL=$<STREQUAL:${RAMAG_CENTROMERE_BLOCK_PLANNER},hierarchical-direction-regions-v4>
  RAMAG_CENTROMERE_REGIONS_ONLY=$<STREQUAL:${RAMAG_CENTROMERE_BLOCK_PLANNER},direction-regions-only-v5>
  RAMAG_CENTROMERE_FAST_DIRECTION=$<STREQUAL:${RAMAG_CENTROMERE_BLOCK_PLANNER},fast-local-direction-v6>
  RAMAG_CENTROMERE_COUNT_CLUSTER=$<STREQUAL:${RAMAG_CENTROMERE_BLOCK_PLANNER},count-cluster-directions-v8>
  RAMAG_CENTROMERE_ORIENTATION_EVIDENCE=$<STREQUAL:${RAMAG_CENTROMERE_BLOCK_PLANNER},orientation-evidence-v9>
  RAMAG_CENTROMERE_ORIENTATION_CONNECTED_GAPS=$<BOOL:${RAMAG_CENTROMERE_ORIENTATION_CONNECTED_GAPS}>
  RAMAG_CENTROMERE_ORIENTATION_BRIDGE_UNKNOWN_GAPS=$<BOOL:${RAMAG_CENTROMERE_ORIENTATION_BRIDGE_UNKNOWN_GAPS}>
  RAMAG_CENTROMERE_ORIENTATION_AUDIT_KNOWN_BLOCKS=$<BOOL:${RAMAG_CENTROMERE_ORIENTATION_AUDIT_KNOWN_BLOCKS}>
  RAMAG_CENTROMERE_RECOVER_LONG_GAP_ENDS=$<BOOL:${RAMAG_CENTROMERE_RECOVER_LONG_GAP_ENDS}>
  RAMAG_CENTROMERE_RECOVER_SPARSE_PIVOTS=$<BOOL:${RAMAG_CENTROMERE_RECOVER_SPARSE_PIVOTS}>
  RAMAG_CENTROMERE_RECOVER_SHORT_GAP_ENDS=$<BOOL:${RAMAG_CENTROMERE_RECOVER_SHORT_GAP_ENDS}>
  RAMAG_CENTROMERE_RECOVER_MIXED_GAP_ENDS=$<BOOL:${RAMAG_CENTROMERE_RECOVER_MIXED_GAP_ENDS}>
  RAMAG_CENTROMERE_RECOVER_SPARSE_GAP_ENDS=$<BOOL:${RAMAG_CENTROMERE_RECOVER_SPARSE_GAP_ENDS}>
  RAMAG_CENTROMERE_PREFER_STRONGER_SPARSE_GAP_END=$<BOOL:${RAMAG_CENTROMERE_PREFER_STRONGER_SPARSE_GAP_END}>
  RAMAG_CENTROMERE_UNIFY_REAL_ANCHOR_ENDS=$<BOOL:${RAMAG_CENTROMERE_UNIFY_REAL_ANCHOR_ENDS}>
  RAMAG_CENTROMERE_REVIEW_LOCAL_COMPETITIVE_GAPS=$<BOOL:${RAMAG_CENTROMERE_REVIEW_LOCAL_COMPETITIVE_GAPS}>
  RAMAG_CENTROMERE_REVIEW_ISOLATED_OPPOSITE_ANCHORS=$<BOOL:${RAMAG_CENTROMERE_REVIEW_ISOLATED_OPPOSITE_ANCHORS}>
  RAMAG_CENTROMERE_REVIEW_COLLINEAR_FRAGMENT_CONNECTIONS=$<BOOL:${RAMAG_CENTROMERE_REVIEW_COLLINEAR_FRAGMENT_CONNECTIONS}>
  RAMAG_CENTROMERE_REVIEW_COLLINEAR_WINDOW_GAPS=$<BOOL:${RAMAG_CENTROMERE_REVIEW_COLLINEAR_WINDOW_GAPS}>
  RAMAG_CENTROMERE_REVIEW_ALTERNATE_ORIGINAL_ENDPOINTS=$<BOOL:${RAMAG_CENTROMERE_REVIEW_ALTERNATE_ORIGINAL_ENDPOINTS}>
  RAMAG_CENTROMERE_REVIEW_NEUTRAL_WINDOW_GAPS=$<BOOL:${RAMAG_CENTROMERE_REVIEW_NEUTRAL_WINDOW_GAPS}>
  RAMAG_CENTROMERE_REVIEW_REGION_CHAINS=$<BOOL:${RAMAG_CENTROMERE_REVIEW_REGION_CHAINS}>
  RAMAG_CENTROMERE_SUPPLEMENT_REGION_ANCHORS=$<BOOL:${RAMAG_CENTROMERE_SUPPLEMENT_REGION_ANCHORS}>
  RAMAG_CENTROMERE_REVIEW_UNIFIED_GAPS=$<BOOL:${RAMAG_CENTROMERE_REVIEW_UNIFIED_GAPS}>
  RAMAG_CENTROMERE_REVIEW_GAP_INTERIOR_OPPOSITION=$<BOOL:${RAMAG_CENTROMERE_REVIEW_GAP_INTERIOR_OPPOSITION}>
  RAMAG_CENTROMERE_REVIEW_KNOWN_BLOCK_COUNTERCHAINS=$<BOOL:${RAMAG_CENTROMERE_REVIEW_KNOWN_BLOCK_COUNTERCHAINS}>
  RAMAG_CENTROMERE_REPAIR_INHERITED_DIRECTIONS=$<BOOL:${RAMAG_CENTROMERE_REPAIR_INHERITED_DIRECTIONS}>
  RAMAG_CENTROMERE_REPAIR_DRIFT_CONSISTENT_GAPS=$<BOOL:${RAMAG_CENTROMERE_REPAIR_DRIFT_CONSISTENT_GAPS}>
  RAMAG_CENTROMERE_UNIFIED_QUERY_LIMIT=${RAMAG_CENTROMERE_UNIFIED_QUERY_LIMIT}
  RAMAG_CENTROMERE_UNIFIED_REFERENCE_LIMIT=${RAMAG_CENTROMERE_UNIFIED_REFERENCE_LIMIT}
  RAMAG_CENTROMERE_UNIFIED_INPUT_BUDGET=${RAMAG_CENTROMERE_UNIFIED_INPUT_BUDGET}
  RAMAG_CENTROMERE_COUNT_WINDOW=${RAMAG_CENTROMERE_COUNT_WINDOW}ULL
  RAMAG_CENTROMERE_COUNT_STEP=${RAMAG_CENTROMERE_COUNT_STEP}ULL
  RAMAG_CENTROMERE_COUNT_MIN_ANCHORS=${RAMAG_CENTROMERE_COUNT_MIN_ANCHORS}U
  RAMAG_CENTROMERE_COUNT_DOMINANCE_RATIO=${RAMAG_CENTROMERE_COUNT_DOMINANCE_RATIO}U
  RAMAG_CENTROMERE_COUNT_CLUSTER_GAP=${RAMAG_CENTROMERE_COUNT_CLUSTER_GAP}ULL
  RAMAG_CENTROMERE_COUNT_INFERRED_GAP=${RAMAG_CENTROMERE_COUNT_INFERRED_GAP}ULL
  RAMAG_CENTROMERE_BLOCK_PLANNER_BYTES=${RAMAG_CENTROMERE_BLOCK_PLANNER_BYTES}ULL
  RAMAG_CENTROMERE_BLOCK_WFA_BYTES=${RAMAG_CENTROMERE_BLOCK_WFA_BYTES}ULL
  RAMAG_CENTROMERE_BLOCK_WFA_SECONDS=${RAMAG_CENTROMERE_BLOCK_WFA_SECONDS}
  RAMAG_CENTROMERE_BLOCK_STAGE_SECONDS=${RAMAG_CENTROMERE_BLOCK_STAGE_SECONDS}
  RAMAG_CENTROMERE_ANCHORS_ONLY=$<BOOL:${RAMAG_CENTROMERE_ANCHORS_ONLY}>
  RAMAG_DIRECTION_LEGACY_SEARCH=$<STREQUAL:${RAMAG_CENTROMERE_DIRECTION_SEARCH},legacy-multiscale-v1>
  RAMAG_CENTROMERE_DIRECTIONS=$<BOOL:${RAMAG_CENTROMERE_DIRECTIONS}>)
target_include_directories(ramag_centromere PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(ramag_centromere PROPERTIES EXPORT_NAME centromere POSITION_INDEPENDENT_CODE ON)
ramag_apply_project_warnings(ramag_centromere)
target_link_libraries(ramag_core PRIVATE ramag_centromere)
if(RAMAG_ENABLE_CENTROMERE)
  if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(FATAL_ERROR "Pinned centromere backend currently requires Linux")
  endif()
  include(FetchContent)
  set(_wfa_commit 99a5879a3a2972eff168cbcc6610b6cf02844ed0)
  if(RAMAG_WFA2_SOURCE_DIR)
    ramag_require_clean_commit("${RAMAG_WFA2_SOURCE_DIR}" "${_wfa_commit}" "WFA2")
    set(_wfa "${RAMAG_WFA2_SOURCE_DIR}")
  else()
    FetchContent_Declare(ramag_wfa2 GIT_REPOSITORY https://github.com/smarco/WFA2-lib.git
      GIT_TAG "${_wfa_commit}" GIT_SHALLOW FALSE)
    FetchContent_GetProperties(ramag_wfa2)
    if(NOT ramag_wfa2_POPULATED)
      FetchContent_Populate(ramag_wfa2)
    endif()
    set(_wfa "${ramag_wfa2_SOURCE_DIR}")
    ramag_require_clean_commit("${_wfa}" "${_wfa_commit}" "WFA2")
  endif()
  add_library(ramag_rama_dependencies OBJECT
    third_party/rama/core/gsacak_guarded.cpp
    third_party/rama/core/rare_match.cpp
    third_party/rama/core/RMQ.cpp
    third_party/rama/core/anchor.cpp
  )
  target_compile_features(ramag_rama_dependencies PRIVATE cxx_std_20)
  target_include_directories(ramag_rama_dependencies PRIVATE third_party/rama/core)
  set_target_properties(ramag_rama_dependencies PROPERTIES POSITION_INDEPENDENT_CODE ON)
  add_library(ramag_wfa_objects OBJECT
  "${_wfa}/wavefront/wavefront_align.c"
  "${_wfa}/wavefront/wavefront_aligner.c"
  "${_wfa}/wavefront/wavefront_attributes.c"
  "${_wfa}/wavefront/wavefront_backtrace_buffer.c"
  third_party/rama/wfa_safety/wavefront_backtrace.c
  "${_wfa}/wavefront/wavefront_backtrace_offload.c"
  "${_wfa}/wavefront/wavefront_bialign.c"
  "${_wfa}/wavefront/wavefront_bialigner.c"
  "${_wfa}/wavefront/wavefront.c"
  "${_wfa}/wavefront/wavefront_components.c"
  "${_wfa}/wavefront/wavefront_compute_affine2p.c"
  "${_wfa}/wavefront/wavefront_compute_affine.c"
  "${_wfa}/wavefront/wavefront_compute.c"
  "${_wfa}/wavefront/wavefront_compute_edit.c"
  "${_wfa}/wavefront/wavefront_compute_linear.c"
  "${_wfa}/wavefront/wavefront_debug.c"
  "${_wfa}/wavefront/wavefront_display.c"
  "${_wfa}/wavefront/wavefront_extend.c"
  "${_wfa}/wavefront/wavefront_heuristic.c"
  third_party/rama/wfa_safety/wavefront_pcigar.c
  "${_wfa}/wavefront/wavefront_penalties.c"
  "${_wfa}/wavefront/wavefront_plot.c"
  "${_wfa}/wavefront/wavefront_sequences.c"
  "${_wfa}/wavefront/wavefront_slab.c"
  "${_wfa}/wavefront/wavefront_unialign.c"
  "${_wfa}/wavefront/wavefront_termination.c"
  third_party/rama/wfa_safety/wavefront_extend_kernels_avx.c
  third_party/rama/wfa_safety/wavefront_extend_kernels.c
  "${_wfa}/system/mm_stack.c"
  "${_wfa}/system/mm_allocator.c"
  "${_wfa}/system/profiler_counter.c"
  "${_wfa}/system/profiler_timer.c"
  "${_wfa}/utils/bitmap.c"
  "${_wfa}/utils/dna_text.c"
  "${_wfa}/utils/sequence_buffer.c"
  "${_wfa}/utils/vector.c"
  "${_wfa}/utils/commons.c"
  "${_wfa}/utils/heatmap.c"
  "${_wfa}/alignment/affine2p_penalties.c"
  "${_wfa}/alignment/affine_penalties.c"
  "${_wfa}/alignment/cigar.c"
  "${_wfa}/alignment/score_matrix.c"
  )
  target_include_directories(ramag_wfa_objects SYSTEM PRIVATE "${_wfa}" "${_wfa}/wavefront")
  target_compile_options(ramag_wfa_objects PRIVATE -fexceptions
    -include "${PROJECT_SOURCE_DIR}/src/rama_wfa_hooks.h")
  set_target_properties(ramag_wfa_objects PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 11)
  target_sources(ramag_centromere PRIVATE src/rama_wfa_memory.cpp
    $<TARGET_OBJECTS:ramag_rama_dependencies> $<TARGET_OBJECTS:ramag_wfa_objects>)
  target_include_directories(ramag_centromere SYSTEM PRIVATE third_party/rama/core "${_wfa}")
  target_compile_definitions(ramag_centromere PRIVATE RAMAG_HAVE_CENTROMERE=1)
  target_link_libraries(ramag_centromere PRIVATE Threads::Threads m)
  if(RAMAG_ENABLE_SANITIZERS)
    foreach(_t ramag_rama_dependencies ramag_wfa_objects)
      target_compile_options(${_t} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    endforeach()
  endif()
else()
  target_compile_definitions(ramag_centromere PRIVATE RAMAG_HAVE_CENTROMERE=0)
endif()
if(RAMAG_ENABLE_SANITIZERS AND NOT MSVC)
  target_compile_options(ramag_centromere PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(ramag_centromere PUBLIC -fsanitize=address,undefined)
endif()
configure_package_config_file(cmake/RaMAGCentromereConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/RaMAGCentromereConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RaMAGCentromere)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/RaMAGCentromereConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(TARGETS ramag_centromere EXPORT RaMAGCentromereTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(EXPORT RaMAGCentromereTargets NAMESPACE RaMAG:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RaMAGCentromere)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/RaMAGCentromereConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/RaMAGCentromereConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RaMAGCentromere)
install(FILES include/ramag/centromere.hpp include/ramag/centromere_direction.hpp include/ramag/centromere_blocks.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ramag)
install(FILES third_party/rama/upstream/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/ramag/third_party/rama)
if(RAMAG_ENABLE_CENTROMERE)
  install(FILES "${_wfa}/LICENSE" DESTINATION ${CMAKE_INSTALL_DATADIR}/ramag/third_party/wfa2)
endif()
