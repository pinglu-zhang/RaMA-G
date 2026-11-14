find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter)
add_test(NAME ramag.public.cli
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_cli.py"
    --binary $<TARGET_FILE:ramag> --root "${CMAKE_CURRENT_BINARY_DIR}/cli")
add_test(NAME ramag.public.formats
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_formats.py" $<TARGET_FILE:ramag>)
if(RAMAG_ENABLE_CENTROMERE AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
  add_test(NAME ramag.public.centromere
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_centromere.py" $<TARGET_FILE:ramag>)
else()
  add_test(NAME ramag.public.centromere
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_disabled.py" $<TARGET_FILE:ramag>)
endif()
add_test(NAME ramag.public.install
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_install.py"
    --build "${PROJECT_BINARY_DIR}" --source "${CMAKE_CURRENT_SOURCE_DIR}"
    --root "${CMAKE_CURRENT_BINARY_DIR}/install" --cmake "${CMAKE_COMMAND}"
    --backend "${RAMAG_ENABLE_CENTROMERE}")
set_tests_properties(ramag.public.cli ramag.public.formats ramag.public.centromere
  ramag.public.install PROPERTIES LABELS "public;synthetic" TIMEOUT 120)
set_tests_properties(ramag.public.install PROPERTIES RUN_SERIAL TRUE TIMEOUT 180)
