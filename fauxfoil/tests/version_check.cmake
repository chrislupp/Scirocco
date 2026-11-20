# Run fauXFOIL with empty input and check the version in its banner.
# Usage: cmake -DFAUXFOIL=<binary> -DEXPECTED=<version> -DWORK=<dir> -P version_check.cmake
file(WRITE "${WORK}/empty_input.txt" "")
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env FAUXFOIL_NOGRAPHICS=1 "${FAUXFOIL}"
    INPUT_FILE "${WORK}/empty_input.txt"
    OUTPUT_VARIABLE output
    RESULT_VARIABLE status
    TIMEOUT 30)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "fauXFOIL exited with ${status}")
endif()
if(NOT output MATCHES "fauXFOIL Version ${EXPECTED}\n")
    message(FATAL_ERROR "banner does not show version ${EXPECTED}:\n${output}")
endif()
