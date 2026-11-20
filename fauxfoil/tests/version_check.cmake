# Run fauXFOIL with empty input and check its banner: the version, and the
# text-only notice exactly when it was built without its plot window.
# Usage: cmake -DFAUXFOIL=<binary> -DEXPECTED=<version> -DTEXT_ONLY=<ON|OFF>
#              -DWORK=<dir> -P version_check.cmake
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
set(notice "Text only: this fauXFOIL was built without its plot window")
string(FIND "${output}" "${notice}" at)
if(TEXT_ONLY AND at EQUAL -1)
    message(FATAL_ERROR "text-only build does not say so:\n${output}")
elseif(NOT TEXT_ONLY AND NOT at EQUAL -1)
    message(FATAL_ERROR "build with the plot window says it is text only")
endif()
