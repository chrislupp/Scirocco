execute_process(COMMAND "${SERVER}" "invalid-address"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error
    TIMEOUT 15)
# WILL_FAIL alone would also accept some abnormal exits. Require our exact code.
if(NOT "${status}" STREQUAL "1" OR NOT error MATCHES "Failed to start server")
    message(FATAL_ERROR "Expected clean startup failure, got ${status}: ${output}${error}")
endif()
