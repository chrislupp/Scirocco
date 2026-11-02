function(run)
    execute_process(COMMAND ${ARGV} RESULT_VARIABLE status)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "Consumer command failed: ${ARGV}")
    endif()
endfunction()
set(root "${BUILD_DIR}/install-smoke")
file(REMOVE_RECURSE "${root}")
run("${CMAKE_COMMAND}" --install "${BUILD_DIR}" --prefix "${root}/original" --config "${CONFIG}")
# Move the prefix to verify the package does not depend on its install location.
file(RENAME "${root}/original" "${root}/relocated")
run("${CMAKE_COMMAND}" -S "${SOURCE_DIR}/tests/install_consumer" -B "${root}/consumer"
    "-DCMAKE_PREFIX_PATH=${root}/relocated" "-DEigen3_DIR=${EIGEN_DIR}"
    "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}" "-DEXPECTED_VERSION=${EXPECTED_VERSION}")
run("${CMAKE_COMMAND}" --build "${root}/consumer" --config "${CONFIG}")
run("${CMAKE_CTEST_COMMAND}" --test-dir "${root}/consumer" -C "${CONFIG}" --output-on-failure)
