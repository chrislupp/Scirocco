cmake_minimum_required(VERSION 3.18)
if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(cppfoil VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(CPPFOIL_BUILD_TESTS "Build unit tests" ON)
option(CPPFOIL_BUILD_SENSITIVITY_DEBUG_TESTS
       "Build chain-piece sensitivity tests (for adjoint debugging)" OFF)
option(CPPFOIL_BUILD_PROBES
       "Build standalone probe diagnostic binaries" OFF)
option(CPPFOIL_BUILD_EXAMPLES "Build example programs" ON)
option(CPPFOIL_NATIVE_ARCH "Use -march=native for Release builds" ON)
option(CPPFOIL_COVERAGE "Instrument library + tests for llvm-cov coverage" OFF)

# Find Eigen3 (3.4+ or 5.x)
find_package(Eigen3 REQUIRED)

# Source files
set(CPPFOIL_SOURCES
    src/spline.cpp
    src/geometry.cpp
    src/paneling.cpp
    src/panel_influence.cpp
    src/panel_influence_deriv.cpp
    src/inviscid.cpp
    src/inviscid_deriv.cpp
    src/wake.cpp
    src/boundary_layer.cpp
    src/transition.cpp
    src/coupled_solver.cpp
    src/forces.cpp
    src/adjoint.cpp
    src/viscous_geom_deriv.cpp
    src/viscous_residual_deriv.cpp
    src/bl_properties.cpp
    src/cppfoil.cpp
)

# Main library
add_library(cppfoil ${CPPFOIL_SOURCES})
target_include_directories(cppfoil PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(cppfoil PUBLIC Eigen3::Eigen)
target_compile_options(cppfoil PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-O3>
    $<$<AND:$<CONFIG:Release>,$<BOOL:${CPPFOIL_NATIVE_ARCH}>>:-march=native>
)

if(CPPFOIL_COVERAGE)
    target_compile_options(cppfoil PUBLIC
        -fprofile-instr-generate -fcoverage-mapping)
    target_link_options(cppfoil PUBLIC
        -fprofile-instr-generate -fcoverage-mapping)
endif()

# Tests
if(CPPFOIL_BUILD_TESTS)
    enable_testing()
    find_package(GTest QUIET)
    find_package(hdf5 CONFIG QUIET)
    if(GTest_FOUND)
        include(GoogleTest)
    endif()
    if(GTest_FOUND AND hdf5_FOUND)
        message(STATUS "Found HDF5: ${hdf5_VERSION}")

        # --- Default suite: primal + end-to-end adjoint coverage ---

        add_executable(cppfoil_tests
            tests/spline_test.cpp
            tests/geometry_test.cpp
            tests/mfoil_reference_test.cpp
        )
        target_link_libraries(cppfoil_tests cppfoil GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(cppfoil_tests PRIVATE
            MFOIL_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/mfoil_reference.h5"
        )
        gtest_discover_tests(cppfoil_tests)

        # Forward reference tests: compare BuildGlobSys R_x/R_U vs mfoil
        add_executable(cppfoil_forward_ref_tests
            tests/forward_reference_test.cpp
        )
        target_link_libraries(cppfoil_forward_ref_tests
            cppfoil GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(cppfoil_forward_ref_tests PRIVATE
            FORWARD_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/forward_reference.h5"
        )
        gtest_discover_tests(cppfoil_forward_ref_tests)

        # Sensitivity reference tests: adjoint vs complex-step HDF5 reference data
        add_executable(cppfoil_sensitivity_ref_tests
            tests/sensitivity_reference_test.cpp
        )
        target_link_libraries(cppfoil_sensitivity_ref_tests
            cppfoil GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(cppfoil_sensitivity_ref_tests PRIVATE
            SENSITIVITY_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/sensitivity_reference.h5"
        )
        gtest_discover_tests(cppfoil_sensitivity_ref_tests)

        # Ping tests: cppfoil analytic vs mfoil-CS Jacobian references.
        add_executable(cppfoil_ping_tests tests/ping_test.cpp)
        target_link_libraries(cppfoil_ping_tests
            cppfoil GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(cppfoil_ping_tests PRIVATE
            PING_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/ping_reference.h5")
        gtest_discover_tests(cppfoil_ping_tests)

        # --- Opt-in: chain-piece sensitivity debug tests ---
        # These isolate which link of the adjoint chain is broken when a
        # sensitivity regression appears. The end-to-end sensitivity_ref test
        # above already exercises the same source lines at machine precision,
        # so these are debug aids rather than regression coverage.
        if(CPPFOIL_BUILD_SENSITIVITY_DEBUG_TESTS)
            message(STATUS "Sensitivity-debug tests enabled")

            # Per-panel influence derivatives
            add_executable(cppfoil_sensitivity_tests tests/sensitivity_ping_test.cpp)
            target_link_libraries(cppfoil_sensitivity_tests
                cppfoil GTest::gtest_main hdf5_cpp-shared)
            target_compile_definitions(cppfoil_sensitivity_tests PRIVATE
                SENSITIVITY_PING_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/sensitivity_ping_reference.h5")
            gtest_discover_tests(cppfoil_sensitivity_tests)

            # Phase 1 geometric chain tests
            add_executable(cppfoil_viscous_chain_tests
                tests/viscous_chain_test.cpp
            )
            target_link_libraries(cppfoil_viscous_chain_tests
                cppfoil GTest::gtest_main hdf5_cpp-shared)
            target_compile_definitions(cppfoil_viscous_chain_tests PRIVATE
                GEOM_CHAIN_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/geom_chain_reference.h5"
            )
            gtest_discover_tests(cppfoil_viscous_chain_tests)

            # Phase 2 residual partial tests
            add_executable(cppfoil_viscous_residual_tests
                tests/viscous_residual_deriv_test.cpp
            )
            target_link_libraries(cppfoil_viscous_residual_tests
                cppfoil GTest::gtest_main hdf5_cpp-shared)
            target_compile_definitions(cppfoil_viscous_residual_tests PRIVATE
                COMPREHENSIVE_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/comprehensive_reference.h5"
            )
            gtest_discover_tests(cppfoil_viscous_residual_tests)

            # Phase 3 viscous output partials test (uses mfoil_cs CS reference)
            add_executable(cppfoil_viscous_output_ping_tests
                tests/viscous_output_ping_test.cpp)
            target_link_libraries(cppfoil_viscous_output_ping_tests
                cppfoil GTest::gtest_main hdf5_cpp-shared)
            target_compile_definitions(cppfoil_viscous_output_ping_tests PRIVATE
                VISCOUS_OUTPUT_PING_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/viscous_output_ping_reference.h5")
            gtest_discover_tests(cppfoil_viscous_output_ping_tests)
        endif()

        # --- Opt-in: probe diagnostic binaries (no gtest) ---
        if(CPPFOIL_BUILD_PROBES)
            message(STATUS "Probe diagnostic binaries enabled")

            add_executable(probe_primal_mismatch tests/probe_primal_mismatch.cpp)
            target_link_libraries(probe_primal_mismatch cppfoil hdf5_cpp-shared)

            add_executable(probe_internal_fd tests/probe_internal_fd.cpp)
            target_link_libraries(probe_internal_fd cppfoil)

            add_executable(probe_lambda_cs tests/probe_lambda_cs.cpp)
            target_link_libraries(probe_lambda_cs cppfoil hdf5_cpp-shared)

            add_executable(probe_rv_compare tests/probe_rv_compare.cpp)
            target_link_libraries(probe_rv_compare cppfoil hdf5_cpp-shared)

            add_executable(probe_sstag_ue tests/probe_sstag_ue.cpp)
            target_link_libraries(probe_sstag_ue cppfoil hdf5_cpp-shared)

            add_executable(probe_rv_decompose tests/probe_rv_decompose.cpp)
            target_link_libraries(probe_rv_decompose cppfoil hdf5_cpp-shared)

            add_executable(probe_drv_dre tests/probe_drv_dre.cpp)
            target_link_libraries(probe_drv_dre cppfoil)

            add_executable(probe_drv_dma tests/probe_drv_dma.cpp)
            target_link_libraries(probe_drv_dma cppfoil)

            add_executable(probe_per_closure_dma tests/probe_per_closure_dma.cpp)
            target_link_libraries(probe_per_closure_dma cppfoil)

            add_executable(probe_lin_ma_per_station tests/probe_lin_ma_per_station.cpp)
            target_link_libraries(probe_lin_ma_per_station cppfoil)

            add_executable(probe_wake_closure_dma tests/probe_wake_closure_dma.cpp)
            target_link_libraries(probe_wake_closure_dma cppfoil)

            add_executable(probe_closure_dma_vs_cs tests/probe_closure_dma_vs_cs.cpp)
            target_link_libraries(probe_closure_dma_vs_cs cppfoil hdf5_cpp-shared)

            add_executable(probe_dJ_dMa tests/probe_dJ_dMa.cpp)
            target_link_libraries(probe_dJ_dMa cppfoil)

            add_executable(probe_drv_dma_cs tests/probe_drv_dma_cs.cpp)
            target_link_libraries(probe_drv_dma_cs cppfoil hdf5_cpp-shared)
        endif()
    elseif(GTest_FOUND)
        message(STATUS "HDF5 not found, reference tests will not be built")
        add_executable(cppfoil_tests
            tests/spline_test.cpp
            tests/geometry_test.cpp
        )
        target_link_libraries(cppfoil_tests cppfoil GTest::gtest_main)
        gtest_discover_tests(cppfoil_tests)
    else()
        message(STATUS "Google Test not found, tests will not be built")
    endif()
endif()

# Examples
if(CPPFOIL_BUILD_EXAMPLES)
    add_executable(naca_example examples/naca_example.cpp)
    target_link_libraries(naca_example cppfoil)

    add_executable(viscous_example examples/viscous_example.cpp)
    target_link_libraries(viscous_example cppfoil)
endif()

# Installation
install(TARGETS cppfoil
    EXPORT cppfoilTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)
install(DIRECTORY include/ DESTINATION include)
install(EXPORT cppfoilTargets
    FILE cppfoilTargets.cmake
    NAMESPACE cppfoil::
    DESTINATION lib/cmake/cppfoil
)
