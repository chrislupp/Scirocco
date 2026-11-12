cmake_minimum_required(VERSION 3.18)
if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(scirocco VERSION 0.1.0 LANGUAGES CXX)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(SCIROCCO_BUILD_TESTS "Build unit tests" ON)
option(SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS
       "Build chain-piece sensitivity tests (for adjoint debugging)" OFF)
option(SCIROCCO_BUILD_PROBES
       "Build standalone probe diagnostic binaries" OFF)
option(SCIROCCO_BUILD_EXAMPLES "Build example programs" ON)
option(SCIROCCO_HDF5_OUTPUT "Build HDF5 output writer library" OFF)
option(SCIROCCO_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(SCIROCCO_NATIVE_ARCH "Use -march=native for Release builds" ON)
# Find Eigen3 (3.4+ or 5.x)
find_package(Eigen3 REQUIRED)

# Source files
set(SCIROCCO_SOURCES
    src/spline.cpp
    src/spline_deriv.cpp
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
    src/scirocco.cpp
)

# Main library
add_library(scirocco ${SCIROCCO_SOURCES})
target_include_directories(scirocco PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(scirocco PUBLIC Eigen3::Eigen)
target_compile_features(scirocco PUBLIC cxx_std_17)
target_compile_options(scirocco PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-O3>
    $<$<AND:$<CONFIG:Release>,$<BOOL:${SCIROCCO_NATIVE_ARCH}>>:-march=native>
)

# Propagate instrumentation and runtime linkage to tests and installed consumers.
if(SCIROCCO_SANITIZERS)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(FATAL_ERROR "SCIROCCO_SANITIZERS requires GCC or Clang")
    endif()
    target_compile_options(scirocco PUBLIC
        -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(scirocco PUBLIC -fsanitize=address,undefined)
endif()

# HDF5 output writer (optional separate library)
if(SCIROCCO_HDF5_OUTPUT)
    find_package(hdf5 CONFIG QUIET)
    if(NOT hdf5_FOUND)
        find_package(HDF5 COMPONENTS CXX REQUIRED)
        if(NOT TARGET hdf5_cpp-shared)
            add_library(hdf5_cpp-shared INTERFACE IMPORTED)
            set_target_properties(hdf5_cpp-shared PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${HDF5_CXX_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${HDF5_CXX_LIBRARIES}"
            )
        endif()
    endif()
    message(STATUS "Building scirocco_hdf5 writer")

    add_library(scirocco_hdf5 src/hdf5_writer.cpp)
    target_include_directories(scirocco_hdf5 PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
    target_link_libraries(scirocco_hdf5 PUBLIC scirocco hdf5_cpp-shared)
    target_compile_options(scirocco_hdf5 PRIVATE
        -Wall -Wextra -Wpedantic
        $<$<CONFIG:Release>:-O3>
        $<$<AND:$<CONFIG:Release>,$<BOOL:${SCIROCCO_NATIVE_ARCH}>>:-march=native>
    )

    if(SCIROCCO_BUILD_EXAMPLES)
        add_executable(hdf5_output_example examples/hdf5_output_example.cpp)
        target_link_libraries(hdf5_output_example scirocco_hdf5)
    endif()
endif()

# Tests
if(SCIROCCO_BUILD_TESTS)
    enable_testing()
    add_test(NAME InstalledConsumer
        COMMAND ${CMAKE_COMMAND} -DBUILD_DIR=${CMAKE_BINARY_DIR}
            -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -DCONFIG=$<CONFIG>
            -DEIGEN_DIR=${Eigen3_DIR} -DCXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DEXPECTED_VERSION=${PROJECT_VERSION}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/install_consumer.cmake)
    find_package(GTest QUIET)
    find_package(hdf5 CONFIG QUIET)
    if(NOT hdf5_FOUND)
        find_package(HDF5 COMPONENTS CXX QUIET)
        if(HDF5_FOUND)
            if(NOT TARGET hdf5_cpp-shared)
                add_library(hdf5_cpp-shared INTERFACE IMPORTED)
                set_target_properties(hdf5_cpp-shared PROPERTIES
                    INTERFACE_INCLUDE_DIRECTORIES "${HDF5_CXX_INCLUDE_DIRS}"
                    INTERFACE_LINK_LIBRARIES "${HDF5_CXX_LIBRARIES}"
                )
            endif()
            set(hdf5_FOUND TRUE)
        endif()
    endif()
    if(GTest_FOUND)
        include(GoogleTest)
        add_executable(scirocco_solver_strategy_tests tests/solver_strategy_test.cpp
            tests/input_validation_test.cpp)
        target_link_libraries(scirocco_solver_strategy_tests scirocco GTest::gtest_main)
        gtest_discover_tests(scirocco_solver_strategy_tests)
        # Timing guard: only meaningful in optimized, uninstrumented builds.
        if(CMAKE_BUILD_TYPE MATCHES "^(Release|RelWithDebInfo)$"
           AND NOT SCIROCCO_SANITIZERS)
            add_executable(scirocco_gradient_cost_tests tests/gradient_cost_test.cpp)
            target_link_libraries(scirocco_gradient_cost_tests scirocco GTest::gtest_main)
            gtest_discover_tests(scirocco_gradient_cost_tests)
        endif()
    endif()
    if(GTest_FOUND AND hdf5_FOUND)
        message(STATUS "Found HDF5")

        # --- Default suite: primal + end-to-end adjoint coverage ---

        add_executable(scirocco_tests
            tests/vortex_stream_test.cpp
            tests/spline_test.cpp
            tests/geometry_test.cpp
            tests/mfoil_reference_test.cpp
        )
        target_link_libraries(scirocco_tests scirocco GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(scirocco_tests PRIVATE
            SOURCE_STREAM_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/source_stream_reference.txt"
            VORTEX_STREAM_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/vortex_stream_reference.txt"
            MFOIL_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/mfoil_reference.h5"
        )
        gtest_discover_tests(scirocco_tests)

        # Forward reference tests: compare BuildGlobSys R_x/R_U vs mfoil
        add_executable(scirocco_forward_ref_tests
            tests/forward_reference_test.cpp
        )
        target_link_libraries(scirocco_forward_ref_tests
            scirocco GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(scirocco_forward_ref_tests PRIVATE
            FORWARD_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/forward_reference.h5"
        )
        gtest_discover_tests(scirocco_forward_ref_tests)

        # Sensitivity reference tests: adjoint vs complex-step HDF5 reference data
        add_executable(scirocco_sensitivity_ref_tests
            tests/sensitivity_reference_test.cpp
        )
        target_link_libraries(scirocco_sensitivity_ref_tests
            scirocco GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(scirocco_sensitivity_ref_tests PRIVATE
            SENSITIVITY_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/sensitivity_reference.h5"
        )
        gtest_discover_tests(scirocco_sensitivity_ref_tests)

        add_executable(scirocco_gradient_contract_tests tests/gradient_contract_test.cpp)
        target_link_libraries(scirocco_gradient_contract_tests
            scirocco GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(scirocco_gradient_contract_tests PRIVATE
            GRADIENT_CONTRACT_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/gradient_contract_reference.h5")
        gtest_discover_tests(scirocco_gradient_contract_tests)

        add_executable(scirocco_remaining_gradient_tests tests/remaining_gradient_test.cpp)
        target_link_libraries(scirocco_remaining_gradient_tests
            scirocco GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(scirocco_remaining_gradient_tests PRIVATE
            REMAINING_GRADIENT_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/remaining_gradient_reference.h5"
            HIGH_PRECISION_GRADIENT_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/high_precision_gradient_reference.txt")
        gtest_discover_tests(scirocco_remaining_gradient_tests)

        # Ping tests: scirocco analytic vs mfoil-CS Jacobian references.
        add_executable(scirocco_ping_tests tests/ping_test.cpp)
        target_link_libraries(scirocco_ping_tests
            scirocco GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(scirocco_ping_tests PRIVATE
            PING_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/ping_reference.h5")
        gtest_discover_tests(scirocco_ping_tests)

        # Geometric chain tests: closed-form forward chains
        # (DFoilSDFoilXCumulative, DSstagDFoilX, DWakeSDFoilX, DXiDFoilX,
        # DWgapDFoilX, DWakeXDFoilX) validated against complex-step references
        # in geom_chain_reference.h5. Promoted into the default suite once all
        # chains validated at 1e-10 rel / 2e-11 abs.
        add_executable(scirocco_viscous_chain_tests
            tests/viscous_chain_test.cpp
        )
        target_link_libraries(scirocco_viscous_chain_tests
            scirocco GTest::gtest_main hdf5_cpp-shared)
        target_compile_definitions(scirocco_viscous_chain_tests PRIVATE
            GEOM_CHAIN_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/geom_chain_reference.h5"
        )
        gtest_discover_tests(scirocco_viscous_chain_tests)

        # --- Opt-in: chain-piece sensitivity debug tests ---
        # These isolate which link of the adjoint chain is broken when a
        # sensitivity regression appears. The end-to-end sensitivity_ref test
        # above already exercises the same source lines at machine precision,
        # so these are debug aids rather than regression coverage.
        if(SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS)
            message(STATUS "Sensitivity-debug tests enabled")

            # Per-panel influence derivatives
            add_executable(scirocco_sensitivity_tests tests/sensitivity_ping_test.cpp)
            target_link_libraries(scirocco_sensitivity_tests
                scirocco GTest::gtest_main hdf5_cpp-shared)
            target_compile_definitions(scirocco_sensitivity_tests PRIVATE
                SENSITIVITY_PING_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/sensitivity_ping_reference.h5")
            gtest_discover_tests(scirocco_sensitivity_tests)

            # Phase 2 residual partial tests
            add_executable(scirocco_viscous_residual_tests
                tests/viscous_residual_deriv_test.cpp
            )
            target_link_libraries(scirocco_viscous_residual_tests
                scirocco GTest::gtest_main hdf5_cpp-shared)
            target_compile_definitions(scirocco_viscous_residual_tests PRIVATE
                COMPREHENSIVE_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/comprehensive_reference.h5"
            )
            gtest_discover_tests(scirocco_viscous_residual_tests)

            # Phase 3 viscous output partials test (uses mfoil_cs CS reference)
            add_executable(scirocco_viscous_output_ping_tests
                tests/viscous_output_ping_test.cpp)
            target_link_libraries(scirocco_viscous_output_ping_tests
                scirocco GTest::gtest_main hdf5_cpp-shared)
            target_compile_definitions(scirocco_viscous_output_ping_tests PRIVATE
                VISCOUS_OUTPUT_PING_REFERENCE_FILE="${CMAKE_CURRENT_SOURCE_DIR}/tests/viscous_output_ping_reference.h5")
            gtest_discover_tests(scirocco_viscous_output_ping_tests)
        endif()

        # --- Opt-in: probe diagnostic binaries (no gtest) ---
        if(SCIROCCO_BUILD_PROBES)
            message(STATUS "Probe diagnostic binaries enabled")

            add_executable(probe_primal_mismatch tests/probe_primal_mismatch.cpp)
            target_link_libraries(probe_primal_mismatch scirocco hdf5_cpp-shared)

            add_executable(probe_internal_fd tests/probe_internal_fd.cpp)
            target_link_libraries(probe_internal_fd scirocco)

            add_executable(probe_lambda_cs tests/probe_lambda_cs.cpp)
            target_link_libraries(probe_lambda_cs scirocco hdf5_cpp-shared)

            add_executable(probe_rv_compare tests/probe_rv_compare.cpp)
            target_link_libraries(probe_rv_compare scirocco hdf5_cpp-shared)

            add_executable(probe_sstag_ue tests/probe_sstag_ue.cpp)
            target_link_libraries(probe_sstag_ue scirocco hdf5_cpp-shared)

            add_executable(probe_rv_decompose tests/probe_rv_decompose.cpp)
            target_link_libraries(probe_rv_decompose scirocco hdf5_cpp-shared)

            add_executable(probe_drv_dre tests/probe_drv_dre.cpp)
            target_link_libraries(probe_drv_dre scirocco)

            add_executable(probe_drv_dma tests/probe_drv_dma.cpp)
            target_link_libraries(probe_drv_dma scirocco)

            add_executable(probe_per_closure_dma tests/probe_per_closure_dma.cpp)
            target_link_libraries(probe_per_closure_dma scirocco)

            add_executable(probe_lin_ma_per_station tests/probe_lin_ma_per_station.cpp)
            target_link_libraries(probe_lin_ma_per_station scirocco)

            add_executable(probe_wake_closure_dma tests/probe_wake_closure_dma.cpp)
            target_link_libraries(probe_wake_closure_dma scirocco)

            add_executable(probe_closure_dma_vs_cs tests/probe_closure_dma_vs_cs.cpp)
            target_link_libraries(probe_closure_dma_vs_cs scirocco hdf5_cpp-shared)

            add_executable(probe_dJ_dMa tests/probe_dJ_dMa.cpp)
            target_link_libraries(probe_dJ_dMa scirocco)

            add_executable(probe_drv_dma_cs tests/probe_drv_dma_cs.cpp)
            target_link_libraries(probe_drv_dma_cs scirocco hdf5_cpp-shared)

            set(SCIROCCO_PROBE_TARGETS
                probe_primal_mismatch
                probe_lambda_cs
                probe_rv_compare
                probe_sstag_ue
                probe_rv_decompose
                probe_closure_dma_vs_cs
                probe_drv_dma_cs
            )
            foreach(probe_target IN LISTS SCIROCCO_PROBE_TARGETS)
                target_compile_definitions(${probe_target} PRIVATE
                    SCIROCCO_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests")
            endforeach()
        endif()
    elseif(GTest_FOUND)
        message(STATUS "HDF5 not found, reference tests will not be built")
        add_executable(scirocco_tests
            tests/spline_test.cpp
            tests/geometry_test.cpp
        )
        target_link_libraries(scirocco_tests scirocco GTest::gtest_main)
        gtest_discover_tests(scirocco_tests)
    else()
        message(STATUS "Google Test not found, tests will not be built")
    endif()
endif()

# Examples
if(SCIROCCO_BUILD_EXAMPLES)
    add_executable(naca_example examples/naca_example.cpp)
    target_link_libraries(naca_example scirocco)

    add_executable(viscous_example examples/viscous_example.cpp)
    target_link_libraries(viscous_example scirocco)

    add_executable(quickstart examples/quickstart.cpp)
    target_link_libraries(quickstart scirocco)

    add_executable(bench examples/bench.cpp)
    target_link_libraries(bench scirocco)

    add_executable(docs_figure_data examples/docs_figure_data.cpp)
    target_link_libraries(docs_figure_data scirocco)
endif()

# Philote-Cpp discipline interface
option(SCIROCCO_BUILD_PHILOTE "Build Philote-Cpp discipline interface" OFF)

if(SCIROCCO_BUILD_PHILOTE)
    find_package(Protobuf CONFIG QUIET)
    if(NOT Protobuf_FOUND)
        find_package(Protobuf REQUIRED)
        if(NOT TARGET protobuf::libprotobuf AND TARGET protobuf::protobuf)
            add_library(protobuf::libprotobuf ALIAS protobuf::protobuf)
        elseif(NOT TARGET protobuf::libprotobuf)
            add_library(protobuf::libprotobuf INTERFACE IMPORTED)
            set_target_properties(protobuf::libprotobuf PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${Protobuf_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${Protobuf_LIBRARIES}"
            )
        endif()
    endif()
    find_package(gRPC CONFIG QUIET)
    if(NOT gRPC_FOUND)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(GRPC REQUIRED grpc++)
        if(NOT TARGET gRPC::grpc++)
            add_library(gRPC::grpc++ INTERFACE IMPORTED)
            set_target_properties(gRPC::grpc++ PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${GRPC_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${GRPC_LINK_LIBRARIES}"
            )
        endif()
    endif()
    find_package(PhiloteCpp REQUIRED)
    add_subdirectory(src/philote)
endif()

# Installation
include(CMakePackageConfigHelpers)
install(TARGETS scirocco
    EXPORT sciroccoTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT sciroccoTargets
    FILE sciroccoTargets.cmake
    NAMESPACE scirocco::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scirocco
)

configure_package_config_file(cmake/sciroccoConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/sciroccoConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scirocco)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/sciroccoConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/sciroccoConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/sciroccoConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scirocco)
