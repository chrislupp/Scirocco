cmake_minimum_required(VERSION 3.18)
project(cppfoil VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(CPPFOIL_BUILD_TESTS "Build unit tests" ON)
option(CPPFOIL_BUILD_EXAMPLES "Build example programs" ON)

# Find Eigen3 (3.4+ or 5.x)
find_package(Eigen3 REQUIRED)

# Source files
set(CPPFOIL_SOURCES
    src/spline.cpp
    src/geometry.cpp
    src/paneling.cpp
    src/panel_influence.cpp
    src/inviscid.cpp
    src/wake.cpp
    src/boundary_layer.cpp
    src/transition.cpp
    src/coupled_solver.cpp
    src/forces.cpp
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
    $<$<CONFIG:Release>:-O3 -march=native>
)

# Tests
if(CPPFOIL_BUILD_TESTS)
    enable_testing()
    find_package(GTest QUIET)
    if(GTest_FOUND)
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
