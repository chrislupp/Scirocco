cmake_minimum_required(VERSION 3.18)
project(xfoil_cpp VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(XFOIL_ENABLE_ENZYME "Enable Enzyme AD support" OFF)
option(XFOIL_BUILD_TESTS "Build unit tests" ON)

# Dependencies
find_package(Eigen3 REQUIRED)

# Main library
add_library(xfoil_cpp
    src/spline.cpp
    src/geometry.cpp
)
target_include_directories(xfoil_cpp PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(xfoil_cpp PUBLIC Eigen3::Eigen)
target_compile_options(xfoil_cpp PRIVATE -Wall -Wextra -Wpedantic)

# Enzyme support (future)
if(XFOIL_ENABLE_ENZYME)
    find_package(Enzyme REQUIRED)
    target_compile_definitions(xfoil_cpp PUBLIC XFOIL_HAS_ENZYME)
endif()

# Tests
if(XFOIL_BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)

    add_executable(xfoil_tests
        tests/spline_test.cpp
        tests/geometry_test.cpp
    )
    target_link_libraries(xfoil_tests xfoil_cpp GTest::gtest_main)
    gtest_discover_tests(xfoil_tests)
endif()
