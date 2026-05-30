add_library(SiroccoPhilote
    philote_discipline.cpp
    geometry_discipline.cpp
)

add_library(Scirocco::Philote ALIAS SiroccoPhilote)

target_include_directories(SiroccoPhilote
    BEFORE PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<TARGET_PROPERTY:PhiloteCpp::PhiloteCpp,INTERFACE_INCLUDE_DIRECTORIES>
)

target_link_libraries(SiroccoPhilote
    PUBLIC
        scirocco
        PhiloteCpp::PhiloteCpp
)

target_compile_features(SiroccoPhilote PUBLIC cxx_std_17)

if(SCIROCCO_HDF5_OUTPUT AND TARGET scirocco_hdf5)
    target_link_libraries(SiroccoPhilote PUBLIC scirocco_hdf5)
    target_compile_definitions(SiroccoPhilote PUBLIC SCIROCCO_HDF5)
endif()

# Server executables
add_executable(scirocco_server server.cpp)
target_link_libraries(scirocco_server PRIVATE SiroccoPhilote)

add_executable(scirocco_geometry_server geometry_server.cpp)
target_link_libraries(scirocco_geometry_server PRIVATE SiroccoPhilote)

# Tests
if(SCIROCCO_BUILD_TESTS AND GTest_FOUND AND hdf5_FOUND)
    add_executable(scirocco_philote_tests
        ${PROJECT_SOURCE_DIR}/tests/philote_discipline_test.cpp
    )
    target_link_libraries(scirocco_philote_tests
        PRIVATE SiroccoPhilote GTest::gtest_main hdf5_cpp-shared
    )
    target_compile_definitions(scirocco_philote_tests PRIVATE
        SENSITIVITY_REFERENCE_FILE="${PROJECT_SOURCE_DIR}/tests/sensitivity_reference.h5"
    )
    gtest_discover_tests(scirocco_philote_tests)
endif()

# Geometry discipline tests (no HDF5 dependency)
if(SCIROCCO_BUILD_TESTS AND GTest_FOUND)
    add_executable(scirocco_geometry_discipline_tests
        ${PROJECT_SOURCE_DIR}/tests/geometry_discipline_test.cpp
    )
    target_link_libraries(scirocco_geometry_discipline_tests
        PRIVATE SiroccoPhilote GTest::gtest_main
    )
    gtest_discover_tests(scirocco_geometry_discipline_tests)
endif()
