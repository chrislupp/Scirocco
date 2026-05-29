cmake_minimum_required(VERSION 3.21)
project(scirocco_viewer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# --- Dependencies ---
find_package(Qt6 REQUIRED COMPONENTS Core Gui Qml Quick QuickControls2 Widgets Charts)
qt_standard_project_setup(REQUIRES 6.5)
qt_policy(SET QTP0004 NEW)

find_package(HDF5 REQUIRED COMPONENTS C CXX)

# --- Sources ---
set(VIEWER_SOURCES
    src/main.cpp
    src/hdf5_reader.h
    src/hdf5_reader.cpp
    src/iteration_data.h
    src/video_exporter.h
    src/video_exporter.cpp
)

# --- Executable (macOS app bundle) ---
qt_add_executable(scirocco_viewer_app MACOSX_BUNDLE ${VIEWER_SOURCES})

set_target_properties(scirocco_viewer_app PROPERTIES
    OUTPUT_NAME "SciroccoViewer"
    MACOSX_BUNDLE TRUE
    MACOSX_BUNDLE_BUNDLE_NAME "Scirocco Viewer"
)

# --- QML module ---
set_source_files_properties(qml/Theme.qml PROPERTIES QT_QML_SINGLETON_TYPE TRUE)

qt_add_qml_module(scirocco_viewer_app
    URI SciroccoViewer
    VERSION 1.0
    QML_FILES
        qml/Main.qml
        qml/Theme.qml
        qml/AirfoilCanvas.qml
        qml/CpChart.qml
        qml/DistributionChart.qml
        qml/ConvergenceChart.qml
        qml/components/Panel.qml
        qml/components/IconRail.qml
        qml/components/StatusBar.qml
        qml/components/InspectorPanel.qml
        qml/components/Chip.qml
        qml/components/SectionHeader.qml
        qml/components/ExportOptions.qml
)

# --- Include & link ---
target_include_directories(scirocco_viewer_app PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(scirocco_viewer_app PRIVATE
    Qt6::Core Qt6::Gui Qt6::Qml Qt6::Quick Qt6::QuickControls2
    Qt6::Widgets Qt6::Charts
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
)
