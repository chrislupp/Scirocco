cmake_minimum_required(VERSION 3.21)
project(scirocco_viewer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# --- Dependencies ---
# 6.6 for the Qt Quick Shapes curve renderer the icons use.
find_package(Qt6 6.6 REQUIRED COMPONENTS Core Gui Qml Quick QuickControls2 Svg Widgets Charts)
qt_standard_project_setup(REQUIRES 6.6)
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

# App icon: resources/app-icon.svg drawn by resources/make_app_icon.py,
# which also writes the macOS icon set.
set(VIEWER_ICNS ${CMAKE_CURRENT_SOURCE_DIR}/resources/AppIcon.icns)
set_source_files_properties(${VIEWER_ICNS} PROPERTIES
    MACOSX_PACKAGE_LOCATION "Resources")
target_sources(scirocco_viewer_app PRIVATE ${VIEWER_ICNS})

set_target_properties(scirocco_viewer_app PROPERTIES
    OUTPUT_NAME "SciroccoViewer"
    MACOSX_BUNDLE TRUE
    MACOSX_BUNDLE_BUNDLE_NAME "Scirocco Viewer"
    MACOSX_BUNDLE_ICON_FILE "AppIcon"
)

# The Scirocco logo comes from the documentation site, so the viewer always
# shows the current one.
set(SCIROCCO_LOGO_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../docs/static/img)
set(VIEWER_RESOURCES resources/app-icon.svg)
foreach(logo logo logo-dark logo-full logo-full-dark)
    set_source_files_properties(${SCIROCCO_LOGO_DIR}/${logo}.svg PROPERTIES
        QT_RESOURCE_ALIAS resources/${logo}.svg)
    list(APPEND VIEWER_RESOURCES ${SCIROCCO_LOGO_DIR}/${logo}.svg)
endforeach()

# --- QML module ---
set_source_files_properties(qml/Theme.qml PROPERTIES QT_QML_SINGLETON_TYPE TRUE)

# The module's build directory would otherwise be named after its URI,
# SciroccoViewer, which clashes with the executable outside a macOS bundle.
qt_add_qml_module(scirocco_viewer_app
    URI SciroccoViewer
    VERSION 1.0
    OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/qml/SciroccoViewer
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
        qml/components/Icon.qml
    RESOURCES
        ${VIEWER_RESOURCES}
)

# --- Include & link ---
target_include_directories(scirocco_viewer_app PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(scirocco_viewer_app PRIVATE
    Qt6::Core Qt6::Gui Qt6::Qml Qt6::Quick Qt6::QuickControls2 Qt6::Svg
    Qt6::Widgets Qt6::Charts
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
)
