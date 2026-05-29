pragma Singleton
import QtQuick

QtObject {
    property bool darkMode: true

    // Surfaces
    readonly property color bgCanvas:    darkMode ? "#0e1116" : "#f5f5f7"
    readonly property color panelBg:     darkMode ? Qt.rgba(0.078, 0.086, 0.110, 0.72) : Qt.rgba(1, 1, 1, 0.85)
    readonly property color panelBgSolid:darkMode ? Qt.rgba(0.078, 0.086, 0.110, 1.0) : Qt.rgba(1, 1, 1, 1.0)
    readonly property color panelBorder: darkMode ? Qt.rgba(1, 1, 1, 0.08) : Qt.rgba(0, 0, 0, 0.10)
    readonly property color railBg:      darkMode ? Qt.rgba(0.047, 0.055, 0.071, 0.85) : Qt.rgba(0.95, 0.95, 0.96, 0.92)
    readonly property color statusBg:    darkMode ? Qt.rgba(0.047, 0.055, 0.071, 0.85) : Qt.rgba(0.95, 0.95, 0.96, 0.92)
    readonly property color hover:       darkMode ? Qt.rgba(1, 1, 1, 0.06) : Qt.rgba(0, 0, 0, 0.05)
    readonly property color active:      darkMode ? Qt.rgba(1, 1, 1, 0.10) : Qt.rgba(0, 0, 0, 0.08)
    readonly property color divider:     darkMode ? Qt.rgba(1, 1, 1, 0.06) : Qt.rgba(0, 0, 0, 0.08)
    readonly property color chartBg:     "transparent"

    // Text
    readonly property color textPrimary: darkMode ? "#e7ecf2" : "#1d1d1f"
    readonly property color textMuted:   darkMode ? "#8b94a3" : "#6e6e73"
    readonly property color textDim:     darkMode ? "#5b6473" : "#aeaeb2"

    // Accent
    readonly property color accent:      darkMode ? "#7cc4ff" : "#0071e3"
    readonly property color accentDim:   darkMode ? Qt.rgba(0.486, 0.769, 1.0, 0.18) : Qt.rgba(0, 0.443, 0.890, 0.12)

    // Lines
    readonly property color foilStroke:  darkMode ? "#e7ecf2" : "#1d1d1f"
    readonly property color wakeStroke:  darkMode ? "#5b6473" : "#aeaeb2"
    readonly property color cpLine:      darkMode ? "#7cc4ff" : "#0071e3"
    readonly property color cpiLine:     darkMode ? Qt.rgba(0.486, 0.769, 1.0, 0.4) : Qt.rgba(0, 0.443, 0.890, 0.35)
    readonly property color arrowStroke: darkMode ? Qt.rgba(0.486, 0.769, 1.0, 0.5) : Qt.rgba(0, 0.443, 0.890, 0.45)
    readonly property color gridLine:    darkMode ? Qt.rgba(1, 1, 1, 0.06) : Qt.rgba(0, 0, 0, 0.08)
    readonly property color axisLabel:   darkMode ? "#8b94a3" : "#6e6e73"

    // Typography
    readonly property string fontMono: "Menlo"
    readonly property int fsCaption:  11
    readonly property int fsBody:     12
    readonly property int fsControl:  13
    readonly property int fsTitle:    14

    // Spacing
    readonly property int s1:  4
    readonly property int s2:  8
    readonly property int s3: 12
    readonly property int s4: 16
    readonly property int s5: 24

    // Radius
    readonly property int rControl:  6
    readonly property int rPanel:   10
    readonly property int rChip:   999

    // Dimensions
    readonly property int railWidth:    52
    readonly property int panelWidth:  300
    readonly property int statusHeight: 36
}
