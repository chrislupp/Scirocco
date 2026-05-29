pragma Singleton
import QtQuick

QtObject {
    // Surfaces
    readonly property color bgCanvas:    "#0e1116"
    readonly property color panelBg:     Qt.rgba(0.078, 0.086, 0.110, 0.72)
    readonly property color panelBgSolid:Qt.rgba(0.078, 0.086, 0.110, 1.0)
    readonly property color panelBorder: Qt.rgba(1, 1, 1, 0.08)
    readonly property color railBg:      Qt.rgba(0.047, 0.055, 0.071, 0.85)
    readonly property color statusBg:    Qt.rgba(0.047, 0.055, 0.071, 0.85)
    readonly property color hover:       Qt.rgba(1, 1, 1, 0.06)
    readonly property color active:      Qt.rgba(1, 1, 1, 0.10)
    readonly property color divider:     Qt.rgba(1, 1, 1, 0.06)
    readonly property color chartBg:     "transparent"

    // Text
    readonly property color textPrimary: "#e7ecf2"
    readonly property color textMuted:   "#8b94a3"
    readonly property color textDim:     "#5b6473"

    // Accent
    readonly property color accent:      "#7cc4ff"
    readonly property color accentDim:   Qt.rgba(0.486, 0.769, 1.0, 0.18)

    // Lines
    readonly property color foilStroke:  "#e7ecf2"
    readonly property color wakeStroke:  "#5b6473"
    readonly property color cpLine:      "#7cc4ff"
    readonly property color cpiLine:     Qt.rgba(0.486, 0.769, 1.0, 0.4)
    readonly property color arrowStroke: Qt.rgba(0.486, 0.769, 1.0, 0.5)
    readonly property color gridLine:    Qt.rgba(1, 1, 1, 0.06)
    readonly property color axisLabel:   "#8b94a3"

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
