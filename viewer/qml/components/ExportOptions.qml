import QtQuick
import QtQuick.Layouts
import SciroccoViewer

Item {
    id: root
    width: 280
    height: col.implicitHeight + 2 * Theme.s4
    visible: false

    property bool exportDark: Theme.darkMode
    property int exportScale: 3
    property int exportFps: 10

    signal exportClicked()

    function show() {
        exportDark = Theme.darkMode
        visible = true
    }

    // Dim background
    Rectangle {
        anchors.fill: parent
        parent: root.parent
        visible: root.visible
        color: Qt.rgba(0, 0, 0, 0.4)
        z: 9

        MouseArea {
            anchors.fill: parent
            onClicked: root.visible = false
        }
    }

    Panel {
        anchors.fill: parent

        ColumnLayout {
            id: col
            anchors.fill: parent
            anchors.margins: Theme.s4
            spacing: Theme.s4

            Text {
                text: "Export Video"
                font.pixelSize: Theme.fsTitle
                font.weight: Font.DemiBold
                color: Theme.textPrimary
            }

            Rectangle { Layout.fillWidth: true; height: 1; color: Theme.divider }

            // Theme
            ColumnLayout {
                spacing: Theme.s2
                SectionHeader { text: "Theme" }
                Row {
                    spacing: Theme.s2
                    Chip {
                        text: "Dark"
                        checked: root.exportDark
                        onClicked: root.exportDark = true
                    }
                    Chip {
                        text: "Light"
                        checked: !root.exportDark
                        onClicked: root.exportDark = false
                    }
                }
            }

            // Resolution
            ColumnLayout {
                spacing: Theme.s2
                SectionHeader { text: "Resolution" }
                Row {
                    spacing: Theme.s2
                    Chip {
                        text: "2x"
                        checked: root.exportScale === 2
                        onClicked: root.exportScale = 2
                    }
                    Chip {
                        text: "3x"
                        checked: root.exportScale === 3
                        onClicked: root.exportScale = 3
                    }
                    Chip {
                        text: "4x"
                        checked: root.exportScale === 4
                        onClicked: root.exportScale = 4
                    }
                }
                Text {
                    property int w: Math.floor(mainArea.width * root.exportScale / 2) * 2
                    property int h: Math.floor(mainArea.height * root.exportScale / 2) * 2
                    text: w + " x " + h + " px"
                    font.family: Theme.fontMono
                    font.pixelSize: Theme.fsCaption
                    color: Theme.textMuted
                }
            }

            // FPS
            ColumnLayout {
                spacing: Theme.s2
                SectionHeader { text: "Frame Rate" }
                Row {
                    spacing: Theme.s2
                    Chip {
                        text: "5 fps"
                        checked: root.exportFps === 5
                        onClicked: root.exportFps = 5
                    }
                    Chip {
                        text: "10 fps"
                        checked: root.exportFps === 10
                        onClicked: root.exportFps = 10
                    }
                    Chip {
                        text: "24 fps"
                        checked: root.exportFps === 24
                        onClicked: root.exportFps = 24
                    }
                }
            }

            Rectangle { Layout.fillWidth: true; height: 1; color: Theme.divider }

            // Export button
            Rectangle {
                Layout.fillWidth: true
                height: 32
                radius: Theme.rControl
                color: exportMa.containsMouse ? Qt.darker(Theme.accent, 1.1) : Theme.accent

                Text {
                    anchors.centerIn: parent
                    text: "Export"
                    font.pixelSize: Theme.fsControl
                    font.weight: Font.DemiBold
                    color: "#ffffff"
                }

                MouseArea {
                    id: exportMa
                    anchors.fill: parent
                    hoverEnabled: true
                    cursorShape: Qt.PointingHandCursor
                    onClicked: {
                        root.visible = false
                        root.exportClicked()
                    }
                }
            }
        }
    }
}
