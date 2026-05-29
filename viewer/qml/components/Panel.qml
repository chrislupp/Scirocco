import QtQuick
import SciroccoViewer

Rectangle {
    id: root
    color: Theme.panelBg
    radius: Theme.rPanel
    border.width: 1
    border.color: Theme.panelBorder

    Rectangle {
        anchors.fill: parent
        anchors.margins: -6
        z: -2
        radius: root.radius + 4
        color: "transparent"
        border.width: 1
        border.color: Qt.rgba(0, 0, 0, 0.18)
        opacity: 0.5
    }

    Rectangle {
        anchors.fill: parent
        anchors.margins: -14
        z: -3
        radius: root.radius + 8
        color: "transparent"
        border.width: 1
        border.color: Qt.rgba(0, 0, 0, 0.10)
        opacity: 0.6
    }
}
