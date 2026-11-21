import QtQuick
import QtQuick.Controls
import SciroccoViewer

AbstractButton {
    id: root
    checkable: true
    implicitHeight: 26
    implicitWidth: contentRow.implicitWidth + Theme.s3 * 2

    contentItem: Row {
        id: contentRow
        spacing: Theme.s1 + 2
        anchors.centerIn: parent

        Rectangle {
            width: 6; height: 6; radius: 3
            anchors.verticalCenter: parent.verticalCenter
            color: root.checked ? Theme.accent : Theme.textDim
        }

        Text {
            text: root.text
            font.pixelSize: Theme.fsBody
            font.weight: root.checked ? Font.DemiBold : Font.Medium
            color: root.checked ? Theme.textPrimary : Theme.textMuted
            anchors.verticalCenter: parent.verticalCenter
        }
    }

    background: Rectangle {
        radius: Theme.rChip
        color: root.checked ? Theme.accentDim
             : root.hovered ? Theme.hover
             : Qt.rgba(1, 1, 1, 0.025)
        border.width: 1
        border.color: root.checked ? Theme.accent : Theme.panelBorder
        opacity: root.enabled ? 1.0 : 0.45
    }
}
