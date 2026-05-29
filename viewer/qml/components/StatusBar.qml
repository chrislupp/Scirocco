import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import SciroccoViewer

Rectangle {
    id: root
    height: Theme.statusHeight
    color: Theme.statusBg

    Rectangle {
        anchors.top: parent.top
        anchors.left: parent.left
        anchors.right: parent.right
        height: 1
        color: Theme.panelBorder
    }

    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: Theme.s4
        anchors.rightMargin: Theme.s4
        spacing: Theme.s3

        // Prev
        Rectangle {
            Layout.preferredWidth: 24; Layout.preferredHeight: 24
            radius: Theme.rControl
            color: prevMa.containsMouse ? Theme.hover : "transparent"
            opacity: reader.currentIteration > 0 ? 1.0 : 0.3

            Text {
                anchors.centerIn: parent
                text: "◀"
                font.pixelSize: 10
                color: Theme.textMuted
            }
            MouseArea {
                id: prevMa
                anchors.fill: parent
                hoverEnabled: true
                cursorShape: Qt.PointingHandCursor
                onClicked: reader.setCurrentIteration(reader.currentIteration - 1)
            }
        }

        // Iteration scrubber
        Item {
            id: scrubber
            Layout.preferredWidth: 160
            Layout.preferredHeight: Theme.statusHeight - 4

            property int maxVal: Math.max(reader.numIterations - 1, 0)
            property real pos: maxVal > 0 ? reader.currentIteration / maxVal : 0

            Rectangle {
                id: track
                anchors.verticalCenter: parent.verticalCenter
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.leftMargin: 7
                anchors.rightMargin: 7
                height: 3
                radius: 2
                color: Theme.divider

                Rectangle {
                    width: scrubber.pos * parent.width
                    height: parent.height
                    radius: 2
                    color: Theme.accent
                }
            }

            Rectangle {
                id: knob
                width: 14; height: 14; radius: 7
                y: (scrubber.height - height) / 2
                x: 7 + scrubber.pos * track.width - 7
                color: scrubDrag.drag.active ? Theme.accent : "#ffffff"
                border.width: 2
                border.color: Theme.accent

                Behavior on color { ColorAnimation { duration: 120 } }
            }

            MouseArea {
                id: scrubDrag
                anchors.fill: parent
                cursorShape: Qt.PointingHandCursor
                enabled: reader.numIterations > 1

                property bool dragging: false

                onPressed: function(mouse) {
                    dragging = true
                    updateFromMouse(mouse.x)
                }
                onReleased: { dragging = false }
                onPositionChanged: function(mouse) {
                    if (dragging) updateFromMouse(mouse.x)
                }

                function updateFromMouse(mx) {
                    var frac = (mx - 7) / track.width
                    frac = Math.max(0, Math.min(1, frac))
                    var idx = Math.round(frac * scrubber.maxVal)
                    reader.setCurrentIteration(idx)
                }
            }
        }

        // Next
        Rectangle {
            Layout.preferredWidth: 24; Layout.preferredHeight: 24
            radius: Theme.rControl
            color: nextMa.containsMouse ? Theme.hover : "transparent"
            opacity: reader.currentIteration < reader.numIterations - 1 ? 1.0 : 0.3

            Text {
                anchors.centerIn: parent
                text: "▶"
                font.pixelSize: 10
                color: Theme.textMuted
            }
            MouseArea {
                id: nextMa
                anchors.fill: parent
                hoverEnabled: true
                cursorShape: Qt.PointingHandCursor
                onClicked: reader.setCurrentIteration(reader.currentIteration + 1)
            }
        }

        // Iteration label
        Text {
            text: reader.numIterations > 0
                  ? reader.currentIteration + " / " + (reader.numIterations - 1)
                  : "—"
            font.family: Theme.fontMono
            font.pixelSize: Theme.fsCaption
            color: Theme.textMuted
        }

        // Divider
        Rectangle { Layout.preferredWidth: 1; Layout.preferredHeight: 14; color: Theme.divider }

        // Key readout
        Row {
            spacing: Theme.s4
            visible: reader.numIterations > 0

            StatusValue { label: "α"; value: reader.data.alpha.toFixed(2) + "°" }
            StatusValue { label: "Cl"; value: reader.data.cl.toFixed(4) }
            StatusValue { label: "Cd"; value: reader.data.cd.toFixed(5) }
        }

        Item { Layout.fillWidth: true }

        // Convergence dot
        Rectangle {
            width: 6; height: 6; radius: 3
            color: reader.numIterations > 0 && reader.data.converged ? Theme.accent : Theme.textDim
            visible: reader.numIterations > 0
        }

        // Divider
        Rectangle { Layout.preferredWidth: 1; Layout.preferredHeight: 14; color: Theme.divider }

        // App name
        Text {
            text: "Scirocco Viewer"
            font.family: Theme.fontMono
            font.pixelSize: Theme.fsCaption
            font.letterSpacing: 1.2
            color: Theme.textDim
        }
    }

    component StatusValue: Row {
        property string label
        property string value
        spacing: Theme.s1

        Text {
            text: label
            font.pixelSize: Theme.fsCaption
            color: Theme.textDim
            anchors.verticalCenter: parent.verticalCenter
        }
        Text {
            text: value
            font.family: Theme.fontMono
            font.pixelSize: Theme.fsCaption
            color: Theme.textPrimary
            anchors.verticalCenter: parent.verticalCenter
        }
    }
}
