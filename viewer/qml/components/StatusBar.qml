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

        // Slider
        Slider {
            id: iterSlider
            Layout.preferredWidth: 160
            from: 0
            to: Math.max(reader.numIterations - 1, 0)
            stepSize: 1
            value: reader.currentIteration
            enabled: reader.numIterations > 1
            onMoved: reader.setCurrentIteration(Math.round(value))

            background: Rectangle {
                x: iterSlider.leftPadding
                y: iterSlider.topPadding + iterSlider.availableHeight / 2 - 1.5
                width: iterSlider.availableWidth
                height: 3
                radius: 2
                color: Theme.divider

                Rectangle {
                    width: iterSlider.visualPosition * parent.width
                    height: parent.height
                    radius: 2
                    color: Theme.accent
                }
            }

            handle: Rectangle {
                x: iterSlider.leftPadding + iterSlider.visualPosition * (iterSlider.availableWidth - width)
                y: iterSlider.topPadding + iterSlider.availableHeight / 2 - height / 2
                width: 14; height: 14
                radius: 7
                color: iterSlider.pressed ? Theme.accent : "#ffffff"
                border.width: 2
                border.color: Theme.accent
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
