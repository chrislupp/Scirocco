import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import SciroccoViewer

Rectangle {
    id: root
    height: Theme.statusHeight
    color: Theme.statusBg

    function togglePlay() {
        if (playTimer.running) {
            playTimer.stop()
        } else {
            if (reader.currentIteration >= reader.numIterations - 1)
                reader.setCurrentIteration(0)
            playTimer.start()
        }
        playIcon.requestPaint()
    }

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

        // Play/Pause
        Rectangle {
            Layout.preferredWidth: 24; Layout.preferredHeight: 24
            radius: Theme.rControl
            color: playMa.containsMouse ? Theme.hover : "transparent"
            opacity: reader.numIterations > 1 ? 1.0 : 0.3

            Canvas {
                id: playIcon
                anchors.centerIn: parent
                width: 12; height: 12
                onPaint: {
                    var ctx = getContext("2d")
                    ctx.clearRect(0, 0, width, height)
                    ctx.fillStyle = playTimer.running ? Theme.accent.toString() : Theme.textMuted.toString()
                    if (playTimer.running) {
                        ctx.fillRect(1, 1, 3.5, 10)
                        ctx.fillRect(7.5, 1, 3.5, 10)
                    } else {
                        ctx.beginPath()
                        ctx.moveTo(1, 0)
                        ctx.lineTo(12, 6)
                        ctx.lineTo(1, 12)
                        ctx.closePath()
                        ctx.fill()
                    }
                }
            }

            MouseArea {
                id: playMa
                anchors.fill: parent
                hoverEnabled: true
                cursorShape: Qt.PointingHandCursor
                onClicked: {
                    if (playTimer.running) {
                        playTimer.stop()
                    } else {
                        if (reader.currentIteration >= reader.numIterations - 1)
                            reader.setCurrentIteration(0)
                        playTimer.start()
                    }
                    playIcon.requestPaint()
                }
            }
        }

        Timer {
            id: playTimer
            interval: 200
            repeat: true
            onTriggered: {
                if (reader.currentIteration < reader.numIterations - 1) {
                    reader.setCurrentIteration(reader.currentIteration + 1)
                } else {
                    stop()
                    playIcon.requestPaint()
                }
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

        // Export progress
        Row {
            spacing: Theme.s2
            visible: videoExporter.exporting

            Text {
                text: "Exporting " + videoExporter.progress + "%"
                font.family: Theme.fontMono
                font.pixelSize: Theme.fsCaption
                color: Theme.accent
                anchors.verticalCenter: parent.verticalCenter
            }
            Rectangle {
                width: 60; height: 3; radius: 2
                color: Theme.divider
                anchors.verticalCenter: parent.verticalCenter
                Rectangle {
                    width: parent.width * videoExporter.progress / 100
                    height: parent.height; radius: 2
                    color: Theme.accent
                }
            }
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
