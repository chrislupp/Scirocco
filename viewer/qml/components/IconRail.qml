import QtQuick
import QtQuick.Layouts
import SciroccoViewer

Rectangle {
    id: root
    width: Theme.railWidth
    color: Theme.railBg

    property int activeSection: -1
    property bool arrowsEnabled: false

    signal sectionClicked(int index)
    signal arrowsToggled(bool on)

    Rectangle {
        anchors.top: parent.top
        anchors.bottom: parent.bottom
        anchors.right: parent.right
        width: 1
        color: Theme.panelBorder
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.topMargin: Theme.s3
        anchors.bottomMargin: Theme.s3
        spacing: Theme.s1

        // Brand mark
        Item {
            Layout.alignment: Qt.AlignHCenter
            Layout.preferredWidth: 32
            Layout.preferredHeight: 32

            Rectangle {
                anchors.centerIn: parent
                width: 22; height: 22
                radius: 4
                color: "transparent"
                border.width: 1.5
                border.color: Theme.accent
            }
            Rectangle {
                anchors.centerIn: parent
                width: 8; height: 8
                radius: 4
                color: Theme.accent
            }
        }

        Item { Layout.preferredHeight: Theme.s3 }

        // Airfoil view
        RailButton {
            icon: "airfoil"
            checked: root.activeSection === 0
            onClicked: root.sectionClicked(root.activeSection === 0 ? -1 : 0)
        }

        // Charts
        RailButton {
            icon: "chart"
            checked: root.activeSection === 1
            onClicked: root.sectionClicked(root.activeSection === 1 ? -1 : 1)
        }

        // Data inspector
        RailButton {
            icon: "data"
            checked: root.activeSection === 2
            onClicked: root.sectionClicked(root.activeSection === 2 ? -1 : 2)
        }

        Item { Layout.fillHeight: true }

        // Arrows toggle
        RailButton {
            icon: "arrows"
            checked: root.arrowsEnabled
            onClicked: {
                root.arrowsEnabled = !root.arrowsEnabled
                root.arrowsToggled(root.arrowsEnabled)
            }
        }
    }

    component RailButton: Rectangle {
        property string icon
        property bool checked: false
        signal clicked()
        onCheckedChanged: iconCanvas.requestPaint()

        Layout.alignment: Qt.AlignHCenter
        Layout.preferredWidth: 36
        Layout.preferredHeight: 36
        radius: Theme.rControl
        color: checked ? Theme.active : mouseArea.containsMouse ? Theme.hover : "transparent"

        Canvas {
            id: iconCanvas
            anchors.centerIn: parent
            width: 20; height: 20
            onPaint: {
                var ctx = getContext("2d")
                ctx.clearRect(0, 0, width, height)
                ctx.strokeStyle = parent.checked ? Theme.accent.toString() : Theme.textMuted.toString()
                ctx.lineWidth = 1.5
                ctx.lineCap = "round"
                ctx.lineJoin = "round"

                if (parent.icon === "airfoil") {
                    // Simplified airfoil shape
                    ctx.beginPath()
                    ctx.moveTo(18, 10)
                    ctx.bezierCurveTo(14, 4, 4, 5, 2, 10)
                    ctx.bezierCurveTo(4, 14, 14, 15, 18, 10)
                    ctx.stroke()
                } else if (parent.icon === "chart") {
                    // Line chart icon
                    ctx.beginPath()
                    ctx.moveTo(2, 16)
                    ctx.lineTo(7, 8)
                    ctx.lineTo(12, 12)
                    ctx.lineTo(18, 4)
                    ctx.stroke()
                    ctx.beginPath()
                    ctx.moveTo(2, 18)
                    ctx.lineTo(18, 18)
                    ctx.stroke()
                } else if (parent.icon === "data") {
                    // Table/grid icon
                    ctx.beginPath()
                    ctx.moveTo(2, 4); ctx.lineTo(18, 4)
                    ctx.moveTo(2, 10); ctx.lineTo(18, 10)
                    ctx.moveTo(2, 16); ctx.lineTo(18, 16)
                    ctx.moveTo(10, 4); ctx.lineTo(10, 16)
                    ctx.stroke()
                } else if (parent.icon === "arrows") {
                    // Arrows icon
                    ctx.beginPath()
                    ctx.moveTo(10, 16); ctx.lineTo(10, 4)
                    ctx.moveTo(6, 7); ctx.lineTo(10, 4); ctx.lineTo(14, 7)
                    ctx.stroke()
                    ctx.beginPath()
                    ctx.moveTo(4, 16); ctx.lineTo(16, 16)
                    ctx.stroke()
                }
            }

            Component.onCompleted: requestPaint()
        }

        MouseArea {
            id: mouseArea
            anchors.fill: parent
            hoverEnabled: true
            cursorShape: Qt.PointingHandCursor
            onClicked: parent.clicked()
        }
    }
}
