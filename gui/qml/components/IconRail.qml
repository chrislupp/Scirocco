import QtQuick
import QtQuick.Layouts
import SciroccoViewer

Rectangle {
    id: root
    width: Theme.railWidth
    color: Theme.railBg

    property int activeSection: -1
    property bool arrowsEnabled: false
    property bool panelsEnabled: false
    property bool hudEnabled: false

    signal sectionClicked(int index)
    signal arrowsMenuRequested(real buttonY)
    signal panelsToggled(bool on)
    signal hudMenuRequested(real buttonY)
    signal exportRequested()

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
        Image {
            Layout.alignment: Qt.AlignHCenter
            Layout.preferredWidth: 34
            Layout.preferredHeight: 26
            source: Theme.darkMode
                    ? "qrc:/qt/qml/SciroccoViewer/resources/logo-dark.svg"
                    : "qrc:/qt/qml/SciroccoViewer/resources/logo.svg"
            sourceSize: Qt.size(34, 26)
            fillMode: Image.PreserveAspectFit
            smooth: true
        }

        Item { Layout.preferredHeight: Theme.s3 }

        // Airfoil view
        RailButton {
            icon: "airfoil"
            tip: "Airfoil only"
            checked: root.activeSection === 0
            onClicked: root.sectionClicked(root.activeSection === 0 ? -1 : 0)
        }

        // Charts
        RailButton {
            icon: "chart"
            tip: "Charts"
            checked: root.activeSection === 1
            onClicked: root.sectionClicked(root.activeSection === 1 ? -1 : 1)
        }

        // Data inspector
        RailButton {
            icon: "results"
            tip: "Results"
            checked: root.activeSection === 2
            onClicked: root.sectionClicked(root.activeSection === 2 ? -1 : 2)
        }

        Item { Layout.fillHeight: true }

        // HUD menu
        RailButton {
            id: hudButton
            icon: "hud"
            tip: "HUD"
            checked: root.hudEnabled
            onClicked: {
                var pos = hudButton.mapToItem(null, 0, 0)
                root.hudMenuRequested(pos.y)
            }
        }

        // Panels toggle
        RailButton {
            icon: "panels"
            tip: "Panel nodes"
            checked: root.panelsEnabled
            onClicked: {
                root.panelsEnabled = !root.panelsEnabled
                root.panelsToggled(root.panelsEnabled)
            }
        }

        // Arrows menu
        RailButton {
            id: arrowsButton
            icon: "arrows"
            tip: "Cp arrows"
            checked: root.arrowsEnabled
            onClicked: {
                var pos = arrowsButton.mapToItem(null, 0, 0)
                root.arrowsMenuRequested(pos.y)
            }
        }

        // Export video
        RailButton {
            icon: "video"
            tip: "Export video"
            checked: false
            onClicked: root.exportRequested()
        }

        Item { Layout.preferredHeight: Theme.s2 }

        // Theme toggle
        RailButton {
            icon: Theme.darkMode ? "sun" : "moon"
            tip: Theme.darkMode ? "Light theme" : "Dark theme"
            checked: false
            onClicked: Theme.darkMode = !Theme.darkMode
        }
    }

    component RailButton: Rectangle {
        id: button
        property string icon
        property string tip
        property bool checked: false
        signal clicked()
        signal rightClicked()

        Layout.alignment: Qt.AlignHCenter
        Layout.preferredWidth: 36
        Layout.preferredHeight: 36
        radius: Theme.rControl
        color: checked ? Theme.active : mouseArea.containsMouse ? Theme.hover : "transparent"

        Behavior on color { ColorAnimation { duration: 120 } }

        // Accent bar marking the active item.
        Rectangle {
            anchors.verticalCenter: parent.verticalCenter
            x: -(Theme.railWidth - parent.width) / 2
            width: 3
            height: button.checked ? 18 : 0
            radius: 1.5
            color: Theme.accent
            Behavior on height { NumberAnimation { duration: 160; easing.type: Easing.OutCubic } }
        }

        Icon {
            anchors.centerIn: parent
            name: button.icon
            size: 20
            color: button.checked ? Theme.accent
                 : mouseArea.containsMouse ? Theme.textPrimary
                 : Theme.textMuted
            Behavior on color { ColorAnimation { duration: 120 } }
        }

        MouseArea {
            id: mouseArea
            anchors.fill: parent
            acceptedButtons: Qt.LeftButton | Qt.RightButton
            hoverEnabled: true
            cursorShape: Qt.PointingHandCursor
            onClicked: function(mouse) {
                if (mouse.button === Qt.RightButton)
                    parent.rightClicked()
                else
                    parent.clicked()
            }
        }

        // Label shown to the right of the rail on hover.
        Rectangle {
            visible: tipTimer.shown && button.tip !== ""
            x: parent.width + (Theme.railWidth - parent.width) / 2 + Theme.s2
            anchors.verticalCenter: parent.verticalCenter
            width: tipText.implicitWidth + 2 * Theme.s2
            height: tipText.implicitHeight + Theme.s1 * 2
            radius: Theme.rControl
            color: Theme.panelBgSolid
            border.width: 1
            border.color: Theme.panelBorder
            z: 100

            Text {
                id: tipText
                anchors.centerIn: parent
                text: button.tip
                font.pixelSize: Theme.fsCaption
                color: Theme.textPrimary
            }
        }

        Timer {
            id: tipTimer
            property bool shown: false
            interval: 450
            running: mouseArea.containsMouse && !shown
            onTriggered: shown = true
        }
        Connections {
            target: mouseArea
            function onContainsMouseChanged() {
                if (!mouseArea.containsMouse) tipTimer.shown = false
            }
        }
    }
}
