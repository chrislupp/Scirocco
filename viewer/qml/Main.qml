import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import SciroccoViewer

ApplicationWindow {
    id: root
    visible: true
    width: 1200
    height: 800
    minimumWidth: 720
    minimumHeight: 480
    title: "Scirocco Viewer"
    color: Theme.bgCanvas

    // -1 = airfoil only, 0 = not used, 1 = charts, 2 = inspector
    property int activeSection: -1
    property bool showCpArrows: false
    property string activeChart: "cp"

    // Keyboard shortcuts
    Shortcut { sequence: "Left";  onActivated: reader.setCurrentIteration(reader.currentIteration - 1) }
    Shortcut { sequence: "Right"; onActivated: reader.setCurrentIteration(reader.currentIteration + 1) }
    Shortcut { sequence: "Home";  onActivated: reader.setCurrentIteration(0) }
    Shortcut { sequence: "End";   onActivated: reader.setCurrentIteration(reader.numIterations - 1) }
    Shortcut { sequences: [StandardKey.Open]; onActivated: fileDialog.open() }

    FileDialog {
        id: fileDialog
        title: "Open HDF5 File"
        nameFilters: ["HDF5 files (*.h5 *.hdf5)", "All files (*)"]
        onAccepted: reader.openFile(selectedFile.toString().replace("file://", ""))
    }

    // --- Layout ---

    IconRail {
        id: rail
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.bottom: statusBar.top
        activeSection: root.activeSection
        arrowsEnabled: root.showCpArrows

        onSectionClicked: function(idx) { root.activeSection = idx }
        onArrowsToggled: function(on) { root.showCpArrows = on }
    }

    // Main content area (right of rail, above status bar)
    Item {
        id: mainArea
        anchors.left: rail.right
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.bottom: statusBar.top

        // Airfoil canvas (always visible, fills available space)
        AirfoilCanvas {
            anchors.left: parent.left
            anchors.right: parent.right
            anchors.top: parent.top
            anchors.bottom: chartPanel.visible ? chartPanel.top : parent.bottom
            showArrows: root.showCpArrows
        }

        // Chart panel (bottom portion, slides up when active)
        Rectangle {
            id: chartPanel
            anchors.left: parent.left
            anchors.right: parent.right
            anchors.bottom: parent.bottom
            height: root.activeSection === 1 ? mainArea.height * 0.45 : 0
            visible: height > 0
            color: "transparent"
            clip: true

            Behavior on height {
                NumberAnimation { duration: 280; easing.type: Easing.OutCubic }
            }

            // Top border
            Rectangle {
                anchors.top: parent.top
                anchors.left: parent.left
                anchors.right: parent.right
                height: 1
                color: Theme.panelBorder
            }

            // Chart chip row
            Flow {
                id: chartChips
                anchors.top: parent.top
                anchors.topMargin: Theme.s2 + 1
                anchors.left: parent.left
                anchors.leftMargin: Theme.s3
                spacing: Theme.s2
                z: 1

                Chip { text: "Cp";      checked: root.activeChart === "cp";      onClicked: root.activeChart = "cp" }
                Chip { text: "Distributions"; checked: root.activeChart === "dist";  onClicked: root.activeChart = "dist" }
                Chip { text: "History"; checked: root.activeChart === "history"; onClicked: root.activeChart = "history" }
            }

            // Chart content
            StackLayout {
                anchors.top: chartChips.bottom
                anchors.topMargin: Theme.s2
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.bottom: parent.bottom
                currentIndex: root.activeChart === "cp" ? 0
                            : root.activeChart === "dist" ? 1 : 2

                CpChart {}
                DistributionChart {}
                ConvergenceChart {}
            }
        }

        // Floating inspector panel
        InspectorPanel {
            id: inspector
            anchors.top: parent.top
            anchors.topMargin: Theme.s3
            anchors.bottom: chartPanel.visible ? chartPanel.top : parent.bottom
            anchors.bottomMargin: Theme.s3
            x: root.activeSection === 2 ? Theme.s3 : -width - Theme.s3
            visible: x > -width
            onCloseRequested: root.activeSection = -1

            Behavior on x {
                NumberAnimation { duration: 280; easing.type: Easing.OutCubic }
            }
        }
    }

    StatusBar {
        id: statusBar
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.bottom: parent.bottom
    }

    // Open file dialog on double-click if no file loaded
    MouseArea {
        anchors.fill: mainArea
        visible: reader.numIterations === 0
        onDoubleClicked: fileDialog.open()
    }
}
