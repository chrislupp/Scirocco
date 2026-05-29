import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs

ApplicationWindow {
    id: root
    visible: true
    width: 1200
    height: 900
    title: reader.filePath ? "Scirocco Viewer — " + reader.filePath : "Scirocco Viewer"

    property bool showCpArrows: false

    menuBar: MenuBar {
        Menu {
            title: "File"
            Action {
                text: "Open..."
                shortcut: StandardKey.Open
                onTriggered: fileDialog.open()
            }
            MenuSeparator {}
            Action {
                text: "Quit"
                shortcut: StandardKey.Quit
                onTriggered: Qt.quit()
            }
        }
        Menu {
            title: "View"
            MenuItem {
                text: "Show Cp Arrows"
                checkable: true
                checked: root.showCpArrows
                onToggled: root.showCpArrows = checked
            }
        }
    }

    header: IterationControl {}

    FileDialog {
        id: fileDialog
        title: "Open HDF5 File"
        nameFilters: ["HDF5 files (*.h5 *.hdf5)", "All files (*)"]
        onAccepted: reader.openFile(selectedFile.toString().replace("file://", ""))
    }

    Shortcut {
        sequence: "Left"
        onActivated: reader.setCurrentIteration(reader.currentIteration - 1)
    }
    Shortcut {
        sequence: "Right"
        onActivated: reader.setCurrentIteration(reader.currentIteration + 1)
    }
    Shortcut {
        sequence: "Home"
        onActivated: reader.setCurrentIteration(0)
    }
    Shortcut {
        sequence: "End"
        onActivated: reader.setCurrentIteration(reader.numIterations - 1)
    }

    SplitView {
        anchors.fill: parent
        orientation: Qt.Vertical

        SplitView {
            SplitView.preferredHeight: parent.height * 0.35
            orientation: Qt.Horizontal

            AirfoilChart {
                SplitView.fillWidth: true
                showArrows: root.showCpArrows
            }

            CoeffPanel {
                SplitView.preferredWidth: 220
                SplitView.minimumWidth: 180
            }
        }

        CpChart {
            SplitView.preferredHeight: parent.height * 0.30
        }

        Item {
            SplitView.preferredHeight: parent.height * 0.35

            TabBar {
                id: bottomTabs
                anchors.top: parent.top
                anchors.left: parent.left
                anchors.right: parent.right

                TabButton { text: "Distributions" }
                TabButton { text: "Convergence" }
            }

            StackLayout {
                anchors.top: bottomTabs.bottom
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.bottom: parent.bottom
                currentIndex: bottomTabs.currentIndex

                DistributionChart {}
                ConvergenceChart {}
            }
        }
    }
}
