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
    property bool showPanels: false
    property string activeChart: "cp"
    property bool showHud: false
    property int hudFontSize: 14
    property bool showHudPanel: false

    // Arrow settings
    property real arrowLineWidth: 1.5
    property real arrowScale: 1.0
    property string arrowColormap: "coolwarm"
    property bool showArrowsPanel: false

    // Keyboard shortcuts
    Shortcut { sequence: "Left";  onActivated: reader.setCurrentIteration(reader.currentIteration - 1) }
    Shortcut { sequence: "Right"; onActivated: reader.setCurrentIteration(reader.currentIteration + 1) }
    Shortcut { sequence: "Home";  onActivated: reader.setCurrentIteration(0) }
    Shortcut { sequence: "End";   onActivated: reader.setCurrentIteration(reader.numIterations - 1) }
    Shortcut { sequences: [StandardKey.Open]; onActivated: fileDialog.open() }
    Shortcut { sequence: "Space"; onActivated: statusBar.togglePlay() }

    FileDialog {
        id: fileDialog
        title: "Open HDF5 File"
        nameFilters: ["HDF5 files (*.h5 *.hdf5)", "All files (*)"]
        onAccepted: reader.openFile(selectedFile.toString().replace("file://", ""))
    }

    property bool savedDarkMode: true

    FileDialog {
        id: exportFileDialog
        title: "Export Video"
        fileMode: FileDialog.SaveFile
        nameFilters: ["MP4 video (*.mp4)"]
        onAccepted: {
            var path = selectedFile.toString().replace("file://", "")
            root.savedDarkMode = Theme.darkMode
            Theme.darkMode = exportOptions.exportDark
            videoExporter.exportVideo(mainArea, path,
                                     exportOptions.exportDark,
                                     exportOptions.exportFps,
                                     exportOptions.exportScale)
        }
    }

    Connections {
        target: videoExporter
        function onFinished(path) { Theme.darkMode = root.savedDarkMode }
        function onError(msg) { Theme.darkMode = root.savedDarkMode }
    }

    Shortcut { sequence: "Ctrl+E"; onActivated: exportOptions.show() }
    Shortcut { sequence: "Meta+E"; onActivated: exportOptions.show() }

    // Export options overlay
    ExportOptions {
        id: exportOptions
        anchors.centerIn: parent
        z: 10
        onExportClicked: exportFileDialog.open()
    }

    // --- Layout ---

    IconRail {
        id: rail
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.bottom: statusBar.top
        activeSection: root.activeSection
        arrowsEnabled: root.showCpArrows
        panelsEnabled: root.showPanels

        hudEnabled: root.showHud

        onSectionClicked: function(idx) { root.activeSection = idx }
        onArrowsMenuRequested: function(buttonY) {
            arrowsPanel.panelY = buttonY
            root.showArrowsPanel = !root.showArrowsPanel
        }
        onPanelsToggled: function(on) { root.showPanels = on }
        onHudMenuRequested: function(buttonY) {
            hudPanel.panelY = buttonY
            root.showHudPanel = !root.showHudPanel
        }
        onExportRequested: exportOptions.show()
    }

    // Main content area (right of rail, above status bar)
    Rectangle {
        id: mainArea
        color: Theme.bgCanvas
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
            showPanels: root.showPanels
            arrowLineWidth: root.arrowLineWidth
            arrowScale: root.arrowScale
            arrowColormap: root.arrowColormap
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

        // HUD overlay (top-right)
        Column {
            id: hudOverlay
            anchors.top: parent.top
            anchors.right: parent.right
            anchors.topMargin: Theme.s3
            anchors.rightMargin: Theme.s3
            spacing: 2
            visible: root.showHud && reader.numIterations > 0
            z: 1

            Text {
                anchors.right: parent.right
                text: "Iteration " + reader.currentIteration
                font.family: Theme.fontMono
                font.pixelSize: root.hudFontSize
                color: Theme.textMuted
            }
            Text {
                anchors.right: parent.right
                text: "α = " + reader.data.alpha.toFixed(2) + "°"
                font.family: Theme.fontMono
                font.pixelSize: root.hudFontSize
                color: Theme.textPrimary
            }
            Text {
                anchors.right: parent.right
                text: "Cl = " + reader.data.cl.toFixed(4)
                font.family: Theme.fontMono
                font.pixelSize: root.hudFontSize
                color: Theme.textPrimary
            }
            Text {
                anchors.right: parent.right
                text: "Cd = " + reader.data.cd.toFixed(5)
                font.family: Theme.fontMono
                font.pixelSize: root.hudFontSize
                color: Theme.textPrimary
            }
        }

        // HUD settings panel (slides out from rail, aligned to icon)
        Panel {
            id: hudPanel
            property real panelY: 0
            width: 200
            height: hudPanelCol.implicitHeight + 2 * Theme.s4
            y: Math.min(panelY, parent.height - height - Theme.s3)
            x: root.showHudPanel ? Theme.s3 : -width - Theme.s3
            visible: x > -width
            z: 2

            Behavior on x {
                NumberAnimation { duration: 280; easing.type: Easing.OutCubic }
            }

            ColumnLayout {
                id: hudPanelCol
                anchors.fill: parent
                anchors.margins: Theme.s4
                spacing: Theme.s3

                RowLayout {
                    Layout.fillWidth: true
                    Text {
                        text: "HUD"
                        font.pixelSize: Theme.fsTitle
                        font.weight: Font.DemiBold
                        color: Theme.textPrimary
                        Layout.fillWidth: true
                    }
                    Rectangle {
                        width: 28; height: 28
                        radius: Theme.rControl
                        color: hudCloseMa.containsMouse ? Theme.hover : "transparent"
                        Text {
                            anchors.centerIn: parent
                            text: "✕"
                            font.pixelSize: 12
                            color: Theme.textMuted
                        }
                        MouseArea {
                            id: hudCloseMa
                            anchors.fill: parent
                            hoverEnabled: true
                            cursorShape: Qt.PointingHandCursor
                            onClicked: root.showHudPanel = false
                        }
                    }
                }

                Rectangle { Layout.fillWidth: true; height: 1; color: Theme.divider }

                // Show HUD toggle
                RowLayout {
                    Layout.fillWidth: true
                    Text {
                        text: "Visible"
                        font.pixelSize: Theme.fsBody
                        color: Theme.textPrimary
                        Layout.fillWidth: true
                    }
                    Chip {
                        text: root.showHud ? "On" : "Off"
                        checked: root.showHud
                        onClicked: root.showHud = !root.showHud
                    }
                }

                Rectangle { Layout.fillWidth: true; height: 1; color: Theme.divider }

                SectionHeader { text: "Text Size" }

                Flow {
                    Layout.fillWidth: true
                    spacing: Theme.s2
                    Chip { text: "S";  checked: root.hudFontSize === 10; onClicked: root.hudFontSize = 10 }
                    Chip { text: "M";  checked: root.hudFontSize === 14; onClicked: root.hudFontSize = 14 }
                    Chip { text: "L";  checked: root.hudFontSize === 20; onClicked: root.hudFontSize = 20 }
                    Chip { text: "XL"; checked: root.hudFontSize === 28; onClicked: root.hudFontSize = 28 }
                }
            }
        }

        // Arrows settings panel (slides out from rail, aligned to icon)
        Panel {
            id: arrowsPanel
            property real panelY: 0
            width: 220
            height: arrowsPanelCol.implicitHeight + 2 * Theme.s4
            y: Math.min(panelY, parent.height - height - Theme.s3)
            x: root.showArrowsPanel ? Theme.s3 : -width - Theme.s3
            visible: x > -width
            z: 2

            Behavior on x {
                NumberAnimation { duration: 280; easing.type: Easing.OutCubic }
            }

            ColumnLayout {
                id: arrowsPanelCol
                anchors.fill: parent
                anchors.margins: Theme.s4
                spacing: Theme.s3

                RowLayout {
                    Layout.fillWidth: true
                    Text {
                        text: "Cp Arrows"
                        font.pixelSize: Theme.fsTitle
                        font.weight: Font.DemiBold
                        color: Theme.textPrimary
                        Layout.fillWidth: true
                    }
                    Rectangle {
                        width: 28; height: 28
                        radius: Theme.rControl
                        color: arrowsCloseMa.containsMouse ? Theme.hover : "transparent"
                        Text {
                            anchors.centerIn: parent
                            text: "✕"
                            font.pixelSize: 12
                            color: Theme.textMuted
                        }
                        MouseArea {
                            id: arrowsCloseMa
                            anchors.fill: parent
                            hoverEnabled: true
                            cursorShape: Qt.PointingHandCursor
                            onClicked: root.showArrowsPanel = false
                        }
                    }
                }

                Rectangle { Layout.fillWidth: true; height: 1; color: Theme.divider }

                // Visible toggle
                RowLayout {
                    Layout.fillWidth: true
                    Text {
                        text: "Visible"
                        font.pixelSize: Theme.fsBody
                        color: Theme.textPrimary
                        Layout.fillWidth: true
                    }
                    Chip {
                        text: root.showCpArrows ? "On" : "Off"
                        checked: root.showCpArrows
                        onClicked: root.showCpArrows = !root.showCpArrows
                    }
                }

                Rectangle { Layout.fillWidth: true; height: 1; color: Theme.divider }

                // Line width
                SectionHeader { text: "Line Width" }
                Flow {
                    Layout.fillWidth: true
                    spacing: Theme.s2
                    Chip { text: "Thin";   checked: root.arrowLineWidth === 1.0; onClicked: root.arrowLineWidth = 1.0 }
                    Chip { text: "Medium"; checked: root.arrowLineWidth === 1.5; onClicked: root.arrowLineWidth = 1.5 }
                    Chip { text: "Thick";  checked: root.arrowLineWidth === 2.5; onClicked: root.arrowLineWidth = 2.5 }
                }

                // Arrow scale
                SectionHeader { text: "Arrow Length" }
                Flow {
                    Layout.fillWidth: true
                    spacing: Theme.s2
                    Chip { text: "0.5x"; checked: root.arrowScale === 0.5; onClicked: root.arrowScale = 0.5 }
                    Chip { text: "1x";   checked: root.arrowScale === 1.0; onClicked: root.arrowScale = 1.0 }
                    Chip { text: "2x";   checked: root.arrowScale === 2.0; onClicked: root.arrowScale = 2.0 }
                    Chip { text: "3x";   checked: root.arrowScale === 3.0; onClicked: root.arrowScale = 3.0 }
                }

                // Colormap
                SectionHeader { text: "Colormap" }
                Flow {
                    Layout.fillWidth: true
                    spacing: Theme.s2
                    Chip { text: "Coolwarm"; checked: root.arrowColormap === "coolwarm"; onClicked: root.arrowColormap = "coolwarm" }
                    Chip { text: "Viridis";  checked: root.arrowColormap === "viridis";  onClicked: root.arrowColormap = "viridis" }
                    Chip { text: "Mono";     checked: root.arrowColormap === "mono";     onClicked: root.arrowColormap = "mono" }
                }
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
