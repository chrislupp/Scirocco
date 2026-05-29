import QtQuick
import QtQuick.Layouts
import SciroccoViewer

Panel {
    id: root
    width: Theme.panelWidth
    signal closeRequested()

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: Theme.s4
        spacing: Theme.s3

        // Header
        RowLayout {
            Layout.fillWidth: true
            Text {
                text: "Results"
                font.pixelSize: Theme.fsTitle
                font.weight: Font.DemiBold
                color: Theme.textPrimary
                Layout.fillWidth: true
            }
            Rectangle {
                width: 28; height: 28
                radius: Theme.rControl
                color: closeMa.containsMouse ? Theme.hover : "transparent"
                Text {
                    anchors.centerIn: parent
                    text: "✕"
                    font.pixelSize: 12
                    color: Theme.textMuted
                }
                MouseArea {
                    id: closeMa
                    anchors.fill: parent
                    hoverEnabled: true
                    cursorShape: Qt.PointingHandCursor
                    onClicked: root.closeRequested()
                }
            }
        }

        // Divider
        Rectangle { Layout.fillWidth: true; height: 1; color: Theme.divider }

        // Scrollable content
        Flickable {
            Layout.fillWidth: true
            Layout.fillHeight: true
            contentHeight: contentCol.implicitHeight
            clip: true
            boundsBehavior: Flickable.StopAtBounds

            ColumnLayout {
                id: contentCol
                width: parent.width
                spacing: Theme.s4

                SectionHeader { text: "Operating Conditions" }

                GridLayout {
                    columns: 2
                    columnSpacing: Theme.s3
                    rowSpacing: Theme.s2
                    Layout.fillWidth: true

                    DataLabel { text: "Name" }
                    DataValue { text: reader.data.name || "—" }
                    DataLabel { text: "α" }
                    DataValue { text: reader.numIterations > 0 ? reader.data.alpha.toFixed(2) + "°" : "—" }
                    DataLabel { text: "Re" }
                    DataValue { text: reader.numIterations > 0 ? reader.data.Re.toExponential(2) : "—" }
                    DataLabel { text: "Ma" }
                    DataValue { text: reader.numIterations > 0 ? reader.data.Ma.toFixed(3) : "—" }
                }

                SectionHeader { text: "Coefficients" }

                GridLayout {
                    columns: 2
                    columnSpacing: Theme.s3
                    rowSpacing: Theme.s2
                    Layout.fillWidth: true

                    DataLabel { text: "Cl" }
                    DataValue { text: reader.numIterations > 0 ? reader.data.cl.toFixed(4) : "—" }
                    DataLabel { text: "Cd" }
                    DataValue { text: reader.numIterations > 0 ? reader.data.cd.toFixed(5) : "—" }
                    DataLabel { text: "Cdf" }
                    DataValue { text: reader.numIterations > 0 ? reader.data.cdf.toFixed(5) : "—" }
                    DataLabel { text: "Cdp" }
                    DataValue { text: reader.numIterations > 0 ? reader.data.cdp.toFixed(5) : "—" }
                    DataLabel { text: "Cm" }
                    DataValue { text: reader.numIterations > 0 ? reader.data.cm.toFixed(4) : "—" }
                }

                SectionHeader { text: "Status" }

                Row {
                    spacing: Theme.s2
                    Rectangle {
                        width: 6; height: 6; radius: 3
                        anchors.verticalCenter: parent.verticalCenter
                        color: reader.numIterations > 0 && reader.data.converged
                               ? Theme.accent : Qt.rgba(1, 0.3, 0.3, 1)
                    }
                    Text {
                        text: reader.numIterations > 0
                              ? (reader.data.converged ? "Converged" : "Not converged")
                              : "—"
                        font.pixelSize: Theme.fsBody
                        color: Theme.textPrimary
                    }
                }
            }
        }
    }

    component DataLabel: Text {
        font.pixelSize: Theme.fsBody
        color: Theme.textMuted
    }

    component DataValue: Text {
        font.family: Theme.fontMono
        font.pixelSize: Theme.fsBody
        color: Theme.textPrimary
    }
}
