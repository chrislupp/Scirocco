import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtCharts
import SciroccoViewer

Item {
    ColumnLayout {
        anchors.fill: parent
        spacing: 0

        // Chip row
        Flow {
            Layout.fillWidth: true
            Layout.leftMargin: Theme.s3
            Layout.topMargin: Theme.s2
            Layout.bottomMargin: Theme.s2
            spacing: Theme.s2

            Chip {
                text: "cf"; checked: true
                ButtonGroup.group: distGroup
            }
            Chip { text: "ue"; ButtonGroup.group: distGroup }
            Chip { text: "θ"; ButtonGroup.group: distGroup }
            Chip { text: "δ*"; ButtonGroup.group: distGroup }
            Chip { text: "Hk"; ButtonGroup.group: distGroup }
            Chip { text: "Reθ"; ButtonGroup.group: distGroup }
            Chip { text: "n/Cτ"; ButtonGroup.group: distGroup }
        }

        ButtonGroup {
            id: distGroup
            exclusive: true
            onCheckedButtonChanged: updateDist()
        }

        ChartView {
            id: chart
            Layout.fillWidth: true
            Layout.fillHeight: true
            antialiasing: true
            backgroundColor: "transparent"
            plotAreaColor: "transparent"
            legend.visible: false
            margins { left: 4; right: 4; top: 4; bottom: 4 }

            ValueAxis {
                id: xAxisDist
                min: 0; max: 1
                color: Theme.divider
                gridLineColor: Theme.gridLine
                labelsColor: Theme.axisLabel
                labelsFont.family: Theme.fontMono
                labelsFont.pixelSize: Theme.fsCaption
                titleText: "x/c"
                titleFont.pixelSize: Theme.fsCaption
                titleBrush: Theme.textDim
            }

            ValueAxis {
                id: yAxisDist
                color: Theme.divider
                gridLineColor: Theme.gridLine
                labelsColor: Theme.axisLabel
                labelsFont.family: Theme.fontMono
                labelsFont.pixelSize: Theme.fsCaption
                titleVisible: false
            }

            LineSeries {
                id: distLine
                axisX: xAxisDist; axisY: yAxisDist
                color: Theme.accent
                width: 1.5
            }
        }
    }

    Connections {
        target: reader
        function onIterationChanged() { updateDist() }
    }

    function updateDist() {
        var name = distGroup.checkedButton ? distGroup.checkedButton.text : "cf"
        var pts = reader.distributionSeries(name)
        distLine.clear()
        for (var i = 0; i < pts.length; i++)
            distLine.append(pts[i].x, pts[i].y)

        if (pts.length > 0) {
            var ymin = Infinity, ymax = -Infinity
            for (var k = 0; k < pts.length; k++) {
                ymin = Math.min(ymin, pts[k].y)
                ymax = Math.max(ymax, pts[k].y)
            }
            var margin = 0.1 * Math.max(ymax - ymin, 1e-10)
            yAxisDist.min = ymin - margin
            yAxisDist.max = ymax + margin
        }
    }
}
