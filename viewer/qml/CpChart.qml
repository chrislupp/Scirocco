import QtQuick
import QtCharts
import SciroccoViewer

Item {
    ChartView {
        id: chart
        anchors.fill: parent
        antialiasing: true
        backgroundColor: "transparent"
        plotAreaColor: "transparent"
        legend.visible: false
        margins { left: 4; right: 4; top: 4; bottom: 4 }

        ValueAxis {
            id: xAxisCp
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
            id: yAxisCp
            reverse: true
            color: Theme.divider
            gridLineColor: Theme.gridLine
            labelsColor: Theme.axisLabel
            labelsFont.family: Theme.fontMono
            labelsFont.pixelSize: Theme.fsCaption
            titleText: "Cp"
            titleFont.pixelSize: Theme.fsCaption
            titleBrush: Theme.textDim
        }

        LineSeries {
            id: cpLine
            axisX: xAxisCp; axisY: yAxisCp
            color: Theme.cpLine
            width: 1.5
        }

        LineSeries {
            id: cpiLine
            axisX: xAxisCp; axisY: yAxisCp
            color: Theme.cpiLine
            width: 1
            style: Qt.DashLine
        }
    }

    Connections {
        target: reader
        function onIterationChanged() { updateCp() }
    }

    function updateCp() {
        var cp = reader.cpSeries
        cpLine.clear()
        for (var i = 0; i < cp.length; i++)
            cpLine.append(cp[i].x, cp[i].y)

        var cpi = reader.cpiSeries
        cpiLine.clear()
        for (var j = 0; j < cpi.length; j++)
            cpiLine.append(cpi[j].x, cpi[j].y)

        if (cp.length > 0) {
            var ymin = Infinity, ymax = -Infinity
            for (var k = 0; k < cp.length; k++) {
                ymin = Math.min(ymin, cp[k].y)
                ymax = Math.max(ymax, cp[k].y)
            }
            var margin = 0.1 * (ymax - ymin)
            yAxisCp.min = ymin - margin
            yAxisCp.max = ymax + margin
        }
    }
}
