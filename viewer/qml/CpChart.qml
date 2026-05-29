import QtQuick
import QtCharts

Item {
    ChartView {
        id: chart
        anchors.fill: parent
        antialiasing: true
        legend.visible: true

        ValueAxis {
            id: xAxisCp
            titleText: "x/c"
            min: 0
            max: 1
        }

        ValueAxis {
            id: yAxisCp
            titleText: "Cp"
            reverse: true
        }

        LineSeries {
            id: cpLine
            axisX: xAxisCp
            axisY: yAxisCp
            name: "Cp"
            color: "blue"
            width: 2
        }

        LineSeries {
            id: cpiLine
            axisX: xAxisCp
            axisY: yAxisCp
            name: "Cp inv"
            color: "blue"
            width: 1
            style: Qt.DashLine
        }
    }

    Connections {
        target: reader
        function onIterationChanged() {
            updateCp()
        }
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

        // Auto-range y axis
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
