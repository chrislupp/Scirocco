import QtQuick
import QtQuick.Controls
import QtCharts

Item {
    property bool showArrows: false

    ChartView {
        id: chart
        anchors.fill: parent
        antialiasing: true
        legend.visible: false

        ValueAxis {
            id: xAxis
            titleText: "x"
            gridVisible: false
        }

        ValueAxis {
            id: yAxis
            titleText: "y"
            gridVisible: false
        }

        LineSeries {
            id: airfoilLine
            axisX: xAxis
            axisY: yAxis
            color: "black"
            width: 2
        }

        LineSeries {
            id: wakeLine
            axisX: xAxis
            axisY: yAxis
            color: "gray"
            width: 1
            style: Qt.DashLine
        }
    }

    Canvas {
        id: arrowCanvas
        anchors.fill: parent
        visible: showArrows && reader.hasArrows
        onPaint: drawArrows()
    }

    onShowArrowsChanged: arrowCanvas.requestPaint()

    Connections {
        target: reader
        function onIterationChanged() {
            updateChart()
        }
    }

    function updateChart() {
        var foil = reader.airfoilSeries
        airfoilLine.clear()
        for (var i = 0; i < foil.length; i++)
            airfoilLine.append(foil[i].x, foil[i].y)

        var wake = reader.wakeSeries
        wakeLine.clear()
        for (var j = 0; j < wake.length; j++)
            wakeLine.append(wake[j].x, wake[j].y)

        // Equal aspect ratio
        if (foil.length > 0) {
            var xmin = Infinity, xmax = -Infinity
            var ymin = Infinity, ymax = -Infinity
            for (var fi = 0; fi < foil.length; fi++) {
                xmin = Math.min(xmin, foil[fi].x)
                xmax = Math.max(xmax, foil[fi].x)
                ymin = Math.min(ymin, foil[fi].y)
                ymax = Math.max(ymax, foil[fi].y)
            }
            // Include arrows in range
            if (showArrows) {
                var tips = reader.arrowTips
                for (var ai = 0; ai < tips.length; ai++) {
                    xmin = Math.min(xmin, tips[ai].x)
                    xmax = Math.max(xmax, tips[ai].x)
                    ymin = Math.min(ymin, tips[ai].y)
                    ymax = Math.max(ymax, tips[ai].y)
                }
            }
            var xRange = xmax - xmin
            var yRange = ymax - ymin
            var margin = 0.1 * Math.max(xRange, yRange)
            var maxRange = Math.max(xRange, yRange) + 2 * margin
            var xCenter = (xmin + xmax) / 2
            var yCenter = (ymin + ymax) / 2
            xAxis.min = xCenter - maxRange / 2
            xAxis.max = xCenter + maxRange / 2
            yAxis.min = yCenter - maxRange / 2
            yAxis.max = yCenter + maxRange / 2
        }

        arrowCanvas.requestPaint()
    }

    function drawArrows() {
        var ctx = arrowCanvas.getContext("2d")
        ctx.clearRect(0, 0, arrowCanvas.width, arrowCanvas.height)
        if (!showArrows || !reader.hasArrows) return

        var bases = reader.arrowBases
        var tips = reader.arrowTips

        ctx.strokeStyle = "rgba(0, 100, 200, 0.6)"
        ctx.lineWidth = 1
        ctx.beginPath()
        for (var i = 0; i < bases.length; i++) {
            var p0 = chart.mapToPosition(Qt.point(bases[i].x, bases[i].y), airfoilLine)
            var p1 = chart.mapToPosition(Qt.point(tips[i].x, tips[i].y), airfoilLine)
            ctx.moveTo(p0.x, p0.y)
            ctx.lineTo(p1.x, p1.y)
        }
        ctx.stroke()
    }
}
