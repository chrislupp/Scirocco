import QtQuick
import SciroccoViewer

Item {
    id: root
    property bool showArrows: false

    property real dataXMin: 0
    property real dataXMax: 1
    property real dataYMin: -0.2
    property real dataYMax: 0.2
    property real pxPerUnit: 1

    Canvas {
        id: canvas
        anchors.fill: parent
        antialiasing: true
        onPaint: draw()
    }

    onShowArrowsChanged: canvas.requestPaint()
    onWidthChanged: canvas.requestPaint()
    onHeightChanged: canvas.requestPaint()

    Connections {
        target: reader
        function onIterationChanged() { computeBounds(); canvas.requestPaint() }
    }

    function computeBounds() {
        var foil = reader.airfoilSeries
        if (foil.length === 0) return

        var xmin = Infinity, xmax = -Infinity
        var ymin = Infinity, ymax = -Infinity
        for (var i = 0; i < foil.length; i++) {
            xmin = Math.min(xmin, foil[i].x)
            xmax = Math.max(xmax, foil[i].x)
            ymin = Math.min(ymin, foil[i].y)
            ymax = Math.max(ymax, foil[i].y)
        }
        if (showArrows && reader.hasArrows) {
            var tips = reader.arrowTips
            for (var j = 0; j < tips.length; j++) {
                xmin = Math.min(xmin, tips[j].x)
                xmax = Math.max(xmax, tips[j].x)
                ymin = Math.min(ymin, tips[j].y)
                ymax = Math.max(ymax, tips[j].y)
            }
        }
        dataXMin = xmin; dataXMax = xmax
        dataYMin = ymin; dataYMax = ymax
    }

    function mapX(x) {
        var w = canvas.width
        var h = canvas.height
        var margin = 0.12
        var xRange = dataXMax - dataXMin
        var yRange = dataYMax - dataYMin
        if (xRange < 1e-12) xRange = 1
        if (yRange < 1e-12) yRange = xRange * 0.3

        var scaleX = w * (1 - 2 * margin) / xRange
        var scaleY = h * (1 - 2 * margin) / yRange
        pxPerUnit = Math.min(scaleX, scaleY)

        var cx = (dataXMin + dataXMax) / 2
        return w / 2 + (x - cx) * pxPerUnit
    }

    function mapY(y) {
        var h = canvas.height
        var cy = (dataYMin + dataYMax) / 2
        return h / 2 - (y - cy) * pxPerUnit
    }

    function draw() {
        var ctx = canvas.getContext("2d")
        ctx.clearRect(0, 0, canvas.width, canvas.height)

        var foil = reader.airfoilSeries
        if (foil.length === 0) {
            // Empty state
            ctx.fillStyle = Theme.textDim.toString()
            ctx.font = Theme.fsControl + "px sans-serif"
            ctx.textAlign = "center"
            ctx.fillText("Open an HDF5 file to begin", canvas.width / 2, canvas.height / 2)
            return
        }

        computeBounds()

        // Draw arrows (behind airfoil)
        if (showArrows && reader.hasArrows) {
            var bases = reader.arrowBases
            var tips = reader.arrowTips
            ctx.strokeStyle = Theme.arrowStroke.toString()
            ctx.lineWidth = 1
            ctx.beginPath()
            for (var a = 0; a < bases.length; a++) {
                ctx.moveTo(mapX(bases[a].x), mapY(bases[a].y))
                ctx.lineTo(mapX(tips[a].x), mapY(tips[a].y))
            }
            ctx.stroke()
        }

        // Draw wake
        var wake = reader.wakeSeries
        if (wake.length > 1) {
            ctx.strokeStyle = Theme.wakeStroke.toString()
            ctx.lineWidth = 1
            ctx.setLineDash([6, 4])
            ctx.beginPath()
            ctx.moveTo(mapX(wake[0].x), mapY(wake[0].y))
            for (var w = 1; w < wake.length; w++)
                ctx.lineTo(mapX(wake[w].x), mapY(wake[w].y))
            ctx.stroke()
            ctx.setLineDash([])
        }

        // Draw airfoil contour
        ctx.strokeStyle = showArrows
            ? Qt.rgba(0.906, 0.925, 0.949, 0.5).toString()
            : Theme.foilStroke.toString()
        ctx.lineWidth = showArrows ? 1.5 : 2
        ctx.beginPath()
        ctx.moveTo(mapX(foil[0].x), mapY(foil[0].y))
        for (var i = 1; i < foil.length; i++)
            ctx.lineTo(mapX(foil[i].x), mapY(foil[i].y))
        ctx.stroke()
    }
}
