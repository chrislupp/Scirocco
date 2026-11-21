import QtQuick
import SciroccoViewer

Item {
    id: root
    property bool showArrows: false
    property bool showPanels: false

    // Arrow display settings (controlled from Main.qml)
    property real arrowLineWidth: 1.5
    property real arrowScale: 1.0
    property string arrowColormap: "coolwarm"   // "coolwarm", "viridis", "mono"

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

    // Shown until a file is open.
    Column {
        id: emptyState
        anchors.centerIn: parent
        spacing: Theme.s4
        visible: reader.numIterations === 0

        Image {
            anchors.horizontalCenter: parent.horizontalCenter
            source: Theme.darkMode
                    ? "qrc:/qt/qml/SciroccoViewer/resources/logo-full-dark.svg"
                    : "qrc:/qt/qml/SciroccoViewer/resources/logo-full.svg"
            sourceSize: Qt.size(220, 151)
            opacity: 0.9
        }

        Text {
            anchors.horizontalCenter: parent.horizontalCenter
            text: "Open an HDF5 file to begin"
            font.pixelSize: Theme.fsControl
            color: Theme.textMuted
        }

        Text {
            anchors.horizontalCenter: parent.horizontalCenter
            text: (Qt.platform.os === "osx" ? "⌘O" : "Ctrl+O") + " or double-click"
            font.pixelSize: Theme.fsCaption
            color: Theme.textDim
        }
    }

    onShowArrowsChanged: canvas.requestPaint()
    onShowPanelsChanged: canvas.requestPaint()
    onArrowLineWidthChanged: canvas.requestPaint()
    onArrowScaleChanged: canvas.requestPaint()
    onArrowColormapChanged: canvas.requestPaint()
    onWidthChanged: canvas.requestPaint()
    onHeightChanged: canvas.requestPaint()

    Connections {
        target: reader
        function onIterationChanged() { computeBounds(); canvas.requestPaint() }
    }

    Connections {
        target: Theme
        function onDarkModeChanged() { canvas.requestPaint() }
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

    function cpColor(cp, cpMin, cpMax) {
        if (arrowColormap === "mono")
            return Theme.arrowStroke.toString()

        var range = cpMax - cpMin
        if (range < 1e-12) range = 1
        var t = (cp - cpMin) / range

        if (arrowColormap === "viridis") {
            // Simplified viridis: purple → teal → yellow
            var r, g, b
            if (t < 0.5) {
                var s = t * 2
                r = Math.round(68 + s * (33 - 68))
                g = Math.round(1 + s * (144 - 1))
                b = Math.round(84 + s * (140 - 84))
            } else {
                var s2 = (t - 0.5) * 2
                r = Math.round(33 + s2 * (253 - 33))
                g = Math.round(144 + s2 * (231 - 144))
                b = Math.round(140 + s2 * (37 - 140))
            }
            return "rgb(" + r + "," + g + "," + b + ")"
        }

        // coolwarm: blue → white → red
        var r2, g2, b2
        if (t < 0.5) {
            var s3 = t * 2
            r2 = Math.round(59 + s3 * (221 - 59))
            g2 = Math.round(76 + s3 * (221 - 76))
            b2 = Math.round(192 + s3 * (221 - 192))
        } else {
            var s4 = (t - 0.5) * 2
            r2 = Math.round(221 + s4 * (180 - 221))
            g2 = Math.round(221 + s4 * (4 - 221))
            b2 = Math.round(221 + s4 * (38 - 221))
        }
        return "rgb(" + r2 + "," + g2 + "," + b2 + ")"
    }

    function draw() {
        var ctx = canvas.getContext("2d")
        ctx.clearRect(0, 0, canvas.width, canvas.height)

        var foil = reader.airfoilSeries
        if (foil.length === 0)
            return  // emptyState shows instead

        computeBounds()

        // Draw arrows (behind airfoil)
        if (showArrows && reader.hasArrows) {
            var bases = reader.arrowBases
            var tips = reader.arrowTips
            var cpVals = reader.arrowCpValues
            var headLen = 5 + root.arrowLineWidth * 2

            // Find Cp range for colormap
            var cpMin = Infinity, cpMax = -Infinity
            for (var c = 0; c < cpVals.length; c++) {
                cpMin = Math.min(cpMin, cpVals[c])
                cpMax = Math.max(cpMax, cpVals[c])
            }

            for (var a = 0; a < bases.length; a++) {
                var bx = mapX(bases[a].x)
                var by = mapY(bases[a].y)
                // Apply arrow scale: interpolate tip between base and full tip
                var fullTx = mapX(tips[a].x)
                var fullTy = mapY(tips[a].y)
                var tx = bx + (fullTx - bx) * root.arrowScale
                var ty = by + (fullTy - by) * root.arrowScale

                var dx = tx - bx
                var dy = ty - by
                var len = Math.sqrt(dx * dx + dy * dy)
                if (len < 2) continue

                var color = cpColor(cpVals[a], cpMin, cpMax)
                ctx.strokeStyle = color
                ctx.fillStyle = color
                ctx.lineWidth = root.arrowLineWidth

                // Arrow shaft
                ctx.beginPath()
                ctx.moveTo(bx, by)
                ctx.lineTo(tx, ty)
                ctx.stroke()

                // Arrowhead direction: positive Cp points toward surface, negative away
                var ux = dx / len
                var uy = dy / len
                var hx, hy, hdx, hdy
                if (cpVals[a] > 0) {
                    // Head at base end (pointing toward airfoil)
                    hx = bx; hy = by
                    hdx = ux; hdy = uy
                } else {
                    // Head at tip end (pointing away from airfoil)
                    hx = tx; hy = ty
                    hdx = -ux; hdy = -uy
                }
                var px1 = hx + headLen * hdx + headLen * 0.4 * hdy
                var py1 = hy + headLen * hdy - headLen * 0.4 * hdx
                var px2 = hx + headLen * hdx - headLen * 0.4 * hdy
                var py2 = hy + headLen * hdy + headLen * 0.4 * hdx

                ctx.beginPath()
                ctx.moveTo(hx, hy)
                ctx.lineTo(px1, py1)
                ctx.lineTo(px2, py2)
                ctx.closePath()
                ctx.fill()
            }
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
        var dimmed = showArrows || showPanels
        ctx.strokeStyle = dimmed
            ? Qt.rgba(Theme.foilStroke.r, Theme.foilStroke.g, Theme.foilStroke.b, 0.5).toString()
            : Theme.foilStroke.toString()
        ctx.lineWidth = dimmed ? 1.5 : 2
        ctx.beginPath()
        ctx.moveTo(mapX(foil[0].x), mapY(foil[0].y))
        for (var i = 1; i < foil.length; i++)
            ctx.lineTo(mapX(foil[i].x), mapY(foil[i].y))
        ctx.stroke()

        // Draw panel node dots
        if (showPanels) {
            ctx.fillStyle = Theme.accent.toString()
            var r = 2.5
            for (var p = 0; p < foil.length; p++) {
                var ppx = mapX(foil[p].x)
                var ppy = mapY(foil[p].y)
                ctx.beginPath()
                ctx.arc(ppx, ppy, r, 0, 2 * Math.PI)
                ctx.fill()
            }
        }
    }
}
