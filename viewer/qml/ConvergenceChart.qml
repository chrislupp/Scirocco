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
        legend.visible: true
        legend.labelColor: Theme.textMuted
        legend.font.pixelSize: Theme.fsCaption
        margins { left: 4; right: 4; top: 4; bottom: 4 }

        ValueAxis {
            id: xAxisConv
            min: 0
            max: Math.max(reader.numIterations - 1, 1)
            color: Theme.divider
            gridLineColor: Theme.gridLine
            labelsColor: Theme.axisLabel
            labelsFont.family: Theme.fontMono
            labelsFont.pixelSize: Theme.fsCaption
            titleText: "Iteration"
            titleFont.pixelSize: Theme.fsCaption
            titleBrush: Theme.textDim
        }

        ValueAxis {
            id: yAxisCl
            color: Theme.divider
            gridLineColor: Theme.gridLine
            labelsColor: Theme.axisLabel
            labelsFont.family: Theme.fontMono
            labelsFont.pixelSize: Theme.fsCaption
            titleVisible: false
        }

        ValueAxis {
            id: yAxisCd
            color: Theme.divider
            gridVisible: false
            labelsColor: Theme.axisLabel
            labelsFont.family: Theme.fontMono
            labelsFont.pixelSize: Theme.fsCaption
            titleVisible: false
        }

        LineSeries {
            id: clLine
            axisX: xAxisConv; axisY: yAxisCl
            name: "Cl"
            color: Theme.cpLine
            width: 1.5
        }

        LineSeries {
            id: cmLine
            axisX: xAxisConv; axisY: yAxisCl
            name: "Cm"
            color: "#5be0a8"
            width: 1.5
        }

        LineSeries {
            id: cdLine
            axisX: xAxisConv; axisYRight: yAxisCd
            name: "Cd"
            color: "#ff7eb6"
            width: 1.5
        }

        ScatterSeries {
            id: currentMarkerCl
            axisX: xAxisConv; axisY: yAxisCl
            markerSize: 8
            color: Theme.cpLine
            borderColor: "transparent"
        }

        ScatterSeries {
            id: currentMarkerCd
            axisX: xAxisConv; axisYRight: yAxisCd
            markerSize: 8
            color: "#ff7eb6"
            borderColor: "transparent"
        }
    }

    Connections {
        target: reader
        function onFileChanged() { loadHistory() }
        function onIterationChanged() { updateMarker() }
    }

    function loadHistory() {
        var cl = reader.clHistory
        var cd = reader.cdHistory
        var cm = reader.cmHistory

        clLine.clear(); cdLine.clear(); cmLine.clear()

        var clMin = Infinity, clMax = -Infinity
        var cdMin = Infinity, cdMax = -Infinity

        for (var i = 0; i < cl.length; i++) {
            clLine.append(cl[i].x, cl[i].y)
            clMin = Math.min(clMin, cl[i].y)
            clMax = Math.max(clMax, cl[i].y)
        }
        for (var j = 0; j < cd.length; j++) {
            cdLine.append(cd[j].x, cd[j].y)
            cdMin = Math.min(cdMin, cd[j].y)
            cdMax = Math.max(cdMax, cd[j].y)
        }
        for (var k = 0; k < cm.length; k++) {
            cmLine.append(cm[k].x, cm[k].y)
            clMin = Math.min(clMin, cm[k].y)
            clMax = Math.max(clMax, cm[k].y)
        }

        if (cl.length > 0) {
            var mCl = 0.1 * (clMax - clMin)
            yAxisCl.min = clMin - mCl
            yAxisCl.max = clMax + mCl
        }
        if (cd.length > 0) {
            var mCd = 0.1 * (cdMax - cdMin)
            yAxisCd.min = cdMin - mCd
            yAxisCd.max = cdMax + mCd
        }

        xAxisConv.max = Math.max(cl.length - 1, 1)
        updateMarker()
    }

    function updateMarker() {
        currentMarkerCl.clear()
        currentMarkerCd.clear()
        var idx = reader.currentIteration
        var cl = reader.clHistory
        var cd = reader.cdHistory
        if (idx < cl.length)
            currentMarkerCl.append(idx, cl[idx].y)
        if (idx < cd.length)
            currentMarkerCd.append(idx, cd[idx].y)
    }
}
