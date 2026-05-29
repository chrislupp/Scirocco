import QtQuick
import QtCharts

Item {
    ChartView {
        id: chart
        anchors.fill: parent
        antialiasing: true
        legend.visible: true

        ValueAxis {
            id: xAxisConv
            titleText: "Iteration"
            min: 0
            max: Math.max(reader.numIterations - 1, 1)
        }

        ValueAxis {
            id: yAxisCl
            titleText: "Cl, Cm"
        }

        ValueAxis {
            id: yAxisCd
            titleText: "Cd"
        }

        LineSeries {
            id: clLine
            axisX: xAxisConv
            axisY: yAxisCl
            name: "Cl"
            color: "blue"
            width: 2
        }

        LineSeries {
            id: cmLine
            axisX: xAxisConv
            axisY: yAxisCl
            name: "Cm"
            color: "green"
            width: 2
        }

        LineSeries {
            id: cdLine
            axisX: xAxisConv
            axisYRight: yAxisCd
            name: "Cd"
            color: "red"
            width: 2
        }

        ScatterSeries {
            id: currentMarkerCl
            axisX: xAxisConv
            axisY: yAxisCl
            markerSize: 10
            color: "blue"
            borderColor: "black"
            borderWidth: 1
        }

        ScatterSeries {
            id: currentMarkerCd
            axisX: xAxisConv
            axisYRight: yAxisCd
            markerSize: 10
            color: "red"
            borderColor: "black"
            borderWidth: 1
        }
    }

    Connections {
        target: reader
        function onFileChanged() {
            loadHistory()
        }
        function onIterationChanged() {
            updateMarker()
        }
    }

    function loadHistory() {
        var cl = reader.clHistory
        var cd = reader.cdHistory
        var cm = reader.cmHistory

        clLine.clear()
        cdLine.clear()
        cmLine.clear()

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
