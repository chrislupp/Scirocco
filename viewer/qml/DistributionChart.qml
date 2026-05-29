import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtCharts

Item {
    ColumnLayout {
        anchors.fill: parent
        spacing: 0

        ComboBox {
            id: distSelector
            Layout.leftMargin: 8
            Layout.topMargin: 4
            Layout.preferredWidth: 200
            model: ["cf", "ue", "ue_inv", "θ", "δ*", "n/Cτ", "Reθ", "Hk"]
            onCurrentTextChanged: updateDist()
        }

        ChartView {
            id: chart
            Layout.fillWidth: true
            Layout.fillHeight: true
            antialiasing: true
            legend.visible: false

            ValueAxis {
                id: xAxisDist
                titleText: "x/c"
                min: 0
                max: 1
            }

            ValueAxis {
                id: yAxisDist
                titleText: distSelector.currentText
            }

            LineSeries {
                id: distLine
                axisX: xAxisDist
                axisY: yAxisDist
                color: "darkgreen"
                width: 2
            }
        }
    }

    Connections {
        target: reader
        function onIterationChanged() {
            updateDist()
        }
    }

    function updateDist() {
        var pts = reader.distributionSeries(distSelector.currentText)
        distLine.clear()
        for (var i = 0; i < pts.length; i++)
            distLine.append(pts[i].x, pts[i].y)

        yAxisDist.titleText = distSelector.currentText

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
