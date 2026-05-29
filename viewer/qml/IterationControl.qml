import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

ToolBar {
    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: 8
        anchors.rightMargin: 8
        spacing: 8

        ToolButton {
            text: "◀"
            enabled: reader.currentIteration > 0
            onClicked: reader.setCurrentIteration(reader.currentIteration - 1)
        }

        Slider {
            id: iterSlider
            Layout.fillWidth: true
            from: 0
            to: Math.max(reader.numIterations - 1, 0)
            stepSize: 1
            value: reader.currentIteration
            enabled: reader.numIterations > 1
            onMoved: reader.setCurrentIteration(Math.round(value))
        }

        ToolButton {
            text: "▶"
            enabled: reader.currentIteration < reader.numIterations - 1
            onClicked: reader.setCurrentIteration(reader.currentIteration + 1)
        }

        Label {
            text: reader.numIterations > 0
                  ? "Iteration " + reader.currentIteration + " / " + (reader.numIterations - 1)
                  : "No file loaded"
            Layout.preferredWidth: 140
        }
    }
}
