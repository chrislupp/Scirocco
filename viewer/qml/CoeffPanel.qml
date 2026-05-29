import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

GroupBox {
    title: "Results"

    GridLayout {
        columns: 2
        columnSpacing: 12
        rowSpacing: 4
        anchors.fill: parent

        Label { text: "Name:"; font.bold: true }
        Label { text: reader.data.name || "—" }

        Label { text: "α:"; font.bold: true }
        Label { text: reader.numIterations > 0 ? reader.data.alpha.toFixed(2) + "°" : "—" }

        Label { text: "Re:"; font.bold: true }
        Label { text: reader.numIterations > 0 ? reader.data.Re.toExponential(2) : "—" }

        Label { text: "Ma:"; font.bold: true }
        Label { text: reader.numIterations > 0 ? reader.data.Ma.toFixed(3) : "—" }

        Rectangle { Layout.columnSpan: 2; Layout.fillWidth: true; height: 1; color: "#ccc" }

        Label { text: "Cl:"; font.bold: true }
        Label { text: reader.numIterations > 0 ? reader.data.cl.toFixed(4) : "—" }

        Label { text: "Cd:"; font.bold: true }
        Label { text: reader.numIterations > 0 ? reader.data.cd.toFixed(5) : "—" }

        Label { text: "Cdf:"; font.bold: true }
        Label { text: reader.numIterations > 0 ? reader.data.cdf.toFixed(5) : "—" }

        Label { text: "Cdp:"; font.bold: true }
        Label { text: reader.numIterations > 0 ? reader.data.cdp.toFixed(5) : "—" }

        Label { text: "Cm:"; font.bold: true }
        Label { text: reader.numIterations > 0 ? reader.data.cm.toFixed(4) : "—" }

        Rectangle { Layout.columnSpan: 2; Layout.fillWidth: true; height: 1; color: "#ccc" }

        Label { text: "Conv:"; font.bold: true }
        Label {
            text: reader.numIterations > 0 ? (reader.data.converged ? "Yes" : "No") : "—"
            color: reader.data.converged ? "green" : "red"
        }
    }
}
