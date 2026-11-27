import QtQuick
import QtQuick.Shapes
import SciroccoViewer

// A line icon drawn from SVG path data on a 24 x 24 grid.
//
// Each entry in `paths` has an optional stroke path `s` (drawn with round
// caps and joins) and an optional fill path `f`. Both take `color`, so an
// icon follows the theme and its state without separate image files.
Item {
    id: root

    property string name
    property color color: Theme.textMuted
    property real size: 20
    // Stroke width in pixels at the drawn size.
    property real strokeWidth: 1.6

    implicitWidth: size
    implicitHeight: size

    readonly property var paths: ({
        // Cambered airfoil at incidence, as in the Scirocco logo.
        "airfoil": {
            s: "M21.50 15.46C21.44 15.40 21.37 15.34 21.11 15.13"
               + "C20.86 14.92 20.46 14.59 19.96 14.21"
               + "C19.46 13.84 18.82 13.36 18.11 12.89"
               + "C17.39 12.42 16.54 11.88 15.66 11.40"
               + "C14.79 10.92 13.81 10.41 12.83 10.01"
               + "C11.86 9.60 10.81 9.20 9.82 8.95C8.84 8.71 7.78 8.56 6.90 8.54"
               + "C6.02 8.52 5.17 8.64 4.52 8.82C3.88 9.00 3.35 9.31 3.01 9.63"
               + "C2.67 9.94 2.51 10.36 2.50 10.72C2.49 11.09 2.68 11.49 2.98 11.82"
               + "C3.28 12.16 3.75 12.46 4.30 12.73C4.86 12.99 5.55 13.21 6.30 13.40"
               + "C7.04 13.60 7.88 13.74 8.77 13.90"
               + "C9.66 14.05 10.66 14.20 11.63 14.33"
               + "C12.60 14.46 13.64 14.57 14.60 14.68"
               + "C15.55 14.80 16.51 14.90 17.33 14.99"
               + "C18.16 15.09 18.94 15.17 19.55 15.24"
               + "C20.16 15.31 20.67 15.37 21.00 15.40"
               + "C21.32 15.44 21.42 15.45 21.50 15.46Z"
        },
        // Plot axes with a curve.
        "chart": {
            s: "M3.5 3.5V20.5H20.5"
               + "M7 16C9.2 16 9.6 8 12.4 8C14.6 8 15 12.5 17 12.5C18.2 12.5 19 11 20 9.5"
        },
        // Label and value rows in a panel.
        "results": {
            s: "M5.5 3.5H18.5A2 2 0 0 1 20.5 5.5V18.5A2 2 0 0 1 18.5 20.5H5.5"
               + "A2 2 0 0 1 3.5 18.5V5.5A2 2 0 0 1 5.5 3.5Z"
               + "M7.5 8.5H10.5M13.5 8.5H16.5M7.5 12H10.5M13.5 12H16.5"
               + "M7.5 15.5H10.5M13.5 15.5H16.5"
        },
        // Readout lines inside viewfinder corners.
        "hud": {
            s: "M3.5 8V5.5A2 2 0 0 1 5.5 3.5H8M16 3.5H18.5A2 2 0 0 1 20.5 5.5V8"
               + "M20.5 16V18.5A2 2 0 0 1 18.5 20.5H16M8 20.5H5.5A2 2 0 0 1 3.5 18.5V16"
               + "M10 9H16.5M12 12H16.5M8.5 15H16.5"
        },
        // An airfoil as six straight panels between nodes.
        "panels": {
            s: "M21.00 12.68L13.96 9.01L6.15 8.41L3.00 12.68L7.05 15.59"
               + "L13.64 14.74Z",
            f: "M22.45 12.68A1.45 1.45 0 1 1 19.55 12.68"
               + "A1.45 1.45 0 1 1 22.45 12.68ZM15.41 9.01"
               + "A1.45 1.45 0 1 1 12.51 9.01A1.45 1.45 0 1 1 15.41 9.01ZM7.60 8.41"
               + "A1.45 1.45 0 1 1 4.70 8.41A1.45 1.45 0 1 1 7.60 8.41ZM4.45 12.68"
               + "A1.45 1.45 0 1 1 1.55 12.68A1.45 1.45 0 1 1 4.45 12.68ZM8.50 15.59"
               + "A1.45 1.45 0 1 1 5.60 15.59A1.45 1.45 0 1 1 8.50 15.59Z"
               + "M15.09 14.74A1.45 1.45 0 1 1 12.19 14.74"
               + "A1.45 1.45 0 1 1 15.09 14.74Z"
        },
        // Pressure arrows off the upper surface of an airfoil.
        "arrows": {
            s: "M21.00 18.16C20.93 18.13 20.85 18.10 20.57 17.99"
               + "C20.30 17.88 19.86 17.71 19.33 17.52"
               + "C18.79 17.34 18.11 17.10 17.37 16.89"
               + "C16.62 16.67 15.74 16.43 14.86 16.23"
               + "C13.97 16.04 12.99 15.84 12.04 15.72"
               + "C11.09 15.60 10.08 15.49 9.17 15.49"
               + "C8.25 15.49 7.31 15.57 6.53 15.72C5.75 15.86 5.03 16.10 4.49 16.36"
               + "C3.95 16.61 3.53 16.94 3.28 17.24C3.03 17.54 2.95 17.88 3.00 18.16"
               + "C3.05 18.44 3.27 18.72 3.60 18.92C3.92 19.12 4.40 19.26 4.95 19.36"
               + "C5.50 19.46 6.17 19.50 6.89 19.51C7.61 19.52 8.43 19.47 9.27 19.42"
               + "C10.12 19.37 11.06 19.29 11.96 19.21"
               + "C12.87 19.12 13.83 19.01 14.71 18.91"
               + "C15.58 18.82 16.46 18.71 17.21 18.62"
               + "C17.97 18.53 18.68 18.45 19.23 18.38"
               + "C19.79 18.31 20.25 18.26 20.54 18.22"
               + "C20.84 18.18 20.92 18.17 21.00 18.16ZM6 14V3.5M4.3 5.2L6 3.5"
               + "L7.7 5.2M12 13.9V5.5M10.3 7.2L12 5.5L13.7 7.2M18 15.3V9M16.3 10.7"
               + "L18 9L19.7 10.7"
        },
        // Film frame with a play mark.
        "video": {
            s: "M5.5 4.5H18.5A2 2 0 0 1 20.5 6.5V17.5A2 2 0 0 1 18.5 19.5H5.5"
               + "A2 2 0 0 1 3.5 17.5V6.5A2 2 0 0 1 5.5 4.5Z"
               + "M7.5 4.5V19.5M16.5 4.5V19.5M3.5 9H7.5M3.5 15H7.5M16.5 9H20.5M16.5 15H20.5",
            f: "M10.5 9.6V14.4C10.5 14.9 11 15.1 11.4 14.9L14.9 12.5"
               + "C15.2 12.3 15.2 11.7 14.9 11.5L11.4 9.1C11 8.9 10.5 9.1 10.5 9.6Z"
        },
        "sun": {
            s: "M16 12A4 4 0 1 1 8 12A4 4 0 1 1 16 12Z"
               + "M12 2.5V4.5M12 19.5V21.5M2.5 12H4.5M19.5 12H21.5"
               + "M5.3 5.3L6.7 6.7M17.3 17.3L18.7 18.7M5.3 18.7L6.7 17.3M17.3 6.7L18.7 5.3"
        },
        "moon": {
            s: "M20.5 13.8A8.5 8.5 0 1 1 10.2 3.5A6.8 6.8 0 0 0 20.5 13.8Z"
        },
        "close": {
            s: "M6.5 6.5L17.5 17.5M17.5 6.5L6.5 17.5"
        },
        "prev": {
            s: "M14.5 6L8.5 12L14.5 18"
        },
        "next": {
            s: "M9.5 6L15.5 12L9.5 18"
        },
        "play": {
            f: "M7.5 5.7V18.3C7.5 19.1 8.4 19.6 9.1 19.1L19 12.8"
               + "C19.6 12.4 19.6 11.6 19 11.2L9.1 4.9C8.4 4.4 7.5 4.9 7.5 5.7Z"
        },
        "pause": {
            f: "M7 5H9.5A1 1 0 0 1 10.5 6V18A1 1 0 0 1 9.5 19H7A1 1 0 0 1 6 18V6A1 1 0 0 1 7 5Z"
               + "M14.5 5H17A1 1 0 0 1 18 6V18A1 1 0 0 1 17 19H14.5A1 1 0 0 1 13.5 18V6"
               + "A1 1 0 0 1 14.5 5Z"
        }
    })

    readonly property var def: paths[name] || ({})
    readonly property real k: size / 24

    Shape {
        anchors.fill: parent
        preferredRendererType: Shape.CurveRenderer

        ShapePath {
            fillColor: root.color
            strokeColor: "transparent"
            strokeWidth: 0
            scale: Qt.size(root.k, root.k)
            PathSvg { path: root.def.f || "" }
        }

        ShapePath {
            fillColor: "transparent"
            strokeColor: root.color
            strokeWidth: root.strokeWidth
            capStyle: ShapePath.RoundCap
            joinStyle: ShapePath.RoundJoin
            scale: Qt.size(root.k, root.k)
            PathSvg { path: root.def.s || "" }
        }
    }
}
