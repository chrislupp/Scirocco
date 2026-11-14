import MDXComponents from "@theme-original/MDXComponents";
import Tabs from "@theme/Tabs";
import TabItem from "@theme/TabItem";
import Next from "@site/src/components/Next";
import BenchTable from "@site/src/components/BenchTable";
import QuickstartSource from "@site/src/components/QuickstartSource";
import N2 from "@site/src/components/N2";
import {
    ConvergencePlot,
    CpPlot,
    E387DragPlot,
    OnsetPlot,
    PanelPlot,
    PolarPlot,
    ShapeGradientPlot,
    TransitionSweepPlot,
} from "@site/src/components/charts";

export default {
    ...MDXComponents,
    Tabs,
    TabItem,
    Next,
    BenchTable,
    QuickstartSource,
    N2,
    CpPlot,
    PolarPlot,
    ConvergencePlot,
    PanelPlot,
    OnsetPlot,
    TransitionSweepPlot,
    E387DragPlot,
    ShapeGradientPlot,
};
