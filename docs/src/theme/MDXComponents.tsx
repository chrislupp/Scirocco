import MDXComponents from "@theme-original/MDXComponents";
import Tabs from "@theme/Tabs";
import TabItem from "@theme/TabItem";
import Next from "@site/src/components/Next";
import BenchTable from "@site/src/components/BenchTable";
import QuickstartSource from "@site/src/components/QuickstartSource";
import N2 from "@site/src/components/N2";
import DockerRun from "@site/src/components/DockerRun";
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
import {
    CodeDistributionPlot,
    CodePolarPlot,
    CodeTransitionPlot,
    ColesWadcockCpPlot,
    ExperimentPolarPlot,
    GregoryCpPlot,
    TrippedPolarPlot,
} from "@site/src/components/charts/vv";

export default {
    ...MDXComponents,
    Tabs,
    TabItem,
    Next,
    BenchTable,
    QuickstartSource,
    N2,
    DockerRun,
    CpPlot,
    PolarPlot,
    ConvergencePlot,
    PanelPlot,
    OnsetPlot,
    TransitionSweepPlot,
    E387DragPlot,
    ShapeGradientPlot,
    CodePolarPlot,
    CodeTransitionPlot,
    CodeDistributionPlot,
    ExperimentPolarPlot,
    GregoryCpPlot,
    TrippedPolarPlot,
    ColesWadcockCpPlot,
};
