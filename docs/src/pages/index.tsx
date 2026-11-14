import React from "react";
import Link from "@docusaurus/Link";
import Layout from "@theme/Layout";
import CodeBlock from "@theme/CodeBlock";
import HeroFoil from "@site/src/components/HeroFoil";
import { useQuickstartCode } from "@site/src/components/QuickstartSource";
import cp from "@site/src/data/cp.json";
import bench from "@site/src/data/bench.json";
import styles from "./index.module.css";

const viscous149 = bench.rows.find((r) => r.panels === 149)!;

/** Format like std::cout's default (six significant digits). */
const g6 = (v: number) => Number(v.toPrecision(6)).toString();

function Hero() {
    return (
        <header className={styles.hero}>
            <div className={styles.heroText}>
                <h1 className={styles.title}>Airfoil analysis for optimization.</h1>
                <p className={styles.lede}>
                    Scirocco couples a panel method to an integral boundary layer with e<sup>N</sup> transition.
                    One solve gives you C<sub>l</sub>, C<sub>d</sub>, and C<sub>m</sub>. One more call gives
                    their derivatives with respect to every coordinate, angle of attack, Reynolds number, and
                    Mach number.
                </p>
                <div className={styles.actions}>
                    <Link className={styles.primary} to="/docs/getting-started/quickstart">
                        Get started
                    </Link>
                    <Link className={styles.secondary} to="/docs/python">
                        Use it from Python
                    </Link>
                </div>
            </div>
            <div className={styles.heroFigure}>
                <HeroFoil />
            </div>
        </header>
    );
}

function Program() {
    const code = useQuickstartCode();
    return (
        <section className={styles.section}>
            <div className={styles.sectionHead}>
                <h2>A viscous solution in three calls</h2>
            </div>
            <div className={styles.program}>
                <div className={styles.code}>
                    <CodeBlock language="cpp" title="quickstart.cpp">
                        {code}
                    </CodeBlock>
                </div>
                <div className={styles.result}>
                    <CodeBlock language="output">
                        {`Cl = ${g6(cp.cl)}\nCd = ${g6(cp.cd)}\nCm = ${g6(cp.cm)}`}
                    </CodeBlock>
                    <p className={styles.timing}>
                        The whole program runs in about a quarter of a second on an Apple M1 Pro. At 149
                        panels, paneling takes {Math.round(viscous149.paneling)} ms and the viscous solve{" "}
                        {Math.round(viscous149.viscous)} ms.
                    </p>
                    <Link className={styles.more} to="/docs/concepts/performance">
                        See all timings
                    </Link>
                </div>
            </div>
        </section>
    );
}

const FEATURES = [
    {
        title: "XFOIL-like viscous solver",
        body: "A linear-vortex panel method coupled by Newton's method to an integral boundary layer and wake, with eᴺ transition and Kármán–Tsien compressibility.",
        to: "/docs/concepts/panel-method",
        link: "How the solver works",
    },
    {
        title: "Verified exact gradients",
        body: "The discrete adjoint returns dCl, dCd, dCm with respect to every panel node, α, Re, and Ma. Each piece of the chain is tested against complex-step reference data from mfoil.",
        to: "/docs/guides/gradients",
        link: "Compute gradients",
    },
    {
        title: "Framework-agnostic MDO interface",
        body: "Connect to Scirocco through Philote, a framework-neutral interface for multidisciplinary design optimization.",
        to: "/docs/python",
        link: "Python and optimization",
    },
];

function Features() {
    return (
        <section className={styles.section}>
            <div className={styles.features}>
                {FEATURES.map((f) => (
                    <div key={f.title} className={styles.feature}>
                        <h3>{f.title}</h3>
                        <p>{f.body}</p>
                        <Link to={f.to}>{f.link}</Link>
                    </div>
                ))}
            </div>
        </section>
    );
}

export default function Home(): React.JSX.Element {
    return (
        <Layout
            title="Airfoil analysis for optimization"
            description="Scirocco is a C++17 library for viscous and inviscid 2D airfoil analysis with adjoint gradients, HDF5 output, and a framework-neutral MDO interface through Philote."
        >
            <main className={styles.main}>
                <Hero />
                <Program />
                <Features />
            </main>
        </Layout>
    );
}
