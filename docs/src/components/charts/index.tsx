/**
 * Documentation figures. Every number comes from src/data/*.json, written by
 * the `docs_figure_data` example program, except onset.json and
 * transition.json, which come from scripts under scirocco/validation (see
 * docs/README.md).
 */
import React from "react";
import Chart from "./Chart";
import styles from "./charts.module.css";
import cp from "@site/src/data/cp.json";
import polar from "@site/src/data/polar.json";
import convergence from "@site/src/data/convergence.json";
import panels from "@site/src/data/panels.json";
import onset from "@site/src/data/onset.json";
import transition from "@site/src/data/transition.json";

const f4 = (v: number) => v.toFixed(4);
const f3 = (v: number) => v.toFixed(3);

function foilPath(x: number[], z: number[], px: (v: number) => number, pz: (v: number) => number) {
    return x.map((xv, i) => `${i ? "L" : "M"}${px(xv).toFixed(1)} ${pz(z[i]).toFixed(1)}`).join("") + "Z";
}

/** Surface pressure on NACA 2412, viscous vs inviscid, with the section drawn in. */
export function CpPlot() {
    const { x, z } = cp.coords;
    const neg = (a: number[]) => a.map((v) => -v);
    return (
        <Chart
            title={`Pressure coefficient on ${cp.airfoil} at ${cp.alpha} degrees, Re 1e6`}
            series={[
                { id: "visc", label: "Viscous, Re = 10⁶", x, y: neg(cp.cp), role: "1" },
                { id: "inv", label: "Inviscid", x, y: neg(cp.cpi), role: "ref", dashed: true },
            ]}
            x={{ domain: [0, 1], ticks: [0, 0.2, 0.4, 0.6, 0.8, 1], label: "x/c" }}
            y={{ domain: [-1.1, 1.5], ticks: [-1, -0.5, 0, 0.5, 1, 1.5], label: "−Cp" }}
            fx={f4}
            fy={f3}
            underlay={(sx, sy) => (
                <path className={styles.foilFill} d={foilPath(x, z, sx, (v) => sy(-0.8 + 2.4 * v))} />
            )}
            caption={
                <>
                    {cp.airfoil} at α = {cp.alpha}°, 160 panels. Viscous: Cl = {cp.cl.toFixed(3)}, Cd ={" "}
                    {cp.cd.toFixed(5)}. The boundary layer lowers the suction peak and relieves the
                    trailing-edge pressure recovery. Section shape drawn to scale in x only.
                </>
            }
            noTable
        />
    );
}

/** Lift curve and drag polar for NACA 2412 at Re 1e6. */
export function PolarPlot() {
    const pts = polar.points.filter((p) => p.converged);
    const alpha = pts.map((p) => p.alpha);
    return (
        <div className={styles.pair}>
            <Chart
                title={`Lift curve for ${polar.airfoil}, Re 1e6, viscous and inviscid`}
                width={360}
                height={300}
                series={[
                    { id: "v", label: "Viscous", x: alpha, y: pts.map((p) => p.cl!), role: "1", markers: true },
                    {
                        id: "i",
                        label: "Inviscid",
                        x: polar.points.map((p) => p.alpha),
                        y: polar.points.map((p) => p.cl_inviscid),
                        role: "ref",
                        dashed: true,
                    },
                ]}
                x={{ domain: [-4, 14], ticks: [-4, 0, 4, 8, 12], label: "α (deg)" }}
                y={{ domain: [-0.5, 1.9], ticks: [-0.5, 0, 0.5, 1, 1.5], label: "Cl" }}
                fy={f4}
            />
            <Chart
                title={`Drag polar for ${polar.airfoil}, Re 1e6`}
                width={360}
                height={300}
                series={[{ id: "p", label: "Viscous", x: pts.map((p) => p.cd!), y: pts.map((p) => p.cl!), role: "1", markers: true }]}
                x={{
                    domain: [0.004, 0.028],
                    ticks: [0.005, 0.01, 0.015, 0.02, 0.025],
                    label: "Cd",
                    format: (v) => v.toFixed(3),
                }}
                y={{ domain: [-0.5, 1.9], ticks: [-0.5, 0, 0.5, 1, 1.5], label: "Cl" }}
                fx={(v) => v.toFixed(5)}
                fy={f4}
            />
        </div>
    );
}

/** Amplification onset ramps against Orr-Sommerfeld envelope slopes. */
export function OnsetPlot() {
    const blasius = onset.stability.reduce((a, b) =>
        Math.abs(a.H - onset.H_model) < Math.abs(b.H - onset.H_model) ? a : b,
    );
    return (
        <Chart
            title="Amplification onset: default and envelope ramps against stability theory"
            width={640}
            height={320}
            series={[
                { id: "os", label: `Orr–Sommerfeld envelope, H = ${blasius.H.toFixed(2)}`, x: blasius.z, y: blasius.r, role: "2", markers: true },
                { id: "env", label: "Envelope ramp (opt-in)", x: onset.z, y: onset.envelope, role: "1" },
                { id: "cub", label: "Cubic ramp (default)", x: onset.z, y: onset.cubic, role: "ref", dashed: true },
            ]}
            x={{ domain: [-0.4, 0.45], ticks: [-0.4, -0.2, 0, 0.2, 0.4], label: "log₁₀ Reθ − lrc(Hk)" }}
            y={{ domain: [0, 1.1], ticks: [0, 0.25, 0.5, 0.75, 1], label: "growth rate / asymptote" }}
            fx={f3}
            fy={f3}
            caption={
                <>
                    Rate of growth of the amplification factor, relative to its asymptotic value, near
                    the critical Reynolds number. Points: spatial Orr–Sommerfeld envelope for the Blasius
                    profile. Growth starts below the cubic ramp and approaches the asymptote slowly; the
                    envelope ramp follows it and keeps the same asymptote, so amplification far
                    downstream is unchanged.
                </>
            }
        />
    );
}

/** Residual history of the coupled Newton solve. */
export function ConvergencePlot({ variant = "typical" }: { variant?: "typical" | "hard" }) {
    const idx = (a: number[]) => a.map((_, i) => i);
    const fmt = (v: number) => v.toExponential(2);
    if (variant === "typical") {
        const r = convergence.ptc.residual;
        return (
            <Chart
                title="Residual history, NACA 2412 at 4 degrees, Re 1e6"
                width={640}
                height={300}
                series={[{ id: "r", label: "Residual", x: idx(r), y: r, role: "1", markers: true }]}
                x={{ domain: [0, r.length - 1], ticks: idx(r), label: "Newton iteration" }}
                y={{ domain: [1e-13, 1e2], ticks: [1e-12, 1e-9, 1e-6, 1e-3, 1], label: "‖R‖", log: true }}
                fx={(v) => `${v}`}
                fy={fmt}
                caption={`NACA 2412, α = 4°, Re = 10⁶, 160 panels. From a cold start the residual falls twelve orders of magnitude in ${r.length - 1} Newton steps.`}
            />
        );
    }
    const a = convergence.hard_ptc.residual;
    const b = convergence.hard_newton.residual;
    const n = Math.max(a.length, b.length);
    return (
        <Chart
            title="Residual history for a difficult start, with and without PTC"
            width={640}
            height={320}
            series={[
                { id: "ptc", label: `PTC (default), ${a.length - 1} steps`, x: idx(a), y: a, role: "1" },
                { id: "newton", label: `Plain Newton, ${b.length - 1} steps`, x: idx(b), y: b, role: "2", dashed: true },
            ]}
            x={{ domain: [0, n - 1], ticks: [0, 10, 20, 30, 40, 50].filter((t) => t < n), label: "Newton iteration" }}
            y={{ domain: [1e-13, 1e2], ticks: [1e-12, 1e-9, 1e-6, 1e-3, 1], label: "‖R‖", log: true }}
            fx={(v) => `${v}`}
            fy={fmt}
            caption="NACA 0012, α = 5°, Re = 10⁶, 199 panels. The residual stalls near 1 while the amplification march struggles at one station, then Newton takes over. PTC reaches 10⁻¹¹ first, then hovers just above the tolerance for ten more steps."
        />
    );
}

/** Panel nodes on NACA 2412 at 60 panels, drawn to scale. */
export function PanelPlot() {
    const { x, z } = panels.coords;
    const W = 640;
    const pad = 14;
    const scale = (W - 2 * pad) / 1.0;
    const zmax = Math.max(...z) + 0.01;
    const zmin = Math.min(...z) - 0.01;
    const H = Math.round((zmax - zmin) * scale + 2 * pad);
    const px = (v: number) => pad + v * scale;
    const pz = (v: number) => pad + (zmax - v) * scale;
    return (
        <figure className={styles.figure}>
            <svg viewBox={`0 0 ${W} ${H}`} className={styles.svg} role="img" aria-label="NACA 2412 paneled with 60 panels; nodes cluster at the leading and trailing edges">
                <path className={styles.panelLine} d={x.map((xv, i) => `${i ? "L" : "M"}${px(xv).toFixed(1)} ${pz(z[i]).toFixed(1)}`).join("")} />
                {x.map((xv, i) => (
                    <circle key={i} className={styles.node} cx={px(xv)} cy={pz(z[i])} r={2.6} />
                ))}
            </svg>
            <figcaption className={styles.caption}>
                NACA 2412 with 60 panels ({x.length} nodes). The paneler clusters nodes where curvature is high,
                at the leading and trailing edges.
            </figcaption>
        </figure>
    );
}

/** Exact dCd/dalpha over a small alpha sweep, cut-cell vs captured transition. */
export function TransitionSweepPlot() {
    const { cutcell, captured } = transition.sweep;
    const e4 = (a: number[]) => a.map((v) => v * 1e4);
    return (
        <Chart
            title="dCd/dα over a small α sweep, NACA 2412, Re 1e6: cut-cell and captured transition"
            width={640}
            height={300}
            series={[
                { id: "capt", label: "Captured transition (opt-in)", x: captured.alpha, y: e4(captured.dcd), role: "1", markers: true },
                { id: "cut", label: "Cut cell (default)", x: cutcell.alpha, y: e4(cutcell.dcd), role: "ref", dashed: true, markers: true },
            ]}
            x={{ domain: [3.9, 4.3], ticks: [3.9, 4.0, 4.1, 4.2, 4.3], label: "α (deg)" }}
            y={{ domain: [3, 10.5], ticks: [4, 6, 8, 10], label: "dCd/dα (10⁻⁴ per degree)" }}
            fx={(v) => v.toFixed(2)}
            fy={(v) => v.toFixed(3)}
            caption={
                <>
                    Complex-step derivative through mfoil, 160 panels, 21 angles. With the cut cell the
                    derivative jumps each time the upper-surface transition point crosses a node. The
                    captured model, at its default transition length of 60θ, gives a smooth curve.
                </>
            }
        />
    );
}

/** Cd error against the UIUC Eppler 387 drag polars at matched Cl. */
export function E387DragPlot() {
    const reLabel = (re: number) => `${(re / 1e5).toFixed(re % 1e5 ? 1 : 0)}×10⁵`;
    return (
        <div className={styles.pair}>
            {transition.polars.map((p) => (
                <Chart
                    key={p.Re}
                    title={`Eppler 387 drag error at Re ${reLabel(p.Re)}`}
                    width={360}
                    height={260}
                    series={[
                        { id: "capt", label: "Captured, L = 60θ", x: p.captured.cl, y: p.captured.err, role: "1", markers: true },
                        { id: "cut", label: "Cut cell", x: p.cutcell.cl, y: p.cutcell.err, role: "ref", dashed: true, markers: true },
                        { id: "exp", label: "Measured", x: [0.1, 1.0], y: [0, 0], role: "2" },
                    ]}
                    x={{ domain: [0, 1.05], ticks: [0, 0.25, 0.5, 0.75, 1], label: "Cl" }}
                    y={{ domain: [-40, 15], ticks: [-40, -30, -20, -10, 0, 10], label: `Cd error (%), Re ${reLabel(p.Re)}` }}
                    fx={(v) => v.toFixed(3)}
                    fy={(v) => v.toFixed(1)}
                />
            ))}
        </div>
    );
}
