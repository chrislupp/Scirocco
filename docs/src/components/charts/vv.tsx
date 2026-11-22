/**
 * Figures for the verification and validation pages. The data come from
 * scirocco/validation/vv/vv.py (see docs/README.md).
 */
import React from "react";
import Chart, { Series } from "./Chart";
import styles from "./charts.module.css";
import codes from "@site/src/data/vv_codes.json";
import experiments from "@site/src/data/vv_experiments.json";

type Code = "scirocco" | "xfoil" | "mfoil" | "flexfoil";
/** The third code drawn next to Scirocco and XFOIL. */
type Third = "mfoil" | "flexfoil";

interface Point {
    alpha: number;
    converged: boolean;
    cl?: number;
    cd?: number;
    cm?: number;
    xtr_upper?: number;
    xtr_lower?: number;
}

interface Distribution {
    alpha: number;
    x: number[];
    y: number[];
    cp: number[];
    cf?: number[];
    theta?: number[];
    dstar?: number[];
}

interface Run {
    points: Point[];
    distributions: Distribution[];
}

interface Case {
    name: string;
    airfoil: string;
    re: number;
    mach: number;
    alphas: number[];
    scirocco: Run;
    xfoil: Run;
    mfoil?: Run;
    flexfoil?: Run;
}

const CODE_STYLE: Record<Code, Pick<Series, "label" | "role" | "dashed" | "markers" | "line" | "open">> = {
    scirocco: { label: "Scirocco", role: "1" },
    xfoil: { label: "XFOIL 6.99", role: "2", dashed: true },
    mfoil: { label: "mfoil", role: "3", markers: true, line: false, open: true },
    flexfoil: { label: `flexfoil ${codes.versions.flexfoil}`, role: "3", markers: true, line: false },
};

const f2 = (v: number) => v.toFixed(2);
const f3 = (v: number) => v.toFixed(3);
const f5 = (v: number) => v.toFixed(5);

function findCase(cases: Case[], name: string): Case {
    const c = cases.find((k) => k.name === name);
    if (!c) throw new Error(`No V&V case ${name}`);
    return c;
}

const codeCase = (name: string) => findCase(codes.cases as unknown as Case[], name);
const experimentCase = (name: string) => findCase(experiments.cases as unknown as Case[], name);

function codesOf(c: Case, third: Third = "mfoil"): Code[] {
    return c[third] ? ["scirocco", "xfoil", third] : ["scirocco", "xfoil"];
}

/** Evenly spaced ticks covering [lo, hi]. */
function ticks(lo: number, hi: number, step: number): number[] {
    const out: number[] = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9; v += step) out.push(Number(v.toFixed(6)));
    return out;
}

/** Axis domain padded to whole steps around the data. */
function domain(values: number[], step: number): [number, number] {
    const lo = Math.min(...values);
    const hi = Math.max(...values);
    return [Math.floor(lo / step + 1e-9) * step, Math.ceil(hi / step - 1e-9) * step];
}

function polarSeries(
    c: Case,
    key: "cl" | "cd" | "xtr_upper",
    xKey: "alpha" | "cd" = "alpha",
    third: Third = "mfoil",
): Series[] {
    return codesOf(c, third).map((code) => {
        const pts = c[code]!.points.filter((p) => p.converged);
        return {
            id: code,
            ...CODE_STYLE[code],
            x: pts.map((p) => (xKey === "alpha" ? p.alpha : p.cd!)),
            y: pts.map((p) => p[key]!),
        };
    });
}

function missing(c: Case, third: Third = "mfoil"): string {
    const parts = codesOf(c, third)
        .map((code) => {
            const bad = c[code]!.points.filter((p) => !p.converged).map((p) => `${p.alpha}°`);
            return bad.length ? `${CODE_STYLE[code].label} at ${bad.join(", ")}` : null;
        })
        .filter(Boolean);
    return parts.length ? ` Not converged: ${parts.join("; ")}.` : " Every point converged in every code.";
}

const reText = (re: number) => {
    const e = Math.floor(Math.log10(re));
    const m = re / 10 ** e;
    return m === 1 ? `10${superscript(e)}` : `${Number(m.toPrecision(3))}×10${superscript(e)}`;
};
const superscript = (n: number) =>
    String(n)
        .split("")
        .map((d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[Number(d)])
        .join("");

function conditions(c: Case) {
    const flow = c.re ? `Re = ${reText(c.re)}` : "inviscid";
    return `${c.airfoil}, ${flow}${c.mach ? `, M = ${c.mach}` : ""}`;
}

/** Lift curve and drag polar from the three codes on shared nodes. */
export function CodePolarPlot({ name, third = "mfoil" }: { name: string; third?: Third }) {
    const c = codeCase(name);
    const lift = polarSeries(c, "cl", "alpha", third);
    const drag = polarSeries(c, "cl", "cd", third);
    const cls = lift.flatMap((s) => s.y);
    const cds = drag.flatMap((s) => s.x);
    const clDomain = domain(cls, 0.25);
    const cdStep = Math.max(...cds) > 0.02 ? 0.005 : 0.0025;
    return (
        <div className={styles.pair}>
            <Chart
                title={`Lift curve, ${conditions(c)}`}
                width={360}
                height={300}
                series={lift}
                x={{ domain: domain(c.alphas, 2), ticks: ticks(...domain(c.alphas, 2), 4), label: "α (deg)" }}
                y={{ domain: clDomain, ticks: ticks(...clDomain, 0.5), label: "Cl" }}
                fx={(v) => `${v}`}
                fy={f3}
                caption={<>{conditions(c)}.{missing(c, third)}</>}
            />
            <Chart
                title={`Drag polar, ${conditions(c)}`}
                width={360}
                height={300}
                series={drag}
                x={{
                    domain: [0, Math.ceil(Math.max(...cds) / cdStep) * cdStep],
                    ticks: ticks(0, Math.ceil(Math.max(...cds) / cdStep) * cdStep, cdStep * 2),
                    label: "Cd",
                }}
                y={{ domain: clDomain, ticks: ticks(...clDomain, 0.5), label: "Cl" }}
                fx={f5}
                fy={f3}
            />
        </div>
    );
}

/** Upper-surface transition location against angle of attack. */
export function CodeTransitionPlot({ name, third = "mfoil" }: { name: string; third?: Third }) {
    const c = codeCase(name);
    return (
        <Chart
            title={`Upper-surface transition, ${conditions(c)}`}
            height={300}
            series={polarSeries(c, "xtr_upper", "alpha", third)}
            x={{ domain: domain(c.alphas, 2), ticks: ticks(...domain(c.alphas, 2), 2), label: "α (deg)" }}
            y={{ domain: [0, 1], ticks: [0, 0.2, 0.4, 0.6, 0.8, 1], label: "x/c transition, upper" }}
            fx={(v) => `${v}`}
            fy={f3}
            caption={<>{conditions(c)}: transition location on the upper surface.</>}
        />
    );
}

/** Index of the leading-edge node (smallest x); nodes run lower TE to upper TE. */
const leIndex = (x: number[]) => x.indexOf(Math.min(...x));

/** Surface distribution of one quantity from the three codes. */
export function CodeDistributionPlot({
    name,
    quantity,
    surface = "both",
    yMax,
    third = "mfoil",
}: {
    name: string;
    quantity: "cp" | "cf";
    surface?: "both" | "upper";
    /** Clip the y axis, e.g. above the leading-edge Cf peak. */
    yMax?: number;
    third?: Third;
}) {
    const c = codeCase(name);
    const dist = (code: Code) => c[code]!.distributions[0];
    const alpha = dist("scirocco").alpha;
    const series: Series[] = codesOf(c, third)
        .filter((code) => dist(code))
        .map((code) => {
            const d = dist(code);
            const start = surface === "upper" ? leIndex(d.x) : 0;
            let x = d.x.slice(start);
            let v = d[quantity]!.slice(start);
            if (code === third) {
                // A marker on every sixth node keeps the lines beneath visible.
                const keep = (_: number, i: number) => i % 6 === 0;
                x = x.filter(keep);
                v = v.filter(keep);
            }
            return {
                id: code,
                ...CODE_STYLE[code],
                x,
                y: quantity === "cp" ? v.map((p) => -p) : v,
            };
        });
    const values = series.flatMap((s) => s.y);
    const isCp = quantity === "cp";
    const top = yMax ?? Math.max(...values);
    const yStep = isCp ? 0.5 : top > 0.008 ? 0.002 : 0.001;
    const yDom: [number, number] = isCp ? domain(values, 0.5) : [Math.floor(Math.min(0, ...values) / yStep) * yStep, top];
    const label = isCp ? "−Cp" : "Cf";
    const where = surface === "upper" ? " on the upper surface" : "";
    return (
        <Chart
            title={`${label}${where}, ${conditions(c)}, α = ${alpha}°`}
            height={330}
            series={series}
            x={{ domain: [0, 1], ticks: [0, 0.2, 0.4, 0.6, 0.8, 1], label: "x/c" }}
            y={{ domain: [yDom[0], Math.ceil(yDom[1] / yStep) * yStep], ticks: ticks(yDom[0], yDom[1], yStep * (isCp ? 1 : 2)), label }}
            fx={f3}
            fy={isCp ? f3 : f5}
            caption={
                <>
                    {conditions(c)}, α = {alpha}°{where}.
                    {yMax !== undefined && " The axis cuts off the leading-edge peak."}
                </>
            }
            noTable
        />
    );
}

const measured = experiments.measured;

/** NACA 0012 lift curve and drag polar at Re 6e6 against measurements. */
export function ExperimentPolarPlot() {
    const c = experimentCase("0012_re6e6");
    const pts = (code: Code) => c[code]!.points.filter((p) => p.converged);
    const code = (k: Code, x: (p: Point) => number, y: (p: Point) => number): Series => ({
        id: k,
        ...CODE_STYLE[k],
        x: pts(k).map(x),
        y: pts(k).map(y),
    });
    const keep = (a: number[]) => a.map((v, i) => [v, i] as const).filter(([v]) => v >= -0.5).map(([, i]) => i);
    const ab = keep(measured.abbott_cl.alpha);
    const ld = keep(measured.ladson["80 grit"].alpha);
    const dragIdx = measured.abbott_cd.cl.map((v, i) => [v, i] as const).filter(([v]) => v >= -0.05).map(([, i]) => i);
    return (
        <div className={styles.pair}>
            <Chart
                title="NACA 0012 lift curve at Re 6 million against measurements"
                width={360}
                height={320}
                series={[
                    code("scirocco", (p) => p.alpha, (p) => p.cl!),
                    code("xfoil", (p) => p.alpha, (p) => p.cl!),
                    {
                        id: "abbott",
                        label: "Abbott & von Doenhoff",
                        role: "ref",
                        line: false,
                        markers: true,
                        open: true,
                        x: ab.map((i) => measured.abbott_cl.alpha[i]),
                        y: ab.map((i) => measured.abbott_cl.cl[i]),
                    },
                    {
                        id: "ladson",
                        label: "Ladson, tripped",
                        role: "ref",
                        line: false,
                        markers: true,
                        x: ld.map((i) => measured.ladson["80 grit"].alpha[i]),
                        y: ld.map((i) => measured.ladson["80 grit"].cl[i]),
                    },
                ]}
                x={{ domain: [0, 20], ticks: [0, 4, 8, 12, 16, 20], label: "α (deg)" }}
                y={{ domain: [0, 2], ticks: [0, 0.5, 1, 1.5, 2], label: "Cl" }}
                fx={f2}
                fy={f3}
                caption={
                    <>
                        NACA 0012, Re = 6×10⁶, M = 0.15. Ladson's data are tripped with 80-grit roughness; lift
                        barely depends on the trip.
                    </>
                }
            />
            <Chart
                title="NACA 0012 drag polar at Re 6 million against measurements"
                width={360}
                height={320}
                series={[
                    code("scirocco", (p) => p.cd!, (p) => p.cl!),
                    code("xfoil", (p) => p.cd!, (p) => p.cl!),
                    {
                        id: "abbott",
                        label: "Abbott & von Doenhoff",
                        role: "ref",
                        line: false,
                        markers: true,
                        open: true,
                        x: dragIdx.map((i) => measured.abbott_cd.cd[i]),
                        y: dragIdx.map((i) => measured.abbott_cd.cl[i]),
                    },
                ]}
                x={{ domain: [0, 0.02], ticks: [0, 0.005, 0.01, 0.015, 0.02], label: "Cd" }}
                y={{ domain: [0, 2], ticks: [0, 0.5, 1, 1.5, 2], label: "Cl" }}
                fx={f5}
                fy={f3}
                caption={<>Free transition in the codes and in Abbott & von Doenhoff's smooth-model data.</>}
            />
        </div>
    );
}

/** Tripped NACA 0012 drag polar against Ladson's roughness-strip data. */
export function TrippedPolarPlot() {
    const c = experimentCase("0012_re6e6_trip");
    const code = (k: Code): Series => {
        const pts = c[k]!.points.filter((p) => p.converged);
        return { id: k, ...CODE_STYLE[k], x: pts.map((p) => p.cd!), y: pts.map((p) => p.cl!) };
    };
    const grit = (name: "80 grit" | "120 grit" | "180 grit", open: boolean): Series => {
        const m = measured.ladson[name];
        const idx = m.alpha.map((v, i) => [v, i] as const).filter(([v, i]) => v >= -0.5 && m.cd[i] <= 0.02).map(([, i]) => i);
        return {
            id: name,
            label: `Ladson, ${name}`,
            role: "ref",
            line: false,
            markers: true,
            open,
            x: idx.map((i) => m.cd[i]),
            y: idx.map((i) => m.cl[i]),
        };
    };
    return (
        <Chart
            title="NACA 0012 drag polar at Re 6 million with transition tripped at 5% chord"
            width={480}
            height={340}
            series={[code("scirocco"), code("xfoil"), grit("80 grit", false), grit("180 grit", true)]}
            x={{ domain: [0, 0.02], ticks: [0, 0.005, 0.01, 0.015, 0.02], label: "Cd" }}
            y={{ domain: [0, 2], ticks: [0, 0.5, 1, 1.5, 2], label: "Cl" }}
            fx={f5}
            fy={f3}
            caption={
                <>
                    NACA 0012, Re = 6×10⁶, M = 0.15, ncrit = 9, 199 panels, trips at x/c = 0.05 on both surfaces in
                    the codes and in Ladson's test (80- and 180-grit strips shown).{missing(c)}
                </>
            }
        />
    );
}

/** Upper-surface pressure on NACA 0012 against Gregory & O'Reilly. */
export function GregoryCpPlot({ alpha }: { alpha: 0 | 10 | 15 }) {
    const c = experimentCase("0012_re288e4");
    const exp = measured.gregory_cp[`alpha=${alpha}` as "alpha=0"];
    const series: Series[] = (["scirocco", "xfoil"] as Code[]).map((code) => {
        const d = c[code]!.distributions.find((k) => k.alpha === alpha)!;
        const le = leIndex(d.x);
        return { id: code, ...CODE_STYLE[code], x: d.x.slice(le), y: d.cp.slice(le).map((p) => -p) };
    });
    series.push({
        id: "gregory",
        label: "Gregory & O'Reilly",
        role: "ref",
        line: false,
        markers: true,
        open: true,
        x: exp.x,
        y: exp.cp.map((p) => -p),
    });
    const yDom = domain(
        series.flatMap((s) => s.y),
        alpha === 0 ? 0.5 : 2,
    );
    return (
        <Chart
            title={`NACA 0012 upper-surface pressure at ${alpha} degrees against Gregory and O'Reilly`}
            height={320}
            series={series}
            x={{ domain: [0, 1], ticks: [0, 0.2, 0.4, 0.6, 0.8, 1], label: "x/c" }}
            y={{ domain: yDom, ticks: ticks(...yDom, alpha === 0 ? 0.5 : 2), label: "−Cp" }}
            fx={f3}
            fy={f3}
            caption={
                <>
                    NACA 0012, α = {alpha}°, Re = 2.88×10⁶, upper surface. Measurements are tripped and
                    digitized from a photocopy.
                </>
            }
        />
    );
}

/** NACA 4412 near maximum lift against Coles & Wadcock. */
export function ColesWadcockCpPlot() {
    const c = experimentCase("4412_re152e4");
    const exp = measured.coles_wadcock_cp;
    const series: Series[] = (["scirocco", "xfoil"] as Code[]).map((code) => {
        const d = c[code]!.distributions[0];
        return { id: code, ...CODE_STYLE[code], x: d.x, y: d.cp.map((p) => -p) };
    });
    series.push({
        id: "coles",
        label: "Coles & Wadcock",
        role: "ref",
        line: false,
        markers: true,
        open: true,
        x: exp.x,
        y: exp.cp.map((p) => -p),
    });
    return (
        <Chart
            title="NACA 4412 pressure at 13.87 degrees against Coles and Wadcock"
            height={340}
            series={series}
            x={{ domain: [0, 1], ticks: [0, 0.2, 0.4, 0.6, 0.8, 1], label: "x/c" }}
            y={{ domain: [-1.5, 7.5], ticks: [-1, 0, 1, 2, 3, 4, 5, 6, 7], label: "−Cp" }}
            fx={f3}
            fy={f3}
            caption={<>NACA 4412, α = 13.87°, Re = 1.52×10⁶, both surfaces.</>}
        />
    );
}
