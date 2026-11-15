import React, { useId, useMemo, useRef, useState } from "react";
import styles from "./charts.module.css";

export type Role = "1" | "2" | "ref";

export interface Series {
    id: string;
    label: string;
    x: number[];
    y: number[];
    role: Role;
    dashed?: boolean;
    /** Draw a marker at every point (for short series). */
    markers?: boolean;
}

export interface Axis {
    domain: [number, number];
    ticks: number[];
    label: string;
    log?: boolean;
    format?: (v: number) => string;
}

export interface ChartProps {
    title: string;
    series: Series[];
    x: Axis;
    y: Axis;
    width?: number;
    height?: number;
    /** Extra SVG drawn in plot coordinates beneath the data. */
    underlay?: (sx: (v: number) => number, sy: (v: number) => number) => React.ReactNode;
    /** Number formatting for the tooltip and data table. */
    fx?: (v: number) => string;
    fy?: (v: number) => string;
    caption?: React.ReactNode;
    /** Hide the "Show data" table (e.g. very long series). */
    noTable?: boolean;
}

const M = { top: 14, right: 18, bottom: 46, left: 58 };
const ROLE_CLASS: Record<Role, string> = {
    "1": styles.role1,
    "2": styles.role2,
    ref: styles.roleRef,
};

const defaultFmt = (v: number) => {
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e5)) return v.toExponential(2);
    return Number(v.toPrecision(4)).toString();
};

export function Legend({ series }: { series: Pick<Series, "id" | "label" | "role" | "dashed">[] }) {
    if (series.length < 2) return null;
    return (
        <ul className={styles.legend}>
            {series.map((s) => (
                <li key={s.id}>
                    <svg width="22" height="10" aria-hidden="true">
                        <line
                            x1="1" y1="5" x2="21" y2="5"
                            className={`${styles.line} ${ROLE_CLASS[s.role]}`}
                            strokeDasharray={s.dashed ? "4 3" : undefined}
                        />
                    </svg>
                    {s.label}
                </li>
            ))}
        </ul>
    );
}

export default function Chart({
    title,
    series,
    x,
    y,
    width = 640,
    height = 360,
    underlay,
    fx = defaultFmt,
    fy = defaultFmt,
    caption,
    noTable,
}: ChartProps) {
    const clipId = useId().replace(/:/g, "");
    const svgRef = useRef<SVGSVGElement>(null);
    const [hover, setHover] = useState<null | { s: number; i: number }>(null);

    const iw = width - M.left - M.right;
    const ih = height - M.top - M.bottom;

    const tx = (v: number) => (x.log ? Math.log10(v) : v);
    const ty = (v: number) => (y.log ? Math.log10(v) : v);
    const [x0, x1] = x.domain.map(tx);
    const [y0, y1] = y.domain.map(ty);
    const sx = (v: number) => M.left + ((tx(v) - x0) / (x1 - x0)) * iw;
    const sy = (v: number) => M.top + ih - ((ty(v) - y0) / (y1 - y0)) * ih;

    const paths = useMemo(
        () =>
            series.map((s) =>
                s.x
                    .map((xv, i) => `${i ? "L" : "M"}${sx(xv).toFixed(1)} ${sy(s.y[i]).toFixed(1)}`)
                    .join(""),
            ),
        // eslint-disable-next-line react-hooks/exhaustive-deps
        [series, width, height],
    );

    function onMove(e: React.PointerEvent<SVGRectElement>) {
        const svg = svgRef.current;
        if (!svg) return;
        const r = svg.getBoundingClientRect();
        const px = ((e.clientX - r.left) / r.width) * width;
        const py = ((e.clientY - r.top) / r.height) * height;
        let best = { d: Infinity, s: 0, i: 0 };
        series.forEach((s, si) =>
            s.x.forEach((xv, i) => {
                const d = (sx(xv) - px) ** 2 + (sy(s.y[i]) - py) ** 2;
                if (d < best.d) best = { d, s: si, i };
            }),
        );
        setHover(best.d < 40 ** 2 ? { s: best.s, i: best.i } : null);
    }

    const hs = hover ? series[hover.s] : null;
    const hx = hs ? sx(hs.x[hover!.i]) : 0;
    const hy = hs ? sy(hs.y[hover!.i]) : 0;
    const fmtTick = (a: Axis) => a.format ?? ((v: number) => (a.log ? `1e${Math.log10(v)}` : `${v}`));

    return (
        <figure className={styles.figure}>
            <Legend series={series} />
            <div className={styles.plot}>
                <svg
                    ref={svgRef}
                    viewBox={`0 0 ${width} ${height}`}
                    role="img"
                    aria-label={title}
                    className={styles.svg}
                >
                    <defs>
                        <clipPath id={clipId}>
                            <rect x={M.left} y={M.top} width={iw} height={ih} />
                        </clipPath>
                    </defs>
                    {y.ticks.map((t) => (
                        <g key={`y${t}`}>
                            <line x1={M.left} x2={M.left + iw} y1={sy(t)} y2={sy(t)} className={styles.grid} />
                            <text x={M.left - 8} y={sy(t)} dy="0.32em" textAnchor="end" className={styles.tick}>
                                {fmtTick(y)(t)}
                            </text>
                        </g>
                    ))}
                    {x.ticks.map((t) => (
                        <text key={`x${t}`} x={sx(t)} y={M.top + ih + 18} textAnchor="middle" className={styles.tick}>
                            {fmtTick(x)(t)}
                        </text>
                    ))}
                    <line x1={M.left} x2={M.left + iw} y1={M.top + ih} y2={M.top + ih} className={styles.axis} />
                    <text x={M.left + iw / 2} y={height - 6} textAnchor="middle" className={styles.axisLabel}>
                        {x.label}
                    </text>
                    <text
                        transform={`translate(14 ${M.top + ih / 2}) rotate(-90)`}
                        textAnchor="middle"
                        className={styles.axisLabel}
                    >
                        {y.label}
                    </text>
                    <g clipPath={`url(#${clipId})`}>
                        {underlay?.(sx, sy)}
                        {series.map((s, si) => (
                            <g key={s.id}>
                                <path
                                    d={paths[si]}
                                    className={`${styles.line} ${ROLE_CLASS[s.role]}`}
                                    strokeDasharray={s.dashed ? "5 4" : undefined}
                                />
                                {s.markers &&
                                    s.x.map((xv, i) => (
                                        <circle
                                            key={i}
                                            cx={sx(xv)}
                                            cy={sy(s.y[i])}
                                            r={3.5}
                                            className={`${styles.marker} ${ROLE_CLASS[s.role]}`}
                                        />
                                    ))}
                            </g>
                        ))}
                    </g>
                    {hs && (
                        <g pointerEvents="none">
                            <line x1={hx} x2={hx} y1={M.top} y2={M.top + ih} className={styles.crosshair} />
                            <circle cx={hx} cy={hy} r={5} className={`${styles.focus} ${ROLE_CLASS[hs.role]}`} />
                        </g>
                    )}
                    <rect
                        x={M.left}
                        y={M.top}
                        width={iw}
                        height={ih}
                        fill="transparent"
                        onPointerMove={onMove}
                        onPointerLeave={() => setHover(null)}
                    />
                </svg>
                {hs && (
                    <div
                        className={styles.tooltip}
                        style={{
                            left: `${(hx / width) * 100}%`,
                            top: `${(hy / height) * 100}%`,
                            transform: `translate(${hx > width * 0.65 ? "calc(-100% - 12px)" : "12px"}, -50%)`,
                        }}
                    >
                        {series.length > 1 && <div className={styles.ttLabel}>{hs.label}</div>}
                        <div>
                            {x.label}: <b>{fx(hs.x[hover!.i])}</b>
                        </div>
                        <div>
                            {y.label}: <b>{fy(hs.y[hover!.i])}</b>
                        </div>
                    </div>
                )}
            </div>
            {caption && <figcaption className={styles.caption}>{caption}</figcaption>}
            {!noTable && (
                <details className={styles.data}>
                    <summary>Show data</summary>
                    <div className={styles.tableWrap}>
                        {series.map((s) => (
                            <table key={s.id}>
                                <caption>{s.label}</caption>
                                <thead>
                                    <tr>
                                        <th>{x.label}</th>
                                        <th>{y.label}</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {s.x.map((xv, i) => (
                                        <tr key={i}>
                                            <td>{fx(xv)}</td>
                                            <td>{fy(s.y[i])}</td>
                                        </tr>
                                    ))}
                                </tbody>
                            </table>
                        ))}
                    </div>
                </details>
            )}
        </figure>
    );
}
