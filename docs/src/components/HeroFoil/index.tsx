import React from "react";
import cp from "@site/src/data/cp.json";
import styles from "./styles.module.css";

/*
 * The homepage figure: a real Scirocco result. Surface pressure (-Cp) over
 * NACA 2412 at 4 degrees, Re 1e6, drawn above the section and its
 * displacement-thickness envelope. All geometry comes from src/data/cp.json.
 */

const W = 720;
const X0 = 36;
const SCALE = 600; // px per chord
const px = (x: number) => X0 + x * SCALE;
// -Cp band at the top
const CP_TOP = 26;
const CP_H = 210;
const CP_MIN = -1.1;
const CP_MAX = 1.45;
const pcp = (v: number) => CP_TOP + ((CP_MAX - v) / (CP_MAX - CP_MIN)) * CP_H;
// Section below, true scale
const FOIL_Y = 330;
const pz = (z: number) => FOIL_Y - z * SCALE;

function line(xs: number[], ys: number[]) {
    return xs.map((x, i) => `${i ? "L" : "M"}${x.toFixed(1)} ${ys[i].toFixed(1)}`).join("");
}

export default function HeroFoil() {
    const { x, z } = cp.coords;
    const n = x.length;

    // Displacement-thickness envelope: offset the surface by delta* along the
    // outward normal of a lower-TE -> LE -> upper-TE loop.
    const ex: number[] = [];
    const ez: number[] = [];
    for (let i = 0; i < n; i++) {
        const a = Math.max(i - 1, 0);
        const b = Math.min(i + 1, n - 1);
        const tx = x[b] - x[a];
        const tz = z[b] - z[a];
        const len = Math.hypot(tx, tz) || 1;
        const d = cp.dstar[i] ?? 0;
        // Outward normal of this clockwise loop is (-tz, tx).
        ex.push(px(x[i] - (tz / len) * d));
        ez.push(pz(z[i] + (tx / len) * d));
    }

    const wakeN = cp.wake.x.findIndex((v: number) => v > 1.13);
    const wx = cp.wake.x.slice(0, wakeN < 0 ? undefined : wakeN);
    const wz = cp.wake.z.slice(0, wx.length);

    const cpPath = line(x.map(px), cp.cp.map((v: number) => pcp(-v)));
    const foil = line(x.map(px), z.map(pz)) + "Z";

    return (
        <figure className={styles.figure}>
            <svg
                viewBox={`0 0 ${W} 400`}
                className={styles.svg}
                role="img"
                aria-label={`Computed pressure distribution over ${cp.airfoil} at ${cp.alpha} degrees and Reynolds number one million`}
            >
                <defs>
                    <linearGradient id="hero-cp" x1="0" x2="1" y1="0" y2="0">
                        <stop offset="0" className={styles.stop0} />
                        <stop offset="0.5" className={styles.stop1} />
                        <stop offset="1" className={styles.stop2} />
                    </linearGradient>
                </defs>

                {/* -Cp axis */}
                {[-1, 0, 1].map((t) => (
                    <g key={t}>
                        <line x1={px(0)} x2={px(1.13)} y1={pcp(t)} y2={pcp(t)} className={t === 0 ? styles.zero : styles.grid} />
                        <text x={px(0) - 8} y={pcp(t)} dy="0.32em" textAnchor="end" className={styles.tick}>
                            {t}
                        </text>
                    </g>
                ))}
                <text x={px(0) - 8} y={CP_TOP - 8} textAnchor="end" className={styles.tick}>
                    −Cp
                </text>

                <path d={cpPath} className={styles.cp} pathLength={1} />

                {/* Section, viscous envelope, and wake */}
                <path d={line(ex, ez)} className={styles.envelope} />
                <path d={line(wx.map(px), wz.map(pz))} className={styles.wake} />
                <path d={foil} className={styles.foil} />

                <text x={px(0.5)} y={pcp(1.18)} className={styles.note} textAnchor="middle">
                    suction side
                </text>
                <text x={px(0.5)} y={pcp(-0.42)} className={styles.note} textAnchor="middle">
                    pressure side
                </text>
            </svg>
            <figcaption className={styles.readout}>
                <span className={styles.case}>
                    {cp.airfoil}, α = {cp.alpha}°, Re = 10⁶
                </span>
                <span>
                    C<sub>l</sub> <b>{cp.cl.toFixed(3)}</b>
                </span>
                <span>
                    C<sub>d</sub> <b>{cp.cd.toFixed(5)}</b>
                </span>
                <span>
                    C<sub>m</sub> <b>{cp.cm.toFixed(4)}</b>
                </span>
            </figcaption>
        </figure>
    );
}
