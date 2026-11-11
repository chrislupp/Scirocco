import React from "react";
import bench from "@site/src/data/bench.json";

const ms = (v: number) => (v >= 1000 ? `${(v / 1000).toFixed(1)} s` : `${v.toFixed(v < 10 ? 1 : 0)} ms`);

/** Phase timings from `scirocco/examples/bench.cpp`. */
export default function BenchTable() {
    return (
        <table>
            <thead>
                <tr>
                    <th>Panels</th>
                    <th>Paneling</th>
                    <th>Inviscid solve</th>
                    <th>Viscous solve</th>
                    <th>Full gradient</th>
                </tr>
            </thead>
            <tbody>
                {bench.rows.map((r) => (
                    <tr key={r.panels}>
                        <td>{r.panels}</td>
                        <td>{ms(r.paneling)}</td>
                        <td>{ms(r.inviscid)}</td>
                        <td>{ms(r.viscous)}</td>
                        <td>{ms(r.gradient)}</td>
                    </tr>
                ))}
            </tbody>
        </table>
    );
}
