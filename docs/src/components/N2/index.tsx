import React, { memo, useCallback, useEffect, useRef, useState } from "react";
import rendered from "@site/src/data/n2/rendered.json";
import styles from "./styles.module.css";

type Diagram = { width: number; height: number; svg: string };
const DIAGRAMS = rendered as Record<string, Diagram>;

/** The SVG itself, never re-rendered, so React leaves its DOM alone. */
const Svg = memo(function Svg({ html }: { html: string }) {
    return <div className={styles.svg} dangerouslySetInnerHTML={{ __html: html }} />;
});

type Tip = {
    from: string;
    to: string;
    vars: [string, string][];
    x: number;
    y: number;
    /** Rendered radius of the visible circle, for the highlight. */
    r: number;
    align: "start" | "center" | "end";
    pinned: boolean;
};

/**
 * An N2 diagram of an example's model, drawn inline. Hovering, focusing, or
 * tapping a coupling circle shows the variables it carries; the tip follows
 * the mouse, and a tap pins it until the next tap elsewhere. The source is a
 * Prelude diagram in src/data/n2/<name>.json; `npm run n2` renders it.
 */
export default function N2({ name, alt }: { name: string; alt: string }): React.JSX.Element {
    const diagram = DIAGRAMS[name];
    const plotRef = useRef<HTMLDivElement>(null);
    const activeRef = useRef<SVGCircleElement | null>(null);
    // Touch has no hover, so a tap pins the tip; a mouse or pen just hovers.
    const pointerRef = useRef<string>("mouse");
    const [tip, setTip] = useState<Tip | null>(null);

    const show = useCallback((hit: SVGCircleElement, pinned: boolean) => {
        const plot = plotRef.current;
        if (!plot) return;
        activeRef.current = hit;
        const p = plot.getBoundingClientRect();
        const c = (hit.previousElementSibling ?? hit).getBoundingClientRect();
        const x = c.left + c.width / 2 - p.left;
        const fraction = x / p.width;
        setTip({
            from: hit.dataset.from ?? "",
            to: hit.dataset.to ?? "",
            vars: JSON.parse(hit.dataset.vars ?? "[]"),
            x,
            y: c.top + c.height / 2 - p.top,
            r: c.width / 2,
            align: fraction < 0.25 ? "start" : fraction > 0.75 ? "end" : "center",
            pinned,
        });
    }, []);

    const hide = useCallback(() => {
        activeRef.current = null;
        setTip(null);
    }, []);

    // Make each coupling reachable by keyboard and named for screen readers.
    useEffect(() => {
        const plot = plotRef.current;
        if (!plot) return;
        for (const hit of Array.from(plot.querySelectorAll<SVGCircleElement>(".n2-hit"))) {
            const vars = (JSON.parse(hit.dataset.vars ?? "[]") as [string, string][]).map(([o]) => o);
            hit.setAttribute("tabindex", "0");
            hit.setAttribute("role", "button");
            hit.setAttribute("aria-label", `${hit.dataset.from} to ${hit.dataset.to}: ${vars.join(", ")}`);
        }
    }, [name]);

    // A pinned tip closes on a click elsewhere or on Escape.
    useEffect(() => {
        if (!tip?.pinned) return;
        const onDown = (e: PointerEvent) => {
            if (!(e.target as Element).closest?.(".n2-hit")) hide();
        };
        const onKey = (e: KeyboardEvent) => e.key === "Escape" && hide();
        document.addEventListener("pointerdown", onDown);
        document.addEventListener("keydown", onKey);
        return () => {
            document.removeEventListener("pointerdown", onDown);
            document.removeEventListener("keydown", onKey);
        };
    }, [tip?.pinned, hide]);

    if (!diagram) throw new Error(`No rendered N2 diagram "${name}"; run npm run n2.`);

    const hitOf = (e: React.SyntheticEvent) => (e.target as Element).closest?.(".n2-hit") as SVGCircleElement | null;

    return (
        <figure className={styles.figure} style={{ maxWidth: diagram.width }}>
            <div
                ref={plotRef}
                className={styles.plot}
                role="group"
                aria-label={alt}
                onPointerDown={(e) => {
                    pointerRef.current = e.pointerType;
                }}
                onPointerOver={(e) => {
                    pointerRef.current = e.pointerType;
                    const hit = hitOf(e);
                    if (hit && e.pointerType !== "touch") show(hit, false);
                }}
                onPointerOut={(e) => {
                    if (hitOf(e) && e.pointerType !== "touch") hide();
                }}
                onFocus={(e) => {
                    const hit = hitOf(e);
                    if (hit) show(hit, false);
                }}
                onBlur={(e) => {
                    if (hitOf(e) && !tip?.pinned) hide();
                }}
                onClick={(e) => {
                    const hit = hitOf(e);
                    if (!hit || pointerRef.current !== "touch") return;
                    if (tip?.pinned && activeRef.current === hit) hide();
                    else show(hit, true);
                }}
            >
                <Svg html={diagram.svg} />
                {tip && (
                    <span
                        className={styles.dot}
                        style={{ left: tip.x, top: tip.y, width: 2 * tip.r, height: 2 * tip.r }}
                    />
                )}
                {tip && (
                    <div
                        className={styles.tip}
                        data-align={tip.align}
                        style={{ left: tip.x, top: tip.y }}
                        role="status"
                    >
                        <div className={styles.tipHead}>
                            {tip.from} → {tip.to}
                        </div>
                        <ul className={styles.tipVars}>
                            {tip.vars.map(([out, inp]) => (
                                <li key={`${out}-${inp}`}>
                                    <code>{out}</code>
                                    {out !== inp && (
                                        <>
                                            {" → "}
                                            <code>{inp}</code>
                                        </>
                                    )}
                                </li>
                            ))}
                        </ul>
                    </div>
                )}
            </div>
        </figure>
    );
}
