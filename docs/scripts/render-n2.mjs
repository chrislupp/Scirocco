/**
 * Render the N2 diagrams in src/data/n2/*.json into src/data/n2/rendered.json,
 * which the <N2> component draws inline.
 *
 * The JSON files are Prelude diagrams: open one in Prelude to edit it, save,
 * and rerun `npm run n2`. Layout and drawing come from mdo-n2, the library
 * Prelude draws with. Colors are left to CSS (src/components/N2), so one SVG
 * serves both themes, and each coupling circle carries the variables it
 * stands for, which the component shows on hover, focus, or tap.
 */
import { JSDOM, VirtualConsole } from "jsdom";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const dom = new JSDOM("", { virtualConsole: new VirtualConsole() });
const { document } = dom.window;
Object.assign(globalThis, { document, XMLSerializer: dom.window.XMLSerializer });

// jsdom has no canvas, and without one mdo-n2 guesses 0.7 em per character,
// which leaves every box far too wide. Measure with Helvetica's character
// widths instead (Arial shares them); the labels are drawn in that font.
// prettier-ignore
const HELVETICA = [ // advance widths for U+0020..U+007E, in 1/1000 em
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278,
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556,
    1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778,
    667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556,
    333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556,
    556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584,
];
dom.window.HTMLCanvasElement.prototype.getContext = function () {
    let size = 13;
    return {
        set font(f) {
            size = Number(/([\d.]+)px/.exec(f)?.[1] ?? size);
        },
        measureText(text) {
            let w = 0;
            for (const ch of text) w += HELVETICA[ch.codePointAt(0) - 32] ?? 556;
            return { width: (w / 1000) * size };
        },
    };
};

const { convert, validateInput, exportSVG } = await import("mdo-n2");

const root = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const srcDir = path.join(root, "src/data/n2");
const outFile = path.join(srcDir, "rendered.json");

// One scale for every diagram: pages show each SVG at its natural size, so
// labels match across figures. The widest (multipoint) must fit the column.
const STYLE = {
    text: { fontFamily: "Helvetica, Arial, sans-serif", fontSize: 13, fontWeight: "500" },
    dimensions: { paddingHorizontal: 28, paddingVertical: 24 },
    connection: { lineWidth: 8 },
    couplingPoint: { radius: 6 },
};

/** Collect every discipline, nested ones included, by id. */
function index(disciplines, out = new Map()) {
    for (const d of disciplines) {
        out.set(d.id, d);
        if (d.disciplines) index(d.disciplines, out);
    }
    return out;
}

/** Strip colors and duplicates from mdo-n2's SVG and annotate the circles. */
function postprocess(svgText, elements) {
    const host = document.createElement("div");
    host.innerHTML = svgText;
    const svg = host.querySelector("svg");
    const width = Number(svg.getAttribute("width"));
    const height = Number(svg.getAttribute("height"));
    svg.removeAttribute("width");
    svg.removeAttribute("height");

    // One polyline per connection: pairs that share a route draw it twice.
    const seen = new Set();
    for (const line of svg.querySelectorAll("#connections polyline")) {
        const key = line.getAttribute("points");
        if (seen.has(key)) line.remove();
        else seen.add(key);
        line.removeAttribute("stroke");
        line.removeAttribute("stroke-width");
    }

    // Boxes and labels: color by element type in CSS.
    for (const rect of svg.querySelectorAll("rect[id]")) {
        const type = elements.get(rect.id)?.type ?? "explicit";
        rect.setAttribute("class", `n2-box n2-${type}`);
        rect.removeAttribute("style");
        rect.removeAttribute("id");
    }
    for (const text of svg.querySelectorAll("text[data-element-id]")) {
        const type = elements.get(text.getAttribute("data-element-id"))?.type ?? "explicit";
        text.setAttribute("class", `n2-label n2-${type}`);
        text.removeAttribute("style");
        text.removeAttribute("data-element-id");
    }

    // One circle per coupling, listing every variable it carries.
    const groups = new Map();
    for (const c of svg.querySelectorAll("#coupling circle")) {
        const key = [c.getAttribute("cx"), c.getAttribute("cy"), c.dataset.source, c.dataset.target].join("|");
        const vars = [c.dataset.output, c.dataset.input];
        if (groups.has(key)) {
            groups.get(key).vars.push(vars);
            c.remove();
        } else {
            groups.set(key, { circle: c, vars: [vars] });
        }
    }
    for (const { circle, vars } of groups.values()) {
        const from = elements.get(circle.dataset.source)?.name ?? circle.dataset.source;
        const to = elements.get(circle.dataset.target)?.name ?? circle.dataset.target;
        for (const a of ["stroke", "stroke-width", "fill", "data-source", "data-target", "data-output", "data-input"]) {
            circle.removeAttribute(a);
        }
        circle.setAttribute("class", "n2-coupling");
        // A larger, invisible circle on top takes the pointer and carries the data.
        const hit = circle.cloneNode();
        hit.setAttribute("r", "14");
        hit.setAttribute("class", "n2-hit");
        hit.setAttribute("data-from", from);
        hit.setAttribute("data-to", to);
        hit.setAttribute("data-vars", JSON.stringify(vars));
        circle.after(hit);
    }

    for (const el of svg.querySelectorAll("[data-element-id]")) el.removeAttribute("data-element-id");

    // Ids would repeat when a page shows two diagrams; the CSS uses classes.
    for (const g of svg.querySelectorAll("g[id]")) {
        g.setAttribute("class", `n2-${g.id}`);
        g.removeAttribute("id");
    }

    return { width, height, svg: host.innerHTML };
}

const rendered = {};
let failed = false;
for (const file of fs.readdirSync(srcDir).filter((f) => f.endsWith(".json") && f !== "rendered.json").sort()) {
    const name = path.basename(file, ".json");
    const input = JSON.parse(fs.readFileSync(path.join(srcDir, file), "utf8"));
    const { valid, errors } = validateInput(input);
    if (!valid) {
        console.error(`${file}: ${errors.join("; ")}`);
        failed = true;
        continue;
    }
    const { data, warnings } = convert(input);
    if (warnings.length) {
        console.error(`${file}: ${warnings.join("; ")}`);
        failed = true;
    }
    data.style = STYLE;
    data.annotate = true;
    rendered[name] = postprocess(exportSVG(data, { padding: 8 }), index(input.disciplines));
    console.log(`rendered ${name} (${rendered[name].width} x ${rendered[name].height})`);
}
fs.writeFileSync(outFile, JSON.stringify(rendered, null, 1) + "\n");
process.exit(failed ? 1 : 0);
