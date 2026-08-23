"""Static HTML report: a downstream consumer of the treeline API (not part of it).

Renders every annotated h5ad in a results directory — per-sample views, the integrated
joint, and any refined within-class views — with the coarse-to-fine slider, hierarchical
colors from `treeline.colors.palette`, per-cluster path tables and NS-Forest substates.
Vanilla JS on canvas; no server, no frameworks.

    .venv/bin/python apps/report.py results/poc report.html
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import scanpy as sc

from treeline.annotate import annotations
from treeline.colors import palette
from treeline.harmonize import observed_paths, path_tree


def build_data(results: Path) -> dict:
    samples, calls, suffixes, anns = {}, {}, {}, []
    for f in sorted(results.glob("*_annotated.h5ad")):
        a = sc.read_h5ad(f)
        name = f.stem.removesuffix("_annotated").replace("refined__", "refined · ")
        ann = annotations(a)
        anns.append(ann)
        keys = ann["cluster_keys"]
        um = a.obsm["X_umap"]
        entry = {
            "x": [round(float(v), 2) for v in um[:, 0]],
            "y": [round(float(v), 2) for v in um[:, 1]],
            "clusters": {k: a.obs[k].astype(str).astype(int).tolist() for k in keys},
        }
        if "sample" in a.obs.columns:  # integrated/refined views: same nuclei as the samples
            names = sorted(a.obs["sample"].unique().tolist())
            entry["sampleNames"] = names
            entry["sample"] = [names.index(s) for s in a.obs["sample"]]
            entry["integrated"] = True
        samples[name] = entry
        calls[name] = ann["calls"]
        suffixes[name] = ann["suffixes"]
    paths = observed_paths(calls)
    panels: dict = {}
    for ann in anns:
        panels.update(ann.get("panels", {}))
    return {
        "samples": samples,
        "calls": calls,
        "suffixes": suffixes,
        "panels": panels,
        "colors": palette(*anns),
        "tree": path_tree(paths),
        "maxDepth": max((len(p) for p in paths), default=1),
    }


def main() -> None:
    results, out = Path(sys.argv[1]), Path(sys.argv[2])
    data = build_data(results)
    html = TEMPLATE.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    out.write_text(html)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")


TEMPLATE = r"""<title>treeline 1619</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,400;12..96,600;12..96,800&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {
  --bg:#F5F7F3; --card:#FFFFFF; --surface:#ECF0E9; --ink:#1B2820; --muted:#5B6B60;
  --line:#D3DCD2; --pine:#2F6B4F; --pine-soft:#DCE9E0; --rock:#75818C; --rock-soft:#E4E8EC;
  --lichen:#A9761F; --lichen-soft:#F1E7D2;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#0F1512; --card:#1A241E; --surface:#151D18; --ink:#E6ECE5; --muted:#93A398;
  --line:#2A3630; --pine:#5FB08A; --pine-soft:#1E3229; --rock:#8B99A6; --rock-soft:#212A31;
  --lichen:#D6A34C; --lichen-soft:#322817;
}}
:root[data-theme="dark"] {
  --bg:#0F1512; --card:#1A241E; --surface:#151D18; --ink:#E6ECE5; --muted:#93A398;
  --line:#2A3630; --pine:#5FB08A; --pine-soft:#1E3229; --rock:#8B99A6; --rock-soft:#212A31;
  --lichen:#D6A34C; --lichen-soft:#322817;
}
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
  font-family:"Bricolage Grotesque",system-ui,sans-serif; font-size:15px; line-height:1.5; }
.mono { font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:0.88em; }
.wrap { max-width:1200px; margin:0 auto; padding:0 20px; }
h1 { font-size:1.9rem; font-weight:800; margin:1.4rem 0 0.1rem; }
h2 { font-size:1.25rem; font-weight:600; margin:2.2rem 0 0.6rem; }
.sub { color:var(--muted); margin:0 0 1rem; }
.controls { position:sticky; top:0; z-index:5; background:var(--bg); border-bottom:1px solid var(--line);
  padding:10px 0; display:flex; gap:1.6rem; align-items:center; flex-wrap:wrap; }
.ctl-label { font-size:0.72rem; font-weight:600; letter-spacing:0.1em; text-transform:uppercase; color:var(--muted); }
.seg { display:inline-flex; border:1px solid var(--line); border-radius:8px; overflow:hidden; }
.seg button { border:0; background:var(--card); color:var(--muted); padding:5px 14px; cursor:pointer;
  font:inherit; font-weight:600; }
.seg button.on { background:var(--pine); color:#fff; }
input[type=range] { accent-color:var(--pine); width:220px; vertical-align:middle; }
#depthVal { font-weight:800; color:var(--pine); min-width:1.2em; display:inline-block; text-align:center; }
.grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:16px; margin:1.2rem 0; }
.panel { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px 14px; }
.panel h3 { margin:0 0 6px; font-size:1rem; font-weight:600; display:flex; justify-content:space-between; align-items:center; }
.panel h3 .n { color:var(--muted); font-weight:400; font-size:0.85rem; }
canvas { width:100%; height:auto; display:block; border-radius:6px; background:var(--surface); }
.legend { display:flex; flex-wrap:wrap; gap:6px 14px; margin:0.6rem 0 0.2rem; }
.legend span { display:inline-flex; align-items:center; gap:6px; font-size:0.85rem; }
.legend i { width:11px; height:11px; border-radius:3px; display:inline-block; }
.legend .cnt { color:var(--muted); font-size:0.78rem; }
table { border-collapse:collapse; width:100%; font-size:0.88rem; }
th { text-align:left; font-size:0.7rem; letter-spacing:0.08em; text-transform:uppercase; color:var(--muted);
  border-bottom:1px solid var(--line); padding:6px 10px 6px 0; }
td { border-bottom:1px solid var(--line); padding:7px 10px 7px 0; vertical-align:top; }
.tablewrap { overflow-x:auto; }
.chip { display:inline-block; padding:1px 8px; border-radius:999px; font-size:0.8rem; margin:1px 2px 1px 0; }
.chip.suffix { background:var(--lichen-soft); color:var(--lichen); border:1px solid var(--lichen);
  font-family:"IBM Plex Mono",monospace; font-size:0.75rem; }
.share { color:var(--muted); font-size:0.78rem; font-variant-numeric:tabular-nums; }
.refused { color:var(--rock); background:var(--rock-soft); border-radius:6px; padding:2px 8px;
  font-size:0.8rem; display:inline-block; margin-top:3px; }
.override { color:var(--lichen); background:var(--lichen-soft); border-radius:6px; padding:2px 8px;
  font-size:0.8rem; display:inline-block; margin-top:3px; max-width:60ch; }
details { margin-top:4px; } summary { cursor:pointer; color:var(--pine); font-size:0.82rem; }
.genes { margin:6px 0 2px; } .genes b { font-weight:600; }
.gene { display:inline-block; background:var(--surface); border:1px solid var(--line); border-radius:4px;
  padding:0 6px; margin:1px 2px; font-family:"IBM Plex Mono",monospace; font-size:0.78rem; }
footer { color:var(--muted); font-size:0.82rem; border-top:1px solid var(--line); margin-top:2.5rem; padding:1rem 0 2rem; }
.tree, .tree ul { list-style:none; margin:0; padding-left:0; }
.tree ul { padding-left:1.35rem; border-left:1px solid var(--line); margin-left:0.42rem; }
.tree li { padding:2px 0; }
.tree .node { display:inline-flex; align-items:center; gap:0.5rem; }
.tree .node i { width:11px; height:11px; border-radius:3px; flex:none; }
.tree .node .nm { font-weight:600; }
.tree .node .ct { color:var(--muted); font-size:0.8rem; font-variant-numeric:tabular-nums; }
.tree .sub-node .nm { font-family:"IBM Plex Mono",monospace; font-weight:500; font-size:0.85rem; }
.tree .sub-node { border:1px solid var(--lichen); border-radius:6px; padding:0 8px; margin:1px 0; }
</style>
<div class="wrap">
  <h1>treeline · 1619 LM / MM</h1>
  <p class="sub">Multi-level annotation, evidence-gated. Slide coarse → fine; labels stop where their vote does.
  Same parent, same hue — deeper labels are darker shades.</p>
  <div class="controls">
    <span><span class="ctl-label">Resolution</span>&nbsp; <span class="seg" id="resSeg"></span></span>
    <span><span class="ctl-label">Depth</span>&nbsp; coarse <input type="range" id="depth" min="1" value="2"> fine
      <span id="depthVal"></span></span>
    <span><span class="ctl-label">Color by</span>&nbsp; <span class="seg" id="colorSeg">
      <button data-v="label" class="on">labels</button><button data-v="cluster">clusters</button></span></span>
    <span><span class="ctl-label">Integrated colors</span>&nbsp; <span class="seg" id="jointSeg">
      <button data-v="label" class="on">labels</button><button data-v="sample">sample</button></span></span>
  </div>
  <div class="grid" id="plots"></div>
  <div class="legend" id="legend"></div>
  <h2>The tree</h2>
  <p class="sub">Every label observed at the current resolution, from the root down; ochre-edged
  entries are NS-Forest substates — what grows above the treeline.</p>
  <div id="treeview"></div>
  <div id="tables"></div>
  <footer>Per-sample views are annotated independently; the integrated view is scANVI on the tree-cut
  consensus prior, reclustered by the driver and re-annotated by the same vote. Gate refusals are
  automated, self-reporting rules; ✍ overrides are signed (SPECS provenance tiers). Generated by
  apps/report.py, a consumer of the treeline API.</footer>
</div>
<script>
const D = __DATA__;
const RES = Object.keys(D.samples[Object.keys(D.samples)[0]].clusters);
const HAS_SUFFIX = Object.keys(D.suffixes || {}).length > 0;
let state = { res: RES[0], depth: 2, jointBy: "label", colorBy: "label" };
const els = {};

function pathOf(sample, res, cl) { const c = D.calls[sample][res][String(cl)]; return c ? c.path : []; }
function suffixOf(sample, res, cl) {
  return ((D.suffixes[sample] || {})[res] || {})[String(cl)] || null; }
function labelAt(path, depth) { return path.length ? path.slice(0, depth).join(" > ") : "Unknown"; }
// above the ontology's max depth, a suffixed cluster is named by its data-driven substate
function labelKey(sample, res, cl, depth) {
  const path = pathOf(sample, res, cl);
  if (depth > D.maxDepth) {
    const e = suffixOf(sample, res, cl);
    if (e && path.length) return path.join(" > ") + " ⊕ " + e.gene;
  }
  return labelAt(path, depth);
}
function colorFor(label) { return D.colors[label] || D.colors.Unknown; }
function resOf(S) { return S.clusters[state.res] ? state.res : Object.keys(S.clusters)[0]; }
function clusterColor(cl) { return `hsl(${(cl * 137.508) % 360} 55% 48%)`; }
function leaf(label) {
  if (label.includes(" ⊕ ")) {
    const [p, g] = label.split(" ⊕ "); const seg = p.split(" > ");
    return g + "+ " + seg[seg.length-1];
  }
  const p = label.split(" > "); return p[p.length-1];
}
function resName(r) { return r.replace(/^leiden_/, ""); }

function makeCanvas(id, title, n) {
  const div = document.createElement("div"); div.className = "panel";
  div.innerHTML = `<h3>${title} <span class="n">${n.toLocaleString()} nuclei</span></h3><canvas id="${id}"></canvas>`;
  document.getElementById("plots").appendChild(div);
  return div.querySelector("canvas");
}

function draw(canvas, xs, ys, colorAt) {
  const dpr = window.devicePixelRatio || 1, w = canvas.clientWidth || 360, h = Math.round(w * 0.82);
  canvas.width = w * dpr; canvas.height = h * dpr; canvas.style.height = h + "px";
  const ctx = canvas.getContext("2d"); ctx.scale(dpr, dpr); ctx.clearRect(0, 0, w, h);
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (let i = 0; i < xs.length; i++) { if (xs[i]<xmin)xmin=xs[i]; if (xs[i]>xmax)xmax=xs[i];
    if (ys[i]<ymin)ymin=ys[i]; if (ys[i]>ymax)ymax=ys[i]; }
  const pad = 8, sx = (w-2*pad)/(xmax-xmin), sy = (h-2*pad)/(ymax-ymin), s = Math.min(sx, sy);
  ctx.globalAlpha = 0.8;
  for (let i = 0; i < xs.length; i++) {
    ctx.fillStyle = colorAt(i);
    ctx.fillRect(pad+(xs[i]-xmin)*s, h-pad-(ys[i]-ymin)*s, 1.8, 1.8);
  }
}

function render() {
  document.getElementById("depthVal").textContent =
    state.depth > D.maxDepth ? D.maxDepth + "+" : state.depth;
  const counts = {};
  const sampleCols = ["#2F6B4F", "#A9761F"];
  for (const [name, S] of Object.entries(D.samples)) {
    const r = resOf(S);  // refined views exist at one resolution only
    const cls = S.clusters[r];
    const labels = cls.map(c => labelKey(name, r, c, state.depth));
    if (!S.integrated) labels.forEach(l => counts[l] = (counts[l]||0)+1);
    draw(els[name], S.x, S.y,
      S.sample && state.jointBy === "sample" ? i => sampleCols[S.sample[i]]
      : state.colorBy === "cluster" ? i => clusterColor(cls[i])
      : i => colorFor(labels[i]));
  }
  // legend, ordered by count
  const lg = document.getElementById("legend"); lg.innerHTML = "";
  Object.entries(counts).sort((a,b) => b[1]-a[1]).forEach(([label, n]) => {
    const s = document.createElement("span");
    s.innerHTML = `<i style="background:${colorFor(label)}"></i>${leaf(label)} <span class="cnt">${n.toLocaleString()}</span>`;
    s.title = label; lg.appendChild(s);
  });
  renderTree();
  renderTables();
}

function renderTree() {
  // nuclei per path prefix, clusters per exact path, substates per (path, gene) — current res
  const nuclei = {}, endClusters = {}, subs = {};
  for (const [name, S] of Object.entries(D.samples)) {
    if (S.integrated) continue;  // same nuclei as the samples — don't double-count
    const sizes = {};
    S.clusters[state.res].forEach(c => sizes[c] = (sizes[c]||0)+1);
    for (const [cl, n] of Object.entries(sizes)) {
      const path = pathOf(name, state.res, cl);
      if (!path.length) { nuclei["Unknown"] = (nuclei["Unknown"]||0)+n; continue; }
      for (let k = 1; k <= path.length; k++) {
        const p = path.slice(0, k).join(" > ");
        nuclei[p] = (nuclei[p]||0)+n;
      }
      const full = path.join(" > ");
      endClusters[full] = (endClusters[full]||0)+1;
      const e = suffixOf(name, state.res, cl);
      if (e) {
        const key = full + " ⊕ " + e.gene;
        subs[full] = subs[full] || {};
        subs[full][key] = subs[full][key] || { gene: e.gene, n: 0, samples: new Set(), fbeta: e.fbeta };
        subs[full][key].n += n; subs[full][key].samples.add(name);
      }
    }
  }
  const node = (color, nm, ct, cls) =>
    `<span class="node ${cls||""}"><i style="background:${color}"></i><span class="nm">${nm}</span>` +
    `<span class="ct">${ct}</span></span>`;
  function walk(sub, prefix) {
    let html = "";
    for (const [label, child] of Object.entries(sub)) {
      const p = [...prefix, label], key = p.join(" > ");
      if (!(key in nuclei)) continue;  // not observed at this resolution
      const kcl = endClusters[key] ? ` · ${endClusters[key]} cluster${endClusters[key]>1?"s":""}` : "";
      let inner = walk(child, p);
      for (const s of Object.values(subs[key] || {}).sort((a,b) => b.n - a.n)) {
        inner += `<li>${node(colorFor(key + " ⊕ " + s.gene), s.gene + "+",
          `${s.n.toLocaleString()} · ${[...s.samples].join(", ")} · Fβ ${s.fbeta}`, "sub-node")}</li>`;
      }
      html += `<li>${node(colorFor(key), label, nuclei[key].toLocaleString() + kcl)}` +
        (inner ? `<ul>${inner}</ul>` : "") + `</li>`;
    }
    return html;
  }
  const total = Object.values(D.samples).filter(S => !S.integrated).reduce((a, S) => a + S.x.length, 0);
  let unknown = nuclei["Unknown"]
    ? `<li>${node(D.colors.Unknown, "Unknown", nuclei["Unknown"].toLocaleString())}</li>` : "";
  document.getElementById("treeview").innerHTML =
    `<ul class="tree"><li>${node("var(--muted)", "cell", total.toLocaleString())}` +
    `<ul>${walk(D.tree, [])}${unknown}</ul></li></ul>`;
}

function renderTables() {
  const host = document.getElementById("tables"); host.innerHTML = "";
  for (const [name, S] of Object.entries(D.samples)) {
    const r = resOf(S);
    const cls = S.clusters[r], sizes = {};
    cls.forEach(c => sizes[c] = (sizes[c]||0)+1);
    const byCl = D.calls[name][r];
    let rows = "";
    Object.keys(byCl).sort((a,b) => sizes[b]-sizes[a]).forEach(cl => {
      const c = byCl[cl];
      const chain = c.levels.map(lv => {
        const lbl = labelAt(c.path, c.levels.indexOf(lv)+1);
        return `<span class="chip" style="background:${colorFor(lbl)};color:#fff">${lv.node}</span>` +
               `<span class="share">${(lv.share*100).toFixed(0)}%</span> `;
      }).join("&rsaquo; ") || `<span class="chip" style="background:${D.colors.Unknown};color:#fff">Unknown</span>`;
      const refusal = c.refused ? `<div class="refused">⛔ ${c.refused}</div>` : "";
      const override = c.overridden ? `<div class="override">✍ ${c.overridden}</div>` : "";
      const e = suffixOf(name, r, cl);
      const suffix = e ? ` <span class="chip suffix" title="NS-Forest markers: ${e.markers.join(" AND ")} · Fbeta=${e.fbeta}">${e.gene}+</span>` : "";
      const genes = c.levels.map(lv => {
        const panel = (D.panels && D.panels[lv.node] || []).map(g => `<span class="gene">${g}</span>`).join("");
        return panel ? `<div class="genes"><b>${lv.node}</b> (${(lv.share*100).toFixed(0)}% agree): ${panel}</div>` : "";
      }).join("");
      const details = genes ? `<details><summary>marker genes per level</summary>${genes}</details>` : "";
      rows += `<tr><td class="mono">${cl}</td><td>${sizes[cl].toLocaleString()}</td>
        <td>${chain}${suffix}${refusal}${override}${details}</td></tr>`;
    });
    const sec = document.createElement("div");
    sec.innerHTML = `<h2>${name} · resolution ${resName(r)} · ${Object.keys(byCl).length} clusters</h2>
      <div class="tablewrap"><table>
      <tr><th>cluster</th><th>n</th><th>label path · vote share per level</th></tr>${rows}</table></div>`;
    host.appendChild(sec);
  }
}

// controls
const resSeg = document.getElementById("resSeg");
RES.forEach(r => {
  const b = document.createElement("button"); b.textContent = resName(r); if (r === state.res) b.className = "on";
  b.onclick = () => { state.res = r; resSeg.querySelectorAll("button").forEach(x => x.className = "");
    b.className = "on"; render(); };
  resSeg.appendChild(b);
});
document.querySelectorAll("#jointSeg button").forEach(b => b.onclick = () => {
  state.jointBy = b.dataset.v;
  document.querySelectorAll("#jointSeg button").forEach(x => x.className = ""); b.className = "on"; render();
});
document.querySelectorAll("#colorSeg button").forEach(b => b.onclick = () => {
  state.colorBy = b.dataset.v;
  document.querySelectorAll("#colorSeg button").forEach(x => x.className = ""); b.className = "on"; render();
});
const depthEl = document.getElementById("depth");
depthEl.max = D.maxDepth + (HAS_SUFFIX ? 1 : 0);  // one stop past the treeline
depthEl.oninput = () => { state.depth = +depthEl.value; render(); };

for (const [name, S] of Object.entries(D.samples)) els[name] = makeCanvas("c_"+name.replace(/[^a-z0-9]/gi,"_"), name, S.x.length);
render();
window.addEventListener("resize", () => render());
</script>
"""


if __name__ == "__main__":
    main()
