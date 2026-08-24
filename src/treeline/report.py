"""The report verb: a static HTML view over any set of annotated h5ads.

A consumer of the other verbs, packaged for users: pass one or more annotated
`.h5ad` files (a single sample, samples side by side, an integrated joint, any mix).
The tree is the navigator — a sticky rail beside the UMAPs; hovering a label highlights
its nuclei in every panel, clicking pins the selection and filters the tables, and the
depth slider draws the treeline as a literal cut across the tree. Substate evidence
opens in a drawer; refusals and overrides are counted in a banner up top. The
resolution toggle is the union of cluster keys across the files; a panel lacking the
selected key falls back to its own first key and says so. Integrated objects are
recognized by the provenance stamp `integrate` writes (`.uns["treeline_integrate"]`),
never guessed. Embeddings are the user's job (like clustering): the report reads the
2-D basis named by --basis (default X_umap) and refuses loudly when it is missing.
Vanilla JS on canvas; no server, no frameworks.

    treeline report a_annotated.h5ad b_annotated.h5ad -o report.html [--basis X_umap]
"""

from __future__ import annotations

import json
from pathlib import Path

from treeline.annotate import annotations
from treeline.colors import palette
from treeline.harmonize import observed_paths, path_tree


def build_data(paths: list[Path], basis: str = "X_umap") -> dict:
    import scanpy as sc

    samples, calls, suffixes, anns = {}, {}, {}, []
    for f in paths:
        a = sc.read_h5ad(f)
        name = Path(f).stem.removesuffix("_annotated")
        if basis not in a.obsm:
            raise ValueError(
                f"{f}: no 2-D embedding at obsm[{basis!r}] — treeline computes no embedding "
                f"(like clustering, it is the user's job). Available obsm keys: "
                f"{list(a.obsm.keys())}. Compute one (e.g. sc.tl.umap) or pass --basis."
            )
        ann = annotations(a)
        anns.append(ann)
        keys = ann["cluster_keys"]
        um = a.obsm[basis]
        entry = {
            "x": [round(float(v), 2) for v in um[:, 0]],
            "y": [round(float(v), 2) for v in um[:, 1]],
            "clusters": {k: a.obs[k].astype(str).astype(int).tolist() for k in keys},
        }
        if "treeline_integrate" in a.uns:  # provenance stamp, not a guess
            entry["integrated"] = True
            if "sample" in a.obs.columns:
                names = sorted(a.obs["sample"].unique().tolist())
                entry["sampleNames"] = names
                entry["sample"] = [names.index(s) for s in a.obs["sample"]]
        samples[name] = entry
        calls[name] = ann["calls"]
        suffixes[name] = ann["suffixes"]
    paths_seen = observed_paths(calls)
    panels: dict = {}
    for ann in anns:
        panels.update(ann.get("panels", {}))
    return {
        "samples": samples,
        "calls": calls,
        "suffixes": suffixes,
        "panels": panels,
        "colors": palette(*anns),
        "tree": path_tree(paths_seen),
        "maxDepth": max((len(p) for p in paths_seen), default=1),
    }


def render_report(paths: list, out, title: str | None = None, basis: str = "X_umap") -> None:
    paths = [Path(p) for p in paths]
    title = title or "treeline · " + ", ".join(p.stem.removesuffix("_annotated") for p in paths)
    data = build_data(paths, basis)
    html = TEMPLATE.replace("__TITLE__", title).replace("__DATA__", json.dumps(data, separators=(",", ":")))
    Path(out).write_text(html)
    print(f"wrote {out} ({Path(out).stat().st_size / 1e6:.1f} MB)")


TEMPLATE = r"""<title>__TITLE__</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,400;12..96,600;12..96,800&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {
  --bg:#F5F7F3; --card:#FFFFFF; --surface:#ECF0E9; --ink:#1B2820; --muted:#5B6B60;
  --line:#D3DCD2; --pine:#2F6B4F; --rock:#75818C; --rock-soft:#E4E8EC;
  --lichen:#A9761F; --lichen-soft:#F1E7D2;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#0F1512; --card:#1A241E; --surface:#151D18; --ink:#E6ECE5; --muted:#93A398;
  --line:#2A3630; --pine:#5FB08A; --rock:#8B99A6; --rock-soft:#212A31;
  --lichen:#D6A34C; --lichen-soft:#322817;
}}
:root[data-theme="dark"] {
  --bg:#0F1512; --card:#1A241E; --surface:#151D18; --ink:#E6ECE5; --muted:#93A398;
  --line:#2A3630; --pine:#5FB08A; --rock:#8B99A6; --rock-soft:#212A31;
  --lichen:#D6A34C; --lichen-soft:#322817;
}
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
  font-family:"Bricolage Grotesque",system-ui,sans-serif; font-size:15px; line-height:1.5; }
.mono { font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:0.88em; }
.wrap { max-width:1360px; margin:0 auto; padding:0 20px; }
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
.seg button:focus-visible { outline:2px solid var(--pine); outline-offset:-2px; }
input[type=range] { accent-color:var(--pine); width:220px; vertical-align:middle; }
#depthVal { font-weight:800; color:var(--pine); min-width:1.2em; display:inline-block; text-align:center; }
#banner { display:none; gap:1.2rem; padding:8px 12px; margin:12px 0 0; border-radius:8px;
  background:var(--rock-soft); font-size:0.85rem; }
#banner span { cursor:pointer; }
#banner span:hover { text-decoration:underline; }
#banner .ov { color:var(--lichen); }
.layout { display:grid; grid-template-columns:260px minmax(0,1fr); gap:20px; align-items:start; margin-top:14px; }
aside#rail { position:sticky; top:58px; max-height:calc(100vh - 74px); overflow:auto;
  background:var(--card); border:1px solid var(--line); border-radius:10px; padding:10px 12px; }
#railCap { font-size:0.7rem; font-weight:600; letter-spacing:0.1em; text-transform:uppercase; color:var(--muted); }
#railCut { font-size:0.72rem; color:var(--lichen); font-weight:600; margin:2px 0 6px; }
#railHint { font-size:0.72rem; color:var(--muted); margin-bottom:6px; }
@media (max-width:900px) { .layout { grid-template-columns:1fr; } aside#rail { position:static; max-height:300px; } }
.grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(300px,1fr)); gap:16px; }
.panel { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px 14px; }
.panel h3 { margin:0 0 6px; font-size:1rem; font-weight:600; display:flex; justify-content:space-between; align-items:center; }
.panel h3 .n { color:var(--muted); font-weight:400; font-size:0.82rem; }
canvas { width:100%; height:auto; display:block; border-radius:6px; background:var(--surface); cursor:crosshair; }
table { border-collapse:collapse; width:100%; font-size:0.88rem; }
th { text-align:left; font-size:0.7rem; letter-spacing:0.08em; text-transform:uppercase; color:var(--muted);
  border-bottom:1px solid var(--line); padding:6px 10px 6px 0; }
td { border-bottom:1px solid var(--line); padding:7px 10px 7px 0; vertical-align:top; }
tr.flash td { background:var(--rock-soft); transition:background 1.2s; }
.tablewrap { overflow-x:auto; }
.chip { display:inline-block; padding:1px 8px; border-radius:999px; font-size:0.8rem; margin:1px 2px 1px 0; }
.chip.suffix { background:var(--lichen-soft); color:var(--lichen); border:1px solid var(--lichen);
  font-family:"IBM Plex Mono",monospace; font-size:0.75rem; cursor:pointer; }
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
/* navigator tree */
.tree, .tree ul { list-style:none; margin:0; padding-left:0; }
.tree ul { padding-left:0.9rem; border-left:1px solid var(--line); margin-left:0.34rem; }
.tree li { padding:1px 0; }
.tree .node { display:inline-flex; align-items:center; gap:0.4rem; cursor:pointer;
  border-radius:5px; padding:1px 5px; font-size:0.82rem; }
.tree .node:hover { background:var(--surface); }
.tree .node.pin { outline:2px solid var(--pine); }
.tree .node i { width:10px; height:10px; border-radius:3px; flex:none; }
.tree .node .nm { font-weight:600; }
.tree .node .ct { color:var(--muted); font-size:0.72rem; font-variant-numeric:tabular-nums; }
.tree .sub-node .nm { font-family:"IBM Plex Mono",monospace; font-weight:500; font-size:0.76rem; }
.tree .sub-node { border:1px solid var(--lichen); border-radius:6px; padding:0 6px; margin:1px 0; }
/* tooltip + drawer */
#tip { position:fixed; display:none; z-index:20; pointer-events:none; background:var(--card);
  border:1px solid var(--line); border-radius:8px; padding:6px 10px; font-size:0.8rem;
  box-shadow:0 4px 14px rgba(0,0,0,0.18); max-width:340px; }
#tip .mono { font-size:0.75rem; color:var(--muted); }
#drawer { position:fixed; right:16px; bottom:16px; z-index:30; display:none; width:min(400px,92vw);
  background:var(--card); border:1px solid var(--lichen); border-radius:12px; padding:14px 16px;
  box-shadow:0 8px 28px rgba(0,0,0,0.25); font-size:0.86rem; }
#drawer h4 { margin:0 0 4px; font-size:0.95rem; }
#drawer .close { float:right; cursor:pointer; border:0; background:none; color:var(--muted);
  font:inherit; font-size:1.1rem; }
#drawer table { font-size:0.8rem; }
#drawer .metrics { color:var(--muted); margin:2px 0 8px; }
.mix { height:12px; border-radius:4px; display:flex; overflow:hidden; margin:6px 0 2px; }
.mixlbl { font-size:0.75rem; color:var(--muted); }
</style>
<div class="wrap">
  <h1>__TITLE__</h1>
  <p class="sub">Multi-level annotation, evidence-gated. Hover the tree to light up nuclei;
  click to pin and filter. Labels stop where their vote does.</p>
  <div class="controls">
    <span><span class="ctl-label">Resolution</span>&nbsp; <span class="seg" id="resSeg"></span></span>
    <span><span class="ctl-label">Depth</span>&nbsp; coarse <input type="range" id="depth" min="1" value="2"> fine
      <span id="depthVal"></span></span>
    <span><span class="ctl-label">Color by</span>&nbsp; <span class="seg" id="colorSeg">
      <button data-v="label" class="on">labels</button><button data-v="cluster">clusters</button></span></span>
    <span id="jointCtl" style="display:none"><span class="ctl-label">Integrated colors</span>&nbsp; <span class="seg" id="jointSeg">
      <button data-v="label" class="on">labels</button><button data-v="sample">sample</button></span></span>
  </div>
  <div id="banner"></div>
  <div class="layout">
    <aside id="rail">
      <div id="railCap">The tree</div>
      <div id="railHint">hover: highlight · click: pin/unpin · esc: clear</div>
      <div id="treeview"></div>
    </aside>
    <main>
      <div class="grid" id="plots"></div>
      <div id="tables"></div>
    </main>
  </div>
  <footer>Per-sample views are annotated independently; an integrated view (recognized by its
  provenance stamp) is scANVI on the tree-cut consensus prior, reclustered by the user and
  re-annotated by the same vote. Gate refusals are automated, self-reporting rules; ✍ overrides
  are signed (SPECS provenance tiers). Generated by treeline report.</footer>
</div>
<div id="tip"></div>
<div id="drawer"></div>
<script>
const D = __DATA__;
const RES = [...new Set(Object.values(D.samples).flatMap(S => Object.keys(S.clusters)))];
const HAS_SUFFIX = Object.values(D.suffixes || {}).some(byRes => Object.values(byRes).some(m => Object.keys(m).length));
let state = { res: RES[0], depth: 2, jointBy: "label", colorBy: "label", hover: null, sel: null };
const els = {};

function pathOf(sample, res, cl) { const c = D.calls[sample][res][String(cl)]; return c ? c.path : []; }
function suffixOf(sample, res, cl) {
  return ((D.suffixes[sample] || {})[res] || {})[String(cl)] || null; }
function labelAt(path, depth) { return path.length ? path.slice(0, depth).join(" > ") : "Unknown"; }
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
const COMBO = {};
for (const [s, byRes] of Object.entries(D.suffixes)) for (const [r, byCl] of Object.entries(byRes))
  for (const [cl, e] of Object.entries(byCl)) {
    const p = ((D.calls[s][r] || {})[cl] || {}).path || [];
    if (p.length) COMBO[p.join(" > ") + " ⊕ " + e.gene] = e.markers.join("+");
  }
function leaf(label) {
  if (label.includes(" ⊕ ")) {
    const [p, g] = label.split(" ⊕ "); const seg = p.split(" > ");
    return (COMBO[label] || g) + "+ " + seg[seg.length-1];
  }
  const p = label.split(" > "); return p[p.length-1];
}
function resName(r) { return r.replace(/^leiden_/, ""); }
const COUNT_ALL = Object.values(D.samples).every(S => S.integrated);
function counted(S) { return COUNT_ALL || !S.integrated; }

// ---- shared selection: {type:"node", path:[...]} | {type:"substate", key} | {type:"unknown"}
function matchesSel(sel, name, r, cl) {
  if (!sel) return true;
  const path = pathOf(name, r, cl);
  if (sel.type === "unknown") return path.length === 0;
  if (sel.type === "node") return sel.path.every((s, i) => path[i] === s);
  if (sel.type === "substate") {
    const e = suffixOf(name, r, cl);
    return !!e && path.length && (path.join(" > ") + " ⊕ " + e.gene) === sel.key;
  }
  return true;
}
function sameSel(a, b) { return JSON.stringify(a) === JSON.stringify(b); }
function activeSel() { return state.sel || state.hover; }

function makeCanvas(id, title) {
  const div = document.createElement("div"); div.className = "panel";
  div.innerHTML = `<h3>${title} <span class="n"></span></h3><canvas id="${id}"></canvas>`;
  document.getElementById("plots").appendChild(div);
  return div.querySelector("canvas");
}

function draw(canvas, xs, ys, colorAt, matchAt) {
  const dpr = window.devicePixelRatio || 1, w = canvas.clientWidth || 360, h = Math.round(w * 0.82);
  canvas.width = w * dpr; canvas.height = h * dpr; canvas.style.height = h + "px";
  const ctx = canvas.getContext("2d"); ctx.scale(dpr, dpr); ctx.clearRect(0, 0, w, h);
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (let i = 0; i < xs.length; i++) { if (xs[i]<xmin)xmin=xs[i]; if (xs[i]>xmax)xmax=xs[i];
    if (ys[i]<ymin)ymin=ys[i]; if (ys[i]>ymax)ymax=ys[i]; }
  const pad = 8, s = Math.min((w-2*pad)/(xmax-xmin), (h-2*pad)/(ymax-ymin));
  const px = new Float32Array(xs.length), py = new Float32Array(xs.length);
  for (let i = 0; i < xs.length; i++) { px[i] = pad+(xs[i]-xmin)*s; py[i] = h-pad-(ys[i]-ymin)*s; }
  if (matchAt) {
    ctx.globalAlpha = 0.06;
    for (let i = 0; i < xs.length; i++) if (!matchAt(i)) { ctx.fillStyle = colorAt(i); ctx.fillRect(px[i], py[i], 1.8, 1.8); }
    ctx.globalAlpha = 0.9;
    for (let i = 0; i < xs.length; i++) if (matchAt(i)) { ctx.fillStyle = colorAt(i); ctx.fillRect(px[i], py[i], 1.8, 1.8); }
  } else {
    ctx.globalAlpha = 0.8;
    for (let i = 0; i < xs.length; i++) { ctx.fillStyle = colorAt(i); ctx.fillRect(px[i], py[i], 1.8, 1.8); }
  }
  canvas._px = px; canvas._py = py;
}

function renderPlots() {
  const sel = activeSel();
  const sampleCols = ["#2F6B4F", "#A9761F", "#75818C", "#5FB08A"];
  for (const [name, S] of Object.entries(D.samples)) {
    const r = resOf(S), cls = S.clusters[r];
    const labels = cls.map(c => labelKey(name, r, c, state.depth));
    const nspan = els[name].parentElement.querySelector(".n");
    nspan.textContent = `${S.x.length.toLocaleString()} nuclei · ${resName(r)}` + (r !== state.res ? " (only)" : "");
    const matchByCl = {};
    if (sel) for (const c of new Set(cls)) matchByCl[c] = matchesSel(sel, name, r, c);
    draw(els[name], S.x, S.y,
      S.sample && state.jointBy === "sample" ? i => sampleCols[S.sample[i]]
      : state.colorBy === "cluster" ? i => clusterColor(cls[i])
      : i => colorFor(labels[i]),
      sel ? i => matchByCl[cls[i]] : null);
  }
}

function render() {
  document.getElementById("depthVal").textContent =
    state.depth > D.maxDepth ? D.maxDepth + "+" : state.depth;
  renderPlots();
  renderBanner();
  renderTree();
  renderTables();
}

function renderBanner() {
  let refs = 0, ovs = 0;
  for (const [name, S] of Object.entries(D.samples)) {
    const r = resOf(S);
    for (const c of Object.values(D.calls[name][r])) { if (c.refused) refs++; if (c.overridden) ovs++; }
  }
  const b = document.getElementById("banner");
  if (!refs && !ovs) { b.style.display = "none"; return; }
  b.style.display = "flex";
  b.innerHTML = (refs ? `<span id="jumpRef">⛔ ${refs} gate refusal${refs>1?"s":""}</span>` : "") +
    (ovs ? `<span class="ov" id="jumpOv">✍ ${ovs} signed override${ovs>1?"s":""}</span>` : "");
  const jump = cls => { const el = document.querySelector("#tables ." + cls);
    if (el) { el.scrollIntoView({behavior:"smooth", block:"center"}); } };
  const jr = document.getElementById("jumpRef"); if (jr) jr.onclick = () => jump("refused");
  const jo = document.getElementById("jumpOv"); if (jo) jo.onclick = () => jump("override");
}

function renderTree() {
  const nuclei = {}, endClusters = {}, subs = {};
  for (const [name, S] of Object.entries(D.samples)) {
    if (!counted(S)) continue;
    const r = resOf(S), sizes = {};
    S.clusters[r].forEach(c => sizes[c] = (sizes[c]||0)+1);
    for (const [cl, n] of Object.entries(sizes)) {
      const path = pathOf(name, r, cl);
      if (!path.length) { nuclei["Unknown"] = (nuclei["Unknown"]||0)+n; continue; }
      for (let k = 1; k <= path.length; k++) {
        const p = path.slice(0, k).join(" > ");
        nuclei[p] = (nuclei[p]||0)+n;
      }
      const full = path.join(" > ");
      endClusters[full] = (endClusters[full]||0)+1;
    }
  }
  // substates come from every view (they may live only on the integrated object)
  for (const [name, S] of Object.entries(D.samples)) {
    const r = resOf(S), sizes = {};
    S.clusters[r].forEach(c => sizes[c] = (sizes[c]||0)+1);
    for (const [cl, e] of Object.entries((D.suffixes[name] || {})[r] || {})) {
      const path = pathOf(name, r, cl);
      if (!path.length) continue;
      const full = path.join(" > ");
      if (!(full in nuclei)) continue;
      const key = full + " ⊕ " + e.gene, n = sizes[cl] || 0;
      subs[full] = subs[full] || {};
      const cur = subs[full][key] || { name: e.markers.join("+"), n: 0, fbeta: e.fbeta, ref: null, refN: 0 };
      cur.n += n;
      if (n > cur.refN) { cur.ref = [name, r, String(cl)]; cur.refN = n; }
      subs[full][key] = cur;
    }
  }
  const pin = state.sel;
  // the tree grows and shrinks with the depth slider: nodes deeper than the current
  // depth stay collapsed; substates unfold only at the above-the-treeline stop
  const maxShow = Math.min(state.depth, D.maxDepth);
  const showSubs = state.depth > D.maxDepth;
  const node = (color, nm, ct, cls, attrs) =>
    `<span class="node ${cls||""}" ${attrs||""}><i style="background:${color}"></i><span class="nm">${nm}</span>` +
    `<span class="ct">${ct}</span></span>`;
  function walk(sub, prefix) {
    let html = "";
    for (const [label, child] of Object.entries(sub)) {
      const p = [...prefix, label], key = p.join(" > ");
      if (!(key in nuclei) || p.length > maxShow) continue;
      const kcl = endClusters[key] ? ` · ${endClusters[key]}cl` : "";
      let cls = pin && pin.type === "node" && pin.path.join(" > ") === key ? "pin" : "";
      let inner = walk(child, p);
      if (showSubs) for (const [skey, s] of Object.entries(subs[key] || {}).sort((a,b) => b[1].n - a[1].n)) {
        let scls = "sub-node";
        if (pin && pin.type === "substate" && pin.key === skey) scls += " pin";
        inner += `<li>${node(colorFor(skey), s.name + "+", `${s.n.toLocaleString()} · Fβ ${s.fbeta}`,
          scls, `data-skey="${encodeURIComponent(skey)}" data-ref="${encodeURIComponent(s.ref.join("|"))}"`)}</li>`;
      }
      html += `<li>${node(colorFor(key), label, nuclei[key].toLocaleString() + kcl, cls,
        `data-path="${encodeURIComponent(key)}"`)}` + (inner ? `<ul>${inner}</ul>` : "") + `</li>`;
    }
    return html;
  }
  const total = Object.values(D.samples).filter(S => counted(S)).reduce((a, S) => a + S.x.length, 0);
  let unknown = nuclei["Unknown"]
    ? `<li>${node(D.colors.Unknown, "Unknown", nuclei["Unknown"].toLocaleString(),
        pin && pin.type === "unknown" ? "pin" : "", 'data-unknown="1"')}</li>` : "";
  document.getElementById("treeview").innerHTML =
    `<ul class="tree"><li>${node("var(--muted)", "cell", total.toLocaleString())}` +
    `<ul>${walk(D.tree, [])}${unknown}</ul></li></ul>`;
}

function renderTables() {
  const host = document.getElementById("tables"); host.innerHTML = "";
  const pin = state.sel;
  for (const [name, S] of Object.entries(D.samples)) {
    const r = resOf(S);
    const cls = S.clusters[r], sizes = {};
    cls.forEach(c => sizes[c] = (sizes[c]||0)+1);
    const byCl = D.calls[name][r];
    let rows = "", shown = 0;
    Object.keys(byCl).sort((a,b) => sizes[b]-sizes[a]).forEach(cl => {
      if (pin && !matchesSel(pin, name, r, cl)) return;
      shown++;
      const c = byCl[cl];
      const chain = c.levels.map(lv => {
        const lbl = labelAt(c.path, c.levels.indexOf(lv)+1);
        return `<span class="chip" style="background:${colorFor(lbl)};color:#fff">${lv.node}</span>` +
               `<span class="share">${(lv.share*100).toFixed(0)}%</span> `;
      }).join("&rsaquo; ") || `<span class="chip" style="background:${D.colors.Unknown};color:#fff">Unknown</span>`;
      const refusal = c.refused ? `<div class="refused">⛔ ${c.refused}</div>` : "";
      const override = c.overridden ? `<div class="override">✍ ${c.overridden}</div>` : "";
      const e = suffixOf(name, r, cl);
      const suffix = e ? ` <span class="chip suffix" data-dref="${encodeURIComponent([name, r, cl].join("|"))}">${e.markers.join("+")}+</span>` : "";
      const genes = c.levels.map(lv => {
        const panel = (D.panels && D.panels[lv.node] || []).map(g => `<span class="gene">${g}</span>`).join("");
        return panel ? `<div class="genes"><b>${lv.node}</b> (${(lv.share*100).toFixed(0)}% agree): ${panel}</div>` : "";
      }).join("");
      const details = genes ? `<details><summary>marker genes per level</summary>${genes}</details>` : "";
      rows += `<tr id="row-${cssId(name)}-${cl}"><td class="mono">${cl}</td><td>${sizes[cl].toLocaleString()}</td>
        <td>${chain}${suffix}${refusal}${override}${details}</td></tr>`;
    });
    const filt = pin ? ` · ${shown} shown (pinned)` : "";
    const sec = document.createElement("div");
    sec.innerHTML = `<h2>${name} · resolution ${resName(r)} · ${Object.keys(byCl).length} clusters${filt}</h2>
      <div class="tablewrap"><table>
      <tr><th>cluster</th><th>n</th><th>label path · vote share per level</th></tr>${rows}</table></div>`;
    host.appendChild(sec);
  }
}

function cssId(s) { return s.replace(/[^a-z0-9]/gi, "_"); }

// ---- tree interactions (delegated)
const treeview = document.getElementById("treeview");
function selFromEl(el) {
  if (el.dataset.unknown) return {type:"unknown"};
  if (el.dataset.skey) return {type:"substate", key: decodeURIComponent(el.dataset.skey)};
  if (el.dataset.path) return {type:"node", path: decodeURIComponent(el.dataset.path).split(" > ")};
  return null;
}
treeview.addEventListener("mouseover", ev => {
  const el = ev.target.closest(".node"); const sel = el && selFromEl(el);
  if (!sameSel(state.hover, sel)) { state.hover = sel; if (!state.sel) renderPlots(); }
});
treeview.addEventListener("mouseleave", () => { state.hover = null; if (!state.sel) renderPlots(); });
treeview.addEventListener("click", ev => {
  const el = ev.target.closest(".node"); if (!el) return;
  const sel = selFromEl(el); if (!sel) return;
  if (sel.type === "substate" && el.dataset.ref) openDrawer(...decodeURIComponent(el.dataset.ref).split("|"));
  state.sel = sameSel(state.sel, sel) ? null : sel;
  render();
});
document.addEventListener("keydown", ev => {
  if (ev.key === "Escape") { state.sel = null; state.hover = null; closeDrawer(); render(); }
});

// ---- hover identity on the plots
const tip = document.getElementById("tip");
function attachTip(name, canvas) {
  canvas.addEventListener("mousemove", ev => {
    const S = D.samples[name], r = resOf(S);
    const rect = canvas.getBoundingClientRect();
    const mx = ev.clientX - rect.left, my = ev.clientY - rect.top;
    const px = canvas._px, py = canvas._py;
    if (!px) return;
    let best = -1, bd = 64; // 8px radius
    for (let i = 0; i < px.length; i++) {
      const d = (px[i]-mx)*(px[i]-mx) + (py[i]-my)*(py[i]-my);
      if (d < bd) { bd = d; best = i; }
    }
    if (best < 0) { tip.style.display = "none"; canvas._hit = null; return; }
    const cl = S.clusters[r][best];
    canvas._hit = cl;
    const path = pathOf(name, r, cl); const e = suffixOf(name, r, cl);
    const n = S.clusters[r].reduce((a, c) => a + (c === cl), 0);
    tip.innerHTML = `<b>cluster ${cl}</b> · ${n.toLocaleString()} nuclei<br>` +
      `<span class="mono">${path.join(" › ") || "Unknown"}${e ? " · " + e.markers.join("+") + "+" : ""}</span>`;
    tip.style.display = "block";
    tip.style.left = Math.min(ev.clientX + 14, innerWidth - 360) + "px";
    tip.style.top = (ev.clientY + 14) + "px";
  });
  canvas.addEventListener("mouseleave", () => { tip.style.display = "none"; canvas._hit = null; });
  canvas.addEventListener("click", () => {
    if (canvas._hit == null) return;
    const row = document.getElementById(`row-${cssId(name)}-${canvas._hit}`);
    if (row) { row.scrollIntoView({behavior:"smooth", block:"center"});
      row.classList.add("flash"); setTimeout(() => row.classList.remove("flash"), 1200); }
  });
}

// ---- evidence drawer
const drawer = document.getElementById("drawer");
function closeDrawer() { drawer.style.display = "none"; }
function openDrawer(name, r, cl) {
  const e = suffixOf(name, r, cl); if (!e) return;
  const S = D.samples[name];
  const path = pathOf(name, r, cl);
  let mix = "";
  if (S.sample) {
    const cls = S.clusters[r]; const counts = {};
    let tot = 0;
    for (let i = 0; i < cls.length; i++) if (String(cls[i]) === String(cl)) { counts[S.sample[i]] = (counts[S.sample[i]]||0)+1; tot++; }
    const cols = ["#2F6B4F", "#A9761F", "#75818C", "#5FB08A"];
    mix = `<div class="mix">` + S.sampleNames.map((sn, i) =>
        `<span style="flex:${counts[i]||0};background:${cols[i]}"></span>`).join("") + `</div>` +
      `<div class="mixlbl">${S.sampleNames.map((sn, i) => `${sn} ${(100*(counts[i]||0)/tot).toFixed(0)}%`).join(" · ")}</div>`;
  }
  drawer.innerHTML = `<button class="close" aria-label="close">×</button>` +
    `<h4>${e.markers.join("+")}+ ${path[path.length-1] || ""}</h4>` +
    `<div class="metrics">${name} · cluster ${cl} · Fβ ${e.fbeta} · PPV ${e.ppv} · recall ${e.recall}<br>` +
    `necessary + sufficient vs siblings under the same label</div>` +
    `<table><tr><th>marker</th><th>binary</th><th>threshold</th><th>on-target</th></tr>` +
    e.marker_stats.map(m => `<tr><td class="mono">${m.gene}</td><td>${m.binary_score}</td>` +
      `<td>${m.threshold}</td><td>${m.on_target}</td></tr>`).join("") + `</table>` + mix;
  drawer.style.display = "block";
  drawer.querySelector(".close").onclick = closeDrawer;
}
document.getElementById("tables").addEventListener("click", ev => {
  const chip = ev.target.closest(".chip.suffix");
  if (chip) openDrawer(...decodeURIComponent(chip.dataset.dref).split("|"));
});

// ---- controls
const resSeg = document.getElementById("resSeg");
RES.forEach(r => {
  const b = document.createElement("button"); b.textContent = resName(r); if (r === state.res) b.className = "on";
  b.onclick = () => { state.res = r; resSeg.querySelectorAll("button").forEach(x => x.className = "");
    b.className = "on"; render(); };
  resSeg.appendChild(b);
});
if (Object.values(D.samples).some(S => S.sample)) document.getElementById("jointCtl").style.display = "";
document.querySelectorAll("#jointSeg button").forEach(b => b.onclick = () => {
  state.jointBy = b.dataset.v;
  document.querySelectorAll("#jointSeg button").forEach(x => x.className = ""); b.className = "on"; render();
});
document.querySelectorAll("#colorSeg button").forEach(b => b.onclick = () => {
  state.colorBy = b.dataset.v;
  document.querySelectorAll("#colorSeg button").forEach(x => x.className = ""); b.className = "on"; render();
});
const depthEl = document.getElementById("depth");
depthEl.max = D.maxDepth + (HAS_SUFFIX ? 1 : 0);
depthEl.oninput = () => { state.depth = +depthEl.value; render(); };

for (const [name, S] of Object.entries(D.samples)) {
  els[name] = makeCanvas("c_" + cssId(name), name);
  attachTip(name, els[name]);
}
render();
window.addEventListener("resize", () => renderPlots());
</script>
"""
