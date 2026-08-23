"""The hierarchical vote: subtree-max descent of the CL DAG, evidence-gated.

Per nucleus, every scored node's gene set yields a score (scanpy score_genes). Descent
is per-node: at the current node, nuclei vote among its children by subtree-max (a
child's score is the best score anywhere in its descendant set — siblings pool rather
than split). DAG semantics per SPECS: a nucleus whose max lies in two siblings' shared
descendants votes for both (shares need not sum to 1); descent goes to the single
highest-share child, only while the share clears DESCEND_AGREE. The root is gated like
any node — a cluster that can't clear it at the first descent is Unknown.

The discriminability gate refuses a descent whose top two children cannot be told apart
on their *exclusive* descendants — BOTH criteria must trip (gene overlap AND score
correlation): measured on 1619, overlap alone wrongly refused the uterine-SMC descent
(5/10 shared genes, but scores discriminate, r~0.6). Shared descendants make correlation
structural, not an evidence failure; a child contained in its sibling is exempt.

Expert overrides (provenance tier 3) are the third influence: signed
{node, decision, justification, author, date} entries; `stop` caps descent at a node,
`prune` removes a node and its exclusive descendants from the candidate space. The
motivating case: the lymphatic-EC panel is generic endothelial genes and wins the EC
descent decisively (share 0.84) — no automated score-side rule can refuse it; knowing
the panel lacks PROX1/CCL21 is domain knowledge, and the contract says domain knowledge
gets signed, not smuggled into thresholds.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from treeline.tree import descendants

DESCEND_AGREE = 0.5  # vote share a child needs before the label descends into it
GATE_OVERLAP = 0.5  # exclusive-panel gene overlap (of the smaller panel) at/above which descent is refused
GATE_SCORE_R = 0.9  # per-nucleus exclusive-score correlation at/above which descent is refused
N_MARKERS = 10  # top genes taken per gene set


def load_panels(csv_path, dag: dict, n_markers: int = N_MARKERS) -> dict[str, list[str]]:
    """Node label -> marker symbols; a node's sets merge (union, order kept)."""
    df = pd.read_csv(csv_path)
    missing = {"Gene Set Name", "Gene Symbol"} - set(df.columns)
    if missing:
        raise ValueError(
            f"gene-set table {csv_path} is missing column(s) {sorted(missing)}; found {list(df.columns)}. "
            "Expected a CellGuide-style CSV: 'Gene Set Name' (named on each set's first row, "
            "blank rows fill down) and 'Gene Symbol' (one gene per row)."
        )
    df["Gene Set Name"] = df["Gene Set Name"].ffill().str.removesuffix(" - marker genes")
    by_set = {name: g["Gene Symbol"].head(n_markers).tolist() for name, g in df.groupby("Gene Set Name", sort=False)}
    panels = {}
    for label, node in dag["nodes"].items():
        genes = list(dict.fromkeys(g for s in node["sets"] for g in by_set[s]))
        if genes:
            panels[label] = genes
    return panels


def score_nodes(adata, panels: dict[str, list[str]]) -> list[str]:
    """Score every panel with >=3 present genes into obs score_{label}; returns labels scored."""
    import scanpy as sc

    scored = []
    for label, panel in panels.items():
        genes = [g for g in panel if g in adata.var_names]
        if len(genes) >= 3:
            sc.tl.score_genes(adata, genes, score_name=f"score_{label}")
            scored.append(label)
    return scored


@dataclass
class Level:
    node: str
    share: float
    runner: str | None = None
    runner_share: float | None = None
    n_nuclei: int = 0


@dataclass
class Call:
    path: list[str] = field(default_factory=list)
    levels: list[Level] = field(default_factory=list)
    refused: str | None = None  # gate refusal message, if descent stopped on the gate
    overridden: str | None = None  # signed expert override that stopped/shaped descent

    @property
    def final(self) -> str:
        return self.path[-1] if self.path else "Unknown"


def load_overrides(path) -> list[dict]:
    """Overrides config: schema-enforced, every entry signed and justified."""
    import json
    from pathlib import Path

    if not Path(path).exists():
        return []
    entries = json.loads(Path(path).read_text())
    for e in entries:
        missing = {"node", "decision", "justification", "author", "date"} - e.keys()
        if missing:
            raise ValueError(f"override entry {e.get('node', '?')!r} missing {sorted(missing)}")
        if e["decision"] not in ("stop", "prune"):  # "relabel" joins when implemented
            raise ValueError(f"override decision {e['decision']!r} not in the closed vocabulary")
    return entries


def assign_cluster(
    dag: dict,
    panels: dict[str, list[str]],
    scores: pd.DataFrame,
    overrides: list[dict] | None = None,
    *,
    descend_agree: float = DESCEND_AGREE,
    gate_overlap: float = GATE_OVERLAP,
    gate_score_r: float = GATE_SCORE_R,
) -> Call:
    """Descend the DAG for one cluster. `scores` is nuclei x score_{label} columns."""
    overrides = overrides or []
    stops = {o["node"]: o for o in overrides if o["decision"] == "stop"}
    pruned = {o["node"] for o in overrides if o["decision"] == "prune"}
    scored = {c.removeprefix("score_") for c in scores.columns} - pruned
    call = Call()
    children = dag["roots"]
    # a set-less single root ("cell") carries no vote information: enter it silently
    while len(children) == 1 and not dag["nodes"][children[0]]["sets"]:
        children = dag["nodes"][children[0]]["children"]
    while children:
        children = [c for c in children if c not in pruned]
        cand = {c: descendants(dag, c) & scored for c in children}
        cand = {c: d for c, d in cand.items() if d}
        if not cand:
            break
        submax = pd.DataFrame({c: scores[[f"score_{n}" for n in d]].max(axis=1) for c, d in cand.items()})
        votes = submax.eq(submax.max(axis=1), axis=0)  # ties/shared descendants vote for both
        shares = votes.mean().sort_values(ascending=False)
        win, share = str(shares.index[0]), float(shares.iloc[0])
        runner = str(shares.index[1]) if len(shares) > 1 else None

        if runner is not None:
            gate = _gate(dag, panels, cand[win], cand[runner], scores, gate_overlap, gate_score_r)
            if gate:
                parent = call.path[-1] if call.path else "root"
                call.refused = f"descent from '{parent}' refused: '{win}' vs '{runner}' {gate}"
                break
        if share < descend_agree:
            break
        call.path.append(win)
        call.levels.append(
            Level(
                node=win,
                share=round(share, 3),
                runner=runner,
                runner_share=round(float(shares.iloc[1]), 3) if runner else None,
                n_nuclei=int(votes[win].sum()),
            )
        )
        if win in stops:
            o = stops[win]
            call.overridden = f"descent stopped at '{win}' by override ({o['author']}, {o['date']}): {o['justification']}"
            break
        scores = scores.loc[votes[win]]
        children = dag["nodes"][win]["children"]
    return call


def _gate(
    dag, panels, desc_a: set[str], desc_b: set[str], scores: pd.DataFrame,
    gate_overlap: float = GATE_OVERLAP, gate_score_r: float = GATE_SCORE_R,
) -> str | None:
    """Refusal reason if the two sibling subtrees are indiscriminable on their exclusive
    descendants; None if the descent may be offered. Containment is exempt (structural)."""
    excl_a, excl_b = desc_a - desc_b, desc_b - desc_a
    if not excl_a or not excl_b:
        return None
    genes_a = {g for n in excl_a for g in panels.get(n, [])}
    genes_b = {g for n in excl_b for g in panels.get(n, [])}
    if not genes_a or not genes_b:
        return None
    shared, smaller = len(genes_a & genes_b), min(len(genes_a), len(genes_b))
    if shared / smaller < gate_overlap:
        return None
    sub_a = scores[[f"score_{n}" for n in excl_a if n in {c.removeprefix('score_') for c in scores.columns}]]
    sub_b = scores[[f"score_{n}" for n in excl_b if n in {c.removeprefix('score_') for c in scores.columns}]]
    if sub_a.empty or sub_b.empty:
        return None
    with np.errstate(invalid="ignore"):
        r = float(np.corrcoef(sub_a.max(axis=1), sub_b.max(axis=1))[0, 1])
    if np.isfinite(r) and r >= gate_score_r:  # both criteria must trip
        return f"siblings share {shared}/{smaller} exclusive genes, score r={r:.2f}"
    return None


def assign_all(
    obs: pd.DataFrame,
    dag: dict,
    panels: dict[str, list[str]],
    cluster_key: str,
    overrides: list[dict] | None = None,
    *,
    descend_agree: float = DESCEND_AGREE,
    gate_overlap: float = GATE_OVERLAP,
    gate_score_r: float = GATE_SCORE_R,
) -> dict[str, Call]:
    score_cols = [c for c in obs.columns if c.startswith("score_")]
    return {
        str(cl): assign_cluster(
            dag, panels, grp[score_cols], overrides,
            descend_agree=descend_agree, gate_overlap=gate_overlap, gate_score_r=gate_score_r,
        )
        for cl, grp in obs.groupby(cluster_key, observed=True)
    }
