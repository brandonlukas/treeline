"""The annotate verb: multi-level annotation written back into the AnnData.

Input contract: `adata.X` log-normalized, cluster labels in the named `.obs` columns
(any clustering method, one or more columns = "resolutions"). Output: per-cluster calls,
shares, refusals and NS-Forest suffixes in `.uns["treeline"]` (one JSON string — h5ad
cannot serialize nested lists of dicts), and a per-nucleus label column
`treeline_<cluster_key>` per clustering. An annotated `.h5ad` is self-contained: it is
the input to `integrate` and to `colors`.
"""

from __future__ import annotations

import dataclasses
import json

from treeline.nsforest import suffixes_for
from treeline.vote import (
    DESCEND_AGREE,
    GATE_OVERLAP,
    GATE_SCORE_R,
    N_MARKERS,
    assign_all,
    load_panels,
    score_nodes,
)

UNS_KEY = "treeline"


def _check_lognorm(adata) -> None:
    """Refuse raw counts loudly: score_genes on counts produces garbage scores and
    therefore confidently wrong labels — the exact failure treeline exists to prevent."""
    mx = adata.X.max()
    if mx > 50:
        raise ValueError(
            f"adata.X max is {mx:.0f} — this looks like raw counts, not log-normalized "
            "expression. annotate scores genes on log-normalized X: keep counts in "
            "layers['counts'] (integrate needs them), then sc.pp.normalize_total + sc.pp.log1p."
        )


def annotate(
    adata,
    dag: dict,
    gene_sets_csv,
    cluster_keys: list[str],
    overrides: list[dict] | None = None,
    *,
    n_markers: int = N_MARKERS,
    descend_agree: float = DESCEND_AGREE,
    gate_overlap: float = GATE_OVERLAP,
    gate_score_r: float = GATE_SCORE_R,
):
    """Annotate every clustering in `cluster_keys` (CL scoring + vote only — cheap, so
    multiple resolutions cost little and all feed the integrate prior); mutates and
    returns `adata`. NS-Forest substate naming is the separate, costly `add_substates`
    step — run it once, on the clustering you settle on."""
    missing = [k for k in cluster_keys if k not in adata.obs.columns]
    if missing:
        raise ValueError(f"cluster keys not in .obs: {missing}; available columns: {list(adata.obs.columns)}")
    _check_lognorm(adata)
    panels = load_panels(gene_sets_csv, dag, n_markers=n_markers)
    score_nodes(adata, panels)
    calls: dict[str, dict] = {}
    for key in cluster_keys:
        cs = assign_all(
            adata.obs, dag, panels, key, overrides,
            descend_agree=descend_agree, gate_overlap=gate_overlap, gate_score_r=gate_score_r,
        )
        calls[key] = {cl: dataclasses.asdict(c) for cl, c in cs.items()}
        label_of = {cl: " > ".join(d["path"]) or "Unknown" for cl, d in calls[key].items()}
        adata.obs[f"treeline_{key}"] = adata.obs[key].astype(str).map(label_of)
    adata.uns[UNS_KEY] = json.dumps(
        {
            "calls": calls,
            "suffixes": {},
            "panels": panels,
            "cluster_keys": list(cluster_keys),
            "dag_fetched": dag.get("fetched"),
            "overrides": overrides or [],
        }
    )
    return adata


def add_substates(adata, cluster_keys: list[str], **nsforest_kwargs):
    """The substates verb: NS-Forest `{gene}+` naming for clusters that collapse under
    one label — on the named clusterings only (costly; typically run once, at the final
    resolution). Requires a previously annotated adata; mutates and returns it.
    `nsforest_kwargs` pass through to `nsforest.suffixes_for`."""
    if UNS_KEY not in adata.uns:
        raise ValueError("not annotated (no .uns['treeline']) — run annotate first")
    ann = json.loads(adata.uns[UNS_KEY])
    missing = [k for k in cluster_keys if k not in ann["calls"]]
    if missing:
        raise ValueError(f"clusterings not annotated: {missing}; annotated: {list(ann['calls'])}")
    for key in cluster_keys:
        clusters = adata.obs[key].astype(str)
        ann["suffixes"][key] = suffixes_for(adata, clusters, ann["calls"][key], **nsforest_kwargs)
    adata.uns[UNS_KEY] = json.dumps(ann)
    return adata


def annotations(adata) -> dict:
    """Decode `.uns['treeline']` from an annotated AnnData."""
    return json.loads(adata.uns[UNS_KEY])
