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


def annotate(
    adata,
    dag: dict,
    gene_sets_csv,
    cluster_keys: list[str],
    overrides: list[dict] | None = None,
    suffixes: bool = True,
    *,
    n_markers: int = N_MARKERS,
    descend_agree: float = DESCEND_AGREE,
    gate_overlap: float = GATE_OVERLAP,
    gate_score_r: float = GATE_SCORE_R,
):
    """Annotate every clustering in `cluster_keys`; mutates and returns `adata`.

    The vote's knobs are keyword parameters with the loud defaults; NS-Forest's own
    knobs (trees, combination size, beta, ...) are on `nsforest.suffixes_for` — call it
    directly if you need them."""
    panels = load_panels(gene_sets_csv, dag, n_markers=n_markers)
    score_nodes(adata, panels)
    calls: dict[str, dict] = {}
    sfx: dict[str, dict] = {}
    for key in cluster_keys:
        cs = assign_all(
            adata.obs, dag, panels, key, overrides,
            descend_agree=descend_agree, gate_overlap=gate_overlap, gate_score_r=gate_score_r,
        )
        calls[key] = {cl: dataclasses.asdict(c) for cl, c in cs.items()}
        if suffixes:
            clusters = adata.obs[key].astype(str)
            sfx[key] = suffixes_for(adata, clusters, calls[key])
        label_of = {cl: " > ".join(d["path"]) or "Unknown" for cl, d in calls[key].items()}
        adata.obs[f"treeline_{key}"] = adata.obs[key].astype(str).map(label_of)
    adata.uns[UNS_KEY] = json.dumps(
        {
            "calls": calls,
            "suffixes": sfx,
            "panels": panels,
            "cluster_keys": list(cluster_keys),
            "dag_fetched": dag.get("fetched"),
            "overrides": overrides or [],
        }
    )
    return adata


def annotations(adata) -> dict:
    """Decode `.uns['treeline']` from an annotated AnnData."""
    return json.loads(adata.uns[UNS_KEY])
