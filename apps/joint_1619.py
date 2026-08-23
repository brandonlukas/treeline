"""POC joint driver: cluster treeline's integrated latent, resubmit for annotation.

The user-side half of the integration loop (SPECS scope: treeline emits the latent and
performs no clustering). This driver plays the user: Leiden on the latent at three
resolutions, then the joint expression + clusterings go back through the same annotate
stage — vote, then NS-Forest per resolution. --refine additionally runs the opt-in
within-class refinement convenience.

    .venv/bin/python apps/joint_1619.py            (after python -m treeline.scanvi)
    .venv/bin/python apps/joint_1619.py --refine
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc

from treeline.nsforest import suffixes_for
from treeline.scanvi import load_joint, refine_classes
from treeline.tree import load
from treeline.vote import assign_all, load_panels, score_nodes

ASSETS = Path("assets")
OUT = Path("results/poc")
RESOLUTIONS = [0.5, 1.0, 2.0]


def main() -> None:
    calls = json.loads((OUT / "calls.json").read_text())
    joint = load_joint(OUT, ASSETS, calls)
    latent = pd.read_parquet(OUT / "integrated_latent.parquet")
    joint.obsm["X_scANVI"] = latent[[c for c in latent.columns if c.startswith("z")]].loc[joint.obs_names].values

    # driver business: cluster the latent (the "user reclusters" step of the loop)
    sc.pp.neighbors(joint, use_rep="X_scANVI")
    sc.tl.umap(joint)
    for res in RESOLUTIONS:
        sc.tl.leiden(joint, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2)

    # resubmit to the annotate stage: same vote, same rules, on the joint clusters
    sc.pp.normalize_total(joint, target_sum=1e4)
    sc.pp.log1p(joint)
    dag = load(ASSETS / "cl_dag.json")
    panels = load_panels(ASSETS / "cellguide_uterus_gene_sets_2026-08-22.csv", dag)
    score_nodes(joint, panels)
    calls["integrated"] = {}
    suffixes = json.loads((OUT / "suffixes.json").read_text())
    suffixes["integrated"] = {}
    for res in RESOLUTIONS:
        cs = assign_all(joint.obs, dag, panels, f"leiden_{res}")
        calls["integrated"][str(res)] = {cl: dataclasses.asdict(c) for cl, c in cs.items()}
        print(f"integrated leiden_{res}: {len(cs)} clusters")
        for cl, c in cs.items():
            print(f"  {cl}: {' > '.join(c.path) or 'Unknown'}")
        clusters = joint.obs[f"leiden_{res}"].astype(str)
        clusters.index = joint.obs_names
        suffixes["integrated"][str(res)] = suffixes_for(joint, clusters, calls["integrated"][str(res)])
        print(f"  {len(suffixes['integrated'][str(res)])} suffixed clusters")

    score_cols = [c for c in joint.obs.columns if c.startswith("score_")]
    out = joint.obs[["sample"] + [f"leiden_{r}" for r in RESOLUTIONS] + score_cols].copy()
    out[["umap1", "umap2"]] = joint.obsm["X_umap"]
    out.to_parquet(OUT / "integrated.parquet")

    if "--refine" in sys.argv:
        refine_classes(joint, calls, suffixes, dag, panels, OUT)

    (OUT / "calls.json").write_text(json.dumps(calls, indent=1))
    (OUT / "suffixes.json").write_text(json.dumps(suffixes, indent=1))
    print("done")


if __name__ == "__main__":
    main()
