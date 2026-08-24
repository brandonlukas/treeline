"""POC joint driver: integrate the annotated samples, recluster the latent, resubmit.

The user-side half of the integration loop (SPECS scope: treeline emits the latent and
performs no clustering). This driver plays the user: Leiden on obsm["X_treeline"] at
three resolutions, then the joint object goes back through the same annotate verb.

    .venv/bin/python apps/joint_1619.py      (after apps/poc_1619.py)
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

from treeline.annotate import add_substates, annotate
from treeline.scanvi import integrate
from treeline.tree import load
from treeline.vote import load_overrides

ASSETS = Path("assets")
OUT = Path("results/poc")
SAMPLES = ["1619LM", "1619MM"]
RESOLUTIONS = [0.5, 1.0, 2.0]
GENE_SETS = ASSETS / "cellguide_uterus_gene_sets_2026-08-22.csv"


def main() -> None:
    joint = integrate({s: sc.read_h5ad(OUT / f"{s}_annotated.h5ad") for s in SAMPLES})

    # driver business: cluster the latent (the "user reclusters" step of the loop)
    sc.pp.neighbors(joint, use_rep="X_treeline")
    sc.tl.umap(joint)
    for res in RESOLUTIONS:
        sc.tl.leiden(joint, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2)

    # resubmit through the annotate verb: same vote, same rules, on the joint clusters
    dag = load(ASSETS / "cl_dag.json")
    annotate(joint, dag, GENE_SETS, [f"leiden_{r}" for r in RESOLUTIONS],
             load_overrides(ASSETS / "overrides.json"))
    add_substates(joint, ["leiden_2.0"])  # costly: final resolution only
    for key in [f"leiden_{r}" for r in RESOLUTIONS]:
        print(f"  {key}: {joint.obs[key].nunique()} clusters -> "
              f"{joint.obs[f'treeline_{key}'].nunique()} labels")
    joint.write_h5ad(OUT / "integrated_annotated.h5ad")
    print("done")


if __name__ == "__main__":
    main()
