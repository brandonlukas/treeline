"""POC per-sample driver: treeline `annotate` over assets/ (patient 1619 LM + MM).

Clustering is upstream of treeline's contract — this driver supplies it (Leiden at
three resolutions) because the assets h5ads ship unclustered, then hands each sample to
the annotate verb and writes the self-contained annotated h5ads the rest of the loop
consumes (apps/joint_1619.py, apps/report.py).

    .venv/bin/python apps/poc_1619.py
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

from treeline.annotate import annotate
from treeline.tree import load
from treeline.vote import load_overrides

ASSETS = Path("assets")
OUT = Path("results/poc")
SAMPLES = ["1619LM", "1619MM"]
RESOLUTIONS = [0.5, 1.0, 2.0]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dag = load(ASSETS / "cl_dag.json")
    overrides = load_overrides(ASSETS / "overrides.json")
    for sample in SAMPLES:
        adata = sc.read_h5ad(ASSETS / f"{sample}_gex.h5ad")
        print(f"{sample}: {adata.shape}")
        adata.layers["counts"] = adata.X.copy()  # integrate trains on counts
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.pp.pca(adata, n_comps=50, mask_var="highly_variable")
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        for res in RESOLUTIONS:
            sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2)
        # no substates per sample: NS-Forest runs once, on the integrated object
        annotate(adata, dag, ASSETS / "cellguide_uterus_gene_sets_2026-08-22.csv",
                 [f"leiden_{r}" for r in RESOLUTIONS], overrides)
        for key in [f"leiden_{r}" for r in RESOLUTIONS]:
            print(f"  {key}: {adata.obs[key].nunique()} clusters -> "
                  f"{adata.obs[f'treeline_{key}'].nunique()} labels")
        adata.write_h5ad(OUT / f"{sample}_annotated.h5ad")
    print("done")


if __name__ == "__main__":
    main()
