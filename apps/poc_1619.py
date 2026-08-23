"""POC driver: treeline over assets/ (patient 1619 LM + MM snMultiome GEX).

Clustering is upstream of treeline's contract — this driver supplies it (Leiden at two
resolutions) because the assets h5ads ship unclustered. Then: score the DAG's panels,
vote per (sample, resolution, cluster), embed a joint UMAP (plain concat, no batch
correction — labeled honestly), and write intermediates for the report step.

    .venv/bin/python apps/poc_1619.py
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pandas as pd
import scanpy as sc

from treeline.tree import load
from treeline.vote import assign_all, load_overrides, load_panels, score_nodes

ASSETS = Path("assets")
OUT = Path("results/poc")
SAMPLES = ["1619LM", "1619MM"]
RESOLUTIONS = [0.5, 1.0, 2.0]


def process(adata) -> None:
    """Standard scanpy: normalize, HVG, PCA, neighbors, UMAP + Leiden per resolution."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata, n_comps=50, mask_var="highly_variable")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    for res in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dag = load(ASSETS / "cl_dag.json")
    panels = load_panels(ASSETS / "cellguide_uterus_gene_sets_2026-08-22.csv", dag)
    overrides = load_overrides(ASSETS / "overrides.json")

    calls: dict[str, dict[str, dict[str, dict]]] = {}
    per_sample = {}
    for sample in SAMPLES:
        adata = sc.read_h5ad(ASSETS / f"{sample}_gex.h5ad")
        print(f"{sample}: {adata.shape}")
        process(adata)
        score_nodes(adata, panels)
        calls[sample] = {}
        for res in RESOLUTIONS:
            key = f"leiden_{res}"
            sample_calls = assign_all(adata.obs, dag, panels, key, overrides)
            calls[sample][str(res)] = {cl: dataclasses.asdict(c) for cl, c in sample_calls.items()}
            print(f"  {key}: {adata.obs[key].nunique()} clusters")
            for cl, c in sample_calls.items():
                note = f"  [{c.refused}]" if c.refused else (f"  [{c.overridden}]" if c.overridden else "")
                print(f"    {cl}: {' > '.join(c.path) or 'Unknown'}{note}")
        score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
        obs = adata.obs[[f"leiden_{r}" for r in RESOLUTIONS] + score_cols].copy()
        obs[["umap1", "umap2"]] = adata.obsm["X_umap"]
        obs.to_parquet(OUT / f"{sample}.parquet")
        per_sample[sample] = adata

    # joint embedding: plain concat, no batch correction — the report labels it as such.
    # Labels come per nucleus from each sample's own clustering; CL terms are the shared
    # vocabulary, so no joint clustering is needed.
    import anndata as ad

    joint = ad.concat(
        {s: sc.read_h5ad(ASSETS / f"{s}_gex.h5ad") for s in SAMPLES}, label="sample", index_unique="-"
    )
    sc.pp.normalize_total(joint, target_sum=1e4)
    sc.pp.log1p(joint)
    sc.pp.highly_variable_genes(joint, n_top_genes=2000)
    sc.pp.pca(joint, n_comps=50, mask_var="highly_variable")
    sc.pp.neighbors(joint)
    sc.tl.umap(joint)
    jobs = joint.obs[["sample"]].copy()
    jobs[["umap1", "umap2"]] = joint.obsm["X_umap"]
    for res in RESOLUTIONS:
        jobs[f"leiden_{res}"] = -1
    for sample in SAMPLES:
        src = per_sample[sample].obs
        mask = jobs["sample"] == sample
        stripped = jobs.index[mask].str.removesuffix(f"-{sample}")
        for res in RESOLUTIONS:
            jobs.loc[mask, f"leiden_{res}"] = src.loc[stripped, f"leiden_{res}"].astype(int).values
    jobs.to_parquet(OUT / "joint.parquet")

    (OUT / "calls.json").write_text(json.dumps(calls, indent=1))
    (OUT / "panels.json").write_text(json.dumps(panels, indent=1))
    print("done")


if __name__ == "__main__":
    main()
