"""Label-aware integration: scANVI semi-supervised by tree-cut treeline labels.

Supervision labels are a PRIOR, not truth: descent paths truncated at SUPERVISE_DEPTH,
kept only where every supplied clustering resolution that resolves agrees
(cross-resolution consensus); nuclei that flip between classes are "Unknown" — scANVI's
let-the-data-decide class. NS-Forest suffixes are never supervision.

This stage integrates and EMITS THE LATENT, nothing more (SPECS scope: treeline performs
no clustering anywhere). The user reclusters the latent at their discretion and resubmits
the joint expression + clusterings to the annotate stage — the POC driver for that loop
is apps/joint_1619.py. `refine_classes` (within-class integration + subcluster + re-vote
+ NS-Forest) is the opt-in convenience that driver exposes as --refine.

    python -m treeline.scanvi results/poc assets    # writes results/poc/integrated_latent.parquet
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from treeline.nsforest import suffixes_for
from treeline.vote import assign_all

SUPERVISE_DEPTH = 2
N_HVG = 2000
N_LATENT = 30
SCVI_EPOCHS = 100
SCANVI_EPOCHS = 50
# label pull vs mixing. Sample==condition here (one LM, one MM), so batch and disease
# effects are unidentifiable — this value is a judgment call, not fittable (see SPECS).
CLASSIFICATION_RATIO = 50.0
REFINE_MIN = 800  # nuclei a coarse class needs before it gets its own integration round
REFINE_RES = 1.0


def cut_label(path: list[str]) -> str:
    return " > ".join(path[:SUPERVISE_DEPTH]) if len(path) >= SUPERVISE_DEPTH else "Unknown"


def consensus_labels(calls_by_res: dict[str, dict], cluster_ids: dict[str, pd.Series]) -> pd.Series:
    """Multi-resolution consensus: the depth-cut label where every resolution that
    resolves agrees; nuclei that flip between classes -> Unknown (labels are a prior,
    not truth — cross-resolution agreement is its strength)."""
    per_res = []
    for res, cl in cluster_ids.items():
        label_of = {c: cut_label(v["path"]) for c, v in calls_by_res[res].items()}
        lab = cl.astype(str).map(label_of)
        per_res.append(lab.where(lab != "Unknown"))
    df = pd.concat(per_res, axis=1)
    n_distinct = df.nunique(axis=1)
    label = df.bfill(axis=1).iloc[:, 0].fillna("Unknown")
    return label.where(n_distinct == 1, "Unknown")


def load_joint(results: Path, assets: Path, calls: dict) -> "ad.AnnData":
    """Concat the per-sample h5ads with consensus prior labels; counts kept as a layer."""
    samples = [s for s in calls if not (s == "integrated" or s.startswith("refined"))]
    parts = {}
    for s in samples:
        a = sc.read_h5ad(assets / f"{s}_gex.h5ad")
        df = pd.read_parquet(results / f"{s}.parquet")
        res_keys = [c.removeprefix("leiden_") for c in df.columns if c.startswith("leiden_")]
        cluster_ids = {r: df[f"leiden_{r}"] for r in res_keys}
        a.obs["cl_label"] = consensus_labels(calls[s], cluster_ids).values
        parts[s] = a
    joint = ad.concat(parts, label="sample", index_unique="-")
    joint.layers["counts"] = joint.X.copy()
    return joint


def integrate(joint: "ad.AnnData") -> pd.DataFrame:
    """SCVI -> SCANVI on the consensus prior; returns the latent as a DataFrame."""
    import scvi

    sub = joint.copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts", batch_key="sample")
    sub = sub[:, sub.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(sub, layer="counts", batch_key="sample", labels_key="cl_label")
    m = scvi.model.SCVI(sub, n_latent=N_LATENT)
    m.train(max_epochs=SCVI_EPOCHS)
    ms = scvi.model.SCANVI.from_scvi_model(m, unlabeled_category="Unknown")
    ms.train(
        max_epochs=SCANVI_EPOCHS,
        n_samples_per_label=100,
        plan_kwargs={"classification_ratio": CLASSIFICATION_RATIO},
    )
    z = ms.get_latent_representation()
    return pd.DataFrame(z, index=joint.obs_names, columns=[f"z{i}" for i in range(z.shape[1])])


def refine_classes(joint, calls: dict, suffixes: dict, dag: dict, panels: dict, results: Path) -> None:
    """Opt-in within-class refinement: per-class HVGs + a fresh per-class integration
    (batch correction acts only inside a type), subcluster, re-vote, NS-Forest.
    `joint` needs: counts layer, log-normalized X, score_ columns, leiden_ columns."""
    import scvi

    res_keys = [c.removeprefix("leiden_") for c in joint.obs.columns if c.startswith("leiden_")]
    joint_ids = {r: joint.obs[f"leiden_{r}"] for r in res_keys}
    class_label = consensus_labels(calls["integrated"], joint_ids)
    for label, n in class_label.value_counts().items():
        if label == "Unknown" or n < REFINE_MIN:
            continue
        leaf = label.split(" > ")[-1]
        key = f"refined · {leaf}"
        print(f"{key}: {n} nuclei")
        sub = joint[(class_label == label).values].copy()
        # seurat flavor on the log-normalized subset: seurat_v3's loess is numerically
        # fragile on small within-class subsets; scvi still trains on the counts layer
        sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG, batch_key="sample")
        subh = sub[:, sub.var["highly_variable"]].copy()
        scvi.model.SCVI.setup_anndata(subh, layer="counts", batch_key="sample")
        mr = scvi.model.SCVI(subh, n_latent=15)
        mr.train(max_epochs=SCVI_EPOCHS)
        sub.obsm["X_refined"] = mr.get_latent_representation()
        sc.pp.neighbors(sub, use_rep="X_refined")
        sc.tl.umap(sub)
        sc.tl.leiden(sub, resolution=REFINE_RES, key_added=f"leiden_{REFINE_RES}", flavor="igraph", n_iterations=2)
        cs = assign_all(sub.obs, dag, panels, f"leiden_{REFINE_RES}")
        calls[key] = {str(REFINE_RES): {cl: dataclasses.asdict(c) for cl, c in cs.items()}}
        for cl, c in cs.items():
            print(f"  {cl}: {' > '.join(c.path) or 'Unknown'}")
        clusters = sub.obs[f"leiden_{REFINE_RES}"].astype(str)
        clusters.index = sub.obs_names
        suffixes[key] = {str(REFINE_RES): suffixes_for(sub, clusters, calls[key][str(REFINE_RES)])}
        rout = sub.obs[["sample", f"leiden_{REFINE_RES}"]].copy()
        rout[["umap1", "umap2"]] = sub.obsm["X_umap"]
        rout.to_parquet(results / f"refined__{leaf}.parquet")


def main() -> None:
    results, assets = Path(sys.argv[1]), Path(sys.argv[2])
    calls = json.loads((results / "calls.json").read_text())
    joint = load_joint(results, assets, calls)
    print("consensus labels:", joint.obs["cl_label"].value_counts().to_dict())
    latent = integrate(joint)
    latent.insert(0, "sample", joint.obs["sample"].values)
    latent.insert(1, "cl_label", joint.obs["cl_label"].values)
    latent.to_parquet(results / "integrated_latent.parquet")
    print(f"wrote {results / 'integrated_latent.parquet'} — recluster it and resubmit (apps/joint_1619.py)")


if __name__ == "__main__":
    main()
