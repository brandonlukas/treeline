"""The integrate verb: scANVI semi-supervised by tree-cut treeline labels.

Input: two or more *annotated* AnnDatas (see annotate.py; each needs a
`layers["counts"]` matrix — training runs on counts). Supervision labels are a PRIOR,
not truth: per-nucleus `treeline_*` labels truncated at SUPERVISE_DEPTH, kept only where
every supplied clustering that resolves agrees (cross-resolution consensus); flips
become "Unknown", scANVI's let-the-data-decide class. NS-Forest suffixes are never
supervision.

`integrate` builds the joint object, trains, attaches the latent as
`obsm["X_treeline"]`, and stops — treeline performs no clustering anywhere (SPECS
scope). The user reclusters the latent and resubmits the joint through `annotate`;
apps/joint_1619.py is the POC driver for that loop. `refine_classes` is the opt-in
within-class convenience it exposes as --refine.
"""

from __future__ import annotations

import json

import anndata as ad
import pandas as pd
import scanpy as sc

from treeline.annotate import annotate

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


def consensus_labels(obs: pd.DataFrame) -> pd.Series:
    """Multi-resolution consensus over the per-nucleus `treeline_*` label columns: the
    depth-cut label where every clustering that resolves agrees; flips -> Unknown
    (labels are a prior, not truth — cross-resolution agreement is its strength)."""
    per = []
    for col in [c for c in obs.columns if c.startswith("treeline_")]:
        cut = obs[col].astype(str).map(
            lambda p: " > ".join(p.split(" > ")[:SUPERVISE_DEPTH])
            if p != "Unknown" and len(p.split(" > ")) >= SUPERVISE_DEPTH
            else None
        )
        per.append(cut)
    df = pd.concat(per, axis=1)
    label = df.bfill(axis=1).iloc[:, 0].fillna("Unknown")
    return label.where(df.nunique(axis=1) == 1, "Unknown")


def integrate(adatas: dict[str, "ad.AnnData"]) -> "ad.AnnData":
    """Annotated AnnDatas -> one joint AnnData with the scANVI latent in
    obsm["X_treeline"] and the consensus prior in obs["cl_prior"]. No clustering."""
    import scvi

    joint = ad.concat(adatas, label="sample", index_unique="-")
    joint.obs["cl_prior"] = consensus_labels(joint.obs).values
    print("consensus prior:", joint.obs["cl_prior"].value_counts().to_dict())

    sub = joint.copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts", batch_key="sample")
    sub = sub[:, sub.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(sub, layer="counts", batch_key="sample", labels_key="cl_prior")
    m = scvi.model.SCVI(sub, n_latent=N_LATENT)
    m.train(max_epochs=SCVI_EPOCHS)
    ms = scvi.model.SCANVI.from_scvi_model(m, unlabeled_category="Unknown")
    ms.train(
        max_epochs=SCANVI_EPOCHS,
        n_samples_per_label=100,
        plan_kwargs={"classification_ratio": CLASSIFICATION_RATIO},
    )
    joint.obsm["X_treeline"] = ms.get_latent_representation()
    return joint


def refine_classes(joint, dag: dict, gene_sets_csv, out_dir) -> None:
    """Opt-in within-class refinement: per-class HVGs + a fresh per-class integration
    (batch correction acts only inside a type), subcluster, re-vote, NS-Forest — one
    annotated h5ad per class. `joint` must already be annotated (post joint-annotate)
    with a counts layer and log-normalized X."""
    import scvi

    class_label = consensus_labels(joint.obs)
    for label, n in class_label.value_counts().items():
        if label == "Unknown" or n < REFINE_MIN:
            continue
        leaf = label.split(" > ")[-1]
        print(f"refined · {leaf}: {n} nuclei")
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
        key = f"leiden_{REFINE_RES}"
        sc.tl.leiden(sub, resolution=REFINE_RES, key_added=key, flavor="igraph", n_iterations=2)
        for c in [c for c in sub.obs.columns if c.startswith("treeline_")]:
            del sub.obs[c]  # stale labels from the parent object
        annotate(sub, dag, gene_sets_csv, [key])
        for cl, d in json.loads(sub.uns["treeline"])["calls"][key].items():
            print(f"  {cl}: {' > '.join(d['path']) or 'Unknown'}")
        sub.write_h5ad(out_dir / f"refined__{leaf}_annotated.h5ad")
