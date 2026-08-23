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
apps/joint_1619.py is the POC driver for that loop.
"""

from __future__ import annotations

import anndata as ad
import pandas as pd
import scanpy as sc

SUPERVISE_DEPTH = 2
N_HVG = 2000
N_LATENT = 30
SCVI_EPOCHS = 100
SCANVI_EPOCHS = 50
# label pull vs mixing. Sample==condition here (one LM, one MM), so batch and disease
# effects are unidentifiable — this value is a judgment call, not fittable (see SPECS).
CLASSIFICATION_RATIO = 50.0


def consensus_labels(obs: pd.DataFrame, supervise_depth: int = SUPERVISE_DEPTH) -> pd.Series:
    """Multi-resolution consensus over the per-nucleus `treeline_*` label columns: the
    depth-cut label where every clustering that resolves agrees; flips -> Unknown
    (labels are a prior, not truth — cross-resolution agreement is its strength)."""
    per = []
    for col in [c for c in obs.columns if c.startswith("treeline_")]:
        cut = obs[col].astype(str).map(
            lambda p: " > ".join(p.split(" > ")[:supervise_depth])
            if p != "Unknown" and len(p.split(" > ")) >= supervise_depth
            else None
        )
        per.append(cut)
    df = pd.concat(per, axis=1)
    label = df.bfill(axis=1).iloc[:, 0].fillna("Unknown")
    return label.where(df.nunique(axis=1) == 1, "Unknown")


def integrate(
    adatas: dict[str, ad.AnnData],
    *,
    supervise_depth: int = SUPERVISE_DEPTH,
    n_hvg: int = N_HVG,
    n_latent: int = N_LATENT,
    scvi_epochs: int = SCVI_EPOCHS,
    scanvi_epochs: int = SCANVI_EPOCHS,
    classification_ratio: float = CLASSIFICATION_RATIO,
) -> ad.AnnData:
    """Annotated AnnDatas -> one joint AnnData with the scANVI latent in
    obsm["X_treeline"] and the consensus prior in obs["cl_prior"]. No clustering."""
    import scvi

    joint = ad.concat(adatas, label="sample", index_unique="-")
    joint.obs["cl_prior"] = consensus_labels(joint.obs, supervise_depth).values
    print("consensus prior:", joint.obs["cl_prior"].value_counts().to_dict())

    sub = joint.copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts", batch_key="sample")
    sub = sub[:, sub.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(sub, layer="counts", batch_key="sample", labels_key="cl_prior")
    m = scvi.model.SCVI(sub, n_latent=n_latent)
    m.train(max_epochs=scvi_epochs)
    ms = scvi.model.SCANVI.from_scvi_model(m, unlabeled_category="Unknown")
    ms.train(
        max_epochs=scanvi_epochs,
        n_samples_per_label=100,
        plan_kwargs={"classification_ratio": classification_ratio},
    )
    joint.obsm["X_treeline"] = ms.get_latent_representation()
    return joint
