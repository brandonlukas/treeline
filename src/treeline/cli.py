"""The treeline CLI: one prep step, three capabilities (SPECS scope).

    treeline derive-tree gene_sets.csv cl_dag.json
    treeline annotate in.h5ad --dag cl_dag.json --gene-sets gene_sets.csv \
        --clusters leiden_0.5 leiden_1.0 [--overrides overrides.json] -o annotated.h5ad
    treeline substates annotated.h5ad --clusters leiden_2.0 -o annotated.h5ad
    treeline integrate a_annotated.h5ad b_annotated.h5ad -o joint.h5ad
    treeline colors annotated.h5ad [more.h5ad ...] -o palette.json
    treeline report annotated.h5ad [more.h5ad ...] -o report.html
    treeline summary annotated.h5ad

`annotate` expects log-normalized X and cluster labels already in .obs (treeline
performs no clustering). `integrate` expects annotated h5ads with a counts layer and
emits the joint object with the latent in obsm["X_treeline"] — recluster it yourself
and resubmit through `annotate`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(prog="treeline")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("derive-tree", help="gene-set CSV -> frozen, dated CL DAG (via EBI OLS4)")
    d.add_argument("gene_sets", help="CL-labeled gene-set CSV (CellGuide export: 'Gene Set Name', 'Gene Symbol')")
    d.add_argument("out", help="frozen DAG JSON, the --dag input to annotate")

    a = sub.add_parser("annotate", help="multi-level annotation, written into the AnnData")
    a.add_argument("h5ad", help="clustered AnnData; X log-normalized, cluster labels in .obs")
    a.add_argument("--dag", required=True, help="DAG JSON from derive-tree")
    a.add_argument("--gene-sets", required=True, help="the same gene-set CSV given to derive-tree")
    a.add_argument("--clusters", nargs="+", required=True, metavar="OBS_KEY",
                   help=".obs columns to annotate (cheap: pass every resolution you have)")
    a.add_argument("--overrides", help="expert overrides JSON: signed {node, decision, justification, author, date} entries")
    a.add_argument("--n-markers", type=int, default=None, help="top genes per gene set (default 10)")
    a.add_argument("--descend-agree", type=float, default=None, help="vote share to descend (default 0.5)")
    a.add_argument("--gate-overlap", type=float, default=None, help="gate: exclusive-gene overlap (default 0.5)")
    a.add_argument("--gate-score-r", type=float, default=None, help="gate: exclusive-score correlation (default 0.9)")
    a.add_argument("-o", "--out", required=True, help="annotated h5ad (self-contained; input to every other verb)")

    n = sub.add_parser("substates", help="NS-Forest {gene}+ substate naming — costly; run on the clustering you settle on")
    n.add_argument("h5ad", help="annotated h5ad from annotate")
    n.add_argument("--clusters", nargs="+", required=True, metavar="OBS_KEY", help="annotated clusterings to name substates in")
    n.add_argument("--min-cluster", type=int, default=None, help="skip clusters smaller than this (default 20)")
    n.add_argument("--n-trees", type=int, default=None, help="random forest size (default 100)")
    n.add_argument("--top-gini", type=int, default=None, help="genes kept by Gini importance (default 15)")
    n.add_argument("--top-combo", type=int, default=None, help="genes entering the combination search (default 6)")
    n.add_argument("--beta", type=float, default=None, help="F-beta weight, <1 favours precision (default 0.5)")
    n.add_argument("-o", "--out", required=True, help="output h5ad (may equal the input)")

    i = sub.add_parser("integrate", help="annotated h5ads -> joint h5ad with the scANVI latent")
    i.add_argument("h5ads", nargs="+", help="two or more annotated h5ads with layers['counts']; sample name = stem minus _annotated")
    i.add_argument("--supervise-depth", type=int, default=None, help="tree cut for the prior (default 2)")
    i.add_argument("--n-hvg", type=int, default=None, help="highly variable genes for training (default 2000)")
    i.add_argument("--n-latent", type=int, default=None, help="latent dimensions (default 30)")
    i.add_argument("--scvi-epochs", type=int, default=None, help="scVI pretraining epochs (default 100)")
    i.add_argument("--scanvi-epochs", type=int, default=None, help="scANVI epochs (default 50)")
    i.add_argument("--classification-ratio", type=float, default=None,
                   help="label pull vs mixing (default 50; a judgment call when sample==condition)")
    i.add_argument("-o", "--out", required=True, help="joint h5ad with the latent in obsm['X_treeline']")

    c = sub.add_parser("colors", help="annotated h5ads -> hierarchical label palette (JSON)")
    c.add_argument("h5ads", nargs="+", help="annotated h5ads")
    c.add_argument("-o", "--out", required=True, help="JSON {label: hex}; siblings share a hue, substates are tints")

    r = sub.add_parser("report", help="annotated h5ads -> static HTML report (slider, tree, tables)")
    r.add_argument("h5ads", nargs="+", help="annotated h5ads with a 2-D embedding; sample name = stem minus _annotated")
    r.add_argument("--title", help="report title (default: derived from filenames)")
    r.add_argument("--basis", default="X_umap", help="2-D obsm embedding to plot (default X_umap)")
    r.add_argument("-o", "--out", required=True, help="self-contained HTML file")

    m = sub.add_parser("summary", help="print an annotated h5ad's labels, per clustering")
    m.add_argument("h5ad", help="annotated h5ad; prints each cluster's path, shares, substate markers, gate refusals and overrides")

    args = p.parse_args()

    def kwargs(*names):  # only knobs the user set; the functions carry the defaults
        return {n: getattr(args, n) for n in names if getattr(args, n) is not None}

    if args.cmd == "derive-tree":
        from treeline.tree import derive_tree

        dag = derive_tree(args.gene_sets)
        Path(args.out).write_text(json.dumps(dag, indent=2))
        print(f"wrote {args.out}: {len(dag['nodes'])} nodes, roots={dag['roots']}")

    elif args.cmd == "annotate":
        import scanpy as sc

        from treeline.annotate import annotate
        from treeline.tree import load
        from treeline.vote import load_overrides

        overrides = load_overrides(args.overrides) if args.overrides else []  # before the slow read
        adata = sc.read_h5ad(args.h5ad)
        annotate(adata, load(args.dag), args.gene_sets, args.clusters, overrides,
                 **kwargs("n_markers", "descend_agree", "gate_overlap", "gate_score_r"))
        adata.write_h5ad(args.out)
        print(f"wrote {args.out}")

    elif args.cmd == "substates":
        import scanpy as sc

        from treeline.annotate import add_substates

        adata = sc.read_h5ad(args.h5ad)
        add_substates(adata, args.clusters, **kwargs("min_cluster", "n_trees", "top_gini", "top_combo", "beta"))
        adata.write_h5ad(args.out)
        print(f"wrote {args.out}")

    elif args.cmd == "integrate":
        import scanpy as sc

        from treeline.annotate import sample_names
        from treeline.scanvi import integrate

        joint = integrate({n: sc.read_h5ad(f) for n, f in sample_names(args.h5ads).items()},
                          **kwargs("supervise_depth", "n_hvg", "n_latent", "scvi_epochs", "scanvi_epochs",
                                   "classification_ratio"))
        joint.write_h5ad(args.out)
        print(f"wrote {args.out} — recluster obsm['X_treeline'] and resubmit through annotate")

    elif args.cmd == "colors":
        import scanpy as sc

        from treeline.annotate import annotations
        from treeline.colors import palette

        anns = [annotations(sc.read_h5ad(f)) for f in args.h5ads]
        Path(args.out).write_text(json.dumps(palette(*anns), indent=1))
        print(f"wrote {args.out}")

    elif args.cmd == "report":
        import scanpy as sc

        from treeline.annotate import sample_names
        from treeline.report import render_report

        render_report({n: sc.read_h5ad(f) for n, f in sample_names(args.h5ads).items()}, args.out, args.title, args.basis)

    elif args.cmd == "summary":
        import scanpy as sc

        from treeline.annotate import summary

        summary(sc.read_h5ad(args.h5ad, backed="r"))  # obs + uns only; X stays on disk


if __name__ == "__main__":
    main()
