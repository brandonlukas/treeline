"""The treeline CLI: one prep step, three capabilities (SPECS scope).

    treeline derive-tree gene_sets.csv cl_dag.json
    treeline annotate in.h5ad --dag cl_dag.json --gene-sets gene_sets.csv \
        --clusters leiden_0.5 leiden_1.0 [--overrides overrides.json] -o annotated.h5ad
    treeline integrate a_annotated.h5ad b_annotated.h5ad -o joint.h5ad
    treeline colors annotated.h5ad [more.h5ad ...] -o palette.json

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
    d.add_argument("gene_sets")
    d.add_argument("out")

    a = sub.add_parser("annotate", help="multi-level annotation, written into the AnnData")
    a.add_argument("h5ad")
    a.add_argument("--dag", required=True)
    a.add_argument("--gene-sets", required=True)
    a.add_argument("--clusters", nargs="+", required=True, metavar="OBS_KEY")
    a.add_argument("--overrides")
    a.add_argument("--no-suffixes", action="store_true")
    a.add_argument("-o", "--out", required=True)

    i = sub.add_parser("integrate", help="annotated h5ads -> joint h5ad with the scANVI latent")
    i.add_argument("h5ads", nargs="+")
    i.add_argument("-o", "--out", required=True)

    c = sub.add_parser("colors", help="annotated h5ads -> hierarchical label palette (JSON)")
    c.add_argument("h5ads", nargs="+")
    c.add_argument("-o", "--out", required=True)

    args = p.parse_args()

    if args.cmd == "derive-tree":
        from treeline.tree import derive, set_names_from_csv

        dag = derive(set_names_from_csv(args.gene_sets))
        Path(args.out).write_text(json.dumps(dag, indent=2))
        print(f"wrote {args.out}: {len(dag['nodes'])} nodes, roots={dag['roots']}")

    elif args.cmd == "annotate":
        import scanpy as sc

        from treeline.annotate import annotate
        from treeline.tree import load
        from treeline.vote import load_overrides

        adata = sc.read_h5ad(args.h5ad)
        overrides = load_overrides(args.overrides) if args.overrides else []
        annotate(adata, load(args.dag), args.gene_sets, args.clusters, overrides, suffixes=not args.no_suffixes)
        adata.write_h5ad(args.out)
        print(f"wrote {args.out}")

    elif args.cmd == "integrate":
        import scanpy as sc

        from treeline.scanvi import integrate

        joint = integrate({Path(f).stem.removesuffix("_annotated"): sc.read_h5ad(f) for f in args.h5ads})
        joint.write_h5ad(args.out)
        print(f"wrote {args.out} — recluster obsm['X_treeline'] and resubmit through annotate")

    elif args.cmd == "colors":
        import scanpy as sc

        from treeline.annotate import annotations
        from treeline.colors import palette

        anns = [annotations(sc.read_h5ad(f)) for f in args.h5ads]
        Path(args.out).write_text(json.dumps(palette(*anns), indent=1))
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
