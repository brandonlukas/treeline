"""The treeline CLI: one prep step, three capabilities (SPECS scope).

    treeline derive-tree gene_sets.csv cl_dag.json
    treeline annotate in.h5ad --dag cl_dag.json --gene-sets gene_sets.csv \
        --clusters leiden_0.5 leiden_1.0 [--overrides overrides.json] -o annotated.h5ad
    treeline substates annotated.h5ad --clusters leiden_2.0 -o annotated.h5ad
    treeline integrate a_annotated.h5ad b_annotated.h5ad -o joint.h5ad
    treeline colors annotated.h5ad [more.h5ad ...] -o palette.json
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
    d.add_argument("gene_sets")
    d.add_argument("out")

    a = sub.add_parser("annotate", help="multi-level annotation, written into the AnnData")
    a.add_argument("h5ad")
    a.add_argument("--dag", required=True)
    a.add_argument("--gene-sets", required=True)
    a.add_argument("--clusters", nargs="+", required=True, metavar="OBS_KEY")
    a.add_argument("--overrides")
    a.add_argument("--n-markers", type=int, default=None, help="top genes per gene set (default 10)")
    a.add_argument("--descend-agree", type=float, default=None, help="vote share to descend (default 0.5)")
    a.add_argument("--gate-overlap", type=float, default=None, help="gate: exclusive-gene overlap (default 0.5)")
    a.add_argument("--gate-score-r", type=float, default=None, help="gate: exclusive-score correlation (default 0.9)")
    a.add_argument("-o", "--out", required=True)

    n = sub.add_parser("substates", help="NS-Forest {gene}+ substate naming — costly; run on the clustering you settle on")
    n.add_argument("h5ad")
    n.add_argument("--clusters", nargs="+", required=True, metavar="OBS_KEY")
    n.add_argument("-o", "--out", required=True)

    i = sub.add_parser("integrate", help="annotated h5ads -> joint h5ad with the scANVI latent")
    i.add_argument("h5ads", nargs="+")
    i.add_argument("--supervise-depth", type=int, default=None, help="tree cut for the prior (default 2)")
    i.add_argument("--n-latent", type=int, default=None, help="latent dimensions (default 30)")
    i.add_argument("--classification-ratio", type=float, default=None,
                   help="label pull vs mixing (default 50; a judgment call when sample==condition)")
    i.add_argument("-o", "--out", required=True)

    c = sub.add_parser("colors", help="annotated h5ads -> hierarchical label palette (JSON)")
    c.add_argument("h5ads", nargs="+")
    c.add_argument("-o", "--out", required=True)

    m = sub.add_parser("summary", help="print an annotated h5ad's labels, per clustering")
    m.add_argument("h5ad")

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
        kw = {k: v for k, v in [("n_markers", args.n_markers), ("descend_agree", args.descend_agree),
                                ("gate_overlap", args.gate_overlap), ("gate_score_r", args.gate_score_r)]
              if v is not None}
        annotate(adata, load(args.dag), args.gene_sets, args.clusters, overrides, **kw)
        adata.write_h5ad(args.out)
        print(f"wrote {args.out}")

    elif args.cmd == "substates":
        import scanpy as sc

        from treeline.annotate import add_substates

        adata = sc.read_h5ad(args.h5ad)
        add_substates(adata, args.clusters)
        adata.write_h5ad(args.out)
        print(f"wrote {args.out}")

    elif args.cmd == "integrate":
        import scanpy as sc

        from treeline.scanvi import integrate

        kw = {k: v for k, v in [("supervise_depth", args.supervise_depth), ("n_latent", args.n_latent),
                                ("classification_ratio", args.classification_ratio)]
              if v is not None}
        joint = integrate({Path(f).stem.removesuffix("_annotated"): sc.read_h5ad(f) for f in args.h5ads}, **kw)
        joint.write_h5ad(args.out)
        print(f"wrote {args.out} — recluster obsm['X_treeline'] and resubmit through annotate")

    elif args.cmd == "colors":
        import scanpy as sc

        from treeline.annotate import annotations
        from treeline.colors import palette

        anns = [annotations(sc.read_h5ad(f)) for f in args.h5ads]
        Path(args.out).write_text(json.dumps(palette(*anns), indent=1))
        print(f"wrote {args.out}")

    elif args.cmd == "summary":
        import scanpy as sc

        from treeline.annotate import annotations

        a = sc.read_h5ad(args.h5ad, backed="r")  # obs + uns only; X stays on disk
        ann = annotations(a)
        for key in ann["cluster_keys"]:
            sizes = a.obs[key].astype(str).value_counts()
            calls = ann["calls"][key]
            print(f"\n{key} — {len(calls)} clusters")
            for cl in sorted(calls, key=lambda c: -sizes.get(c, 0)):
                c = calls[cl]
                chain = " > ".join(f"{lv['node']} {lv['share'] * 100:.0f}%" for lv in c["levels"]) or "Unknown"
                sfx = ann["suffixes"].get(key, {}).get(cl)
                extra = ""
                if sfx:
                    combo = "+".join(sfx["markers"])
                    extra = f"  [{combo}+ Fbeta {sfx['fbeta']} PPV {sfx['ppv']} recall {sfx['recall']}]"
                if c["refused"]:
                    extra += f"\n{'':>13}gate: {c['refused']}"
                if c["overridden"]:
                    extra += f"\n{'':>13}override: {c['overridden']}"
                print(f"  {cl:>3} {sizes.get(cl, 0):>7,}  {chain}{extra}")


if __name__ == "__main__":
    main()
