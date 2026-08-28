"""treeline — multi-level cell type annotation, evidence-gated descent of the Cell Ontology.

One function per verb, all on this namespace:

    import treeline as tln
    dag = tln.derive_tree("gene_sets.csv")            # once; network
    tln.annotate(adata, dag, "gene_sets.csv", ["leiden_1.0"])
    tln.summary(adata)
    tln.add_substates(adata, ["leiden_1.0"])
    joint = tln.integrate({"a": a, "b": b})
    tln.report({"a": a, "b": b}, "report.html")
"""

from treeline.annotate import add_substates, annotate, annotations, coarse_labels, summary
from treeline.colors import palette
from treeline.report import render_report as report
from treeline.scanvi import integrate
from treeline.tree import derive_tree, subdag
from treeline.tree import load as load_tree
from treeline.vote import load_overrides

__all__ = [
    "add_substates", "annotate", "annotations", "coarse_labels", "derive_tree", "integrate",
    "load_overrides", "load_tree", "palette", "report", "subdag", "summary",
]
