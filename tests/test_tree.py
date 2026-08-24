"""Toy checks for the DAG/label helpers: subdag restriction and the coarse-label depth cut."""

import anndata as ad
import numpy as np
import pandas as pd

from treeline.annotate import coarse_labels
from treeline.tree import subdag

# root -> A, B; A -> A1; M is a multi-parent child of both A and B
DAG = {
    "roots": ["A", "B"],
    "fetched": "2026-01-01",
    "nodes": {
        "A": {"sets": ["A"], "parents": [], "children": ["A1", "M"]},
        "B": {"sets": ["B"], "parents": [], "children": ["M"]},
        "A1": {"sets": ["A1"], "parents": ["A"], "children": []},
        "M": {"sets": ["M"], "parents": ["A", "B"], "children": []},
    },
}


def test_subdag_keeps_descendants_and_trims_outside_parents():
    sub = subdag(DAG, "A")
    assert sub["roots"] == ["A"]
    assert set(sub["nodes"]) == {"A", "A1", "M"}
    assert sub["nodes"]["M"]["parents"] == ["A"]  # B is outside the sub-DAG
    assert sub["fetched"] == DAG["fetched"]  # non-node fields carried through
    assert DAG["nodes"]["M"]["parents"] == ["A", "B"]  # original untouched


def test_subdag_of_leaf_is_singleton():
    sub = subdag(DAG, "A1")
    assert sub["roots"] == ["A1"]
    assert set(sub["nodes"]) == {"A1"}


def test_coarse_labels_cuts_below_root_and_keeps_root_stops():
    a = ad.AnnData(np.zeros((3, 1)))
    a.obs["treeline_leiden"] = pd.Categorical(["cell > A > A1", "cell", "Unknown"])
    assert list(coarse_labels(a, "leiden")) == ["A", "cell", "Unknown"]
