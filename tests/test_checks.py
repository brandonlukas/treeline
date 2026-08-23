"""Input-contract checks fail loudly, not deep in a stack."""

import anndata as ad
import numpy as np
import pytest

from treeline.annotate import _check_lognorm
from treeline.scanvi import _check_annotated_counts
from treeline.vote import load_panels


def test_raw_counts_refused_by_annotate():
    raw = ad.AnnData(np.array([[0.0, 500.0], [3.0, 1.0]]))
    with pytest.raises(ValueError, match="raw counts"):
        _check_lognorm(raw)
    _check_lognorm(ad.AnnData(np.log1p(raw.X)))  # log-normalized passes


def test_integrate_refuses_unprepared_inputs():
    a = ad.AnnData(np.log1p(np.array([[0.0, 5.0], [3.0, 1.0]])))
    with pytest.raises(ValueError, match="not annotated"):
        _check_annotated_counts("s1", a)
    a.uns["treeline"] = "{}"
    a.obs["treeline_leiden"] = ["x", "y"]
    with pytest.raises(ValueError, match="layers\\['counts'\\]"):
        _check_annotated_counts("s1", a)
    a.layers["counts"] = np.array([[0.0, 5.0], [3.0, 1.5]])  # fractional
    with pytest.raises(ValueError, match="not raw counts"):
        _check_annotated_counts("s1", a)
    a.layers["counts"] = np.array([[0.0, 5.0], [3.0, 1.0]])
    _check_annotated_counts("s1", a)  # passes


def test_gene_set_table_schema_error(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("set,gene\nfoo,ACTA2\n")
    with pytest.raises(ValueError, match="Gene Set Name"):
        load_panels(bad, {"nodes": {}})
