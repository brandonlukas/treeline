"""Toy-DAG checks for the vote: descent, shared-descendant votes, gate, Unknown."""

import numpy as np
import pandas as pd

from treeline.vote import Call, assign_cluster

# root -> A, B, C; A -> A1, A2; M is a multi-parent child of both A and B
DAG = {
    "roots": ["A", "B", "C"],
    "nodes": {
        "A": {"sets": ["A"], "parents": [], "children": ["A1", "A2", "M"]},
        "B": {"sets": ["B"], "parents": [], "children": ["M"]},
        "C": {"sets": ["C"], "parents": [], "children": []},
        "A1": {"sets": ["A1"], "parents": ["A"], "children": []},
        "A2": {"sets": ["A2"], "parents": ["A"], "children": []},
        "M": {"sets": ["M"], "parents": ["A", "B"], "children": []},
    },
}
PANELS = {
    "A": ["g1", "g2"],
    "B": ["g3", "g4"],
    "C": ["ga", "gb"],
    "A1": ["g5", "g6"],
    "A2": ["g7", "g8"],
    "M": ["g9", "g10"],
}


def frame(**cols):
    return pd.DataFrame({f"score_{k}": v for k, v in cols.items()})


def test_descends_to_clear_winner():
    n = 20
    call = assign_cluster(
        DAG, PANELS, frame(A=np.ones(n), B=np.zeros(n), A1=np.full(n, 2.0), A2=np.zeros(n), M=np.zeros(n))
    )
    assert call.path == ["A", "A1"]
    assert call.levels[0].share == 1.0  # A1's 2.0 pools up into A's subtree-max
    assert call.final == "A1"


def test_shared_descendant_votes_for_both_parents():
    n = 10  # M dominates: every nucleus's max is M, which sits under both A and B
    call = assign_cluster(
        DAG, PANELS, frame(A=np.zeros(n), B=np.zeros(n), A1=np.zeros(n), A2=np.zeros(n), M=np.ones(n))
    )
    assert call.levels[0].share == 1.0 and call.levels[0].runner_share == 1.0  # shares sum > 1
    assert call.path[-1] == "M"


def test_root_below_threshold_is_unknown():
    # three-way split, 4/4/4: no root child reaches DESCEND_AGREE
    z = [0.0] * 12
    scores = frame(
        A=[1.0] * 4 + [0.0] * 8,
        B=[0.0] * 4 + [1.0] * 4 + [0.0] * 4,
        C=[0.0] * 8 + [1.0] * 4,
        A1=z,
        A2=z,
        M=[-1.0] * 12,
    )
    call = assign_cluster(DAG, PANELS, scores)
    assert call.final == "Unknown" and call.path == []


def test_gate_refuses_indiscriminable_siblings():
    # both criteria trip: identical exclusive panels AND correlated per-nucleus scores
    panels = dict(PANELS, A1=["g5", "g6"], A2=["g5", "g6"])
    n = 20
    a1 = np.linspace(1.0, 2.0, n)
    call = assign_cluster(
        DAG, panels, frame(A=np.ones(n), B=np.zeros(n), A1=a1, A2=a1 * 0.95, M=np.zeros(n))
    )
    assert call.path == ["A"]
    assert call.refused and "share 2/2 exclusive genes" in call.refused


def test_gate_needs_both_criteria():
    # overlapping panels but uncorrelated scores -> descent allowed (the SMC lesson)
    panels = dict(PANELS, A1=["g5", "g6"], A2=["g5", "g6"])
    rng = np.random.default_rng(0)
    n = 200
    a1 = np.concatenate([np.full(150, 2.0), np.full(50, 0.0)]) + rng.normal(0, 0.01, n)
    a2 = rng.normal(1.0, 0.5, n)
    call = assign_cluster(DAG, panels, frame(A=np.ones(n), B=np.zeros(n), A1=a1, A2=a2, M=np.zeros(n)))
    assert call.refused is None and len(call.path) == 2


def test_override_stop():
    n = 20
    ov = [{"node": "A", "decision": "stop", "justification": "j", "author": "a", "date": "d"}]
    call = assign_cluster(
        DAG, PANELS, frame(A=np.ones(n), B=np.zeros(n), A1=np.full(n, 2.0), A2=np.zeros(n), M=np.zeros(n)), ov
    )
    assert call.path == ["A"]
    assert call.overridden and "override (a, d): j" in call.overridden


def test_containment_exempt_from_gate():
    # at the root, A's descendants contain B's shared child M -> no gate between A and B
    n = 20
    call = assign_cluster(
        DAG, PANELS, frame(A=np.ones(n), B=np.full(n, 0.5), A1=np.zeros(n), A2=np.zeros(n), M=np.full(n, 0.9))
    )
    assert isinstance(call, Call) and call.path[0] == "A"
