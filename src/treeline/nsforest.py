"""Name what grows above the treeline: {gene}+ suffixes for subclusters that collapse
under one deepest-supportable CL label.

This is NS-Forest v4.0 (Liu et al., BMC Methods 2024, s44330-024-00015-2) run
one-vs-rest *within* each collapsed group, reimplemented on sklearn:

  1. candidate genes = expressed in >50% of nuclei of >=1 group cluster (positive medians)
  2. Binary Expression Score per (gene, target):  sum_i!=T max(0, 1 - m_gi/m_gT) / (n-1)
  3. BinaryFirst_high pre-selection: score >= mean + 2*std of the group's score distribution
  4. random forest one-vs-rest -> top genes by Gini importance, re-ranked by Binary Score
  5. per-gene single-split decision tree -> optimal expression threshold
  6. all AND-combinations of the top candidates -> best F-beta (beta=0.5, precision-weighted)

The chosen combination is the cluster's minimal marker set; the suffix is its most
binary member. Computed per clustering, one-vs-rest against siblings only — the
classification task a substate name actually has to win. Called by `annotate`; call
`suffixes_for` directly for the full set of knobs.

# ponytail: 100 trees, top 15 by Gini, top 6 into the combination search — smaller than
# the paper's defaults; raise if marker quality on bigger groups disappoints.
"""

from __future__ import annotations

import re
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import fbeta_score
from sklearn.tree import DecisionTreeClassifier

N_TREES = 100
TOP_GINI = 15
TOP_COMBO = 6
BETA = 0.5
MIN_CLUSTER = 20

# structural/housekeeping genes and unnamed ids make bad marker names; the minimal set
# IS the substate's name now, so they are excluded from candidacy entirely
BAD_NAME = re.compile(r"^(MT-|RPL|RPS|MRPL|MRPS|ENSG|MALAT1$|NEAT1$|ACTB$|TMSB4X$|B2M$|EEF1A1$)")


def binary_scores(medians: pd.DataFrame) -> pd.DataFrame:
    """Binary Expression Score, clusters x genes: 1 = expressed in the target only."""
    n = len(medians)
    out = pd.DataFrame(0.0, index=medians.index, columns=medians.columns)
    for target in medians.index:
        m_t = medians.loc[target]
        others = medians.drop(index=target)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = 1 - others.div(m_t, axis=1)
        score = ratio.clip(lower=0).sum(axis=0) / (n - 1)
        out.loc[target] = score.where(m_t > 0, 0.0)
    return out


def nsforest_group(
    x: np.ndarray, genes: np.ndarray, clusters: pd.Series, *,
    n_trees: int = N_TREES, top_gini: int = TOP_GINI, top_combo: int = TOP_COMBO, beta: float = BETA,
) -> dict[str, dict]:
    """NS-Forest over one collapsed group. x: nuclei x candidate genes (dense)."""
    ids = sorted(clusters.unique())
    medians = pd.DataFrame(
        {g: 0.0 for g in genes} | {},
        index=ids,
    )
    for cl in ids:
        medians.loc[cl] = np.median(x[(clusters == cl).values], axis=0)
    scores = binary_scores(medians)
    flat = scores.values.ravel()
    cutoff = flat.mean() + 2 * flat.std()  # BinaryFirst_high

    results: dict[str, dict] = {}
    for cl in ids:
        y = (clusters == cl).values
        s = scores.loc[cl]
        cand = np.where(s.values >= cutoff)[0]
        if len(cand) < top_gini:  # small groups can have coarse score distributions
            cand = np.argsort(s.values)[::-1][:top_gini]
        rf = RandomForestClassifier(n_estimators=n_trees, n_jobs=-1, random_state=0, class_weight="balanced")
        rf.fit(x[:, cand], y)
        by_gini = cand[np.argsort(rf.feature_importances_)[::-1][:top_gini]]
        top = sorted(by_gini, key=lambda j: -s.values[j])[:top_combo]

        thresholds = {}
        for j in top:
            dt = DecisionTreeClassifier(max_depth=1, random_state=0).fit(x[:, [j]], y)
            thresholds[j] = dt.tree_.threshold[0] if dt.tree_.node_count > 1 else 0.0
        best, best_f = (), -1.0
        for k in range(1, len(top) + 1):
            for combo in combinations(top, k):  # ties prefer the smaller combo
                pred = np.all([x[:, j] > thresholds[j] for j in combo], axis=0)
                f = fbeta_score(y, pred, beta=beta, zero_division=0)
                if f > best_f:
                    best, best_f = combo, f
        # the evidence behind the name, NS-Forest v4's own metric set — never secret
        pred = np.all([x[:, j] > thresholds[j] for j in best], axis=0)
        tp = float((pred & y).sum())
        ordered = sorted(best, key=lambda j: -s.values[j])
        results[str(cl)] = {
            "markers": [str(genes[j]) for j in ordered],
            "fbeta": round(float(best_f), 3),
            "ppv": round(tp / pred.sum(), 3) if pred.sum() else 0.0,
            "recall": round(tp / y.sum(), 3) if y.sum() else 0.0,
            "marker_stats": [
                {
                    "gene": str(genes[j]),
                    "binary_score": round(float(s.values[j]), 3),
                    "threshold": round(float(thresholds[j]), 3),
                    # on-target fraction: how much of this gene's expression sits in the
                    # target cluster (1 = expressed there only) — within the group
                    "on_target": round(float(x[y, j].sum() / x[:, j].sum()), 3) if x[:, j].sum() else 0.0,
                }
                for j in ordered
            ],
        }
    return results


def suffixes_for(
    adata, clusters: pd.Series, calls: dict[str, dict], *,
    min_cluster: int = MIN_CLUSTER, n_trees: int = N_TREES, top_gini: int = TOP_GINI,
    top_combo: int = TOP_COMBO, beta: float = BETA,
) -> dict[str, dict]:
    """cluster id -> {gene, markers, fbeta} for clusters sharing a path with a sibling."""
    by_path: dict[tuple, list[str]] = defaultdict(list)
    for cl, call in calls.items():
        if call["path"]:
            by_path[tuple(call["path"])].append(cl)
    out: dict[str, dict] = {}
    for group in by_path.values():
        counts = clusters.value_counts()
        keep = [c for c in group if counts.get(c, 0) >= min_cluster]
        if len(keep) < 2:
            continue
        mask = clusters.isin(keep).values
        sub, subcl = adata[mask], clusters[mask]
        # candidates: nonzero in >50% of nuclei of some member cluster (positive median),
        # excluding genes that make bad names (the minimal set is the name)
        nz = pd.DataFrame.sparse.from_spmatrix(sub.X > 0, columns=sub.var_names).groupby(subcl.values).mean()
        namable = np.array([not BAD_NAME.match(g) for g in sub.var_names])
        cand_mask = (nz > 0.5).any(axis=0).values & namable
        x = np.asarray(sub[:, cand_mask].X.todense())
        genes = sub.var_names[cand_mask].to_numpy()
        group_res = nsforest_group(
            x, genes, subcl.reset_index(drop=True),
            n_trees=n_trees, top_gini=top_gini, top_combo=top_combo, beta=beta,
        )
        used: set[str] = set()
        for cl in keep:
            r = group_res[cl]
            name = next((g for g in r["markers"] if not BAD_NAME.match(g) and g not in used), None)
            if name is None:
                continue
            used.add(name)
            out[cl] = {"gene": name, **r}
    return out
