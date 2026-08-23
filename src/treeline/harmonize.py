"""Cross-sample harmonization: CL descent paths ARE the shared vocabulary.

Per-sample cluster IDs never need to match — their paths do. Harmonization is prefix
grouping: at display depth k, a cluster's label is its path truncated to k (a cluster
that stopped shallower keeps its deepest label). The union of observed paths, as a
tree, is the harmonized label space.
"""

from __future__ import annotations


def observed_paths(calls: dict) -> list[tuple[str, ...]]:
    """Every distinct path across {sample: {res: {cluster: {'path': [...]}}}}."""
    seen: dict[tuple[str, ...], None] = {}
    for by_res in calls.values():
        for by_cluster in by_res.values():
            for call in by_cluster.values():
                seen[tuple(call["path"])] = None
    return list(seen)


def path_tree(paths: list[tuple[str, ...]]) -> dict:
    """Union of paths as nested {label: subtree} — the harmonized label space."""
    root: dict = {}
    for p in paths:
        node = root
        for label in p:
            node = node.setdefault(label, {})
    return root
