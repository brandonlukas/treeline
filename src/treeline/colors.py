"""Hierarchical colors: labels under the same parent share a hue family.

The observed descent paths form a tree; each top-level branch owns an arc of the hue
wheel, children subdivide their parent's arc proportionally to leaf count, and depth
darkens the shade — so subclusters of one cell type read as shades of one hue.
Pass-through levels (single child) inherit the full arc without burning hue space.
Colors are keyed by full path (" > "-joined), so a multi-parent CL term reached through
two different parents gets each parent's hue — the duplication stays visible, as the
DAG-on-a-slider rule in SPECS requires.
"""

from __future__ import annotations

import colorsys

UNKNOWN = "#9a9a9a"
SATURATION = 0.58
L_START, L_STEP, L_MIN = 0.66, 0.08, 0.30  # lightness darkens with path depth
SHRINK = 0.4  # below the first branching, children stay within this fraction of the parent's arc


def observed_paths(calls: dict) -> list[tuple[str, ...]]:
    """Every distinct path across {sample: {res: {cluster: {'path': [...]}}}}. Per-sample
    cluster ids never need to match — their CL paths are the shared vocabulary."""
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


def _leaves(node: dict) -> int:
    return sum(_leaves(c) for c in node.values()) if node else 1


def _hex(h: float, lightness: float) -> str:
    r, g, b = colorsys.hls_to_rgb((h % 360) / 360, lightness, SATURATION)
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"


def palette(*annotation_dicts: dict) -> dict[str, str]:
    """The colors verb: one or more `annotate` outputs (`annotations(adata)`) -> label
    hex map, covering every observed path prefix and every `{path} ⊕ {gene}` substate
    (substates ramp within the parent hue: poline-style saturation/lightness curve)."""
    paths: dict[tuple[str, ...], None] = {}
    for ann in annotation_dicts:
        for by_cl in ann["calls"].values():
            for call in by_cl.values():
                paths[tuple(call["path"])] = None
    out = assign(list(paths))
    by_parent: dict[str, list[str]] = {}
    for ann in annotation_dicts:
        for key, sfx in ann.get("suffixes", {}).items():
            for cl, e in sfx.items():
                parent = " > ".join(ann["calls"][key][cl]["path"])
                k = f"{parent} ⊕ {e['gene']}"
                by_parent.setdefault(parent, [])
                if k not in by_parent[parent]:
                    by_parent[parent].append(k)
    for parent, keys in by_parent.items():
        for i, k in enumerate(sorted(keys)):
            out[k] = _substate_hex(out[parent], i, len(keys))
    return out


def _substate_hex(parent_hex: str, i: int, k: int) -> str:
    r, g, b = (int(parent_hex[j : j + 2], 16) / 255 for j in (1, 3, 5))
    h, _, _ = colorsys.rgb_to_hls(r, g, b)
    t = i / (k - 1) if k > 1 else 0.5
    r, g, b = colorsys.hls_to_rgb(h, 0.64 - 0.40 * t, 0.42 + 0.38 * t)
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"


def assign(paths: list[tuple[str, ...]]) -> dict[str, str]:
    """Path (and every prefix of it) -> hex color. Siblings spread across the wheel at
    the first branching; below it children keep close to the parent's hue, and lightness
    darkens one step per level — subclusters of one type read as shades of one hue."""
    out: dict[str, str] = {"Unknown": UNKNOWN}
    tree = path_tree([p for p in paths if p])

    def walk(node: dict, prefix: tuple[str, ...], mid: float, width: float, spread: int) -> None:
        total = sum(_leaves(c) for c in node.values()) or 1
        use = width if spread > 0 else width * SHRINK
        pos = mid - use / 2
        branching = len(node) > 1
        for label, child in node.items():
            w = use * _leaves(child) / total
            cmid = pos + w / 2
            p = prefix + (label,)
            lightness = max(L_MIN, L_START - L_STEP * (len(p) - 1))
            out[" > ".join(p)] = _hex(cmid, lightness)
            walk(child, p, cmid, w, spread - 1 if branching else spread)
            pos += w

    # full hue spread for the first two real branchings (the major families), then
    # children hold close to the parent hue and depth darkens the shade
    walk(tree, (), 180.0, 360.0, 2)
    return out
