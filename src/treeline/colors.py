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


def _tree_of(paths: list[tuple[str, ...]]) -> dict:
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
    return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))


def assign(paths: list[tuple[str, ...]]) -> dict[str, str]:
    """Path (and every prefix of it) -> hex color. Siblings spread across the wheel at
    the first branching; below it children keep close to the parent's hue, and lightness
    darkens one step per level — subclusters of one type read as shades of one hue."""
    out: dict[str, str] = {"Unknown": UNKNOWN}
    tree = _tree_of([p for p in paths if p])

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
