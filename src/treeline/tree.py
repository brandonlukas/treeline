"""Derive the Cell Ontology DAG over a gene-set table's CL term labels, via EBI OLS4.

Unlike the matchafinn-apps prototype (derive_cl_tree.py): the DAG is preserved — every
minimal is_a parent within the node set is kept, no pins, no tie-breaks. Multi-parent
nodes (pericyte) legitimately appear under all their parents; the vote handles them
(subtree-max is well-defined over descendant sets). Single-child chains still collapse,
which auto-merges near-synonyms (myometrial cell folds into uterine smooth muscle cell).

Live API at derivation, frozen JSON at inference:

    treeline derive-tree <gene_sets.csv> <out.json>
"""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path
from typing import Any

OLS = "https://www.ebi.ac.uk/ols4/api"


def _get(url: str) -> Any:
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.load(r)


def resolve(label: str) -> tuple[str, str]:
    """CL term label -> (IRI, CL id): exact label match, else exact synonym match
    (gene-set tables carry older CL labels that the ontology has since renamed)."""
    q = urllib.parse.urlencode(
        {"q": label, "ontology": "cl", "exact": "true", "rows": 10, "fieldList": "label,obo_id,iri,synonym"}
    )
    docs = [d for d in _get(f"{OLS}/search?{q}")["response"]["docs"] if d.get("obo_id", "").startswith("CL:")]
    for doc in docs:
        if doc["label"].lower() == label.lower():
            return str(doc["iri"]), str(doc["obo_id"])
    for doc in docs:
        if label.lower() in (s.lower() for s in doc.get("synonym", [])):
            print(f"[resolve] {label!r} -> current CL label {doc['label']!r} (synonym match)", file=sys.stderr)
            return str(doc["iri"]), str(doc["obo_id"])
    raise ValueError(f"no exact CL match for {label!r}")


def ancestors(iri: str) -> list[str]:
    """is_a ancestor labels (lowercase), CL terms only."""
    enc = urllib.parse.quote(urllib.parse.quote(iri, safe=""), safe="")
    terms = _get(f"{OLS}/ontologies/cl/terms/{enc}/ancestors?size=500").get("_embedded", {}).get("terms", [])
    return [str(t["label"]).lower() for t in terms if str(t.get("obo_id") or "").startswith("CL:")]


def set_names_from_csv(csv_path: str | Path) -> list[str]:
    """Ordered unique gene-set names ('<label> - marker genes' -> '<label>')."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    if "Gene Set Name" not in df.columns:
        raise ValueError(
            f"gene-set table {csv_path} has no 'Gene Set Name' column; found {list(df.columns)}. "
            "Expected a CellGuide-style CSV (set name on each set's first row, blank rows fill down)."
        )
    names = df["Gene Set Name"].ffill().str.removesuffix(" - marker genes")
    return list(dict.fromkeys(names))


def derive(set_names: list[str]) -> dict:
    """Build the DAG: nodes = set-bearing terms + CL ancestors shared by >=2 of them."""
    anc: dict[str, set[str]] = {}
    cl_id: dict[str, str] = {}
    dropped: list[str] = []
    for name in set_names:
        try:
            iri, cl_id[name.lower()] = resolve(name)
        except ValueError as e:
            print(f"[derive] set dropped, {e}", file=sys.stderr)
            dropped.append(name)
            continue
        anc[name.lower()] = set(ancestors(iri))
        print(f"{name}: {len(anc[name.lower()])} ancestors", file=sys.stderr)
    set_names = [s for s in set_names if s not in dropped]
    names = set(anc)

    shared = {a for a in set().union(*anc.values()) if sum(a in v for v in anc.values()) >= 2}
    for n in sorted(shared - names):
        iri, cl_id[n] = resolve(n)
        anc[n] = set(ancestors(iri))
        print(f"[grouping] {n}: {len(anc[n])} ancestors", file=sys.stderr)
    nodes = names | shared

    def parents_of(n: str) -> list[str]:
        cand = [a for a in anc[n] if a in nodes]
        return sorted(a for a in cand if not any(a in anc[b] for b in cand if b != a))

    parents = {n: parents_of(n) for n in nodes}

    def children_of() -> dict[str, list[str]]:
        ch: dict[str, list[str]] = {n: [] for n in nodes}
        for n, ps in parents.items():
            for p in ps:
                ch[p].append(n)
        return {n: sorted(c) for n, c in ch.items()}

    sets = {n: [s for s in set_names if s.lower() == n] for n in nodes}

    # splice barren grouping nodes (set-less, <2 children): parents adopt the children
    changed = True
    while changed:
        changed = False
        children = children_of()
        for n in sorted(nodes):
            if not sets[n] and len(children[n]) < 2:
                nodes.discard(n)
                for c in children[n]:
                    parents[c] = sorted((set(parents[c]) - {n}) | set(parents[n]))
                del parents[n], sets[n]
                changed = True
                break  # recompute children before the next splice

    # collapse single-child chains (child has this sole parent): parent absorbs the
    # child's sets and children — the near-synonym auto-merge; keeps the label with
    # sets, preferring the deeper (child) label when the parent is set-less
    changed = True
    while changed:
        changed = False
        children = children_of()
        for n in sorted(nodes):
            kids = children.get(n, [])
            if len(kids) == 1 and parents[kids[0]] == [n]:
                c = kids[0]
                label = n if sets[n] else c
                merged_sets = sets[n] + sets[c]
                grandkids = children_of().get(c, [])
                nodes.discard(n)
                nodes.discard(c)
                nodes.add(label)
                sets.pop(n), sets.pop(c, None)
                keep_parents = parents.pop(n)
                parents.pop(c, None)
                sets[label] = merged_sets
                parents[label] = keep_parents
                for g in grandkids:
                    parents[g] = sorted((set(parents[g]) - {c}) | {label})
                changed = True
                break

    # transitive reduction: splicing unions parents without re-checking minimality,
    # so drop any parent that is an ancestor of another parent (redundant edge)
    def dag_ancestors(n: str, seen: set[str] | None = None) -> set[str]:
        seen = seen if seen is not None else set()
        for p in parents.get(n, []):
            if p not in seen:
                seen.add(p)
                dag_ancestors(p, seen)
        return seen

    for n in sorted(nodes):
        ps = parents[n]
        parents[n] = sorted(p for p in ps if not any(p in dag_ancestors(q) for q in ps if q != p))

    children = children_of()
    roots = sorted(n for n in nodes if not parents[n])
    return {
        "source": "EBI OLS4, ontology=cl",
        "fetched": date.today().isoformat(),  # noqa: DTZ011 — human-readable provenance date
        "derived_from_sets": set_names,
        "dropped_sets": dropped,
        "roots": roots,
        "nodes": {
            n: {"sets": sets[n], "parents": parents[n], "children": children[n], "cl_id": cl_id.get(n)}
            for n in sorted(nodes)
        },
    }


def load(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())


def descendants(dag: dict, label: str) -> set[str]:
    """A node's descendant set, itself included. DAG-safe."""
    out, stack = set(), [label]
    while stack:
        n = stack.pop()
        if n not in out:
            out.add(n)
            stack += dag["nodes"][n]["children"]
    return out


def subdag(dag: dict, root: str) -> dict:
    """The DAG restricted to `root` and its descendants (parents trimmed to match)."""
    keep = descendants(dag, root)
    nodes = {
        n: {**dag["nodes"][n], "parents": [p for p in dag["nodes"][n]["parents"] if p in keep]}
        for n in keep
    }
    return {**dag, "roots": [root], "nodes": nodes}
