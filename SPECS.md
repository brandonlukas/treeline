# treeline — technical specification (seed)

*The treeline is the elevation past which conditions stop supporting growth. This tool
descends the Cell Ontology only as far as the evidence supports, then stops.*

Multi-level cell type annotation + label-aware batch correction for single-cell/nucleus
RNA data. Annotation is treated as **structured prediction over the Cell Ontology (CL)**,
not flat classification: a cluster's label is a *path* through the ontology whose depth is
earned level by level, with a slider from coarse to fine. Data-driven suffixes (NS-Forest
minimal markers) name what grows *above* the treeline — substates the ontology has no
term for yet.

Status: **proof of concept**. Target: make it work on the one dataset in `assets/`.
Comprehensive benchmarks, other datasets, other ontologies: explicit non-goals for now.

## Scope: one prep step, three capabilities

1. **`derive-tree`** (prep, once per gene-set table): CL-labeled gene sets -> the frozen,
   dated CL DAG via OLS4. Explicit — never folded into annotate — because the dated
   artifact pins provenance (live API at derivation, frozen artifact at inference).
2. **`annotate`**: scRNA-seq expression + cluster labels (one or more resolutions, any
   clustering method) + gene sets + frozen DAG -> multi-level annotation, as descriptive
   as the ontology and marker sets support, falling back to NS-Forest `{gene}+` suffixes
   above the treeline. Results are written back into the AnnData (per-cluster calls,
   shares, refusals, suffixes in `.uns["treeline"]`; per-nucleus labels in `.obs`), so an
   annotated `.h5ad` is self-contained.
3. **`integrate`**: two or more *annotated* AnnDatas -> scANVI cell-type-aware
   semi-supervised integration. Supervision is the tree-cut consensus prior (multiple
   supplied resolutions strengthen it; cross-resolution disagreement at the cut ->
   `Unknown`, scANVI's unlabeled class). Emits the joint latent and stops.
4. **`colors`**: annotation output -> the hierarchical palette (same parent, same hue;
   depth darkens; substates ramp within the parent hue).

**Treeline performs no clustering anywhere.** The user's loop: cluster (their method) ->
annotate per sample -> integrate -> recluster the latent (their method) -> annotate the
joint -> colors at every step. Not in the public API: upstream QC/ambient correction,
within-type reclustering (`--refine` remains an opt-in convenience), and the HTML report —
the report is a downstream *consumer* of these verbs (its natural views: one sample,
samples side by side, integrated, integrated + relabels), kept in `apps/` for development
and to be rebuilt against the stable API later. **Not treeline's job by default:** upstream QC and
ambient correction, per-sample clustering, and within-cell-type reclustering — subtype
discovery beyond the supplied clusterings is the user's discretion (`--refine` ships as
an opt-in convenience because the latent and labels are already in hand).

## Why (the gap)

Existing hierarchical/ontology-aware annotators (CellO, OnClass, scANVI-based transfer,
SOCAM, GPTAnno) are all **reference/supervision-dependent** — the ontology structures
their output space but a trained classifier or reference atlas supplies the evidence.
When no trusted reference exists for a tissue (the motivating case: uterine
leiomyoma/myometrium snMultiome), you fall back to marker sets — and flat marker voting
fails in a characteristic way: **sibling panels split votes** (pericyte vs SMC, T vs NK
vs macrophage), producing confidently-wrong leaf labels. treeline's bet: marker-set
evidence + ontology structure + evidence-gated descent is enough for trustworthy
annotation without any trained reference, and the labels are then good enough to
semi-supervise integration (scANVI-style).

## Inputs (the tool's whole contract)

1. **Per-sample AnnData** (`.h5ad`), ambient-corrected counts, QC'd, with one or more
   categorical `.obs` columns naming clusterings. Clustering is upstream, like QC:
   treeline is agnostic to the algorithm — the vote consumes cluster assignments and
   nothing else. (The POC driver runs its own Leiden at two resolutions because the
   assets h5ads ship unclustered; that is the driver's business, not the tool's
   contract.)
2. **Gene set table** (CSV), user-supplied at inference time: gene sets labeled with CL
   term names. CellGuide's per-tissue "export gene set annotations" satisfies this
   directly; any custom table with CL-labeled sets works. Supplied as a *file*, never
   fetched live — a dated artifact pins provenance; API corpora drift.
3. **CL structure**, fetched once via the EBI OLS4 API from the gene-set labels and
   frozen to JSON (`treeline derive-tree`). Live API at derivation, frozen artifact at
   inference: offline runs, reproducible labels, versionable curation.
4. **Overrides config** (see the provenance contract below) — expert decisions only,
   every entry `{node, decision, justification, author, date}`, justification enforced
   by schema and carried verbatim into every output artifact. `decision` is a closed
   vocabulary — `prune` (drop node + subtree from the candidate space), `stop` (cap
   descent at this node), `relabel` (override a cluster's assigned path) — so overrides
   stay domain knowledge, not an escape hatch; extending the vocabulary means editing
   the schema, on purpose.

## The provenance contract (core design principle)

Every influence on a label belongs to exactly one of three tiers; nothing silent:

1. **Ontology facts, verbatim — including the ambiguity.** CL is a DAG; treeline keeps
   it a DAG rather than tie-breaking it into a tree. The vote is DAG-native: subtree-max
   is well-defined over a node's descendant set regardless of parent count, so a node
   with two parents legitimately contributes evidence under both, and descent follows
   whichever pools stronger. (This dissolves most of the "pins" the matchafinn-apps
   prototype needed — they existed only to force a tree out of a DAG.)
2. **Automated, self-reporting rules.** Data-driven decisions are parameterized rules
   that print their reasoning, never hand edits. Chief example, the **discriminability
   gate**: refuse to offer a descent whose sibling sets cannot be told apart — gene-set
   overlap ≥ `GATE_OVERLAP` (POC 0.5) **and** per-nucleus score correlation ≥
   `GATE_SCORE_R` (POC 0.9); both must trip, and the refusal reports both numbers.
   (Measured on 1619: overlap alone wrongly refused the uterine-SMC descent — 5/10
   shared genes but score r≈0.6, the scores discriminate fine.) The vote-share descent
   threshold is the same species of rule; all three live as loud module-level constants
   next to `DESCEND_AGREE`.
3. **Expert overrides, signed.** Whatever remains is domain knowledge in the overrides
   config, with mandatory justification, propagated into outputs. No anonymous curation.
   A signature means a *human domain expert*. An override proposed by an AI assistant
   (or carried from AI-assisted prototype code) is not expert knowledge: it must say so
   in its author field, be marked PROVISIONAL in its justification, separate measured
   facts from asserted premises, and stay flagged until a human expert reviews it. The
   data can verify a panel's contents; only a biologist can verify what the panel's
   absences mean.

The matchafinn-apps prototype fails this standard in two known ways — silent tie-breaks
in tree induction, and pins compensating for algorithm weakness rather than expressing
domain knowledge. Fixing both is a primary goal of the POC, not a v2 nicety.

## Pipeline

```
per-sample AnnData ─┬─> user-supplied clusterings from .obs (any algorithm; the POC
                    │   driver adds Leiden at R=2 resolutions, e.g. 0.5 / 2.0)
gene sets + tree ───┴─> per-cluster hierarchical vote  ──> multi-level labels per
                        (subtree-max scores, descend        (sample, resolution, cluster)
                         while vote ≥ threshold)                    │
                                                                    v
                        cross-sample harmonization: CL terms ARE the shared vocabulary —
                        per-sample cluster IDs never need to match, their paths do
                                                                    │
                                                                    v
                        (optional) one scANVI round: coarse labels as semi-supervision,
                        "Unknown" for below-threshold clusters ──> integrated latent
                                                                    │
                                                                    v
                        (optional) NS-Forest minimal markers per subcluster within a
                        coarse label ──> "{gene}+ <coarse label>" suffixes
                                                                    │
                                                                    v
                        static HTML report: the slider (coarse -> fine), per-level labels,
                        agreements, DE evidence, sample composition per cluster
```

### The vote (validated in matchafinn-apps `apps/fig4_scatac/04_annotate.py`)

Per nucleus, score every node's gene set (`scanpy.tl.score_genes`). Descent is per-node,
not per-"level" (a DAG has no global levels): at the current node, nuclei vote among its
child subtrees by **subtree-max** score (siblings pool rather than split); the cluster's
label descends into the winning child only while its vote share clears `DESCEND_AGREE`
(POC default 0.5). Report the full path + per-level agreement, never just the deepest
label.

**DAG semantics of the vote (pinned, not implied):**

- A nucleus whose max-scoring node lies in the *shared* descendant set of two siblings
  votes for **both**. Vote shares are per-child fractions of nuclei and therefore need
  not sum to 1 when sibling subtrees overlap — expected, reported, not an error. Descent
  goes to the single highest-share child, and only if it clears `DESCEND_AGREE`; if two
  children clear it, descend into the higher and report the runner-up's share.
- A cluster's **label is its descent history** — the sequence of nodes actually descended
  through — not "an ontology path" (a multi-parent node has several). Harmonization and
  the report slider operate on descent histories.
- The root is gated like any other node: a cluster that cannot clear `DESCEND_AGREE` at
  the first descent is labeled `Unknown` (the same token the scANVI round consumes).
  The sibling vote always produces a winner; the threshold is what stops
  confidently-wrong-at-the-root.

### Integration (label-aware, tree-cut supervision)

scANVI's label head is a flat categorical — ontology-blind, so it would separate
`fibroblast` from `stromal cell of ovary` as hard as from `uterine smooth muscle cell`.
The tree fixes this at the input: supervise with descent paths **truncated at
`SUPERVISE_DEPTH`** (POC: 2, taken from the fine clustering), where sibling labels are
as ontologically parallel as CL gets and uncertain deep splits collapse into their
common parent (fibroblast and stromal-of-ovary both become `connective tissue cell`).
Clusters stalled above the cut are `Unknown` — scANVI's native unlabeled class.
Tree-distance-weighted label losses are v2 research. NS-Forest suffixes are **never**
supervision: they are artifacts of one clustering, and a flat classifier would carve
the latent along per-sample Leiden boundaries.

**Labels are a prior, not truth.** scANVI pins labeled cells to their label; treeline
softens this with a **multi-resolution consensus**: a nucleus is supervised only when
its depth-cut label agrees across all clustering resolutions where it resolves;
nuclei that flip between classes get `Unknown` (scANVI's let-the-data-decide class),
and clusters that stall at one resolution inherit their consensus from the others.
Prior strength = cross-resolution agreement. A true soft-label ELBO (per-nucleus label
distributions the model can overrule) is v2 research.

**Integration strength is a chosen prior, not a fittable parameter, on this data.**
Sample ≡ condition here (one LM, one MM), so batch effect and disease effect are the
same axis — unidentifiable. `CLASSIFICATION_RATIO` is the exposed knob (label pull vs
mixing); with n=1 per condition its value is a judgment call. More patients
(batch=patient, condition=tissue) make it identifiable; partial within-type mixing may
be real disease states and is not, by itself, a failure.

After training, treeline emits the latent and stops; the user (POC: the driver) clusters
it and resubmits, and the same vote runs on the joint clusters (marker scores are
expression-derived, so re-annotation is not circular through the latent), then NS-Forest
once on joint clusters — a shared
cross-sample substate vocabulary. Then **within-class refinement (opt-in, `--refine`)**:
for each coarse class with enough nuclei, per-class HVGs and a fresh per-class
integration (scVI, batch=sample) on just those nuclei — within-type substructure hides
in genes global HVG selection never keeps — subclustered and re-voted, NS-Forest per
class; batch correction acts only within a type, so it cannot smear classes together.
Off by default per the scope statement. Caveat stands (see design
decisions): joint clusters are label-influenced; downstream differential analysis on
them inherits the supervision.

### Structure derivation (prototype: `derive_cl_tree.py` in matchafinn-apps — supersede it)

Resolve set labels -> CL IRIs (OLS4 search, exact), fetch is_a ancestor closures
(CL-only), induce the hierarchy over set-bearing terms + shared ancestors, prune barren
generic ancestors, collapse single-child chains (which auto-merges near-synonyms:
myometrial ⊂ uterine SMC). Unlike the prototype: **preserve the DAG** (no pins, no
tie-breaks — see the provenance contract); multi-parent nodes are kept and the vote
handles them. The frozen JSON records the DAG, the OLS4 fetch date, and nothing else.

## Design decisions already paid for (do not relearn)

- **Whole-cell marker sets are weak in nuclei** — cytoplasm-abundant markers (ACTA2-type)
  under-detect; immune panels especially. Pooling up the tree is the mitigation; expect
  coarse labels to be far more reliable than leaves.
- **Subtree-max has a bias**: a node with several correlated sets beats a single-set
  sibling (max of k correlated scores > one score). Known; acceptable at POC; a mean- or
  size-corrected variant is a v2 experiment.
- **Labels must not feed the embedding that gates downstream measurement** when the
  downstream analysis measures the modality used for gating (see matchafinn-apps fig4:
  ATAC held out of cell-state definition). scANVI integration is a *product* of treeline,
  but consumers doing differential analysis on gated populations should be warned in docs.
- **Ambient correction is upstream but decisive** — on uncorrected nuclei counts, ambient
  dominant-cell-type transcripts flip marker votes in minority populations. Assume
  CellBender-corrected input; document loudly.

## Known-hard problems (the actual research content)

1. **Resolutions don't nest.** Clusterings at two resolutions (Leiden or otherwise)
   are not a hierarchy (clusters
   split AND merge). POC: present resolutions independently, no fake nesting. v2 options:
   over-cluster once and merge upward; cluster-the-clusters dendrogram.
2. **Cross-sample harmonization** beyond exact path match (partial-depth matches, one
   sample resolving deeper than another). POC: group by longest common path prefix.
3. **NS-Forest cost per slider position** — markers are per-clustering. POC: compute at
   the fine resolution only, lazily.
4. **Displaying a DAG on a tree slider** — a multi-parent node (pericyte) appears under
   more than one parent. POC: allow the duplication in the report and mark it; do not
   fake a tree.
5. **Generic-restatement sibling panels are invisible to the gate.** Measured on 1619:
   the lymphatic-EC set is generic endothelial genes and *wins* the EC descent
   decisively (share 0.84, overlap 4/10, score r 0.3–0.5 — every automated measure
   reads as a legitimate descent). Refusing it requires knowing the panel lacks the
   markers that define the type (PROX1/CCL21) — domain knowledge. A tier-3 `stop` was
   drafted for it, but its premise is AI-asserted and unreviewed, so it sits in
   `assets/overrides.proposed.json` awaiting a biologist's signature and is NOT active:
   the POC labels follow the inputs (all EC clusters descend to the lymphatic leaf),
   and the disagreement stays visible — one EC subcluster is CCL21+/PROX1+ by measurement,
   the rest are not. The clean resolution is fixing the *input* (a real lymphatic panel
   in the gene-set table), not curation. An automated rule is open research, not assumed.

## POC acceptance test (the one check that matters)

On `assets/` data: per-sample clustering + tree vote reproduces the fig4 findings
without per-dataset code — SMC clusters resolve to `uterine smooth muscle cell` with
high agreement; the pericyte-DE cluster stops at a coarse contractile/SMC level rather
than mislabeling; the mixed immune cluster resolves through `leukocyte`; endothelial
stops at `endothelial cell` (no unsupported lymphatic descent) — NOTE: this last
criterion encodes fig4's hand prune, whose premise is unreviewed (known-hard #5). It is
met only if the proposed override is expert-signed or the input gains a real lymphatic
panel; until then the tool honestly descends, and that deviation is expected, not a bug.

Secondary smoke checks (not evidence for the bet — the vote is; these just prove the
optional stages wire up): scANVI round runs; NS-Forest emits `{gene}+` suffixes for SMC
subclusters; HTML report renders both resolutions. The core test passing with the smoke
checks unbuilt is still a successful POC.

## Repo layout (proposed)

```
treeline/
  SPECS.md                this file
  assets/                 dev/test data (gitignored except small files — see assets/README.md)
  src/treeline/           tree.py (derive/load), vote.py, harmonize.py, scanvi.py,
                          nsforest.py, report.py — flat, no subpackages
  apps/poc_1619.py        the POC driver over assets/
```

Conventions inherited from the sibling repos: `src/`-layout, ruff + mypy + pytest,
helpers-first/`main()`-last in apps, public-first in `src/`.

## Dependencies

scanpy, anndata, scvi-tools (scANVI), stdlib urllib for OLS4. No new frameworks; the
slider report is static HTML + vanilla JS, no server. NS-Forest is reimplemented in
`nsforest.py` on sklearn (already a scanpy dependency) following the v4.0 paper
(Liu et al., BMC Methods 2024): Binary Score -> BinaryFirst pre-selection -> RF Gini ->
per-gene decision-tree thresholds -> F-beta(0.5) combination search — ~80 lines beats
depending on the thin `nsforest` PyPI package, and runs one-vs-rest within a collapsed
group, which the package's whole-dataset workflow doesn't expose.
