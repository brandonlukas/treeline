# assets/ — POC development data

Dev/test inputs for the proof of concept. The h5ads are **gitignored** (too large for
git; 1619MM exceeds GitHub's 100MB limit) — this README is the committed record.
Provenance: all produced by `../matchafinn-apps` `apps/fig4_scatac/` (see that repo for
the exact pipeline and versions).

| file | what it is |
|---|---|
| `1619LM_gex.h5ad` | Patient 1619 leiomyoma snMultiome GEX. CellBender 0.4.0 ambient-corrected counts, QC'd (per-sample MAD filters, both-modality barcode intersection): 7,139 nuclei. From `results/fig4_scatac/qc/`. |
| `1619MM_gex.h5ad` | Same, matched myometrium: 15,883 nuclei. |
| `cellguide_uterus_gene_sets_2026-08-22.csv` | CellGuide computational marker gene sets for uterus (CL-term-labeled), exported by hand from the CELLxGENE CellGuide uterus explorer 2026-08-22. The user-supplied gene-set input. |
| `cl_tree.json` | The *prototype's* tree (matchafinn-apps `derive_cl_tree.py`): OLS4-derived but tie-broken into a tree by pins, endothelial descent pruned by hand. Kept as a comparison baseline — treeline's own `derive-tree` must NOT reproduce it (see SPECS provenance contract: preserve the DAG, no pins). The endothelial prune turned out NOT to be automatable from scores (SPECS known-hard #5) — it lives on as the signed `stop` override in `overrides.json`. |
| `cl_dag.json` | treeline's derived DAG (`python -m treeline.tree`): 29 nodes, multi-parent nodes kept (pericyte under 4 parents), no pins. |
| `overrides.json` | Tier-3 expert overrides, schema-enforced `{node, decision, justification, author, date}`. Currently EMPTY — no human-signed overrides exist. |
| `overrides.proposed.json` | Drafted overrides awaiting expert review; NOT loaded by the pipeline. One entry: AI-proposed `stop` at `endothelial cell` (see SPECS known-hard #5) — move into `overrides.json` only with a real biologist's signature. |

To regenerate the h5ads: run matchafinn-apps `apps/fig4_scatac/00`–`02` (raw FASTQs on
Box `nu_fibroids_snmultiome`; cellranger-arc + CellBender + QC, ~1 day of compute).
