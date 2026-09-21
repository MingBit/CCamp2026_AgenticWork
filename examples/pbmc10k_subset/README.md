# PBMC 10k Multiome subset: clustering to SCENIC+ eRegulons

Outputs of a complete `run` of the workflow (inspection -> QC -> normalization/PCA ->
graph -> clustering with local-LLM annotation -> SCENIC+ regulons -> validation ->
report) on a 3,000-cell subset of the 10x PBMC 10k Multiome data. Unlike
`examples/scenicplus_smoke/`, this is a genome-wide run with real marker-based cell
types, so the results are interpretable, with the caveats below.

## How the run was made

- **Data:** `public_data/pbmc10k_multiome/` (healthy donor, granulocytes removed);
  3,000 of 11,898 cells drawn at random (seed 17), all 36,601 genes, all chromosomes.
  The ATAC fragments were filtered to those cells only. Cells are the only subsetting.
- **Cell types:** Leiden clusters annotated against the bundled human PBMC panel
  (`scrna_workflow/knowledge/human_pbmc_markers.yaml`) reviewed by a local Ollama model
  (`qwen2.5:7b`, `annotation_backend: ollama`). A label is accepted only when at least two
  panel markers support it; 3 of 9 clusters remained `unknown` and form no pseudobulk.
- **Configuration:** `run_scripts/pbmc10k_subset.yaml`; reduced topic modelling
  (3 models, 150 iterations) for run time. Launched with `run_scripts/run_subset.sbatch`
  on a SLURM cluster (40 CPUs, 128 GB); see `run_scripts/slurm-6506341.out`.
- **Run recorded here:** SLURM job 6506341, 2026-09-21, about 2 hours wall clock.
  Code as committed with these files.

## Results in brief

45 direct eRegulons (24 TFs, 13,412 TF-region-gene links) and 25 extended eRegulons;
1,558 of 3,000 cells scored (the rest failed ATAC QC). Networks were built for the five
labelled cell types. The most cell-type-specific eRegulons (RSS) are lineage regulators:
SPIB and POU2AF1 (B cells), TBX21 and RUNX3 (NK cells), SPI1 and KLF4 (FCGR3A+
monocytes), LEF1 and TCF7 (T cells), BCL11A (plasmacytoid dendritic cells).

## Contents

| Path | What it is |
|---|---|
| `regulon/eregulon_triplets_{direct,extended}.tsv` | TF-region-gene links per eRegulon |
| `regulon/eregulon_activity_*_{direct,extended}.tsv.gz` | Cells x eRegulons AUC (gene- and region-based); cells failing ATAC QC are empty |
| `regulon/activity_group_<n>[_<modality>_<kind>].tsv` | Mean AUC per group over scored cells, with `n_cells`/`n_cells_scored` (0 = cell type, 1 = cluster) |
| `regulon/rss_group_<n>[...].tsv`, `rss_ranks*.png`, `rss_heatmap*.png` | Regulon specificity scores, rank plots and heatmaps per cell type |
| `regulon/networks/` | Per-cell-type TF -> region -> gene views: PNG, GraphML (Cytoscape/Gephi), node/edge TSVs, `summary.tsv` |
| `regulon/diagnostics.json` | Stage fingerprints, metrics, warnings, package versions |
| `clustering/annotations.json`, `llm_annotation_audit.json` | Cluster labels, and the model's proposal, accepted label and reason per cluster |
| `clustering/embedding.png`, `resolution_diagnostics.csv` | UMAP by label; resolution choice |
| `report/` | `scientific_report.md`, evidence ledger and figures |
| `stages/<stage>.json` | Per-stage SCENIC+ results (peaks, cisTopic/MALLET, database checks, eGRN) |
| `run_state.json`, `events.jsonl` | Task statuses, metrics, warnings and the task timeline |
| `run_scripts/` | Config, SLURM script and job log |

Not included (large, regenerate by rerunning): `clustering/annotated.h5ad` (159 MB),
`clustering/markers.csv`, `clustering/report.html`, the per-stage SCENIC+ intermediates
(cisTopic object, topic models, Snakemake outputs, bigWigs) and the input subset.

## Reading these results

- eRegulons are inferred associations from one donor's cells: TF-target co-variation,
  region-gene correlation and motif enrichment measured in the same cells. They are not
  TF binding, perturbation or causal evidence.
- Networks are cell-type *views* of one global network, filtered by RSS, that type's
  pseudobulk peaks and detection in its cells; they are not separately inferred networks.
- No CD14+ monocyte cluster passed the two-marker rule, so a major PBMC population is
  probably inside one of the `unknown` clusters and has no pseudobulk, region sets or
  network here.
- 3,000 cells is a quarter of the dataset; rare populations are thin, and AUC and RSS
  values depend on the cell types present.
