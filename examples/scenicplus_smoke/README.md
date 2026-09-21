# SCENIC+ smoke-run outputs (example only)

Outputs of the SCENIC+ regulon step (`regulon_method: scenicplus`) run on a deliberately
small subset of the 10x PBMC 10k Multiome data. They show what the step produces and let
you check formats and downstream code without a server run. **They are not biological
results**: the subset is restricted to two chromosomes, cell labels are crude single-marker
calls, and several parameters are relaxed.

## How the run was made

- **Data:** 1,200 random cells from `public_data/pbmc10k_multiome/` (seed 1); genes on
  chr21/chr22 plus all TFs in `allTFs_hg38.txt` (2,474 genes); ATAC fragments of those cells on
  chr21, chr22 and unplaced GL/KI contigs.
- **Cell labels (`annotation`):** `T cell` if CD3E > 0, else `B cell` if MS4A1 > 0, else
  `Myeloid` if LYZ > 0, else `unknown` (unknown cells are kept but form no pseudobulk).
- **Settings:** `configs/pbmc10k_multiome.yaml` with the overrides in
  `run_scripts/smoke_regulon.py`: `keep_chromosomes` chr21/chr22, topic models with 5 and 10
  topics and 50 iterations, relaxed ATAC QC minima without Otsu thresholds, `gsea_n_perm` 100,
  `min_target_genes` 3, MALLET memory 8 GB.
- **Execution:** SLURM, 8 CPUs; the regulon step called directly (no inspection/QC/clustering
  tasks), see `run_scripts/`. `run_scripts/smoke_regulon.py` is kept as it ran; with the current
  package layout import `_result` from `scrna_workflow.tools._downstream_common` and `PATH_KEYS`
  from `scrna_workflow.cli` instead of `scrna_workflow.downstream` / `scrna_workflow.runner`. Steps to reproduce are in the SCENIC+ section of the main README.
- **Run recorded here:** SLURM job 6504131 on eodd, 2026-09-17, stages 1-3 from job 6504017
  (2026-09-16). Code as committed with these files, except that stage 4 ran before the
  direct-only fallback (`scenicplus_allow_direct_only`) was added; both eRegulon sets were built,
  so the fallback did not apply. Update this line when copying outputs of a new run.

## Contents

| Path | What it is |
|---|---|
| `regulon/eregulon_triplets_direct.tsv` | TF-region-gene links of eRegulons whose TF-motif link has direct evidence |
| `regulon/eregulon_triplets_extended.tsv` | Same for motifs linked to the TF through orthology (weaker evidence); absent if none were built |
| `regulon/eregulon_activity_gene_based_direct.tsv.gz` | Cells x eRegulons, AUC of each eRegulon's target genes; cells failing ATAC QC are empty |
| `regulon/eregulon_activity_region_based_direct.tsv.gz` | Same, scored on target regions |
| `regulon/eregulon_activity_*_extended.tsv.gz` | Extended-eRegulon AUC. In this example taken from stage 4, so only scored cells are listed; newer runs write it aligned to all cells like the direct tables |
| `regulon/activity_group_<n>.tsv` | Mean direct gene-based AUC per group over scored cells; `<n>` follows `summaries` in `diagnostics.json` (here 0 = `annotation`, 1 = `cluster`, identical labels) |
| `regulon/activity_group_<n>_<modality>_<kind>.tsv` | Same for region-based and extended AUC. Newer runs only; every summary also starts with `n_cells` and `n_cells_scored` columns. This example predates both |
| `regulon/rss_group_<n>[_<modality>_<kind>].tsv`, `rss_ranks*.png`, `rss_heatmap*.png` | Regulon specificity scores per group and, for the cell-type column, rank plots and heatmaps. Newer runs only; this example predates them |
| `regulon/networks/<cell type>.png`, `.graphml`, `_nodes.tsv`, `_edges.tsv`, `summary.tsv` | Cell-type views of the direct eRegulon network (TF -> region -> gene). Newer runs only; this example predates them |
| `regulon/diagnostics.json` | Hand-off numbers and every stage's fingerprint, metrics, warnings, versions and output paths |
| `stages/<stage>.json` | Each stage's result record; `stage2_otsu_thresholds.tsv` holds pycisTopic's automatic QC thresholds |
| `scenicplus_reports/ctx_results.html`, `dem_results.html` | SCENIC+ motif enrichment reports (cisTarget and DEM) |
| `scenicplus_reports/snakemake_config.yaml` | The SCENIC+ Snakemake config generated from our settings |
| `run_scripts/` | Driver script, SLURM script and SLURM log of the run |

Column meanings of the triplet tables and interpretation caveats are described in the SCENIC+
section and "Scientific boundaries" of the main README. Paths inside the JSON files point to the
server run directory and are kept as a record only.

Not included (large, regenerate with a run): `cistopic_object.pkl`, `topic_models.pkl`,
`gex.h5ad`, `ACC_GEX.h5mu`, `AUCell_*.h5mu`, `scplusmdata.h5mu`, pseudobulk BED/bigWig files,
`tf_to_gene_adj.tsv`, `region_to_gene_adj.tsv`, `search_space.tsv`, consensus peaks and region sets.

## Updating these files

From a machine that can see the smoke directory:

```bash
examples/scenicplus_smoke/collect_outputs.sh /path/to/Smoke_Data <slurm job id>
```

The script only copies the files listed above (about 1-2 MB) and overwrites previous copies.
