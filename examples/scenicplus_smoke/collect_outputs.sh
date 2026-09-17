#!/usr/bin/env bash
# Copy the small, shareable outputs of a SCENIC+ smoke run into this folder.
# Usage: examples/scenicplus_smoke/collect_outputs.sh <smoke dir> [slurm job id]
#   <smoke dir>  the $SMOKE directory used for the run (contains run/, work/, smoke_regulon.py)
# Large intermediates (pickles, h5ad/h5mu, bigWigs, adjacency tables) are deliberately not copied.
set -euo pipefail

usage="usage: collect_outputs.sh <smoke dir> [slurm job id]  (outputs are written next to this script)"
if [ $# -lt 1 ] || [ $# -gt 2 ]; then echo "$usage" >&2; exit 2; fi
SRC=$1
JOB=${2:-}
if [ ! -d "$SRC/run/regulon" ] || [ ! -d "$SRC/work" ]; then
  echo "Not a smoke run directory (expected run/regulon/ and work/ inside): $SRC" >&2
  echo "$usage" >&2; exit 2
fi
if [ -n "$JOB" ] && ! [[ "$JOB" =~ ^[0-9]+$ ]]; then echo "SLURM job id must be a number, got: $JOB" >&2; exit 2; fi
DEST=$(cd "$(dirname "$0")" && pwd)
WORK="$SRC/work"
mkdir -p "$DEST/regulon" "$DEST/stages" "$DEST/scenicplus_reports" "$DEST/run_scripts"

copy_if_present() {  # copy_if_present <source> <destination>
  if [ -f "$1" ]; then cp "$1" "$2"; else echo "skipped (not found): $1"; fi
}

# Final regulon outputs, as produced by the workflow's regulon task.
cp "$SRC"/run/regulon/* "$DEST/regulon/"
# Extended activity scores exist only when extended eRegulons were built.
for f in auc_gene_based_extended.tsv.gz auc_region_based_extended.tsv.gz; do
  copy_if_present "$WORK/stage4_scenicplus/$f" "$DEST/regulon/eregulon_activity_${f#auc_}"
done

# Per-stage metrics, warnings, package versions and fingerprints.
for stage in stage1_peaks stage2_cistopic stage3_motif_databases stage4_scenicplus; do
  copy_if_present "$WORK/$stage/stage_result.json" "$DEST/stages/$stage.json"
done
copy_if_present "$WORK/stage2_cistopic/qc/pbmc10k.otsu_thresholds.tsv" "$DEST/stages/stage2_otsu_thresholds.tsv"

# SCENIC+ motif enrichment reports (open in a browser).
for f in ctx_results.html dem_results.html; do
  copy_if_present "$WORK/stage4_scenicplus/Snakemake/$f" "$DEST/scenicplus_reports/$f"
done
copy_if_present "$WORK/stage4_scenicplus/Snakemake/config/config.yaml" "$DEST/scenicplus_reports/snakemake_config.yaml"

# How the run was launched.
copy_if_present "$SRC/smoke_regulon.py" "$DEST/run_scripts/smoke_regulon.py"
copy_if_present "$SRC/smoke_regulon.sbatch" "$DEST/run_scripts/smoke_regulon.sbatch"
if [ -n "$JOB" ]; then copy_if_present "$SRC/slurm-$JOB.out" "$DEST/run_scripts/slurm-$JOB.out"; fi

du -sh "$DEST"
