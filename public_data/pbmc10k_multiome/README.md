# PBMC 10k Multiome (scRNA-seq + scATAC-seq)

10x Genomics public dataset "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)",
Chromium Single Cell Multiome ATAC + Gene Expression, processed with Cell Ranger ARC 2.0.0 (GRCh38).
Source: https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0

Data files are git-ignored. To re-download:

```bash
B=https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k
for s in summary.csv web_summary.html atac_peaks.bed atac_peak_annotation.tsv atac_fragments.tsv.gz.tbi \
         per_barcode_metrics.csv filtered_feature_bc_matrix.h5 atac_fragments.tsv.gz; do
  curl -fL -C - -o pbmc_granulocyte_sorted_10k_$s "${B}_$s"
done
```

| File | Contents |
|---|---|
| `*_filtered_feature_bc_matrix.h5` | Cell x feature counts; both "Gene Expression" and "Peaks" feature types |
| `*_atac_fragments.tsv.gz` (+ `.tbi`) | Per-cell ATAC fragments (for Signac/ArchR/SnapATAC2) |
| `*_atac_peaks.bed` | Called ATAC peaks |
| `*_atac_peak_annotation.tsv` | Peak-to-gene annotation |
| `*_per_barcode_metrics.csv` | Per-barcode QC for both modalities |
| `*_summary.csv`, `*_web_summary.html` | Run-level QC summary |
