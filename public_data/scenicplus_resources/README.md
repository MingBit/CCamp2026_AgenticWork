# SCENIC+ reference resources (hg38)

Local copies of everything pycisTopic/SCENIC+ would otherwise fetch at run time, so the
workflow runs offline. Data files are git-ignored (~49 GB); rebuild with the commands below
from this directory.

## `hg38/cistarget/` (aertslab cisTarget resources)

| File | Use |
|---|---|
| `hg38_screen_v10_clust.regions_vs_motifs.rankings.feather` (35 GB) | cisTarget motif enrichment (`ctx_db_fname`) |
| `hg38_screen_v10_clust.regions_vs_motifs.scores.feather` (14 GB) | DEM motif enrichment (`dem_db_fname`) |
| `*.feather.sha1sum.txt` | Upstream SHA-1 checksums |
| `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` | Motif-to-TF annotation (`path_to_motif_annotations`) |
| `allTFs_hg38.txt` | Human TF list |

```bash
B=https://resources.aertslab.org/cistarget
DB=$B/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs
mkdir -p hg38/cistarget && cd hg38/cistarget
for s in rankings scores; do
  curl -fL -C - -O $DB.$s.feather && curl -fLO $DB.$s.feather.sha1sum.txt
done
sha1sum -c *.sha1sum.txt
curl -fLO $B/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
curl -fLO $B/tf_lists/allTFs_hg38.txt
```

The SCREEN databases score predefined ENCODE SCREEN regions. Regions in the consensus
peaks are matched to them by overlap (`fraction_overlap_w_*_database`). A custom database
built from the consensus peaks (`envs/create_cistarget_databases.yml`) is the alternative.

## `hg38/genome/`

| File | Use |
|---|---|
| `hg38-blacklist.v2.bed` | ENCODE blacklist v2 (copy shipped in pycisTopic) for peak calling/QC |
| `hg38.chrom.sizes` | UCSC chromosome sizes |
| `Homo_sapiens.GRCh38.98.gtf.gz` | Ensembl 98 GTF (= GENCODE v32, annotation of 10x GRCh38-2020-A) |
| `tss_ensembl98_protein_coding.bed` | pycisTopic TSS annotation (QC, TSS enrichment) |
| `genome_annotation.tsv`, `chromsizes.tsv` | SCENIC+ search space (replace its BioMart/NCBI download step) |

```bash
mkdir -p hg38/genome && cd hg38/genome
curl -fLO https://raw.githubusercontent.com/aertslab/pycisTopic/787ce422a37f5975b0ebb9e7b19eeaed44847501/blacklist/hg38-blacklist.v2.bed
curl -fLO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
curl -fLO https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
curl -fL -o ensembl98_CHECKSUMS https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/CHECKSUMS
sum -r Homo_sapiens.GRCh38.98.gtf.gz; grep GRCh38.98.gtf.gz ensembl98_CHECKSUMS   # must match
python3 make_annotations.py
```

`make_annotations.py` reproduces the protein-coding BioMart queries of
`pycistopic tss get_tss` and `scenicplus prepare_data download_genome_annotations`,
keeping chromosomes 1–22, X, Y and M with UCSC names. The Ensembl 98 BioMart archive
is retired, so the release GTF is used instead.
