"""Build offline gene annotations for pycisTopic and SCENIC+ from the Ensembl 98 GTF.

Ensembl 98 (GENCODE v32) is the annotation of the 10x GRCh38-2020-A reference used by
Cell Ranger ARC for the PBMC 10k Multiome data; its BioMart archive is retired, so the
BioMart queries both tools would run are reproduced here from the release GTF.

Outputs (UCSC chromosome names, assembled chromosomes only):
  tss_ensembl98_protein_coding.bed  pycisTopic TSS BED, as `pycistopic tss get_tss`
  genome_annotation.tsv             SCENIC+ `download_genome_annotations` gene annotation
  chromsizes.tsv                    SCENIC+ `download_genome_annotations` chromosome sizes

Usage: python3 make_annotations.py  (stdlib only; run from this directory)
"""
import csv
import gzip
import re

GTF = "Homo_sapiens.GRCh38.98.gtf.gz"
UCSC_SIZES = "hg38.chrom.sizes"
ASSEMBLED = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
TO_UCSC = {c: "chrM" if c == "MT" else f"chr{c}" for c in ASSEMBLED}
ATTR = re.compile(r'(\S+) "([^"]*)"')

genes, transcripts = {}, []
with gzip.open(GTF, "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        chrom, _, feature, start, end, _, strand, _, attributes = line.rstrip("\n").split("\t")
        if chrom not in TO_UCSC or feature not in ("gene", "transcript"):
            continue
        attrs = dict(ATTR.findall(attributes))
        if feature == "gene":
            genes[attrs["gene_id"]] = (int(start), int(end))
        elif attrs.get("transcript_biotype") == "protein_coding":
            tss = int(start) if strand == "+" else int(end)
            name = attrs.get("gene_name", attrs["gene_id"])
            transcripts.append((TO_UCSC[chrom], tss, strand, name, attrs["gene_id"]))

# pycisTopic: chromosome, TSS (0-based BED), gene name, score, strand, transcript type.
tss_rows = sorted({(c, tss - 1, tss, name, ".", strand, "protein_coding") for c, tss, strand, name, _ in transcripts})
with open("tss_ensembl98_protein_coding.bed", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["# Chromosome", "Start", "End", "Gene", "Score", "Strand", "Transcript_type"])
    writer.writerows(tss_rows)

# SCENIC+: gene start/end (1-based, BioMart style) with one row per distinct transcript TSS.
annotation_rows = sorted({(c, *genes[gene_id], strand, name, tss, "protein_coding") for c, tss, strand, name, gene_id in transcripts})
with open("genome_annotation.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["Chromosome", "Start", "End", "Strand", "Gene", "Transcription_Start_Site", "Transcript_type"])
    writer.writerows(annotation_rows)

sizes = dict(line.split("\t") for line in open(UCSC_SIZES).read().splitlines())
with open("chromsizes.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["Chromosome", "Start", "End"])
    writer.writerows((TO_UCSC[c], 0, int(sizes[TO_UCSC[c]])) for c in ASSEMBLED)

print(f"TSS rows: {len(tss_rows)}; annotation rows: {len(annotation_rows)}; "
      f"genes: {len({row[4] for row in annotation_rows})}; chromosomes: {len(ASSEMBLED)}")
