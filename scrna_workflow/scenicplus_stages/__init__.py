"""SCENIC+ stage scripts, executed in the separate `scenicplus` conda environment.

stage1_peaks            cell-type pseudobulks, MACS2 peaks, consensus peaks, contig removal
stage2_cistopic         ATAC QC, cisTopic object, MALLET topics, topic/DAR region sets
stage3_motif_databases  cisTarget database integrity and consensus-peak coverage
stage4_scenicplus       SCENIC+ Snakemake eGRN inference and TSV export

Only `common` is imported by the workflow environment; stage modules import
pycisTopic/SCENIC+ inside `run()`.
"""
