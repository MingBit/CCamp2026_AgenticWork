"""Stage 1 (scenicplus env): cell-type pseudobulks, MACS2 peak calling, consensus peaks.

Fragments of cells in eligible cell types are split into one pseudobulk per type
(chromosomes restricted to `keep_chromosomes`), peaks are called per pseudobulk with
MACS2 and merged into fixed-width consensus peaks (Corces et al. 2018, as implemented
in pycisTopic). Consensus peaks outside `keep_chromosomes` (unplaced/unlocalised
contigs, chrM) are removed before any downstream use.
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402


def run(p):
    import pandas as pd
    import pyranges as pr
    from pycisTopic.iterative_peak_calling import get_consensus_peaks
    from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling

    stage = Path(p["stage_dir"])
    sample_id, keep = p["sample_id"], list(p["keep_chromosomes"])
    macs2 = shutil.which("macs2")
    if macs2 is None:
        raise RuntimeError("macs2 is not on PATH in the SCENIC+ environment.")

    cells = pd.read_table(p["cells_tsv"], dtype=str, keep_default_na=False)
    eligible = cells.loc[cells["pseudobulk_group"] != ""]
    cell_data = pd.DataFrame(
        {"barcode": eligible["barcode"].to_numpy(), "sample_id": sample_id,
         "pseudobulk_group": eligible["pseudobulk_group"].to_numpy()},
        index=[common.cell_name(b, sample_id) for b in eligible["barcode"]])

    chromsizes = pd.read_table(p["chromsizes"])
    chromsizes = chromsizes.loc[chromsizes["Chromosome"].isin(keep), ["Chromosome", "Start", "End"]]
    absent = sorted(set(keep) - set(chromsizes["Chromosome"]))
    if absent:
        raise ValueError(f"keep_chromosomes not present in chromsizes: {absent}")

    _, bed_paths = export_pseudobulk(
        input_data=cell_data, variable="pseudobulk_group", sample_id_col="sample_id",
        chromsizes=chromsizes, bed_path=str(stage / "pseudobulk_bed"),
        bigwig_path=str(stage / "pseudobulk_bigwig"), path_to_fragments={sample_id: p["fragments"]},
        n_cpu=int(p["n_cpu"]), normalize_bigwig=True, split_pattern=common.SPLIT_PATTERN,
        temp_dir=p["temp_dir"])
    missing = sorted(set(cell_data["pseudobulk_group"]) - set(bed_paths))
    if missing:
        raise RuntimeError(f"No pseudobulk fragments were written for: {missing}")

    narrow_peaks = peak_calling(
        macs_path=macs2, bed_paths=bed_paths, outdir=str(stage / "macs2"),
        genome_size=p["macs_genome_size"], n_cpu=int(p["n_cpu"]), input_format="BEDPE",
        shift=int(p["macs_shift"]), ext_size=int(p["macs_extsize"]), keep_dup=p["macs_keep_dup"],
        q_value=float(p["macs_qvalue"]), _temp_dir=p["temp_dir"])
    without_peaks = sorted(set(bed_paths) - set(narrow_peaks))

    consensus = get_consensus_peaks(
        narrow_peaks_dict=narrow_peaks, peak_half_width=int(p["peak_half_width"]),
        chromsizes=chromsizes, path_to_blacklist=p["blacklist"])
    peaks = consensus.df
    kept = peaks.loc[peaks["Chromosome"].astype(str).isin(keep)]
    bed = stage / "consensus_regions.bed"
    pr.PyRanges(kept).to_bed(path=str(bed), keep=True, compression="infer", chain=False)

    warnings = []
    if without_peaks:
        warnings.append(f"MACS2 returned no peaks for pseudobulks: {without_peaks}")
    return {
        "outputs": {"consensus_regions": str(bed), "pseudobulk_bed_dir": str(stage / "pseudobulk_bed"),
                    "pseudobulk_bigwig_dir": str(stage / "pseudobulk_bigwig"), "macs2_dir": str(stage / "macs2")},
        "metrics": {"pseudobulks": len(bed_paths), "cells_in_pseudobulks": int(len(cell_data)),
                    "peaks_per_pseudobulk": {k: int(len(v)) for k, v in narrow_peaks.items()},
                    "consensus_peaks_before_chromosome_filter": int(len(peaks)),
                    "consensus_peaks_removed_outside_keep_chromosomes": int(len(peaks) - len(kept)),
                    "consensus_peaks": int(len(kept))},
        "warnings": warnings,
        "versions": common.package_versions(["pycisTopic", "MACS2", "pyranges"]),
    }


if __name__ == "__main__":
    common.stage_main(run)
