"""Stage 3 (scenicplus env): validate precomputed cisTarget motif databases.

The hg38 SCREEN databases score predefined regions, not our consensus peaks. SCENIC+
maps region sets onto database regions by overlap, so this stage verifies database
integrity (upstream SHA-1), consistency between the rankings/scores databases and the
motif annotation, and reports how much of the consensus peak set the databases cover.
"""
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

MIN_COVERED_FRACTION = 0.5


def sha1(path, block=16 * 1024 * 1024):
    digest = hashlib.sha1()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(block), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_checksum(path):
    """Compare against `<db>.sha1sum.txt` when present; returns (checked, ok)."""
    sums = Path(str(path) + ".sha1sum.txt")
    if not sums.is_file():
        return False, None
    expected = sums.read_text().split()[0]
    return True, sha1(path) == expected


def database_names(path):
    """Region names (columns) and motif names (index column) of a cisTarget feather DB.

    The databases are compressed Arrow IPC files stored as one record batch; only the
    index column is decoded, reading a whole batch needs tens of GB of memory.
    """
    import pyarrow as pa
    import pyarrow.ipc as ipc
    with pa.memory_map(str(path), "r") as source:
        names = ipc.open_file(source).schema.names
        index_column = next((n for n in ("motifs", "tracks", "features") if n in names), None)
        motifs = []
        if index_column is not None:
            options = ipc.IpcReadOptions(included_fields=[names.index(index_column)])
            reader = ipc.open_file(source, options=options)
            for i in range(reader.num_record_batches):
                motifs.extend(reader.get_batch(i).column(0).to_pylist())
    return [n for n in names if n != index_column], motifs


def run(p):
    import pandas as pd
    import pyranges as pr

    warnings, checksums = [], {}
    for key in ("ctx_db", "dem_db"):
        checked, ok = verify_checksum(p[key])
        checksums[key] = "verified" if ok else ("not_available" if not checked else "MISMATCH")
        if checked and not ok:
            raise RuntimeError(f"SHA-1 mismatch for {p[key]}; re-download the database.")
        if not checked:
            warnings.append(f"No upstream checksum found for {p[key]}; integrity not verified.")

    ctx_regions, ctx_motifs = database_names(p["ctx_db"])
    dem_regions, dem_motifs = database_names(p["dem_db"])
    if set(ctx_regions) != set(dem_regions):
        raise RuntimeError("Rankings and scores databases contain different regions.")

    annotation = pd.read_table(p["motif_annotations"], dtype=str)
    motif_column = next(c for c in ("#motif_id", "motif_id") if c in annotation.columns)
    annotated = set(annotation[motif_column].astype(str))
    motif_coverage = len(annotated & set(ctx_motifs)) / max(len(set(ctx_motifs)), 1)
    if motif_coverage == 0:
        raise RuntimeError("No database motif has an entry in the motif annotation table.")

    database = pr.PyRanges(pd.DataFrame([common.parse_region(r) for r in ctx_regions if ":" in r],
                                        columns=["Chromosome", "Start", "End"]))
    consensus = pr.read_bed(p["consensus_regions"]).df[["Chromosome", "Start", "End"]]
    consensus = pr.PyRanges(consensus.assign(Name=range(len(consensus))))
    joined = consensus.join(database, report_overlap=True).df
    minimum = float(p["fraction_overlap_w_ctx_database"])
    if len(joined):
        # Same rule as pycistarget.utils.target_to_query: overlap fraction of either region.
        fraction = (joined["Overlap"] / (joined["End"] - joined["Start"])).combine(
            joined["Overlap"] / (joined["End_b"] - joined["Start_b"]), max)
        covered = joined.loc[fraction > minimum, "Name"].nunique()
    else:
        covered = 0
    covered_fraction = covered / max(len(consensus), 1)
    if covered_fraction < MIN_COVERED_FRACTION:
        warnings.append(f"Only {covered_fraction:.1%} of consensus peaks overlap database regions by > {minimum:.0%} of either region; "
                        "consider a custom cisTarget database built from the consensus peaks.")
    return {
        "outputs": {"ctx_db": p["ctx_db"], "dem_db": p["dem_db"], "motif_annotations": p["motif_annotations"]},
        "metrics": {"checksums": checksums, "database_regions": len(ctx_regions), "database_motifs": len(ctx_motifs),
                    "annotated_motif_fraction": round(motif_coverage, 4), "consensus_peaks": int(len(consensus)),
                    "consensus_peaks_covered_by_database": int(covered),
                    "consensus_peak_coverage_fraction": round(covered_fraction, 4),
                    "dem_database_motifs": len(dem_motifs)},
        "warnings": warnings,
        "versions": common.package_versions(["pyarrow", "pyranges", "pycistarget"]),
    }


if __name__ == "__main__":
    common.stage_main(run)
