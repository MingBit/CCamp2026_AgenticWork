"""Stage 2 (scenicplus env): ATAC QC, cisTopic object, MALLET topic models, region sets.

Cells must pass pycisTopic ATAC QC on the consensus peaks and be retained by the
workflow's RNA QC. Topic models are trained with MALLET; the selected model yields
region sets as binarized topics (Otsu and top-n) and cell-type differentially
accessible regions (DARs), written as BED files for SCENIC+ motif enrichment.
"""
import os
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

EXCLUDED_GROUP = "__not_in_pseudobulk__"


def qc_thresholds(otsu_file, minima, automatic):
    """Per-metric threshold: user minimum, raised to the Otsu threshold when automatic."""
    header, values = Path(otsu_file).read_text().splitlines()[:2]
    otsu = dict(zip(header.split("\t"), (float(v) for v in values.split("\t"))))
    unique = float(minima["unique_fragments_in_peaks"])
    tss = float(minima["tss_enrichment"])
    if automatic:
        unique = max(unique, otsu["unique_fragments_in_peaks_count_otsu_threshold"])
        tss = max(tss, otsu["tss_enrichment_otsu_threshold"])
    return {"unique_fragments_threshold": unique, "tss_enrichment_threshold": tss,
            "frip_threshold": float(minima["frip"])}, otsu


def run(p):
    import numpy as np
    import pandas as pd
    from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
    from pycisTopic.cli.subcommand.qc import qc
    from pycisTopic.diff_features import (find_diff_features, find_highly_variable_features,
                                          impute_accessibility, normalize_scores)
    from pycisTopic.lda_models import evaluate_models, run_cgs_models_mallet
    from pycisTopic.qc import get_barcodes_passing_qc_for_sample
    from pycisTopic.topic_binarization import binarize_topics

    stage = Path(p["stage_dir"])
    sample_id, n_cpu = p["sample_id"], int(p["n_cpu"])
    cells = pd.read_table(p["cells_tsv"], dtype=str, keep_default_na=False)

    # ATAC QC on consensus peaks, thresholds bounded below by configured minima.
    qc_dir = stage / "qc"
    qc_dir.mkdir(exist_ok=True)
    qc(fragments_tsv_filename=p["fragments"], regions_bed_filename=p["consensus_regions"],
       tss_annotation_bed_filename=p["tss_annotation"], output_prefix=str(qc_dir / sample_id), no_threads=n_cpu)
    thresholds, otsu = qc_thresholds(qc_dir / f"{sample_id}.otsu_thresholds.tsv", p["qc_minima"],
                                     bool(p["qc_use_automatic_thresholds"]))
    passing, _ = get_barcodes_passing_qc_for_sample(
        sample_id=sample_id, pycistopic_qc_output_dir=str(qc_dir), use_automatic_thresholds=False, **thresholds)
    passing = set(np.atleast_1d(passing).astype(str))
    rna = set(cells["barcode"])
    valid = sorted(passing & rna)
    if len(valid) < int(p["min_cells"]):
        raise RuntimeError(f"Only {len(valid)} cells pass both ATAC QC and RNA QC (minimum {p['min_cells']}).")
    pd.Series(valid, name="barcode").to_csv(stage / "cells_passing_atac_and_rna_qc.tsv", sep="\t", index=False)

    # cisTopic object: consensus peaks x cells passing both modalities' QC.
    cistopic = create_cistopic_object_from_fragments(
        path_to_fragments=p["fragments"], path_to_regions=p["consensus_regions"],
        path_to_blacklist=p["blacklist"], valid_bc=valid, n_cpu=n_cpu, project=sample_id,
        split_pattern=common.SPLIT_PATTERN)
    meta = cells.set_index("barcode").loc[valid, ["cell_type", "pseudobulk_group"]].copy()
    meta["dar_group"] = meta["pseudobulk_group"].replace("", EXCLUDED_GROUP)
    meta.index = [common.cell_name(b, sample_id) for b in meta.index]
    cistopic.add_cell_data(meta, split_pattern=common.SPLIT_PATTERN)

    # Topic modelling with MALLET; one model per topic number.
    os.environ["MALLET_MEMORY"] = f"{int(p['mallet_memory_gb'])}G"
    models = run_cgs_models_mallet(
        cistopic, n_topics=[int(n) for n in p["n_topics"]], n_cpu=n_cpu, n_iter=int(p["topic_n_iter"]),
        random_state=int(p["seed"]), alpha=float(p["topic_alpha"]), alpha_by_topic=bool(p["topic_alpha_by_topic"]),
        eta=float(p["topic_eta"]), eta_by_topic=bool(p["topic_eta_by_topic"]), tmp_path=p["temp_dir"],
        save_path=None, mallet_path=p["mallet_path"])
    with open(stage / "topic_models.pkl", "wb") as handle:
        pickle.dump(models, handle)
    selection = p["topic_selection"]
    model = evaluate_models(models, select_model=None if selection == "auto" else int(selection),
                            return_model=True, plot=False)
    cistopic.add_LDA_model(model)

    # Region sets: binarized topics and cell-type DARs.
    region_sets = stage / "region_sets"
    counts = {}
    ntop = int(p["topic_binarization_ntop"])
    binarizations = {f"Topics_top_{ntop}": binarize_topics(cistopic, method="ntop", ntop=ntop, plot=False)}
    if p["topic_binarization_otsu"]:
        binarizations["Topics_otsu"] = binarize_topics(cistopic, method="otsu", plot=False)
    for folder, topics in binarizations.items():
        for topic, table in topics.items():
            counts[f"{folder}/{topic}"] = common.write_region_bed(region_sets / folder / f"{topic}.bed", table.index)

    imputed = impute_accessibility(cistopic, selected_cells=None, selected_regions=None, scale_factor=10**6)
    normalized = normalize_scores(imputed, scale_factor=10**4)
    variable_regions = find_highly_variable_features(
        normalized, min_disp=0.05, min_mean=0.0125, max_mean=3, max_disp=np.inf, n_bins=20,
        n_top_features=None, plot=False)
    dars = find_diff_features(
        cistopic, imputed, variable="dar_group", var_features=variable_regions, contrasts=None,
        adjpval_thr=float(p["dar_adjpval_threshold"]), log2fc_thr=float(p["dar_log2fc_threshold"]),
        split_pattern=common.SPLIT_PATTERN, n_cpu=n_cpu, _temp_dir=p["temp_dir"])
    for group, table in dars.items():
        if group != EXCLUDED_GROUP and len(table):
            counts[f"DARs_cell_type/{group}"] = common.write_region_bed(
                region_sets / "DARs_cell_type" / f"{common.safe_label(group)}.bed", table.index)
    if not counts:
        raise RuntimeError("No region sets were produced from topics or DARs.")

    cistopic_file = stage / "cistopic_object.pkl"
    with open(cistopic_file, "wb") as handle:
        pickle.dump(cistopic, handle)
    warnings = []
    empty_dars = sorted(g for g, t in dars.items() if g != EXCLUDED_GROUP and not len(t))
    if empty_dars:
        warnings.append(f"No DARs passed thresholds for: {empty_dars}")
    return {
        "outputs": {"cistopic_object": str(cistopic_file), "region_set_folder": str(region_sets),
                    "topic_models": str(stage / "topic_models.pkl"), "qc_dir": str(qc_dir),
                    "cells": str(stage / "cells_passing_atac_and_rna_qc.tsv")},
        "metrics": {"barcodes_passing_atac_qc": len(passing), "rna_qc_cells": len(rna), "cells_used": len(valid),
                    "qc_thresholds": thresholds, "otsu_thresholds": otsu,
                    "regions": int(len(cistopic.region_names)), "selected_topics": int(model.n_topic),
                    "variable_regions": int(len(variable_regions)), "region_set_sizes": counts},
        "warnings": warnings,
        "versions": common.package_versions(["pycisTopic", "gensim", "lda"]),
    }


if __name__ == "__main__":
    common.stage_main(run)
