"""SCENIC+ regulon specialist: orchestrates the stage scripts from the workflow environment.

SCENIC+ pins pandas 1.5/scanpy 1.8, so its code runs as subprocesses in its own conda
environment (envs/scenicplus.yml). This module checks prerequisites, hands the RNA
counts and cell-type labels over in version-neutral files, runs the four stages with
per-stage fingerprints (so a failed or repeated run reuses completed stages from
`scenicplus_work_dir`) and converts SCENIC+ results into the regulon task's outputs.

eRegulons are inferred TF-region-gene associations supported by motif enrichment and
correlation in the same cells; they are not binding or causal evidence.
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

from .scenicplus_stages import common

STAGE_DIR = Path(__file__).resolve().parent / "scenicplus_stages"
RESOURCE_KEYS = {
    "fragments": "scenicplus_fragments_path", "ctx_db": "scenicplus_ctx_db_path",
    "dem_db": "scenicplus_dem_db_path", "motif_annotations": "scenicplus_motif_annotations_path",
    "blacklist": "scenicplus_blacklist_path", "tss_annotation": "scenicplus_tss_annotation_path",
    "genome_annotation": "scenicplus_genome_annotation_path", "chromsizes": "scenicplus_chromsizes_path",
    "mallet_path": "scenicplus_mallet_path",
}
DEFAULT_EXCLUDED_LABELS = ["unknown", "ambiguous"]
SNAKEMAKE_KEYS = [
    "is_multiome", "nr_cells_per_metacells", "search_space_upstream", "search_space_downstream",
    "search_space_extend_tss", "motif_similarity_fdr", "orthologous_identity_threshold", "annotations_to_use",
    "fraction_overlap_w_dem_database", "dem_max_bg_regions", "dem_balance_number_of_promoters",
    "dem_promoter_space", "dem_adj_pval_thr", "dem_log2fc_thr", "dem_mean_fg_thr", "dem_motif_hit_thr",
    "fraction_overlap_w_ctx_database", "ctx_auc_threshold", "ctx_nes_threshold", "ctx_rank_threshold",
    "tf_to_gene_importance_method", "region_to_gene_importance_method", "region_to_gene_correlation_method",
    "order_regions_to_genes_by", "order_tfs_to_genes_by", "gsea_n_perm", "quantile_thresholds_region_to_gene",
    "top_n_region_to_genes_per_gene", "top_n_region_to_genes_per_region", "min_regions_per_gene",
    "rho_threshold", "min_target_genes",
]
STAGE1_KEYS = ["keep_chromosomes", "macs_genome_size", "macs_shift", "macs_extsize", "macs_qvalue",
               "macs_keep_dup", "peak_half_width"]
STAGE2_KEYS = ["qc_use_automatic_thresholds", "n_topics", "topic_n_iter", "topic_alpha", "topic_alpha_by_topic",
               "topic_eta", "topic_eta_by_topic", "topic_selection", "topic_binarization_otsu",
               "topic_binarization_ntop", "dar_adjpval_threshold", "dar_log2fc_threshold", "mallet_memory_gb"]
STAGE_RESOURCES = {
    "stage1_peaks": ["fragments", "chromsizes", "blacklist"],
    "stage2_cistopic": ["fragments", "tss_annotation", "blacklist", "mallet_path"],
    "stage3_motif_databases": ["ctx_db", "dem_db", "motif_annotations"],
    "stage4_scenicplus": ["ctx_db", "dem_db", "motif_annotations", "genome_annotation", "chromsizes"],
}
# Execution settings that do not change results; excluded from stage fingerprints.
EXECUTION_KEYS = {"stage_dir", "n_cpu", "temp_dir", "mallet_memory_gb", "mallet_path"}
INTERPRETATION = ("eRegulons are inferred TF-region-gene associations (TF and target co-variation, region "
                  "accessibility-gene correlation and motif enrichment in the same cells); they are not binding, "
                  "perturbation or causal evidence.")


def _setting(cfg, key):
    value = cfg.get(f"scenicplus_{key}")
    if value is None:
        raise KeyError(f"scenicplus_{key}")
    return value


def prerequisites(cfg, source):
    """Return a list of unmet prerequisites (empty when the SCENIC+ step can run)."""
    problems = []
    if not source:
        problems.append("A validated upstream AnnData dataset is required.")
    for key in RESOURCE_KEYS.values():
        path = cfg.get(key)
        if not path or not Path(path).is_file():
            problems.append(f"Local resource missing: {key} ({path}).")
    mallet = cfg.get("scenicplus_mallet_path")
    if mallet and Path(mallet).is_file() and not os.access(mallet, os.X_OK):
        problems.append(f"MALLET is not executable: {mallet}.")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", str(cfg.get("scenicplus_sample_id") or "")):
        problems.append("scenicplus_sample_id must be set and match [A-Za-z0-9_.-]+.")
    if not cfg.get("scenicplus_cell_type_column"):
        problems.append("scenicplus_cell_type_column must name the obs column with cell-type labels.")
    missing = [f"scenicplus_{k}" for k in STAGE1_KEYS + STAGE2_KEYS + SNAKEMAKE_KEYS if cfg.get(f"scenicplus_{k}") is None]
    missing += [k for k in ("scenicplus_species", "scenicplus_motif_annotation_version", "scenicplus_n_cpu")
                if cfg.get(k) is None]
    if missing:
        problems.append("Missing SCENIC+ parameters: " + ", ".join(missing))
    python = cfg.get("scenicplus_python")
    if python and not Path(python).is_file():
        problems.append(f"scenicplus_python not found: {python}.")
    if not python and shutil.which("conda") is None:
        problems.append("Neither scenicplus_python nor a `conda` executable is available to run the SCENIC+ environment.")
    return problems


def _labels(a, cfg):
    """Cell-type labels and pseudobulk groups (types with enough cells, excluding placeholders)."""
    import pandas as pd
    column = cfg["scenicplus_cell_type_column"]
    if column not in a.obs:
        raise KeyError(column)
    labels = a.obs[column].astype("string").fillna("").astype(str)
    excluded = {str(x) for x in cfg.get("scenicplus_excluded_cell_types", DEFAULT_EXCLUDED_LABELS)} | {""}
    sizes = labels[~labels.isin(excluded)].value_counts()
    eligible = sizes[sizes >= int(cfg.get("scenicplus_min_cells_per_cell_type", 10))].index
    mapping = common.unique_safe_labels(eligible)
    groups = labels.map(lambda x: mapping.get(x, "")).astype(str)
    return labels, groups, pd.Series(sizes, dtype=int)


def export_handoff(a, cfg, directory):
    """Write counts (npz), cells (barcode, labels, groups) and genes; return files and a content fingerprint."""
    import numpy as np
    import pandas as pd
    from scipy import sparse

    directory.mkdir(parents=True, exist_ok=True)
    labels, groups, sizes = _labels(a, cfg)
    barcode_column = cfg.get("scenicplus_barcode_column")
    barcodes = a.obs[barcode_column].astype(str) if barcode_column else pd.Series(a.obs_names.astype(str), index=a.obs_names)
    if not barcodes.is_unique:
        raise ValueError("Cell barcodes used for ATAC matching are not unique.")
    counts = sparse.csr_matrix(a.layers["counts"])
    counts.sort_indices()
    h = hashlib.sha256()
    for part in (counts.data, counts.indices, counts.indptr, np.asarray(counts.shape)):
        h.update(np.ascontiguousarray(part).tobytes())
    cells = pd.DataFrame({"barcode": barcodes.to_numpy(), "cell_type": labels.to_numpy(), "pseudobulk_group": groups.to_numpy()})
    genes = pd.DataFrame({"gene": a.var_names.astype(str)})
    for frame in (cells, genes):
        h.update(frame.to_csv(sep="\t", index=False).encode())
    fingerprint = h.hexdigest()
    files = {"counts_npz": directory / "counts.npz", "cells_tsv": directory / "cells.tsv", "genes_tsv": directory / "genes.tsv"}
    marker = directory / "handoff.json"
    if not (marker.is_file() and common.read_json(marker).get("fingerprint") == fingerprint and all(f.is_file() for f in files.values())):
        sparse.save_npz(files["counts_npz"], counts)
        cells.to_csv(files["cells_tsv"], sep="\t", index=False)
        genes.to_csv(files["genes_tsv"], sep="\t", index=False)
        common.write_json(marker, {"fingerprint": fingerprint})
    groups_used = sorted(set(groups) - {""})
    return {k: str(v) for k, v in files.items()}, fingerprint, {
        "cells": int(a.n_obs), "genes": int(a.n_vars), "pseudobulk_groups": len(groups_used),
        "cells_in_pseudobulk_groups": int((groups != "").sum()), "cell_type_sizes": {str(k): int(v) for k, v in sizes.items()}}


def launcher(cfg):
    """Command prefix and environment for running Python in the SCENIC+ environment."""
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env.update(PYTHONNOUSERSITE="1", MPLBACKEND="Agg")
    python = cfg.get("scenicplus_python")
    if python:
        env["PATH"] = str(Path(python).resolve().parent) + os.pathsep + env.get("PATH", "")
        return [str(python)], env
    return [shutil.which("conda"), "run", "--no-capture-output", "-n", cfg.get("scenicplus_conda_env") or "scenicplus", "python"], env


def stage_fingerprint(name, params, upstream):
    """Parameters, resource file stamps, upstream fingerprints and stage code determine reuse."""
    code = [(f.name, hashlib.sha256(f.read_bytes()).hexdigest()) for f in (STAGE_DIR / f"{name}.py", STAGE_DIR / "common.py")]
    stamps = {k: common.file_stamp(params[k]) for k in RESOURCE_KEYS
              if k in params and k not in EXECUTION_KEYS and Path(str(params[k])).is_file()}
    payload = {"stage": name, "params": {k: v for k, v in params.items() if k not in EXECUTION_KEYS},
               "resources": stamps, "upstream": upstream, "code": code}
    return common.sha256_text(json.dumps(payload, sort_keys=True, default=str))


def run_stage(name, params, upstream, work, cfg, runner=None):
    """Run one stage script unless a completed result with the same fingerprint exists.

    Files left by an unfinished attempt are kept only when that attempt had the same
    fingerprint (e.g. Snakemake resumes after a crash); otherwise the stage directory is
    cleared, because timestamp-based reuse inside a stage cannot see parameter changes.
    """
    stage_dir = work / name
    stage_dir.mkdir(parents=True, exist_ok=True)
    params = {**params, "stage_dir": str(stage_dir)}
    fingerprint = stage_fingerprint(name, params, upstream)
    marker = stage_dir / "stage_result.json"
    if marker.is_file():
        previous = common.read_json(marker)
        outputs_present = all(Path(p).exists() for p in previous.get("outputs", {}).values())
        if previous.get("fingerprint") == fingerprint and previous.get("status") == "completed" and outputs_present:
            return {**previous, "reused": True}
    attempt = stage_dir / "attempt.json"
    same_attempt = attempt.is_file() and common.read_json(attempt).get("fingerprint") == fingerprint
    if not same_attempt and any(stage_dir.iterdir()):
        shutil.rmtree(stage_dir)
        stage_dir.mkdir(parents=True)
    common.write_json(attempt, {"fingerprint": fingerprint})
    marker.unlink(missing_ok=True)
    (stage_dir / common.STAGE_OUTPUT).unlink(missing_ok=True)
    params_file = stage_dir / "params.json"
    common.write_json(params_file, params)
    log = stage_dir / "stage.log"
    command, env = launcher(cfg)
    if runner is not None:
        returncode = runner(name, command + [str(STAGE_DIR / f"{name}.py"), str(params_file)], env, log)
    else:
        with log.open("w") as handle:
            returncode = subprocess.run(command + [str(STAGE_DIR / f"{name}.py"), str(params_file)],
                                        stdout=handle, stderr=subprocess.STDOUT, env=env).returncode
    if returncode != 0 or not (stage_dir / common.STAGE_OUTPUT).is_file():
        tail = log.read_text(errors="replace").splitlines()[-25:] if log.is_file() else []
        raise RuntimeError(f"SCENIC+ {name} failed (exit {returncode}); log: {log}\n" + "\n".join(tail))
    result = {**common.read_json(stage_dir / common.STAGE_OUTPUT), "status": "completed", "fingerprint": fingerprint}
    common.write_json(marker, result)
    return {**result, "reused": False}


def stage_parameters(cfg, handoff, work):
    """Parameters shared by all stages plus per-stage settings from the configuration."""
    temp_dir = cfg.get("scenicplus_temp_dir") or str(work / "tmp")
    Path(temp_dir).mkdir(parents=True, exist_ok=True)
    shared = {"sample_id": cfg["scenicplus_sample_id"], "n_cpu": int(_setting(cfg, "n_cpu")), "temp_dir": temp_dir,
              "seed": int(cfg.get("seed", 0))}
    resources = {stage: {k: cfg[RESOURCE_KEYS[k]] for k in keys} for stage, keys in STAGE_RESOURCES.items()}
    stage1 = {**shared, **resources["stage1_peaks"], "cells_tsv": handoff["cells_tsv"], **{k: _setting(cfg, k) for k in STAGE1_KEYS}}
    stage2 = {**shared, **resources["stage2_cistopic"], "cells_tsv": handoff["cells_tsv"], **{k: _setting(cfg, k) for k in STAGE2_KEYS},
              "min_cells": int(cfg.get("scenicplus_min_cells", 50)),
              "qc_minima": {"unique_fragments_in_peaks": cfg.get("scenicplus_qc_min_unique_fragments_in_peaks", 0),
                            "tss_enrichment": cfg.get("scenicplus_qc_min_tss_enrichment", 0),
                            "frip": cfg.get("scenicplus_qc_min_frip", 0)}}
    stage3 = {**shared, **resources["stage3_motif_databases"],
              "fraction_overlap_w_ctx_database": _setting(cfg, "fraction_overlap_w_ctx_database")}
    stage4 = {**shared, **resources["stage4_scenicplus"], **handoff, "species": _setting(cfg, "species"),
              "motif_annotation_version": _setting(cfg, "motif_annotation_version"),
              "allow_direct_only": bool(cfg.get("scenicplus_allow_direct_only", True)),
              **{k: _setting(cfg, k) for k in SNAKEMAKE_KEYS}}
    return stage1, stage2, stage3, stage4, temp_dir


def _aligned_activity(path, index, destination):
    import pandas as pd
    values = pd.read_csv(path, sep="\t", index_col=0)
    values.index = values.index.astype(str)
    aligned = values.reindex(index)
    aligned.index.name = "cell_id"
    aligned.to_csv(destination, sep="\t")
    return aligned, int(aligned.notna().any(axis=1).sum())


def run(ctx, out, cfg, source, result, runner=None):
    """Regulon task body for `regulon_method: scenicplus`; `result` builds the task result dict."""
    problems = prerequisites(cfg, source)
    if problems:
        return result("skipped", [source] if source else [], warnings=problems, actions=problems)
    import anndata as ad
    a = ad.read_h5ad(source)
    if "counts" not in a.layers or not a.uns.get("workflow_has_counts", False):
        return result("skipped", [source], warnings=["Verified raw counts are unavailable."],
                      actions=["Provide verified counts for regulatory preprocessing."])
    if cfg["scenicplus_cell_type_column"] not in a.obs:
        message = f"Cell-type column '{cfg['scenicplus_cell_type_column']}' is absent from the clustered dataset."
        return result("skipped", [source], warnings=[message], actions=["Provide cell-type labels (e.g. curated markers) before SCENIC+."])

    work = Path(cfg.get("scenicplus_work_dir") or out / "scenicplus_work")
    handoff, handoff_fingerprint, handoff_metrics = export_handoff(a, cfg, work / "handoff")
    if handoff_metrics["pseudobulk_groups"] < 2:
        message = ("At least two cell types with >= scenicplus_min_cells_per_cell_type cells (excluding "
                   f"{cfg.get('scenicplus_excluded_cell_types', DEFAULT_EXCLUDED_LABELS)}) are required for pseudobulk peak calling.")
        return result("skipped", [source], metrics={"handoff": handoff_metrics}, warnings=[message],
                      actions=["Annotate cell types before running SCENIC+."])

    stage1, stage2, stage3, stage4, temp_dir = stage_parameters(cfg, handoff, work)
    stages = {}
    stages["stage1_peaks"] = run_stage("stage1_peaks", stage1, [handoff_fingerprint], work, cfg, runner)
    peaks = stages["stage1_peaks"]
    stage2["consensus_regions"] = stage3["consensus_regions"] = peaks["outputs"]["consensus_regions"]
    stages["stage2_cistopic"] = run_stage("stage2_cistopic", stage2, [handoff_fingerprint, peaks["fingerprint"]], work, cfg, runner)
    stages["stage3_motif_databases"] = run_stage("stage3_motif_databases", stage3, [peaks["fingerprint"]], work, cfg, runner)
    topics = stages["stage2_cistopic"]
    stage4.update(cistopic_object=topics["outputs"]["cistopic_object"], region_set_folder=topics["outputs"]["region_set_folder"])
    stages["stage4_scenicplus"] = run_stage(
        "stage4_scenicplus", stage4,
        [handoff_fingerprint, topics["fingerprint"], stages["stage3_motif_databases"]["fingerprint"]], work, cfg, runner)
    final = stages["stage4_scenicplus"]["outputs"]

    outputs = {}
    for kind in ("direct", "extended"):
        if f"triplets_{kind}" not in final:
            continue
        edges = out / f"eregulon_triplets_{kind}.tsv"
        shutil.copyfile(final[f"triplets_{kind}"], edges)
        outputs["edges" if kind == "direct" else "edges_extended"] = str(edges)
    activity, scored = _aligned_activity(final["auc_gene_based_direct"], a.obs_names.astype(str), out / "eregulon_activity_gene_based_direct.tsv.gz")
    outputs["activity"] = str(out / "eregulon_activity_gene_based_direct.tsv.gz")
    _aligned_activity(final["auc_region_based_direct"], a.obs_names.astype(str), out / "eregulon_activity_region_based_direct.tsv.gz")
    outputs["activity_region_based"] = str(out / "eregulon_activity_region_based_direct.tsv.gz")

    summaries = []
    labels = a.obs[cfg["scenicplus_cell_type_column"]].astype(str)
    groups = [c for c in dict.fromkeys([cfg["scenicplus_cell_type_column"], "cluster", cfg.get("sample_column"),
                                        cfg.get("donor_column"), cfg.get("condition_column")]) if c and c in a.obs]
    for number, column in enumerate(groups):
        f = out / f"activity_group_{number}.tsv"
        activity.groupby(a.obs[column].astype(str).to_numpy(), observed=True).mean().to_csv(f, sep="\t")
        outputs[f"summary_{number}"] = str(f)
        summaries.append({"group_column": column, "path": str(f)})
    warnings = [w for s in stages.values() for w in s.get("warnings", [])]
    if scored < a.n_obs:
        warnings.append(f"{a.n_obs - scored} of {a.n_obs} RNA-QC cells have no eRegulon activity (failed ATAC QC); their activity is empty.")
    diagnostics = {"method": "scenicplus", "work_dir": str(work), "temp_dir": temp_dir, "handoff": handoff_metrics,
                   "stages": {name: {k: s.get(k) for k in ("fingerprint", "reused", "metrics", "warnings", "versions", "outputs")}
                              for name, s in stages.items()},
                   "summaries": summaries, "cell_type_labels": sorted(set(labels)), "interpretation": INTERPRETATION}
    outputs["diagnostics"] = common.write_json(out / "diagnostics.json", diagnostics)
    direct = stages["stage4_scenicplus"]["metrics"]["direct"]
    metrics = {"method": "scenicplus", "eregulons": direct["eregulons"], "tfs": direct["tfs"], "target_genes": direct["target_genes"],
               "regions": direct["regions"], "triplets": direct["triplets"], "cells_with_activity": scored,
               "consensus_peaks": peaks["metrics"]["consensus_peaks"],
               "contig_peaks_removed": peaks["metrics"]["consensus_peaks_removed_outside_keep_chromosomes"],
               "atac_rna_cells": topics["metrics"]["cells_used"], "selected_topics": topics["metrics"]["selected_topics"],
               "stages_reused": [n for n, s in stages.items() if s.get("reused")]}
    if direct["eregulons"] == 0:
        return result("inconclusive", [source], outputs, metrics, warnings + ["No direct eRegulons were inferred."])
    return result(inputs=[source] + [cfg[k] for k in RESOURCE_KEYS.values()], outputs=outputs, metrics=metrics,
                  warnings=[INTERPRETATION, "AUC activity is computed on the same cells used for inference and is not independent evidence of a cell state."] + warnings)
