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
# Ray creates unix sockets under the temp directory: <temp_dir>/session_<40 chars>/sockets/plasma_store
# must fit AF_UNIX's 107-byte limit, which leaves ~45 characters for the temp directory itself.
MAX_TEMP_DIR_LENGTH = 45
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


# (modality, kind) -> (stage 4 output key, run output key, file stem); the first entry is the primary activity.
ACTIVITY_TABLES = [
    (("gene_based", "direct"), "auc_gene_based_direct", "activity", "eregulon_activity_gene_based_direct"),
    (("region_based", "direct"), "auc_region_based_direct", "activity_region_based", "eregulon_activity_region_based_direct"),
    (("gene_based", "extended"), "auc_gene_based_extended", "activity_extended", "eregulon_activity_gene_based_extended"),
    (("region_based", "extended"), "auc_region_based_extended", "activity_region_based_extended", "eregulon_activity_region_based_extended"),
]


def _aligned_activity(path, index, destination):
    import pandas as pd
    values = pd.read_csv(path, sep="\t", index_col=0)
    values.index = values.index.astype(str)
    aligned = values.reindex(index)
    aligned.index.name = "cell_id"
    aligned.to_csv(destination, sep="\t")
    return aligned, int(aligned.notna().any(axis=1).sum())


def group_activity_summary(activity, groups):
    """Mean AUC per group over scored cells, with group size and number of scored cells first."""
    import pandas as pd
    groups = pd.Series(groups, index=activity.index, name="group").astype(str)
    scored = activity.notna().any(axis=1)
    counts = pd.DataFrame({"n_cells": groups.value_counts(), "n_cells_scored": groups[scored].value_counts()})
    counts = counts.fillna(0).astype(int)
    means = activity.groupby(groups.to_numpy(), observed=True).mean()
    summary = counts.join(means)
    summary.index.name = "group"
    return summary.sort_index()


def regulon_specificity_scores(activity, groups):
    """Regulon specificity scores (Suo et al. 2018), as in `scenicplus.RSS.regulon_specificity_scores_df`.

    RSS = 1 - Jensen-Shannon distance between an eRegulon's AUC distribution over scored
    cells (normalized to sum 1) and a group's membership indicator (normalized to sum 1).
    Cells without AUC (failed ATAC QC) are ignored; rows are groups, columns eRegulons.
    """
    import numpy as np
    import pandas as pd
    from scipy.spatial.distance import jensenshannon

    groups = pd.Series(groups, index=activity.index).astype(str)
    scored = activity.notna().any(axis=1)
    labels = sorted(groups[scored].unique())
    scores = pd.DataFrame(np.nan, index=pd.Index(labels, name="group"), columns=activity.columns)
    for regulon in activity.columns:
        valid = activity[regulon].notna().to_numpy()
        values = activity[regulon].to_numpy(dtype=float)[valid]
        members = groups.to_numpy()[valid]
        if values.sum() <= 0:
            continue
        for label in labels:
            indicator = (members == label).astype(float)
            if indicator.sum() > 0:
                scores.loc[label, regulon] = 1.0 - jensenshannon(values / values.sum(), indicator / indicator.sum())
    return scores


def _top_regulons(scores, top_n):
    return list(dict.fromkeys(name for group in scores.index for name in scores.loc[group].dropna().nlargest(top_n).index))


def plot_rss_ranks(scores, destination, top_n=5, title="eRegulon specificity"):
    """One panel per group: RSS of all eRegulons by rank, top eRegulons labelled."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    ncols = min(3, len(scores.index))
    nrows = int(np.ceil(len(scores.index) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 3.6 * nrows), squeeze=False)
    for ax, group in zip(axes.ravel(), scores.index):
        values = scores.loc[group].dropna().sort_values(ascending=False)
        ranks = np.arange(1, len(values) + 1)
        ax.scatter(ranks, values.to_numpy(), s=8, color="#8497ad", rasterized=True)
        top = values.head(top_n)
        ax.scatter(ranks[:len(top)], top.to_numpy(), s=18, color="#c0504d", zorder=3)
        # Labels in a column at the upper right, in rank order, so leader lines do not cross.
        for i, (name, score) in enumerate(top.items()):
            ax.annotate(str(name), (i + 1, score), xytext=(0.42, 0.94 - 0.085 * i), textcoords="axes fraction",
                        fontsize=7, va="center", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
                        arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5))
        ax.set(title=str(group), xlabel="eRegulon rank", ylabel="RSS")
    for ax in axes.ravel()[len(scores.index):]:
        ax.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(destination, dpi=150)
    plt.close(fig)
    return str(destination)


def plot_rss_heatmap(scores, destination, top_n=5, title="eRegulon specificity (RSS)"):
    """Groups x union of each group's top eRegulons by RSS."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chosen = _top_regulons(scores, top_n)
    data = scores[chosen]
    fig, ax = plt.subplots(figsize=(max(5, 0.32 * len(chosen) + 2), max(2.8, 0.38 * len(data.index) + 1.8)))
    image = ax.imshow(data.to_numpy(dtype=float), cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(chosen)), chosen, rotation=90, fontsize=7)
    ax.set_yticks(range(len(data.index)), data.index)
    ax.set(title=title)
    fig.colorbar(image, ax=ax, label="RSS")
    fig.tight_layout()
    fig.savefig(destination, dpi=150)
    plt.close(fig)
    return str(destination)


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
    if len(temp_dir) > MAX_TEMP_DIR_LENGTH:
        message = (f"scenicplus_temp_dir is {len(temp_dir)} characters ({temp_dir}); Ray's unix sockets need "
                   f"at most {MAX_TEMP_DIR_LENGTH}. Set a short local path, e.g. /tmp/scenicplus.")
        return result("skipped", [source], metrics={"handoff": handoff_metrics}, warnings=[message], actions=[message])
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
    tables = {}
    for (modality, kind), stage_key, output_key, stem in ACTIVITY_TABLES:
        if stage_key not in final:  # extended tables are absent when only direct eRegulons were built
            continue
        destination = out / f"{stem}.tsv.gz"
        tables[(modality, kind)], n_scored = _aligned_activity(final[stage_key], a.obs_names.astype(str), destination)
        outputs[output_key] = str(destination)
        if (modality, kind) == ("gene_based", "direct"):
            scored = n_scored

    # Per-group summaries for every activity table. Direct gene-based keeps the original
    # activity_group_<n>.tsv names; the others add "_<modality>_<kind>".
    summaries = []
    labels = a.obs[cfg["scenicplus_cell_type_column"]].astype(str)
    groups = [c for c in dict.fromkeys([cfg["scenicplus_cell_type_column"], "cluster", cfg.get("sample_column"),
                                        cfg.get("donor_column"), cfg.get("condition_column")]) if c and c in a.obs]
    for number, column in enumerate(groups):
        for (modality, kind), activity in tables.items():
            primary = (modality, kind) == ("gene_based", "direct")
            suffix = "" if primary else f"_{modality}_{kind}"
            f = out / f"activity_group_{number}{suffix}.tsv"
            group_activity_summary(activity, a.obs[column].to_numpy()).to_csv(f, sep="\t")
            outputs[f"summary_{number}{suffix}"] = str(f)
            summaries.append({"group_column": column, "modality": modality, "eregulon_kind": kind, "path": str(f)})

    # Regulon specificity scores per group, with rank plots and heatmaps for the cell-type column.
    rss_files, top_rss, cell_type_rss = [], {}, None
    top_n = int(cfg.get("scenicplus_rss_top_n", 5))
    for number, column in enumerate(groups):
        for (modality, kind), activity in tables.items():
            primary = (modality, kind) == ("gene_based", "direct")
            suffix = "" if primary else f"_{modality}_{kind}"
            scores = regulon_specificity_scores(activity, a.obs[column].to_numpy())
            f = out / f"rss_group_{number}{suffix}.tsv"
            scores.to_csv(f, sep="\t")
            outputs[f"rss_{number}{suffix}"] = str(f)
            entry = {"group_column": column, "modality": modality, "eregulon_kind": kind, "path": str(f)}
            if column == cfg["scenicplus_cell_type_column"] and scores.size:
                label = f"{modality.replace('_', '-')} {kind} eRegulons"
                entry["rank_plot"] = outputs[f"rss_rank_plot{suffix}"] = plot_rss_ranks(
                    scores, out / f"rss_ranks{suffix}.png", top_n, f"eRegulon specificity per {column} ({label})")
                entry["heatmap"] = outputs[f"rss_heatmap{suffix}"] = plot_rss_heatmap(
                    scores, out / f"rss_heatmap{suffix}.png", top_n, f"RSS per {column}, top {top_n} per group ({label})")
                if primary:
                    top_rss = {group: scores.loc[group].dropna().nlargest(3).round(4).to_dict() for group in scores.index}
                    cell_type_rss = scores
            rss_files.append(entry)
    # Cell-type views of the direct eRegulon network (TF -> region -> target gene).
    network_summary = []
    warnings = [w for s in stages.values() for w in s.get("warnings", [])]
    if cell_type_rss is not None and "edges" in outputs and peaks["outputs"].get("macs2_dir"):
        from .regulon_networks import build_cell_type_networks
        network_outputs, network_summary, network_warnings = build_cell_type_networks(
            outputs["edges"], cell_type_rss, handoff["cells_tsv"], peaks["outputs"]["macs2_dir"], a,
            cfg["scenicplus_cell_type_column"], out / "networks",
            top_eregulons=int(cfg.get("scenicplus_network_top_eregulons", 5)),
            min_gene_fraction=float(cfg.get("scenicplus_network_min_gene_fraction", 0.1)),
            max_targets_per_tf=int(cfg.get("scenicplus_network_max_targets_per_tf", 20)))
        outputs.update(network_outputs)
        warnings += network_warnings
    if scored < a.n_obs:
        warnings.append(f"{a.n_obs - scored} of {a.n_obs} RNA-QC cells have no eRegulon activity (failed ATAC QC); their activity is empty.")
    diagnostics = {"method": "scenicplus", "work_dir": str(work), "temp_dir": temp_dir, "handoff": handoff_metrics,
                   "stages": {name: {k: s.get(k) for k in ("fingerprint", "reused", "metrics", "warnings", "versions", "outputs")}
                              for name, s in stages.items()},
                   "summaries": summaries, "regulon_specificity": rss_files, "cell_type_networks": network_summary,
                   "cell_type_labels": sorted(set(labels)),
                   "interpretation": INTERPRETATION}
    outputs["diagnostics"] = common.write_json(out / "diagnostics.json", diagnostics)
    direct = stages["stage4_scenicplus"]["metrics"]["direct"]
    extended = stages["stage4_scenicplus"]["metrics"].get("extended") or {}
    metrics = {"method": "scenicplus", "eregulons": direct["eregulons"], "tfs": direct["tfs"], "target_genes": direct["target_genes"],
               "regions": direct["regions"], "triplets": direct["triplets"], "cells_with_activity": scored,
               "extended_eregulons": extended.get("eregulons"), "extended_triplets": extended.get("triplets"),
               "top_rss_per_cell_type": top_rss,
               "cell_type_networks": {row["cell_type"]: {k: row[k] for k in ("tfs", "regions", "target_genes", "edges")}
                                      for row in network_summary},
               "consensus_peaks": peaks["metrics"]["consensus_peaks"],
               "contig_peaks_removed": peaks["metrics"]["consensus_peaks_removed_outside_keep_chromosomes"],
               "atac_rna_cells": topics["metrics"]["cells_used"], "selected_topics": topics["metrics"]["selected_topics"],
               "stages_reused": [n for n, s in stages.items() if s.get("reused")]}
    if direct["eregulons"] == 0:
        return result("inconclusive", [source], outputs, metrics, warnings + ["No direct eRegulons were inferred."])
    return result(inputs=[source] + [cfg[k] for k in RESOURCE_KEYS.values()], outputs=outputs, metrics=metrics,
                  warnings=[INTERPRETATION, "AUC activity is computed on the same cells used for inference and is not independent evidence of a cell state.",
                            "RSS ranks eRegulons by how concentrated their activity is in a group relative to all scored cells; it is descriptive, depends on the group definitions, and is not a statistical test."] + warnings)
