"""SCENIC+ orchestration: handoff, stage reuse, result alignment and report wiring.

Stage scripts need the separate SCENIC+ environment, so a fake stage runner writes the
files each stage would produce; the tests check the workflow-side contract only.
"""
import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scrna_workflow import regulon_scenicplus as rs
from scrna_workflow.downstream import regulon, report
from scrna_workflow.tools._downstream_common import _result
from scrna_workflow.scenicplus_stages import common
from scrna_workflow.scenicplus_stages.stage2_cistopic import qc_thresholds
from scrna_workflow.scenicplus_stages.stage4_scenicplus import empty_extended_egrn_failure

LABELS = ["B cell"] * 20 + ["T/NK"] * 20 + ["unknown"] * 15 + ["rare"] * 3


def make_dataset(tmp_path):
    rng = np.random.default_rng(0)
    n = len(LABELS)
    counts = sparse.csr_matrix(rng.poisson(2, size=(n, 12)).astype(float))
    obs = pd.DataFrame({"annotation": LABELS, "cluster": [str(i % 3) for i in range(n)]},
                       index=[f"BC{i:03d}-1" for i in range(n)])
    a = ad.AnnData(counts.copy(), obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(12)]))
    a.layers["counts"] = counts
    a.uns["workflow_has_counts"] = True
    f = tmp_path / "annotated.h5ad"
    a.write_h5ad(f)
    return a, f


def make_config(tmp_path):
    resources = tmp_path / "resources"
    resources.mkdir()
    cfg = {"regulon_method": "scenicplus", "seed": 3, "scenicplus_python": str(tmp_path / "python"),
           "scenicplus_work_dir": str(tmp_path / "work"), "scenicplus_sample_id": "s1",
           "scenicplus_cell_type_column": "annotation", "scenicplus_min_cells_per_cell_type": 10,
           "scenicplus_n_cpu": 2, "scenicplus_species": "homo_sapiens", "scenicplus_motif_annotation_version": "v10nr_clust"}
    Path(cfg["scenicplus_python"]).write_text("")
    for key in rs.RESOURCE_KEYS.values():
        f = resources / key
        f.write_text(key)
        cfg[key] = str(f)
    Path(cfg["scenicplus_mallet_path"]).chmod(0o755)
    for key in rs.STAGE1_KEYS + rs.STAGE2_KEYS + rs.SNAKEMAKE_KEYS:
        cfg[f"scenicplus_{key}"] = [1, 2] if key.startswith(("search_space", "n_topics", "quantile", "top_n", "annotations")) else 1
    cfg["scenicplus_keep_chromosomes"] = ["chr1"]
    return cfg


class FakeStages:
    """Writes the outputs each stage script would produce; scored cells exclude the last two."""

    def __init__(self, barcodes, extended=True):
        self.calls, self.barcodes, self.extended = [], list(barcodes), extended

    def __call__(self, name, command, env, log):
        self.calls.append(name)
        params = json.loads(Path(command[-1]).read_text())
        stage = Path(params["stage_dir"])
        log.write_text("fake stage log\n")
        outputs, metrics = {}, {}
        if name == "stage1_peaks":
            outputs["consensus_regions"] = str(stage / "consensus_regions.bed")
            Path(outputs["consensus_regions"]).write_text("chr1\t100\t600\n")
            metrics = {"consensus_peaks": 1, "consensus_peaks_removed_outside_keep_chromosomes": 2}
        elif name == "stage2_cistopic":
            assert params["consensus_regions"].endswith("consensus_regions.bed")
            outputs = {"cistopic_object": str(stage / "cistopic_object.pkl"), "region_set_folder": str(stage / "region_sets")}
            Path(outputs["cistopic_object"]).write_text("pickle")
            Path(outputs["region_set_folder"]).mkdir(exist_ok=True)
            metrics = {"cells_used": len(self.barcodes) - 2, "selected_topics": 5}
        elif name == "stage4_scenicplus":
            assert Path(params["counts_npz"]).is_file() and Path(params["cistopic_object"]).is_file()
            scored = self.barcodes[:-2]
            for kind in ("direct", "extended") if self.extended else ("direct",):
                triplets = pd.DataFrame({"tf": ["TF1", "TF1", "TF2"], "target": ["G1", "G2", "G3"],
                                         "region": ["chr1:100-600"] * 3, "eregulon": ["TF1_+/+", "TF1_+/+", "TF2_+/+"],
                                         "importance_tf2g": [3.0, 2.0, 1.0], "rho_tf2g": [0.4, 0.3, 0.2]})
                outputs[f"triplets_{kind}"] = str(stage / f"eregulon_triplets_{kind}.tsv")
                triplets.to_csv(outputs[f"triplets_{kind}"], sep="\t", index=False)
                for modality in ("gene_based", "region_based"):
                    f = stage / f"auc_{modality}_{kind}.tsv.gz"
                    frame = pd.DataFrame({"TF1_+/+_(2g)": np.linspace(0, 1, len(scored))}, index=pd.Index(scored, name="cell_id"))
                    frame.iloc[::-1].to_csv(f, sep="\t")  # reversed order must be realigned
                    outputs[f"auc_{modality}_{kind}"] = str(f)
            metrics = {kind: {"eregulons": 2, "tfs": 2, "target_genes": 3, "regions": 1, "triplets": 3, "cells_scored": len(scored)}
                       for kind in ("direct", "extended") if self.extended or kind == "direct"}
        warnings = [] if self.extended or name != "stage4_scenicplus" else ["only direct-annotation eRegulons were built"]
        common.write_json(stage / common.STAGE_OUTPUT, {"outputs": outputs, "metrics": metrics, "warnings": warnings, "versions": {}})
        return 0


def run_regulon(tmp_path, cfg, a, source, fake):
    out = tmp_path / "run" / "regulon"
    out.mkdir(parents=True, exist_ok=True)
    return rs.run({"config": cfg}, out, cfg, str(source), _result, runner=fake)


def test_snakemake_overrides_and_missing_template_keys():
    cfg = make_config(Path(tempfile.mkdtemp()))
    params = {"sample_id": "s1", "n_cpu": 4, "seed": 1, "temp_dir": "/tmp/x", "cistopic_object": "c", "gex_anndata": "g",
              "region_set_folder": "r", "ctx_db": "d", "dem_db": "e", "motif_annotations": "m", "species": "homo_sapiens",
              "motif_annotation_version": "v10nr_clust",
              **{k: cfg[f"scenicplus_{k}"] for k in rs.SNAKEMAKE_KEYS}}
    overrides = common.snakemake_overrides(params)
    assert overrides[("params_data_preparation", "bc_transform_func")] == "\"lambda x: f'{x}___s1'\""
    assert overrides[("params_data_preparation", "species")] == "hsapiens"
    assert overrides[("params_inference", "quantile_thresholds_region_to_gene")] == "1 2"
    template = {section: {key: None for (s, key) in overrides if s == section} for section, _ in overrides}
    assert common.apply_overrides(template, overrides)["input_data"]["ctx_db_fname"] == "d"
    del template["params_inference"]["rho_threshold"]
    with pytest.raises(KeyError, match="rho_threshold"):
        common.apply_overrides(template, overrides)
    with pytest.raises(ValueError):
        common.snakemake_overrides({**params, "sample_id": "bad id"})


def test_regions_labels_and_qc_thresholds(tmp_path):
    assert common.parse_region("chr1:10-20") == ("chr1", 10, 20)
    assert common.write_region_bed(tmp_path / "s" / "r.bed", ["chr2:5-9", "chr1:10-20", "chr1:10-20"]) == 2
    assert (tmp_path / "s" / "r.bed").read_text() == "chr1\t10\t20\nchr2\t5\t9\n"
    assert common.unique_safe_labels(["T/NK", "T NK", "B"]) == {"B": "B", "T NK": "T_NK", "T/NK": "T_NK_2"}
    otsu = tmp_path / "otsu.tsv"
    otsu.write_text("unique_fragments_in_peaks_count_otsu_threshold\ttss_enrichment_otsu_threshold\n800.5\t7.2\n")
    minima = {"unique_fragments_in_peaks": 1000, "tss_enrichment": 5, "frip": 0}
    assert qc_thresholds(otsu, minima, True)[0] == {"unique_fragments_threshold": 1000.0, "tss_enrichment_threshold": 7.2, "frip_threshold": 0.0}
    assert qc_thresholds(otsu, minima, False)[0]["tss_enrichment_threshold"] == 5.0


def test_missing_resources_skip_without_running(tmp_path):
    a, source = make_dataset(tmp_path)
    ctx = {"root": tmp_path / "run", "config": {"regulon_method": "scenicplus"},
           "artifacts": {"clustering": {"outputs": {"dataset": str(source)}}}}
    result = regulon(ctx)
    assert result["status"] == "skipped"
    assert any("scenicplus_fragments_path" in w for w in result["warnings"])
    assert regulon({**ctx, "config": {"regulon_method": "other"}})["status"] == "skipped"


def test_handoff_stages_reuse_and_aligned_outputs(tmp_path):
    a, source = make_dataset(tmp_path)
    cfg = make_config(tmp_path)
    fake = FakeStages(a.obs_names)
    result = run_regulon(tmp_path, cfg, a, source, fake)
    assert result["status"] == "completed", result["warnings"]
    assert fake.calls == ["stage1_peaks", "stage2_cistopic", "stage3_motif_databases", "stage4_scenicplus"]

    cells = pd.read_table(Path(cfg["scenicplus_work_dir"]) / "handoff" / "cells.tsv", dtype=str, keep_default_na=False)
    assert dict(zip(cells["cell_type"], cells["pseudobulk_group"])) == {"B cell": "B_cell", "T/NK": "T_NK", "unknown": "", "rare": ""}
    counts = sparse.load_npz(Path(cfg["scenicplus_work_dir"]) / "handoff" / "counts.npz")
    assert (counts != sparse.csr_matrix(a.layers["counts"])).nnz == 0

    activity = pd.read_csv(result["outputs"]["activity"], sep="\t", index_col=0)
    assert list(activity.index) == list(a.obs_names)
    assert activity.iloc[:-2].notna().all().all() and activity.iloc[-2:].isna().all().all()
    assert activity.iloc[0, 0] == 0 and activity.iloc[-3, 0] == 1
    assert result["metrics"]["method"] == "scenicplus" and result["metrics"]["cells_with_activity"] == len(a) - 2
    assert any("no eRegulon activity" in w for w in result["warnings"])
    assert pd.read_csv(result["outputs"]["edges"], sep="\t")["tf"].tolist() == ["TF1", "TF1", "TF2"]

    # Identical inputs: every stage is reused. A peak-calling change reruns stage 1 and its dependents.
    fake.calls.clear()
    again = run_regulon(tmp_path, cfg, a, source, fake)
    assert fake.calls == [] and len(again["metrics"]["stages_reused"]) == 4
    run_regulon(tmp_path, {**cfg, "scenicplus_n_cpu": 64, "scenicplus_temp_dir": str(tmp_path / "other_tmp")}, a, source, fake)
    assert fake.calls == []
    stale = Path(cfg["scenicplus_work_dir"]) / "stage4_scenicplus" / "eRegulon_direct.tsv"
    stale.write_text("built with previous parameters")
    run_regulon(tmp_path, {**cfg, "scenicplus_peak_half_width": 300}, a, source, fake)
    assert fake.calls == ["stage1_peaks", "stage2_cistopic", "stage3_motif_databases", "stage4_scenicplus"]
    assert not stale.exists()  # outputs from a different fingerprint never leak into a rerun


def test_failed_attempt_keeps_files_only_for_same_parameters(tmp_path):
    a, source = make_dataset(tmp_path)
    cfg = make_config(tmp_path)
    fake = FakeStages(a.obs_names)
    partial = Path(cfg["scenicplus_work_dir"]) / "stage4_scenicplus" / "Snakemake" / "tf_to_gene_adj.tsv"

    def crash_in_stage4(name, command, env, log):
        if name != "stage4_scenicplus":
            return fake(name, command, env, log)
        partial.parent.mkdir(parents=True, exist_ok=True)
        partial.write_text("finished Snakemake rule")
        log.write_text("Error in rule eGRN_extended\n")
        return 1
    with pytest.raises(RuntimeError, match="stage4_scenicplus failed"):
        run_regulon(tmp_path, cfg, a, source, crash_in_stage4)
    seen = []

    def check_partial(name, command, env, log):
        seen.append(partial.exists())
        return fake(name, command, env, log)
    run_regulon(tmp_path, cfg, a, source, check_partial)  # same parameters: Snakemake may resume
    assert seen == [True]
    with pytest.raises(RuntimeError):  # new parameters: stage 4 reruns and crashes again
        run_regulon(tmp_path, {**cfg, "scenicplus_min_target_genes": 3}, a, source, crash_in_stage4)
    seen.clear()
    run_regulon(tmp_path, {**cfg, "scenicplus_min_target_genes": 4}, a, source, check_partial)
    assert seen == [False]  # different parameters than the crashed attempt: stale Snakemake outputs removed


# Tail of the SCENIC+ v1.0a2 failure seen on the smoke subset (job 6504017).
EMPTY_EXTENDED_LOG = """2026-09-16 18:04:54,673 SCENIC+      INFO     Formatting eGRN as table.
Traceback (most recent call last):
  File ".../scenicplus/cli/commands.py", line 790, in _format_egrns
    eRegulon_metadata = pd.concat(eRegulons_formatted)
    raise ValueError("No objects to concatenate")
ValueError: No objects to concatenate
[Wed Sep 16 18:04:57 2026]
Error in rule eGRN_extended:
    jobid: 13
    output: eRegulons_extended.tsv
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
WorkflowError:
At least one job did not complete successfully.""".splitlines()


def test_direct_only_fallback_detection_and_results(tmp_path):
    assert empty_extended_egrn_failure(EMPTY_EXTENDED_LOG)
    assert not empty_extended_egrn_failure([l.replace("eGRN_extended", "eGRN_direct") for l in EMPTY_EXTENDED_LOG])
    assert not empty_extended_egrn_failure([l.replace("No objects to concatenate", "MemoryError") for l in EMPTY_EXTENDED_LOG])
    assert not empty_extended_egrn_failure(EMPTY_EXTENDED_LOG + ["Error in rule AUCell_direct:"])

    a, source = make_dataset(tmp_path)
    cfg = make_config(tmp_path)
    result = run_regulon(tmp_path, cfg, a, source, FakeStages(a.obs_names, extended=False))
    assert result["status"] == "completed"
    assert "edges" in result["outputs"] and "edges_extended" not in result["outputs"]
    assert any("only direct-annotation" in w for w in result["warnings"])
    params = json.loads((Path(cfg["scenicplus_work_dir"]) / "stage4_scenicplus" / "params.json").read_text())
    assert params["allow_direct_only"] is True


def test_single_cell_type_and_stage_failure(tmp_path):
    a, source = make_dataset(tmp_path)
    cfg = make_config(tmp_path)
    one_type = a.copy()
    one_type.obs["annotation"] = "unknown"
    one_type.obs.loc[one_type.obs_names[:20], "annotation"] = "B cell"
    single = tmp_path / "single.h5ad"
    one_type.write_h5ad(single)
    skipped = run_regulon(tmp_path, cfg, one_type, single, FakeStages(a.obs_names))
    assert skipped["status"] == "skipped" and "two cell types" in skipped["warnings"][0]

    def failing(name, command, env, log):
        log.write_text("Traceback: MACS2 exploded\n")
        return 1
    with pytest.raises(RuntimeError, match="(?s)stage1_peaks failed.*MACS2 exploded"):
        run_regulon(tmp_path, cfg, a, source, failing)
    assert not (Path(cfg["scenicplus_work_dir"]) / "stage1_peaks" / "stage_result.json").exists()


def test_report_uses_eregulon_outputs(tmp_path):
    a, source = make_dataset(tmp_path)
    cfg = make_config(tmp_path)
    regulon_result = run_regulon(tmp_path, cfg, a, source, FakeStages(a.obs_names))
    ctx = {"root": tmp_path / "run", "config": {"organism": "human"},
           "artifacts": {"clustering": {"status": "completed", "outputs": {"dataset": str(source)}}, "regulon": regulon_result}}
    result = report(ctx)
    names = {Path(p).name for k, p in result["outputs"].items() if k.startswith("figure_")}
    assert {"eregulon_activity_heatmap.png", "eregulon_tf_network.png"} <= names
    ledger = pd.read_csv(result["outputs"]["evidence_ledger"], sep="\t")
    assert ledger["finding"].str.contains("SCENIC\\+ eRegulons").any()


class ScenicPlusTests(unittest.TestCase):
    def test_scenicplus_contract(self):
        for check in (test_regions_labels_and_qc_thresholds, test_missing_resources_skip_without_running,
                      test_handoff_stages_reuse_and_aligned_outputs, test_failed_attempt_keeps_files_only_for_same_parameters,
                      test_direct_only_fallback_detection_and_results, test_single_cell_type_and_stage_failure,
                      test_report_uses_eregulon_outputs):
            with self.subTest(check=check.__name__), tempfile.TemporaryDirectory() as directory:
                check(Path(directory))
        test_snakemake_overrides_and_missing_template_keys()
