"""Integrity checks at the cell-to-sample boundary; no biological claims."""
import json
import tempfile
import unittest
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scrna_workflow.downstream import discovery, validation


def make_context(tmp_path, verified=True):
    x = sparse.csr_matrix(np.arange(1, 49).reshape(12, 4))
    obs = pd.DataFrame({"sample": np.repeat(["s1", "s2", "s3"], 4),
                        "donor": np.repeat(["d1", "d2", "d3"], 4),
                        "condition": np.repeat(["A", "A", "B"], 4),
                        "cluster": ["0", "1"] * 6}, index=[f"c{i}" for i in range(12)])
    a = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=["g1", "g2", "g3", "g4"]))
    a.layers["counts"] = x.copy()
    a.uns["workflow_has_counts"] = verified
    f = tmp_path / "data.h5ad"
    a.write_h5ad(f)
    ctx = {"root": tmp_path / "run", "config": {"sample_column": "sample", "donor_column": "donor", "condition_column": "condition"},
           "artifacts": {"inspection": {"outputs": {"dataset": str(f)}}, "qc": {"outputs": {"dataset": str(f)}}, "clustering": {"outputs": {"dataset": str(f)}}}}
    return ctx, a


def test_pseudobulk_preserves_each_sample_population_sum(tmp_path):
    ctx, a = make_context(tmp_path)
    result = discovery(ctx)
    pb = ad.read_h5ad(result["outputs"]["pseudobulk"])
    assert pb.X.sum() == a.layers["counts"].sum()
    for name, row in pb.obs.iterrows():
        cells = (a.obs["sample"] == row["sample"]) & (a.obs["cluster"] == row["population"])
        np.testing.assert_array_equal(pb[name].X.toarray().ravel(), np.asarray(a.layers["counts"][cells.to_numpy()].sum(axis=0)).ravel())
    props = pd.read_csv(result["outputs"]["abundance"], sep="\t", index_col=0)
    np.testing.assert_allclose(props.sum(axis=1), 1)
    assert "differential_expression" not in result["outputs"]
    assert json.loads(Path(result["outputs"]["de_diagnostics"]).read_text())["status"] == "skipped"


def test_unverified_counts_are_not_exported(tmp_path):
    ctx, _ = make_context(tmp_path, verified=False)
    result = discovery(ctx)
    assert "pseudobulk" not in result["outputs"]


def test_validation_catches_count_corruption(tmp_path):
    ctx, a = make_context(tmp_path)
    a.layers["counts"][0, 0] += 1
    f = tmp_path / "corrupted.h5ad"
    a.write_h5ad(f)
    ctx["artifacts"]["clustering"]["outputs"]["dataset"] = str(f)
    result = validation(ctx)
    assert result["status"] == "failed"
    assert not result["metrics"]["checks"]["counts_preserved_since_qc"]
    assert not result["metrics"]["checks"]["retained_counts_preserved_since_inspection"]


class DownstreamTests(unittest.TestCase):
    def test_integrity_and_design_guards(self):
        for check in (test_pseudobulk_preserves_each_sample_population_sum,
                      test_unverified_counts_are_not_exported, test_validation_catches_count_corruption):
            with self.subTest(check=check.__name__), tempfile.TemporaryDirectory() as directory:
                check(Path(directory))
