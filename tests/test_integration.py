"""End-to-end regression using the public PBMC3k count dataset when available."""
import json
import os
from pathlib import Path
import tempfile
import unittest

import anndata as ad
from scipy import sparse
import yaml

from scrna_workflow.orchestrator import digest, execute


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_INPUT = ROOT / "results/pbmc3k/pbmc3k_counts_annotated.h5ad"
PUBLIC_CONFIG = ROOT / "results/pbmc3k/config.yaml"


@unittest.skipUnless(
    os.environ.get("SCRNA_RUN_PBMC_INTEGRATION") == "1"
    and PUBLIC_INPUT.is_file() and PUBLIC_CONFIG.is_file(),
    "Set SCRNA_RUN_PBMC_INTEGRATION=1 with the public PBMC3k files to rerun the full pipeline",
)
class PublicPbmcIntegrationTests(unittest.TestCase):
    def test_public_pbmc_run_and_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            config = yaml.safe_load(PUBLIC_CONFIG.read_text())
            config["input_path"] = str(PUBLIC_INPUT)
            config["output_dir"] = str(Path(directory) / "run")
            config["annotation_backend"] = "markers"
            self.assertEqual(execute(config), 0)
            root = Path(config["output_dir"])
            state = json.loads((root / "run_state.json").read_text())
            self.assertEqual(state["status"], "completed")
            cluster = state["tasks"]["clustering"]
            self.assertTrue(Path(cluster["outputs"]["html_report"]).is_file())
            self.assertTrue(Path(cluster["outputs"]["annotations"]).is_file())
            final = ad.read_h5ad(cluster["outputs"]["dataset"])
            original = ad.read_h5ad(PUBLIC_INPUT)
            aligned = original[final.obs_names, final.var_names]
            self.assertEqual((sparse.csr_matrix(aligned.X) -
                              sparse.csr_matrix(final.layers["counts"])).nnz, 0)
            self.assertTrue(all(matrix.shape == (final.n_obs, final.n_obs)
                                for matrix in final.obsp.values()))
            hashes = {path: digest(path) for result in state["tasks"].values()
                      for path in result.get("artifact_hashes", {})}
            self.assertEqual(execute(config, resume=True), 0)
            self.assertEqual(hashes, {path: digest(path) for path in hashes})
            self.assertIn("checkpoint_reused", (root / "events.jsonl").read_text())


if __name__ == "__main__":
    unittest.main()
