"""Safety gates for reference-free, provisional marker annotation."""
import unittest

import anndata as ad
import pandas as pd

from scrna_workflow.agents.annotation import OllamaAnnotationAgent


class AnnotationAgentTests(unittest.TestCase):
    def setUp(self):
        self.adata = ad.AnnData(
            obs=pd.DataFrame({"cluster": pd.Categorical(["0", "0"])},
                             index=["cell1", "cell2"]),
            var=pd.DataFrame(index=["CD3D", "CD3E", "MS4A1"]),
        )
        self.markers = pd.DataFrame({
            "group": ["0", "0", "0"],
            "names": ["CD3D", "CD3E", "MS4A1"],
            "logfoldchanges": [1.8, 1.2, -0.4],
        })

    def test_accepts_only_observed_multi_gene_support(self):
        response = {"annotations": [{"cluster": "0", "label": "T cells",
                                     "supporting_genes": ["CD3D", "CD3E"],
                                     "reason": "Two T-lineage markers"}]}
        agent = OllamaAnnotationAgent(chat=lambda *args, **kwargs: response)
        annotations, audit = agent.annotate_from_markers(self.adata, self.markers)
        self.assertEqual(annotations[0]["label"], "T cells")
        self.assertTrue(audit[0]["accepted"])
        self.assertLessEqual(annotations[0]["confidence"], 0.45)

    def test_rejects_unobserved_or_single_gene_support(self):
        response = {"annotations": [{"cluster": "0", "label": "B cells",
                                     "supporting_genes": ["MS4A1", "CD79A"],
                                     "reason": "unsupported"}]}
        agent = OllamaAnnotationAgent(chat=lambda *args, **kwargs: response)
        annotations, audit = agent.annotate_from_markers(self.adata, self.markers)
        self.assertEqual(annotations[0]["label"], "unknown")
        self.assertFalse(audit[0]["accepted"])


if __name__ == "__main__":
    unittest.main()
