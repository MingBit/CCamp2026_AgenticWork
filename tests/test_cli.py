"""Exercise the public CLI and marker annotation contract."""
import json
from io import BytesIO
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import anndata as ad
import pandas as pd

from scrna_workflow.cli import _marker_sets
from scrna_workflow.question import prepare_question
from scrna_workflow.llm_planner import FIELDS, plan_question as plan_llm_question
from scrna_workflow.planner import plan_question as plan_rule_question
from scrna_workflow.orchestrator import selected_tasks
from scrna_workflow.tools.annotation import annotate_clusters


class CliTests(unittest.TestCase):
    def test_plan_and_tools_expose_all_stages(self):
        for command in ("plan", "tools"):
            result = subprocess.run(
                [sys.executable, "-m", "scrna_workflow", command],
                check=True, capture_output=True, text=True,
            )
            rows = json.loads(result.stdout)
            self.assertEqual(len(rows), 9)
            self.assertEqual({row["task"] for row in rows}, {
                "inspection", "qc", "representation", "graph", "clustering",
                "regulon", "discovery", "validation", "report",
            })

    def test_marker_file_requires_multiple_genes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "markers.yaml"
            path.write_text("markers:\n  T cells: [CD3D, CD3E]\n")
            self.assertEqual(_marker_sets(path), {"T cells": ["CD3D", "CD3E"]})
            path.write_text("markers:\n  T cells: [CD3D, CD3D]\n")
            with self.assertRaisesRegex(ValueError, "at least two"):
                _marker_sets(path)

    def test_annotation_resolves_gene_symbols(self):
        data = ad.AnnData(
            obs=pd.DataFrame({"cluster": pd.Categorical(["0", "0"])}),
            var=pd.DataFrame({"gene_symbols": ["CD3D", "CD3E", "MS4A1"]},
                             index=["ENSG1", "ENSG2", "ENSG3"]),
        )
        ranked = pd.DataFrame({"group": ["0", "0", "0"],
                               "names": ["ENSG1", "ENSG2", "ENSG3"],
                               "logfoldchanges": [2.0, 1.0, -1.0]})
        result = annotate_clusters(data, ranked,
                                   {"T cells": ["CD3D", "CD3E"],
                                    "B cells": ["MS4A1", "CD79A"]})
        self.assertEqual(result[0]["label"], "T cells")
        self.assertEqual(data.obs["cell_type"].unique().tolist(), ["T cells"])

    def test_question_extracts_local_input_and_marker_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "cells.h5ad"
            markers = root / "markers.yaml"
            data.touch()
            markers.write_text("markers:\n  T cells: [CD3D, CD3E]\n")
            question = f"Please analyze {data} and annotate cells using {markers}"
            config = prepare_question(question, {"output_dir": str(root / "run")})
            self.assertEqual(config["input_path"], str(data.resolve()))
            self.assertEqual(config["marker_path"], str(markers.resolve()))
            self.assertEqual(config["biological_question"], question)

    def test_question_requires_unambiguous_dataset(self):
        with self.assertRaisesRegex(ValueError, "one dataset path"):
            prepare_question("Please analyze my scRNA-seq data", {})

    def test_local_llm_plan_selects_bounded_stages(self):
        with tempfile.TemporaryDirectory() as directory:
            data = Path(directory) / "cells.h5ad"
            data.touch()
            question = f"Please only do QC on {data}"
            proposal = {name: None for name in FIELDS}
            proposal.update(input_path=str(data), primary_comparison=None,
                            requested_analyses=["qc", "report"],
                            unsupported_requests=["differential_expression"])
            config = plan_llm_question(question, {}, chat=lambda messages, model: proposal)
            self.assertEqual(config["requested_tasks"], ["qc"])
            self.assertEqual(config["question_plan"]["unsupported_requests"], [])
            self.assertEqual(config["input_path"], str(data.resolve()))

    def test_local_llm_cannot_invent_input_path(self):
        with tempfile.TemporaryDirectory() as directory:
            data = Path(directory) / "invented.h5ad"
            data.touch()
            proposal = {name: None for name in FIELDS}
            proposal.update(input_path=str(data), primary_comparison=None,
                            requested_analyses=["qc"], unsupported_requests=[])
            with self.assertRaisesRegex(ValueError, "absent from the question"):
                plan_llm_question("Please do QC", {}, chat=lambda messages, model: proposal)

    def test_local_llm_does_not_reassign_configured_knowledge_base(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "cells.h5ad"
            knowledge = root / "human_pbmc_markers.yaml"
            data.touch()
            knowledge.touch()
            proposal = {name: None for name in FIELDS}
            proposal.update(input_path=str(data), marker_path=str(knowledge),
                            primary_comparison=None, requested_analyses=["annotation"],
                            unsupported_requests=[])
            config = plan_llm_question(
                "Analyze human PBMC cells and annotate clusters",
                {"input_path": str(data), "knowledge_base_path": str(knowledge)},
                chat=lambda messages, model: proposal,
            )
            self.assertEqual(config["knowledge_base_path"], str(knowledge))
            self.assertNotIn("marker_path", config)
            self.assertEqual(config["requested_tasks"], ["clustering"])

    def test_qc_only_plan_excludes_downstream_stages(self):
        plan = plan_rule_question("Analyze cells.h5ad and only do the QC", planner="rules")
        self.assertEqual(plan["targets"], ["qc"])
        self.assertEqual(list(selected_tasks(plan["targets"])), ["inspection", "qc"])

    def test_local_llm_can_only_choose_registered_tasks(self):
        payload = {"message": {"content": json.dumps({"targets": ["clustering"]})}}
        with patch("scrna_workflow.planner.request.urlopen", return_value=BytesIO(json.dumps(payload).encode())) as call:
            plan = plan_rule_question("Identify cell types", planner="ollama", model="local-model")
        self.assertEqual(plan["targets"], ["clustering"])
        self.assertEqual(list(selected_tasks(plan["targets"])),
                         ["inspection", "qc", "representation", "graph", "clustering"])
        request_body = json.loads(call.call_args.args[0].data)
        self.assertFalse(request_body["stream"])
        self.assertEqual(request_body["format"]["properties"]["targets"]["items"]["enum"],
                         list(selected_tasks(["report"])))
        with self.assertRaisesRegex(RuntimeError, "Cloud Ollama"):
            plan_rule_question("Identify cell types", planner="ollama", model="remote:cloud")


if __name__ == "__main__":
    unittest.main()
