"""Interactive answers use saved evidence and do not rerun clustering."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scrna_workflow.cli import main
from scrna_workflow.interactive import answer_question, run_chat


class InteractiveTests(unittest.TestCase):
    def test_follow_up_sends_compact_evidence_to_local_model(self):
        evidence = {"available_stages": {"clustering": "completed"},
                    "annotations": [{"cluster": "0", "label": "T cells"}],
                    "top_markers": {"0": [{"gene": "CD3D", "logfoldchange": 2.0}]}}
        observed = {}

        def local_chat(messages, model, format_schema):
            observed["messages"] = messages
            observed["model"] = model
            return {"answer": "Cluster 0 is labeled T cells, supported by CD3D."}

        answer = answer_question("What is cluster 0?", evidence,
                                 {"annotations": "/tmp/annotations.json"}, chat=local_chat)
        self.assertIn("T cells", answer)
        self.assertEqual(observed["model"], "qwen2.5:7b")
        self.assertIn("CD3D", observed["messages"][1]["content"])

    def test_unsupported_condition_question_does_not_call_model(self):
        evidence = {"available_stages": {"clustering": "completed"}}
        answer = answer_question("Was there differential expression?", evidence, {},
                                 chat=lambda *args, **kwargs: self.fail("Model should not be called"))
        self.assertIn("no completed discovery analysis", answer)

    def test_reopened_chat_records_follow_up(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            evidence = {"available_stages": {"clustering": "completed"},
                        "annotations": [{"cluster": "0", "label": "unknown"}]}
            responses = iter(["What is cluster 0?", "exit"])
            output = []
            with patch("scrna_workflow.interactive.load_evidence",
                       return_value=(root, evidence, {"annotations": str(root / "annotations.json")})):
                self.assertEqual(run_chat(root, input_fn=lambda prompt: next(responses),
                                          output_fn=output.append,
                                          chat=lambda *args, **kwargs: {"answer": "Cluster 0 is unknown."}), 0)
            transcript = [json.loads(line) for line in (root / "conversation.jsonl").read_text().splitlines()]
            self.assertEqual(len(transcript), 1)
            self.assertEqual(transcript[0]["question"], "What is cluster 0?")
            self.assertIn("unknown", output[-1])

    def test_fresh_chat_clusters_once_and_reopen_skips_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "counts.h5ad"
            data.touch()
            run_dir = root / "run"
            with patch("scrna_workflow.cli.execute", return_value=0) as execute, \
                    patch("scrna_workflow.interactive.run_chat", return_value=0) as chat:
                self.assertEqual(main(["chat", "--input", str(data),
                                       "--output", str(run_dir)]), 0)
                config = execute.call_args.args[0]
                self.assertEqual(config["requested_tasks"], ["clustering"])
                self.assertEqual(config["input_path"], str(data.resolve()))
                self.assertEqual(config["annotation_backend"], "ollama")
                self.assertFalse(config.get("knowledge_base_path"))
                self.assertEqual(chat.call_count, 1)
                execute.reset_mock()
                self.assertEqual(main(["chat", "--run-dir", str(run_dir)]), 0)
                execute.assert_not_called()


if __name__ == "__main__":
    unittest.main()
