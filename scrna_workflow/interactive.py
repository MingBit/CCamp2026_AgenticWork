"""Evidence-grounded follow-up questions for an existing clustering run."""
from datetime import datetime, timezone
import json
from pathlib import Path

import pandas as pd

from .llm_planner import _ollama_chat
from .orchestrator import valid_checkpoint


ANSWER_SCHEMA = {
    "type": "object",
    "properties": {"answer": {"type": "string"}},
    "required": ["answer"],
    "additionalProperties": False,
}


def load_evidence(run_dir):
    """Read compact, versioned results without sending the expression matrix to a model."""
    root = Path(run_dir).expanduser().resolve()
    state_file = root / "run_state.json"
    if not state_file.is_file():
        raise ValueError(f"No workflow run exists at {root}")
    state = json.loads(state_file.read_text(encoding="utf-8"))
    tasks = state.get("tasks", {})
    cluster = tasks.get("clustering", {})
    if cluster.get("status") != "completed" or not valid_checkpoint(cluster):
        raise ValueError("A completed, intact clustering checkpoint is required for chat")
    outputs = cluster["outputs"]
    annotations = json.loads(Path(outputs["annotations"]).read_text(encoding="utf-8"))
    markers = pd.read_csv(outputs["markers"])
    top_markers = {}
    if not markers.empty:
        for group, rows in markers.groupby("group", sort=False):
            top_markers[str(group)] = [
                {"gene": str(row.names),
                 "logfoldchange": round(float(row.logfoldchanges), 3)}
                for row in rows.head(10).itertuples()
            ]
    import anndata as ad
    dataset = ad.read_h5ad(outputs["dataset"], backed="r")
    try:
        sizes = dataset.obs["cluster"].astype(str).value_counts().to_dict()
    finally:
        dataset.file.close()
    audit_path = outputs.get("llm_annotation_audit")
    audit = json.loads(Path(audit_path).read_text(encoding="utf-8")) if audit_path else None
    evidence = {
        "run_status": state.get("status"),
        "available_stages": {name: result.get("status") for name, result in tasks.items()},
        "input": tasks.get("inspection", {}).get("metrics", {}),
        "qc": tasks.get("qc", {}).get("metrics", {}),
        "clustering": cluster.get("metrics", {}),
        "cluster_sizes": {str(key): int(value) for key, value in sizes.items()},
        "annotations": annotations,
        "top_markers": top_markers,
        "annotation_audit": audit,
        "limitations": ["Cluster markers are exploratory cell-level rankings.",
                        "Cell identity and model explanations are computational inferences.",
                        "No condition-level inference is supported by clustering alone."],
    }
    sources = {key: outputs[key] for key in
               ("dataset", "markers", "annotations", "html_report", "llm_annotation_audit")
               if key in outputs}
    return root, evidence, sources


def answer_question(question, evidence, sources, model="qwen2.5:7b", history=None, chat=None):
    """Answer from recorded evidence; never request an unrun analysis implicitly."""
    question = question.strip()
    if not question:
        raise ValueError("Enter a question about the completed run")
    if model.lower().endswith(":cloud"):
        raise ValueError("Cloud Ollama models are disabled for local analysis chat")
    lower = question.lower()
    if any(term in lower for term in ("differential expression", "deg", "condition difference",
                                      "trajectory", "pseudotime", "ligand receptor")):
        if evidence["available_stages"].get("discovery") != "completed":
            return ("This session has clustering results but no completed discovery analysis. "
                    "That question needs appropriate sample and donor metadata and a separate "
                    "validated analysis run.")
    chat = chat or _ollama_chat
    messages = [
        {"role": "system", "content": (
            "Answer the user's follow-up about one completed scRNA-seq run using only "
            "the evidence JSON in this message. Prior chat is conversational context, "
            "not scientific evidence. Cite cluster IDs and supporting gene names when relevant. "
            "State when evidence is absent or ambiguous. Do not invent cell identities, "
            "condition effects, statistical significance, causal links, or new analyses. "
            "Annotation confidence is a heuristic marker-overlap score, not a calibrated probability. "
            "If annotation_backend is ollama_provisional, call every assigned identity provisional. "
            "Return a concise JSON object with one answer string."
        )},
        {"role": "user", "content": json.dumps({"evidence": evidence,
                                                "artifact_paths": sources})},
    ]
    for turn in (history or [])[-3:]:
        messages.append({"role": "user", "content": turn["question"]})
        messages.append({"role": "assistant", "content": turn["answer"]})
    messages.append({"role": "user", "content": question})
    result = chat(messages, model, format_schema=ANSWER_SCHEMA)
    if not isinstance(result, dict) or not isinstance(result.get("answer"), str) or not result["answer"].strip():
        raise ValueError("Local Ollama returned no usable answer")
    return result["answer"].strip()


def run_chat(run_dir, model="qwen2.5:7b", input_fn=input, output_fn=print, chat=None):
    """Open or resume a terminal conversation without rerunning analysis stages."""
    root, evidence, sources = load_evidence(run_dir)
    transcript = root / "conversation.jsonl"
    history = []
    if transcript.is_file():
        history = [json.loads(line) for line in transcript.read_text(encoding="utf-8").splitlines()]
    output_fn(f"Clustering ready: {len(evidence['annotations'])} clusters. Ask about this run; type exit to leave.")
    while True:
        try:
            question = input_fn("scRNA> ").strip()
        except (EOFError, KeyboardInterrupt):
            output_fn("")
            break
        if question.lower() in {"exit", "quit", ":q"}:
            break
        if not question:
            continue
        try:
            answer = answer_question(question, evidence, sources, model, history, chat)
        except (OSError, ValueError, RuntimeError) as exc:
            output_fn(f"Could not answer: {exc}")
            continue
        output_fn(answer)
        record = {"time": datetime.now(timezone.utc).isoformat(),
                  "question": question, "answer": answer, "model": model}
        with transcript.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")
        history.append(record)
    return 0
