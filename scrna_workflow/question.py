"""Conservative local interpretation and evidence-grounded summary of a CLI question.

This is a deterministic interface, not an unrestricted language model. Paths and
scientific settings are taken from explicit input, never guessed from biology.
"""
from datetime import datetime, timezone
import json
from pathlib import Path
import re


DATA_FILE = re.compile(r"(?<!\w)(?:~?/|\.\.?/)?[^\s'\";,]+\.(?:h5ad|h5|hdf5)\b", re.I)
MARKER_FILE = re.compile(r"(?<!\w)(?:~?/|\.\.?/)?[^\s'\";,]+\.(?:yaml|yml|json)\b", re.I)


def prepare_question(question, config):
    """Fill only unambiguous local paths from a question; retain CLI/config priority."""
    question = question.strip()
    if not question:
        raise ValueError("Write a scRNA-seq analysis question after 'ask'.")
    config = dict(config)
    if not config.get("input_path"):
        paths = [Path(hit).expanduser().resolve() for hit in DATA_FILE.findall(question)]
        paths = list(dict.fromkeys(paths))
        if len(paths) != 1:
            raise ValueError("Specify one dataset path in the question or with --input (.h5ad or 10x .h5).")
        config["input_path"] = str(paths[0])
    if not Path(config["input_path"]).exists():
        raise ValueError(f"Dataset does not exist: {config['input_path']}")
    if not config.get("marker_path") and not config.get("markers"):
        candidates = [Path(hit).expanduser().resolve() for hit in MARKER_FILE.findall(question)]
        candidates = [path for path in candidates if "marker" in path.name.lower()]
        if len(candidates) == 1:
            config["marker_path"] = str(candidates[0])
    if not config.get("output_dir"):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        config["output_dir"] = str((Path.cwd() / "results" / f"question_{stamp}").resolve())
    config["biological_question"] = question
    return config


def summarize_run(output_dir, question):
    """Build a compact summary solely from this run's recorded results."""
    root = Path(output_dir)
    state = json.loads((root / "run_state.json").read_text(encoding="utf-8"))
    tasks = state.get("tasks", {})
    get = lambda task: tasks.get(task, {})
    inspection = get("inspection").get("metrics", {})
    qc = get("qc").get("metrics", {})
    cluster = get("clustering")
    discovery = get("discovery")
    lines = ["# scRNA-seq analysis summary", "", f"Question: {question}",
             f"Run status: {state.get('status', 'unknown')}",
             "Executed stages: " + ", ".join(tasks), ""]
    plan = state.get("config", {}).get("question_plan", {})
    if plan.get("provider"):
        lines.append(f"Question planner: {plan['provider']}" +
                     (f" ({plan['model']})" if plan.get("model") else "") + ".")
    if plan.get("unsupported_requests"):
        lines.append("Unavailable requested analyses: " +
                     ", ".join(plan["unsupported_requests"]) + ".")
    if inspection:
        lines.append(f"Input: {inspection.get('cells', 'unknown')} cells and {inspection.get('genes', 'unknown')} genes; expression scale: {inspection.get('matrix_kind', 'unknown')}.")
    if qc:
        lines.append(f"QC retained {qc.get('retained_cells', 'unknown')} of {qc.get('input_cells', 'unknown')} cells.")
        qc_outputs = get("qc").get("outputs", {})
        if qc_outputs.get("thresholds"):
            lines.append(f"QC thresholds: {qc_outputs['thresholds']}")
        if qc_outputs.get("cell_qc"):
            lines.append(f"Per-cell QC results: {qc_outputs['cell_qc']}")
    selected = cluster.get("metrics", {}).get("selected", {})
    if selected:
        lines.append(f"Clustering: {selected.get('n_clusters', 'unknown')} clusters at Leiden resolution {selected.get('resolution', 'unknown')}.")
    annotation_file = cluster.get("outputs", {}).get("annotations")
    if annotation_file and Path(annotation_file).is_file():
        annotations = json.loads(Path(annotation_file).read_text(encoding="utf-8"))
        provisional = cluster.get("metrics", {}).get("annotation_backend") == "ollama_provisional"
        lines.append(("Provisional LLM marker hypotheses: " if provisional else "Marker-based cluster labels: ") + "; ".join(
            f"{item['cluster']} = {item['label']} (heuristic confidence {item['confidence']:.2f})"
            for item in annotations) + ".")
        if provisional:
            lines.append("No compatible annotation reference was used; review these identities before biological inference.")
        if all(item["label"] in {"unknown", "ambiguous"} for item in annotations):
            lines.append("No supported cell identities were assigned. Supply suitable multi-gene marker sets or a validated reference.")
    elif cluster.get("status") in {"failed", "blocked", "skipped"}:
        lines.append(f"Clustering was {cluster['status']}.")
    if "discovery" in tasks:
        if discovery.get("outputs", {}).get("differential_expression"):
            lines.append("A donor-aware pseudobulk differential-expression table was produced; review its model diagnostics and adjusted p-values before interpretation.")
        else:
            reason = ""
            skipped_file = discovery.get("outputs", {}).get("skipped_analyses")
            if skipped_file and Path(skipped_file).is_file():
                reason = json.loads(Path(skipped_file).read_text(encoding="utf-8")).get("differential_expression", "")
            lines.append("Condition-level differential expression was not established." + (f" Reason: {reason}" if reason else ""))
    failed = [name for name, result in tasks.items() if result.get("status") in {"failed", "blocked"}]
    if failed:
        lines.append("Failed or blocked stages: " + ", ".join(failed) + ".")
    report = get("report").get("outputs", {}).get("report")
    if report:
        lines.append(f"Full scientific report: {report}")
    html = cluster.get("outputs", {}).get("html_report")
    if html:
        lines.append(f"QC and clustering HTML report: {html}")
    dataset = cluster.get("outputs", {}).get("dataset")
    if dataset:
        lines.append(f"Annotated dataset: {dataset}")
    if report:
        lines.append("Interpretation is limited to recorded computational evidence; review warnings and validation checks in the full report.")
    else:
        lines.append("Only the requested stages and their dependencies ran; review their result.json files for warnings.")
    summary = "\n".join(lines) + "\n"
    (root / "question_summary.md").write_text(summary, encoding="utf-8")
    return summary
