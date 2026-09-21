"""Local LLM planning for natural-language analysis requests.

The model proposes structured configuration; it never executes shell commands or
supplies biological marker evidence. Only explicitly named local paths are used.
"""
import json
import os
from pathlib import Path
import re
from urllib.error import URLError
from urllib.request import Request, urlopen

from .question import prepare_question


FIELDS = ("input_path", "metadata_path", "marker_path", "knowledge_base_path", "output_dir",
          "sample_column", "donor_column", "condition_column", "organism", "tissue")
ANALYSES = ("qc", "clustering", "annotation", "graph", "regulon",
            "differential_expression", "discovery", "report")
TASK_FOR_ANALYSIS = {"qc": "qc", "clustering": "clustering",
                     "annotation": "clustering", "graph": "graph", "regulon": "regulon",
                     "differential_expression": "discovery", "discovery": "discovery",
                     "report": "report"}
SCHEMA = {
    "type": "object",
    "properties": {
        **{name: {"type": ["string", "null"]} for name in FIELDS},
        "primary_comparison": {"type": ["array", "null"], "items": {"type": "string"}},
        "requested_analyses": {"type": "array", "items": {"type": "string"}},
        "unsupported_requests": {"type": "array", "items": {"type": "string"}},
    },
    "required": [*FIELDS, "primary_comparison", "requested_analyses", "unsupported_requests"],
    "additionalProperties": False,
}


def _ollama_setting(name, default):
    """Read a positive integer override from the environment; keep the default otherwise.

    CPU-only hosts need a longer timeout than the 120 s default (SCRNA_OLLAMA_TIMEOUT), and
    Ollama otherwise starts one thread per host core, which oversubscribes a job that was
    allocated fewer cores (SCRNA_OLLAMA_NUM_THREADS).
    """
    value = os.environ.get(name)
    if value is None or not value.strip():
        return default
    try:
        number = int(value)
    except ValueError as exc:
        raise ValueError(f"{name} must be a positive whole number, got {value!r}") from exc
    if number <= 0:
        raise ValueError(f"{name} must be a positive whole number, got {value!r}")
    return number


def _ollama_chat(messages, model, timeout=None, format_schema=None, max_output_tokens=700):
    """Call only the loopback Ollama API; no dataset matrix is transmitted."""
    options = {"temperature": 0, "num_predict": max_output_tokens}
    threads = _ollama_setting("SCRNA_OLLAMA_NUM_THREADS", None)
    if threads:
        options["num_thread"] = threads
    timeout = timeout if timeout is not None else _ollama_setting("SCRNA_OLLAMA_TIMEOUT", 120)
    payload = json.dumps({"model": model, "messages": messages,
                          "format": format_schema or SCHEMA, "stream": False,
                          "options": options}).encode()
    request = Request("http://127.0.0.1:11434/api/chat", data=payload,
                      headers={"Content-Type": "application/json"})
    try:
        with urlopen(request, timeout=timeout) as response:
            result = json.load(response)
    except (OSError, URLError) as exc:
        raise RuntimeError(
            f"Local Ollama call failed after {timeout} s ({type(exc).__name__}). Start Ollama and install "
            "the selected model, raise SCRNA_OLLAMA_TIMEOUT for slow CPU-only hosts, choose a smaller "
            "model, or use --planner rules / annotation_backend: markers."
        ) from exc
    try:
        return json.loads(result["message"]["content"])
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError("The local model did not return a valid structured response.") from exc


def _explicit_path(value, question, label):
    """Reject paths that were invented by the model."""
    if not isinstance(value, str) or not value.strip():
        return None
    raw = value.strip().strip("\"'")
    if raw not in question:
        raise ValueError(f"The model proposed a {label} path absent from the question: {raw}")
    return str(Path(raw).expanduser().resolve())


def plan_question(question, base_config, model="qwen2.5:7b", chat=None):
    """Translate a question to validated workflow configuration and an audit plan."""
    if not question.strip():
        raise ValueError("Write a scRNA-seq analysis question after 'ask'.")
    chat = chat or _ollama_chat
    if model.lower().endswith(":cloud"):
        raise ValueError("Cloud Ollama models are disabled for local data planning.")
    context = {
        "question": question,
        "already_configured": {key: base_config.get(key) for key in
                               (*FIELDS, "primary_comparison") if base_config.get(key)},
        "available_analyses": ANALYSES,
    }
    messages = [
        {"role": "system", "content": (
            "Plan a local scRNA-seq analysis. Return only JSON matching the supplied schema. "
            "Extract paths and metadata column names only when stated in the user's question. "
            "Use null for missing values. Do not invent input data, markers, comparisons, "
            "results, shell commands, or biological annotations. The existing workflow "
            "will do QC, clustering, marker-based annotation, graph, candidate regulon, "
            "guarded donor-aware discovery, validation, and reporting. Identify requested "
            "analyses from the allowed list; list unsupported requests only when they "
            "are explicitly present in the user's words. Respect 'only' and 'just' scope limits. "
            "For a broad request to analyze the dataset, choose report."
        )},
        {"role": "user", "content": json.dumps(context)},
    ]
    proposed = chat(messages, model)
    if not isinstance(proposed, dict) or set(proposed) != set(SCHEMA["required"]):
        raise ValueError("The local model returned an incomplete analysis plan.")
    if not isinstance(proposed["requested_analyses"], list) or not isinstance(proposed["unsupported_requests"], list):
        raise ValueError("The local model returned invalid analysis lists.")
    if any(item not in ANALYSES for item in proposed["requested_analyses"]):
        raise ValueError("The local model requested an unavailable analysis tool.")
    if any(not isinstance(item, str) for item in proposed["unsupported_requests"]):
        raise ValueError("The local model returned invalid unsupported requests.")
    comparison = proposed["primary_comparison"]
    if comparison is not None and (not isinstance(comparison, list) or len(comparison) != 2
                                   or any(not isinstance(item, str) for item in comparison)
                                   or (not base_config.get("primary_comparison")
                                       and any(item not in question for item in comparison))):
        raise ValueError("The model proposed a comparison that was not stated in the question.")
    config = dict(base_config)
    accepted = {}
    path_fields = ("input_path", "metadata_path", "marker_path",
                   "knowledge_base_path", "output_dir")
    configured_paths = {
        str(Path(config[field]).expanduser().resolve())
        for field in path_fields if config.get(field)
    }
    for field in path_fields:
        if config.get(field):
            continue
        value = proposed[field]
        if isinstance(value, str) and value.strip():
            candidate = str(Path(value.strip().strip("\"'")).expanduser().resolve())
            if candidate in configured_paths:
                # The model may copy an existing knowledge-base path into marker_path.
                # Keep the user's configured role instead of treating it as a new path.
                continue
        path = _explicit_path(value, question, field)
        if path:
            config[field] = accepted[field] = path
    for field in ("sample_column", "donor_column", "condition_column", "organism", "tissue"):
        value = proposed[field]
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ValueError(f"The model returned an invalid {field}.")
        if value and not config.get(field):
            config[field] = accepted[field] = value.strip()
    if comparison and not config.get("primary_comparison"):
        config["primary_comparison"] = accepted["primary_comparison"] = comparison
    config = prepare_question(question, config)
    analyses = proposed["requested_analyses"]
    targets = list(dict.fromkeys(TASK_FOR_ANALYSIS[item] for item in analyses)) or ["report"]
    if re.search(r"\b(?:only|just)\b.{0,35}\b(?:qc|quality[ -]?control)\b", question, re.I):
        targets = ["qc"]
    if "report" in targets:
        targets = ["report"]  # The dependency closure includes every analysis stage.
    unsupported = [item for item in proposed["unsupported_requests"]
                   if item.lower().replace("_", " ") in question.lower()]
    config["requested_tasks"] = targets
    config["question_plan"] = {
        "provider": "local_ollama", "model": model,
        "requested_analyses": analyses, "targets": targets,
        "unsupported_requests": unsupported,
        "accepted_fields": accepted,
    }
    return config
