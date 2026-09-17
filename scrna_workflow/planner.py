"""Route a question to a bounded set of workflow stages.

Optional Ollama planning returns data, never executable code or shell commands.
The dependency scheduler validates and expands the chosen task targets.
"""
import json
import os
import re
from urllib import request

from .orchestrator import DEPS


def rule_targets(question):
    """Handle common intents, including explicit scope limits, without an LLM."""
    text = question.lower()
    found = []
    patterns = (
        ("qc", r"\b(?:qc|quality[ -]?control|filter(?:ing)? cells)\b"),
        ("clustering", r"\b(?:annotat\w*|cell[ -]?typ\w*|cluster\w*|marker genes?)\b"),
        ("graph", r"\b(?:neighbor(?:hood)? graph|knn graph|cell[ -]?cell graph)\b"),
        ("regulon", r"\b(?:regulon\w*|tf[ -]?target|transcription factor network)\b"),
        ("discovery", r"\b(?:differential expression|deg\b|pseudobulk|condition\w*|abundance)\b"),
    )
    for task, pattern in patterns:
        if re.search(pattern, text):
            found.append(task)
    if found:
        return found
    return ["report"]


def _ollama_targets(question, model):
    if not model:
        raise ValueError("Provide --llm-model or SCRNA_LLM_MODEL for Ollama planning.")
    if model.lower().endswith(":cloud"):
        raise ValueError("Cloud Ollama models are disabled for local data planning.")
    schema = {
        "type": "object",
        "properties": {"targets": {"type": "array", "items": {"type": "string", "enum": list(DEPS)},
                                   "minItems": 1, "maxItems": 4}},
        "required": ["targets"], "additionalProperties": False,
    }
    prompt = (
        "Select only the requested scRNA-seq workflow stages. Return JSON matching the schema. "
        "Use qc for quality control only; clustering for clustering or marker annotation; "
        "graph for neighborhood graphs; regulon for TF-target modules; discovery for "
        "sample abundance or differential expression; report for a full analysis. "
        "Dependencies will be added automatically. Respect words such as 'only' and 'just'. "
        "Do not generate shell commands, paths, biological claims, or parameters.\n"
        f"Question: {question}"
    )
    body = json.dumps({"model": model, "messages": [{"role": "user", "content": prompt}],
                       "stream": False, "format": schema, "options": {"temperature": 0}}).encode()
    outgoing = request.Request("http://127.0.0.1:11434/api/chat", data=body,
                               headers={"Content-Type": "application/json"}, method="POST")
    with request.urlopen(outgoing, timeout=30) as response:
        payload = json.load(response)
    parsed = json.loads(payload["message"]["content"])
    targets = parsed.get("targets")
    if not isinstance(targets, list) or not targets or len(targets) > 4 or any(t not in DEPS for t in targets):
        raise ValueError("Ollama returned workflow tasks outside the allowed schema.")
    return list(dict.fromkeys(targets))


def plan_question(question, planner="auto", model=None):
    """Choose task targets; auto uses Ollama only when a model is configured."""
    if planner not in {"auto", "rules", "ollama"}:
        raise ValueError("Planner must be auto, rules, or ollama.")
    model = model or os.environ.get("SCRNA_LLM_MODEL")
    rules = rule_targets(question)
    # Explicit QC-only wording is a hard scope limit even if an LLM disagrees.
    qc_only = bool(re.search(r"\b(?:only|just)\b.{0,35}\b(?:qc|quality[ -]?control)\b", question, re.I))
    if qc_only:
        return {"targets": ["qc"], "planner": "rules", "reason": "explicit QC-only request"}
    if planner == "rules" or (planner == "auto" and not model):
        return {"targets": rules, "planner": "rules"}
    try:
        targets = _ollama_targets(question, model)
        return {"targets": targets, "planner": "ollama", "model": model}
    except (OSError, TimeoutError, ValueError, KeyError, json.JSONDecodeError) as exc:
        if planner == "ollama":
            raise RuntimeError(f"Local Ollama planning failed: {exc}") from exc
        return {"targets": rules, "planner": "rules", "fallback_reason": str(exc)}
