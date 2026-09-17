"""Grounded local-model review of cluster marker evidence."""
import json

from ..llm_planner import _ollama_chat
from ..tools.annotation import annotate_clusters


class OllamaAnnotationAgent:
    """Review only labels supported by a supplied marker knowledge base."""

    def __init__(self, model="qwen2.5:7b", chat=None):
        if model.lower().endswith(":cloud"):
            raise ValueError("Cloud Ollama models are disabled for cell annotation.")
        self.model = model
        self.chat = chat or _ollama_chat

    def annotate(self, adata, ranked_markers, marker_sets, symbol_column=None):
        baseline = annotate_clusters(adata, ranked_markers, marker_sets, symbol_column)
        supported = {
            row["cluster"]: {name: genes for name, genes, _ in row["evidence"] if len(genes) >= 2}
            for row in baseline
        }
        schema = {
            "type": "object",
            "properties": {"annotations": {"type": "array", "items": {"type": "object",
                "properties": {"cluster": {"type": "string"},
                               "label": {"type": "string", "enum": [*marker_sets, "unknown", "ambiguous"]},
                               "reason": {"type": "string"}},
                "required": ["cluster", "label", "reason"], "additionalProperties": False}}},
            "required": ["annotations"], "additionalProperties": False,
        }
        evidence = [{"cluster": row["cluster"], "candidate_markers": supported[row["cluster"]],
                     "overlap_label": row["label"]} for row in baseline]
        messages = [
            {"role": "system", "content": (
                "You review scRNA-seq cluster identities using only the supplied marker "
                "knowledge base and observed positive marker overlap. Return JSON only. "
                "Choose a cell type only when at least two listed markers support it. "
                "If evidence conflicts, choose ambiguous; if insufficient, unknown. "
                "Do not infer activation, exhaustion, senescence, spatial contact, or causality."
            )},
            {"role": "user", "content": json.dumps({"knowledge_base": marker_sets,
                                                    "cluster_evidence": evidence})},
        ]
        response = self.chat(messages, self.model, format_schema=schema)
        if not isinstance(response, dict) or not isinstance(response.get("annotations"), list):
            raise ValueError("Ollama annotation response is not a list of cluster labels")
        predictions = response["annotations"]
        if len(predictions) != len(baseline) or len({p.get("cluster") for p in predictions}) != len(baseline):
            raise ValueError("Ollama must return exactly one label per cluster")
        by_cluster = {p["cluster"]: p for p in predictions}
        audit = []
        for row in baseline:
            cluster = row["cluster"]
            if cluster not in by_cluster:
                raise ValueError(f"Ollama omitted cluster {cluster}")
            proposed = by_cluster[cluster]
            label = proposed.get("label")
            if label not in {*marker_sets, "unknown", "ambiguous"}:
                raise ValueError(f"Ollama proposed an unknown cell type: {label}")
            accepted = label in {"unknown", "ambiguous"} or label in supported[cluster]
            final_label = label if accepted else "unknown"
            matching = supported[cluster].get(final_label, [])
            confidence = (min(0.8, len(matching) / max(len(marker_sets[final_label]), 1))
                          if matching else 0.25 if final_label == "ambiguous" else 0.0)
            mask = adata.obs["cluster"].astype(str) == cluster
            adata.obs.loc[mask, "cell_type"] = final_label
            adata.obs.loc[mask, "annotation"] = final_label
            adata.obs.loc[mask, "annotation_confidence"] = confidence
            row["label"] = final_label
            row["confidence"] = confidence
            audit.append({"cluster": cluster, "proposed_label": label,
                          "accepted_label": final_label, "accepted": accepted,
                          "supporting_markers": matching,
                          "model_reason": str(proposed.get("reason", ""))[:500]})
        return baseline, audit
