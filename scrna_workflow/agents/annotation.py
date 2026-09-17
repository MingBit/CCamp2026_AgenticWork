"""Grounded local-model review of cluster marker evidence."""
import json
import re

from ..llm_planner import _ollama_chat
from ..tools.annotation import annotate_clusters


class OllamaAnnotationAgent:
    """Review marker evidence with a local reference, or make provisional calls."""

    def __init__(self, model="qwen2.5:7b", chat=None):
        if model.lower().endswith(":cloud"):
            raise ValueError("Cloud Ollama models are disabled for cell annotation.")
        self.model = model
        self.chat = chat or _ollama_chat

    def annotate(self, adata, ranked_markers, marker_sets, symbol_column=None):
        baseline = annotate_clusters(adata, ranked_markers, marker_sets, symbol_column)
        adata.obs["annotation_source"] = "local_ollama_knowledge_review"
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
        response = self.chat(messages, self.model, format_schema=schema,
                             max_output_tokens=min(2400, max(700, 100 * len(evidence))))
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

    def annotate_from_markers(self, adata, ranked_markers, organism=None, tissue=None,
                              symbol_column=None):
        """Propose broad labels without a reference; verify every cited marker."""
        symbol_column = symbol_column or next(
            (name for name in ("gene_symbols", "gene_symbol", "gene_name")
             if name in adata.var), None
        )
        symbols = (adata.var[symbol_column].astype(str).to_numpy()
                   if symbol_column else adata.var_names.astype(str).to_numpy())
        feature_symbols = {str(feature).upper(): str(symbol).upper()
                           for feature, symbol in zip(adata.var_names, symbols)}
        evidence = []
        allowed = {}
        for cluster in adata.obs["cluster"].cat.categories:
            key = str(cluster)
            rows = ranked_markers.loc[
                (ranked_markers["group"].astype(str) == key)
                & (ranked_markers["logfoldchanges"] > 0.25)
            ].head(20)
            genes = list(dict.fromkeys(
                feature_symbols.get(str(feature).upper(), str(feature).upper())
                for feature in rows["names"]
            ))
            evidence.append({"cluster": key, "positive_top_markers": genes})
            allowed[key] = set(genes)
        schema = {
            "type": "object",
            "properties": {"annotations": {"type": "array", "items": {"type": "object",
                "properties": {"cluster": {"type": "string"}, "label": {"type": "string"},
                               "supporting_genes": {"type": "array", "items": {"type": "string"}}},
                "required": ["cluster", "label", "supporting_genes"],
                "additionalProperties": False}}},
            "required": ["annotations"], "additionalProperties": False,
        }
        messages = [
            {"role": "system", "content": (
                "You are reviewing scRNA-seq clusters without a reference database. "
                "Use only the observed positive top marker genes in the input. "
                "Propose the most likely broad cell identity when at least two distinctive "
                "genes support it. Use unknown only when the pattern genuinely does not "
                "support a broad identity; use ambiguous for conflicting lineages. "
                "Do not invent genes. "
                "Separate identity from activation, cycling, stress, exhaustion, "
                "senescence and disease state. Your labels are provisional hypotheses, "
                "not validated annotations. Return short labels and gene lists only, "
                "without explanations. Return JSON only."
            )},
            {"role": "user", "content": json.dumps({"organism": organism or "unknown",
                                                    "tissue": tissue or "unknown",
                                                    "clusters": evidence})},
        ]
        response = (self.chat(messages, self.model, format_schema=schema,
                              max_output_tokens=min(2400, max(700, 100 * len(evidence))))
                    if any(item["positive_top_markers"] for item in evidence)
                    else {"annotations": [
                        {"cluster": item["cluster"], "label": "unknown",
                         "supporting_genes": []}
                        for item in evidence
                    ]})
        predictions = response.get("annotations") if isinstance(response, dict) else None
        # A model can over-apply the caution above and return unknown for every
        # cluster. Retry that specific inconclusive outcome once with a clearer
        # request; retain unknown if the second pass is also inconclusive.
        if (isinstance(predictions, list) and len(predictions) == len(evidence)
                and predictions and all(str(p.get("label", "")).lower() == "unknown"
                                        for p in predictions)
                and any(len(item["positive_top_markers"]) >= 2 for item in evidence)):
            retry_messages = [messages[0], {"role": "user", "content": json.dumps({
                "organism": organism or "unknown", "tissue": tissue or "unknown",
                "clusters": evidence,
                "instruction": "The first pass called every cluster unknown. Reassess each "
                               "marker pattern and propose a provisional broad identity "
                               "where two distinctive observed genes support it. Keep "
                               "unknown for genuinely unresolvable clusters.",
            })}]
            retry = self.chat(retry_messages, self.model, format_schema=schema,
                              max_output_tokens=min(2400, max(700, 100 * len(evidence))))
            if isinstance(retry, dict) and isinstance(retry.get("annotations"), list):
                predictions = retry["annotations"]
        if not isinstance(predictions, list) or len(predictions) != len(evidence):
            raise ValueError("Ollama must return one provisional annotation per cluster")
        by_cluster = {item.get("cluster"): item for item in predictions}
        if set(by_cluster) != set(allowed) or len(by_cluster) != len(predictions):
            raise ValueError("Ollama returned duplicate or unknown cluster identifiers")
        adata.obs["cell_type"] = "unknown"
        adata.obs["annotation"] = "unknown"
        adata.obs["annotation_confidence"] = 0.0
        adata.obs["annotation_source"] = "provisional_llm_no_reference"
        annotations, audit = [], []
        for item in evidence:
            cluster = item["cluster"]
            proposal = by_cluster[cluster]
            label = str(proposal.get("label", "unknown")).strip()
            genes = proposal.get("supporting_genes", [])
            if not isinstance(genes, list) or any(not isinstance(gene, str) for gene in genes):
                raise ValueError("Ollama returned invalid supporting genes")
            matched = sorted({gene.upper() for gene in genes} & allowed[cluster])
            complex_state = bool(re.search(
                r"activat|exhaust|senesc|stress|cycling|proliferat|malignan|diseas",
                label, re.I,
            ))
            accepted = (label.lower() not in {"unknown", "ambiguous"}
                        and not complex_state and len(matched) >= 2)
            final_label = label if accepted else (
                "ambiguous" if label.lower() == "ambiguous" else "unknown"
            )
            confidence = min(0.45, 0.25 + 0.05 * len(matched)) if accepted else 0.0
            mask = adata.obs["cluster"].astype(str) == cluster
            adata.obs.loc[mask, "cell_type"] = final_label
            adata.obs.loc[mask, "annotation"] = final_label
            adata.obs.loc[mask, "annotation_confidence"] = confidence
            annotations.append({"cluster": cluster, "label": final_label,
                                "confidence": confidence,
                                "evidence": [("observed positive markers", matched, len(allowed[cluster]))]})
            audit.append({"cluster": cluster, "proposed_label": label,
                          "accepted_label": final_label, "accepted": accepted,
                          "supporting_markers": matched,
                          "model_reason": "No external reference was supplied; label is a marker-based hypothesis.",
                          "evidence_type": "provisional local-LLM marker hypothesis"})
        return annotations, audit
