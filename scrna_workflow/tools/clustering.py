"""Checkpointed clustering stage using focused scientific tools."""
import pandas as pd
from ._core_common import _result, _load_step, _json, _save
from .annotation import annotate_clusters
from .embedding import compute_umap, plot_embedding
from .markers import rank_cluster_markers
from .cluster_report import render_cluster_report
from .knowledge import load_knowledge_base


def _cluster(adata, resolutions, seed):
    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score, silhouette_score

    if not resolutions:
        raise ValueError("Provide at least one clustering resolution")
    candidates, diagnostics = [], []
    for resolution in resolutions[:5]:
        labels = []
        for random_state in (seed, seed + 1):
            sc.tl.leiden(adata, resolution=float(resolution), key_added="candidate", random_state=random_state)
            labels.append(adata.obs["candidate"].astype(str).to_numpy())
        n_clusters = len(set(labels[0]))
        silhouette = silhouette_score(adata.obsm["X_pca"], labels[0], sample_size=min(2000, adata.n_obs), random_state=seed) if 1 < n_clusters < adata.n_obs else -1
        diagnostics.append({"resolution": float(resolution), "n_clusters": n_clusters, "silhouette": float(silhouette), "seed_ARI": float(adjusted_rand_score(*labels))})
        candidates.append(labels[0])
    selected = max(range(len(diagnostics)), key=lambda i: diagnostics[i]["silhouette"])
    adata.obs["cluster"] = pd.Categorical(candidates[selected])
    del adata.obs["candidate"]
    return diagnostics, diagnostics[selected]


def clustering(ctx):
    """Select a Leiden resolution, rank markers, annotate, and plot embedding."""
    adata, out, source = _load_step(ctx, "graph", "clustering")
    config, seed, warnings = ctx["config"], ctx.get("seed", 0), []
    diagnostics, selected = _cluster(adata, config.get("resolutions", [0.3, 0.6, 1.0]), seed)
    pd.DataFrame(diagnostics).to_csv(out / "resolution_diagnostics.csv", index=False)
    markers = rank_cluster_markers(adata)
    markers.to_csv(out / "markers.csv", index=False)
    if not markers.empty:
        warnings.append("Marker p-values are exploratory cell-level rankings, not condition-level inference.")
    marker_sets = config.get("markers") or {}
    knowledge = {"path": config.get("marker_path"), "source": "user supplied marker sets"}
    if config.get("knowledge_base_path"):
        marker_sets, knowledge = load_knowledge_base(config["knowledge_base_path"])
        organism = str(config.get("organism") or "").lower()
        if knowledge.get("organism") and organism and organism not in {
                str(knowledge["organism"]).lower(), "homo sapiens" if knowledge["organism"] == "human" else ""}:
            raise ValueError("Knowledge-base organism conflicts with workflow organism")
    backend = config.get("annotation_backend", "markers")
    annotation_audit = []
    if backend == "ollama" and marker_sets:
        from ..agents.annotation import OllamaAnnotationAgent
        annotations, annotation_audit = OllamaAnnotationAgent(
            config.get("annotation_model", "qwen2.5:7b")
        ).annotate(adata, markers, marker_sets, config.get("gene_symbol_column"))
    elif backend == "markers" or not marker_sets:
        annotations = annotate_clusters(adata, markers, marker_sets,
                                        config.get("gene_symbol_column"))
    else:
        raise ValueError("annotation_backend must be markers or ollama")
    if not marker_sets:
        warnings.append("No compatible marker knowledge base supplied: biological annotations remain unknown.")
    warning = compute_umap(adata, seed)
    if warning:
        warnings.append(warning)
    plot_embedding(adata, out / "embedding.png")
    outputs = {
        "dataset": _save(adata, out / "annotated.h5ad"), "markers": str(out / "markers.csv"),
        "annotations": _json(out / "annotations.json", annotations), "embedding": str(out / "embedding.png"),
        "resolution_diagnostics": str(out / "resolution_diagnostics.csv"),
    }
    if annotation_audit:
        outputs["llm_annotation_audit"] = _json(out / "llm_annotation_audit.json", {
            "backend": "local_ollama", "model": config.get("annotation_model", "qwen2.5:7b"),
            "knowledge_base": knowledge, "clusters": annotation_audit,
        })
    for role in ("sample", "donor", "condition"):
        column = config.get(f"{role}_column")
        if column:
            path = out / f"cluster_by_{role}.csv"
            pd.crosstab(adata.obs["cluster"], adata.obs[column]).to_csv(path)
            outputs[f"cluster_by_{role}"] = str(path)
    outputs.update(render_cluster_report(
        adata, out, ctx["artifacts"].get("qc", {}).get("outputs", {}),
        markers, diagnostics, annotations, warnings,
    ))
    metrics = {"method": "leiden", "selected": selected, "cells": adata.n_obs,
               "annotation_backend": backend if marker_sets else "none",
               "knowledge_base": knowledge,
               "annotation_confidence_type": "heuristic marker overlap, not calibrated probability"}
    return _result(outputs, metrics, warnings, inputs=[source])
