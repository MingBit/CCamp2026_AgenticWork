"""Local HTML report for QC, PCA, clustering, and marker evidence."""
from html import escape
from pathlib import Path
import json

import numpy as np
import pandas as pd


def _table(frame):
    return frame.to_html(index=False, escape=True, classes="data", border=0) if not frame.empty else "<p>No results available.</p>"


def render_cluster_report(adata, out, qc_outputs, markers, diagnostics, annotations, warnings):
    """Write a report using only artifacts from this run; return its file paths."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out)
    qc_path = qc_outputs.get("cell_qc")
    qc = pd.read_csv(qc_path, index_col=0) if qc_path and Path(qc_path).exists() else pd.DataFrame()
    passed = qc["qc_pass"].astype(str).str.lower().eq("true") if "qc_pass" in qc else pd.Series(dtype=bool)
    qc_summary = pd.DataFrame([
        {"measure": "Input cells", "value": len(qc) if len(qc) else adata.n_obs},
        {"measure": "QC passed", "value": int(passed.sum()) if len(passed) else adata.n_obs},
        {"measure": "Cells after filtering", "value": adata.n_obs},
        {"measure": "Genes after filtering", "value": adata.n_vars},
        {"measure": "RPL/RPS genes removed", "value": adata.uns.get("workflow_gene_filter", {}).get("removed_genes", 0)},
    ])
    for column, name in (("total_counts", "Median library size"),
                         ("n_genes_by_counts", "Median detected genes"),
                         ("pct_counts_mt", "Median mitochondrial %")):
        if column in qc and len(qc):
            values = pd.to_numeric(qc[column], errors="coerce")
            qc_summary.loc[len(qc_summary)] = [name, f"{values.median():.2f}"]

    figures = {}
    qc_columns = [name for name in ("total_counts", "n_genes_by_counts", "pct_counts_mt") if name in qc]
    if qc_columns:
        fig, axes = plt.subplots(1, len(qc_columns), figsize=(4 * len(qc_columns), 3.2), squeeze=False)
        for ax, column in zip(axes.ravel(), qc_columns):
            values = pd.to_numeric(qc[column], errors="coerce").dropna()
            ax.hist(values, bins=40, color="#367a9c")
            ax.set(xlabel=column, ylabel="Input cells")
        fig.tight_layout()
        path = out / "cluster_qc.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        figures["qc_figure"] = str(path)

    variance = np.asarray(adata.uns.get("pca", {}).get("variance_ratio", []), dtype=float)
    n_pcs = adata.obsm["X_pca"].shape[1]
    pc_table = pd.DataFrame({
        "PC": np.arange(1, min(len(variance), n_pcs) + 1),
        "Variance explained (%)": np.round(100 * variance[:n_pcs], 2),
    })
    if len(variance):
        fig, ax = plt.subplots(figsize=(7, 3.2))
        ax.bar(np.arange(1, len(variance) + 1), 100 * variance, color="#496ca6")
        ax.set(xlabel="Principal component", ylabel="Variance explained (%)",
               title=f"Selected PCs: {n_pcs}")
        fig.tight_layout()
        path = out / "selected_pcs.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        figures["pc_figure"] = str(path)

    cluster_sizes = adata.obs["cluster"].astype(str).value_counts().sort_index()
    annotation_by_cluster = {row["cluster"]: row for row in annotations}
    cluster_table = pd.DataFrame([
        {"Cluster": cluster, "Cells": int(size),
         "Cell type": annotation_by_cluster[cluster]["label"],
         "Confidence": round(float(annotation_by_cluster[cluster]["confidence"]), 2),
         "Marker evidence": "; ".join(
             f"{name}: {', '.join(genes)}" for name, genes, _ in annotation_by_cluster[cluster]["evidence"] if genes
         ) or "none"}
        for cluster, size in cluster_sizes.items()
    ])
    marker_columns = [c for c in ("group", "names", "logfoldchanges", "pvals_adj") if c in markers]
    top_markers = (
        markers.sort_values(["group", "pvals_adj"], na_position="last")
        .groupby("group", observed=True, sort=False).head(10)[marker_columns].copy()
        if not markers.empty else pd.DataFrame(columns=marker_columns)
    )
    if "logfoldchanges" in top_markers:
        top_markers["logfoldchanges"] = top_markers["logfoldchanges"].round(3)
    if "pvals_adj" in top_markers:
        top_markers["pvals_adj"] = top_markers["pvals_adj"].map(lambda x: f"{x:.3g}" if pd.notna(x) else "")
    thresholds_path = qc_outputs.get("thresholds")
    thresholds = pd.DataFrame(json.loads(Path(thresholds_path).read_text())) if thresholds_path and Path(thresholds_path).exists() else pd.DataFrame()
    resolution = pd.DataFrame(diagnostics)
    warning_items = "".join(f"<li>{escape(str(w))}</li>" for w in warnings)
    qc_image = '<img src="cluster_qc.png" alt="QC distributions">' if "qc_figure" in figures else ""
    pc_image = '<img src="selected_pcs.png" alt="Variance explained by selected PCs">' if "pc_figure" in figures else ""
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Clustering report</title>
<style>body{{font:15px/1.5 system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#172536}}
h1,h2{{color:#163b5e}} .data{{border-collapse:collapse;width:100%;margin:1rem 0 2rem}}
.data th,.data td{{padding:.45rem .6rem;border-bottom:1px solid #dce4eb;text-align:left}}
.data th{{background:#eef4f8}}img{{max-width:100%;height:auto;margin:1rem 0}}
.note{{background:#f3f7fa;padding:.7rem 1rem;border-left:4px solid #367a9c}}</style></head>
<body><h1>Clustering and QC report</h1>
<p class="note">Counts and QC are observations. PCA, clusters, marker rankings and cell labels are computational inferences.
Marker p-values are exploratory cell-level rankings, not donor-aware condition tests.</p>
<h2>QC results</h2>{_table(qc_summary)}{qc_image}
<h3>Sample-aware QC thresholds</h3>{_table(thresholds)}
<h2>Selected principal components</h2><p>{n_pcs} PCs were used for the expression graph.</p>{pc_image}{_table(pc_table)}
<h2>Clusters and marker-based annotation</h2><p>Embedding: <a href="embedding.png">view annotated embedding</a>.</p>
<img src="embedding.png" alt="Cluster embedding">{_table(cluster_table)}
<h3>Resolution diagnostics</h3>{_table(resolution)}
<h2>Top marker genes by cluster</h2>{_table(top_markers)}
<p>Complete markers: <a href="markers.csv">markers.csv</a>. Labels need multiple supporting markers;
unsupported or tied labels remain unknown or ambiguous.</p>
<h2>Warnings and limits</h2><ul>{warning_items}</ul></body></html>"""
    path = out / "report.html"
    path.write_text(html, encoding="utf-8")
    return {"html_report": str(path), **figures}
