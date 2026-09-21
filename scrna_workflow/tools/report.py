"""Scientific report, evidence ledger, and figures."""
from pathlib import Path
import json
from ._downstream_common import _get, _base, _result, _dataset, _json

def report(ctx):
    """Write an artifact-linked scientific report and figures from actual data only."""
    out, cfg, artifacts = _base(ctx, "report")
    source = _dataset(artifacts)
    lines = ["# Single-cell analysis report", "", "This report separates observed outputs, computational inferences, and untested hypotheses.", "", "## Scope", "",
             f"Organism: {cfg.get('organism') or 'unspecified'}. Biological question: {cfg.get('biological_question') or 'unspecified'}.", ""]
    figures = []
    if source:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        a = ad.read_h5ad(source)
        lines += [f"Processed dataset contains {a.n_obs:,} cells and {a.n_vars:,} features. Source artifact: `{source}`.", ""]
        qc_cols = [c for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt") if c in a.obs]
        if qc_cols:
            fig, axes = plt.subplots(1, len(qc_cols), figsize=(4 * len(qc_cols), 3.5), squeeze=False)
            for ax, col in zip(axes.ravel(), qc_cols):
                ax.hist(a.obs[col].dropna(), bins=40, color="#367a9c")
                ax.set(xlabel=col, ylabel="Cells")
            fig.tight_layout()
            f = out / "qc_distributions.png"
            fig.savefig(f, dpi=300)
            plt.close(fig)
            figures.append(str(f))
        embedding = next((key for key in ("X_umap", "X_pca") if key in a.obsm), None)
        if embedding and a.obsm[embedding].shape[1] >= 2:
            label = next((c for c in ("annotation", "cluster") if c in a.obs), None)
            coords = a.obsm[embedding]
            fig, ax = plt.subplots(figsize=(7, 5))
            if label:
                groups = a.obs[label].astype(str)
                for name in sorted(groups.unique()):
                    ix = np.asarray(groups == name)
                    ax.scatter(coords[ix, 0], coords[ix, 1], s=3, label=name, rasterized=True)
                ax.legend(loc="upper left", bbox_to_anchor=(1, 1), markerscale=3, fontsize=7)
            else:
                ax.scatter(coords[:, 0], coords[:, 1], s=3, rasterized=True)
            ax.set(xlabel=embedding + " 1", ylabel=embedding + " 2", title="Computational expression embedding")
            fig.tight_layout()
            f = out / "embedding.png"
            fig.savefig(f, dpi=300, bbox_inches="tight")
            plt.close(fig)
            figures.append(str(f))
            lines += ["Embedding distances represent a computational expression projection and do not establish physical proximity.", ""]
        abundance = artifacts.get("discovery", {}).get("outputs", {}).get("abundance")
        if abundance and Path(abundance).is_file():
            table = pd.read_csv(abundance, sep="\t", index_col=0)
            if len(table) and len(table.columns):
                fig, ax = plt.subplots(figsize=(max(6, min(14, len(table) * .4)), 4))
                table.plot.bar(stacked=True, ax=ax, width=.8)
                ax.set(ylabel="Fraction of retained cells", title="Descriptive sample-level population abundance")
                ax.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=7)
                fig.tight_layout()
                f = out / "sample_abundance.png"
                fig.savefig(f, dpi=300, bbox_inches="tight")
                plt.close(fig)
                figures.append(str(f))
        regulon_method = artifacts.get("regulon", {}).get("metrics", {}).get("method", "coexpression")
        scenicplus = regulon_method == "scenicplus"
        activities = artifacts.get("regulon", {}).get("outputs", {}).get("activity")
        if activities and Path(activities).is_file():
            values = pd.read_csv(activities, sep="\t", index_col=0)
            if not values.index.astype(str).equals(a.obs_names.astype(str)):
                raise ValueError("Regulon activity cell identifiers/order do not align with report dataset.")
            groups = a.obs["cluster"] if "cluster" in a.obs else pd.Series("all", index=a.obs_names)
            means = values.groupby(groups, observed=True).mean()
            chosen = means.var(axis=0).nlargest(min(20, means.shape[1])).index
            means = means[chosen]
            if means.size:
                fig, ax = plt.subplots(figsize=(max(5, len(chosen) * .35), max(3, len(means) * .3)))
                plot = ax.imshow(means, cmap="viridis", aspect="auto")
                ax.set_xticks(range(len(chosen)), chosen, rotation=90)
                ax.set_yticks(range(len(means)), means.index)
                ax.set(title="SCENIC+ eRegulon activity (direct, gene-based)" if scenicplus else "Candidate coexpression module activity", ylabel="Cluster")
                fig.colorbar(plot, ax=ax, label="Mean AUC (cells with ATAC and RNA QC)" if scenicplus else "Mean target log-normalized expression")
                fig.tight_layout()
                f = out / ("eregulon_activity_heatmap.png" if scenicplus else "candidate_module_heatmap.png")
                fig.savefig(f, dpi=300)
                plt.close(fig)
                figures.append(str(f))
        rssfile = artifacts.get("regulon", {}).get("outputs", {}).get("rss_0")
        if scenicplus and rssfile and Path(rssfile).is_file():
            from ..regulon_scenicplus import plot_rss_heatmap
            scores = pd.read_csv(rssfile, sep="\t", index_col=0)
            if scores.notna().any().any():
                figures.append(plot_rss_heatmap(scores, out / "eregulon_rss_heatmap.png", int(cfg.get("scenicplus_rss_top_n", 5)),
                                                "SCENIC+ eRegulon specificity (RSS, direct gene-based) per cell type"))
                lines += ["Regulon specificity scores (RSS) describe how concentrated eRegulon activity is in each cell type; they are descriptive, not a statistical test.", ""]
        edgesfile = artifacts.get("regulon", {}).get("outputs", {}).get("edges")
        if edgesfile and Path(edgesfile).is_file():
            edges = pd.read_csv(edgesfile, sep="\t")
            weight = "rho_tf2g" if scenicplus else "spearman_r"
            if len(edges) and weight in edges:
                if scenicplus:
                    edges = edges.sort_values("importance_tf2g", ascending=False).drop_duplicates(["tf", "target"])
                    tf = edges.groupby("tf")["target"].nunique().idxmax()
                    selected = edges[edges["tf"] == tf].nlargest(15, "importance_tf2g")
                else:
                    tf = edges.groupby("tf")["spearman_r"].mean().idxmax()
                    selected = edges[edges["tf"] == tf].nlargest(15, "spearman_r")
                fig, ax = plt.subplots(figsize=(7, 6))
                angle = np.linspace(0, 2 * np.pi, len(selected), endpoint=False)
                for theta, (_, row) in zip(angle, selected.iterrows()):
                    xx, yy = np.cos(theta), np.sin(theta)
                    ax.plot([0, xx], [0, yy], color="#8497ad", lw=1 + abs(row[weight]))
                    ax.scatter(xx, yy, s=250, color="#cbe3e8", zorder=2)
                    ax.text(xx * 1.13, yy * 1.13, str(row.target), ha="center", va="center", fontsize=8)
                ax.scatter([0], [0], s=500, color="#efb86b", zorder=3)
                ax.text(0, 0, str(tf), ha="center", va="center", fontsize=9)
                title = ("Selected SCENIC+ eRegulon targets (top 15 by TF-gene importance)\nMotif-supported association, not causal evidence"
                         if scenicplus else "Selected TF–target coexpression candidates\nNo motif or causal support")
                ax.set(title=title, xlim=(-1.4, 1.4), ylim=(-1.4, 1.4), aspect="equal")
                ax.axis("off")
                fig.tight_layout()
                f = out / ("eregulon_tf_network.png" if scenicplus else "candidate_tf_network.png")
                fig.savefig(f, dpi=300)
                plt.close(fig)
                figures.append(str(f))
        markerfile = artifacts.get("clustering", {}).get("outputs", {}).get("markers")
        if markerfile and Path(markerfile).is_file() and "cluster" in a.obs:
            markers = pd.read_csv(markerfile)
            if len(markers) and {"group", "names"}.issubset(markers):
                selected = markers.groupby("group", observed=True).head(3)["names"].astype(str)
                genes = list(dict.fromkeys(g for g in selected if g in a.var_names))[:30]
                if genes:
                    from scipy import sparse
                    matrix = a[:, genes].X
                    values = matrix.toarray() if sparse.issparse(matrix) else matrix
                    means = pd.DataFrame(values, index=a.obs_names, columns=genes).groupby(a.obs["cluster"], observed=True).mean()
                    fig, ax = plt.subplots(figsize=(max(6, len(genes) * .3), max(3, len(means) * .3)))
                    plot = ax.imshow(means, aspect="auto", cmap="magma")
                    ax.set_xticks(range(len(genes)), genes, rotation=90)
                    ax.set_yticks(range(len(means)), means.index)
                    ax.set(title="Exploratory cluster marker expression", ylabel="Cluster")
                    fig.colorbar(plot, ax=ax, label="Mean normalized expression")
                    fig.tight_layout()
                    f = out / "cluster_marker_heatmap.png"
                    fig.savefig(f, dpi=300)
                    plt.close(fig)
                    figures.append(str(f))
        defile = artifacts.get("discovery", {}).get("outputs", {}).get("differential_expression")
        if defile and Path(defile).is_file():
            de = pd.read_csv(defile, sep="\t")
            if len(de):
                fig, ax = plt.subplots(figsize=(6, 4))
                finite = np.isfinite(de["log2FoldChange"]) & np.isfinite(de["padj_global"])
                for population, rows in de.loc[finite].groupby("population"):
                    ax.scatter(rows["log2FoldChange"], -np.log10(np.maximum(rows["padj_global"], 1e-300)), s=7, alpha=.65, label=str(population))
                ax.axhline(-np.log10(.05), color="gray", linestyle="--", linewidth=.8)
                ax.set(xlabel="log2 fold change (test / reference)", ylabel="−log10 global BH adjusted p", title="Independent-donor pseudobulk contrasts")
                ax.legend(title="Population", fontsize=7)
                fig.tight_layout()
                f = out / "pseudobulk_de_summary.png"
                fig.savefig(f, dpi=300)
                plt.close(fig)
                figures.append(str(f))
    else:
        lines += ["No supplied dataset was available for biological analysis. No biological result has been fabricated.", ""]
    findings = []
    scenicplus = artifacts.get("regulon", {}).get("metrics", {}).get("method") == "scenicplus"
    regulon_finding = (("SCENIC+ eRegulons (TF-region-gene triplets)", "pycisTopic topics/DARs, cisTarget/DEM motif enrichment, GBM TF-gene and region-gene importance with correlation (SCENIC+)",
                        "Exploratory; motif and accessibility supported association in the same cells, no binding or causal support",
                        "Validate key TF-region-target links with perturbation, TF ChIP/CUT&RUN, or held-out donors.")
                       if scenicplus else
                       ("Candidate TF-target coexpression modules", "Spearman correlation and same-data bootstrap sign support",
                        "Exploratory; no motif, binding or causal support", "Validate candidate TFs with perturbation and compatible motif/binding evidence."))
    if source:
        findings.append({"type": "observation", "finding": f"Retained dataset: {a.n_obs} cells and {a.n_vars} features",
                         "artifact": source, "method": "AnnData dimensions after QC", "confidence": "High for recorded dimensions; depends on QC integrity",
                         "follow_up": "Review per-sample cell losses and sample metadata."})
    for task, artifact_key, finding, method, confidence, follow_up in (
        ("clustering", "dataset", "Expression-based cluster assignments", "PCA neighborhood graph and Leiden clustering", "Exploratory; review stability and marker evidence", "Validate cell identities with independent marker panels or a compatible reference."),
        ("regulon", "edges", *regulon_finding),
        ("discovery", "abundance", "Sample-level population proportions", "Retained-cell counts divided by each sample total", "Descriptive; capture and filtering can alter proportions", "Test in independently sampled biological replicates with the correct design."),
        ("discovery", "differential_expression", "Donor-aware expression contrasts", "PyDESeq2 independent-donor pseudobulk; BH adjusted p-values", "Conditional on specified design and model assumptions", "Replicate leading effects in held-out donors and orthogonal assays."),
    ):
        artifact = artifacts.get(task, {}).get("outputs", {}).get(artifact_key)
        if artifact:
            findings.append({"type": "observation" if artifact_key == "abundance" else "computational inference", "finding": finding,
                             "artifact": artifact, "method": method, "confidence": confidence, "follow_up": follow_up})
    import csv
    ledger = out / "evidence_ledger.tsv"
    with ledger.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["type", "finding", "artifact", "method", "confidence", "follow_up"], delimiter="\t")
        writer.writeheader()
        writer.writerows(findings)
    lines += ["## Principal results and evidence", "", "Full evidence ledger: [evidence_ledger.tsv](evidence_ledger.tsv).", ""]
    for finding in findings:
        lines += [f"- **{finding['finding']}** ({finding['type']}). Method: {finding['method']}. Confidence: {finding['confidence']}. Follow-up: {finding['follow_up']}"]
    lines += [""]
    for task, result in artifacts.items():
        lines += ["## " + task.replace("_", " ").title(), "", "Status: **" + result.get("status", "unknown") + "**.", ""]
        if result.get("metrics"):
            lines += ["Observed metrics:", "", "```json", json.dumps(result["metrics"], indent=2, default=str), "```", ""]
        for key, path in result.get("outputs", {}).items():
            lines += [f"- Artifact ({key}): `{path}`"]
        for warning in result.get("warnings", []):
            lines += ["- Limitation: " + str(warning)]
        for action in result.get("recommended_next_actions", []):
            lines += ["- Follow-up: " + str(action)]
        lines += [""]
    lines += ["## Evidence interpretation", "", "Dataset counts and QC summaries are observations. Clusters, annotations and regulatory modules (coexpression candidates or SCENIC+ eRegulons) are computational inferences. Regulatory mechanisms remain hypotheses requiring independent perturbation or binding evidence. Sample-level descriptive tables do not establish statistically significant condition effects.", "", "Confidence is limited by the checks and unresolved warnings recorded above; successful execution is not biological validation.", "", "## Figures", ""]
    for f in figures:
        lines += [f"![{Path(f).stem}]({Path(f).name})", ""]
    f = out / "scientific_report.md"
    f.write_text("\n".join(lines))
    return _result(inputs=[source] if source else [], outputs={"report": str(f), "evidence_ledger": str(ledger), **{f"figure_{i}": p for i, p in enumerate(figures)}}, metrics={"figures": len(figures)}, warnings=[] if source else ["Biological reporting blocked by missing input data."])
