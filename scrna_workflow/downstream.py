"""File-backed downstream specialists. Regulatory edges are hypotheses, not causality.

Each public function accepts the orchestrator context and returns a JSON-serializable
specialist result. No private data leave the host. Output ownership is per task.
"""
from pathlib import Path
import json


def _get(ctx, key, default=None):
    return ctx.get(key, default) if isinstance(ctx, dict) else getattr(ctx, key, default)


def _base(ctx, name):
    directory = Path(_get(ctx, "root")) / name
    directory.mkdir(parents=True, exist_ok=True)
    return directory, _get(ctx, "config", {}), _get(ctx, "artifacts", {})


def _result(status="completed", inputs=None, outputs=None, metrics=None, warnings=None, actions=None):
    return dict(status=status, input_references=inputs or [], outputs=outputs or {},
                metrics=metrics or {}, warnings=warnings or [], recommended_next_actions=actions or [])


def _dataset(artifacts):
    for task in ("clustering", "qc", "inspection"):
        outputs = artifacts.get(task, {}).get("outputs", {})
        for key in ("dataset", "annotated_dataset", "filtered_dataset", "adata"):
            candidate = outputs.get(key)
            if candidate and Path(candidate).suffix == ".h5ad" and Path(candidate).is_file():
                return str(candidate)
    return None


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False, default=str))
    return str(path)


def regulon(ctx):
    """Infer regulatory candidates with the configured `regulon_method`.

    `coexpression` (default, RNA only): bounded, bootstrap-tested TF coexpression
    candidates from a local TF list. `scenicplus` (paired RNA + ATAC): SCENIC+ eRegulons
    run in the separate SCENIC+ environment (see regulon_scenicplus).
    """
    out, cfg, artifacts = _base(ctx, "regulon")
    source = _dataset(artifacts)
    method = cfg.get("regulon_method") or "coexpression"
    if method == "scenicplus":
        from . import regulon_scenicplus
        return regulon_scenicplus.run(ctx, out, cfg, source, _result)
    if method != "coexpression":
        message = f"Unknown regulon_method '{method}'; use 'coexpression' or 'scenicplus'."
        return _result("skipped", warnings=[message], actions=[message])
    return _coexpression_regulon(out, cfg, source, _get(ctx, "seed", 0))


def _coexpression_regulon(out, cfg, source, seed):
    """Bounded, bootstrap-tested TF coexpression candidates from local resources.

    Requires a local one-column TF-symbol list with matching declared organism and
    verified counts. Rank association is computed on log library-normalized counts.
    Activity is a target-expression score, not a binding or causal estimate.
    """
    tf_file = cfg.get("tf_list")
    prerequisites = []
    if not source:
        prerequisites.append("A validated upstream AnnData dataset is required.")
    if not tf_file or not Path(tf_file).is_file():
        prerequisites.append("A local, one-column TF list is required (tf_list).")
    if not cfg.get("organism") or cfg.get("tf_resource_organism") != cfg.get("organism"):
        prerequisites.append("Declare organism and matching tf_resource_organism.")
    if prerequisites:
        return _result("skipped", warnings=prerequisites, actions=prerequisites)
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse, stats
    a = ad.read_h5ad(source)
    if "counts" not in a.layers or not a.uns.get("workflow_has_counts", False):
        return _result("skipped", [source], warnings=["Verified raw counts are unavailable."], actions=["Provide verified counts for regulatory preprocessing."])
    if a.n_obs < 30:
        return _result("skipped", [source], warnings=["Fewer than 30 cells; candidate-module stability is not assessable."])
    tf_names = set(line.strip().split("\t")[0] for line in Path(tf_file).read_text().splitlines() if line.strip() and not line.startswith("#"))
    names = np.asarray(a.var_names.astype(str))
    tf_idx = np.flatnonzero(np.isin(names, list(tf_names)))
    if not len(tf_idx):
        return _result("skipped", [source, str(tf_file)], warnings=["No TF identifiers match dataset var_names; supply an explicit identifier mapping."])
    raw = sparse.csr_matrix(a.layers["counts"], dtype=float)
    totals = np.asarray(raw.sum(axis=1)).ravel()
    x = (sparse.diags(1e4 / np.maximum(totals, 1)) @ raw).tocsr()
    x.data = np.log1p(x.data)
    mean = np.asarray(x.mean(axis=0)).ravel()
    var = np.asarray(x.multiply(x).mean(axis=0)).ravel() - mean ** 2
    target_idx = np.argsort(var)[-min(int(cfg.get("regulon_max_genes", 2000)), a.n_vars):]
    rng = np.random.default_rng(seed)
    sample = np.sort(rng.choice(a.n_obs, min(a.n_obs, int(cfg.get("regulon_max_cells", 2000))), replace=False))
    target = x[sample][:, target_idx].toarray()
    ranks = stats.rankdata(target, axis=0)
    ranks -= ranks.mean(axis=0)
    ranks /= np.maximum(np.linalg.norm(ranks, axis=0), 1e-12)
    edges, activities = [], {}
    max_targets = min(50, max(1, int(cfg.get("regulon_max_targets", 25))))
    for idx in tf_idx[:min(200, int(cfg.get("regulon_max_tfs", 100)))]:
        tf = x[sample, idx].toarray().ravel()
        if np.std(tf) == 0:
            continue
        r = stats.rankdata(tf)
        r -= r.mean()
        r /= max(np.linalg.norm(r), 1e-12)
        corr = r @ ranks
        eligible = np.flatnonzero((corr >= float(cfg.get("regulon_min_correlation", .3))) & (target_idx != idx))
        chosen = eligible[np.argsort(corr[eligible])[-max_targets:]]
        if len(chosen) < 3:
            continue
        boot = []
        for _ in range(10):
            b = rng.choice(len(sample), len(sample), replace=True)
            boot.append([float(stats.spearmanr(tf[b], target[b, j]).statistic) for j in chosen])
        boot = np.asarray(boot)
        for k, j in enumerate(chosen):
            edges.append(dict(tf=names[idx], target=names[target_idx[j]], spearman_r=float(corr[j]),
                              bootstrap_positive_fraction=float(np.mean(boot[:, k] > 0)),
                              evidence_type="coexpression_only", motif_support=False,
                              confidence="exploratory; same-data bootstrap; not independent validation"))
        activities[names[idx]] = np.asarray(x[:, target_idx[chosen]].mean(axis=1)).ravel()
    if not edges:
        return _result("inconclusive", [source, str(tf_file)], warnings=["No TF has at least three targets passing the candidate correlation threshold."])
    edgefile = out / "tf_target_candidates.tsv"
    pd.DataFrame(edges).to_csv(edgefile, sep="\t", index=False)
    activity = pd.DataFrame(activities, index=a.obs_names)
    activity.index.name = "cell_id"
    afile = out / "candidate_module_activity.tsv.gz"
    activity.to_csv(afile, sep="\t")
    diagnostics = []
    for tf in activity:
        rho = stats.spearmanr(activity[tf], totals).statistic
        diagnostics.append({"tf": tf, "library_size_spearman": None if not np.isfinite(rho) else float(rho)})
    groups = [c for c in ["annotation", "cluster", cfg.get("sample_column"), cfg.get("donor_column"), cfg.get("condition_column")] if c and c in a.obs]
    summary_files = []
    for number, column in enumerate(dict.fromkeys(groups)):
        f = out / f"activity_group_{number}.tsv"
        activity.groupby(a.obs[column], observed=True).mean().to_csv(f, sep="\t")
        summary_files.append({"group_column": column, "path": str(f)})
    info = _json(out / "diagnostics.json", {"depth_associations": diagnostics, "summaries": summary_files,
                 "bootstrap_replicates": 10, "cell_cycle_adjusted": False, "batch_adjusted": False,
                 "sampling_unit": "cells; bootstrap is internal consistency only"})
    return _result(inputs=[source, str(tf_file)], outputs={"edges": str(edgefile), "activity": str(afile), "diagnostics": info,
                   **{f"summary_{i}": s["path"] for i, s in enumerate(summary_files)}},
                   metrics={"method": "coexpression", "candidate_tfs": len(activity.columns), "edges": len(edges), "cells_used_for_inference": len(sample)},
                   warnings=["Coexpression-only candidates: no motif or causal support; cell identity, batch, and cell cycle can confound associations.", "Activity uses the same data as inference and is not independent evidence of a cell state."])


def _differential_expression(pb, cfg, out):
    """Guarded independent-donor pseudobulk DE; paired designs require explicit work.

    Contrast is [reference, test]. BH is supplied per population by PyDESeq2 and
    also across every exported gene/population p-value. Filtering is count-only.
    """
    import numpy as np
    import pandas as pd
    from scipy import stats
    from scipy import sparse
    contrast = cfg.get("primary_comparison")
    condition, donor = cfg.get("condition_column"), cfg.get("donor_column")
    reason = None
    if not isinstance(contrast, (list, tuple)) or len(contrast) != 2 or contrast[0] == contrast[1]:
        reason = "Provide primary_comparison: [reference, test] for an explicit condition contrast."
    elif not condition or not donor or condition not in pb.obs or donor not in pb.obs:
        reason = "DE requires condition_column and donor_column with complete sample metadata."
    elif cfg.get("paired", False) or cfg.get("longitudinal", False):
        reason = "Paired/longitudinal DE is skipped: implement and review its subject-aware design explicitly."
    if reason:
        return {}, {"status": "skipped", "reason": reason}
    sample_meta = pb.obs.drop_duplicates("sample")
    sample_meta = sample_meta[sample_meta[condition].isin(contrast)]
    if sample_meta[[condition, donor]].isna().any().any():
        return {}, {"status": "skipped", "reason": "Missing donor or condition metadata."}
    if sample_meta[donor].duplicated().any():
        return {}, {"status": "skipped", "reason": "Donors contribute multiple samples or conditions; independent-sample design is invalid. Review pairing or aggregate technical replicates."}
    if any((sample_meta[condition] == level).sum() < 3 for level in contrast):
        return {}, {"status": "skipped", "reason": "At least three independent donors in each selected condition are required by this workflow."}
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        return {}, {"status": "skipped", "reason": "Install the optional pydeseq2 backend to run count-based pseudobulk DE."}
    batch = cfg.get("batch_column")
    if batch and batch not in pb.obs:
        return {}, {"status": "skipped", "reason": "Configured batch metadata unavailable or inconsistent within samples."}
    tables, details = [], []
    min_cells = max(1, int(cfg.get("de_min_cells_per_sample", 10)))
    min_count = max(1, int(cfg.get("de_min_total_count", 10)))
    for population in sorted(pb.obs["population"].unique()):
        subset = pb[(pb.obs["population"] == population) & pb.obs[condition].isin(contrast) & (pb.obs["n_cells"] >= min_cells)].copy()
        if any((subset.obs[condition] == level).sum() < 3 for level in contrast):
            details.append({"population": population, "status": "skipped", "reason": "Fewer than three donors per condition after population cell-count filtering."})
            continue
        counts = pd.DataFrame(subset.X.toarray() if sparse.issparse(subset.X) else subset.X, index=subset.obs_names, columns=subset.var_names)
        counts = counts.loc[:, counts.sum(axis=0) >= min_count]
        if counts.shape[1] < 10 or (counts.sum(axis=1) == 0).any():
            details.append({"population": population, "status": "skipped", "reason": "Too few expressed genes or an empty retained pseudobulk sample."})
            continue
        metadata = pd.DataFrame({"condition": subset.obs[condition].astype(str)}, index=subset.obs_names)
        design = "~condition"
        if batch:
            if subset.obs[batch].isna().any():
                details.append({"population": population, "status": "skipped", "reason": "Missing batch identifiers."})
                continue
            if subset.obs[batch].nunique() > 1:
                metadata["batch"] = subset.obs[batch].astype(str)
                design = "~batch + condition"
        model_matrix = pd.get_dummies(metadata, drop_first=True, dtype=float)
        model_matrix.insert(0, "intercept", 1.)
        if np.linalg.matrix_rank(model_matrix.to_numpy()) < model_matrix.shape[1] or len(model_matrix) - model_matrix.shape[1] < 2:
            details.append({"population": population, "status": "skipped", "reason": "Confounded or saturated design matrix."})
            continue
        try:
            import warnings as pywarnings
            with pywarnings.catch_warnings(record=True) as model_warnings:
                pywarnings.simplefilter("always")
                dds = DeseqDataSet(counts=counts.astype("int64"), metadata=metadata, design=design,
                                  n_cpus=1, quiet=True, size_factors_fit_type="poscounts")
                dds.deseq2()
                result = DeseqStats(dds, contrast=["condition", str(contrast[1]), str(contrast[0])], n_cpus=1, quiet=True)
                result.summary()
            table = result.results_df.copy().rename_axis("gene").reset_index()
            table.insert(0, "population", population)
            table["contrast"] = str(contrast[1]) + " vs " + str(contrast[0])
            tables.append(table)
            details.append({"population": population, "status": "completed", "donors": len(subset), "genes": len(table), "design": design,
                            "model_warnings": list(dict.fromkeys(str(w.message) for w in model_warnings))})
        except Exception as exc:
            details.append({"population": population, "status": "failed", "reason": type(exc).__name__ + ": " + str(exc), "retry": "No automatic retry; review model and data."})
    diagnostics = {"status": "completed" if tables else "inconclusive", "populations": details,
                   "min_cells_per_sample_population": min_cells, "min_total_gene_count": min_count,
                   "independent_units": "donors, one sample per donor", "size_factor_method": "poscounts",
                   "multiple_testing": "padj: within-population BH with independent filtering; padj_global: BH across exported finite p-values",
                   "log2FoldChange_direction": str(contrast[1]) + " / " + str(contrast[0])}
    outputs = {}
    if tables:
        table = pd.concat(tables, ignore_index=True)
        valid = np.isfinite(table["pvalue"])
        table["padj_global"] = np.nan
        table.loc[valid, "padj_global"] = stats.false_discovery_control(table.loc[valid, "pvalue"].to_numpy(), method="bh")
        file = out / "pseudobulk_differential_expression.tsv"
        table.to_csv(file, sep="\t", index=False)
        outputs["differential_expression"] = str(file)
    return outputs, diagnostics


def discovery(ctx):
    """Export sample-level descriptive abundance and count-preserving pseudobulk.

    Independent-donor DE runs only with an explicit contrast, adequate replication,
    and PyDESeq2. Paired or repeated-donor designs are deliberately deferred.
    """
    out, cfg, artifacts = _base(ctx, "discovery")
    source = _dataset(artifacts)
    if not source:
        return _result("skipped", warnings=["No processed dataset is available."])
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    a = ad.read_h5ad(source)
    sample = cfg.get("sample_column")
    condition = cfg.get("condition_column")
    donor = cfg.get("donor_column")
    label = "annotation" if "annotation" in a.obs and (a.obs["annotation"].astype(str) != "unknown").any() else "cluster"
    if label not in a.obs or not sample or sample not in a.obs:
        return _result("skipped", [source], warnings=["Sample-level discovery requires sample_column and clustering labels."], actions=["Supply validated biological sample identifiers; do not substitute individual cells."])
    if a.obs[sample].isna().any():
        return _result("skipped", [source], warnings=["Missing sample identifiers prevent defensible sample-level aggregation."])
    obs = a.obs.copy()
    warnings = []
    columns = list(dict.fromkeys(c for c in (condition, donor, cfg.get("batch_column")) if c and c in obs))
    metadata = obs.groupby(sample, observed=True)[columns].first() if columns else pd.DataFrame(index=pd.Index(obs[sample].unique(), name=sample))
    inconsistent = [c for c in columns if (obs.groupby(sample, observed=True)[c].nunique(dropna=False) > 1).any()]
    if inconsistent:
        warnings.append("Sample metadata are inconsistent within sample: " + ", ".join(inconsistent))
        metadata = metadata.drop(columns=inconsistent)
    table = pd.crosstab(obs[sample], obs[label])
    proportions = table.div(table.sum(axis=1), axis=0)
    abundance = out / "sample_abundance.tsv"
    proportions.to_csv(abundance, sep="\t")
    nfile = out / "sample_cell_counts.tsv"
    table.to_csv(nfile, sep="\t")
    mfile = out / "sample_metadata.tsv"
    metadata.to_csv(mfile, sep="\t")
    outputs = {"abundance": str(abundance), "cell_counts": str(nfile), "sample_metadata": str(mfile)}
    has_counts = "counts" in a.layers and bool(a.uns.get("workflow_has_counts", False))
    if has_counts:
        groups = obs.groupby([sample, label], observed=True, sort=True).indices
        rows, records = [], []
        for number, ((s, celltype), indices) in enumerate(groups.items()):
            rows.append(sparse.csr_matrix(a.layers["counts"])[indices].sum(axis=0))
            records.append({"pseudobulk_id": f"pb_{number}", "sample": str(s), "population": str(celltype), "n_cells": len(indices)})
        matrix = sparse.csr_matrix(np.vstack([np.asarray(r).ravel() for r in rows]))
        pb = ad.AnnData(matrix, obs=pd.DataFrame(records).set_index("pseudobulk_id"), var=a.var.copy())
        for c in metadata.columns:
            pb.obs[c] = pb.obs["sample"].map({str(k): v for k, v in metadata[c].items()})
        f = out / "pseudobulk_counts.h5ad"
        pb.write_h5ad(f)
        outputs["pseudobulk"] = str(f)
        de_outputs, de_info = _differential_expression(pb, cfg, out)
        outputs.update(de_outputs)
        outputs["de_diagnostics"] = _json(out / "de_diagnostics.json", de_info)
        warnings.extend("DE population " + str(p["population"]) + ": " + warning
                        for p in de_info.get("populations", []) for warning in p.get("model_warnings", []))
    else:
        warnings.append("Pseudobulk counts skipped: verified raw count layer unavailable.")
    prerequisites = {"differential_expression": "Reviewed primary contrast, independent biological replicates, pairing/covariates and a count-based DE backend are required. No cell-level condition tests were run.",
                     "condition_abundance_tests": "Confirm independent donor/sample units and paired/longitudinal design before inference.",
                     "pathway_comparison": "Compatible local gene-set resource and reviewed donor-aware design required.",
                     "ligand_receptor": "Organism-compatible local ligand-receptor resource and replicated design required.",
                     "trajectory": "A biological trajectory question and supported root/direction evidence are required."}
    if "differential_expression" in outputs:
        prerequisites.pop("differential_expression")
    elif has_counts:
        prerequisites["differential_expression"] = de_info.get("reason", "Review population-specific DE diagnostics.")
    outputs["skipped_analyses"] = _json(out / "skipped_analyses.json", prerequisites)
    warnings.append("Abundances are descriptive; no abundance p-values or causal mechanisms are inferred.")
    return _result(inputs=[source], outputs=outputs, metrics={"samples": len(table), "populations": len(table.columns), "population_column": label}, warnings=warnings, actions=list(prerequisites.values()))


def validation(ctx):
    """Check alignment and preservation; record limits rather than independent proof."""
    out, cfg, artifacts = _base(ctx, "validation")
    source = _dataset(artifacts)
    if not source:
        return _result("skipped", warnings=["No dataset is available for validation."])
    import anndata as ad
    import numpy as np
    from scipy import sparse
    a = ad.read_h5ad(source)
    checks = {"cell_identifiers_unique": bool(a.obs_names.is_unique), "gene_identifiers_unique": bool(a.var_names.is_unique),
              "metadata_aligned": len(a.obs) == a.n_obs and len(a.var) == a.n_vars}
    for key, matrix in a.obsp.items():
        checks["graph_dimensions_" + key] = matrix.shape == (a.n_obs, a.n_obs)
    if "counts" in a.layers and a.uns.get("workflow_has_counts", False):
        counts = a.layers["counts"]
        values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
        checks["count_shape"] = counts.shape == a.shape
        checks["counts_nonnegative_integer"] = bool(np.isfinite(values).all() and (values >= 0).all() and np.allclose(values, np.rint(values)))
    qcsource = artifacts.get("qc", {}).get("outputs", {}).get("dataset")
    if qcsource and Path(qcsource).is_file():
        q = ad.read_h5ad(qcsource)
        checks["qc_cell_order_preserved"] = bool(a.obs_names.equals(q.obs_names))
        checks["qc_gene_order_preserved"] = bool(a.var_names.equals(q.var_names))
        if "counts" in a.layers and "counts" in q.layers and a.shape == q.shape and checks["qc_cell_order_preserved"] and checks["qc_gene_order_preserved"]:
            delta = sparse.csr_matrix(a.layers["counts"]) - sparse.csr_matrix(q.layers["counts"])
            checks["counts_preserved_since_qc"] = delta.nnz == 0
    inspected = artifacts.get("inspection", {}).get("outputs", {}).get("dataset")
    if inspected and Path(inspected).is_file() and a.uns.get("workflow_has_counts", False) and "counts" in a.layers:
        original = ad.read_h5ad(inspected)
        identifiers_present = bool(a.obs_names.isin(original.obs_names).all() and a.var_names.isin(original.var_names).all())
        checks["retained_identifiers_in_inspected_input"] = identifiers_present
        if identifiers_present and "counts" in original.layers:
            retained = original[a.obs_names, a.var_names]
            delta = sparse.csr_matrix(a.layers["counts"]) - sparse.csr_matrix(retained.layers["counts"])
            checks["retained_counts_preserved_since_inspection"] = delta.nnz == 0
    issues = []
    for task, result in artifacts.items():
        for warning in result.get("warnings", []):
            issues.append({"owner": task, "issue": warning, "revision_cycles": 0, "resolution": "human review or targeted new evidence required"})
    concerns = ["These checks establish internal integrity, not independent biological validation.", "No external dataset or held-out donor validation was performed.", "Regulons defining a state cannot independently validate that state."]
    f = _json(out / "validation.json", {"checks": checks, "issues": issues, "limitations": concerns, "maximum_revision_cycles_per_issue": 2})
    return _result("completed" if all(checks.values()) else "failed", [source], {"validation": f}, {"checks": checks, "failed_checks": sum(not v for v in checks.values())}, concerns,
                   ["Review unresolved issues and use targeted donor-held-out or external validation where available."])


def report(ctx):
    """Write an artifact-linked scientific report and figures from actual data only."""
    out, cfg, artifacts = _base(ctx, "report")
    source = _dataset(artifacts)
    lines = ["# Single-cell analysis report", "", "This report separates observed outputs, computational inferences, and untested hypotheses.", "", "## Scope", "",
             f"Organism: {cfg.get('organism') or 'unspecified'}. Biological question: {cfg.get('biological_question') or 'unspecified'}.", ""]
    if cfg.get("organism") == "synthetic" or cfg.get("synthetic", False):
        lines += ["**SYNTHETIC TEST DATA: all outputs demonstrate software behavior only. They are not findings about biological samples or an organism.**", ""]
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
