"""Guarded independent-donor pseudobulk differential expression tool."""
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


