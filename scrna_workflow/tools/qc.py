"""Sample-aware count QC and preservation checks."""
import numpy as np
import pandas as pd
from scipy import sparse
from ._core_common import _result, _load_step, _json, _save

def _qc_metrics(adata, config):
    counts = adata.layers["counts"]
    totals = np.asarray(counts.sum(axis=1)).ravel()
    adata.obs["total_counts"] = totals
    adata.obs["n_genes_by_counts"] = np.asarray((counts > 0).sum(axis=1)).ravel()
    prefixes = {"human": "MT-", "homo sapiens": "MT-", "mouse": "mt-", "mus musculus": "mt-"}
    prefix = config.get("mitochondrial_prefix") or prefixes.get(str(config.get("organism", "")).lower())
    mitochondrial = adata.var["source_gene_id"].str.startswith(prefix).to_numpy() if prefix else np.zeros(adata.n_vars, dtype=bool)
    adata.obs["pct_counts_mt"] = (
        100 * np.asarray(counts[:, mitochondrial].sum(axis=1)).ravel() / np.maximum(totals, 1)
        if mitochondrial.any() else np.nan
    )
    return bool(mitochondrial.any())


def _qc_thresholds(adata, sample_column):
    """Use within-sample median ± 3 scaled MAD; only flag the high count tail."""
    groups = adata.obs[sample_column].astype(str) if sample_column else pd.Series("pooled_unknown", index=adata.obs_names)
    keep = (adata.obs["total_counts"].to_numpy() > 0) & (adata.obs["n_genes_by_counts"].to_numpy() > 0)
    high = np.zeros(adata.n_obs, dtype=bool)
    thresholds = []
    for sample in groups.unique():
        indices = np.flatnonzero(groups.to_numpy() == sample)
        row = {"sample": sample, "n_cells": len(indices), "method": "median +/- 3 scaled MAD; count metrics log1p-transformed"}
        for column in ("total_counts", "n_genes_by_counts", "pct_counts_mt"):
            values = adata.obs[column].to_numpy()[indices]
            is_count = column != "pct_counts_mt"
            transformed = np.log1p(values) if is_count else values
            usable = len(indices) >= 20 and np.isfinite(transformed).all()
            median = np.median(transformed) if usable else 0
            mad = 1.4826 * np.median(np.abs(transformed - median)) if usable else 0
            low, upper = (median - 3 * mad, median + 3 * mad) if mad > 0 else (-np.inf, np.inf)
            if is_count:
                low, upper = max(0, np.expm1(low)), np.expm1(upper)
                keep[indices] &= values >= low
                high[indices] |= values > upper
                row[f"{column}_lower"] = float(low)
                row[f"{column}_upper_flag"] = float(upper) if np.isfinite(upper) else None
            else:
                if np.isfinite(upper):
                    keep[indices] &= values <= upper
                row["mitochondrial_upper"] = float(upper) if np.isfinite(upper) else None
        thresholds.append(row)
    adata.obs["qc_pass"] = keep
    adata.obs["high_library_complexity_flag"] = high
    return thresholds


def qc(ctx):
    """Calculate sample-aware QC, filter cells/empty genes, and verify counts."""
    adata, out, source = _load_step(ctx, "inspection", "qc")
    config, warnings, outputs = ctx["config"], [], {}
    has_counts = adata.uns["workflow_has_counts"]
    if has_counts:
        has_mito = _qc_metrics(adata, config)
        sample = config.get("sample_column") or config.get("donor_column")
        outputs["thresholds"] = _json(out / "thresholds.json", _qc_thresholds(adata, sample))
        if not sample:
            warnings.append("No sample identifier: pooled QC may mask sample effects.")
        if not has_mito:
            warnings.append("Mitochondrial fraction unavailable: provide organism-aware identifiers/prefix.")
        warnings.append("High counts flag potential doublets only. Dedicated doublet/ambient RNA inference requires suitable inputs and configuration.")
    else:
        adata.obs["qc_pass"] = True
        warnings.append("No raw counts: count-based QC, filtering and doublet detection skipped.")
    adata.obs.to_csv(out / "cell_qc.csv")
    filtered = adata[adata.obs["qc_pass"].to_numpy()].copy()
    if has_counts:
        expressed = np.asarray(filtered.layers["counts"].sum(axis=0)).ravel() > 0
        filtered = filtered[:, expressed].copy()
        expected = adata[filtered.obs_names, filtered.var_names].layers["counts"]
        if (sparse.csr_matrix(expected) != sparse.csr_matrix(filtered.layers["counts"])).nnz:
            raise ValueError("Raw count preservation check failed")
    if min(filtered.shape) < 3:
        raise ValueError("QC leaves fewer than three cells or genes; review thresholds")

    condition, batch = config.get("condition_column"), config.get("batch_column")
    if condition and batch and batch in filtered.obs:
        table = pd.crosstab(filtered.obs[batch], filtered.obs[condition])
        table.to_csv(out / "batch_condition.csv")
        outputs["batch_condition"] = str(out / "batch_condition.csv")
        if (table.gt(0).sum(axis=1) == 1).all():
            warnings.append("Batch is nested in condition: batch and condition may be confounded.")
    outputs.update(dataset=_save(filtered, out / "filtered.h5ad"), cell_qc=str(out / "cell_qc.csv"))
    metrics = {
        "input_cells": adata.n_obs, "retained_cells": filtered.n_obs,
        "retained_genes": filtered.n_vars, "unique_cell_ids": True,
        "retained_counts_exactly_preserved": bool(has_counts),
    }
    return _result(outputs, metrics, warnings, inputs=[source])
