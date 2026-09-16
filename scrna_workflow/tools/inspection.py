"""Input inspection, metadata alignment, and immutable source snapshots."""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from scipy import sparse
from ._core_common import _result, _dir, _json, _save, _digest, _countlike

def _snapshot(source, target):
    """Copy once and reject changes to an existing source snapshot."""
    target.parent.mkdir(parents=True, exist_ok=True)
    checksum = _digest(source)
    if target.exists() and _digest(target) != checksum:
        raise ValueError("Immutable source snapshot differs from input; start a new run")
    if not target.exists():
        shutil.copy2(source, target)
        target.chmod(0o444)
    return {"path": str(source), "snapshot": str(target), "sha256": checksum}


def _read_input(source):
    import scanpy as sc

    if source.is_dir():
        return sc.read_10x_mtx(source, var_names="gene_symbols", make_unique=False)
    if source.suffix.lower() == ".h5ad":
        return sc.read_h5ad(source)
    if source.suffix.lower() in {".h5", ".hdf5"}:
        return sc.read_10x_h5(source)
    raise ValueError("Supported inputs: .h5ad, 10x .h5, or a 10x Matrix Market directory")


def _merge_metadata(adata, path):
    metadata = pd.read_csv(
        path, sep="\t" if path.suffix in {".tsv", ".txt"} else ",", index_col=0
    )
    metadata.index = metadata.index.astype(str)
    if (not metadata.index.is_unique or len(metadata) != adata.n_obs
            or not adata.obs_names.isin(metadata.index).all()):
        raise ValueError("Metadata requires unique cell IDs and coverage of every expression barcode")
    metadata = metadata.loc[adata.obs_names]
    for column in metadata:
        if column in adata.obs and not adata.obs[column].astype(str).equals(metadata[column].astype(str)):
            raise ValueError(f"Conflicting embedded and supplied metadata column: {column}")
        adata.obs[column] = metadata[column]


def _set_expression_kind(adata, config, warnings):
    kind = config.get("matrix_kind", "auto")
    if kind not in {"auto", "counts", "normalized", "log1p"}:
        raise ValueError("matrix_kind must be auto, counts, normalized, or log1p")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all():
        raise ValueError("Expression contains nonfinite values")
    layer = config.get("counts_layer") or ("counts" if "counts" in adata.layers else None)
    if layer and layer not in adata.layers:
        raise ValueError(f"Configured counts layer is missing: {layer}")
    counts = adata.layers[layer] if layer else adata.X
    countlike = _countlike(counts)
    if (kind == "counts" or layer) and not countlike:
        raise ValueError("Declared raw counts contain negative or noninteger values")
    has_counts = bool(countlike and (kind in {"auto", "counts"} or layer))
    if has_counts:
        adata.layers["counts"] = counts.copy()
    if kind == "auto":
        kind = "counts" if _countlike(adata.X) else "log1p" if "log1p" in adata.uns else "unknown"
        warnings.append("Expression scale is inferred; integer nonnegative values do not prove raw-count provenance.")
    adata.uns["workflow_matrix_kind"] = kind
    adata.uns["workflow_has_counts"] = has_counts


def inspect_data(ctx):
    """Read and inventory data, align metadata, and preserve source/counts."""
    config = ctx["config"]
    source = config.get("input_path")
    if not source or not Path(source).exists():
        result = _result(status="blocked", warnings=["Provide an existing input_path: h5ad or 10x counts."])
        result["recommended_next_actions"] = ["Supply the dataset path and any barcode-indexed metadata."]
        return result

    source = Path(source).resolve()
    out = _dir(ctx, "inspection")
    snapshot = out / "source"
    files = sorted(p for p in source.rglob("*") if p.is_file()) if source.is_dir() else [source]
    manifest = [
        _snapshot(path, snapshot / (path.relative_to(source) if source.is_dir() else path.name))
        for path in files
    ]
    adata = _read_input(snapshot if source.is_dir() else snapshot / source.name)
    if min(adata.shape) < 3:
        raise ValueError("At least three cells and three genes are required")
    if not adata.obs_names.is_unique:
        raise ValueError("Duplicate cell identifiers require upstream resolution")

    warnings = []
    duplicates = int(adata.var_names.duplicated().sum())
    adata.var["source_gene_id"] = adata.var_names.astype(str)
    if duplicates:
        adata.var_names_make_unique()
        warnings.append("Duplicate feature names made unique; original labels retained in source_gene_id.")
    if config.get("metadata_path"):
        metadata = Path(config["metadata_path"]).resolve()
        copied = snapshot / ("external_metadata" + metadata.suffix)
        manifest.append(_snapshot(metadata, copied))
        _merge_metadata(adata, copied)
    _set_expression_kind(adata, config, warnings)

    inventory = {
        "cells": adata.n_obs, "genes": adata.n_vars, "duplicate_genes": duplicates,
        "matrix_kind": adata.uns["workflow_matrix_kind"],
        "counts_available": adata.uns["workflow_has_counts"],
        "sparse": sparse.issparse(adata.X), "metadata_columns": list(adata.obs),
        "metadata_missing": adata.obs.isna().sum().to_dict(),
        "existing_obsm": list(adata.obsm), "existing_layers": list(adata.layers),
        "existing_graphs": list(adata.obsp),
    }
    for role in ("sample", "donor", "condition"):
        column = config.get(f"{role}_column")
        if column and column not in adata.obs:
            raise ValueError(f"Configured {role} column absent: {column}")
        if column:
            inventory[f"{role}_counts"] = adata.obs[column].value_counts(dropna=False).to_dict()
    outputs = {
        "dataset": _save(adata, out / "inspected.h5ad"),
        "inventory": _json(out / "inventory.json", inventory),
        "source_manifest": _json(out / "source_manifest.json", manifest),
    }
    return _result(outputs, inventory, warnings, inputs=[str(source)])
