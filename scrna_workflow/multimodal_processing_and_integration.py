"""Xenium and multimodal processing built around the existing core specialists."""
from pathlib import Path

from .graph import spatial_cell_graph, plot_spatial_cells, plot_spatial_gene_density


def _result(outputs=None, metrics=None, warnings=None, status="completed", inputs=None):
    return dict(status=status, input_references=inputs or [], outputs=outputs or {},
                metrics=metrics or {}, warnings=warnings or [], recommended_next_actions=[])


def process_xenium(ctx):
    """Prepare Xenium cell-by-gene counts and reuse core QC/PCA/cell graph steps."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse

    source = Path(ctx["config"]["input_path"]).expanduser().resolve()
    out = Path(ctx["root"]) / "multimodal analysis output"
    out.mkdir(parents=True, exist_ok=True)
    
    matrix_path = source / "cell_feature_matrix.h5"
    cells_path = source / "cells.csv.gz"
    if not matrix_path.exists() or not cells_path.exists():
        return _result(status="blocked", warnings=["Xenium input requires cell_feature_matrix.h5 and cells.csv.gz in input directory."])

    data = sc.read_10x_h5(matrix_path)
    cells = pd.read_csv(cells_path)
    cells["cell_id"] = cells["cell_id"].astype(str)
    data.obs_names = data.obs_names.astype(str)
    cells = cells.set_index("cell_id").reindex(data.obs_names)
    
    if cells.isna().all(axis=1).any():
        raise ValueError("Xenium cell metadata does not align with feature-matrix barcodes.")
        
    data.obs = cells
    data.obs["x_centroid"] = pd.to_numeric(data.obs["x_centroid"])
    data.obs["y_centroid"] = pd.to_numeric(data.obs["y_centroid"])
    data.obsm["spatial"] = data.obs[["x_centroid", "y_centroid"]].to_numpy()
    data.uns["workflow_matrix_kind"] = "counts"
    data.uns["workflow_has_counts"] = True
    data.layers["counts"] = sparse.csr_matrix(data.X, dtype=np.int32)
    data.var["source_gene_id"] = data.var_names.astype(str)
    
    dataset = out / "xenium_prepared.h5ad"
    data.write_h5ad(dataset)

    # Execute spatial graph & plotting
    spatial_res = spatial_cell_graph(ctx, dataset_path=dataset)
    plot_res = plot_spatial_cells(ctx, dataset_path=spatial_res["outputs"]["dataset"])
    
    return _result(
        outputs={
            "dataset": str(dataset),
            "spatial_graph": spatial_res["outputs"]["dataset"],
            "spatial_overview": plot_res["outputs"]["overview"]
        },
        metrics={"cells": data.n_obs, "genes": data.n_vars},
        inputs=[str(source)]
    )


def run(ctx):
    """Process Xenium when configured; return an explicit skip for other inputs."""
    source = Path(ctx["config"].get("input_path", ""))
    if source.is_dir() and (source / "cell_feature_matrix.h5").exists():
        return process_xenium(ctx)
    return _result(status="skipped", warnings=["No Xenium cell_feature_matrix.h5 directory was supplied."])
