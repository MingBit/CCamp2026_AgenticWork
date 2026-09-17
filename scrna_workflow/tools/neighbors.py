"""Versioned expression-neighborhood graph tool."""
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import sparse
from ._core_common import _result, _load_step, _json, _save

def build_neighbors(ctx):
    """Build a standard Scanpy neighbor graph and export its edges/diagnostics."""
    import scanpy as sc
    from scipy.sparse.csgraph import connected_components

    adata, out, source = _load_step(ctx, "representation", "graph")
    n_neighbors = min(int(ctx["config"].get("n_neighbors", 15)), adata.n_obs - 1)
    if n_neighbors < 2:
        raise ValueError("n_neighbors must be at least two")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca", metric="euclidean", random_state=ctx.get("seed", 0))
    graph = adata.obsp["connectivities"].tocsr()
    if graph.shape != (adata.n_obs, adata.n_obs) or not np.isfinite(graph.data).all():
        raise ValueError("Invalid neighbor graph dimensions or weights")
    if (graph - graph.T).nnz or graph.diagonal().any():
        raise ValueError("Neighbor graph must be symmetric and have no self-edges")
    spec = json.dumps({"cells": adata.obs_names.tolist(), "genes": adata.var_names.tolist(), "params": adata.uns["neighbors"]["params"]}, sort_keys=True)
    version = hashlib.sha256(adata.obsm["X_pca"].tobytes() + spec.encode()).hexdigest()[:16]
    adata.uns["expression_graph_version"] = version
    edges = sparse.triu(graph).tocoo()
    pd.DataFrame({"source": adata.obs_names[edges.row], "target": adata.obs_names[edges.col], "weight": edges.data, "graph_version": version}).to_csv(out / "edges.csv", index=False)
    sparse.save_npz(out / "expression_graph.npz", graph)
    diagnostics = {
        "type": "expression_similarity", "representation": "X_pca", "metric": "euclidean",
        "method": "scanpy.pp.neighbors", "n_neighbors": n_neighbors, "version": version,
        "components": int(connected_components(graph)[0]), "cells": adata.n_obs, "edges": len(edges.data),
    }
    outputs = {
        "dataset": _save(adata, out / "graph.h5ad"), "graph": str(out / "expression_graph.npz"),
        "edges": str(out / "edges.csv"), "diagnostics": _json(out / "diagnostics.json", diagnostics),
    }
    return _result(outputs, diagnostics, ["Expression-similarity edges do not imply physical proximity."], inputs=[source])
