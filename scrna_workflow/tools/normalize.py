"""Prepare the expression layer used for clustering without changing raw counts."""

def normalize_expression(adata):
    """Normalize counts, log normalized values, or accept declared log1p input."""
    import scanpy as sc

    kind = adata.uns["workflow_matrix_kind"]
    if kind == "counts":
        adata.X = adata.layers["counts"].copy()
        adata.uns.pop("log1p", None)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif kind == "normalized":
        sc.pp.log1p(adata)
    elif kind != "log1p":
        return False
    adata.layers["log_expression"] = adata.X.copy()
    return True
