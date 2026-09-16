"""PCA on selected genes."""
import numpy as np


def compute_pca(adata, n_pcs, seed):
    """Compute PCA and return the effective number of components."""
    import scanpy as sc

    selected = int(adata.var["highly_variable"].sum())
    components = min(int(n_pcs), adata.n_obs - 1, selected - 1)
    if components < 1:
        raise ValueError("PCA needs at least two cells and two selected genes")
    sc.pp.pca(adata, n_comps=components, mask_var="highly_variable", random_state=seed)
    if not np.isfinite(adata.obsm["X_pca"]).all():
        raise ValueError("PCA produced nonfinite values")
    return components
