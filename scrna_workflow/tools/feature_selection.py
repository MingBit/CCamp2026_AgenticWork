"""Select highly variable genes for dimensionality reduction."""

def select_hvg(adata, n_top_genes):
    """Return the number of selected genes and any fallback warning."""
    import scanpy as sc

    warning = None
    try:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=min(int(n_top_genes), adata.n_vars), flavor="seurat"
        )
    except (ValueError, IndexError):
        adata.var["highly_variable"] = True
        warning = "HVG selection failed on small/degenerate data; using all genes."
    if adata.var["highly_variable"].sum() < 3:
        adata.var["highly_variable"] = True
    return int(adata.var["highly_variable"].sum()), warning
