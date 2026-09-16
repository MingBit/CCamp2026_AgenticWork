"""Exploratory marker ranking by cluster."""
import pandas as pd


def rank_cluster_markers(adata):
    import scanpy as sc

    markers = pd.DataFrame(columns=["group", "names", "scores", "logfoldchanges", "pvals", "pvals_adj"])
    sizes = adata.obs["cluster"].value_counts()
    if len(sizes) > 1 and sizes.min() > 1:
        sc.tl.rank_genes_groups(adata, "cluster", method="wilcoxon", use_raw=False)
        markers = sc.get.rank_genes_groups_df(adata, group=None)
    return markers
