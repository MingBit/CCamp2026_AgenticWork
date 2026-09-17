"""Compute and plot expression embeddings."""

def plot_embedding(adata, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    key = "X_umap" if "X_umap" in adata.obsm else "X_pca"
    coordinates = adata.obsm[key]
    fig, ax = plt.subplots(figsize=(7, 5))
    for cluster in adata.obs["cluster"].cat.categories:
        mask = (adata.obs["cluster"] == cluster).to_numpy()
        label = adata.obs.loc[mask, "cell_type"].iloc[0]
        ax.scatter(coordinates[mask, 0], coordinates[mask, 1], s=4, alpha=0.7, label=f"{cluster}: {label}", rasterized=True)
    axis = "UMAP" if key == "X_umap" else "PC"
    ax.set(xlabel=f"{axis} 1", ylabel=f"{axis} 2")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=2)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def compute_umap(adata, seed):
    """Attempt UMAP; retain PCA if the embedding is unavailable."""
    import scanpy as sc
    try:
        sc.tl.umap(adata, random_state=seed)
        return None
    except (ValueError, TypeError, RuntimeError) as error:
        adata.obsm.pop("X_umap", None)
        return f"UMAP unavailable: {error}; PCA embedding retained."
