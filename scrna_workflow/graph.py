"""PCA-independent gene expression graph specialist."""
import hashlib
import json
from pathlib import Path
import os

import pandas as pd

from .core import _dir, _json, _previous, _result


def _multimodal_dir(ctx):
    root = Path(ctx["root"])
    out = root if root.name == "multimodal analysis output" else root / "multimodal analysis output"
    out.mkdir(parents=True, exist_ok=True)
    return out


def spatial_cell_graph(ctx, dataset_path=None, n_neighbors=8):
    """Build a spatial KNN graph from Xenium cell centroids."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from sklearn.neighbors import NearestNeighbors

    path = dataset_path or _previous(ctx, "graph")
    data = ad.read_h5ad(path)
    if "spatial" not in data.obsm:
        raise ValueError("Spatial graph requires obsm['spatial'] with x/y coordinates.")
    coordinates = np.asarray(data.obsm["spatial"], dtype=float)[:, :2]
    k = min(int(n_neighbors), data.n_obs - 1)
    if k < 1:
        raise ValueError("At least three cells are required for a spatial graph.")
    distances, indices = NearestNeighbors(n_neighbors=k + 1).fit(coordinates).kneighbors(coordinates)
    rows = np.repeat(np.arange(data.n_obs), k)
    columns = np.concatenate([ids[ids != row][:k] for row, ids in enumerate(indices)])
    values = np.concatenate([dist[ids != row][:k] for row, (ids, dist) in enumerate(zip(indices, distances))])
    scale = max(float(np.median(values[values > 0])) if (values > 0).any() else 1.0, 1e-12)
    weights = sparse.csr_matrix((np.exp(-values / scale), (rows, columns)), shape=(data.n_obs, data.n_obs))
    weights = weights.maximum(weights.T)
    weights.setdiag(0)
    weights.eliminate_zeros()
    upper = sparse.triu(weights).tocoo()
    out = _multimodal_dir(ctx)
    edges = pd.DataFrame({"source": data.obs_names[upper.row], "target": data.obs_names[upper.col], "weight": upper.data})
    edge_path = out / "spatial_cell_edges.csv"
    edges.to_csv(edge_path, index=False)
    data.obsp["spatial_connectivities"] = weights
    data.uns["spatial_graph"] = {"n_neighbors": k, "metric": "euclidean", "kernel": "exp(-distance / median_positive_distance)"}
    dataset = out / "spatial_graph.h5ad"
    data.write_h5ad(dataset)
    return _result(outputs={"dataset": str(dataset), "edges": str(edge_path)},
                   metrics={"nodes": data.n_obs, "edges": len(edges), "n_neighbors": k},
                   warnings=["Spatial edges indicate centroid proximity, not physical or molecular interaction."], inputs=[str(path)])


def plot_spatial_cells(ctx, dataset_path=None):
    """Plot whole-slide transcript density alongside a zoomed-in spatial KNN graph ROI."""
    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import sparse

    path = dataset_path or _previous(ctx, "spatial_graph")
    data = ad.read_h5ad(path)
    
    coords = np.asarray(data.obsm["spatial"], dtype=float)[:, :2]
    matrix = sparse.csr_matrix(data.obsp["spatial_connectivities"])
    edges = sparse.triu(matrix).tocoo()
    
    k_neighbors = data.uns.get("spatial_graph", {}).get("n_neighbors", 8)
    counts = np.log1p(np.asarray(data.X.sum(axis=1)).ravel())
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    
    # 1. Whole-slide transcript density
    sc1 = axes[0].scatter(
        coords[:, 0], coords[:, 1], 
        c=counts, s=1.0, cmap="viridis", linewidths=0, alpha=0.8
    )
    axes[0].set_title("Whole Slide Cell Transcript Density", fontsize=12, fontweight="bold")
    axes[0].set_xlabel("X Centroid Coordinate (µm)", fontsize=10)
    axes[0].set_ylabel("Y Centroid Coordinate (µm)", fontsize=10)
    axes[0].set_aspect("equal")
    axes[0].invert_yaxis()
    
    cb1 = fig.colorbar(sc1, ax=axes[0], fraction=0.046, pad=0.04)
    cb1.set_label("Total Expression [log1p(counts)]", fontsize=9)
    
    # 2. Zoomed-in ROI (700 µm x 700 µm patch at tissue center)
    x_mid, y_mid = np.median(coords[:, 0]), np.median(coords[:, 1])
    roi_radius = 350.0  # µm
    
    roi_mask = (
        (coords[:, 0] >= x_mid - roi_radius) & (coords[:, 0] <= x_mid + roi_radius) &
        (coords[:, 1] >= y_mid - roi_radius) & (coords[:, 1] <= y_mid + roi_radius)
    )
    roi_indices = set(np.where(roi_mask)[0])
    
    # Draw graph edges inside ROI
    for row, col in zip(edges.row, edges.col):
        if row in roi_indices and col in roi_indices:
            axes[1].plot(
                coords[[row, col], 0], coords[[row, col], 1],
                color="#e74c3c", alpha=0.7, linewidth=0.8, zorder=1
            )
            
    # Draw ROI cell nodes
    sc2 = axes[1].scatter(
        coords[roi_mask, 0], coords[roi_mask, 1],
        c=counts[roi_mask], s=20, cmap="viridis", edgecolors="black", linewidths=0.3, zorder=2
    )
    
    axes[1].set_title(f"Spatial {k_neighbors}-NN Graph Topology (700 µm ROI)", fontsize=12, fontweight="bold")
    axes[1].set_xlabel("X Centroid Coordinate (µm)", fontsize=10)
    axes[1].set_ylabel("Y Centroid Coordinate (µm)", fontsize=10)
    axes[1].set_aspect("equal")
    axes[1].set_xlim(x_mid - roi_radius, x_mid + roi_radius)
    axes[1].set_ylim(y_mid + roi_radius, y_mid - roi_radius)  # Match inverted Y orientation
    
    cb2 = fig.colorbar(sc2, ax=axes[1], fraction=0.046, pad=0.04)
    cb2.set_label("Total Expression [log1p(counts)]", fontsize=9)
    
    # Draw ROI indicator box on full slide
    rect = plt.Rectangle(
        (x_mid - roi_radius, y_mid - roi_radius), roi_radius * 2, roi_radius * 2,
        fill=False, edgecolor="#e74c3c", linewidth=1.5, linestyle="--"
    )
    axes[0].add_patch(rect)

    # Extract cell IDs inside ROI mask
    roi_nodes = data.obs_names[roi_mask].tolist()

    out_dir = _multimodal_dir(ctx)
    fig_out = out_dir / "spatial_cells_overview.png"

    # Define metadata columns you wish to display (e.g. total counts and cell area)
    # You can change this list to match any dataset format (e.g., ["cell_type", "total_counts"])

    edge_path = out_dir / "spatial_cell_edges.csv"
    if edge_path.exists():
        edges_df = pd.read_csv(edge_path)
        export_graph_metrics_md(
            edges_df,
            output_dir=str(out_dir),
            roi_nodes=roi_nodes,
            node_metadata=data.obs,
            label_cols=["total_counts", "cell_area"]  # data.obs.columns
        )

    fig.savefig(fig_out, dpi=300)
    plt.close(fig)
    return _result(outputs={"overview": str(fig_out)})


def plot_spatial_gene_density(ctx, dataset_path=None, gene="EPCAM"):
    """Plot gene expression sorted by intensity so expressing cells render on top."""
    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import sparse

    path = dataset_path or _previous(ctx, "spatial_graph")
    data = ad.read_h5ad(path)
    
    values = data[:, gene].X
    values = np.asarray(values.toarray()).ravel() if sparse.issparse(values) else np.asarray(values).ravel()
    coords = np.asarray(data.obsm["spatial"], dtype=float)[:, :2]
    
    # Sort cells so high-expressing cells are drawn ON TOP of zero-expressing cells
    sort_idx = np.argsort(values)
    coords_sorted = coords[sort_idx]
    values_sorted = np.log1p(values[sort_idx])
    
    fig, ax = plt.subplots(figsize=(8, 6.5), constrained_layout=True)
    
    # Pass 1D X and Y coordinate vectors explicitly
    sc = ax.scatter(
        coords_sorted[:, 0], 
        coords_sorted[:, 1], 
        c=values_sorted, 
        s=2.0, 
        cmap="magma", 
        linewidths=0, 
        alpha=0.85
    )
    
    ax.set_title(f"Spatial Expression: {gene}", fontsize=13, fontweight="bold")
    ax.set_xlabel("X Centroid Coordinate (µm)", fontsize=10)
    ax.set_ylabel("Y Centroid Coordinate (µm)", fontsize=10)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    
    cbar = fig.colorbar(sc, ax=ax, pad=0.03)
    cbar.set_label("Expression Level [log1p(counts)]", fontsize=10)
    
    out = _multimodal_dir(ctx) / f"spatial_gene_{gene}.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return _result(outputs={"overview": str(out)})


def gene_graph(ctx):
    """Build a KNN graph whose nodes are genes and profiles are cells."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components
    from sklearn.neighbors import NearestNeighbors
    from sklearn.decomposition import TruncatedSVD

    p = _previous(ctx, "qc")
    a = sc.read_h5ad(p)
    out = _dir(ctx, "graph") / "gene_graph"
    out.mkdir(parents=True, exist_ok=True)
    k = min(int(ctx["config"].get("n_neighbors", 15)), a.n_vars - 1)
    if k < 1:
        raise ValueError("n_neighbors must be positive")

    if a.uns.get("workflow_has_counts", False):
        x = sparse.csr_matrix(a.layers["counts"], dtype=float)
        totals = np.asarray(x.sum(axis=1)).ravel()
        x = sparse.diags(1e4 / np.maximum(totals, 1)) @ x
        x.data = np.log1p(x.data)
        expression_source = "counts normalized per cell to 10,000 and log1p transformed"
    else:
        x = sparse.csr_matrix(a.X, dtype=float)
        expression_source = "QC-filtered expression matrix without PCA"

    if x.shape[1] > 2000:
        variance = np.asarray(x.multiply(x).mean(axis=0)).ravel() - np.asarray(x.mean(axis=0)).ravel() ** 2
        feature_idx = np.argsort(variance)[-2000:]
        x = x[:, feature_idx]
        gene_names = a.var_names[feature_idx]
        feature_selection = "top 2,000 variance-ranked genes"
    else:
        gene_names = a.var_names
        feature_selection = "all available genes"

    # Transpose so each row is a gene profile across cells.
    x = x.T.tocsr()

    def build(n):
        ds, ix = NearestNeighbors(n_neighbors=n + 1, metric="euclidean").fit(x).kneighbors(x)
        clean = [(ids[ids != row][:n], dist[ids != row][:n]) for row, (ids, dist) in enumerate(zip(ix, ds))]
        rows = np.repeat(np.arange(x.shape[0]), n)
        cols = np.concatenate([value[0] for value in clean])
        distances = np.concatenate([value[1] for value in clean])
        scale = max(float(np.median(distances[distances > 0])) if (distances > 0).any() else 1.0, 1e-12)
        weights = sparse.csr_matrix((np.exp(-distances / scale), (rows, cols)), shape=(x.shape[0], x.shape[0]))
        weights = weights.maximum(weights.T)
        weights.setdiag(0)
        weights.eliminate_zeros()
        distance_matrix = sparse.csr_matrix((distances, (rows, cols)), shape=weights.shape)
        distance_matrix = distance_matrix.maximum(distance_matrix.T)
        return weights, distance_matrix

    weights, distances = build(k)
    if weights.shape != (x.shape[0], x.shape[0]) or (weights - weights.T).nnz or weights.diagonal().any():
        raise ValueError("Graph integrity check failed")

    spec = json.dumps(dict(k=k, metric="euclidean", representation="log_expression",
                           genes=gene_names.tolist(), cells=a.obs_names.tolist(),
                           feature_selection=feature_selection), sort_keys=True)
    version = hashlib.sha256(x.data.tobytes() + x.indices.tobytes() + spec.encode()).hexdigest()[:16]
    gene_data = ad.AnnData(x, obs=a.var.loc[gene_names].copy(), var=a.obs.copy())
    gene_data.obs_names = gene_names
    gene_data.var_names = a.obs_names
    gene_data.obsp["connectivities"] = weights
    gene_data.obsp["distances"] = distances
    gene_data.uns["neighbors"] = dict(connectivities_key="connectivities", distances_key="distances",
                               params=dict(n_neighbors=k, method="custom_exponential", metric="euclidean",
                                           use_rep="gene_expression_profile"))
    gene_data.uns["expression_graph_version"] = version

    diag = dict(type="expression_similarity", representation="log_expression", expression_source=expression_source,
                feature_selection=feature_selection, metric="euclidean", n_neighbors=k, version=version,
                kernel="exp(-distance / median_positive_neighbor_distance); scale=1 if all distances zero",
                symmetrization="maximum of directed neighbor weights; diagonal removed",
                nodes="genes", genes=x.shape[0], cells=a.n_obs,
                components=int(connected_components(weights)[0]),
                edges=int(weights.nnz // 2), sensitivity={})
    diag["metadata_mixing_descriptive"] = {}
    edge_rows, edge_cols = weights.nonzero()
    for neighbor_count in sorted(set([max(1, k // 2), min(x.shape[0] - 1, k * 2)])):
        alternate, _ = build(neighbor_count)
        union = (weights + alternate).astype(bool).nnz
        diag["sensitivity"][str(neighbor_count)] = dict(
            components=int(connected_components(alternate)[0]),
            edge_jaccard=float(weights.astype(bool).multiply(alternate.astype(bool)).nnz / max(union, 1)))

    upper = sparse.triu(weights).tocoo()
    edges = pd.DataFrame(dict(source=gene_names[upper.row], target=gene_names[upper.col],
                              weight=upper.data, graph_version=version))
    edges.to_csv(out / "edges.csv", index=False)
    degree = np.asarray(weights.getnnz(axis=1)).ravel()
    strength = np.asarray(weights.sum(axis=1)).ravel()
    diag["degree"] = dict(min=int(degree.min()), median=float(np.median(degree)), mean=float(degree.mean()),
                           max=int(degree.max()), isolated_genes=int(np.sum(degree == 0)))
    diag["weight"] = dict(min=float(upper.data.min()), median=float(np.median(upper.data)),
                           mean=float(upper.data.mean()), max=float(upper.data.max()))

    summary = out / "summary.md"
    summary.write_text("\n".join([
        "# Expression graph overview", "",
        f"- Genes: {x.shape[0]:,}", f"- Undirected edges: {len(edges):,}",
        f"- Neighbors per gene: {k}", "- Input: QC-filtered expression, independent of PCA",
        f"- Cells used as gene-profile dimensions: {a.n_obs:,}",
        f"- Connected components: {diag['components']}", f"- Genes used: {x.shape[0]:,} ({feature_selection})",
        f"- Median node degree: {np.median(degree):.1f}",
        f"- Edge weight range: {upper.data.min():.3f} to {upper.data.max():.3f}", "",
        "## Interpretation", "",
        "Edges connect genes with similar expression profiles across cells. Edge weights are relative similarity "
        "scores, not probabilities or physical interactions. The overview image uses a graph-derived "
        "layout and a deterministic sample of edges so dense graphs remain readable.", "",
        "## Files", "", "- `graph.h5ad`: gene-oriented AnnData dataset with graph matrices",
        "- `expression_graph.npz`: sparse weighted adjacency matrix",
        "- `edges.csv`: one row per undirected edge", "- `diagnostics.json`: machine-readable metrics",
        "- `overview.png`: visual graph summary", ""
    ]))

    rng = np.random.default_rng(ctx.get("seed", 0))
    edge_limit = min(len(edges), 5000)
    chosen = rng.choice(len(edges), edge_limit, replace=False) if edge_limit else np.array([], dtype=int)
    layout = TruncatedSVD(n_components=2, random_state=ctx.get("seed", 0)).fit_transform(x)
    layout -= layout.mean(axis=0, keepdims=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    axes[0].scatter(layout[:, 0], layout[:, 1], c=degree, s=7, cmap="viridis", alpha=.75, linewidths=0)
    axes[0].set_title("Gene graph layout colored by degree")
    axes[0].set_xlabel("SVD 1 of log-expression gene profiles")
    axes[0].set_ylabel("SVD 2 of log-expression gene profiles")
    degree_plot = axes[0].collections[-1]
    fig.colorbar(degree_plot, ax=axes[0], fraction=.046, pad=.04, label="Number of gene neighbors")
    for row in edges.iloc[chosen].itertuples():
        source = gene_data.obs_names.get_loc(row.source)
        target = gene_data.obs_names.get_loc(row.target)
        axes[1].plot(layout[[source, target], 0], layout[[source, target], 1],
                     color="0.55", alpha=.08, linewidth=.4)
    axes[1].scatter(layout[:, 0], layout[:, 1], c=strength, s=7, cmap="plasma", alpha=.8, linewidths=0)
    axes[1].set_title(f"Sampled gene edges ({edge_limit:,} of {len(edges):,})")
    axes[1].set_xlabel("SVD 1 of log-expression gene profiles")
    axes[1].set_ylabel("SVD 2 of log-expression gene profiles")
    strength_plot = axes[1].collections[-1]
    fig.colorbar(strength_plot, ax=axes[1], fraction=.046, pad=.04, label="Sum of edge similarities")
    fig.suptitle("Gene expression similarity graph", fontsize=14)
    fig.savefig(out / "overview.png", dpi=180)
    plt.close(fig)

    sparse.save_npz(out / "expression_graph.npz", weights)
    file = out / "graph.h5ad"
    gene_data.write_h5ad(file)
    return _result(dict(dataset=str(file), graph=str(out / "expression_graph.npz"),
                        edges=str(out / "edges.csv"), diagnostics=_json(out / "diagnostics.json", diag),
                        summary=str(summary), overview=str(out / "overview.png")), diag,
                   ["Gene expression edges indicate similarity, not regulatory, physical, or causal interactions."],
                   inputs=[p])


def graph(ctx):
    """Backward-compatible name for the normalized-expression gene graph."""
    return gene_graph(ctx)


def interpret_pca(ctx, pca_dataset=None, gene_graph_dataset=None, top_n=20, top_pcs=6):
    """Interpret PCA loadings using neighborhoods in a gene-level graph.

    The PCA dataset must contain ``varm['PCs']`` and the gene graph dataset must
    contain gene nodes in ``obs_names`` and ``obsp['connectivities']``. Only
    genes present in both artifacts are compared.
    """
    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy import sparse

    def artifact_path(task, explicit):
        if explicit:
            return str(explicit)
        item = ctx.get("artifacts", {}).get(task, {})
        outputs = item.get("outputs", item)
        return outputs.get("dataset")

    pca_path = artifact_path("representation", pca_dataset)
    gene_graph_path = gene_graph_dataset or artifact_path("gene_graph", None) or artifact_path("graph", None)
    if not pca_path or not gene_graph_path:
        raise ValueError("Provide PCA and gene-graph datasets, or matching representation/graph artifacts in ctx.")

    pca = ad.read_h5ad(pca_path)
    gene_graph = ad.read_h5ad(gene_graph_path)
    if "PCs" not in pca.varm or "X_pca" not in pca.obsm:
        raise ValueError("PCA dataset must contain varm['PCs'] and obsm['X_pca'].")
    if "connectivities" not in gene_graph.obsp:
        raise ValueError("Gene graph dataset must contain obsp['connectivities'].")

    loading_genes = pd.Index(pca.var_names.astype(str))
    graph_genes = pd.Index(gene_graph.obs_names.astype(str))
    shared = loading_genes.intersection(graph_genes)
    if len(shared) < 3:
        raise ValueError("PCA and gene graph share fewer than three genes.")
    top_n = max(1, int(top_n))
    loading_index = loading_genes.get_indexer(shared)
    graph_index = graph_genes.get_indexer(shared)
    loadings = np.asarray(pca.varm["PCs"])[loading_index]
    adjacency = sparse.csr_matrix(gene_graph.obsp["connectivities"])[graph_index][:, graph_index]
    degree = np.asarray(adjacency.getnnz(axis=1)).ravel()
    graph_rows, graph_cols = adjacency.nonzero()
    variance_ratio = np.asarray(pca.uns.get("pca", {}).get("variance_ratio", []), dtype=float)
    total_pcs = loadings.shape[1]
    n_pcs = min(max(1, int(top_pcs)), total_pcs)
    loadings = loadings[:, :n_pcs]
    variance_ratio = variance_ratio[:n_pcs]
    records = []
    pc_metrics = []
    for pc_index in range(n_pcs):
        values = loadings[:, pc_index]
        order = np.argsort(np.abs(values))[::-1]
        positive = [i for i in order if values[i] > 0][:top_n]
        negative = [i for i in order if values[i] < 0][:top_n]
        for side, indices in (("positive", positive), ("negative", negative)):
            for rank, gene_index in enumerate(indices, 1):
                neighbors = adjacency.getrow(gene_index).indices
                neighbor_values = values[neighbors] if len(neighbors) else np.array([])
                records.append({
                    "pc": pc_index + 1,
                    "explained_variance_ratio": float(variance_ratio[pc_index]) if pc_index < len(variance_ratio) else None,
                    "side": side,
                    "rank": rank,
                    "gene": str(shared[gene_index]),
                    "loading": float(values[gene_index]),
                    "graph_degree": int(degree[gene_index]),
                    "neighbor_mean_loading": float(np.mean(neighbor_values)) if len(neighbor_values) else None,
                    "neighbor_same_sign_fraction": float(np.mean(np.sign(neighbor_values) == np.sign(values[gene_index]))) if len(neighbor_values) else None,
                })
        edge_products = values[graph_rows] * values[graph_cols]
        edge_scales = np.abs(values[graph_rows]) * np.abs(values[graph_cols])
        pc_metrics.append({
            "pc": pc_index + 1,
            "explained_variance_ratio": float(variance_ratio[pc_index]) if pc_index < len(variance_ratio) else None,
            "shared_genes": len(shared),
            "positive_genes": int(np.sum(values > 0)),
            "negative_genes": int(np.sum(values < 0)),
            "graph_edges_within_shared_genes": int(len(graph_rows) // 2),
            "edge_same_sign_fraction": float(np.mean(edge_products >= 0)) if len(edge_products) else None,
            "weighted_edge_loading_agreement": float(np.average(np.sign(edge_products), weights=edge_scales)) if len(edge_products) and edge_scales.sum() else None,
        })

    out = _dir(ctx, "graph") / "pca_interpretation"
    out.mkdir(parents=True, exist_ok=True)
    table = out / "pca_gene_loadings.csv"
    pd.DataFrame(records).to_csv(table, index=False)
    metrics_path = out / "pca_graph_metrics.json"
    _json(metrics_path, {"pca_dataset": str(pca_path), "gene_graph_dataset": str(gene_graph_path),
                         "shared_genes": len(shared), "total_pcs": total_pcs,
                         "interpreted_pcs": n_pcs, "top_n_per_side": top_n, "pcs": pc_metrics})
    metrics_frame = pd.DataFrame(pc_metrics)
    overview = out / "overview.png"
    leading_pcs = min(6, n_pcs)
    figure = plt.figure(figsize=(15, 10), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, width_ratios=(1, 1.35), height_ratios=(1, 1.25))

    variance_axis = figure.add_subplot(grid[0, 0])
    variance_values = variance_ratio[:n_pcs] if len(variance_ratio) else np.zeros(n_pcs)
    variance_axis.bar(np.arange(1, n_pcs + 1), variance_values, color="#2a6f97", width=.8)
    variance_axis.set_title("Variance explained by each PC")
    variance_axis.set_xlabel("Principal component")
    variance_axis.set_ylabel("Variance ratio")
    variance_axis.set_xlim(.4, n_pcs + .6)

    coherence_axis = figure.add_subplot(grid[0, 1])
    coherence = metrics_frame["edge_same_sign_fraction"].fillna(0).to_numpy()
    coherence_axis.bar(np.arange(1, n_pcs + 1), coherence, color="#d95f02", width=.8)
    coherence_axis.axhline(.5, color="0.45", linestyle="--", linewidth=1)
    coherence_axis.set_title("Graph coherence of PC loadings")
    coherence_axis.set_xlabel("Principal component")
    coherence_axis.set_ylabel("Same-sign edge fraction")
    coherence_axis.set_ylim(0, 1)
    coherence_axis.set_xlim(.4, n_pcs + .6)

    signed_axis = figure.add_subplot(grid[1, 0])
    first_values = loadings[:, 0]
    first_order = np.argsort(np.abs(first_values))[-min(top_n, len(first_values)):][::-1]
    labels = [str(shared[i]) for i in first_order][::-1]
    signed_values = first_values[first_order][::-1]
    colors = ["#b2182b" if value < 0 else "#2166ac" for value in signed_values]
    signed_axis.barh(np.arange(len(labels)), signed_values, color=colors)
    signed_axis.set_yticks(np.arange(len(labels)), labels, fontsize=8)
    signed_axis.set_ylabel("Gene")
    signed_axis.set_title("Largest absolute loadings: PC1")
    signed_axis.set_xlabel("Signed loading")
    signed_axis.axvline(0, color="0.25", linewidth=.8)

    heatmap_axis = figure.add_subplot(grid[1, 1])
    heatmap_genes = []
    for pc_index in range(leading_pcs):
        order = np.argsort(np.abs(loadings[:, pc_index]))[::-1]
        heatmap_genes.extend(order[:min(5, len(order))].tolist())
    heatmap_genes = list(dict.fromkeys(heatmap_genes))
    heatmap_values = loadings[heatmap_genes, :leading_pcs]
    image = heatmap_axis.imshow(heatmap_values, aspect="auto", cmap="RdBu_r", vmin=-np.max(np.abs(heatmap_values)),
                                vmax=np.max(np.abs(heatmap_values)))
    heatmap_axis.set_yticks(np.arange(len(heatmap_genes)), [str(shared[i]) for i in heatmap_genes], fontsize=8)
    heatmap_axis.set_xticks(np.arange(leading_pcs), [f"PC{i}" for i in range(1, leading_pcs + 1)])
    heatmap_axis.set_title("Top loading genes across leading PCs")
    heatmap_axis.set_ylabel("Gene")
    heatmap_axis.set_xlabel("Blue = negative, red = positive")
    figure.colorbar(image, ax=heatmap_axis, fraction=.046, pad=.04, label="Loading")
    figure.suptitle("PCA interpretation through the gene expression graph", fontsize=16)
    figure.savefig(overview, dpi=180)
    plt.close(figure)
    summary = out / "summary.md"
    lines = ["# PCA interpretation using the gene graph", "",
             f"- Shared genes: {len(shared):,}", f"- Interpreted leading PCs: {n_pcs} of {total_pcs}",
             f"- Genes reported per PC and sign: {top_n}", "",
             "For each PC, the CSV lists genes with the largest positive and negative loadings. "
             "Graph coherence measures whether connected genes tend to have loadings with the same sign. "
             "This supports interpretation of coordinated expression programs; it is not evidence of regulation or causality.", "",
             "## Most graph-coherent components", ""]
    ranked_metrics = sorted(pc_metrics, key=lambda item: item["edge_same_sign_fraction"] if item["edge_same_sign_fraction"] is not None else -1, reverse=True)
    for item in ranked_metrics[:min(5, len(ranked_metrics))]:
        ratio = item["explained_variance_ratio"]
        ratio_text = f"; variance ratio {ratio:.4f}" if ratio is not None else ""
        lines.append(f"- PC{item['pc']}: same-sign edge fraction {item['edge_same_sign_fraction']:.3f}{ratio_text}")
    lines.extend(["", "## Leading PC genes", ""])
    for pc_index in range(min(5, n_pcs)):
        values = loadings[:, pc_index]
        order = np.argsort(np.abs(values))[::-1]
        positive_genes = [str(shared[i]) for i in order if values[i] > 0][:top_n]
        negative_genes = [str(shared[i]) for i in order if values[i] < 0][:top_n]
        positive = ", ".join(positive_genes)
        negative = ", ".join(negative_genes)
        lines.append(f"- PC{pc_index + 1}: positive loading genes are {positive}; negative loading genes are {negative}")
    lines.extend(["", "## Files", "", "- `overview.png`: visual PCA and gene-graph interpretation",
                  "- `pca_gene_loadings.csv`: ranked genes and local graph neighborhood statistics",
                  "- `pca_graph_metrics.json`: per-PC graph coherence metrics", "- `summary.md`: this interpretation"])
    summary.write_text("\n".join(lines) + "\n")
    return _result(outputs={"loadings": str(table), "metrics": str(metrics_path), "summary": str(summary), "overview": str(overview)},
                   metrics={"shared_genes": len(shared), "pcs": n_pcs, "top_n_per_side": top_n},
                   warnings=["Graph coherence describes co-expression structure around PCA loadings; it does not establish regulation or causality."],
                   inputs=[str(pca_path), str(gene_graph_path)])


def integrate_graph_evidence(ctx):
    """Compare gene-graph neighborhoods with PCA, markers, and regulon candidates."""
    import anndata as ad
    import pandas as pd
    from scipy import sparse

    artifacts = ctx.get("artifacts", {})
    root = _dir(ctx, "graph") / "integration"
    root.mkdir(parents=True, exist_ok=True)
    pca_path = artifacts.get("representation", {}).get("outputs", {}).get("dataset")
    gene_graph_path = artifacts.get("gene_graph", {}).get("outputs", {}).get("dataset")
    if not gene_graph_path:
        gene_graph_path = artifacts.get("graph", {}).get("outputs", {}).get("dataset")
    if not pca_path or not gene_graph_path:
        return _result(status="skipped", warnings=["PCA and gene-graph artifacts are required for integration."])

    gene_graph = ad.read_h5ad(gene_graph_path)
    graph_genes = pd.Index(gene_graph.obs_names.astype(str))
    adjacency = sparse.csr_matrix(gene_graph.obsp["connectivities"])

    def evidence_metrics(label, genes):
        unique = list(dict.fromkeys(str(gene) for gene in genes if pd.notna(gene)))
        indices = graph_genes.get_indexer(unique)
        indices = indices[indices >= 0]
        subgraph = adjacency[indices][:, indices] if len(indices) else sparse.csr_matrix((0, 0))
        return {"source": label, "genes_considered": len(unique), "genes_in_gene_graph": int(len(indices)),
                "within_set_edges": int(subgraph.nnz // 2),
                "within_set_edge_density": float((subgraph.nnz // 2) / max(len(indices) * (len(indices) - 1) / 2, 1))}

    interpretation = interpret_pca(ctx, pca_dataset=pca_path, gene_graph_dataset=gene_graph_path)
    outputs = dict(interpretation["outputs"])
    outputs["pca_summary"] = outputs.pop("summary")
    outputs["pca_overview"] = outputs.pop("overview")
    marker_path = artifacts.get("clustering", {}).get("outputs", {}).get("markers")
    regulon_path = artifacts.get("regulon", {}).get("outputs", {}).get("edges")
    comparison = []
    if marker_path and Path(marker_path).is_file():
        markers = pd.read_csv(marker_path)
        if {"group", "names"}.issubset(markers.columns):
            for group, table in markers.groupby("group", observed=True):
                genes = table["names"].astype(str).dropna().drop_duplicates().head(50).tolist()
                item = evidence_metrics("cluster_markers", genes)
                item["group"] = str(group)
                item["note"] = "Marker genes are descriptive; graph density measures co-expression structure, not identity proof."
                comparison.append(item)
    if regulon_path and Path(regulon_path).is_file():
        regulons = pd.read_csv(regulon_path, sep="\t")
        if {"tf", "target"}.issubset(regulons.columns):
            for tf, table in regulons.groupby("tf", observed=True):
                genes = [tf] + table["target"].astype(str).dropna().drop_duplicates().tolist()
                item = evidence_metrics("regulon_candidates", genes)
                item["group"] = str(tf)
                item["note"] = "Candidate TF-target edges are hypotheses; graph density is not regulatory validation."
                comparison.append(item)
    comparison_path = root / "evidence_comparison.csv"
    pd.DataFrame(comparison, columns=["source", "group", "genes_considered", "genes_in_gene_graph",
                                      "within_set_edges", "within_set_edge_density", "note"]).to_csv(comparison_path, index=False)
    summary = root / "summary.md"
    summary.write_text("\n".join([
        "# Integrated graph evidence", "",
        "This task keeps the cell graph, gene graph, PCA loadings, cluster markers, and candidate regulons as distinct evidence types.",
        "The PCA interpretation links gene neighborhoods to component loadings; the comparison table measures marker and regulon overlap with the gene graph.",
        "These overlaps describe consistency of evidence and do not establish regulation, interaction, or causality.", "",
        "## Files", "", "- `pca_summary`: PCA and gene-graph interpretation",
        "- `evidence_comparison.csv`: overlap and within-set graph density for marker/regulon gene sets", ""
    ]))
    outputs["comparison"] = str(comparison_path)
    outputs["summary"] = str(summary)
    return _result(outputs=outputs, metrics={"comparison_rows": len(comparison)},
                   warnings=["Integrated graph evidence is exploratory and does not independently validate biological mechanisms."],
                   inputs=[str(pca_path), str(gene_graph_path)] + ([str(marker_path)] if marker_path else []) + ([str(regulon_path)] if regulon_path else []))


def export_graph_metrics_md(
    graph_input,
    output_dir="multimodal analysis output",
    roi_nodes=None,
    top_k=5,
    node_metadata=None,
    label_cols=None
):
    """
    Computes graph theory metrics using python-igraph and exports a Markdown report.

    Parameters:
    -----------
    graph_input : igraph.Graph, networkx.Graph, or pd.DataFrame
        Input graph structure.
    output_dir : str
        Target directory where 'graph_metrics.md' will be saved.
    roi_nodes : iterable, optional
        List/set of ROI node identifiers to compute metrics on the ROI subgraph.
    top_k : int, optional
        Number of top nodes to display for each centrality ranking.
    node_metadata : pd.DataFrame, optional
        DataFrame containing node metadata indexed by cell/node ID (e.g., data.obs).
    label_cols : list[str] or str, optional
        List of column names from node_metadata to include as descriptive node attributes.
    """
    import igraph as ig
    import numpy as np

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "graph_metrics.md")

    # Format label_cols parameter
    if isinstance(label_cols, str):
        label_cols = [label_cols]
    elif label_cols is None:
        label_cols = []

    # 1. Input conversion to igraph.Graph
    if isinstance(graph_input, pd.DataFrame):
        has_weight = "weight" in graph_input.columns
        cols = ["source", "target"] + (["weight"] if has_weight else [])
        G = ig.Graph.TupleList(graph_input[cols].itertuples(index=False), directed=False, weights=has_weight)
    elif hasattr(graph_input, "to_undirected"):
        G = ig.Graph.from_networkx(graph_input)
    elif isinstance(graph_input, ig.Graph):
        G = graph_input
    else:
        raise ValueError("Unsupported input type. Provide an igraph.Graph, networkx.Graph, or pandas DataFrame.")

    # Retrieve node names
    if "name" in G.vs.attributes():
        node_names = G.vs["name"]
    elif "_nx_name" in G.vs.attributes():
        node_names = G.vs["_nx_name"]
    else:
        node_names = [str(i) for i in range(G.vcount())]

    G.vs["node_id"] = [str(n) for n in node_names]

    # 2. Isolate ROI if provided
    if roi_nodes is not None:
        roi_set = set(str(n) for n in roi_nodes)
        valid_indices = [v.index for v in G.vs if v["node_id"] in roi_set]
        target_G = G.subgraph(valid_indices)
        analysis_scope = f"ROI Analysis ({len(valid_indices)} valid nodes)"
    else:
        target_G = G
        analysis_scope = "Global Graph Analysis"

    num_nodes = target_G.vcount()
    num_edges = target_G.ecount()

    if num_nodes == 0:
        print("Warning: Graph/ROI contains no valid nodes. Report generation skipped.")
        return

    # 3. Global Topological Statistics
    density = target_G.density()
    is_directed = target_G.is_directed()

    components = target_G.connected_components(mode="weak" if is_directed else "strong")
    largest_comp_subgraph = components.giant()

    try:
        avg_path_length = largest_comp_subgraph.average_path_length(directed=is_directed)
        diameter = largest_comp_subgraph.diameter(directed=is_directed)
    except Exception:
        avg_path_length = "N/A"
        diameter = "N/A"

    avg_clustering = target_G.transitivity_avglocal_undirected(mode="zero")

    # 4. Centrality Calculations
    denom = max(1, num_nodes - 1)
    degrees = target_G.degree()
    degree_cent = {v["node_id"]: d / denom for v, d in zip(target_G.vs, degrees)}

    betweenness_scores = target_G.betweenness(directed=is_directed, normalized=True)
    betweenness_cent = {v["node_id"]: b for v, b in zip(target_G.vs, betweenness_scores)}

    closeness_scores = target_G.closeness(normalized=True)
    closeness_cent = {v["node_id"]: c for v, c in zip(target_G.vs, closeness_scores)}

    def get_top(metric_dict, k):
        return sorted(metric_dict.items(), key=lambda x: x[1], reverse=True)[:k]

    top_hubs = get_top(degree_cent, top_k)
    top_bridges = get_top(betweenness_cent, top_k)
    top_closeness = get_top(closeness_cent, top_k)

    # 5. Helper function to build dynamic Markdown ranking tables
    def build_ranking_table(top_list, metric_title):
        headers = ["Rank", "Node"] + [c.replace("_", " ").title() for c in label_cols] + ["Score"]
        table_lines = [
            f"\n### Top {top_k} {metric_title}",
            "| " + " | ".join(headers) + " |",
            "| " + " | ".join([":---"] * len(headers)) + " |"
        ]

        for idx, (node, score) in enumerate(top_list, 1):
            row = [str(idx), f"`{node}`"]

            # Lookup metadata row for current node
            meta_row = None
            if node_metadata is not None and not node_metadata.empty:
                if node in node_metadata.index:
                    meta_row = node_metadata.loc[node]
                elif str(node) in node_metadata.index:
                    meta_row = node_metadata.loc[str(node)]
                elif str(node).isdigit() and int(node) in node_metadata.index:
                    meta_row = node_metadata.loc[int(node)]

            for col in label_cols:
                if meta_row is not None and col in meta_row:
                    val = meta_row[col]
                    row.append(f"`{val:.2f}`" if isinstance(val, (float, np.floating)) else f"`{val}`")
                else:
                    row.append("`N/A`")

            row.append(f"`{score:.4f}`")
            table_lines.append("| " + " | ".join(row) + " |")

        return table_lines

    # 6. Build Markdown Content
    largest_clique_size = target_G.clique_number()
    try:
        num_cliques = len(target_G.maximal_cliques())
    except Exception:
        num_cliques = "N/A"

    md = [
        "# Graph Theory Metrics Summary\n",
        f"**Scope:** `{analysis_scope}`  ",
        f"**Graph Type:** `{'Directed' if is_directed else 'Undirected'}`\n",
        "---",
        "## 1. Global Topology & Connectivity",
        f"- **Total Nodes:** `{num_nodes}`",
        f"- **Total Edges:** `{num_edges}`",
        f"- **Graph Density:** `{density:.4f}`",
        f"- **Connected Components:** `{len(components)}`",
        f"- **Average Clustering Coefficient:** `{avg_clustering:.4f}`",
        f"- **Largest Component Diameter:** `{diameter}`",
        f"- **Largest Component Avg Path Length:** `{avg_path_length if isinstance(avg_path_length, str) else f'{avg_path_length:.4f}'}`\n",
        "## 2. Key Centrality Rankings"
    ]

    md.extend(build_ranking_table(top_hubs, "Hubs (Degree Centrality)"))
    md.extend(build_ranking_table(top_bridges, "Bridges (Betweenness Centrality)"))
    md.extend(build_ranking_table(top_closeness, "Closest Nodes (Closeness Centrality)"))

    md.extend([
        "\n## 3. Substructures & Cliques",
        f"- **Total Maximal Cliques:** `{num_cliques}`",
        f"- **Largest Clique Size:** `{largest_clique_size}` nodes"
    ])

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(md))

    print(f"[+] Metrics successfully written to {output_file}")