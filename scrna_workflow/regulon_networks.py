"""Cell-type views of the SCENIC+ eRegulon network: TF -> region -> target gene.

SCENIC+ infers one network from all cells. A cell-type view keeps the part of it
with evidence in that cell type:

- eRegulons: the `top_eregulons` direct gene-based eRegulons with the highest regulon
  specificity score (RSS) for the cell type;
- regions: triplets whose region overlaps a MACS2 peak called on that cell type's
  pseudobulk (accessible in the cell type);
- target genes: detected (count > 0) in at least `min_gene_fraction` of the cell type's
  cells, keeping per TF the `max_targets_per_tf` targets with the highest TF-gene
  importance.

The views are for navigation and hypothesis generation; they are not separate
cell-type-specific network inferences.
"""
from pathlib import Path

NODE_STYLE = {
    "tf": {"color": "#e08a3c", "marker": "D", "size": 230, "label": "Transcription factor"},
    "region": {"color": "#5b8cc0", "marker": "s", "size": 34, "label": "Regulatory region"},
    "gene": {"color": "#6aa56a", "marker": "o", "size": 90, "label": "Target gene"},
}
POSITIVE, NEGATIVE = "#c0504d", "#4f81bd"


def read_peaks(path):
    """narrowPeak/BED intervals per chromosome as (starts, running max of ends), sorted by start."""
    import numpy as np
    import pandas as pd
    peaks = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"],
                        dtype={"chrom": str})
    index = {}
    for chrom, rows in peaks.sort_values(["chrom", "start"]).groupby("chrom"):
        index[chrom] = (rows["start"].to_numpy(), np.maximum.accumulate(rows["end"].to_numpy()))
    return index


def overlapping_regions(region_names, peak_index):
    """Region names ("chrom:start-end", half-open) overlapping at least one peak."""
    import numpy as np
    from .scenicplus_stages.common import parse_region
    hits = set()
    for name in set(region_names):
        chrom, start, end = parse_region(name)
        if chrom not in peak_index:
            continue
        starts, max_ends = peak_index[chrom]
        n = int(np.searchsorted(starts, end, side="left"))  # peaks starting before the region ends
        if n and max_ends[n - 1] > start:
            hits.add(name)
    return hits


def detected_fraction(counts, cell_mask, genes, var_names):
    """Fraction of selected cells with a nonzero count, for the requested genes present in var_names."""
    import numpy as np
    from scipy import sparse
    positions = {g: i for i, g in enumerate(var_names)}
    present = [g for g in genes if g in positions]
    if not present or not np.any(cell_mask):
        return {}
    matrix = sparse.csr_matrix(counts)[np.flatnonzero(cell_mask)][:, [positions[g] for g in present]]
    fractions = np.asarray((matrix > 0).mean(axis=0)).ravel()
    return dict(zip(present, fractions.astype(float)))


def cell_type_network(triplets, rss_row, accessible, gene_fraction, top_eregulons=5, min_gene_fraction=0.1,
                      max_targets_per_tf=20):
    """Select triplets for one cell type; returns (nodes, edges, selected eRegulon signatures)."""
    import pandas as pd

    selected = list(rss_row.dropna().nlargest(top_eregulons).index)
    rows = triplets[triplets["gene_signature"].isin(selected)]
    rows = rows[rows["region"].isin(accessible)]
    rows = rows[rows["target"].map(lambda g: gene_fraction.get(g, 0.0) >= min_gene_fraction)]
    if len(rows):
        strength = rows.groupby(["tf", "target"])["importance_tf2g"].max().reset_index()
        strength = strength.sort_values(["tf", "importance_tf2g"], ascending=[True, False])
        keep = strength.groupby("tf").head(max_targets_per_tf)[["tf", "target"]]
        rows = rows.merge(keep, on=["tf", "target"])

    def joined(values):
        return ";".join(sorted(set(map(str, values))))

    edges = []
    for (tf, region), group in rows.groupby(["tf", "region"]):
        edges.append({"source": f"TF:{tf}", "target": f"region:{region}", "interaction": "tf_binds_region",
                      "eregulons": joined(group["eregulon"]), "importance": None, "rho": None})
    for (region, gene), group in rows.groupby(["region", "target"]):
        best = group.loc[group["importance_r2g"].idxmax()]
        edges.append({"source": f"region:{region}", "target": f"gene:{gene}", "interaction": "region_regulates_gene",
                      "eregulons": joined(group["eregulon"]), "importance": float(best["importance_r2g"]),
                      "rho": float(best["rho_r2g"])})
    for (tf, gene), group in rows.groupby(["tf", "target"]):
        best = group.loc[group["importance_tf2g"].idxmax()]
        edges.append({"source": f"TF:{tf}", "target": f"gene:{gene}", "interaction": "tf_regulates_gene",
                      "eregulons": joined(group["eregulon"]), "importance": float(best["importance_tf2g"]),
                      "rho": float(best["rho_tf2g"])})
    edges = pd.DataFrame(edges, columns=["source", "target", "interaction", "eregulons", "importance", "rho"])

    nodes = []
    for tf, group in rows.groupby("tf"):
        nodes.append({"id": f"TF:{tf}", "label": tf, "node_type": "tf", "eregulons": joined(group["eregulon"]),
                      "detected_fraction": gene_fraction.get(tf), "n_targets": int(group["target"].nunique())})
    for region, group in rows.groupby("region"):
        nodes.append({"id": f"region:{region}", "label": region, "node_type": "region",
                      "eregulons": joined(group["eregulon"]), "detected_fraction": None,
                      "n_targets": int(group["target"].nunique())})
    for gene, group in rows.groupby("target"):
        nodes.append({"id": f"gene:{gene}", "label": gene, "node_type": "gene", "eregulons": joined(group["eregulon"]),
                      "detected_fraction": gene_fraction.get(gene), "n_targets": 0})
    nodes = pd.DataFrame(nodes, columns=["id", "label", "node_type", "eregulons", "detected_fraction", "n_targets"])
    return nodes, edges, selected


def write_graphml(nodes, edges, path):
    """GraphML for Cytoscape/Gephi; missing numeric attributes are omitted."""
    import math
    import networkx as nx

    def clean(record):
        return {k: v for k, v in record.items() if v is not None and not (isinstance(v, float) and math.isnan(v))}
    graph = nx.MultiDiGraph()
    for record in nodes.to_dict("records"):
        graph.add_node(record.pop("id"), **clean(record))
    for record in edges.to_dict("records"):
        graph.add_edge(record.pop("source"), record.pop("target"), **clean(record))
    nx.write_graphml(graph, path)
    return str(path)


def draw_network(nodes, edges, destination, title):
    """Layered drawing: TFs (left) -> regions (middle) -> target genes (right).

    Genes are grouped by their strongest TF and regions follow their genes, which keeps
    each TF's module together. Region-gene edges are coloured by correlation sign.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np

    tf_gene = edges[edges["interaction"] == "tf_regulates_gene"]
    region_gene = edges[edges["interaction"] == "region_regulates_gene"]
    tf_region = edges[edges["interaction"] == "tf_binds_region"]
    tfs = nodes[nodes["node_type"] == "tf"].sort_values(["n_targets", "label"], ascending=[False, True])["id"].tolist()
    tf_rank = {tf: i for i, tf in enumerate(tfs)}
    owner = (tf_gene.assign(rank=tf_gene["source"].map(tf_rank))
             .sort_values(["target", "importance", "rank"], ascending=[True, False, True])
             .drop_duplicates("target").set_index("target")["rank"])
    genes = sorted(nodes.loc[nodes["node_type"] == "gene", "id"], key=lambda g: (owner.get(g, len(tfs)), g))
    gene_rank = {g: i for i, g in enumerate(genes)}
    region_order = region_gene.assign(rank=region_gene["target"].map(gene_rank)).groupby("source")["rank"].mean()
    regions = sorted(nodes.loc[nodes["node_type"] == "region", "id"], key=lambda r: (region_order.get(r, 0), r))

    def spread(ids):
        return {node: (1.0 - (i + 0.5) / len(ids)) for i, node in enumerate(ids)} if ids else {}
    x_of = {"tf": 0.0, "region": 1.0, "gene": 2.0}
    position = {}
    for kind, ids in (("tf", tfs), ("region", regions), ("gene", genes)):
        position.update({node: (x_of[kind], y) for node, y in spread(ids).items()})

    tall = max(len(regions), len(genes), len(tfs) * 3)
    fig, ax = plt.subplots(figsize=(11, min(40, max(5, 0.16 * tall + 1.5))))
    for _, edge in tf_region.iterrows():
        (x0, y0), (x1, y1) = position[edge["source"]], position[edge["target"]]
        ax.plot([x0, x1], [y0, y1], color="#b0b0b0", lw=0.4, alpha=0.6, zorder=1)
    scale = max(region_gene["rho"].abs().max(), 1e-9) if len(region_gene) else 1.0
    for _, edge in region_gene.iterrows():
        (x0, y0), (x1, y1) = position[edge["source"]], position[edge["target"]]
        ax.plot([x0, x1], [y0, y1], color=POSITIVE if edge["rho"] >= 0 else NEGATIVE,
                lw=0.4 + 1.6 * abs(edge["rho"]) / scale, alpha=0.7, zorder=1)
    for kind, ids in (("tf", tfs), ("region", regions), ("gene", genes)):
        style = NODE_STYLE[kind]
        if ids:
            ax.scatter([position[i][0] for i in ids], [position[i][1] for i in ids], s=style["size"], c=style["color"],
                       marker=style["marker"], edgecolors="white", linewidths=0.6, zorder=3)
    labels = nodes.set_index("id")["label"]
    for node in tfs:
        ax.text(position[node][0] - 0.06, position[node][1], labels[node], ha="right", va="center", fontsize=9,
                fontweight="bold")
    for node in genes:
        ax.text(position[node][0] + 0.05, position[node][1], labels[node], ha="left", va="center",
                fontsize=6 if len(genes) > 40 else 8)
    if len(regions) <= 40:
        for node in regions:
            ax.text(position[node][0], position[node][1] + 0.35 / max(len(regions), 1), labels[node], ha="center",
                    va="bottom", fontsize=5, color="#555555")
    handles = [Line2D([], [], marker=s["marker"], color="none", markerfacecolor=s["color"], markeredgecolor="white",
                      markersize=9, label=s["label"]) for s in NODE_STYLE.values()]
    handles += [Line2D([], [], color="#b0b0b0", lw=1, label="TF motif in region"),
                Line2D([], [], color=POSITIVE, lw=1.5, label="Region-gene, positive correlation"),
                Line2D([], [], color=NEGATIVE, lw=1.5, label="Region-gene, negative correlation")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=3, fontsize=7, frameon=False)
    ax.set_xlim(-0.75, 2.75)
    ax.set_ylim(-0.02, 1.02)
    ax.axis("off")
    ax.set_title(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(destination, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(destination)


def build_cell_type_networks(triplets_file, rss, cells_tsv, macs2_dir, adata, cell_type_column, out_dir,
                             top_eregulons=5, min_gene_fraction=0.1, max_targets_per_tf=20):
    """Write nodes/edges TSV, GraphML and PNG per pseudobulk cell type; returns (outputs, summary rows, warnings)."""
    import numpy as np
    import pandas as pd
    from .scenicplus_stages.common import safe_label

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    triplets = pd.read_csv(triplets_file, sep="\t")
    missing = {"tf", "target", "region", "eregulon", "gene_signature", "importance_tf2g", "rho_tf2g",
               "importance_r2g", "rho_r2g"} - set(triplets.columns)
    if missing:
        return {}, [], [f"eRegulon networks skipped: triplet table lacks columns {sorted(missing)}."]
    cells = pd.read_csv(cells_tsv, sep="\t", dtype=str, keep_default_na=False)
    groups = cells.loc[cells["pseudobulk_group"] != "", ["cell_type", "pseudobulk_group"]].drop_duplicates()
    labels = adata.obs[cell_type_column].astype(str).to_numpy()
    genes = pd.unique(pd.concat([triplets["tf"], triplets["target"]]).astype(str))
    var_names = adata.var_names.astype(str)
    outputs, summary, warnings = {}, [], []
    for cell_type, group in sorted(groups.itertuples(index=False, name=None)):
        peaks_file = Path(macs2_dir) / f"{group}_peaks.narrowPeak"
        if cell_type not in rss.index or not peaks_file.is_file():
            warnings.append(f"eRegulon network for {cell_type} skipped: no RSS row or pseudobulk peaks ({peaks_file.name}).")
            continue
        accessible = overlapping_regions(triplets["region"], read_peaks(peaks_file))
        fraction = detected_fraction(adata.layers["counts"], labels == cell_type, genes, var_names)
        nodes, edges, selected = cell_type_network(triplets, rss.loc[cell_type], accessible, fraction,
                                                   top_eregulons, min_gene_fraction, max_targets_per_tf)
        stem = safe_label(cell_type)
        row = {"cell_type": cell_type, "eregulons_selected": ";".join(selected),
               "tfs": int((nodes["node_type"] == "tf").sum()), "regions": int((nodes["node_type"] == "region").sum()),
               "target_genes": int((nodes["node_type"] == "gene").sum()), "edges": int(len(edges)),
               "cells": int(np.sum(labels == cell_type))}
        summary.append(row)
        if not len(nodes):
            warnings.append(f"eRegulon network for {cell_type} is empty after accessibility and expression filters.")
            continue
        nodes.to_csv(out_dir / f"{stem}_nodes.tsv", sep="\t", index=False)
        edges.to_csv(out_dir / f"{stem}_edges.tsv", sep="\t", index=False)
        outputs[f"network_nodes_{stem}"] = str(out_dir / f"{stem}_nodes.tsv")
        outputs[f"network_edges_{stem}"] = str(out_dir / f"{stem}_edges.tsv")
        outputs[f"network_graphml_{stem}"] = write_graphml(nodes, edges, out_dir / f"{stem}.graphml")
        outputs[f"network_plot_{stem}"] = draw_network(
            nodes, edges, out_dir / f"{stem}.png",
            f"{cell_type}: top {len(selected)} eRegulons by RSS ({row['tfs']} TFs, {row['regions']} regions, "
            f"{row['target_genes']} genes)")
    if summary:
        pd.DataFrame(summary).to_csv(out_dir / "summary.tsv", sep="\t", index=False)
        outputs["network_summary"] = str(out_dir / "summary.tsv")
    return outputs, summary, warnings
