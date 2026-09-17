"""Local HTML report for QC, PCA, clustering, and marker evidence."""
from html import escape
from pathlib import Path
import json

import numpy as np
import pandas as pd


def _table(frame):
    return frame.to_html(index=False, escape=True, classes="data", border=0) if not frame.empty else "<p>No results available.</p>"


def _interactive_browser(adata, markers):
    """Embed a self-contained cluster and gene browser in the HTML report."""
    from plotly.offline import get_plotlyjs
    from scipy import sparse

    key = "X_umap" if "X_umap" in adata.obsm else "X_pca"
    xy = np.asarray(adata.obsm[key][:, :2], dtype=float)
    if xy.shape != (adata.n_obs, 2):
        raise ValueError("Embedding and cell identifiers are misaligned")
    clusters = adata.obs["cluster"].astype(str).tolist()
    labels = adata.obs["cell_type"].astype(str).tolist()
    if len(clusters) != len(labels) or any(not label for label in labels):
        raise ValueError("Cell annotations and embedding are misaligned")
    positive = markers.loc[markers["logfoldchanges"] > 0].copy() if not markers.empty else markers
    if not positive.empty:
        positive = positive.sort_values(["group", "pvals_adj"], na_position="last")
        genes = list(dict.fromkeys(positive.groupby("group", sort=False)["names"]
                                   .head(5).astype(str)))[:50]
    else:
        genes = []
    genes = [gene for gene in genes if gene in adata.var_names]
    expression = {}
    matrix = adata.layers["log_expression"] if "log_expression" in adata.layers else adata.X
    for gene in genes:
        column = matrix[:, adata.var_names.get_loc(gene)]
        values = column.toarray().ravel() if sparse.issparse(column) else np.asarray(column).ravel()
        expression[gene] = np.round(values, 3).tolist()
    payload = {"x": np.round(xy[:, 0], 3).tolist(),
               "y": np.round(xy[:, 1], 3).tolist(),
               "cell": adata.obs_names.astype(str).tolist(),
               "cluster": clusters, "label": labels, "expression": expression,
               "axis": "UMAP" if key == "X_umap" else "PC"}
    data = json.dumps(payload, separators=(",", ":"), ensure_ascii=True).replace("</", "<\\/")
    script = """<script>__PLOTLY__</script>
<script>
const browserData = __DATA__;
const clusterSelect = document.getElementById('cluster-select');
const geneSelect = document.getElementById('gene-select');
const colors = ['#286fb4','#e57832','#2a9d79','#a565b8','#c4515c','#8e7831','#2a9fb0','#777e91'];
const clusters = [...new Set(browserData.cluster)].sort((a,b)=>Number(a)-Number(b));
for (const id of clusters) {
  const option = document.createElement('option'); option.value=id;
  option.textContent=`Cluster ${id}: ${browserData.label[browserData.cluster.indexOf(id)]}`;
  clusterSelect.append(option);
}
for (const gene of Object.keys(browserData.expression).sort()) {
  const option=document.createElement('option'); option.value=gene; option.textContent=gene;
  geneSelect.append(option);
}
function redraw() {
  const selected=clusterSelect.value, gene=geneSelect.value;
  const traces=[];
  if (gene) {
    const indices=browserData.cell.map((_,i)=>i).filter(i=>selected==='all'||browserData.cluster[i]===selected);
    traces.push({type:'scattergl',mode:'markers',name:gene,
      x:indices.map(i=>browserData.x[i]),y:indices.map(i=>browserData.y[i]),
      customdata:indices.map(i=>[browserData.cell[i],browserData.cluster[i],browserData.label[i],browserData.expression[gene][i]]),
      hovertemplate:'Cell %{customdata[0]}<br>Cluster %{customdata[1]}: %{customdata[2]}<br>'+gene+': %{customdata[3]}<extra></extra>',
      marker:{size:6,opacity:.8,color:indices.map(i=>browserData.expression[gene][i]),colorscale:'Viridis',
        showscale:true,colorbar:{title:gene}}});
  } else {
    clusters.filter(id=>selected==='all'||id===selected).forEach((id,k)=>{
      const indices=browserData.cell.map((_,i)=>i).filter(i=>browserData.cluster[i]===id);
      traces.push({type:'scattergl',mode:'markers',name:`${id}: ${browserData.label[indices[0]]}`,
        x:indices.map(i=>browserData.x[i]),y:indices.map(i=>browserData.y[i]),
        customdata:indices.map(i=>browserData.cell[i]),
        hovertemplate:'Cell %{customdata}<br>Cluster '+id+': '+browserData.label[indices[0]]+'<extra></extra>',
        marker:{size:6,opacity:.75,color:colors[clusters.indexOf(id)%colors.length]}});
    });
  }
  Plotly.react('cell-browser',traces,{xaxis:{title:browserData.axis+' 1',zeroline:false},
    yaxis:{title:browserData.axis+' 2',zeroline:false},legend:{orientation:'h'},
    margin:{l:55,r:35,t:20,b:55},hovermode:'closest'},
    {responsive:true,displaylogo:false,modeBarButtonsToRemove:['sendDataToCloud']});
  const summary = gene ? clusters.filter(id=>selected==='all'||id===selected).map(id=>({
    type:'violin',name:`Cluster ${id}`,box:{visible:true},meanline:{visible:true},
    points:false,y:browserData.expression[gene].filter((_,i)=>browserData.cluster[i]===id),
    line:{color:colors[clusters.indexOf(id)%colors.length]}
  })) : [{type:'bar',x:clusters.filter(id=>selected==='all'||id===selected),
    y:clusters.filter(id=>selected==='all'||id===selected).map(id=>browserData.cluster.filter(x=>x===id).length),
    marker:{color:clusters.filter(id=>selected==='all'||id===selected).map(id=>colors[clusters.indexOf(id)%colors.length])}}];
  Plotly.react('gene-summary',summary,{title:gene ? `${gene} expression by cluster` : 'Cells by cluster',
    xaxis:{title:'Cluster'},yaxis:{title:gene ? 'Log-normalized expression' : 'Cells'},
    showlegend:false,margin:{l:65,r:20,t:55,b:55}},
    {responsive:true,displaylogo:false,modeBarButtonsToRemove:['sendDataToCloud']});
  document.querySelectorAll('#cluster-table tbody tr,#marker-table tbody tr').forEach(row=>{
    row.hidden=selected!=='all' && row.cells[0].textContent.trim()!==selected;
  });
}
clusterSelect.addEventListener('change',redraw); geneSelect.addEventListener('change',redraw);
redraw();
</script>"""
    return script.replace("__PLOTLY__", get_plotlyjs()).replace("__DATA__", data)


def render_cluster_report(adata, out, qc_outputs, markers, diagnostics, annotations, warnings):
    """Write a report using only artifacts from this run; return its file paths."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out)
    qc_path = qc_outputs.get("cell_qc")
    qc = pd.read_csv(qc_path, index_col=0) if qc_path and Path(qc_path).exists() else pd.DataFrame()
    passed = qc["qc_pass"].astype(str).str.lower().eq("true") if "qc_pass" in qc else pd.Series(dtype=bool)
    qc_summary = pd.DataFrame([
        {"measure": "Input cells", "value": len(qc) if len(qc) else adata.n_obs},
        {"measure": "QC passed", "value": int(passed.sum()) if len(passed) else adata.n_obs},
        {"measure": "Cells after filtering", "value": adata.n_obs},
        {"measure": "Genes after filtering", "value": adata.n_vars},
        {"measure": "RPL/RPS genes removed", "value": adata.uns.get("workflow_gene_filter", {}).get("removed_genes", 0)},
    ])
    for column, name in (("total_counts", "Median library size"),
                         ("n_genes_by_counts", "Median detected genes"),
                         ("pct_counts_mt", "Median mitochondrial %")):
        if column in qc and len(qc):
            values = pd.to_numeric(qc[column], errors="coerce")
            qc_summary.loc[len(qc_summary)] = [name, f"{values.median():.2f}"]

    figures = {}
    qc_columns = [name for name in ("total_counts", "n_genes_by_counts", "pct_counts_mt") if name in qc]
    if qc_columns:
        fig, axes = plt.subplots(1, len(qc_columns), figsize=(4 * len(qc_columns), 3.2), squeeze=False)
        for ax, column in zip(axes.ravel(), qc_columns):
            values = pd.to_numeric(qc[column], errors="coerce").dropna()
            ax.hist(values, bins=40, color="#367a9c")
            ax.set(xlabel=column, ylabel="Input cells")
        fig.tight_layout()
        path = out / "cluster_qc.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        figures["qc_figure"] = str(path)

    variance = np.asarray(adata.uns.get("pca", {}).get("variance_ratio", []), dtype=float)
    n_pcs = adata.obsm["X_pca"].shape[1]
    pc_table = pd.DataFrame({
        "PC": np.arange(1, min(len(variance), n_pcs) + 1),
        "Variance explained (%)": np.round(100 * variance[:n_pcs], 2),
    })
    if len(variance):
        fig, ax = plt.subplots(figsize=(7, 3.2))
        ax.bar(np.arange(1, len(variance) + 1), 100 * variance, color="#496ca6")
        ax.set(xlabel="Principal component", ylabel="Variance explained (%)",
               title=f"Selected PCs: {n_pcs}")
        fig.tight_layout()
        path = out / "selected_pcs.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        figures["pc_figure"] = str(path)

    cluster_sizes = adata.obs["cluster"].astype(str).value_counts().sort_index()
    annotation_by_cluster = {row["cluster"]: row for row in annotations}
    for cluster in cluster_sizes.index:
        observed_labels = adata.obs.loc[adata.obs["cluster"].astype(str) == cluster,
                                         "cell_type"].astype(str).unique().tolist()
        if (cluster not in annotation_by_cluster or len(observed_labels) != 1
                or observed_labels[0] != annotation_by_cluster[cluster]["label"]):
            raise ValueError(f"Cluster {cluster} annotation disagrees with cell metadata")
    cluster_table = pd.DataFrame([
        {"Cluster": cluster, "Cells": int(size),
         "Cell type": annotation_by_cluster[cluster]["label"],
         "Confidence": round(float(annotation_by_cluster[cluster]["confidence"]), 2),
         "Marker evidence": "; ".join(
             f"{name}: {', '.join(genes)}" for name, genes, _ in annotation_by_cluster[cluster]["evidence"] if genes
         ) or "none"}
        for cluster, size in cluster_sizes.items()
    ])
    marker_columns = [c for c in ("group", "names", "logfoldchanges", "pvals_adj") if c in markers]
    top_markers = (
        markers.sort_values(["group", "pvals_adj"], na_position="last")
        .groupby("group", observed=True, sort=False).head(10)[marker_columns].copy()
        if not markers.empty else pd.DataFrame(columns=marker_columns)
    )
    if "logfoldchanges" in top_markers:
        top_markers["logfoldchanges"] = top_markers["logfoldchanges"].round(3)
    if "pvals_adj" in top_markers:
        top_markers["pvals_adj"] = top_markers["pvals_adj"].map(lambda x: f"{x:.3g}" if pd.notna(x) else "")
    thresholds_path = qc_outputs.get("thresholds")
    thresholds = pd.DataFrame(json.loads(Path(thresholds_path).read_text())) if thresholds_path and Path(thresholds_path).exists() else pd.DataFrame()
    resolution = pd.DataFrame(diagnostics)
    warning_items = "".join(f"<li>{escape(str(w))}</li>" for w in warnings)
    qc_image = '<img src="cluster_qc.png" alt="QC distributions">' if "qc_figure" in figures else ""
    pc_image = '<img src="selected_pcs.png" alt="Variance explained by selected PCs">' if "pc_figure" in figures else ""
    browser = _interactive_browser(adata, markers)
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Clustering report</title>
<style>body{{font:15px/1.5 system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#172536}}
h1,h2{{color:#163b5e}} .data{{border-collapse:collapse;width:100%;margin:1rem 0 2rem}}
.data th,.data td{{padding:.45rem .6rem;border-bottom:1px solid #dce4eb;text-align:left}}
.data th{{background:#eef4f8}}img{{max-width:100%;height:auto;margin:1rem 0}}
.note{{background:#f3f7fa;padding:.7rem 1rem;border-left:4px solid #367a9c}}
.controls{{display:flex;flex-wrap:wrap;gap:1rem;margin:1rem 0}}
.controls label{{display:flex;flex-direction:column;gap:.25rem;font-weight:600}}
.controls select{{font:inherit;padding:.4rem;min-width:15rem}}
#cell-browser{{width:100%;height:610px}}#gene-summary{{width:100%;height:370px}}
tr[hidden]{{display:none}}</style></head>
<body><h1>Clustering and QC report</h1>
<p class="note">Counts and QC are observations. PCA, clusters, marker rankings and cell labels are computational inferences.
Marker p-values are exploratory cell-level rankings, not donor-aware condition tests.</p>
<h2>QC results</h2>{_table(qc_summary)}{qc_image}
<h3>Sample-aware QC thresholds</h3>{_table(thresholds)}
<h2>Selected principal components</h2><p>{n_pcs} PCs were used for the expression graph.</p>{pc_image}{_table(pc_table)}
<h2>Explore cells and markers</h2>
<p>Filter by cluster, color cells by a top marker gene, and hover to inspect cell IDs and annotations. The embedding shows expression similarity, not physical proximity. Gene colors show log-normalized expression where available.</p>
<div class="controls"><label>Cluster<select id="cluster-select"><option value="all">All clusters</option></select></label>
<label>Color by gene<select id="gene-select"><option value="">Cell annotation</option></select></label></div>
<div id="cell-browser" role="img" aria-label="Interactive cell embedding"></div>
<div id="gene-summary" role="img" aria-label="Interactive gene expression or cluster size summary"></div>
<p>Static export: <a href="embedding.png">annotated embedding</a>.</p>
<h2>Clusters and marker-based annotation</h2><div id="cluster-table">{_table(cluster_table)}</div>
<h3>Resolution diagnostics</h3>{_table(resolution)}
<h2>Top marker genes by cluster</h2><div id="marker-table">{_table(top_markers)}</div>
<p>Complete markers: <a href="markers.csv">markers.csv</a>. Labels need multiple supporting markers;
unsupported or tied labels remain unknown or ambiguous.</p>
<h2>Warnings and limits</h2><ul>{warning_items}</ul>{browser}</body></html>"""
    path = out / "report.html"
    path.write_text(html, encoding="utf-8")
    return {"html_report": str(path), **figures}
