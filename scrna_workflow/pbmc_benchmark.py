"""Descriptive same-data PBMC reference-label comparison and cell-type marker test.

This benchmark does not perform a disease/condition contrast and never treats
individual cells as biological replicates. Published labels are an internal
comparison on these same cells, not independent biological validation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
from datetime import datetime, timezone


def _hash(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def _broad(label):
    """Explicit lineage-only harmonization; no subtype/state assumptions."""
    s=str(label).strip().lower().replace("_"," ")
    if s in ("unknown","ambiguous","unavailable","nan",""):
        return None
    if s in ("b", "b cells", "b cell"):
        return "B cells"
    if s in ("t", "t cells", "t cell", "cd4 t cells", "cd8 t cells", "cd4+ t cells", "cd8+ t cells"):
        return "T cells"
    if s in ("nk", "nk cells", "nk cell", "natural killer cells"):
        return "NK cells"
    if s in ("monocytes", "monocyte", "cd14+ monocytes", "fcgr3a+ monocytes", "cd14 monocytes", "fcgr3a monocytes"):
        return "Monocytes"
    if s in ("dendritic cells", "dendritic", "dc"):
        return "Dendritic cells"
    if s in ("megakaryocytes", "megakaryocyte", "platelets", "platelet"):
        return "Platelet/megakaryocyte"
    return None


def benchmark(run_dir, output, group="B cells", reference="CD4 T cells", seed=0):
    run_dir=Path(run_dir).resolve(); output=Path(output).resolve()
    source=run_dir/"clustering"/"annotated.h5ad"
    if not source.exists():
        raise FileNotFoundError(f"Missing clustering artifact: {source}")
    if output==source.parent or output in source.parents:
        raise ValueError("Benchmark output must not overwrite parent workflow artifact directories")
    output.mkdir(parents=True,exist_ok=True)
    os.environ.setdefault("NUMBA_CACHE_DIR",str(output/".cache"/"numba"))
    os.environ.setdefault("MPLCONFIGDIR",str(output/".cache"/"matplotlib"))
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a=sc.read_h5ad(source)
    for col in ("reference_cell_type","cluster"):
        if col not in a.obs:
            raise ValueError(f"Required metadata missing: {col}")
    valid=a.obs["reference_cell_type"].notna() & ~a.obs["reference_cell_type"].astype(str).str.lower().isin(["unavailable","unknown", "nan", ""])
    b=a[valid.to_numpy()].copy()
    if not b.n_obs:
        raise ValueError("No cells have available reference annotations")
    labels=b.obs["reference_cell_type"].astype(str)
    groups=labels.value_counts()
    if group not in groups or reference not in groups or min(groups[group],groups[reference])<2:
        raise ValueError(f"DEG comparison needs at least two cells per label. Available labels: {groups.to_dict()}")
    if "log_expression" not in b.layers:
        raise ValueError("Benchmark needs workflow log_expression layer; raw scale must not be guessed")
    contingency=pd.crosstab(labels,b.obs["cluster"])
    contingency.to_csv(output/"reference_by_cluster.tsv",sep="\t")
    metrics=dict(total_workflow_cells=a.n_obs,reference_labeled_cells=b.n_obs,
                 excluded_without_reference=int((~valid).sum()),
                 reference_counts={str(k):int(v) for k,v in groups.items()},
                 cluster_reference_ARI=float(adjusted_rand_score(labels,b.obs["cluster"].astype(str))),
                 cluster_reference_NMI=float(normalized_mutual_info_score(labels,b.obs["cluster"].astype(str))),
                 reference_scope="Same dataset annotations; not independent validation",seed=seed)
    inferred_col="annotation" if "annotation" in b.obs else "cell_type" if "cell_type" in b.obs else None
    mapping={str(v):_broad(v) for v in labels.unique()}
    if inferred_col:
        ref=labels.map(_broad); inferred=b.obs[inferred_col].astype(str).map(_broad)
        eligible=ref.notna() & inferred.notna()
        metrics["annotation_broad_mapping"]={"reference":mapping,"inferred":{str(v):_broad(v) for v in b.obs[inferred_col].unique()}}
        metrics["broad_annotation_coverage"]=float(eligible.mean())
        metrics["broad_annotation_accuracy_on_mapped_cells"]=float((ref[eligible]==inferred[eligible]).mean()) if eligible.any() else None
        pd.crosstab(ref.fillna("unmapped"),inferred.fillna("unmapped")).to_csv(output/"reference_by_annotation_broad.tsv",sep="\t")
    b.obs["reference_cell_type"]=pd.Categorical(labels)
    sc.tl.rank_genes_groups(b,groupby="reference_cell_type",groups=[group],reference=reference,
                           method="wilcoxon",tie_correct=True,use_raw=False,layer="log_expression",
                           corr_method="benjamini-hochberg",pts=True,key_added="benchmark_markers")
    deg=sc.get.rank_genes_groups_df(b,group=group,key="benchmark_markers")
    expr=b.layers["log_expression"]
    group_mask=(labels==group).to_numpy(); ref_mask=(labels==reference).to_numpy()
    def fractions(mask):
        return np.asarray((expr[mask]>0).mean(axis=0)).ravel()
    gf=pd.Series(fractions(group_mask),index=b.var_names); rf=pd.Series(fractions(ref_mask),index=b.var_names)
    deg["fraction_detected_group"]=deg["names"].map(gf)
    deg["fraction_detected_reference"]=deg["names"].map(rf)
    deg["group"]=group; deg["reference"]=reference
    deg["inference_scope"]="exploratory cell-level cell-type marker comparison; no biological replication"
    deg.to_csv(output/"exploratory_cell_type_deg.tsv",sep="\t",index=False)
    metrics["deg"]=dict(group=group,reference=reference,group_cells=int(group_mask.sum()),reference_cells=int(ref_mask.sum()),
                        tested_genes=len(deg),fdr_below_005=int((deg.pvals_adj<.05).sum()),
                        method="Scanpy Wilcoxon rank-sum with tie correction; BH across tested genes",
                        expression="log1p library-size-normalized expression",
                        logfoldchange="Scanpy approximation from means of log-normalized expression",
                        limitation="Cell-level marker statistics; cells are not donor replicates; no condition-level inference")

    fig,ax=plt.subplots(figsize=(7,5))
    x=deg.logfoldchanges.to_numpy(); y=-np.log10(np.maximum(deg.pvals_adj.to_numpy(),np.finfo(float).tiny))
    sig=(deg.pvals_adj<.05).to_numpy()
    ax.scatter(x[~sig],y[~sig],s=5,c="#bbbbbb",rasterized=True)
    ax.scatter(x[sig],y[sig],s=6,c="#2563a5",rasterized=True)
    for _,row in deg.reindex(deg.logfoldchanges.abs().sort_values(ascending=False).index).head(8).iterrows():
        ax.annotate(row["names"],(row.logfoldchanges,-np.log10(max(row.pvals_adj,np.finfo(float).tiny))),fontsize=7)
    ax.set(xlabel=f"Approximate log2 fold change: {group} / {reference}",ylabel="−log10 BH-adjusted cell-level p-value",
           title="Exploratory cell-type markers — same dataset")
    fig.tight_layout(); fig.savefig(output/"exploratory_deg_volcano.png",dpi=300); fig.savefig(output/"exploratory_deg_volcano.pdf"); plt.close(fig)

    ranked=deg.sort_values("scores",ascending=False)
    selected=list(dict.fromkeys(ranked.head(8)["names"].tolist()+ranked.tail(8)["names"].tolist()))
    celltypes=sorted(labels.unique())
    matrix=b[:,selected].layers["log_expression"]
    means=np.vstack([np.asarray(matrix[(labels==ct).to_numpy()].mean(axis=0)).ravel() for ct in celltypes])
    scale=means.std(axis=0); scale[scale==0]=1
    standardized=(means-means.mean(axis=0))/scale
    pd.DataFrame(means,index=celltypes,columns=selected).to_csv(output/"marker_mean_log_expression.tsv",sep="\t")
    fig,ax=plt.subplots(figsize=(11,max(3, .5*len(celltypes))))
    im=ax.imshow(standardized,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
    ax.set_xticks(range(len(selected)),selected,rotation=65,ha="right",fontsize=8)
    ax.set_yticks(range(len(celltypes)),celltypes,fontsize=9)
    ax.set_title("Selected contrast markers across reference cell types")
    fig.colorbar(im,ax=ax,label="Per-gene z score of mean log expression"); fig.tight_layout()
    fig.savefig(output/"marker_heatmap.png",dpi=300); fig.savefig(output/"marker_heatmap.pdf"); plt.close(fig)

    embedding="X_umap" if "X_umap" in b.obsm else "X_pca"
    z=b.obsm[embedding]
    fig,axes=plt.subplots(1,2,figsize=(14,6))
    fields=[("reference_cell_type","Published same-data reference"),(inferred_col or "cluster","Inferred annotation" if inferred_col else "Inferred cluster")]
    for ax,(column,title) in zip(axes,fields):
        cats=b.obs[column].astype(str)
        for i,label in enumerate(sorted(cats.unique())):
            mask=(cats==label).to_numpy()
            ax.scatter(z[mask,0],z[mask,1],s=4,alpha=.75,label=label,color=plt.get_cmap("tab20")(i%20),rasterized=True)
        ax.set(title=title,xlabel="UMAP 1" if embedding=="X_umap" else "PC 1",ylabel="UMAP 2" if embedding=="X_umap" else "PC 2")
        ax.legend(loc="upper center",bbox_to_anchor=(.5,-.12),fontsize=7,ncol=2,markerscale=3)
    fig.tight_layout(); fig.savefig(output/"reference_vs_inferred_embedding.png",dpi=300,bbox_inches="tight")
    fig.savefig(output/"reference_vs_inferred_embedding.pdf",bbox_inches="tight"); plt.close(fig)

    (output/"metrics.json").write_text(json.dumps(metrics,indent=2))
    top=ranked.head(10)
    bullets="\n".join(f"- {row['names']}: approximate log2 fold change {row.logfoldchanges:.2f}; detection {row.fraction_detected_group:.1%} vs {row.fraction_detected_reference:.1%}." for _,row in top.iterrows())
    accuracy=metrics.get("broad_annotation_accuracy_on_mapped_cells")
    accuracy_text=f"{accuracy:.1%} among mapped cells (coverage {metrics['broad_annotation_coverage']:.1%})" if accuracy is not None else "not estimable: no mapped annotation pairs"
    report=f"""# PBMC same-data annotation and marker benchmark

Analyzed {b.n_obs:,} cells with available reference labels; {int((~valid).sum()):,} retained workflow cells without reference labels were excluded from this comparison. Reference labels were never supplied to clustering by this benchmark.

Cluster/reference adjusted Rand index: **{metrics['cluster_reference_ARI']:.3f}**; normalized mutual information: **{metrics['cluster_reference_NMI']:.3f}**. Broad-lineage annotation agreement is {accuracy_text}. The exact conservative lineage mappings and denominators are in `metrics.json`; unmatched and unknown labels are excluded from agreement and retained in the contingency table. These are internal consistency measures against labels from the same expression dataset, not independent validation. Annotation tuning using these labels would invalidate any claim of held-out performance.

## Exploratory cell-type differential-expression test

Contrast: **{group}** ({group_mask.sum():,} cells) versus **{reference}** ({ref_mask.sum():,} cells). Wilcoxon rank-sum tests use tie correction on log1p library-size-normalized expression; Benjamini–Hochberg correction covers {len(deg):,} tested genes. {metrics['deg']['fdr_below_005']:,} genes have adjusted cell-level p < 0.05. Approximate log2 fold changes are Scanpy estimates derived from mean log expression. Detection fractions accompany every gene.

The strongest positive ranked genes are:

{bullets}

These statistics describe cell-type marker separation within one dataset. Cells are not independent biological samples. There is no replicated donor/condition comparison, and small adjusted p-values do not establish population-level differential expression. The reference cell types themselves were inferred from expression, so marker enrichment also has selection circularity.

## Artifacts and follow-up

- `reference_by_cluster.tsv`: cluster/reference contingency; ARI/NMI in `metrics.json`.
- `exploratory_cell_type_deg.tsv`: all tested genes, scores, approximate effects, adjusted p-values, and detection fractions.
- `exploratory_deg_volcano.png` and `.pdf`: exploratory marker effects and cell-level significance.
- `marker_heatmap.png` and `.pdf`: mean log-expression profiles of selected markers, standardized across reference cell types.
- `reference_vs_inferred_embedding.png` and `.pdf`: reference and inferred annotations on the same embedding. Embedding distances do not represent physical proximity.

Principal observation: reference cell types show expression separation quantified by the marker table. Confidence is descriptive within this dataset. Follow up with independent donors and predefined marker panels; use donor-level pseudobulk for any condition comparison. Biological annotation agreement should be tested on held-out samples or an independent reference. No causal, disease-associated, or regulatory mechanism is established by this benchmark.
"""
    (output/"REPORT.md").write_text(report)
    manifest=dict(created_at=datetime.now(timezone.utc).isoformat(),input_path=str(source),input_sha256=_hash(source),
                  parameters=dict(group=group,reference=reference,seed=seed),scanpy_version=sc.__version__,
                  benchmark_source_sha256=_hash(Path(__file__)),artifacts={p.name:_hash(p) for p in sorted(output.iterdir()) if p.is_file() and p.name not in ("manifest.json","tool_calls.jsonl")})
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2))
    calls=[dict(tool="scanpy.read_h5ad",input=str(source)),dict(tool="sklearn.metrics",methods=["adjusted_rand_score","normalized_mutual_info_score"]),
           dict(tool="scanpy.tl.rank_genes_groups",groupby="reference_cell_type",groups=[group],reference=reference,method="wilcoxon",tie_correct=True,layer="log_expression",use_raw=False,corr_method="benjamini-hochberg")]
    (output/"tool_calls.jsonl").write_text("".join(json.dumps(c)+"\n" for c in calls))
    return metrics


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir",required=True); p.add_argument("--output",required=True)
    p.add_argument("--group",default="B cells"); p.add_argument("--reference",default="CD4 T cells")
    p.add_argument("--seed",type=int,default=0)
    args=p.parse_args()
    print(json.dumps(benchmark(args.run_dir,args.output,args.group,args.reference,args.seed),indent=2))


if __name__=="__main__":
    main()
