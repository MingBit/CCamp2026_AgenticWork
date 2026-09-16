"""File-backed scientific specialists. No network access or implicit reference downloads."""
from pathlib import Path
import json
import hashlib
import shutil


def _result(outputs=None, metrics=None, warnings=None, status="completed", inputs=None):
    return dict(status=status, input_references=inputs or [], outputs=outputs or {},
                metrics=metrics or {}, warnings=warnings or [], recommended_next_actions=[])


def _dir(ctx, name):
    p = Path(ctx["root"]) / name
    p.mkdir(parents=True, exist_ok=True)
    return p


def _previous(ctx, task, key="dataset"):
    item = ctx["artifacts"][task]
    return item.get("outputs", item)[key]


def _json(path, obj):
    path.write_text(json.dumps(obj, indent=2, default=str))
    return str(path)


def _digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _countlike(x):
    import numpy as np
    from scipy import sparse
    values = x.data if sparse.issparse(x) else np.asarray(x).ravel()
    return bool(np.isfinite(values).all() and (values >= 0).all()
                and np.allclose(values, np.round(values), atol=1e-6, rtol=0))


def inspect_data(ctx):
    """Inventory an h5ad, 10x HDF5, or 10x Matrix Market directory; snapshot sources."""
    cfg = ctx["config"]
    source = cfg.get("input_path")
    if not source or not Path(source).exists():
        r = _result(status="blocked", warnings=["No existing input_path: provide h5ad or a 10x count matrix."])
        r["recommended_next_actions"] = ["Supply a dataset path; optionally a barcode-indexed metadata CSV/TSV."]
        return r
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse
    source = Path(source).resolve()
    out = _dir(ctx, "inspection")
    snap = out / "source"
    snap.mkdir(exist_ok=True)
    manifest = []
    sources = sorted(source.rglob("*")) if source.is_dir() else [source]
    for f in sources:
        if f.is_file():
            rel = f.relative_to(source) if source.is_dir() else Path(f.name)
            target = snap / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            digest = _digest(f)
            if target.exists() and _digest(target) != digest:
                raise ValueError("Immutable source snapshot differs from current input; start a new run")
            if not target.exists():
                shutil.copy2(f, target)
                target.chmod(0o444)
            manifest.append(dict(path=str(f), snapshot=str(target), sha256=digest, bytes=f.stat().st_size))
    source_copy = snap if source.is_dir() else snap / source.name
    if source.is_dir():
        a = sc.read_10x_mtx(source_copy, var_names="gene_symbols", make_unique=False)
    elif source.suffix.lower() == ".h5ad":
        a = sc.read_h5ad(source_copy)
    elif source.suffix.lower() in (".h5", ".hdf5"):
        a = sc.read_10x_h5(source_copy)
    else:
        raise ValueError("Supported formats: .h5ad, 10x .h5, and 10x Matrix Market directory")
    if a.n_obs < 3 or a.n_vars < 3:
        raise ValueError("At least three cells and three genes are required")
    if not a.obs_names.is_unique:
        raise ValueError("Duplicate cell identifiers require upstream resolution")
    warnings = []
    duplicate_genes = int(a.var_names.duplicated().sum())
    a.var["source_gene_id"] = a.var_names.astype(str)
    if duplicate_genes:
        a.var_names_make_unique()
        warnings.append("Duplicate gene labels made unique; source_gene_id preserves originals; no genes merged.")
    if cfg.get("metadata_path"):
        mp = Path(cfg["metadata_path"])
        m = pd.read_csv(mp, sep="\t" if mp.suffix in (".tsv", ".txt") else ",", index_col=0)
        m.index = m.index.astype(str)
        if not m.index.is_unique or not a.obs_names.isin(m.index).all() or len(m) != a.n_obs:
            raise ValueError("Metadata requires unique cell IDs and coverage of every expression barcode")
        m = m.loc[a.obs_names]
        for col in m:
            if col in a.obs and not a.obs[col].astype(str).equals(m[col].astype(str)):
                raise ValueError(f"Conflicting embedded and supplied metadata column: {col}")
            a.obs[col] = m[col]
        metadata_copy = snap / ("external_metadata" + mp.suffix)
        metadata_digest = _digest(mp)
        if metadata_copy.exists() and _digest(metadata_copy) != metadata_digest:
            raise ValueError("Immutable metadata snapshot differs from current input; start a new run")
        if not metadata_copy.exists():
            shutil.copy2(mp, metadata_copy); metadata_copy.chmod(0o444)
        manifest.append(dict(path=str(mp.resolve()), snapshot=str(metadata_copy), sha256=metadata_digest))
    x = a.X
    values = x.data if sparse.issparse(x) else np.asarray(x).ravel()
    if not np.isfinite(values).all():
        raise ValueError("Expression contains nonfinite values")
    layer = cfg.get("counts_layer") or ("counts" if "counts" in a.layers else None)
    kind = cfg.get("matrix_kind", "auto")
    if kind not in ("auto", "counts", "log1p", "normalized"):
        raise ValueError("matrix_kind must be auto, counts, log1p, or normalized")
    if layer and layer not in a.layers:
        raise ValueError(f"Configured counts layer is missing: {layer}")
    candidate = a.layers[layer] if layer else a.X
    countlike = _countlike(candidate)
    if (kind == "counts" or layer) and not countlike:
        raise ValueError("Declared raw counts contain negative or noninteger values")
    has_counts = bool(countlike and (kind in ("auto", "counts") or layer))
    if has_counts:
        a.layers["counts"] = candidate.copy()
    if kind == "auto":
        kind = "counts" if _countlike(a.X) else ("log1p" if "log1p" in a.uns else "unknown")
        warnings.append("Matrix scale inferred heuristically; integer nonnegative values suggest counts but do not establish provenance.")
    a.uns["workflow_matrix_kind"] = kind
    a.uns["workflow_has_counts"] = has_counts
    inventory = dict(cells=a.n_obs, genes=a.n_vars, duplicate_genes=duplicate_genes,
                     matrix_kind=kind, counts_available=has_counts, sparse=bool(sparse.issparse(a.X)),
                     metadata_columns=list(a.obs.columns), metadata_missing=a.obs.isna().sum().to_dict(),
                     existing_obsm=list(a.obsm), existing_layers=list(a.layers), existing_graphs=list(a.obsp))
    for role in ("sample", "donor", "condition"):
        col = cfg.get(role + "_column")
        if col and col not in a.obs:
            raise ValueError(f"Configured {role} column absent: {col}")
        if col:
            inventory[role + "_counts"] = a.obs[col].value_counts(dropna=False).to_dict()
    file = out / "inspected.h5ad"
    a.write_h5ad(file)
    return _result(dict(dataset=str(file), inventory=_json(out / "inventory.json", inventory),
                        source_manifest=_json(out / "source_manifest.json", manifest)), inventory, warnings, inputs=[str(source)])


def qc(ctx):
    """Sample-specific robust lower-tail QC; upper tails flagged, never silently deleted."""
    import numpy as np
    import pandas as pd
    import scanpy as sc
    p = _previous(ctx, "inspection")
    a = sc.read_h5ad(p)
    out = _dir(ctx, "qc")
    cfg = ctx["config"]
    warnings = []
    if not a.uns["workflow_has_counts"]:
        warnings.append("No raw counts: library-size/mitochondrial QC, doublet detection and count filtering skipped.")
        a.obs["qc_pass"] = True
    else:
        x = a.layers["counts"]
        total = np.asarray(x.sum(axis=1)).ravel()
        detected = np.asarray((x > 0).sum(axis=1)).ravel()
        prefix = cfg.get("mitochondrial_prefix")
        if not prefix:
            org = str(cfg.get("organism", "")).lower()
            prefix = "MT-" if org in ("human", "homo sapiens") else "mt-" if org in ("mouse", "mus musculus") else None
        mt = a.var["source_gene_id"].str.startswith(prefix).to_numpy() if prefix else np.zeros(a.n_vars, dtype=bool)
        mt_fraction = np.asarray(x[:, mt].sum(axis=1)).ravel() / np.maximum(total, 1) if mt.any() else np.full(a.n_obs, np.nan)
        a.obs["total_counts"] = total
        a.obs["n_genes_by_counts"] = detected
        a.obs["pct_counts_mt"] = 100 * mt_fraction
        sample = cfg.get("sample_column") or cfg.get("donor_column")
        groups = a.obs[sample].astype(str) if sample else pd.Series("pooled_unknown", index=a.obs_names)
        if not sample:
            warnings.append("No sample identifier: QC uses pooled distribution, which may mask sample effects.")
        keep = (total > 0) & (detected > 0)
        thresholds = []
        upper = np.zeros(a.n_obs, dtype=bool)
        for group in groups.unique():
            ix = np.flatnonzero(groups.to_numpy() == group)
            row = dict(sample=group, n_cells=len(ix), method="median +/- 3 scaled MAD of log1p count metrics")
            for name, values in (("total_counts", total), ("n_genes_by_counts", detected)):
                v = np.log1p(values[ix]); med = np.median(v); mad = 1.4826 * np.median(abs(v-med))
                lo = max(0, np.expm1(med - 3 * mad)) if len(ix) >= 20 and mad > 0 else 0
                hi = np.expm1(med + 3 * mad) if len(ix) >= 20 and mad > 0 else float("inf")
                keep[ix] &= values[ix] >= lo
                upper[ix] |= values[ix] > hi
                row[name + "_lower"] = float(lo)
                row[name + "_upper_flag"] = float(hi) if np.isfinite(hi) else None
            if mt.any() and len(ix) >= 20:
                v = 100 * mt_fraction[ix]; med = np.median(v); mad = 1.4826 * np.median(abs(v-med))
                cutoff = med + 3 * mad if mad > 0 else None
                if cutoff is not None:
                    keep[ix] &= v <= cutoff
                row["mitochondrial_upper"] = cutoff
            thresholds.append(row)
        a.obs["qc_pass"] = keep
        a.obs["high_library_complexity_flag"] = upper
        _json(out / "thresholds.json", thresholds)
        warnings.append("Adaptive QC is exploratory; high-count cells flagged as potential doublets, not classified as doublets. Dedicated doublet and ambient RNA inference skipped without validated configuration/raw droplets.")
        if not mt.any():
            warnings.append("Mitochondrial fraction unavailable: supply organism-aware identifiers/prefix.")
    a.obs.to_csv(out / "cell_qc.csv")
    n = a.n_obs
    source_counts = a.layers["counts"].copy() if a.uns["workflow_has_counts"] else None
    original_cells = a.obs_names.copy(); original_genes = a.var_names.copy()
    a = a[a.obs["qc_pass"].to_numpy()].copy()
    if a.n_obs < 3:
        raise ValueError("QC leaves fewer than three cells; review thresholds")
    if a.uns["workflow_has_counts"]:
        a = a[:, np.asarray(a.layers["counts"].sum(axis=0)).ravel() > 0].copy()
    if a.n_vars < 3:
        raise ValueError("QC leaves fewer than three expressed genes")
    if source_counts is not None:
        reference = source_counts[original_cells.get_indexer(a.obs_names), :][:, original_genes.get_indexer(a.var_names)]
        difference = reference != a.layers["counts"]
        if bool(difference.sum()):
            raise ValueError("Raw count preservation check failed")
    condition, sample = cfg.get("condition_column"), cfg.get("sample_column")
    batch = cfg.get("batch_column")
    if condition and batch and batch in a.obs:
        tab = pd.crosstab(a.obs[batch], a.obs[condition]); tab.to_csv(out / "batch_condition.csv")
        if (tab.gt(0).sum(axis=1) == 1).all():
            warnings.append("Batch is nested in condition: batch/condition effects may be confounded.")
    file = out / "filtered.h5ad"; a.write_h5ad(file)
    outputs = dict(dataset=str(file), cell_qc=str(out / "cell_qc.csv"))
    for extra in ("thresholds.json", "batch_condition.csv"):
        if (out / extra).exists():
            outputs[Path(extra).stem] = str(out / extra)
    return _result(outputs,
                   dict(input_cells=n, retained_cells=a.n_obs, retained_genes=a.n_vars, unique_cell_ids=True,
                        retained_counts_exactly_preserved=source_counts is not None), warnings, inputs=[p])


def representation(ctx):
    import scanpy as sc
    import numpy as np
    p = _previous(ctx, "qc"); a = sc.read_h5ad(p); out = _dir(ctx, "representation")
    kind = a.uns["workflow_matrix_kind"]
    warnings = ["No automatic batch integration: biological preservation cannot be assessed without study design."]
    if kind == "counts":
        a.X = a.layers["counts"].copy()
        a.uns.pop("log1p", None)
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    elif kind == "normalized":
        sc.pp.log1p(a)
    elif kind != "log1p":
        return _result(status="blocked", warnings=["Unknown expression scale: explicitly set matrix_kind before normalization/PCA."], inputs=[p])
    a.layers["log_expression"] = a.X.copy()
    try:
        sc.pp.highly_variable_genes(a, n_top_genes=min(int(ctx["config"].get("n_top_genes",2000)), a.n_vars), flavor="seurat")
    except (ValueError, IndexError):
        a.var["highly_variable"] = True
        warnings.append("HVG selection failed on small/degenerate data; all genes used.")
    n_hvg = int(a.var["highly_variable"].sum())
    if n_hvg < 3:
        a.var["highly_variable"] = True; n_hvg = a.n_vars
    n_pc = min(int(ctx["config"].get("n_pcs",30)), a.n_obs-1, n_hvg-1)
    sc.pp.pca(a, n_comps=n_pc, mask_var="highly_variable", random_state=ctx.get("seed", 0))
    if not np.isfinite(a.obsm["X_pca"]).all():
        raise ValueError("PCA produced nonfinite values")
    file = out / "representation.h5ad"; a.write_h5ad(file)
    return _result(dict(dataset=str(file)), dict(n_pcs=n_pc, highly_variable_genes=n_hvg), warnings, inputs=[p])


def graph(ctx):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components
    from sklearn.neighbors import NearestNeighbors
    
    p = _previous(ctx, "representation"); a = sc.read_h5ad(p); out = _dir(ctx, "graph")
    k = min(int(ctx["config"].get("n_neighbors", 15)), a.n_obs-1)
    if k < 1:
        raise ValueError("n_neighbors must be positive")
    z = a.obsm["X_pca"]
    def build(n):
        ds, ix = NearestNeighbors(n_neighbors=n+1).fit(z).kneighbors(z)
        clean = [(ids[ids != row][:n], dist[ids != row][:n]) for row,(ids,dist) in enumerate(zip(ix,ds))]
        rows = np.repeat(np.arange(len(z)), n)
        cols = np.concatenate([v[0] for v in clean]); d = np.concatenate([v[1] for v in clean])
        scale = max(float(np.median(d[d > 0])) if (d > 0).any() else 1., 1e-12)
        w = sparse.csr_matrix((np.exp(-d/scale), (rows, cols)), shape=(len(z), len(z)))
        w = w.maximum(w.T); w.setdiag(0); w.eliminate_zeros()
        distances = sparse.csr_matrix((d, (rows, cols)), shape=w.shape).maximum(sparse.csr_matrix((d, (rows, cols)), shape=w.shape).T)
        return w, distances
    w, d = build(k)
    if w.shape != (a.n_obs, a.n_obs) or (w-w.T).nnz or w.diagonal().any():
        raise ValueError("Graph integrity check failed")
    spec=json.dumps(dict(k=k,metric="euclidean",representation="X_pca",cells=a.obs_names.tolist(),
                         genes=a.var_names.tolist(),n_pcs=z.shape[1],n_top_genes=ctx["config"].get("n_top_genes",2000)),sort_keys=True)
    version = hashlib.sha256(z.tobytes() + spec.encode()).hexdigest()[:16]
    a.obsp["connectivities"] = w; a.obsp["distances"] = d
    a.uns["neighbors"] = dict(connectivities_key="connectivities", distances_key="distances", params=dict(n_neighbors=k, method="custom_exponential", metric="euclidean", use_rep="X_pca"))
    a.uns["expression_graph_version"] = version
    diag = dict(type="expression_similarity", representation="X_pca", metric="euclidean", n_neighbors=k, version=version,
                kernel="exp(-distance / median_positive_neighbor_distance); scale=1 if all distances zero",
                symmetrization="maximum of directed neighbor weights; diagonal removed",
                components=int(connected_components(w)[0]), cells=a.n_obs, edges=int(w.nnz//2), sensitivity={})
    diag["metadata_mixing_descriptive"] = {}
    edge_rows, edge_cols = w.nonzero()
    for role in ("sample", "batch", "condition"):
        col = ctx["config"].get(role+"_column")
        if col and col in a.obs:
            v=a.obs[col].astype(str).to_numpy()
            diag["metadata_mixing_descriptive"][role] = dict(column=col,
                fraction_edges_same_group=float(np.mean(v[edge_rows] == v[edge_cols])),
                random_label_baseline=float((a.obs[col].value_counts(normalize=True)**2).sum()))
    for nk in sorted(set([max(1, k//2), min(a.n_obs-1, k*2)])):
        ww, _ = build(nk); union=(w+ww).astype(bool).nnz
        diag["sensitivity"][str(nk)] = dict(components=int(connected_components(ww)[0]), edge_jaccard=float(w.astype(bool).multiply(ww.astype(bool)).nnz/max(union,1)))
    upper=sparse.triu(w).tocoo()
    pd.DataFrame(dict(source=a.obs_names[upper.row], target=a.obs_names[upper.col], weight=upper.data, graph_version=version)).to_csv(out / "edges.csv", index=False)
    sparse.save_npz(out / "expression_graph.npz", w)
    file=out / "graph.h5ad"; a.write_h5ad(file)
    return _result(dict(dataset=str(file), graph=str(out / "expression_graph.npz"), edges=str(out / "edges.csv"), diagnostics=_json(out / "diagnostics.json", diag)), diag,
                   ["Expression edges indicate similarity, not physical proximity; no spatial coordinates or interaction network inferred."], inputs=[p])


def clustering(ctx):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score, silhouette_score
    p=_previous(ctx,"graph"); a=sc.read_h5ad(p); out=_dir(ctx,"clustering")
    seed=ctx.get("seed",0); warnings=[]; rows=[]; candidates=[]
    try:
        import leidenalg
        method="leiden"
    except ImportError:
        method="kmeans"
        warnings.append("Leiden unavailable: explicit PCA k-means fallback; resolution is number of clusters, not graph modularity.")
    resolutions=ctx["config"].get("resolutions",[0.3,0.6,1.0])[:5]
    for j,r in enumerate(resolutions):
        if method=="leiden":
            sc.tl.leiden(a,resolution=float(r),key_added="candidate",random_state=seed)
            labels=a.obs["candidate"].astype(str).to_numpy()
            sc.tl.leiden(a,resolution=float(r),key_added="candidate_repeat",random_state=seed+1)
            repeat=a.obs["candidate_repeat"].astype(str).to_numpy()
        else:
            k=min(max(2,round(float(r)*5)), a.n_obs-1)
            labels=KMeans(n_clusters=k,n_init=10,random_state=seed).fit_predict(a.obsm["X_pca"]).astype(str)
            repeat=KMeans(n_clusters=k,n_init=10,random_state=seed+1).fit_predict(a.obsm["X_pca"]).astype(str)
        n=len(set(labels)); sil=float(silhouette_score(a.obsm["X_pca"],labels,sample_size=min(2000,a.n_obs),random_state=seed)) if 1<n<a.n_obs else -1.
        rows.append(dict(resolution=float(r),n_clusters=n,silhouette=sil,seed_ARI=float(adjusted_rand_score(labels,repeat))))
        candidates.append(labels)
    chosen=max(range(len(rows)),key=lambda i: rows[i]["silhouette"])
    a.obs["cluster"]=pd.Categorical(candidates[chosen]); pd.DataFrame(rows).to_csv(out / "resolution_diagnostics.csv",index=False)
    markers=pd.DataFrame(columns=["group","names","scores","logfoldchanges","pvals","pvals_adj"])
    sizes=a.obs["cluster"].value_counts()
    if len(sizes)>1 and sizes.min()>1:
        sc.tl.rank_genes_groups(a,"cluster",method="wilcoxon",use_raw=False)
        markers=sc.get.rank_genes_groups_df(a,group=None)
        warnings.append("Marker p-values are exploratory cell-level ranking statistics; they are not condition-level inference.")
    markers.to_csv(out / "markers.csv",index=False)
    a.obs["cell_type"]="unknown"; a.obs["annotation_confidence"]=0.
    marker_sets=ctx["config"].get("markers") or {}
    annotation=[]
    for c in a.obs["cluster"].cat.categories:
        tops=markers.loc[markers["group"].astype(str)==str(c)].head(50) if len(markers) else markers
        positive=set(tops.loc[tops["logfoldchanges"]>0,"names"].astype(str)) if len(tops) else set()
        evidence=[(name,sorted(set(genes)&positive),len(set(genes)&set(a.var_names))) for name,genes in marker_sets.items()]
        ranked=sorted(evidence,key=lambda t:len(t[1]),reverse=True)
        label="unknown"; confidence=0.
        if ranked and len(ranked[0][1])>=2:
            label=ranked[0][0] if len(ranked)==1 or len(ranked[0][1])>len(ranked[1][1]) else "ambiguous"
            confidence=min(.8,len(ranked[0][1])/max(ranked[0][2],1)) if label!="ambiguous" else .25
        a.obs.loc[a.obs["cluster"]==c,"cell_type"]=label
        a.obs.loc[a.obs["cluster"]==c,"annotation_confidence"]=confidence
        annotation.append(dict(cluster=str(c),label=label,confidence=confidence,evidence=evidence))
    if not marker_sets:
        warnings.append("No curated organism/tissue marker sets supplied: all biological annotations remain unknown.")
    a.obs["annotation"] = a.obs["cell_type"].copy()
    try:
        sc.tl.umap(a,random_state=seed)
    except (ValueError,TypeError,RuntimeError) as exc:
        warnings.append(f"UMAP unavailable on this dataset: {exc}; PCA embedding retained.")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    emb=a.obsm.get("X_umap",a.obsm["X_pca"])
    fig,ax=plt.subplots(figsize=(7,5))
    for c in a.obs["cluster"].cat.categories:
        mask=(a.obs["cluster"]==c).to_numpy(); label=a.obs.loc[mask,"cell_type"].iloc[0]
        ax.scatter(emb[mask,0],emb[mask,1],s=4,alpha=.7,label=f"{c}: {label}",rasterized=True)
    ax.set(xlabel="UMAP 1" if "X_umap" in a.obsm else "PC 1",ylabel="UMAP 2" if "X_umap" in a.obsm else "PC 2")
    ax.legend(bbox_to_anchor=(1.02,1),loc="upper left",markerscale=2); fig.tight_layout(); fig.savefig(out / "embedding.png",dpi=300); plt.close(fig)
    for role in ("sample","donor","condition"):
        col=ctx["config"].get(role+"_column")
        if col:
            pd.crosstab(a.obs["cluster"],a.obs[col]).to_csv(out / f"cluster_by_{role}.csv")
    file=out / "annotated.h5ad"; a.write_h5ad(file)
    outputs=dict(dataset=str(file),markers=str(out / "markers.csv"),annotations=_json(out / "annotations.json",annotation),embedding=str(out / "embedding.png"),resolution_diagnostics=str(out / "resolution_diagnostics.csv"))
    for role in ("sample","donor","condition"):
        table=out / f"cluster_by_{role}.csv"
        if table.exists():
            outputs[f"cluster_by_{role}"]=str(table)
    return _result(outputs,
                   dict(method=method,selected=rows[chosen],cells=a.n_obs,annotation_confidence_type="heuristic marker overlap, not calibrated probability"),warnings,inputs=[p])
