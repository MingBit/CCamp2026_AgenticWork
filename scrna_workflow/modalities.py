"""Local modality detection and optional multimodal preprocessing specialists."""
from pathlib import Path
import json

from .core import _dir, _json, _result


def _text(value):
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return str(value)


def detect_modalities(input_path):
    """Return modality names detected in an h5, h5ad, or h5mu input."""
    path = Path(input_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    if path.is_dir() and (path / "cell_feature_matrix.h5").exists() and (path / "cells.csv.gz").exists():
        return ["rna"]
    suffix = path.suffix.lower()
    if suffix not in {".h5", ".hdf5", ".h5ad", ".h5mu"}:
        raise ValueError("Modality detection supports .h5, .h5ad, and .h5mu files.")

    if suffix == ".h5ad":
        import anndata as ad
        data = ad.read_h5ad(path, backed="r")
        labels = set()
        for key in ("modality", "modality_type"):
            if key in data.uns:
                labels.add(_text(data.uns[key]).lower())
        if "feature_types" in data.var:
            labels.update(_text(value).lower() for value in data.var["feature_types"].dropna().unique())
        if "modality" in data.var:
            labels.update(_text(value).lower() for value in data.var["modality"].dropna().unique())
        modalities = set()
        for label in labels:
            if "atac" in label or "peak" in label or "chromatin" in label:
                modalities.add("atac")
            elif "protein" in label or "antibody" in label or "adt" in label:
                modalities.add("protein")
            elif "rna" in label or "gene expression" in label or "gene" in label:
                modalities.add("rna")
        if not modalities:
            modalities.add("rna")
        return sorted(modalities)

    import h5py
    with h5py.File(path, "r") as handle:
        names = []
        def visit(name, obj):
            parts = name.lower().split("/")
            if len(parts) > 1 and parts[0] == "mod" and parts[1] in {"rna", "atac", "protein"}:
                names.append(parts[1])
            if hasattr(obj, "shape") and name.lower().endswith(("feature_type", "feature_types")):
                names.extend(_text(value).lower() for value in obj[()])
        handle.visititems(visit)
        keys = [name.lower() for name in handle.keys()]
    modalities = set()
    for label in names + keys:
        if "atac" in label or "peak" in label or "chromatin" in label:
            modalities.add("atac")
        elif "protein" in label or "antibody" in label or "adt" in label:
            modalities.add("protein")
        elif "rna" in label or "gene expression" in label:
            modalities.add("rna")
    if suffix == ".h5mu":
        for label in keys:
            if label in {"rna", "atac", "protein"}:
                modalities.add(label)
    if not modalities:
        modalities.add("rna")
    return sorted(modalities)


def detect_modalities_task(ctx):
    source = ctx["config"].get("input_path")
    if not source:
        return _result(status="blocked", warnings=["No input_path was supplied for modality detection."])
    modalities = detect_modalities(source)
    out = _dir(ctx, "inspection")
    result_path = out / "modalities.json"
    _json(result_path, {"input_path": str(Path(source).resolve()), "modalities": modalities})
    return _result(outputs={"modalities": str(result_path)}, metrics={"modalities": modalities}, inputs=[str(source)])


def _get_modalities(ctx):
    item = ctx.get("artifacts", {}).get("modality_detection", {})
    output = item.get("outputs", {}).get("modalities")
    if output and Path(output).exists():
        return json.loads(Path(output).read_text()).get("modalities", [])
    return []


def _read_input(ctx):
    source = Path(ctx["config"]["input_path"]).resolve()
    if source.suffix.lower() == ".h5mu":
        try:
            import mudata as md
        except ImportError as exc:
            raise RuntimeError(".h5mu input requires the optional mudata package.") from exc
        data = md.read_h5mu(source)
        return data
    import scanpy as sc
    if source.suffix.lower() == ".h5ad":
        return sc.read_h5ad(source)
    return sc.read_10x_h5(source)


def _select_modality(data, modality):
    """Subset a loaded AnnData object when feature-type metadata is available."""
    for column in ("feature_types", "modality"):
        if column in data.var:
            values = data.var[column].astype(str).str.lower()
            if modality == "atac":
                mask = values.str.contains("atac|peak|chromatin", regex=True)
            else:
                mask = values.str.contains("protein|antibody|adt", regex=True)
            if mask.any():
                return data[:, mask.to_numpy()].copy()
    return data


def atac_lsi(ctx):
    """Run TF-IDF and truncated-SVD LSI when ATAC is present."""
    import numpy as np
    from scipy import sparse
    from sklearn.decomposition import TruncatedSVD
    modalities = _get_modalities(ctx)
    if "atac" not in modalities:
        return _result(status="skipped", warnings=["ATAC modality is absent; TF-IDF/LSI was not run."])
    data = _read_input(ctx)
    if hasattr(data, "mod"):
        data = data.mod["atac"]
    else:
        data = _select_modality(data, "atac")
    matrix = sparse.csr_matrix(data.X, dtype=float)
    tf = matrix.multiply(1 / np.maximum(np.asarray(matrix.sum(axis=1)).ravel(), 1)[:, None])
    document_frequency = np.asarray((matrix > 0).sum(axis=0)).ravel()
    idf = np.log1p(matrix.shape[0] / np.maximum(document_frequency, 1))
    tfidf = tf.multiply(idf)
    n_components = min(int(ctx["config"].get("lsi_components", 30)), tfidf.shape[0] - 1, tfidf.shape[1] - 1)
    if n_components < 2:
        return _result(status="inconclusive", warnings=["ATAC matrix is too small for LSI."])
    embedding = TruncatedSVD(n_components=n_components, random_state=ctx.get("seed", 0)).fit_transform(tfidf)
    out = _dir(ctx, "atac")
    sparse.save_npz(out / "tfidf.npz", tfidf.tocsr())
    import anndata as ad
    result_data = ad.AnnData(embedding, obs=data.obs.copy())
    result_data.obsm["X_lsi"] = embedding
    result_data.uns["lsi"] = {"n_components": n_components, "method": "TF-IDF followed by TruncatedSVD"}
    dataset = out / "lsi.h5ad"
    result_data.write_h5ad(dataset)
    return _result(outputs={"dataset": str(dataset)}, metrics={"modality": "atac", "n_components": n_components}, inputs=[str(ctx["config"]["input_path"])])


def protein_clr(ctx):
    """Run centered log-ratio normalization when protein is present."""
    import numpy as np
    from scipy import sparse
    modalities = _get_modalities(ctx)
    if "protein" not in modalities:
        return _result(status="skipped", warnings=["Protein modality is absent; CLR was not run."])
    data = _read_input(ctx)
    if hasattr(data, "mod"):
        data = data.mod["protein"]
    else:
        data = _select_modality(data, "protein")
    values = data.X.toarray() if sparse.issparse(data.X) else np.asarray(data.X, dtype=float)
    log_values = np.log1p(values)
    clr = log_values - log_values.mean(axis=1, keepdims=True)
    import anndata as ad
    out = _dir(ctx, "protein")
    result_data = ad.AnnData(clr, obs=data.obs.copy(), var=data.var.copy())
    result_data.uns["clr"] = {"method": "centered log-ratio on log1p protein counts"}
    dataset = out / "clr.h5ad"
    result_data.write_h5ad(dataset)
    return _result(outputs={"dataset": str(dataset)}, metrics={"modality": "protein", "features": clr.shape[1]}, inputs=[str(ctx["config"]["input_path"])])


def joint_integration(ctx):
    """Select WNN or MultiVI for multi-modal inputs when its backend is installed."""
    modalities = _get_modalities(ctx)
    if len(modalities) <= 1:
        return _result(status="skipped", metrics={"modalities": modalities}, warnings=["Fewer than two modalities are present; joint integration was not run."])
    backend = None
    try:
        import muon  # noqa: F401
        backend = "WNN"
    except ImportError:
        try:
            import scvi  # noqa: F401
            backend = "MultiVI"
        except ImportError:
            return _result(status="skipped", metrics={"modalities": modalities}, warnings=["Multiple modalities detected, but neither muon (WNN) nor scvi-tools (MultiVI) is installed."])
    return _result(status="inconclusive", metrics={"modalities": modalities, "backend": backend}, warnings=[f"{backend} backend detected; joint integration requires modality-specific preparation and review before execution."])


def plot_multimodal_diagnostics(ctx):
    """Plot diagnostic summaries for available modalities (ATAC / Protein / Joint)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import anndata as ad
    import numpy as np
    
    out = _dir(ctx, "multimodal_plots")
    outputs = {}
    
    # Protein CLR Plotting
    protein_file = Path(ctx["root"]) / "protein" / "clr.h5ad"
    if protein_file.exists():
        pdata = ad.read_h5ad(protein_file)
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        means = np.asarray(pdata.X.mean(axis=0)).ravel()
        features = pdata.var_names[:30] # Top 30 proteins
        ax.bar(range(len(features)), means[:len(features)], color="#2b5c8f")
        ax.set_xticks(range(len(features)))
        ax.set_xticklabels(features, rotation=90, fontsize=8)
        ax.set_ylabel("Mean CLR Intensity", fontsize=10)
        ax.set_xlabel("Protein Surface Marker", fontsize=10)
        ax.set_title("Surface Protein Abundance Profile (CLR Normalized)", fontsize=12, fontweight="bold")
        p_path = out / "protein_clr_summary.png"
        fig.savefig(p_path, dpi=300)
        plt.close(fig)
        outputs["protein_summary"] = str(p_path)
        
    # ATAC LSI Plotting
    atac_file = Path(ctx["root"]) / "atac" / "lsi.h5ad"
    if atac_file.exists():
        adata = ad.read_h5ad(atac_file)
        lsi = adata.obsm["X_lsi"]
        fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        scatter = ax.scatter(lsi[:, 0], lsi[:, 1], s=4, c="teal", alpha=0.6)
        ax.set_xlabel("LSI Dimension 1", fontsize=10)
        ax.set_ylabel("LSI Dimension 2", fontsize=10)
        ax.set_title("Chromatin Accessibility Projection (ATAC LSI)", fontsize=12, fontweight="bold")
        a_path = out / "atac_lsi_summary.png"
        fig.savefig(a_path, dpi=300)
        plt.close(fig)
        outputs["atac_summary"] = str(a_path)
        
    return _result(outputs=outputs)