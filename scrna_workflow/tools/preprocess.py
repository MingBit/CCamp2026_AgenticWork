"""Normalization, feature selection, and PCA tool."""
import numpy as np
from ._core_common import _result, _load_step, _save

from .normalize import normalize_expression
from .feature_selection import select_hvg
from .pca import compute_pca
from .gene_filter import remove_ribosomal_genes


def preprocess(ctx):
    """Run normalization, HVG selection and PCA as one checkpointed stage."""
    adata, out, source = _load_step(ctx, "qc", "representation")
    config = ctx["config"]
    warnings = ["Batch integration is not applied automatically; it requires a justified study design."]
    adata, removed = remove_ribosomal_genes(adata, config.get("gene_symbol_column"))
    removed_path = out / "removed_rpl_rps_genes.tsv"
    removed.to_csv(removed_path, sep="\t", index=False)
    if removed.empty:
        warnings.append("No RPL/RPS symbols matched; verify gene identifiers if symbols are unavailable.")
    if not normalize_expression(adata):
        return _result(status="blocked", warnings=["Set matrix_kind before preprocessing an unknown expression scale."], inputs=[source])
    n_genes, warning = select_hvg(adata, config.get("n_top_genes", 2000))
    if warning:
        warnings.append(warning)
    n_pcs = compute_pca(adata, config.get("n_pcs", 30), ctx.get("seed", 0))
    return _result(
        {"dataset": _save(adata, out / "representation.h5ad"),
         "removed_rpl_rps_genes": str(removed_path)},
        {"n_pcs": n_pcs, "highly_variable_genes": n_genes,
         "removed_rpl_rps_genes": len(removed)}, warnings, inputs=[source],
    )
