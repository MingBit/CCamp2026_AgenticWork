"""Independent artifact-integrity checks and critique."""
from pathlib import Path
import json
from ._downstream_common import _get, _base, _result, _dataset, _json

def validation(ctx):
    """Check alignment and preservation; record limits rather than independent proof."""
    out, cfg, artifacts = _base(ctx, "validation")
    source = _dataset(artifacts)
    if not source:
        return _result("skipped", warnings=["No dataset is available for validation."])
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    a = ad.read_h5ad(source)
    checks = {"cell_identifiers_unique": bool(a.obs_names.is_unique), "gene_identifiers_unique": bool(a.var_names.is_unique),
              "metadata_aligned": len(a.obs) == a.n_obs and len(a.var) == a.n_vars}
    for key, matrix in a.obsp.items():
        checks["graph_dimensions_" + key] = matrix.shape == (a.n_obs, a.n_obs)
    if "counts" in a.layers and a.uns.get("workflow_has_counts", False):
        counts = a.layers["counts"]
        values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
        checks["count_shape"] = counts.shape == a.shape
        checks["counts_nonnegative_integer"] = bool(np.isfinite(values).all() and (values >= 0).all() and np.allclose(values, np.rint(values)))
    qcsource = artifacts.get("qc", {}).get("outputs", {}).get("dataset")
    if qcsource and Path(qcsource).is_file():
        q = ad.read_h5ad(qcsource)
        checks["qc_cell_order_preserved"] = bool(a.obs_names.equals(q.obs_names))
        positions = q.var_names.get_indexer(a.var_names)
        checks["qc_gene_subset_order_preserved"] = bool((positions >= 0).all() and np.all(np.diff(positions) > 0))
        gene_filter = a.uns.get("workflow_gene_filter", {})
        excluded = q.var_names[~q.var_names.isin(a.var_names)]
        checks["gene_filter_accounted_for"] = bool(
            len(excluded) == int(gene_filter.get("removed_genes", 0))
            and q.n_vars - a.n_vars == len(excluded)
        )
        symbol_source = gene_filter.get("symbol_source", "var_names")
        symbols = q.var[symbol_source] if symbol_source in q.var else q.var_names
        dropped_symbols = pd.Index(symbols)[~q.var_names.isin(a.var_names)]
        checks["excluded_genes_match_rpl_rps_rule"] = bool(
            dropped_symbols.astype(str).str.upper().str.startswith(("RPL", "RPS")).all()
        ) if len(excluded) else True
        if "counts" in a.layers and "counts" in q.layers and checks["qc_cell_order_preserved"] and checks["qc_gene_subset_order_preserved"]:
            delta = sparse.csr_matrix(a.layers["counts"]) - sparse.csr_matrix(q[:, a.var_names].layers["counts"])
            checks["counts_preserved_since_qc"] = delta.nnz == 0
    inspected = artifacts.get("inspection", {}).get("outputs", {}).get("dataset")
    if inspected and Path(inspected).is_file() and a.uns.get("workflow_has_counts", False) and "counts" in a.layers:
        original = ad.read_h5ad(inspected)
        identifiers_present = bool(a.obs_names.isin(original.obs_names).all() and a.var_names.isin(original.var_names).all())
        checks["retained_identifiers_in_inspected_input"] = identifiers_present
        if identifiers_present and "counts" in original.layers:
            retained = original[a.obs_names, a.var_names]
            delta = sparse.csr_matrix(a.layers["counts"]) - sparse.csr_matrix(retained.layers["counts"])
            checks["retained_counts_preserved_since_inspection"] = delta.nnz == 0
    issues = []
    for task, result in artifacts.items():
        for warning in result.get("warnings", []):
            issues.append({"owner": task, "issue": warning, "revision_cycles": 0, "resolution": "human review or targeted new evidence required"})
    concerns = ["These checks establish internal integrity, not independent biological validation.", "No external dataset or held-out donor validation was performed.", "Regulons defining a state cannot independently validate that state."]
    f = _json(out / "validation.json", {"checks": checks, "issues": issues, "limitations": concerns, "maximum_revision_cycles_per_issue": 2})
    return _result("completed" if all(checks.values()) else "failed", [source], {"validation": f}, {"checks": checks, "failed_checks": sum(not v for v in checks.values())}, concerns,
                   ["Review unresolved issues and use targeted donor-held-out or external validation where available."])
