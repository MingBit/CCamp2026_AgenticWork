"""Filter ribosomal protein genes from the analysis representation."""
import pandas as pd


def remove_ribosomal_genes(adata, symbol_column=None):
    """Return a gene-filtered copy and the removed-feature table.

    The caller retains its QC artifact and immutable input. Both X and layers in
    this analysis copy are subset together, so count/feature alignment is kept.
    """
    if symbol_column and symbol_column not in adata.var:
        raise ValueError(f"Gene symbol column is missing: {symbol_column}")
    if symbol_column is None:
        symbol_column = next(
            (name for name in ("gene_symbols", "gene_symbol", "gene_name") if name in adata.var),
            None,
        )
    symbols = (
        adata.var[symbol_column].astype("string")
        if symbol_column else pd.Series(adata.var_names, index=adata.var_names, dtype="string")
    )
    names = symbols.fillna("").str.upper()
    ribosomal = names.str.startswith(("RPL", "RPS")).fillna(False).to_numpy(dtype=bool)
    removed = pd.DataFrame({
        "feature_id": adata.var_names[ribosomal].astype(str),
        "gene_symbol": symbols.iloc[ribosomal].astype(str).to_numpy(),
    })
    filtered = adata[:, ~ribosomal].copy()
    if filtered.n_vars < 3:
        raise ValueError("RPL/RPS filtering leaves fewer than three genes")
    filtered.uns["workflow_gene_filter"] = {
        "rule": "case-insensitive gene symbol prefix RPL or RPS",
        "symbol_source": symbol_column or "var_names",
        "removed_genes": len(removed),
        "input_genes": adata.n_vars,
    }
    return filtered, removed
