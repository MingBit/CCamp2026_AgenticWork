"""Versioned local marker knowledge bases for grounded annotation."""
from pathlib import Path
import yaml


def load_knowledge_base(path):
    path = Path(path).expanduser().resolve()
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or not isinstance(data.get("markers"), dict):
        raise ValueError("Knowledge base must contain a marker mapping")
    markers = {}
    for label, genes in data["markers"].items():
        if not isinstance(label, str) or not label.strip() or not isinstance(genes, list):
            raise ValueError("Knowledge base labels need nonempty gene lists")
        cleaned = list(dict.fromkeys(str(gene).strip().upper() for gene in genes))
        if len(cleaned) < 2 or any(not gene for gene in cleaned):
            raise ValueError(f"Knowledge base label {label!r} needs two distinct marker genes")
        markers[label] = cleaned
    return markers, {"path": str(path), "name": data.get("name"),
                     "organism": data.get("organism"), "tissue": data.get("tissue"),
                     "source": data.get("source")}


def built_in_pbmc_knowledge(config, question):
    """Offer the bundled human PBMC panel only when species and tissue are explicit."""
    organism = str(config.get("organism") or "").lower()
    tissue = str(config.get("tissue") or "").lower()
    text = question.lower()
    human = organism in {"human", "homo sapiens"} or "human" in text
    pbmc = any(term in tissue or term in text for term in
               ("pbmc", "peripheral blood mononuclear"))
    if human and pbmc:
        return str((Path(__file__).resolve().parents[1] / "knowledge" /
                    "human_pbmc_markers.yaml").resolve())
    return None


def detected_immune_panel(adata, ranked_markers, organism=None, symbol_column=None):
    """Return the bundled panel only for a clear multi-lineage immune signature."""
    if str(organism or "").lower() not in {"", "human", "homo sapiens"}:
        return None
    column = symbol_column or next((name for name in ("gene_symbols", "gene_symbol", "gene_name")
                                    if name in adata.var), None)
    symbols = (adata.var[column].astype(str).to_numpy() if column
               else adata.var_names.astype(str).to_numpy())
    feature_symbols = {str(feature).upper(): str(symbol).upper()
                       for feature, symbol in zip(adata.var_names, symbols)}
    signatures = ({"CD3D", "CD3E", "IL7R", "CCR7"},
                  {"MS4A1", "CD79A", "CD79B"},
                  {"GNLY", "NKG7", "PRF1"},
                  {"CD14", "LYZ", "S100A8", "S100A9"})
    signature_genes = set().union(*signatures)
    matching_symbols = [str(symbol) for symbol in symbols
                        if str(symbol).upper() in signature_genes]
    if (len(matching_symbols) < 8 or
            sum(symbol.isupper() for symbol in matching_symbols) / len(matching_symbols) < 0.8):
        return None
    supported = set()
    for cluster in adata.obs["cluster"].cat.categories:
        rows = ranked_markers.loc[(ranked_markers["group"].astype(str) == str(cluster))
                                  & (ranked_markers["logfoldchanges"] > 0)].head(50)
        positive = {feature_symbols.get(str(name).upper(), str(name).upper())
                    for name in rows["names"]}
        supported.update(i for i, genes in enumerate(signatures)
                         if len(positive & genes) >= 2)
    if len(supported) >= 3:
        return str((Path(__file__).resolve().parents[1] / "knowledge" /
                    "human_pbmc_markers.yaml").resolve())
    return None
