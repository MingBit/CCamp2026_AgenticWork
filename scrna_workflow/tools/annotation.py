"""Evidence-based cluster annotation tool."""
import pandas as pd

def annotate_clusters(adata, markers, marker_sets, gene_symbol_column=None):
    """Annotate from positive cluster markers, matching feature IDs or gene symbols."""
    symbol_column = gene_symbol_column or next(
        (name for name in ("gene_symbols", "gene_symbol", "gene_name") if name in adata.var), None
    )
    symbols = (adata.var[symbol_column].astype(str).to_numpy()
               if symbol_column in adata.var else adata.var_names.astype(str).to_numpy())
    feature_symbols = {str(feature).upper(): str(symbol).upper()
                       for feature, symbol in zip(adata.var_names, symbols)}
    available = set(feature_symbols.values())
    adata.obs["cell_type"] = "unknown"
    adata.obs["annotation_confidence"] = 0.0
    adata.obs["annotation_source"] = "marker_set_overlap"
    annotations = []
    for cluster in adata.obs["cluster"].cat.categories:
        top = markers.loc[markers["group"].astype(str) == str(cluster)].head(50)
        positive = {feature_symbols.get(str(feature).upper(), str(feature).upper())
                    for feature in top.loc[top["logfoldchanges"] > 0, "names"]}
        evidence = [(name, sorted({gene.upper() for gene in genes} & positive),
                     len({gene.upper() for gene in genes} & available))
                    for name, genes in marker_sets.items()]
        ranked = sorted(evidence, key=lambda item: len(item[1]), reverse=True)
        label, confidence = "unknown", 0.0
        if ranked and len(ranked[0][1]) >= 2:
            tied = len(ranked) > 1 and len(ranked[0][1]) == len(ranked[1][1])
            label = "ambiguous" if tied else ranked[0][0]
            confidence = 0.25 if tied else min(0.8, len(ranked[0][1]) / max(ranked[0][2], 1))
        mask = adata.obs["cluster"] == cluster
        adata.obs.loc[mask, "cell_type"] = label
        adata.obs.loc[mask, "annotation_confidence"] = confidence
        annotations.append({"cluster": str(cluster), "label": label, "confidence": confidence, "evidence": evidence})
    adata.obs["annotation"] = adata.obs["cell_type"].copy()
    return annotations
