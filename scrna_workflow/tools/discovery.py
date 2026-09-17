"""Sample summaries and guarded donor-level pseudobulk DE tools."""
from pathlib import Path
import json
from ._downstream_common import _get, _base, _result, _dataset, _json

from .pseudobulk_de import _differential_expression


def discovery(ctx):
    """Export sample-level descriptive abundance and count-preserving pseudobulk.

    Independent-donor DE runs only with an explicit contrast, adequate replication,
    and PyDESeq2. Paired or repeated-donor designs are deliberately deferred.
    """
    out, cfg, artifacts = _base(ctx, "discovery")
    source = _dataset(artifacts)
    if not source:
        return _result("skipped", warnings=["No processed dataset is available."])
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    a = ad.read_h5ad(source)
    sample = cfg.get("sample_column")
    condition = cfg.get("condition_column")
    donor = cfg.get("donor_column")
    provisional = ("annotation_source" in a.obs and
                   a.obs["annotation_source"].astype(str).eq("provisional_llm_no_reference").all())
    label = ("annotation" if not provisional and "annotation" in a.obs
             and (a.obs["annotation"].astype(str) != "unknown").any() else "cluster")
    if label not in a.obs or not sample or sample not in a.obs:
        return _result("skipped", [source], warnings=["Sample-level discovery requires sample_column and clustering labels."], actions=["Supply validated biological sample identifiers; do not substitute individual cells."])
    if a.obs[sample].isna().any():
        return _result("skipped", [source], warnings=["Missing sample identifiers prevent defensible sample-level aggregation."])
    obs = a.obs.copy()
    warnings = []
    if provisional:
        warnings.append("Reference-free LLM labels are provisional; discovery groups cells by cluster IDs.")
    columns = list(dict.fromkeys(c for c in (condition, donor, cfg.get("batch_column")) if c and c in obs))
    metadata = obs.groupby(sample, observed=True)[columns].first() if columns else pd.DataFrame(index=pd.Index(obs[sample].unique(), name=sample))
    inconsistent = [c for c in columns if (obs.groupby(sample, observed=True)[c].nunique(dropna=False) > 1).any()]
    if inconsistent:
        warnings.append("Sample metadata are inconsistent within sample: " + ", ".join(inconsistent))
        metadata = metadata.drop(columns=inconsistent)
    table = pd.crosstab(obs[sample], obs[label])
    proportions = table.div(table.sum(axis=1), axis=0)
    abundance = out / "sample_abundance.tsv"
    proportions.to_csv(abundance, sep="\t")
    nfile = out / "sample_cell_counts.tsv"
    table.to_csv(nfile, sep="\t")
    mfile = out / "sample_metadata.tsv"
    metadata.to_csv(mfile, sep="\t")
    outputs = {"abundance": str(abundance), "cell_counts": str(nfile), "sample_metadata": str(mfile)}
    has_counts = "counts" in a.layers and bool(a.uns.get("workflow_has_counts", False))
    if has_counts:
        groups = obs.groupby([sample, label], observed=True, sort=True).indices
        rows, records = [], []
        for number, ((s, celltype), indices) in enumerate(groups.items()):
            rows.append(sparse.csr_matrix(a.layers["counts"])[indices].sum(axis=0))
            records.append({"pseudobulk_id": f"pb_{number}", "sample": str(s), "population": str(celltype), "n_cells": len(indices)})
        matrix = sparse.csr_matrix(np.vstack([np.asarray(r).ravel() for r in rows]))
        pb = ad.AnnData(matrix, obs=pd.DataFrame(records).set_index("pseudobulk_id"), var=a.var.copy())
        for c in metadata.columns:
            pb.obs[c] = pb.obs["sample"].map({str(k): v for k, v in metadata[c].items()})
        f = out / "pseudobulk_counts.h5ad"
        pb.write_h5ad(f)
        outputs["pseudobulk"] = str(f)
        de_outputs, de_info = _differential_expression(pb, cfg, out)
        outputs.update(de_outputs)
        outputs["de_diagnostics"] = _json(out / "de_diagnostics.json", de_info)
        warnings.extend("DE population " + str(p["population"]) + ": " + warning
                        for p in de_info.get("populations", []) for warning in p.get("model_warnings", []))
    else:
        warnings.append("Pseudobulk counts skipped: verified raw count layer unavailable.")
    prerequisites = {"differential_expression": "Reviewed primary contrast, independent biological replicates, pairing/covariates and a count-based DE backend are required. No cell-level condition tests were run.",
                     "condition_abundance_tests": "Confirm independent donor/sample units and paired/longitudinal design before inference.",
                     "pathway_comparison": "Compatible local gene-set resource and reviewed donor-aware design required.",
                     "ligand_receptor": "Organism-compatible local ligand-receptor resource and replicated design required.",
                     "trajectory": "A biological trajectory question and supported root/direction evidence are required."}
    if "differential_expression" in outputs:
        prerequisites.pop("differential_expression")
    elif has_counts:
        prerequisites["differential_expression"] = de_info.get("reason", "Review population-specific DE diagnostics.")
    outputs["skipped_analyses"] = _json(out / "skipped_analyses.json", prerequisites)
    warnings.append("Abundances are descriptive; no abundance p-values or causal mechanisms are inferred.")
    return _result(inputs=[source], outputs=outputs, metrics={"samples": len(table), "populations": len(table.columns), "population_column": label}, warnings=warnings, actions=list(prerequisites.values()))
