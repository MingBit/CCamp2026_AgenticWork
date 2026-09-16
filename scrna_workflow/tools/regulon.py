"""Candidate TF-target coexpression tool."""
from pathlib import Path
import json
from ._downstream_common import _get, _base, _result, _dataset, _json

def regulon(ctx):
    """Infer bounded, bootstrap-tested TF coexpression candidates from local resources.

    Requires a local one-column TF-symbol list with matching declared organism and
    verified counts. Rank association is computed on log library-normalized counts.
    Activity is a target-expression score, not a binding or causal estimate.
    """
    out, cfg, artifacts = _base(ctx, "regulon")
    source = _dataset(artifacts)
    tf_file = cfg.get("tf_list")
    prerequisites = []
    if not source:
        prerequisites.append("A validated upstream AnnData dataset is required.")
    if not tf_file or not Path(tf_file).is_file():
        prerequisites.append("A local, one-column TF list is required (tf_list).")
    if not cfg.get("organism") or cfg.get("tf_resource_organism") != cfg.get("organism"):
        prerequisites.append("Declare organism and matching tf_resource_organism.")
    if prerequisites:
        return _result("skipped", warnings=prerequisites, actions=prerequisites)
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse, stats
    a = ad.read_h5ad(source)
    if "counts" not in a.layers or not a.uns.get("workflow_has_counts", False):
        return _result("skipped", [source], warnings=["Verified raw counts are unavailable."], actions=["Provide verified counts for regulatory preprocessing."])
    if a.n_obs < 30:
        return _result("skipped", [source], warnings=["Fewer than 30 cells; candidate-module stability is not assessable."])
    tf_names = set(line.strip().split("\t")[0] for line in Path(tf_file).read_text().splitlines() if line.strip() and not line.startswith("#"))
    names = np.asarray(a.var_names.astype(str))
    tf_idx = np.flatnonzero(np.isin(names, list(tf_names)))
    if not len(tf_idx):
        return _result("skipped", [source, str(tf_file)], warnings=["No TF identifiers match dataset var_names; supply an explicit identifier mapping."])
    raw = sparse.csr_matrix(a.layers["counts"], dtype=float)
    totals = np.asarray(raw.sum(axis=1)).ravel()
    x = (sparse.diags(1e4 / np.maximum(totals, 1)) @ raw).tocsr()
    x.data = np.log1p(x.data)
    mean = np.asarray(x.mean(axis=0)).ravel()
    var = np.asarray(x.multiply(x).mean(axis=0)).ravel() - mean ** 2
    target_idx = np.argsort(var)[-min(int(cfg.get("regulon_max_genes", 2000)), a.n_vars):]
    rng = np.random.default_rng(_get(ctx, "seed", 0))
    sample = np.sort(rng.choice(a.n_obs, min(a.n_obs, int(cfg.get("regulon_max_cells", 2000))), replace=False))
    target = x[sample][:, target_idx].toarray()
    ranks = stats.rankdata(target, axis=0)
    ranks -= ranks.mean(axis=0)
    ranks /= np.maximum(np.linalg.norm(ranks, axis=0), 1e-12)
    edges, activities = [], {}
    max_targets = min(50, max(1, int(cfg.get("regulon_max_targets", 25))))
    for idx in tf_idx[:min(200, int(cfg.get("regulon_max_tfs", 100)))]:
        tf = x[sample, idx].toarray().ravel()
        if np.std(tf) == 0:
            continue
        r = stats.rankdata(tf)
        r -= r.mean()
        r /= max(np.linalg.norm(r), 1e-12)
        corr = r @ ranks
        eligible = np.flatnonzero((corr >= float(cfg.get("regulon_min_correlation", .3))) & (target_idx != idx))
        chosen = eligible[np.argsort(corr[eligible])[-max_targets:]]
        if len(chosen) < 3:
            continue
        boot = []
        for _ in range(10):
            b = rng.choice(len(sample), len(sample), replace=True)
            boot.append([float(stats.spearmanr(tf[b], target[b, j]).statistic) for j in chosen])
        boot = np.asarray(boot)
        for k, j in enumerate(chosen):
            edges.append(dict(tf=names[idx], target=names[target_idx[j]], spearman_r=float(corr[j]),
                              bootstrap_positive_fraction=float(np.mean(boot[:, k] > 0)),
                              evidence_type="coexpression_only", motif_support=False,
                              confidence="exploratory; same-data bootstrap; not independent validation"))
        activities[names[idx]] = np.asarray(x[:, target_idx[chosen]].mean(axis=1)).ravel()
    if not edges:
        return _result("inconclusive", [source, str(tf_file)], warnings=["No TF has at least three targets passing the candidate correlation threshold."])
    edgefile = out / "tf_target_candidates.tsv"
    pd.DataFrame(edges).to_csv(edgefile, sep="\t", index=False)
    activity = pd.DataFrame(activities, index=a.obs_names)
    activity.index.name = "cell_id"
    afile = out / "candidate_module_activity.tsv.gz"
    activity.to_csv(afile, sep="\t")
    diagnostics = []
    for tf in activity:
        rho = stats.spearmanr(activity[tf], totals).statistic
        diagnostics.append({"tf": tf, "library_size_spearman": None if not np.isfinite(rho) else float(rho)})
    groups = [c for c in ["annotation", "cluster", cfg.get("sample_column"), cfg.get("donor_column"), cfg.get("condition_column")] if c and c in a.obs]
    summary_files = []
    for number, column in enumerate(dict.fromkeys(groups)):
        f = out / f"activity_group_{number}.tsv"
        activity.groupby(a.obs[column], observed=True).mean().to_csv(f, sep="\t")
        summary_files.append({"group_column": column, "path": str(f)})
    info = _json(out / "diagnostics.json", {"depth_associations": diagnostics, "summaries": summary_files,
                 "bootstrap_replicates": 10, "cell_cycle_adjusted": False, "batch_adjusted": False,
                 "sampling_unit": "cells; bootstrap is internal consistency only"})
    return _result(inputs=[source, str(tf_file)], outputs={"edges": str(edgefile), "activity": str(afile), "diagnostics": info,
                   **{f"summary_{i}": s["path"] for i, s in enumerate(summary_files)}},
                   metrics={"candidate_tfs": len(activity.columns), "edges": len(edges), "cells_used_for_inference": len(sample)},
                   warnings=["Coexpression-only candidates: no motif or causal support; cell identity, batch, and cell cycle can confound associations.", "Activity uses the same data as inference and is not independent evidence of a cell state."])
