"""Stage 4 (scenicplus env): SCENIC+ eGRN inference via its Snakemake pipeline.

Builds the gene-expression AnnData inside this environment (the workflow's anndata
version writes files SCENIC+'s older anndata may not read), fills the upstream
Snakemake config from stage parameters, pre-stages the offline genome annotation and
chromosome sizes (replacing the BioMart/NCBI download rule), runs Snakemake and exports
eRegulon triplets and AUC scores as plain TSV files for the workflow environment.

SCENIC+ v1.0a2 crashes (`No objects to concatenate`) when no extended eRegulon passes
its filters, which happens on small inputs. When `allow_direct_only` is set, that one
failure is recovered by building only the direct-annotation targets; results are then
direct-only and a warning is recorded. Any other failure still fails the stage.
"""
import shutil
import subprocess
import sys
from collections import deque
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

EMPTY_EGRN_ERROR = "No objects to concatenate"
DIRECT_ONLY_TARGETS = ["eRegulon_direct.tsv", "AUCell_direct.h5mu"]
EREGULON_COLUMNS = {
    "TF": "tf", "Gene": "target", "Region": "region", "eRegulon_name": "eregulon",
    "Gene_signature_name": "gene_signature", "Region_signature_name": "region_signature",
    "is_extended": "is_extended", "regulation": "regulation",
    "importance_TF2G": "importance_tf2g", "rho_TF2G": "rho_tf2g",
    "importance_x_abs_rho_TF2G": "importance_x_abs_rho_tf2g",
    "importance_R2G": "importance_r2g", "rho_R2G": "rho_r2g",
    "importance_x_abs_rho_R2G": "importance_x_abs_rho_r2g", "triplet_rank": "triplet_rank",
}


def build_gex_anndata(p, out):
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    from scipy import sparse

    counts = sparse.load_npz(p["counts_npz"]).tocsr().astype("float32")
    cells = pd.read_table(p["cells_tsv"], dtype=str, keep_default_na=False).set_index("barcode")
    genes = pd.read_table(p["genes_tsv"], dtype=str, keep_default_na=False).set_index("gene")
    if counts.shape != (len(cells), len(genes)):
        raise ValueError(f"Count matrix {counts.shape} does not match {len(cells)} cells x {len(genes)} genes.")
    adata = ad.AnnData(X=counts, obs=cells[["cell_type"]].copy(), var=pd.DataFrame(index=genes.index))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # prepare_GEX_ACC reads .raw (all retained genes, log-normalized)
    adata.write_h5ad(out)
    return adata.n_obs, adata.n_vars


def run_snakemake(command, cwd, tail_lines=400):
    """Run Snakemake, streaming its output to the stage log; returns (exit code, last lines)."""
    print("Running:", " ".join(command), flush=True)
    tail = deque(maxlen=tail_lines)
    with subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, errors="replace", bufsize=1) as process:
        for line in process.stdout:
            sys.stdout.write(line)
            tail.append(line.rstrip("\n"))
        sys.stdout.flush()
    return process.returncode, list(tail)


def empty_extended_egrn_failure(lines):
    """True only if the single failing Snakemake job was eGRN_extended with SCENIC+'s empty-result error."""
    failed_rules = [line for line in lines if line.startswith("Error in rule ")]
    return failed_rules == ["Error in rule eGRN_extended:"] and any(EMPTY_EGRN_ERROR in line for line in lines)


def export_results(snakemake_dir, stage, kinds=("direct", "extended")):
    import mudata
    import pandas as pd

    outputs, metrics = {}, {}
    files = {"direct": ("eRegulon_direct.tsv", "AUCell_direct.h5mu"), "extended": ("eRegulons_extended.tsv", "AUCell_extended.h5mu")}
    for kind in kinds:
        eregulon_file, auc_file = files[kind]
        table = pd.read_table(snakemake_dir / eregulon_file)
        missing = {"TF", "Gene", "Region", "eRegulon_name"} - set(table.columns)
        if missing:
            raise RuntimeError(f"{eregulon_file} lacks columns {sorted(missing)}")
        edges = table[[c for c in EREGULON_COLUMNS if c in table.columns]].rename(columns=EREGULON_COLUMNS)
        edges["evidence_type"] = f"scenicplus_eregulon_{kind}"
        edges_file = stage / f"eregulon_triplets_{kind}.tsv"
        edges.to_csv(edges_file, sep="\t", index=False)
        outputs[f"triplets_{kind}"] = str(edges_file)
        auc = mudata.read(str(snakemake_dir / auc_file))
        for modality in ("Gene_based", "Region_based"):
            frame = auc[modality].to_df()
            frame.index = [common.barcode_from_cell_name(c) for c in frame.index]
            frame.index.name = "cell_id"
            f = stage / f"auc_{modality.lower()}_{kind}.tsv.gz"
            frame.to_csv(f, sep="\t")
            outputs[f"auc_{modality.lower()}_{kind}"] = str(f)
        metrics[kind] = {"eregulons": int(edges["eregulon"].nunique()), "tfs": int(edges["tf"].nunique()),
                         "target_genes": int(edges["target"].nunique()), "regions": int(edges["region"].nunique()),
                         "triplets": int(len(edges)), "cells_scored": int(auc.n_obs)}
    if (snakemake_dir / "scplusmdata.h5mu").is_file():
        outputs["scplus_mudata"] = str(snakemake_dir / "scplusmdata.h5mu")
    return outputs, metrics


def run(p):
    from importlib.resources import files

    import yaml

    stage = Path(p["stage_dir"])
    gex = stage / "gex.h5ad"
    n_cells, n_genes = build_gex_anndata(p, gex)

    # Equivalent to `scenicplus init_snakemake`, but idempotent so an interrupted run resumes.
    snakemake_dir = stage / "Snakemake"
    (snakemake_dir / "config").mkdir(parents=True, exist_ok=True)
    (snakemake_dir / "workflow").mkdir(parents=True, exist_ok=True)
    package = files("scenicplus.snakemake")
    shutil.copyfile(str(package.joinpath("Snakefile")), snakemake_dir / "workflow" / "Snakefile")
    template = yaml.safe_load(Path(str(package.joinpath("config.yaml"))).read_text())
    config = common.apply_overrides(template, common.snakemake_overrides({**p, "gex_anndata": str(gex)}))
    (snakemake_dir / "config" / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))

    # Offline annotation: outputs of rule download_genome_annotations are provided up front.
    for key, source in (("genome_annotation", p["genome_annotation"]), ("chromsizes", p["chromsizes"])):
        target = snakemake_dir / config["output_data"][key]
        if not target.exists() or target.read_bytes() != Path(source).read_bytes():
            shutil.copy2(source, target)

    snakemake = shutil.which("snakemake")
    if snakemake is None:
        raise RuntimeError("snakemake is not on PATH in the SCENIC+ environment.")
    command = [snakemake, "--cores", str(int(p["n_cpu"])), "--rerun-triggers", "mtime", "--rerun-incomplete",
               "--snakefile", "workflow/Snakefile"]
    returncode, lines = run_snakemake(command, snakemake_dir)
    kinds, warnings = ("direct", "extended"), []
    if returncode != 0:
        if not (p.get("allow_direct_only") and empty_extended_egrn_failure(lines)):
            raise RuntimeError(f"Snakemake failed (exit {returncode}); see the log above.")
        print("No extended eRegulons passed the SCENIC+ filters; building direct-annotation results only.", flush=True)
        returncode, _ = run_snakemake(command + DIRECT_ONLY_TARGETS, snakemake_dir)
        if returncode != 0:
            raise RuntimeError(f"Snakemake failed building direct-only targets (exit {returncode}).")
        kinds = ("direct",)
        warnings.append("No extended (orthology-annotated) eRegulon passed the SCENIC+ filters "
                        f"(min_target_genes={p['min_target_genes']}); only direct-annotation eRegulons were built "
                        "and the combined scplusmdata.h5mu was not created.")

    outputs, metrics = export_results(snakemake_dir, stage, kinds)
    outputs.update(gex_anndata=str(gex), snakemake_dir=str(snakemake_dir))
    if metrics["direct"]["eregulons"] == 0:
        warnings.append("SCENIC+ produced no direct-annotation eRegulons.")
    return {"outputs": outputs,
            "metrics": {"gex_cells": n_cells, "gex_genes": n_genes, "extended_eregulons_built": "extended" in kinds, **metrics},
            "warnings": warnings,
            "versions": common.package_versions(["scenicplus", "pycisTopic", "pycistarget", "snakemake", "anndata", "mudata"])}


if __name__ == "__main__":
    common.stage_main(run)
