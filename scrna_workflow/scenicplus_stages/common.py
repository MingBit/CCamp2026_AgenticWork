"""Standard-library helpers shared by the workflow and the SCENIC+ stage scripts.

The stage scripts run inside the separate `scenicplus` conda environment (see
envs/scenicplus.yml), which cannot import the workflow package or its dependencies.
They load this module by path, so it must only use the standard library (PyYAML is
imported lazily; both environments provide it).
"""
import hashlib
import json
import re
import sys
from pathlib import Path

SPLIT_PATTERN = "___"
STAGE_OUTPUT = "stage_output.json"
_SAFE = re.compile(r"[^A-Za-z0-9_.-]+")


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    tmp.replace(path)
    return str(path)


def sha256_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def file_stamp(path):
    """Cheap identity of a large, immutable resource: resolved path, size and mtime."""
    stat = Path(path).stat()
    return [str(Path(path).resolve()), stat.st_size, stat.st_mtime_ns]


def safe_label(label):
    """Filesystem/shell-safe token for a cell-type label (MACS2 names, BED file names)."""
    token = _SAFE.sub("_", str(label)).strip("._")
    return token or "unnamed"


def unique_safe_labels(labels):
    """Map each distinct label to a unique safe token, stable under label order."""
    mapping, used = {}, set()
    for label in sorted({str(x) for x in labels}):
        token, n = safe_label(label), 2
        base = token
        while token in used:
            token, n = f"{base}_{n}", n + 1
        mapping[label] = token
        used.add(token)
    return mapping


def cell_name(barcode, sample_id):
    """pycisTopic cell name: <barcode>___<sample_id>."""
    return f"{barcode}{SPLIT_PATTERN}{sample_id}"


def barcode_from_cell_name(name):
    return str(name).split(SPLIT_PATTERN)[0]


def parse_region(name):
    """Split a "chrom:start-end" region name into (chrom, start, end)."""
    chrom, span = str(name).rsplit(":", 1)
    start, end = span.split("-")
    return chrom, int(start), int(end)


def write_region_bed(path, region_names):
    """Write region names as a 3-column BED file; returns the number of regions."""
    rows = sorted({parse_region(r) for r in region_names})
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows))
    return len(rows)


def spaced(values):
    """SCENIC+ Snakemake passes multi-value parameters as space-separated strings."""
    return " ".join(str(v) for v in values)


def snakemake_overrides(p):
    """SCENIC+ v1.0a2 Snakemake config values derived from the stage parameters.

    Keys are (section, key) of the upstream config.yaml; every key must already exist
    in the template so that a renamed upstream parameter fails loudly.
    """
    sample_id = p["sample_id"]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample_id):
        raise ValueError(f"scenicplus_sample_id must match [A-Za-z0-9_.-]+, got {sample_id!r}")
    return {
        ("input_data", "cisTopic_obj_fname"): p["cistopic_object"],
        ("input_data", "GEX_anndata_fname"): p["gex_anndata"],
        ("input_data", "region_set_folder"): p["region_set_folder"],
        ("input_data", "ctx_db_fname"): p["ctx_db"],
        ("input_data", "dem_db_fname"): p["dem_db"],
        ("input_data", "path_to_motif_annotations"): p["motif_annotations"],
        ("params_general", "temp_dir"): p["temp_dir"],
        ("params_general", "n_cpu"): int(p["n_cpu"]),
        ("params_general", "seed"): int(p["seed"]),
        # Evaluated by the SCENIC+ CLI: maps RNA barcodes onto cisTopic cell names.
        ("params_data_preparation", "bc_transform_func"): f"\"lambda x: f'{{x}}{SPLIT_PATTERN}{sample_id}'\"",
        ("params_data_preparation", "is_multiome"): bool(p["is_multiome"]),
        ("params_data_preparation", "key_to_group_by"): p.get("key_to_group_by") or "",
        ("params_data_preparation", "nr_cells_per_metacells"): int(p["nr_cells_per_metacells"]),
        ("params_data_preparation", "species"): "hsapiens" if p["species"] == "homo_sapiens" else p["species"],
        ("params_data_preparation", "search_space_upstream"): spaced(p["search_space_upstream"]),
        ("params_data_preparation", "search_space_downstream"): spaced(p["search_space_downstream"]),
        ("params_data_preparation", "search_space_extend_tss"): spaced(p["search_space_extend_tss"]),
        ("params_motif_enrichment", "species"): p["species"],
        ("params_motif_enrichment", "annotation_version"): p["motif_annotation_version"],
        ("params_motif_enrichment", "motif_similarity_fdr"): p["motif_similarity_fdr"],
        ("params_motif_enrichment", "orthologous_identity_threshold"): p["orthologous_identity_threshold"],
        ("params_motif_enrichment", "annotations_to_use"): spaced(p["annotations_to_use"]),
        ("params_motif_enrichment", "fraction_overlap_w_dem_database"): p["fraction_overlap_w_dem_database"],
        ("params_motif_enrichment", "dem_max_bg_regions"): int(p["dem_max_bg_regions"]),
        ("params_motif_enrichment", "dem_balance_number_of_promoters"): bool(p["dem_balance_number_of_promoters"]),
        ("params_motif_enrichment", "dem_promoter_space"): int(p["dem_promoter_space"]),
        ("params_motif_enrichment", "dem_adj_pval_thr"): p["dem_adj_pval_thr"],
        ("params_motif_enrichment", "dem_log2fc_thr"): p["dem_log2fc_thr"],
        ("params_motif_enrichment", "dem_mean_fg_thr"): p["dem_mean_fg_thr"],
        ("params_motif_enrichment", "dem_motif_hit_thr"): p["dem_motif_hit_thr"],
        ("params_motif_enrichment", "fraction_overlap_w_ctx_database"): p["fraction_overlap_w_ctx_database"],
        ("params_motif_enrichment", "ctx_auc_threshold"): p["ctx_auc_threshold"],
        ("params_motif_enrichment", "ctx_nes_threshold"): p["ctx_nes_threshold"],
        ("params_motif_enrichment", "ctx_rank_threshold"): p["ctx_rank_threshold"],
        ("params_inference", "tf_to_gene_importance_method"): p["tf_to_gene_importance_method"],
        ("params_inference", "region_to_gene_importance_method"): p["region_to_gene_importance_method"],
        ("params_inference", "region_to_gene_correlation_method"): p["region_to_gene_correlation_method"],
        ("params_inference", "order_regions_to_genes_by"): p["order_regions_to_genes_by"],
        ("params_inference", "order_TFs_to_genes_by"): p["order_tfs_to_genes_by"],
        ("params_inference", "gsea_n_perm"): int(p["gsea_n_perm"]),
        ("params_inference", "quantile_thresholds_region_to_gene"): spaced(p["quantile_thresholds_region_to_gene"]),
        ("params_inference", "top_n_regionTogenes_per_gene"): spaced(p["top_n_region_to_genes_per_gene"]),
        ("params_inference", "top_n_regionTogenes_per_region"): spaced(p["top_n_region_to_genes_per_region"]),
        ("params_inference", "min_regions_per_gene"): int(p["min_regions_per_gene"]),
        ("params_inference", "rho_threshold"): p["rho_threshold"],
        ("params_inference", "min_target_genes"): int(p["min_target_genes"]),
    }


def apply_overrides(template, overrides):
    """Return a copy of the Snakemake config with overrides applied to existing keys only."""
    config = json.loads(json.dumps(template))
    missing = [f"{section}.{key}" for section, key in overrides if key not in config.get(section, {})]
    if missing:
        raise KeyError("SCENIC+ Snakemake config lacks expected keys (version mismatch?): " + ", ".join(missing))
    for (section, key), value in overrides.items():
        config[section][key] = value
    return config


def stage_main(run):
    """Entry point for stage scripts: `python stageN_x.py params.json`."""
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} params.json")
    params = read_json(sys.argv[1])
    result = run(params)
    write_json(Path(params["stage_dir"]) / STAGE_OUTPUT, result)


def package_versions(names):
    from importlib import metadata
    versions = {}
    for name in names:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions
