import json, sys, yaml
from pathlib import Path
repo, smoke, spy, ncpu = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], int(sys.argv[4])
sys.path.insert(0, str(repo))
from scrna_workflow import regulon_scenicplus as rs
from scrna_workflow.downstream import _result
from scrna_workflow.runner import PATH_KEYS
cfg = yaml.safe_load((repo / "configs/pbmc10k_multiome.yaml").read_text())
for k in PATH_KEYS:
    if cfg.get(k): cfg[k] = str((repo / "configs" / cfg[k]).resolve())
cfg.update(scenicplus_python=spy, scenicplus_min_target_genes=3, scenicplus_work_dir=str(smoke / "work"), scenicplus_temp_dir=str(smoke / "t"),
           scenicplus_fragments_path=str(smoke / "fragments.tsv.gz"), scenicplus_n_cpu=ncpu,
           scenicplus_mallet_memory_gb=8, scenicplus_keep_chromosomes=["chr21", "chr22"],
           scenicplus_n_topics=[5, 10], scenicplus_topic_n_iter=50, scenicplus_min_cells=50,
           scenicplus_qc_min_unique_fragments_in_peaks=10, scenicplus_qc_min_tss_enrichment=1,
           scenicplus_qc_use_automatic_thresholds=False, scenicplus_gsea_n_perm=100)
out = smoke / "run" / "regulon"; out.mkdir(parents=True, exist_ok=True)
print("prerequisites:", rs.prerequisites(cfg, str(smoke / "smoke.h5ad")) or "ok")
r = rs.run({"config": cfg}, out, cfg, str(smoke / "smoke.h5ad"), _result)
print(json.dumps({k: r[k] for k in ("status", "metrics", "warnings")}, indent=1, default=str))
