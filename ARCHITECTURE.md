# Agent interfaces and evidence ownership

| Specialist | Owner | Validated input | Exclusive outputs |
|---|---|---|---|
| Inspection | Orchestrator/QC | Local h5ad or 10x, optional metadata | Source snapshot, inventory, aligned dataset |
| QC | QC specialist | Inspection dataset | QC metrics/thresholds, filtered counts |
| Representation | Member 1 | Filtered dataset | Log representation, selected genes, PCA |
| Graph | Member 2 | PCA and cell IDs | Versioned neighbors, edges, diagnostics |
| Clustering | Member 1 | Validated graph | Labels, marker evidence, annotation confidence |
| Regulon | Member 3 | Counts + compatible local TF resource; for SCENIC+ also cell-type labels, ATAC fragments, motif databases (separate conda env) | Candidate TF–target edges or SCENIC+ eRegulon triplets, cell activities |
| Discovery | Integration | Annotated counts + design | Sample summaries, pseudobulk, eligible DE |
| Validation | Independent critique | All preceding task results | Integrity checks, issues and review requests |
| Reporting | Report specialist | Validation + evidence files | Scientific report, figures |

Modality routing occurs before inspection. RNA follows QC, PCA, the cell graph, and clustering. ATAC branches to TF-IDF/LSI, protein branches to CLR normalization, and multi-modal inputs request optional WNN/MultiVI integration. Unsupported or unavailable branches are recorded as explicit skips.

```mermaid
flowchart LR
  input[Read-only input] --> inspection --> qc --> representation --> graph --> clustering
  clustering --> regulon
  clustering --> discovery
  regulon --> validation
  discovery --> validation
  qc --> validation
  graph --> validation
  validation --> report
  validation -. review issue; max 2 revisions .-> orchestrator
  regulon -. candidate state evidence; review only .-> clustering
  orchestrator <--> workspace[(run_state + manifest + JSONL + artifacts)]
```

Runtime specialists are deterministic Python functions, not autonomous remote language-model calls. The orchestrator schedules independent functions concurrently. Human/LLM review can use the compact result JSON without loading count matrices into a prompt. Any future LLM integration must remain local or receive explicit authorization before sharing private data.

A task is eligible only after its parents return; failures block downstream numerical work. Critique/reporting still run to explain failures. Checkpoints bind code, configuration, software and input hashes; output hashes are verified. Distinct task subdirectories prohibit concurrent file ownership. Expression neighborhoods, TF–target networks and any future ligand–receptor/spatial networks have separate schemas and storage.

Review feedback does not silently modify annotation or validate its own discovery evidence. Automatic scientific revision is disabled (zero cycles); targeted revisions require a changed assumption in a new run, and the issue ledger specifies a maximum of two. There is no autonomous selection of biological hypotheses from unsupported inputs.
