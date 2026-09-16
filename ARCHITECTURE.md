# Agent interfaces and evidence ownership

| Specialist | Owner | Validated input | Exclusive outputs |
|---|---|---|---|
| Inspection | Orchestrator/QC | Local h5ad or 10x, optional metadata | Source snapshot, inventory, aligned dataset |
| QC | QC specialist | Inspection dataset | QC metrics/thresholds, filtered counts |
| Representation | Member 1 | Filtered dataset | Log representation, selected genes, PCA |
| Cell graph | Member 2 | PCA and cell IDs | Versioned cell neighbors, edges, diagnostics |
| Gene graph | Graph specialist | QC-filtered normalized expression | Gene neighborhoods, edges, diagnostics |
| Graph interpretation | Graph specialist | PCA, gene graph, markers, candidate regulons | Cross-evidence interpretation |
| Clustering | Member 1 | Validated graph | Labels, marker evidence, annotation confidence |
| Regulon | Member 3 | Counts + compatible local TF resource | Candidate TF–target edges, cell activities |
| Discovery | Integration | Annotated counts + design | Sample summaries, pseudobulk, eligible DE |
| Validation | Independent critique | All preceding task results | Integrity checks, issues and review requests |
| Reporting | Report specialist | Validation + evidence files | Scientific report, figures |

```mermaid
flowchart LR
  input[Read-only input] --> inspection --> qc
  qc --> representation --> graph[Cell graph] --> clustering
  qc --> gene_graph[Gene graph]
  clustering --> regulon
  clustering --> discovery
  representation --> graph_interpretation[Compare PCA, markers, regulon]
  gene_graph --> graph_interpretation
  regulon --> graph_interpretation
  clustering --> graph_interpretation
  regulon --> validation
  discovery --> validation
  qc --> validation
  graph --> validation
  gene_graph --> validation
  validation --> report
  validation -. review issue; max 2 revisions .-> orchestrator
  regulon -. candidate state evidence; review only .-> clustering
  orchestrator <--> workspace[(run_state + manifest + JSONL + artifacts)]
```

Runtime specialists are deterministic Python functions, not autonomous remote language-model calls. The orchestrator schedules independent functions concurrently. Human/LLM review can use the compact result JSON without loading count matrices into a prompt. Any future LLM integration must remain local or receive explicit authorization before sharing private data.

A task is eligible only after its parents return; failures block downstream numerical work. Critique/reporting still run to explain failures. Checkpoints bind code, configuration, software and input hashes; output hashes are verified. Distinct task subdirectories prohibit concurrent file ownership. Expression neighborhoods, TF–target networks and any future ligand–receptor/spatial networks have separate schemas and storage.

Review feedback does not silently modify annotation or validate its own discovery evidence. Automatic scientific revision is disabled (zero cycles); targeted revisions require a changed assumption in a new run, and the issue ledger specifies a maximum of two. There is no autonomous selection of biological hypotheses from unsupported inputs.
