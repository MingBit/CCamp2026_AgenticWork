#!/bin/bash
#SBATCH --job-name=pbmc10k-subset
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --output=/mnt/storageBig8/work/muratli/PBMC10k_subset/slurm-%j.out
set -euo pipefail

# Ollama must listen on 127.0.0.1:11434 on this node; the workflow calls that address.
export SCRNA_OLLAMA_TIMEOUT=1800 SCRNA_OLLAMA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SCRNA_OLLAMA_TIMEOUT=86400 SCRNA_OLLAMA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OLLAMA_MODELS=/mnt/storageBig8/work/muratli/ollama/models OLLAMA_HOST=127.0.0.1:11434 OLLAMA_KEEP_ALIVE=-1
/mnt/storageBig8/work/muratli/ollama/bin/ollama serve > /mnt/storageBig8/work/muratli/PBMC10k_subset/ollama-$SLURM_JOB_ID.log 2>&1 &
OLLAMA_PID=$!
trap 'kill $OLLAMA_PID 2>/dev/null || true' EXIT
for i in $(seq 60); do curl -sf http://127.0.0.1:11434/api/tags > /dev/null && break; sleep 2; done
time /mnt/storageBig8/work/muratli/ollama/bin/ollama run qwen2.5:7b "Reply with OK."   # load the model before annotation

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
cd /mnt/storageBig8/work/muratli/Coding_Camps/2026
/mnt/storageBig8/work/muratli/Smoke_Data/wfenv/bin/python -u -m scrna_workflow run --config configs/pbmc10k_subset.yaml --output /mnt/storageBig8/work/muratli/PBMC10k_subset/results_run4
