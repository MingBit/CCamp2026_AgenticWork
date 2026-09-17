#!/bin/bash
#SBATCH --job-name=scenicplus-smoke
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --output=/mnt/storageBig8/work/muratli/Smoke_Data/slurm-%j.out
set -euo pipefail

# module load java                    # only if Java comes from a module on compute nodes
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # the stages parallelise via n_cpu

cd /mnt/storageBig8/work/muratli/Coding_Camps/2026
echo "node: $(hostname)  cpus: $SLURM_CPUS_PER_TASK"
java -version 2>&1 | head -1
/mnt/storageBig8/work/muratli/Smoke_Data/wfenv/bin/python -u /mnt/storageBig8/work/muratli/Smoke_Data/smoke_regulon.py /mnt/storageBig8/work/muratli/Coding_Camps/2026 /mnt/storageBig8/work/muratli/Smoke_Data /home/muratli/.conda/envs/multiome/bin/python $SLURM_CPUS_PER_TASK
