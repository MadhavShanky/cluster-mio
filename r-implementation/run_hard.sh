#!/bin/bash
#SBATCH --job-name=mio_hard
#SBATCH --partition=mit_normal
#SBATCH --array=1-24%24
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=/home/dengly/mio-sim/logs/hard_%a.out
#SBATCH --error=/home/dengly/mio-sim/logs/hard_%a.err
# Harder-regime robustness study: one array task per cell (cid = SLURM_ARRAY_TASK_ID).
module load miniforge/25.11.0-0 gcc/12.2.0
eval "$(conda shell.bash hook)"
conda activate mio_r
unset SIM_QUICK
export MIO_DIR=/home/dengly/mio-sim
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
cid=${SLURM_ARRAY_TASK_ID}
export SIM_CELLS="${cid}-${cid}"
cd "$MIO_DIR"
echo "task ${SLURM_ARRAY_TASK_ID}: cid=${cid}"
Rscript sim/run_hard.R
