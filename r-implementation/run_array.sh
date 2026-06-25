#!/bin/bash
#SBATCH --job-name=mio_sim
#SBATCH --partition=mit_normal
#SBATCH --array=1-120%30
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/home/dengly/mio-sim/logs/cell_%a.out
#SBATCH --error=/home/dengly/mio-sim/logs/cell_%a.err
module load miniforge/25.11.0-0 gcc/12.2.0
eval "$(conda shell.bash hook)"
conda activate mio_r
unset SIM_QUICK
export MIO_DIR=/home/dengly/mio-sim
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export SIM_CELLS="${SLURM_ARRAY_TASK_ID}-${SLURM_ARRAY_TASK_ID}"
cd "$MIO_DIR"
Rscript sim/run_sims.R
