#!/bin/bash
#SBATCH --job-name=mio_bench
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/home/dengly/mio-sim/logs/bench.out
#SBATCH --error=/home/dengly/mio-sim/logs/bench.err
# Head-to-head vs L0Learn + brute force. Single job (small grid; no sharding).
module load miniforge/25.11.0-0 gcc/12.2.0
eval "$(conda shell.bash hook)"
conda activate mio_r
unset SIM_QUICK
export MIO_DIR=/home/dengly/mio-sim
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
cd "$MIO_DIR"
Rscript sim/bench_l0learn.R
