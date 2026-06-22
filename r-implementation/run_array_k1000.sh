#!/bin/bash
#SBATCH --job-name=mio_k1000
#SBATCH --partition=mit_normal
#SBATCH --array=1-45%30
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/home/dengly/mio-sim/logs/k1000_%a.out
#SBATCH --error=/home/dengly/mio-sim/logs/k1000_%a.err
# Re-run the 9 timed-out K=1000 cells (cids 82-90) as rep-shards: 5 shards x 20 reps = 100 reps/cid total,
# pooled by scenario keys in figures.R/tables.R. taskid 1..45 -> (cid, rep-range).
module load miniforge/25.11.0-0 gcc/12.2.0
eval "$(conda shell.bash hook)"
conda activate mio_r
unset SIM_QUICK
export MIO_DIR=/home/dengly/mio-sim
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
g=$(( SLURM_ARRAY_TASK_ID - 1 ))
cid=$(( 82 + g / 5 ))
shard=$(( g % 5 ))
rstart=$(( shard * 20 + 1 ))
rend=$(( rstart + 19 ))
export SIM_CELLS="${cid}-${cid}"
export SIM_REPS_RANGE="${rstart}-${rend}"
cd "$MIO_DIR"
echo "task ${SLURM_ARRAY_TASK_ID}: cid=${cid} reps ${rstart}-${rend}"
Rscript sim/run_sims.R
