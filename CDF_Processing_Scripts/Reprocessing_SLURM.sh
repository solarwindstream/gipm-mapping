#!/bin/bash
#SBATCH -J cluster_gipm_reprocessing   # job name
#SBATCH --partition computeshort
#SBATCH -n 1      # request 1 core/task
#SBATCH --mem-per-cpu=20GB   # request 1GB RAM (can also use MB or decimal values)
#SBATCH --time=1:00:00  # request 1 hour(s) runtime - n.b. is this per task?
#SBATCH --array=0	# array task range (each task sets SLURM_ARRAY_TASK_ID var)
#SBATCH -o /gpfs/scratch/apx059/%x_%A_%a.log	# log output

echo "Starting job #${SLURM_JOBID}"
echo "Starting task #${SLURM_ARRAY_TASK_ID}"
module purge    # purge any loaded modules

module load python/3.11
source /data/home/apx059/GIPM_new_env/bin/activate

set -e
python CDF_Processing_Script.py