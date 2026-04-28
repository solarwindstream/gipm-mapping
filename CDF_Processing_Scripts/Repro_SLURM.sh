#!/bin/bash
#SBATCH --partition=compute
#SBATCH -J cl_reproc   # script name
#SBATCH --ntasks=1      # request 1 core(s)
#SBATCH --mem-per-cpu=20G   # request 5GB RAM
#SBATCH --array=0-95        # task ID (called using SLURM_ARRAY_TASK_ID)
#SBATCH --time=10:0:0  # request 240.0 hour(s) runtime
#SBATCH -o /data/home/apx059/CDF_Processing_Scripts/logs/%x_%A_%a.log  # set output dir for logs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=r.l.atkinson@qmul.ac.uk

echo "Starting task #${SLURM_JOB_ID}"
echo "Starting task #${SLURM_ARRAY_TASK_ID}"

module purge

module load python/3.11
echo 'Loaded Python!'
source /data/home/apx059/GIPM_new_env/bin/activate
echo 'loaded env'

set -e
echo 'whatever the fuck this does'
python CDF_Processing_Script.py