#!/bin/bash
#SBATCH --job-name=hspike_cluster
#SBATCH --partition=normal
#SBATCH --time=99:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=24
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm/%j_%A-%a-%x-output.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm/%j_%A-%a-%x-error.txt
#SBATCH --array=1-7

module load MATLAB/R2020b
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "hspike_spike_cluster($SLURM_ARRAY_TASK_ID);"
sleep 5;
