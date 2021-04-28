#!/bin/bash
#SBATCH --job-name=pnh_cluster
#SBATCH --partition=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=24
#SBATCH --time=99:00:00
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/vn_pnh/slurm/%A-%a-%x-output.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/vn_pnh/slurm/%A-%a-%x-error.txt
#SBATCH --array=1

module load MATLAB/R2020b
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "pnh_cluster($SLURM_ARRAY_TASK_ID);"
sleep 5;

###SBATCH --mem=220G
###SBATCH --cpus-per-task=24