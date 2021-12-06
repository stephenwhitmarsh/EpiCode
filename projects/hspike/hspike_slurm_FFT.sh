#!/bin/bash
#SBATCH --job-name=hspike_FFT
#SBATCH --partition=normal
#SBATCH --time=99:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --chdir=.
#SBATCH --error=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm/hspike_FFT_%A_%a_%j-%x_error.txt
#SBATCH --output=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm/hspike_FFT_%A_%a_%j-%x_output.txt
#SBATCH --array=1-8

module load MATLAB/R2020b
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "hspike_cluster_FFT($SLURM_ARRAY_TASK_ID);"
sleep 5;
