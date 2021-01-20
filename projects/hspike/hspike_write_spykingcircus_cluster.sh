#!/bin/bash
#SBATCH --job-name=hspike_write_spykingcircus_cluster
#SBATCH --partition=bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=220G
#SBATCH --cpus-per-task=1
#SBATCH --chdir=.
#SBATCH --error=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/%A_%a_%j-%x_error.txt
#SBATCH --output=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/%A_%a_%j-%x_output.txt
#SBATCH --array=2

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "hspike_write_spykingcircus_cluster($SLURM_ARRAY_TASK_ID);"

sleep 5;
