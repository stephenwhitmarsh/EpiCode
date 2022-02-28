#!/bin/bash
#SBATCH --job-name=preictal_project
#SBATCH --partition=bigmem,normal
#SBATCH --time=99:99:99
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/vn_peictal/scripts/slurm_output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/vn_peictal/scripts/slurm_output/error-%j_%a-%x.txt
#SBATCH --array=1-3


module load MATLAB/R2020b
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "slurm_preictal_project($SLURM_ARRAY_TASK_ID);"
sleep 5;
