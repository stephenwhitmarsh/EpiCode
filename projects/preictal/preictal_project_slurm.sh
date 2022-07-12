#!/bin/bash
#SBATCH --job-name=preictal
#SBATCH --partition=normal,bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/slurm_output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/slurm_error/error-%j_%a-%x.txt
#SBATCH --array=1

module load MATLAB/R2020b
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "preictal_project($SLURM_ARRAY_TASK_ID);"
sleep 5;
