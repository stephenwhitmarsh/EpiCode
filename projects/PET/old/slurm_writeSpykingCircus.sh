#!/bin/bash
#SBATCH --job-name=write_sc
#SBATCH --partition=bigmem,normal
#SBATCH --time=99:99:99
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/slurm_output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/slurm_output/error-%j_%a-%x.txt
#SBATCH --array=1


module load MATLAB/R2019b
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "preictal_writeSpykingCircus($SLURM_ARRAY_TASK_ID);"
sleep 5;
