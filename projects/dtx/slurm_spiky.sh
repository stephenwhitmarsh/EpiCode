#!/bin/bash
#SBATCH --job-name=Spiky
#SBATCH --partition=bigmem,normal
#SBATCH --time=99:99:99
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j_%a-%x.txt
#SBATCH --array=1-5

module load MATLAB/R2019b
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "dtx_Spiky($SLURM_ARRAY_TASK_ID);"
sleep 5;
