#!/bin/bash
#SBATCH --job-name=spikes
#SBATCH --partition=bigmem,normal
#SBATCH --time=99:99:99
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j_%a-%x.txt
#SBATCH --mail-user=paul.baudin@icm-institute.org
#SBATCH --mail-type=ALL
#SBATCH --array=1-7

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "dtx_project_spikes($SLURM_ARRAY_TASK_ID);"

sleep 5;

