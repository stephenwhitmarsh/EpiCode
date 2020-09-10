#!/bin/bash

#SBATCH --job-name=dtx_spikes

#SBATCH --partition=normal,bigmem

#SBATCH --time=99:99:99

#SBATCH --mem=120G

#SBATCH --cpus-per-task=2

#SBATCH --chdir=.

#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j_%a-%x.txt

#SBATCH --array=0

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "wod_project($SLURM_ARRAY_TASK_ID);"

sleep 5;

