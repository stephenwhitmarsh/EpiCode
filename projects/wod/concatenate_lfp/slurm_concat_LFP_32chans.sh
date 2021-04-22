#!/bin/bash

#SBATCH --job-name=concat_lfp

#SBATCH --partition=normal,bigmem

#SBATCH --time=99:99:99

#SBATCH --mem=64G

#SBATCH --cpus-per-task=2

#SBATCH --chdir=.

#SBATCH --output=/network/lustre/iss01/charpier/analyses/wod/Sofia/slurm-output/output-%j_%a-%x.txt

#SBATCH --error=/network/lustre/iss01/charpier/analyses/wod/Sofia/slurm-output/error-%j_%a-%x.txt

#SBATCH --array=3-6

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "wod_concatenateLFP($SLURM_ARRAY_TASK_ID, 'wod_setparams_32chan');"

sleep 5;

