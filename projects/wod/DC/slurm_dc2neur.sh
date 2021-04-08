#!/bin/bash

#SBATCH --job-name=wod_dc2neur

#SBATCH --partition=normal

#SBATCH --time=99:99:99

#SBATCH --mem=20G

#SBATCH --cpus-per-task=2

#SBATCH --chdir=.

#SBATCH --output=/network/lustre/iss01/charpier/analyses/wod/slurm-output/output-%j_%a-%x.txt

#SBATCH --error=/network/lustre/iss01/charpier/analyses/wod/slurm-output/error-%j_%a-%x.txt

#SBATCH --array=1-12

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "ced_to_neuralynx($SLURM_ARRAY_TASK_ID);"

sleep 5;

