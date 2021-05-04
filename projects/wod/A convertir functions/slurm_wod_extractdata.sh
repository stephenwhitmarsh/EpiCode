#!/bin/bash

#SBATCH --job-name=wod-freqdata

#SBATCH --partition=normal

#SBATCH --time=99:99:99

#SBATCH --mem=20G

#SBATCH --cpus-per-task=1

#SBATCH --chdir=.

#SBATCH --output=/network/lustre/iss01/charpier/analyses/wod/Antoine/slurm-output/output-%j_%a-%x.txt

#SBATCH --error=/network/lustre/iss01/charpier/analyses/wod/Antoine/slurm-output/error-%j_%a-%x.txt


module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "wod_extractdata('wod_setparams');"

sleep 5;

