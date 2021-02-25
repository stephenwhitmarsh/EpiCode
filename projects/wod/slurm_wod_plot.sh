#!/bin/bash

#SBATCH --job-name=wod_plot

#SBATCH --partition=normal

#SBATCH --time=99:99:99

#SBATCH --mem=40G

#SBATCH --cpus-per-task=2

#SBATCH --chdir=.

#SBATCH --output=/network/lustre/iss01/charpier/analyses/wod/slurm-output/output-%j_-%x.txt

#SBATCH --error=/network/lustre/iss01/charpier/analyses/wod/slurm-output/error-%j_-%x.txt


module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "wod_plot;"

sleep 5;

