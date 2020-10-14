#!/bin/bash

#SBATCH --job-name=readN

#SBATCH --partition=bigmem,normal

#SBATCH --time=99:99:99

#SBATCH --mem=120G

#SBATCH --cpus-per-task=2

#SBATCH --chdir=.

#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j_%a-%x.txt

#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j_%a-%x.txt

#SBATCH --mail-user=paul.baudin@icm-institute.org

#SBATCH --mail-type=ALL

#SBATCH --array=1



module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "testreadnatus();"

sleep 5;

