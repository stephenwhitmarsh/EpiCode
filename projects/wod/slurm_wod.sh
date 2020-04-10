#!/bin/bash
#SBATCH --job-name=testwod
#SBATCH --partition=normal
#SBATCH --time=99:99:99
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j-%x.txt
#SBATCH --mail-user=paul.baudin@icm-institute.org
#SBATCH --mail-type=ALL

module load MATLAB/R2019b
matlab -nodesktop -softwareopengl -nosplash -nodisplay < testfunction.m
sleep 5;
