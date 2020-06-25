#!/bin/bash
#SBATCH --job-name=Spiky
#SBATCH --partition=bigmem,normal
#SBATCH --time=99:99:99
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j_%a.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j_%a.txt

module load MATLAB
matlab -nodesktop -softwareopengl -nosplash -nodisplay < /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/test_Spiky.m
sleep 5;
