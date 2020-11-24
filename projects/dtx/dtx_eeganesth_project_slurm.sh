#!/bin/bash
#SBATCH --job-name=batch
#SBATCH --partition=normal,bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j_%a-%x.txt
#SBATCH --mail-user=paul.baudin@icm-institute.org
#SBATCH --mail-type=ALL

module load MATLAB/R2019b
matlab -nodesktop -softwareopengl -nosplash -nodisplay < dtx_project_anesth_eeg.m
sleep 5;
