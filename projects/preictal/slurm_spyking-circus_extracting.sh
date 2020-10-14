#!/bin/bash
#SBATCH --job-name=spyking-circus-extracting
#SBATCH --partition=normal,bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=120G
#SBATCH --cpus-per-task=28
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/slurm_output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/slurm_output/error-%j_%a-%x.txt
#SBATCH --array=1

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "preictal_launch_SpykingCircus_extracting($SLURM_ARRAY_TASK_ID);"

sleep 5;

