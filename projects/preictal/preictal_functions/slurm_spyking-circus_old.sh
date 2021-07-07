#!/bin/bash
#SBATCH --job-name=spyking-circus
#SBATCH --partition=normal,bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=120G
#SBATCH --cpus-per-task=28
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/slurm_output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/slurm_output/error-%j_%a-%x.txt
#SBATCH --array=2

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "preictal_launch_SpykingCircus($SLURM_ARRAY_TASK_ID,'extracting');"

sleep 5;

