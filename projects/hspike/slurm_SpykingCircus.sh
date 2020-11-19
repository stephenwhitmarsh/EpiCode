#!/bin/bash
#SBATCH --job-name=HspikeSC
#SBATCH --partition=normal,bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=120G
#SBATCH --cpus-per-task=28
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/error-%j_%a-%x.txt
#SBATCH --array=1

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "hspike_cluster_SpykingCircus($SLURM_ARRAY_TASK_ID);"

sleep 5;

