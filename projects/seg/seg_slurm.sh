#!/bin/bash
#SBATCH --job-name=Seg
#SBATCH --partition=normal
#SBATCH --time=99:99:99
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/output-%j.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/error-%j.txt

module load MATLAB
matlab -nodesktop -softwareopengl -nosplash -nodisplay < /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epilepsy/seg/seg_cluster.m
sleep 5;
