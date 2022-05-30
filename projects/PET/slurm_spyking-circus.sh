#!/bin/bash
#SBATCH --job-name=SC
#SBATCH --partition=normal,bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=120G
#SBATCH --cpus-per-task=28
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/slurm_output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/slurm_output/error-%j_%a-%x.txt
#SBATCH --array=15

eval $(sed -n "$SLURM_ARRAY_TASK_ID"p /network/lustre/iss01/charpier/analyses/vn_preictal/data/preictal_slurm_job_list.txt)