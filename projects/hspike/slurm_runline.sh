#!/bin/bash
#SBATCH --job-name=SC_runline
#SBATCH --partition=bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=128G
#SBATCH --cpus-per-task=24
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/%j_%A-%a-%x-output.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/%j_%A-%a-%x-error.txt
#SBATCH --array=7-12

eval $(sed -n "$SLURM_ARRAY_TASK_ID"p /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data/hspike/slurm_job_list.txt)

