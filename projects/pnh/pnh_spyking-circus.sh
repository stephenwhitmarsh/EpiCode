#!/bin/bash
#SBATCH --job-name=spyking-circus
#SBATCH --partition=bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=128G
#SBATCH --cpus-per-task=24
#SBATCH --output=/network/lustre/iss01/charpier/analyses/vn_pnh/slurm/%j_%A-%a-%x-output.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/vn_pnh/slurm/%j_%A-%a-%x-error.txt
#SBATCH --array=4

eval $(sed -n "$SLURM_ARRAY_TASK_ID"p /network/lustre/iss01/charpier/analyses/vn_pnh/data/pnh/slurm_job_list.txt)
