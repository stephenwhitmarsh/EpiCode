#!/bin/bash
#SBATCH --job-name=SC
#SBATCH --partition=normal,bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=60G
#SBATCH --cpus-per-task=14
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j_%a-%x.txt
#SBATCH --array=1-10

eval $(sed -n "$SLURM_ARRAY_TASK_ID"p /network/lustre/iss01/charpier/analyses/lgi1/DTX-PROBE/data/dtx_slurm_job_list.txt)
