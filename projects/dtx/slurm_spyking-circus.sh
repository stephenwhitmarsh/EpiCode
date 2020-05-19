#!/bin/bash

#SBATCH --job-name=SC_end

#SBATCH --partition=normal

#SBATCH --time=99:99:99

#SBATCH --mem=120G

#SBATCH --cpus-per-task=28

#SBATCH --chdir=.

#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j-%x.txt

#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j-%x.txt

#SBATCH --mail-user=paul.baudin@icm-institute.org

#SBATCH --mail-type=ALL



module load spyking-circus/0.9.1





cd /network/lustre/iss01/charpier/analyses/lgi1/DTX-PROBE/data/spike/DTX7/p1

spyking-circus DTX7-p1-multifile-E06.ncs -c 28

spyking-circus DTX7-p1-multifile-E06.ncs -m converting -c 28

sleep 5;



