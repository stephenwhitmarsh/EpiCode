#!/bin/bash
#SBATCH --job-name=avg
#SBATCH --partition=normal,bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=50G
#SBATCH --cpus-per-task=2
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/slurm-error/error-%j_%a-%x.txt
#SBATCH --mail-user=paul.baudin@icm-institute.org
#SBATCH --mail-type=ALL

module load MATLAB/R2020b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "dtx_spikes_grandaverage;"

sleep 5;

