#!/bin/bash
#SBATCH --job-name=spyking-circus
#SBATCH --partition=normal,bigmem
#SBATCH --time=99:99:99
#SBATCH --mem=120G
#SBATCH --cpus-per-task=28
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/slurm_output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/slurm_output/error-%j_%a-%x.txt


module load spyking-circus/0.9.9

cd /network/lustre/iss01/charpier/analyses/vn_preictal/data/pat_02256_0700_Crise2_m2mCi/p1
spyking-circus pat_02256_0700_Crise2_m2mCi-p1-multifile-m2mCi_2.ncs -c 28
spyking-circus pat_02256_0700_Crise2_m2mCi-p1-multifile-m2mCi_2.ncs -m converting -c 28

cd /network/lustre/iss01/charpier/analyses/vn_preictal/data/pat_02379_0828_Crise1_mHa2d/p1
spyking-circus pat_02379_0828_Crise1_mHa2d-p1-multifile-mHa2d_2.ncs -c 28
spyking-circus pat_02379_0828_Crise1_mHa2d-p1-multifile-mHa2d_2.ncs -m converting -c 28

cd /network/lustre/iss01/charpier/analyses/vn_preictal/data/pat_02379_0828_Crise1_mHaBg/p1
spyking-circus pat_02379_0828_Crise1_mHaBg-p1-multifile-mHaBg_6.ncs -c 28
spyking-circus pat_02379_0828_Crise1_mHaBg-p1-multifile-mHaBg_6.ncs -m converting -c 28

cd /network/lustre/iss01/charpier/analyses/vn_preictal/data/pat_02379_0828_Crise2_mHa2d/p1
spyking-circus pat_02379_0828_Crise2_mHa2d-p1-multifile-mHa2d_2.ncs -c 28
spyking-circus pat_02379_0828_Crise2_mHa2d-p1-multifile-mHa2d_2.ncs -m converting -c 28

cd /network/lustre/iss01/charpier/analyses/vn_preictal/data/pat_02379_0828_Crise2_mHaBg/p1
spyking-circus pat_02379_0828_Crise2_mHaBg-p1-multifile-mHaBg_4.ncs -c 28
spyking-circus pat_02379_0828_Crise2_mHaBg-p1-multifile-mHaBg_4.ncs -m converting -c 28

cd /network/lustre/iss01/charpier/analyses/vn_preictal/data/pat_02599_1057_Crise1_mHaT2/p1
spyking-circus pat_02599_1057_Crise1_mHaT2-p1-multifile-mHaT2_6.ncs -c 28
spyking-circus pat_02599_1057_Crise1_mHaT2-p1-multifile-mHaT2_6.ncs -m converting -c 28

cd /network/lustre/iss01/charpier/analyses/vn_preictal/data/pat_02599_1057_Crise2_mHaT2/p1
spyking-circus pat_02599_1057_Crise2_mHaT2-p1-multifile-mHaT2_6.ncs -c 28
spyking-circus pat_02599_1057_Crise2_mHaT2-p1-multifile-mHaT2_6.ncs -m converting -c 28

cd /network/lustre/iss01/charpier/analyses/vn_preictal/data/pat_02599_1057_Crise3_mHaT2/p1
spyking-circus pat_02599_1057_Crise3_mHaT2-p1-multifile-mHaT2_6.ncs -c 28
spyking-circus pat_02599_1057_Crise3_mHaT2-p1-multifile-mHaT2_6.ncs -m converting -c 28

sleep 5;