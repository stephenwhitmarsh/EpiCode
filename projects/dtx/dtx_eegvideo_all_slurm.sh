#!/bin/bash
JOBID1=$(sbatch --parsable slurm_eegvideo_cluster.sh)
JOBID2=$(sbatch --dependency afterok:${JOBID1} slurm_eegvideo_project.sh)
sbatch --dependency afterok:${JOBID2} slurm_eegvideo_plotallseizures.sh                                                                              