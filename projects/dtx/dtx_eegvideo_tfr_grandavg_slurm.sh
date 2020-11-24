#!/bin/bash
JOBID=$(sbatch --parsable slurm_eegvideo_tfr.sh)
sbatch --dependency afterok:${JOBID} --array=0 slurm_eegvideo_tfr.sh     