#!/bin/bash
JOBID=$(sbatch --array=1-15 --parsable slurm_eeganesth_tfr.sh)
sbatch --dependency afterok:${JOBID} --array=0 slurm_eeganesth_tfr.sh     