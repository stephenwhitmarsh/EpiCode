#!/bin/bash
JOBID=$(sbatch --array=4-11 --parsable slurm_wod_project.sh)
sbatch --dependency afterok:${JOBID} --array=0 slurm_wod_project.sh                                                                                