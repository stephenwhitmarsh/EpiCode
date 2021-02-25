#!/bin/bash
JOBID1=$(sbatch --parsable slurm_eeganesth_cluster.sh)
JOBID2=$(sbatch --dependency afterok:${JOBID1} slurm_eeganesth_project.sh)
