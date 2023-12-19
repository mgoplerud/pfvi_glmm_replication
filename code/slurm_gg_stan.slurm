#!/bin/bash

#SBATCH --array=1-9
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=mhg23@pitt.edu

module purge
module load gcc/8.2.0
module load lapack/3.8.0
module load r/4.1.0

Rscript --vanilla code/run_stan_gg.R ${1} $SLURM_ARRAY_TASK_ID ${2} $SLURM_CPUS_PER_TASK