#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cluster=smp
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=mhg23@pitt.edu
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --array=1-18
#SBATCH --output=slurm_files/gg_%a_out.txt
#SBATCH --error=slurm_files/gg_%a_err.txt

module purge
module load gcc/8.2.0
module load r/4.1.0

# sbatch code/slurm_gg_mrp.slurm
# reportseff 11951315 | more

Rscript --vanilla code/summarize_gg_mrp.R ${SLURM_ARRAY_TASK_ID}