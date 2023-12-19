#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cluster=smp
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=mhg23@pitt.edu
#SBATCH --time=23:00:00
#SBATCH --mem=8G
#SBATCH --output=slurm_files/vglmer_gg_out.txt
#SBATCH --error=slurm_files/vglmer_gg_err.txt

module purge
module load gcc/8.2.0
module load r/4.1.0

# reportseff 11947172 | more

Rscript --vanilla code/run_vglmer_gg.R