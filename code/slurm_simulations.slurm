#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cluster=smp
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=mhg23@pitt.edu

module purge
module load gcc/8.2.0 
module load lapack/3.8.0
module load r/4.1.0
module load python/anaconda3.9-2021.11

echo "FAMILY :" ${1}
echo "SIM    :" ${SLURM_ARRAY_TASK_ID}
echo "SIZE   :" ${3}
Rscript --vanilla code/run_simulations.R ${1} ${SLURM_ARRAY_TASK_ID} ${2}