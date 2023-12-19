#!/bin/bash

run_time=23:00:00             # set time to run
cores=4
# begin loop 

for yr in 1 2; do
	for outcome in "binomial"; do
	   echo ${yr}-${outcome};
	   sbatch --job-name=gold_${yr}_${outcome} \
	    -o slurm_files/stan_gg_${outcome}_${yr}_%a_out.txt \
	    -e slurm_files/stan_gg_${outcome}_${yr}_%a_err.txt \
	    -t ${run_time} \
	    --mem 24G \
	    --cpus-per-task=${cores} \
	    --cluster=smp \
	    code/slurm_gg_stan.slurm ${yr} ${outcome}
	    sleep 2
	done
done

# reportseff 11952444 11952445 | more