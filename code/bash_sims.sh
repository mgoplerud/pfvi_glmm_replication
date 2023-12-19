#!/bin/bash

set -e

array_size=100

for fmly in "linear"; do
	for lvl in {6..6}; do
		echo "$fmly + $lvl"
		# Configure parameters
		if (($lvl < 5)); then
			# Definitely too much time but better be safe
			run_time=10:59:00
			memory=16G
		elif (( $lvl == 6 )); then
			run_time=23:59:00
			memory=16G
		else
			run_time=2-23:00:00
			memory=64G
		fi
		# Run simulations
		sbatch --job-name=sim_${fmly}_${lvl} \
		    -o slurm_files/array_${fmly}_${lvl}_%a_out.txt \
		    -e slurm_files/array_${fmly}_${lvl}_%a_err.txt \
		    -t ${run_time} --array=1-${array_size} \
		    --mem=${memory} \
		    code/slurm_simulations.slurm ${fmly} ${lvl}
		# Pause to keep scheduler happy
		sleep 5	
	done
done

# binomial
# reportseff 11848757 11848867 11848980 11849093 11849208 11849223 | more
# linear
# reportseff 11849314 11849819 11849932 11850024 11850086 11850099 | more