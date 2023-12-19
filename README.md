# Replication Data For "Partially Factorized Variational Inference for High-Dimensional Mixed Models" 

Version 1.0; 19 December 2023

Authors: Max Goplerud, Omiros Papaspiliopoulos, Giacomo Zanella

## Overview

This `README.md` explains the replication archive for the analyses conducted in our paper. To use PF-VI for your own data, please skip to the "Software" section for installation instructions.

## File Listing

The files contained in this archive are listed below. Their usage is described in more detail in the following sections. To simply re-create all of the figures in the paper and examine the output of the analyses, please run the `code/create_figures.R` script.

- `build_gg`: The replicator should add the raw data from the `data.zip` from Ghitza and Gelman's (2013) replication archive (`https://doi.org/10.7910/DVN/PZAOO6`) to this folder. This is required to run the `prepare_GG.R` script that creates the `data/gg_data.RDS` object that is needed for all scripts involving the Ghitza-Gelman data.
- `code`: This contains the scripts needed to replicate the analyses.
  - `aux_functions.R`: Auxiliary functions needed in the simulations, e.g. calculating the accuracy of the variational approximation.
  - `bash_condaenv.sh`: Create conda environment to call Python from R.
  - `bash_sims.sh`: Send the simulations to the SLURM scheduler; see "Simulations" section for more information.
  - `bash_stan_gg.sh`: Send the HMC models in STAN for the Ghitza-Gelman analysis to the SLURM scheduler; see "Ghitza-Gelman Analysis" section for more information.
  - `create_figures.R`: This will create all tables and figures used in the manuscript. It relies on saved versions of the simulation and analysis output that were run on a high-performance computing environment.
  - `collect_sims.R`: Pulls together the simulation results into a single `.RDS` file.
  - `gen_data_PSFZ.py`: Python script used to sample from the collapsed Gibbs sampler from Papaspiliopoulos, Stumpf-Fetizon, and Zanella (2023) in the simulations. It is lightly adapted from `https://github.com/timsf/crossed-effects`
  - `prepare_GG.R`: This will create the relevant data needed for the Ghitza-Gelman analysis. This script directly taken from an author's previous work (Goplerud; `https://doi.org/10.7910/DVN/DI19IB`)
  - `run_simulations.R`: This script will run the simulations; see "Simulation" section for more information about its usage.
  - `run_stan_gg.R`: This script will estimate the HMC models for Ghitza-Gelman; see "Ghitza-Gelman Analysis" section for more information about its usage.
  - `run_vglmer_gg.R`: This script will estimate the VI models for the Ghitza-Gelman analysis; see "Ghitza-Gelman Analysis" section for information about its usage.
  - `slurm_.*_.slurm`: Five SLURM files (`slurm_gg_mrp`, `slurm_gg_stan`, `slurm_gg_uqf`, `slurm_gg_vglmer`, `slurm_simulations`) are used to schedule the numerical experiments. See the respective sections for information on their usage.
  - `summarize_gg_mrp.R`: Estimates the post-stratified quantities and the linear predictor for each type of respondent; see "Ghitza-Gelman Analysis" section for more information about its usage. 
  - `summarize_gg_uqf.R`: Estimates the UQF from the Ghitza-Gelman CAVI models; see "Ghitza-Gelman Analysis" section for more information about its usage.
- `crossed_effects_python`: This contains Python code to implement the samplers outlined in Papaspiliopoulos, Stumpf-Fetizon, and Zanella (2023); it is downloaded from here `https://github.com/timsf/crossed-effects`
- `data`: This contains data needed for the analyses. It contains a number of files.
  - `formula_list.RDS`: A list of formula for the Ghitza-Gelman models.
  - `gg_data.RDS`, `gg_poststrat.RDS` and `gg_statevalues.RDS`: These are needed for the Ghitza-Gelman analyses; they are not provided with the replication data but can be created as discussed above using the `prepare_GG.R` script.
  - `gg_table_description.xlsx`: This contains the text used in Table 1
- `figures`: This contains the tables and figures in the manuscript.
  - `gg_table.tex` can only be produced if `output_vglmer` contains all of the estimated models. Otherwise, the GitHub will use the pre-produced version.
- `output`: This contains the output of the simulations and the Ghitza-Gelman analysis.
  - `all_simulations.RDS`: This contains the results for all of the simulations.
  - `output_mrp`: This contains the state-level post-stratified estimates and the linear predictor estimates for the Ghitza-Gelman analyses.
  - `output_uqf`: This contains the estimated UQF for the Ghitza-Gelman analyses.
  - `output_stan`: This contains the HMC models fit using Stan for the Ghitza-Gelman analyses. These are too large to include on GitHub, and thus the folder is empty, but can be produced as noted below.
  - `output_vglmer`: This contains many different variational models estimated for the Ghitza-Gelman analyses with different factorization schemes. These are too large to include on GitHub, and thus the folder is empty, but can be produced as noted below.
  - `vglmer_disagg_timing.RDS`: This contains detailed information on the estimation time needed for each variational model.
  - `vglmer_fit.RDS`: This contains estimated information on the variational models.
- `slurm_files`: This folder contains the SLURM files (error and output logs). It is empty but it is included to ensure that the SLURM scripts would run as expected.

## Software

The software used in this analysis is a branch of the `vglmer` package that is hosted on GitHub at `https://github.com/mgoplerud/vglmer/tree/collapsed`. This corresponds to the `collapsed` branch of the `vglmer` package. It can be installed in R as follows:

```
remotes::install_github('mgoplerud/vglmer', ref = 'collapsed')
```

This branch will eventually be integrated into the main `vglmer` branch and then put onto CRAN.

# Creating Tables and Figures

To create the tables and figures in the manuscript, you may run the following script that pulls together the output from the underlying simulations that has been collected into a small number of files.

```
Rscript --vanilla code/create_figures.R
```

If you are interested in running the simulations or Ghitza-Gelman data analysis directly, please see the following sections.

# Simulations

First, to run the simulations, you must create a conda environment to call the Python scripts. This can be done as follows:

```
bash code/bash_condaenv.sh
```

We use a high-performance computing environment to parallelize the simulations. We rely on the University of Pittsburgh's Center for Research Computing that uses a SLURM task scheduler. Multiple calls to the scheduler to run all simulations can be done using the following `bash` script:

```
bash code/bash_sims.sh
```

An individual simulation can be run using the following R command where the arguments indicate, respectively: the type of model to be run ("linear", i.e., Gaussian, or "binomial"), the simulation number, and the size of the model (i.e., the number of levels of the random effect---see `run_simulations.R` for exactly how this argument maps onto the number of levels).

```
Rscript --vanilla code/run_simulations.R linear 1 1
```

The output from the full simulations are individual files that are stored in `output_simulation` folder. As this contains many files, we collect them into the zipped `output_simulation.zip` for this replication archive. These individual files are put into a single object (`output/all_simulations.RDS`) by running `Rscript --vanilla code/collect_sims.R`. If you wish to run this on your own machine, please extract the contents of `output_simulation.zip` into the `output_simulation` folder and then run `Rscript --vanilla code/collect_sims.R`.

The relevant tables and figures are created as discussed in the "Creating Tables and Figures" section.

# Ghitza-Gelman Analysis

We compare our results against an MCMC baseline estimated using STAN. The models are provided with the replication code (in `output_stan`) but can be re-estimated, if desired, using the following command. This again relies on a SLURM scheduler to parallelize the process. To do this, one must create the Ghitza-Gelman data using the steps noted above to download the data and then run `Rscript --vanilla code/prepare_GG.R`

```
bash code/bash_stan_gg.sh
```

The variational models are run using the `run_vglmer_gg.R` script as follows. This is run sequentially but a SLURM script is used to simplify the workflow. The output of each model is saved in `output_vglmer`.

```
sbatch code/slurm_gg_vglmer.slurm
```

Additional quantities after model estimation are computing, including the UQF and the accuracy of state-level predictions---see the paper for more details. These can be estimated using the following two scripts:

```
# To estimate the UQF
sbatch code/slurm_gg_uqf.slurm
# To estimate the MRP quantities and compute accuracy
sbatch code/slurm_gg_mrp.slurm
```

The output is saved, respectively, in the `output_uqf` and `output_mrp` folders and processed for the tables and figures by the code in the "Creating Tables and Figures" section.
