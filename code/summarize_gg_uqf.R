rm(list=ls())
library(tidyverse)
library(glue)
library(blme)
library(reshape2)
library(brms)
library(Matrix)
library(vglmer)

source('code/aux_functions.R')

mdls_ps <- c('vglmer_weak', 'vglmer_strong',  'vglmer_collapse_none',
             'vglmer_collapse_FE', 'vglmer_collapse_main', 
             'vglmer_collapse_all')

all_dist_eigen <- all_uqf <- data.frame()
id <- as.numeric(commandArgs(trailingOnly = TRUE))
id <- expand.grid(mdl = 1:9, year = 1:2)[id,]
year <- id$year
mdl <- id$mdl
outcome <- 'binomial'
print(c(year, mdl))

m_HMC <- readRDS(glue(
  'output/output_stan/brms_y_{year}_m_{mdl}_o_{outcome}.RDS'
))

samples_HMC <- as.matrix(m_HMC$model$orig_brm)
samples_HMC <- samples_HMC[,!grepl(colnames(samples_HMC), pattern='^sd_|^var_|^lp_|^chol_|^S_[0-9]+\\[|^lprior$')]
samples_HMC <- rename_stan(samples_HMC)
covar_HMC <- cov(samples_HMC)

out_dist_eigen <- out_uqf <- data.frame()

for (type in mdls_ps){

  print(c(type, year, mdl))
  
  m_vglmer <- readRDS(glue('output/output_vglmer/{type}_y_{year}_m_{mdl}_o_{outcome}.RDS'))
  
  m_uqf <- get_uqf(model_vi = m_vglmer, 
    MCMC_covariance = covar_HMC, 
    MCMC_samples = samples_HMC, allow_unequal = TRUE)

  m_dist_eigen <- get_dist_eigen(model_vi = m_vglmer, 
    MCMC_covariance = covar_HMC, MCMC_samples = samples_HMC, allow_unequal = TRUE)

  tp <- type
  out_uqf <- bind_rows(out_uqf, m_uqf %>% mutate(factor_method = type))
  out_dist_eigen <- bind_rows(out_dist_eigen, m_dist_eigen %>% mutate(method = tp))

}

all_uqf <- bind_rows(all_uqf, 
  out_uqf %>% mutate(model = mdl, year = year, outcome = outcome)
)
all_dist_eigen <- bind_rows(all_dist_eigen, 
  out_dist_eigen %>% mutate(model = mdl, year = year, outcome = outcome)
)
path <- glue('output/output_uqf/out_{year}_{mdl}.RDS')
saveRDS(list(uqf = all_uqf, dist = all_dist_eigen), path)