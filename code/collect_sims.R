rm(list=ls())
library(dplyr)
library(stringr)

file <- dir('output/output_simulation')
all_output <- lapply(file, FUN=function(i){readRDS(paste0('output/output_simulation/', i))})
names(all_output) <- file

out_uqf <- lapply(all_output, FUN=function(i){i$uqf}) %>% bind_rows(.id = 'file') %>% 
	mutate(family = str_extract(file, pattern='binomial|linear'))

out_time <- lapply(all_output, FUN=function(i){i$time}) %>% bind_rows(.id = 'file') %>% 
	mutate(family = str_extract(file, pattern='binomial|linear'))

out_ELBO <- lapply(all_output, FUN=function(i){i$ELBO}) %>% bind_rows(.id = 'file') %>% 
	mutate(family = str_extract(file, pattern='binomial|linear'))

out_eigen <- lapply(all_output, FUN=function(i){i$eigen}) %>% bind_rows(.id = 'file') %>% 
	mutate(family = str_extract(file, pattern='binomial|linear'))

out_acc <- lapply(all_output, FUN=function(i){i$acc}) %>% bind_rows(.id = 'file') %>% 
	mutate(family = str_extract(file, pattern='binomial|linear'))


saveRDS(list(uqf = out_uqf, 
	time = out_time, acc = out_acc,
	ELBO = out_ELBO, eigen = out_eigen), 'output/all_simulations.RDS')
