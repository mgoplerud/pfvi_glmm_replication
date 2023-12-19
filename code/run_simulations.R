rm(list=ls())
library(reticulate)
library(tidyverse)
library(RSpectra)
library(rstan)
library(glue)
library(vglmer)
library(Matrix)
library(caret)
library(lme4)

source('code/aux_functions.R')
args <- commandArgs(trailingOnly = T)
fmly <- as.character(args[1])
sim <- as.numeric(args[2])
v <- as.numeric(args[3])
v <- seq(5, 15)[v]

use_condaenv('python_modern', required = TRUE)
source_python('code/gen_data_PSFZ.py')

message('Simulation Information')
message(paste(c(fmly, sim, v), collapse = ' '))

offset <- 0
seed <- as.integer(offset + sim * 100 + v * 10^3)
n_samples <- 10^5
n_warmup <- 10^4

time_MCMC <- Sys.time()
print(time_MCMC)
if (fmly == 'linear'){
	sim_data <- generate_PSFZ_data(
		as.integer(2^v), as.integer(seed),
		as.integer(n_samples), as.integer(n_warmup), 'gaussian')
}else{
	sim_data <- generate_PSFZ_data(
		as.integer(2^v), as.integer(seed),
		as.integer(n_samples), as.integer(n_warmup), 'binomial')
}

print(time_MCMC - Sys.time())
PSFZ_samples <- sim_data[[2]]
PSFZ_samples <- t(sapply(PSFZ_samples, FUN=function(i){unlist(i[[1]])}))
colnames(PSFZ_samples) <- c('(Intercept)', 
	paste0('X1 @ (Intercept) @ ', 1:2^v),
	paste0('X2 @ (Intercept) @ ', 1:2^v)
)

print(dim(PSFZ_samples))
sim_data <- sim_data[[1]]
if (fmly == 'linear'){
	names(sim_data) <- c('y', 'junk1', 'junk2', 'levels', 'Z')
	vglmer_data <- data.frame(y = sim_data[[1]], sim_data[[5]] + 1)
}else{
	names(sim_data) <- c('y', 'junk1', 'levels', 'Z')
	vglmer_data <- data.frame(y = sim_data[[1]], sim_data[[4]] + 1)
}

vglmer_data$X1 <- factor(vglmer_data$X1, levels = 1:2^v)
vglmer_data$X2 <- factor(vglmer_data$X2, levels = 1:2^v)
levels_min <- pmin(length(unique(vglmer_data$X1)), length(unique(vglmer_data$X2)))

N <- length(sim_data$y)
sim_data$Z <- sim_data$Z + 1

data_stan <- list(
	N = N, Y = sim_data$y, 
	N_1 = sim_data$levels[1], N_2 = sim_data$levels[2],
	J_1 = sim_data$Z[,1], J_2 = sim_data$Z[,2],
	nu_1 = 1, nu_2 = 1,
	prior_S_1 = 1/2, prior_S_2 = 1/2,
	Z_1_1 = rep(1, N), Z_2_1 = rep(1, N),
	prior_only = FALSE
)

if (fmly == 'binomial'){
	data_stan$trials <- rep(1, N)
}

gc()


param_names_format <- c('(Intercept)',
	paste0('X1 @ (Intercept) @ ', 1:2^v),
	paste0('X2 @ (Intercept) @ ', 1:2^v))

# Get covariance from PSFZ collapsed Gibbs sampler 
PSFZ_covar <- cov(PSFZ_samples)

fmt_PSFZ <- data.frame(
	name = colnames(PSFZ_samples),
	mean = colMeans(PSFZ_samples),
	var = apply(PSFZ_samples, MARGIN = 2, var)
)
rownames(fmt_PSFZ) <- NULL

###########
# Variational Analyses
###########

model_configurations <- list(
  'lme4' = NA,
  'CAVI' = vglmer_control(
    drop.unused.levels = FALSE,
    factorization_method = 'strong',
    do_timing = TRUE, quiet_rho = FALSE,
    linpred_method = 'cyclical',
    do_SQUAREM = FALSE,
    tolerance_elbo = 1e-06, tolerance_parameters = 0,
    parameter_expansion = 'none',
    return_data = TRUE,
    prior_variance = 'mean_exists'),
  'Collapse [FE]' = vglmer_control(
    drop.unused.levels = FALSE,
    factorization_method = 'partially_factorized', collapse_set = 'FE',
    do_timing = TRUE, quiet_rho = FALSE,
    linpred_method = 'cyclical',
    tolerance_elbo = 1e-06, tolerance_parameters = 0,
    do_SQUAREM = FALSE,
    parameter_expansion = 'none',
    prior_variance = 'mean_exists'),
  'Collapse [0]' = vglmer_control(
    drop.unused.levels = FALSE,
    factorization_method = 'partially_factorized', collapse_set = 0,
    do_timing = TRUE, quiet_rho = FALSE,
    linpred_method = 'cyclical',
    tolerance_elbo = 1e-06, tolerance_parameters = 0,
    do_SQUAREM = FALSE,
    parameter_expansion = 'none',
    prior_variance = 'mean_exists'),
  'Collapse [Inf]' = vglmer_control(
    drop.unused.levels = FALSE,
    factorization_method = 'partially_factorized', collapse_set = Inf,
    do_timing = TRUE, quiet_rho = FALSE,
    linpred_method = 'cyclical',
    tolerance_elbo = 1e-06, tolerance_parameters = 0,
    do_SQUAREM = FALSE,
    parameter_expansion = 'none',
    prior_variance = 'mean_exists'),
  'CAVI Unfactorized' = vglmer_control(
    drop.unused.levels = FALSE,
    factorization_method = 'weak',
    do_timing = TRUE, quiet_rho = FALSE,
    linpred_method = 'cyclical',
    tolerance_elbo = 1e-06, tolerance_parameters = 0,
    do_SQUAREM = FALSE,
    parameter_expansion = 'none',
    return_data = TRUE,
    prior_variance = 'mean_exists')
)

out_ELBO <- out_time <- out_acc <- out_uqf <- out_eigen <- list()
for (m in names(model_configurations)){
  print(m)
  m_config <- model_configurations[[m]]
  if (m == 'lme4'){
    time_start <- Sys.time()
    if (fmly == 'linear'){
      fit_lme4 <- suppressWarnings(
          lmer(y ~ (1 | X1) + (1 | X2), 
          data = vglmer_data)
      )
    }else{
      fit_lme4 <- suppressWarnings(
          glmer(y ~ (1 | X1) + (1 | X2), 
          data = vglmer_data)
      )
    }
    time_end <- Sys.time() - time_start
    print(time_end)

    out_time[[m]] <- data.frame(
      time = as.double(time_end, units = 'mins')
    )

    rm(fit_lme4); gc()
    next
  }else{
    m_config$return_data <- TRUE
    time_start <- Sys.time()
    fit_vglmer <- suppressWarnings(
        vglmer(y ~ (1 | X1) + (1 | X2), 
        control = m_config,
        data = vglmer_data, family = fmly)
    )
    time_end <- Sys.time() - time_start
    print(time_end)

    out_ELBO[[m]] <- data.frame(ELBO = ELBO(fit_vglmer, 'trajectory')) %>% 
      mutate(it = 1:n())    
    out_time[[m]] <- data.frame(time = as.double(time_end, units = 'mins'),
      time_CAVI_only = fit_vglmer$internal_parameters$CAVI_time / 60,
      iterations = fit_vglmer$internal_parameters$it_used)
  }

  out_dist_eigen <- get_dist_eigen(model_vi = fit_vglmer, 
    MCMC_covariance = PSFZ_covar, 
    MCMC_samples = PSFZ_samples)

  out_summary <- data.frame(
    get_uqf(model_vi = fit_vglmer, MCMC_covariance = PSFZ_covar, MCMC_samples = PSFZ_samples)
  ) %>% mutate(method = 'PSFZ')

  out_marginal <- data.frame(
    get_marginal(model_vi = fit_vglmer, MCMC_covariance = PSFZ_covar, MCMC_samples = PSFZ_samples)
  )

  out_summary <- cbind(out_summary, out_marginal)

	if (m_config$factorization_method == 'partially_factorized'){
		chol_pre <- Cholesky(drop0(forceSymmetric(
			vglmer:::extract_precision(fit_vglmer)
		)))
		chol_pre <- expand(chol_pre)
		z <- matrix(rnorm(nrow(chol_pre$L) * 4000), ncol = 4000)
		rawsamples <- t(chol_pre$P) %*% solve(t(chol_pre$L), z)
		rawsamples <- t(rawsamples) 
		rm(z, chol_pre); gc()
		rawsamples <- sweep(rawsamples, MARGIN = 2, STATS = c(fit_vglmer$beta$mean, fit_vglmer$alpha$mean), FUN = '+')
		colnames(rawsamples) <- c(rownames(fit_vglmer$beta$mean), rownames(fit_vglmer$alpha$mean))
	}else{
		#Get samples from HMC and glmer
		#aligned with vglmer; use unexported ::: function to do so
		#simply puts them in the same order as the existing columns
		rawsamples <- predict(fit_vglmer,
      newdata = vglmer_data[1:5,], samples = 4000, 
      samples_only = T, allow_missing_levels = FALSE)
		colnames(rawsamples) <- c(rownames(fit_vglmer$beta$mean), rownames(fit_vglmer$alpha$mean))
	}

  fmt_VI <- format_vglmer(fit_vglmer)

  acc_PSFZ <- accuracy_bayesian(rawsamples, PSFZ_samples, colnames(rawsamples))
  out_summary$acc <- c(min(acc_PSFZ$accuracy))
  out_summary$acc_RE <- c(min(acc_PSFZ$accuracy[-1]))
  out_summary$acc_FE <- c(min(acc_PSFZ$accuracy[1]))

	sum_PSFZ_marg <- left_join(fmt_VI, 
		fmt_PSFZ %>% rename(mcmc_mean = mean, mcmc_var = var), by = 'name') %>% 
		summarize(min_ratio = min(var/mcmc_var), 
			bias_mean = mean( (mean - mcmc_mean)/abs(mcmc_mean) ) * 100, 
			bias_var = mean( (var - mcmc_var)/mcmc_var ) * 100)

  out_summary <- cbind(out_summary, sum_PSFZ_marg)    

  out_acc[[m]] <- acc_PSFZ
	out_summary <- out_summary %>% relocate(method)
  out_uqf[[m]] <- out_summary
  out_eigen[[m]] <- out_dist_eigen
  rm(time_start, time_end, fit_vglmer, rawsamples); gc()
}

out_acc <- bind_rows(out_acc, .id = 'model')
out_eigen <- bind_rows(out_eigen, .id = 'model')
out_ELBO <- bind_rows(out_ELBO, .id = 'model')
out_time <- bind_rows(out_time, .id = 'model')
out_uqf <- bind_rows(out_uqf, .id = 'model')

out_acc <- out_acc %>% mutate(levels = v, sim = sim, seed = seed)
out_eigen <- out_eigen %>% mutate(levels = v, sim = sim, seed = seed)
out_ELBO <- out_ELBO %>% mutate(levels = v, sim = sim, seed = seed)
out_time <- out_time %>% mutate(levels = v, sim = sim, seed = seed)
out_uqf <- out_uqf %>% mutate(levels = v, sim = sim, seed = seed)
out_uqf$levels_min <- levels_min

path <- glue('output/output_simulation/out_{sim}_{v}_{fmly}.RDS')
saveRDS(list(ELBO = out_ELBO, time = out_time, 
              acc = out_acc,
              uqf = out_uqf, eigen = out_eigen), path)