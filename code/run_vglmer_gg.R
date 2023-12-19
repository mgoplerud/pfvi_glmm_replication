rm(list=ls())

library(brms)
library(lme4)
library(vglmer)
library(glue)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RSpectra)

draw_samples <- TRUE

source('code/aux_functions.R')

store_all <- list()

data <- readRDS('data/gg_data.RDS')
data$fake_outcome <- as.numeric(data$success > data$failure)
data$frac_success <- data$success
data$frac_failure <- data$failure

data$success <- round(data$success)
data$failure <- round(data$failure)

year_list <- c(1,2)
outcome_list <- c('binomial')

fmt_stan <- function(object, useSigma = FALSE) {
  post_stan <- as.matrix(object)
  if (!useSigma) {
    post_stan <- post_stan[, !grepl(colnames(post_stan), pattern = "^Sigma")]
  }
  if (any(grepl(colnames(post_stan), pattern='^z_[0-9]+\\[[0-9,]+\\]$'))){
    post_stan <- post_stan[, !grepl(colnames(post_stan), pattern='^z_[0-9]+\\[[0-9,]+\\]$')]
  }
  if (any(grepl(colnames(post_stan), pattern='^chol_[0-9]+\\['))){
    post_stan <- post_stan[, !grepl(colnames(post_stan), pattern='^chol_[0-9]+\\[')]
  }
  if (any(grepl(colnames(post_stan), pattern='^(S|var)_[0-9]+\\[[0-9]+(,[0-9]+)?\\]'))){
    post_stan <- post_stan[, !grepl(colnames(post_stan), pattern='^(S|var)_[0-9]+\\[[0-9]+(,[0-9]+)?\\]')]
  }
  
  post_stan <- post_stan[, !grepl(colnames(post_stan), pattern='^lp__$|^var_[0-9]+\\[[0-9]+\\]$')]
  post_stan <- post_stan[, !grepl(colnames(post_stan), pattern='^sd_')]
  colnames(post_stan) <- gsub(colnames(post_stan), pattern='^b_', replacement = '')
  
  parse_stan_names <- strsplit(x = colnames(post_stan),
                               split = '^r_|^b\\[| |\\[|\\]')
  # parse_stan_names <- strsplit(x = colnames(post_stan), split = "^b\\[| |\\]", perl = T)
  
  fmt_stan_names <- sapply(parse_stan_names, FUN = function(i) {
    if (length(i) == 1) {
      return(i)
    } else {
      i_one <- unlist(strsplit(i[3], split = ":|,"))
      if (any(grepl(i[3], pattern=':'))){
        return(paste(i_one[1], i[2], i_one[2], sep = " @ "))
      }else{
        return(paste(i[2], i_one[2], i_one[1], sep = " @ "))
      }
    }
  })
  colnames(post_stan) <- fmt_stan_names
  colnames(post_stan) <- gsub(colnames(post_stan), pattern='(?!<=\\))Intercept(?!\\))', perl = T, replacement ='(Intercept)')
  output <- data.frame(var = apply(post_stan, MARGIN = 2, var),
                       mean = colMeans(post_stan))
  output$name <- colnames(post_stan)
  rownames(output) <- NULL
  return(output)
}

rename_stan <- function(post_stan){
  colnames(post_stan) <- gsub(colnames(post_stan), pattern='^b_', replacement = '')
  
  parse_stan_names <- strsplit(x = colnames(post_stan),
                               split = '^r_|^b\\[| |\\[|\\]')
  # parse_stan_names <- strsplit(x = colnames(post_stan), split = "^b\\[| |\\]", perl = T)
  
  fmt_stan_names <- sapply(parse_stan_names, FUN = function(i) {
    if (length(i) == 1) {
      return(i)
    } else {
      i_one <- unlist(strsplit(i[3], split = ":|,"))
      if (any(grepl(i[3], pattern=':'))){
        return(paste(i_one[1], i[2], i_one[2], sep = " @ "))
      }else{
        return(paste(i[2], i_one[2], i_one[1], sep = " @ "))
      }
    }
  })
  colnames(post_stan) <- fmt_stan_names
  colnames(post_stan) <- gsub(colnames(post_stan), pattern='(?!<=\\))Intercept(?!\\))', perl = T, replacement ='(Intercept)')
  return(post_stan)
}

track_uqf <- disagg_timing <- data.frame()

colVar <- vglmer:::colVar

# Set seed for the MAVB / posterior draws to ensure replicability
# vglmer itself is deterministic given its fixed initialization

set.seed(15260)

for (y in 1:2){
  for (model in 1:9){
    for (outcome in outcome_list){
        message(glue('Model: {model}; Year: {y}; Outcome: {outcome}'))
      sub_data <- data %>% filter(year == y)
      if (nrow(sub_data) == nrow(data)){stop()}
      
      est_HMC <- readRDS(glue('output/output_stan/brms_y_{y}_m_{model}_o_{outcome}.RDS'))
      formula <- as.formula(formula(est_HMC$model$orig_brm))
      if (outcome == 'binomial'){
        formula <- update(formula, 'cbind(success, failure) ~ .')
      }
      sub_data <- sub_data %>% mutate(id = 1:n())
      sub_data <- sub_data %>% mutate(n_trials = success + failure)
      sub_data <- sub_data %>% mutate(pr = frac_success/(frac_success + frac_failure))
      
      options_for_vglmer <- expand.grid(
        factor_method = c('strong', 'weak',
          'collapse_none', 
          'collapse_all', 
          'collapse_FE',
          'collapse_main'
        ), stringsAsFactors = F)
      rownames(options_for_vglmer) <- NULL
      
      if (outcome %in% 'binomial'){
        options_for_vglmer$vi_family <- 'binomial'
      }else{
        options_for_vglmer$vi_family <- 'linear'
      }
      
      store_lp <- store_timing <- store_parameters <- data.frame()
      store_rspectra_vector <- store_ELBO <- data.frame()

      HMC_samples <- as.matrix(est_HMC$model$orig_brm)        
      HMC_samples <- HMC_samples[,!grepl(colnames(HMC_samples), pattern='^sd_|^var_|^lp_|^chol_|^S_[0-9]+\\[|^lprior$')]
      HMC_samples <- rename_stan(HMC_samples)
      
      HMC_covariance <- cov(HMC_samples)
      
      for (v in 1:nrow(options_for_vglmer)){
        options_v <- options_for_vglmer[v,,drop=F]
        print(options_v)
        rownames(options_v) <- NULL
        
        vi_family <- options_v$vi_family
        lp_method <- NULL # Use default except for "strong"/"weak"

        if (options_v$factor_method == 'collapse_FE'){
          factor_method <- 'partially_factorized'
          collapse <- 'FE'
        }else if (options_v$factor_method == 'collapse_none'){
          factor_method <- 'partially_factorized'
          collapse <- 0
        }else if (options_v$factor_method == 'collapse_all'){
          factor_method <- 'partially_factorized'
          collapse <- Inf
        }else if (options_v$factor_method == 'collapse_main'){
          factor_method <- 'partially_factorized'
          re_fmla <- sapply(lme4::findbars(formula), FUN=function(i){as.character(i)[3]})
          int_re <- re_fmla[grepl(re_fmla, pattern='\\.')]
          collapse <- unique(unlist(strsplit(int_re, split='\\.')))
          collapse <- c('FE', collapse)
        }else{
          factor_method <- options_v$factor_method
          collapse <- NA
          lp_method <- 'joint'
        }
        #Estimate vglmer
        time_vglmer <- proc.time()
        #If fractional, need to set "force_whole = TRUE"
        est_vglmer <- vglmer(formula = formula, 
            data = sub_data, family = vi_family, 
            control = vglmer_control(
              factorization_method = factor_method, 
              collapse_set = collapse,
              prior_variance = 'mean_exists', iterations = 5000,
              parameter_expansion = 'none', do_SQUAREM = FALSE,
              linpred_method = lp_method,
              tolerance_parameters = 1e-30, 
              tolerance_elbo = 1e-6,
              do_timing = TRUE, return_data = TRUE)
        )
        time_vglmer <- time_vglmer - proc.time()

        
        disagg_timing <- bind_rows(disagg_timing,
         bind_rows(est_vglmer$timing, data.frame(stage = 'Total', total = -as.numeric(time_vglmer[3]))) %>%
           mutate(year = y, formula = model, factor_method = options_v$factor_method, outcome = outcome)
        )
        
        out_vglmer <- cbind(format_vglmer(est_vglmer), options_v %>% mutate(MAVB = FALSE))

        saveRDS(est_vglmer, glue('output/output_vglmer/vglmer_{options_v$factor_method}_y_{y}_m_{model}_o_{outcome}.RDS'))
        
        param_names <- c(rownames(est_vglmer$beta$mean), rownames(est_vglmer$alpha$mean))
        stopifnot(all.equal(param_names[match(colnames(HMC_covariance), param_names)], colnames(HMC_covariance)))
        
        full_precision_VI <- vglmer:::extract_precision(est_vglmer)
        precision_VI <- full_precision_VI[
          match(colnames(HMC_covariance), param_names),
          match(colnames(HMC_covariance), param_names)
        ]
        
        product_matrix <- precision_VI %*% HMC_covariance
        
        uqf_ev <- NA
        # all_ev <- eigen(product_matrix)$values
        # uqf_ev <- 1/max(all_ev)
        # stopifnot(all(all_ev > - sqrt(.Machine$double.eps)))
        
        rspectra_fit <- eigs(product_matrix, k = 1, which = 'LM')
        rspectra_ev <- rspectra_fit$values
        stopifnot(all(Im(rspectra_ev) == 0))
        rspectra_ev <- Re(rspectra_ev)
        uqf_rspectra <- 1/max(rspectra_ev)
        
        rspectra_vector <- data.frame(vector = Re(rspectra_fit$vectors[,1]), 
            name = colnames(HMC_covariance),
            stringsAsFactors = F)
        
        if (draw_samples){
          
          if (grepl(factor_method, pattern='partially_factorized')){
            chol_pre <- Cholesky(drop0(forceSymmetric(full_precision_VI)))
            chol_pre <- expand(chol_pre)
            z <- matrix(rnorm(nrow(full_precision_VI) * 4000), ncol = 4000)
            rawsamples <- t(chol_pre$P) %*% solve(t(chol_pre$L), z)
            rawsamples <- t(rawsamples) 
            rm(z, chol_pre); gc()
            rawsamples <- sweep(rawsamples, MARGIN = 2, STATS = c(est_vglmer$beta$mean, est_vglmer$alpha$mean), FUN = '+')
            colnames(rawsamples) <- c(rownames(est_vglmer$beta$mean), rownames(est_vglmer$alpha$mean))
          }else{
            #Get samples from HMC and glmer
            #aligned with vglmer; use unexported ::: function to do so
            #simply puts them in the same order as the existing columns
            
            rawsamples <- predict(est_vglmer, newdata = sub_data,
             samples = 4000, 
             samples_only = T, 
             allow_missing_levels = TRUE)
            colnames(rawsamples) <- c(rownames(est_vglmer$beta$mean), rownames(est_vglmer$alpha$mean))
          }
          
          #stopifnot(identical(dim(HMC_samples), dim(rawsamples)))
          #stopifnot(identical(dim(HMC_samples), dim(lMAVB)))
          acc_raw <- accuracy_bayesian(HMC_samples, rawsamples, name = colnames(HMC_samples))
          out_vglmer <- full_join(out_vglmer, acc_raw, by ='name')
          
        }
        
        store_ELBO <- bind_rows(store_ELBO, cbind(est_vglmer$ELBO_trajectory, options_v))
        store_parameters <- bind_rows(store_parameters,  out_vglmer)
        store_rspectra_vector <- bind_rows(
          store_rspectra_vector, cbind(rspectra_vector, options_v)
        )
        store_timing <- bind_rows(store_timing,
          cbind(data.frame(
            iterations = nrow(est_vglmer$ELBO_trajectory),
            uqf_rspectra = uqf_rspectra, uqf_ev = uqf_ev,
            change_param = max(est_vglmer$parameter.change),
            change_mean = max(est_vglmer$parameter.change[1:2]),
            change_ELBO = diff(tail(est_vglmer$ELBO_trajectory, 2)$ELBO),
            time = -time_vglmer[3]/60, est_vglmer$ELBO[1:3]), options_v))
        
        # glmer_samples <- vglmer:::custom_glmer_samples(est_glmer$model, ordering= c(rownames(est_vglmer$beta$mean), rownames(est_vglmer$alpha$mean)), samples = 4000)
        # #stopifnot(identical(dim(HMC_samples), dim(glmer_samples)))
        # acc_glmer <- accuracy_bayesian(HMC_samples, glmer_samples, name = colnames(HMC_samples))
        # output_glmer <- full_join(output_glmer, acc_glmer, by = 'name')
        
      }

      store_parameters <- bind_rows(store_parameters ,
        fmt_stan(est_HMC$model$orig_brm) %>% mutate(factor_method = 'hmc'))
      
      store_timing <- bind_rows(store_timing, data.frame(factor_method = c('hmc', 'glmer'), time_MAVB = rep(0,2), time = c(est_HMC$time[3], NA) * -1/60, stringsAsFactors = F))
      print(store_timing %>% select(factor_method, time, time_MAVB, ELBO, iterations, change_param, uqf_ev, uqf_rspectra))
      
      track_uqf <- bind_rows(track_uqf, store_timing %>% select(factor_method, uqf_ev, uqf_rspectra) %>% mutate(model = model, family = outcome))
      
      print(
        ggplot(track_uqf %>% filter(!is.na(uqf_rspectra))) + 
          geom_point(aes(x=model,y=uqf_rspectra, col = factor_method)) +
          geom_line(aes(x=model,y=uqf_rspectra, col = factor_method)) +
          facet_wrap(~ family)
      )
      
      store_all[[paste0(y, '-', model, '-', outcome)]] <- list(
        ELBO = store_ELBO, lp = store_lp,
        formula = formula,
        rspectra_eigen = store_rspectra_vector,
        parameters = store_parameters, timing = store_timing
      )
      
      saveRDS(store_all, 'output/vglmer_fit.RDS', version = 2)
      saveRDS(disagg_timing, 'output/vglmer_disagg_timing.RDS', version = 2)
    }
  }
}

