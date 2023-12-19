rm(list=ls())

# Must use R 4.1.0 to run on CRC

library(tidyverse)
library(glue)
library(rstan)
library(brms)


if (!(packageVersion('brms') >= '2.18.0')){
  stop('brms must be above 2.18.0')
}

args <- commandArgs(trailingOnly=TRUE)

yr_list <- as.numeric(args[1])
mdl_list <- as.numeric(args[2])
outcome_list <- args[3]
cores <- as.numeric(args[4])

print(args)

options(mc.cores = cores)

convert_brm <- function(formula = NULL, data = NULL, brms_args = list(),
  compile = TRUE, prior_method = 'default_iw'){
  message('Compiling original BRM model')
  
  brms_args$formula <- formula
  brms_args$data <- data
  brms_args$chains <- 0
  brm_object <- do.call('brm', brms_args)
  
  if (family(brm_object)$family == 'gaussian'){
    is_gaussian <- TRUE
  }else{
    is_gaussian <- FALSE
  }
  
  data_brm <- standata(brm_object)
  
  #Get the code and split by whitespace
  code_brm <- stancode(brm_object)
  code_brm <- strsplit(code_brm, split='\\n')[[1]]
  copy_brm <<- code_brm
  message('Rewriting Code')
  re_id <- names(data_brm)[grepl(names(data_brm), pattern='^M_[0-9]+$')]
  
  for (v in re_id){
    id <- gsub(v, pattern='^M_', replacement = '')
    size_v <- data_brm[[v]]
    if (size_v == 1){
      if (prior_method == 'default_iw'){
        #Default inverse wishart (d_j + 1) and I
        #so (d_j + 1)/2 = 1
        # 1/2
        data_brm[[paste0('prior_S_', id)]] <- 1/2
        data_brm[[paste0('nu_', id)]] <- 1
      }
      insert_data <- c(glue("  // ----------------- [Manual Prior for IW {id}]"),
                       glue("  real nu_{id};"),
                       glue("  real prior_S_{id};"),
                       "  // -----------------")
      pos_data <- grep(code_brm, pattern = glue("vector\\[N\\] Z_{id}_1;"))
      if (length(pos_data) > 1){stop('---')}
      code_brm <- c(code_brm[1:(pos_data-1)], insert_data, code_brm[(pos_data):length(code_brm)])  
      
      pos_parameter <- grep(code_brm, pattern=glue("vector<lower=0>\\[M_{id}\\] sd_{id};  // group-level standard deviations"))
      insert_param <- c(glue("  vector<lower=0>[M_{id}] var_{id}; // group-level variance"))
      code_brm <- c(code_brm[1:(pos_parameter - 1)], insert_param,
                    code_brm[(pos_parameter + 1):length(code_brm)])      
      
      pos_model <- grep(code_brm, pattern=glue("vector\\[N_{id}\\] r_{id}_1;  // actual group-level effects"))
      
      code_brm <- c(code_brm[1:pos_model],
                    glue("  vector[M_{id}] sd_{id} = sqrt(var_{id});"),
                    code_brm[(pos_model+1):length(code_brm)])
      
      pos_model <- grep(code_brm, pattern=glue("lprior \\+= student_t_lpdf\\(sd_{id} \\| [0-9,\\. ]+\\)"))
      if (length(pos_model) != 1){stop('..')}
      if (grepl(code_brm[pos_model + 1], pattern='- 1 \\* student_t_lccdf(0 \\| [0-9,\\. ]+);')){
        stop('..')
      }
      insert_model <- c(glue("  lprior += inv_gamma_lpdf(var_{id}[1] | nu_{id}, prior_S_{id});"))
      code_brm <- c(code_brm[1:(pos_model[1]-1)], 
                    insert_model, 
                    code_brm[(pos_model[1] + 2):length(code_brm)])
      if (is_gaussian){
        pos_gen <- grep(code_brm, pattern=glue('^ *r_{id}.*'))
        if (length(pos_gen) != 1){stop('..')}
        code_brm[pos_gen] <- gsub(code_brm[pos_gen], pattern='= \\(', replacement = '= \\(sigma * ')
      }
      
    }else{
      
      if (prior_method == 'default_iw'){
        data_brm[[paste0('prior_S_', id)]] <- diag(size_v)
        data_brm[[paste0('nu_', id)]] <- size_v + 1
      }
      
      insert_data <- c(glue("  // ----------------- [Manual Prior for IW {id}]"),
                       glue("  cov_matrix[M_1] prior_S_{id};"),
                       glue("  real nu_{id};"),
                       "  // -----------------")
      pos_data <- grep(code_brm, pattern = glue("int<lower=1> NC_{id};"))
      if (length(pos_data) > 1){stop('---')}
      code_brm <- c(code_brm[1:(pos_data-1)], insert_data, code_brm[(pos_data + 1):length(code_brm)])  
      
      pos_parameter <- grep(code_brm, pattern=glue("cholesky_factor_corr\\[M_{id}\\] L_{id};|vector<lower=0>\\[M_{id}\\] sd_{id};"))
      if (diff(pos_parameter) != 2){stop('Misaligned')}
      
      insert_param <- c(glue("  cov_matrix[M_{id}] S_{id};  // covariance matrix"))
      code_brm <- c(code_brm[1:pos_parameter[1]-1], 
                    code_brm[pos_parameter[1] + 1],
                    insert_param,
                    code_brm[pos_parameter[2]:length(code_brm)]
      )
      
      
      pos_transform <- grep(code_brm, pattern=glue("r_{id} = \\(diag_pre_multiply\\(sd_{id}, L_{id}\\) \\* z_{id}\\)';"))
      if (length(pos_transform) == 0){
        pos_transform <- grep(code_brm, pattern=glue('scale_r_cor\\(z_{id}, sd_{id}, L_{id}\\)'))
        if (length(pos_transform) != 1){
          stop('Did not remove diag_pre_multiply/etc.')
        }
      }
      insert_transform <- c(glue("  cholesky_factor_cov[M_{id}] chol_{id} = cholesky_decompose(S_{id});"),
                            glue("  r_{id} = (chol_{id} * z_{id})';"))
      code_brm <- c(
        code_brm[1:(pos_transform-1)],
        insert_transform,
        code_brm[(pos_transform+1):length(code_brm)]
      )
      pat_model <- c(glue("lprior \\+= student_t_lpdf\\(sd_{id} \\| [0-9,\\. ]+\\)"),
                     glue("lprior \\+= lkj_corr_cholesky_lpdf\\(L_{id} \\| 1);"))
      pos_model <- grep(code_brm, pattern=paste(pat_model, collapse = '|'))
      insert_model <- c(glue("  lprior += inv_wishart_lpdf(S_{id} | nu_{id}, prior_S_{id});"))
      code_brm <- c(code_brm[1:(pos_model[1]-1)], code_brm[pos_model[1] + 2],
                    insert_model, code_brm[(pos_model[2]+1):length(code_brm)])
      
      pos_gen <- grep(code_brm, pattern=glue("cor_{id}\\[choose\\(k - 1, 2\\) \\+ j\\] = Cor_{id}\\[j, k\\];"))
      if (length(pos_gen) == 0){stop('posgen missing (cor)')}
      code_brm <- code_brm[-(pos_gen + -3:2)]
      pos_gen <- grep(code_brm, pattern = glue("corr_matrix\\[M_{id}\\] Cor_{id} = multiply_lower_tri_self_transpose\\(L_{id}\\);"))
      if (length(pos_gen) == 0){stop('posgen missing (corr)')}
      code_brm <- code_brm[-c(pos_gen + -1:1)]
      
      if (is_gaussian){
        pos_gen <- grep(code_brm, pattern=glue('r_{id} = \\(chol_{id}'))
        if (length(pos_gen) != 1){stop('..')}
        code_brm[pos_gen] <- gsub(code_brm[pos_gen], pattern='= \\(', replacement = '= sigma * (')
      }
    }
  }  
  
  pos_chol <- grep(code_brm, pattern=glue("cholesky_factor_cov\\[M_[0-9]+\\] chol_[0-9]+ = cholesky_decompose\\(S_[0-9]+\\);"))
  if (length(pos_chol) > 0){
    select_chol <- code_brm[pos_chol]
    code_brm <- code_brm[-pos_chol]
    code_brm <- c(code_brm[1:(pos_chol[1]-1)], select_chol,
                  code_brm[(pos_chol[1] ):length(code_brm)])
  }
  
  
  
  fix_z <- grep(code_brm, pattern= 'target \\+= std_normal_lpdf\\(z_[0-9]+\\[1\\])')
  if (length(fix_z) > 0){
    code_brm[fix_z] <- gsub(code_brm[fix_z], pattern='target \\+= std_normal_lpdf\\((z_([0-9]+))(\\[1\\]\\));', 
                            replacement = 'for (j in 1:M_\\2){\n    target += std_normal_lpdf(z_\\2[j]);\n}')
  }
  
  # Flat prior on intercept
  code_intercept <- grep(code_brm, pattern='lprior.*Intercept')
  code_brm <- code_brm[-code_intercept]
  if (is_gaussian){
    # Adjust prior for sigma to be 1/sigma^2
    code_sigma <- grep(code_brm, pattern='lprior \\+= student_t_lpdf\\(sigma')
    if (!grepl(code_brm[code_sigma + 1], pattern='student_t_lccdf')){
      stop('...')
    }
    code_brm <- code_brm[-(code_sigma + 1)]
    code_brm[code_sigma] <- '  lprior += -2 * log(sigma);'
    
  }
  final_brm <<- code_brm
  code_brm <- paste(code_brm, collapse='\n')
  
  data_brm <- data_brm[!(names(data_brm) %in% paste0('NC', gsub(re_id, pattern='^M', replacement = '')))]
  
  if (compile){
    message('Compiling Again')
    compile_adjusted <- stan_model(model_code = code_brm, save_dso = TRUE)
    #Check convergence
    convg <- suppressMessages(sampling(compile_adjusted, data = data_brm, chains = 0))
  }else{
    compile_adjusted <- NULL
  }
  
  output <- list(model = compile_adjusted, data = data_brm,
                 code = code_brm, orig_brm = brm_object)
  return(output)
}

data <- readRDS('data/gg_data.RDS')
data$fake_outcome <- as.numeric(data$success > data$failure)
data$frac_success <- data$success
data$frac_failure <- data$failure

data$success <- round(data$success)
data$failure <- round(data$failure)
data$n_trials <- data$success + data$failure

formula_list <- readRDS('data/formula_list.RDS')

for (y in yr_list){
  for (model in mdl_list){
    for (outcome in outcome_list){
      message(c(y, '-', outcome, '-', model))
      fmla <- formula_list[[model]]
      
      reg_data <- data %>% filter(year == y) %>% filter(n_trials > 0)
      fmla <- update.formula(fmla, 'success | trials(n_trials) ~ .')
      reg_family <- list(family = binomial())
      
      #Run using STAN
      adjust_brm <- convert_brm(
        formula = fmla, 
        data = reg_data, 
        brms_args = reg_family)
      
      custom_exclude <- adjust_brm$orig_brm$exclude
      if (is.null(custom_exclude)){
        custom_exclude <- brms:::exclude_pars(adjust_brm$orig_brm)
        if (is.null(custom_exclude)){stop('No exclusion??')}
        custom_exclude <- c(custom_exclude, c('chol_', 'sigma'))
      }else{
        custom_exclude <- c(custom_exclude, c('chol_', 'sigma'))
      }
      print(custom_exclude)

      time_start <- proc.time()
      mcmc_iw <- sampling(object = adjust_brm$model, data = adjust_brm$data,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 15), 
            refresh = 100, 
            iter = 10^4,
            chains = 4, include = FALSE, pars = custom_exclude)
      time_stan <- proc.time() - time_start
      
      adjust_brm$orig_brm$fit <- mcmc_iw
      adjust_brm$orig_brm <- brms:::rename_pars(adjust_brm$orig_brm)
      
      saveRDS(list(model = adjust_brm, time = time_stan, year = y, outcome = outcome),
              glue('output/output_stan/brms_y_{y}_m_{model}_o_{outcome}.RDS'),
              version = 2)
      
      rm(mcmc_iw, adjust_brm)
      gc()
      
    }
  }  
}