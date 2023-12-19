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

accuracy_bayesian <- function(model_1, model_2, names){
  require(KernSmooth)
  output <- sapply(names, FUN=function(k){
    samples_1 <- model_1[,k]
    samples_2 <- model_2[,k]
    range_s <- range(c(samples_1, samples_2))
    kern_1 <- suppressWarnings(tryCatch(list(warning = 0, kern = bkde(samples_1, range.x = range_s)), warning = function(w){list(warning = 1, kern = bkde(samples_1, range.x = range_s))}))
    kern_2 <- suppressWarnings(tryCatch(list(warning = 0, kern = bkde(samples_2, range.x = range_s)), warning = function(w){list(warning = 1, kern = bkde(samples_2, range.x = range_s))}))
    stopifnot(isTRUE(all.equal(kern_1$kern$x, kern_2$kern$x)))
    #Get the difference of the heights
    iae <- sum(abs(kern_1$kern$y - kern_2$kern$y))
    #Normalize by the bandwidths
    iae <- mean(diff(kern_1$kern$x)) * iae
    if (is.na(iae)){
      iae <- 2
      kern_1$warning <- TRUE
      kern_2$warning <- TRUE
    }else{
      if (iae < 0 | iae > 2){stop('IAE in invalid range')}
    }
    acc <- 1 - 1/2 * iae    
    return(c(acc, kern_1$warning, kern_2$warning))
  })
  output <- data.frame(t(output)) 
  names(output) <- c('accuracy', 'warn_1', 'warn_2')
  output <- output %>% mutate(name = rownames(.))
  return(output)
}

get_uqf <- function(model_vi, MCMC_covariance, MCMC_samples, allow_unequal = FALSE){
  require(RSpectra); require(caret)
  full_precision_VI <- vglmer:::extract_precision(model_vi)
  param_names <- c(rownames(model_vi$beta$mean), rownames(model_vi$alpha$mean))
  both_names <- intersect(colnames(MCMC_covariance), param_names)

  if (!allow_unequal){
    stopifnot(length(both_names) == length(param_names))
    stopifnot(length(both_names) == ncol(MCMC_covariance))
  }

  # precision_VI <- full_precision_VI[
  #   match(both_names, param_names), 
  #   match(both_names, param_names)
  # ]
  cov_VI <- solve(full_precision_VI)[
    match(both_names, param_names), match(both_names, param_names)
  ]
  precision_VI <- solve(cov_VI)

  product_matrix <- precision_VI %*% MCMC_covariance[both_names,both_names] 
  # Get the leading eigenvector and eigenvalue 
  rspectra_fit <- RSpectra::eigs(product_matrix, k = 1, which = 'LM')
  rspectra_ev <- rspectra_fit$values
  stopifnot(all(Im(rspectra_ev) == 0))
  rspectra_ev <- Re(rspectra_ev)
  uqf_rspectra <- 1/max(rspectra_ev)
  # Now, calculate it a different way
  eigenvector <- Re(rspectra_fit$vector)

  # "both_names" refers to the parameters that appear in both the MCMC
  # and VI; (levels that are not in the estimation data are not estimated with VI)
  # Proposal 2: 
  var_MCMC <- var(MCMC_samples[,both_names] %*% as.vector(eigenvector))
  cov_VI <- solve(full_precision_VI)[match(both_names, param_names), match(both_names, param_names)]
  var_VI <- t(eigenvector) %*% cov_VI %*% eigenvector
  uqf_ratio <- as.numeric(var_VI)/as.numeric(var_MCMC)

  # Proposal 3: Split-Sample MCMC
  folds <- caret::createFolds(y = 1:nrow(MCMC_samples), k = 5)

  uqf_split_sample <- sapply(folds, FUN=function(fold){
    # Get covariance on the K-1 folds
    MCMC_cov_sample <- cov(MCMC_samples[-fold,])[both_names, both_names]
    # Get top 5 eigenvectors
    product_matrix_sample <- precision_VI %*% MCMC_cov_sample
    eigen_prod_sample <- Re(RSpectra::eigs(product_matrix_sample, k = 50, which = 'LM')$vector)
    # Get the covariance of held out fold (K) (e.g., 20%) 
    # against the eigenvectors (P x 5 matrix)
    var_MCMC <- cov(MCMC_samples[fold, both_names] %*% eigen_prod_sample)
    # Get the precision of the VI against this; estimate covariance the covariance from VI
    # and then invert...
    prec_VI <- solve(t(eigen_prod_sample) %*% cov_VI %*% eigen_prod_sample)
    uqf_split_sample <- 1/eigen(var_MCMC %*% prec_VI)$values[1]
    return(uqf_split_sample)
  })
  uqf_split_first <- uqf_split_sample[1]
  uqf_split_sample <- mean(uqf_split_sample)

  return(data.frame(uqf_rspectra = uqf_rspectra,
    uqf_split_first = uqf_split_first,
    uqf_ratio = uqf_ratio, 
    uqf_split_sample = uqf_split_sample))
}

get_dist_eigen <- function(
    model_vi, 
    MCMC_covariance, 
    MCMC_samples, allow_unequal = FALSE){
  require(RSpectra); require(caret)
  full_precision_VI <- vglmer:::extract_precision(model_vi)
  param_names <- c(rownames(model_vi$beta$mean), rownames(model_vi$alpha$mean))
  both_names <- intersect(colnames(MCMC_covariance), param_names)
  
  if (!allow_unequal){
    stopifnot(length(both_names) == length(param_names))
    stopifnot(length(both_names) == ncol(MCMC_covariance))    
  }

  # precision_VI <- full_precision_VI[
  #   match(both_names, param_names), 
  #   match(both_names, param_names)
  # ]
  cov_VI <- solve(full_precision_VI)[
    match(both_names, param_names), match(both_names, param_names)
  ]
  precision_VI <- solve(cov_VI)

  # Proposal 3: Split-Sample MCMC
  folds <- caret::createFolds(y = 1:nrow(MCMC_samples), k = 5)

  uqf_split_sample <- lapply(folds, FUN=function(fold){
    # Get covariance on the K-1 folds
    store_out <- data.frame()
    for (type in c('LM', 'SM')){
      MCMC_cov_sample <- cov(MCMC_samples[-fold,])[both_names, both_names]
      # Get top 5 eigenvectors
      product_matrix_sample <- precision_VI %*% MCMC_cov_sample
      ncv_size <- min(c(250, ncol(product_matrix_sample)))
      eigen_prod_sample <- Re(RSpectra::eigs(product_matrix_sample,
       k = 25, opts = list(ncv = ncv_size),
       which = type)$vector)
      if (ncol(eigen_prod_sample) < 25){message('FAILED TO OBTAIN 25 EIGENVECTORS')}
      # Get the covariance of held out fold (K) (e.g., 20%) 
      # against the eigenvectors (P x 5 matrix)
      var_MCMC <- cov(MCMC_samples[fold, both_names] %*% eigen_prod_sample)
      # Get the precision of the VI against this; estimate covariance the covariance from VI
      # and then invert...
      prec_VI <- solve(t(eigen_prod_sample) %*% cov_VI %*% eigen_prod_sample)      
      dist_ev <- 1/eigen(var_MCMC %*% prec_VI)$values
      out <- data.frame(eigenvalue = dist_ev, id = 1:length(dist_ev))
      out$type <- type
      store_out <- rbind(store_out, out)
    }
    return(store_out)
  })
  uqf_split_sample <- bind_rows(uqf_split_sample, .id = 'fold')

  return(uqf_split_sample)
}

get_marginal <- function(model_vi, MCMC_covariance, MCMC_samples){
  
  require(caret)

  full_precision_VI <- vglmer:::extract_precision(model_vi)
  param_names <- c(rownames(model_vi$beta$mean), rownames(model_vi$alpha$mean))
  both_names <- intersect(colnames(MCMC_covariance), param_names)

  stopifnot(length(both_names) == length(param_names))
  stopifnot(length(both_names) == ncol(MCMC_covariance))

  # precision_VI <- full_precision_VI[match(both_names, param_names), 
  #   match(both_names, param_names)]
  cov_VI <- solve(full_precision_VI)[
    match(both_names, param_names), match(both_names, param_names)
  ]
  precision_VI <- solve(cov_VI)

  stopifnot(all.equal(diag(cov_VI), diag(solve(precision_VI))))

  ratio_fullsample <- min(diag(cov_VI)/diag(MCMC_covariance)[both_names])

  folds <- caret::createFolds(y = 1:nrow(MCMC_samples), k = 5)
  min_ratio_split <- sapply(folds, FUN=function(fold){
    MCMC_cov_sample <- cov(MCMC_samples[-fold,])[both_names, both_names]
    ratio_variance <- diag(cov_VI)/diag(MCMC_cov_sample)
    # Get the five worst ratios
    worst_k <- order(ratio_variance, decreasing = FALSE)[1:50]
    # Get the variance of the five worst variables on the held out sample
    MCMC_split_var <- diag(cov(MCMC_samples[fold, both_names][, worst_k]))
    VI_split_var <- diag(cov_VI)[worst_k]
    # Get the ratio of the worst variable across those five
    ratio_worst <- min(VI_split_var/MCMC_split_var)
    return(ratio_worst)
  })

  mean_split_first <- min_ratio_split[1]
  mean_split_sample <- mean(min_ratio_split)

  return(data.frame(
    ratio_fullsample = ratio_fullsample,
    ratio_split_first = mean_split_first,
    ratio_split = mean_split_sample
  ))
}
