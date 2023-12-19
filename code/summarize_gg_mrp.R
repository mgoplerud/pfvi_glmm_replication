rm(list=ls())
library(tidyverse)
library(glue)
library(blme)
library(reshape2)
library(brms)
library(Matrix)
library(vglmer)

args <- commandArgs(trailingOnly=TRUE)
id <- as.numeric(args[1])

year <- id %% 2 + 1
mdl <- id %% 9 + 1

print(c(year, mdl))

id_list <- list(
  c('stt'), 
  c('stt', 'eth'),
  c('stt', 'eth', 'inc'),
  c('stt', 'eth', 'inc', 'age')
)
id_list <- id_list[c(1,4)]

source('code/aux_functions.R')

gg.data <- readRDS('data/gg_data.RDS')
gg.data <- gg.data[gg.data$year == year,]
gg.data$total_survey <- round(gg.data$success) + round(gg.data$failure)
gg.data <- gg.data[, c('stt', 'eth', 'age', 'inc', 'total_survey')]

sum.pop <- readRDS('data/gg_poststrat.RDS')
sum.state <- readRDS('data/gg_statevalues.RDS')


dat.out <- left_join(sum.pop, 
                     sum.state %>% select(stt, matches('^z')), by = 'stt')
dat.out <- dat.out %>% filter(!(stt == -1 | eth == -1 | inc == -1 | age == -1))

L.z.incstt <- list('z.inc2004', 'z.inc2007')
L.z.repprv <- list('z.rep2000', 'z.rep2004')
L.z.trnprv <- list('z.trn2000', 'z.trn2004')
L.pop <- list('pop2004', 'pop2008')

reg.lkup <- c(3,4,4,3,4,4,1,1,5,3,3,4,4,2,2,2,2,3,3,1,1,1,2,2,3,2,4,2,4,1,1,4,1,3,2,2,3,4,1,1,3,2,3,3,4,1,3,4,1,2,4)
dat.out$reg <- reg.lkup[dat.out$stt]

for (r in 2:4){
  combo <- combn(c('reg', 'stt', 'eth', 'inc', 'age'), r)
  for (v in array_branch(combo, margin = 2)){
    if ('stt' %in% v & 'reg' %in% v){
      #Skip if both state and region.
      next
    }
    dat.out[,paste(v, collapse='.')] <- apply(dat.out[, v,drop=F], MARGIN = 1, FUN=function(i){paste(i, collapse='-')})
  }
}

mdls_ps <- c('vglmer_weak', 'vglmer_strong', 'vglmer_collapse_FE', 'vglmer_collapse_main', 
             'vglmer_collapse_all')

m_HMC <- readRDS(glue(
    'output/output_stan/brms_y_{year}_m_{mdl}_o_binomial.RDS'
))

copy_dat <- dat.out

mf <- model.frame(m_HMC$model$orig_brm)
unique_inc <- mf %>% select(inc, z.inc)  %>% group_by(inc, z.inc) %>% tally() %>% select(-n)
stopifnot(nrow(unique_inc) == 5)
copy_dat$z.incstt <- copy_dat[[L.z.incstt[[year]]]]
copy_dat$z.trnprv <- copy_dat[[L.z.trnprv[[year]]]]
copy_dat$pop_year <- copy_dat[[L.pop[[year]]]]
copy_dat <- left_join(copy_dat, unique_inc, by = 'inc')
copy_dat <- left_join(copy_dat, gg.data, by = c("stt", "eth", "age", "inc"))
stopifnot(nrow(copy_dat) == 4080)

pred_HMC <- posterior_linpred(
  m_HMC$model$orig_brm, newdata = copy_dat %>% mutate(n_trials = 1), allow_new_levels = TRUE
)
pred_HMC <- plogis(pred_HMC)

all_acc <- data.frame()
for (type in mdls_ps){
  print(c(type, year, mdl))
  
  m_vglmer <- readRDS(glue('output/output_vglmer/{type}_y_{year}_m_{mdl}_o_binomial.RDS'))
  
  pred_vglmer <- predict(
    m_vglmer,
    newdata = copy_dat, samples = 4000, summary = F
  )
  pred_vglmer <- plogis(as.matrix(pred_vglmer))
  colnames(pred_HMC) <- colnames(pred_vglmer) <- 1:ncol(pred_HMC)
  
  all_id <- data.frame()
  for (id in id_list){
    print(id)
    id_name <- paste(id, collapse='_')
    jid <- apply(copy_dat[,id], MARGIN = 1, paste, collapse = '_')
    uid <- unique(jid)
    if (length(uid) == 4080){
      M <- sparseMatrix(i = 1:4080, j = 1:4080, x = 1)
      colnames(M) <- uid
    }else{
      M <- sparseMatrix(i = 1:nrow(copy_dat), j = match(jid, uid), x = copy_dat$pop_year)
      M <- M %*% Diagonal(x = 1/colSums(M))
      colnames(M) <- uid
    }
    
    zero_id <- which(apply(pred_HMC %*% M, MARGIN = 2, var) %in% c(0, NA))
    est_id <- accuracy_bayesian(pred_HMC %*% M, pred_vglmer %*% M, names = setdiff(uid, uid[zero_id]))
    est_id$combo <- id_name
    
    est_id$HMC_var <- apply(pred_HMC %*% M, MARGIN = 2, var)
    est_id$vglmer_var <- apply(pred_vglmer %*% M, MARGIN = 2, var)

    # Change into a matrix for summation
    M@x <- rep(1, length(M@x))
    pop_by_cell <- as.vector(t(M) %*% copy_dat$pop_year)
    survey_by_cell <- as.vector(t(M) %*% (copy_dat$total_survey))
    est_id$poststrat_population <- pop_by_cell
    est_id$survey_population <- survey_by_cell 
    all_id <- bind_rows(all_id, est_id)
  }
  
  all_acc <- bind_rows(
    all_acc, all_id %>% mutate(method = type, year = year, model = mdl)
  )
  rm(pred_vglmer);
}

path <- glue('output/output_mrp/out_{year}_{mdl}.RDS')
saveRDS(all_acc, path)   
