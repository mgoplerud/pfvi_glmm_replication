# This file is taken from previous work by one of the authors 
# (Goplerud; https://doi.org/10.7910/DVN/DI19IB)

# The first few lines of the replication data for Ghitza and Gelman (2013).
# Used to prepare the dataset used in the analysis.

library(arm)
library(car)
# no longer available on CRAN
# Can install with this. Or easily replaced with dplyr
# remotes::install_version('mapReduce', version = '1.2.6')
library(mapReduce)
library(maps)
library(purrr)
library(dplyr)

#Mostly verbatim from G&G replication R script
rm(list = ls())

#TO begin, download data from replication archive and install into folder
#such as "build_gg"
#https://doi.org/10.7910/DVN/PZAOO6

########################################################################################
### helper functions

WMean <- function(a, w=rep(1,length(a)), subset=rep(TRUE,length(a))) {
  keep <- !is.na(a) & !is.na(w) & !is.na(subset) & subset
  return(sum((w*a)[keep])/sum(w[keep]))
}

FindDelta <- function(delta, a, w, x0)
  abs(x0-sum(invlogit(logit(a) + delta)*w))

CorrectWeighted <- function(a, w, x0) {
  delta <- optimize(FindDelta, interval=c(-5,5), a, w, x0)$minimum
  corrected <- invlogit(logit(a) + delta)
  return(list(delta=delta, corrected=corrected))
}

logit <- function (a) log(a/(1-a))

########################################################################################
### labels

reg.lkup <- c(3,4,4,3,4,4,1,1,5,3,3,4,4,2,2,2,2,3,3,1,1,1,2,2,3,2,4,2,4,1,1,4,1,3,2,2,3,4,1,1,3,2,3,3,4,1,3,4,1,2,4)

reg.label <- c("Northeast", "Midwest", "South", "West", "DC")
stt.label <- c(state.abb[1:8], "DC", state.abb[9:50])
eth.label <- c("White", "Black", "Hispanic", "Other")
inc.label <- c("$0-20k", "$20-40k", "$40-75k", "$75-150k", "$150k+")
age.label <- c("18-29", "30-44", "45-64", "65+")
sex.label <- c("Male", "Female")
edu.label <- c("< HS", "HS", "Some Coll", "Coll", "Post-Grad")
mar.label <- c("Married", "Single")
rel.label <- c("Protestant (not born-again)","Born-again Protestant","Catholic","Mormon","Jewish","Other religion","No religion")
chu.label <- c("Nonattenders","Rare attenders","Occasional\nattenders","Frequent attenders","Very frequent\nchurch attenders")
kid.label <- c("Kids", "No Kids")
cit.label <- c("Citizen", "Not Citizen")

n.reg <- length(reg.label)
n.stt <- length(stt.label)
n.eth <- length(eth.label)
n.inc <- length(inc.label)
n.age <- length(age.label)
n.sex <- length(sex.label)
n.edu <- length(edu.label)
n.mar <- length(mar.label)
n.rel <- length(rel.label)
n.chu <- length(chu.label)
n.kid <- length(kid.label)

########################################################################################
### get data

### state-level data
dat.stt <- read.table("build_gg/state-stats.dat", header=T, sep="\t")
dat.stt$z.inc2000 <- rescale(dat.stt$inc2000)
dat.stt$z.inc2004 <- rescale(dat.stt$inc2004)
dat.stt$z.inc2007 <- rescale(dat.stt$inc2007)
dat.stt$z.rep1996 <- rescale(dat.stt$rep1996)
dat.stt$z.rep2000 <- rescale(dat.stt$rep2000)
dat.stt$z.rep2004 <- rescale(dat.stt$rep2004)
dat.stt$z.rep2008 <- rescale(dat.stt$rep2008)
dat.stt$z.trn1996 <- rescale(dat.stt$vote1996/dat.stt$pop1996)
dat.stt$z.trn2000 <- rescale(dat.stt$vote2000/dat.stt$pop2000)
dat.stt$z.trn2004 <- rescale(dat.stt$vote2004/dat.stt$pop2004)
dat.stt$z.trn2008 <- rescale(dat.stt$vote2008/dat.stt$pop2007)
dat.stt$stt <- 1:nrow(dat.stt)

### CPS data for modeling turnout
dat.cps <- read.table("build_gg/cps2000-04-08-DKs.dat", header=T, sep="\t")
ok <- apply(is.na(dat.cps[1:8]) , 1, sum)==0 & (dat.cps$cit==1 | is.na(dat.cps$cit))
dat.cps <- dat.cps[ok,]
dat.cps$reg <- reg.lkup[dat.cps$stt]
dat.cps$year <- as.numeric(gsub(as.character(dat.cps$file), pattern='[^0-9]', replacement = ''))
dat.cps <- dat.cps[dat.cps$year != 2000,]

### annenberg / pew data for modeling vote choice
dat.vot <- read.table("build_gg/votechoice2000-04-08.dat", header=T, sep="\t")
dat.vot$weight[is.na(dat.vot$weight)] <- 1
ok <- apply(is.na(dat.vot[1:8]), 1, sum)==0 & 
      (dat.vot$cit==1 | is.na(dat.vot$cit)) & 
      (dat.vot$regist=="registered" | is.na(dat.vot$regist))
dat.vot <- dat.vot[ok,]
dat.vot$reg <- reg.lkup[dat.vot$stt]
dat.vot$year <- as.numeric(gsub(as.character(dat.vot$file), pattern='[^0-9]', replacement = ''))
dat.vot <- dat.vot[dat.vot$year != 2000,]

### census data from PUMS for population cell sizes
dat.pop <- read.table("build_gg/census-pums-pop-2000-04-08.dat", header=T, sep="\t")

### prepare data for looping through years
years <- c(2004, 2008)
L.z.incstt <- list(dat.stt$z.inc2004, dat.stt$z.inc2007)
L.z.repprv <- list(dat.stt$z.rep2000, dat.stt$z.rep2004)
L.z.trnprv <- list(dat.stt$z.trn2000, dat.stt$z.trn2004)

run.full.models <- FALSE

########################################################################################
### Run MRP for stt, eth, inc, age

D <- as.data.frame(expand.grid(1:n.stt, 1:n.eth, 1:n.inc, 1:n.age))
colnames(D) <- c("stt", "eth", "inc", "age")
D$grp <- apply(D, 1, paste, collapse="_")
D$ix <- 1:nrow(D)

### multilevel models
M.cps <- M.vot <- list()

store_data <- data.frame()
for (i.year in 1:2) {
  cat(paste("***** Multilevel Models for", years[i.year], "\n"))
  
  ### covariates
  stt <- D$stt
  eth <- D$eth
  inc <- D$inc
  age <- D$age
  reg <- reg.lkup[stt]
  reg.eth <- 10*reg + eth
  reg.inc <- 10*reg + inc
  reg.age <- 10*reg + age
  stt.eth <- 10*stt + eth
  stt.inc <- 10*stt + inc
  stt.age <- 10*stt + age
  eth.inc <- 10*eth + inc
  eth.age <- 10*eth + age
  inc.age <- 10*inc + age
  reg.eth.inc <- 100*reg + 10*eth + inc
  stt.eth.inc <- 100*stt + 10*eth + inc
  reg.eth.age <- 100*reg + 10*eth + age
  reg.inc.age <- 100*reg + 10*inc + age
  stt.eth.age <- 100*stt + 10*eth + age
  stt.inc.age <- 100*stt + 10*inc + age
  eth.inc.age <- 100*eth + 10*inc + age
  reg.eth.inc.age <- 1000*reg + 100*eth + 10*inc + age
  stt.eth.inc.age <- 1000*stt + 100*eth + 10*inc + age
  z.inc <- rescale(inc)
  z.incstt <- L.z.incstt[[i.year]][stt]
  z.trnprv <- L.z.trnprv[[i.year]][stt]
  z.repprv <- L.z.repprv[[i.year]][stt]
  

  tmp <- dat.cps[dat.cps$year==years[i.year],]
  tmp$grp <- apply(tmp[, c("stt", "eth", "inc", "age")], 1, paste, collapse="_")
  tmp$ones <- 1

  ### turnout model
  cat("*****   CPS Turnout Model\n")
  mr <- as.data.frame(mapReduce(data=tmp, map=grp,
                                n=sum(ones),
                                ybar.wt=sum(vote*weight)/sum(weight),
                                des.eff.cell=1 + var(weight/mean(weight))))
  mr$grp <- rownames(mr)
  D.tmp <- merge(x=D, y=mr, by="grp", all.x=T)  
  D.tmp <- D.tmp[order(D.tmp$ix),]
  D.tmp$n[is.na(D.tmp$n)] <- 0
  des.eff <- WMean(D.tmp$des.eff.cell, D.tmp$n, D.tmp$n > 1)
  D.tmp$n.eff <- D.tmp$n/des.eff
  D.tmp$ybar.wt[D.tmp$n.eff==0] <- 0.5
  y <- cbind(D.tmp$ybar.wt * D.tmp$n.eff, (1 - D.tmp$ybar.wt) * D.tmp$n.eff)

  
  stt <- D$stt
  eth <- D$eth
  inc <- D$inc
  age <- D$age
  reg <- reg.lkup[stt]
  
  ###
  ### Extra additions to create relevant columns.
  ###
  fmt_data <- data.frame(success = y[,1], failure = y[,2], 
                         stt, reg, eth, inc, age, z.inc,
             z.incstt,
             z.trnprv,
             z.repprv)
  
  for (r in 2:4){
    combo <- combn(c('reg', 'stt', 'eth', 'inc', 'age'), r)
    for (v in array_branch(combo, margin = 2)){
      if ('stt' %in% v & 'reg' %in% v){
        #Skip if both state and region.
        next
      }
      fmt_data[,paste(v, collapse='.')] <- apply(fmt_data[, v,drop=F], MARGIN = 1, FUN=function(i){paste(i, collapse='-')})
    }
  }
  fmt_data$year <- i.year
  store_data <- rbind(store_data, fmt_data)
}

saveRDS(store_data, 'data/gg_data.RDS')

# Get the postratification matrix, i.e. number of people
# of each demographic-type in each state
store_poststrat <- dat.pop %>%
  group_by(stt, eth, inc, age) %>%
  summarize(pop2004 = sum(wtd2004), pop2008 = sum(wtd2008)) %>%
  ungroup
saveRDS(store_poststrat, 'data/gg_poststrat.RDS')
# Save out state-level covariates
saveRDS(dat.stt, 'data/gg_statevalues.RDS')
