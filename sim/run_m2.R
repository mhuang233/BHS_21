### ::: Run ::: ###
{
  # rm(list = ls())
  library(loo)
  library(rstan)
  # library(gtools) # for combinations
  library(rstanarm)
  library(LaplacesDemon)
  library(tidyverse)
  
  source("bhs_fun.R")
  options(scipen = 999)
  options(mc.cores = 4)
}

# - DON'T RUN EXCEPT CHTC #################### CHTC - Starts ################# #

args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
}

arguments = commandArgs(trailingOnly=TRUE)
rep = as.numeric(arguments[[1]])
c = as.numeric(arguments[[2]])
ni = as.numeric(arguments[[3]])
nj = as.numeric(arguments[[4]])
icc = as.numeric(arguments[[5]])

# ##################### CHTC - Ends ######################### - DON'T RUN ENDS #

set.seed(c)

# dgf
{
  ni = ni
  nj = nj
  gamma00 = 400
  gamma01 = -20
  gamma02 = 6
  gamma03 = -2
  gamma04 = 8
  gamma05 = -4
  u_0 = sqrt(icc)
  u_1 = sqrt(icc)
  u_2 = sqrt(icc)
  w_0 = sqrt(1-icc)
  w_1 = sqrt(1-icc)
  w_2 = sqrt(1-icc)
  sigma = 1

  x = matrix(99, nrow = ni*nj, ncol = 5)
  y = matrix(99, nrow = ni*nj, ncol = 1) 
  r = matrix(99, nrow = ni*nj, ncol = 1) 
  U = matrix(99, nrow = ni*nj, ncol = 3)
  W = matrix(99, nrow = ni*nj, ncol = 3)
  sim = matrix(99, nrow = ni*nj, ncol = 12) 
}


df <- dgf(ni, nj, gamma00, gamma01, gamma02, gamma03, gamma04, gamma05,
          u_0, u_1, u_2, w_0, w_1, w_2, sigma)

nr <- nrow(df)
df1 <- df[1:(nr/2), ]
df2 <- df[(nr/2+1):nr, ]

# combinations for covariates
covariates <- df %>% select(x1, x2, x3, x4, x5) %>% names()

#library(gtools)
cov1 <- combinations(4, 1, covariates, repeats.allowed = F)
cov2 <- combinations(4, 2, covariates, repeats.allowed = F)
cov3 <- combinations(4, 3, covariates, repeats.allowed = F)
cov4 <- combinations(4, 4, covariates, repeats.allowed = F)

cov1 <- matrix(c("x1", "x2", "x3", "x4"), ncol = 1, nrow = 4, byrow = T)
cov3 <- matrix(c("x1", "x2", "x3",
                 "x1", "x2", "x4",
                 "x1", "x2", "x5",
                 "x1", "x3", "x4",
                  "x1", "x3", "x5",
                  "x1", "x4", "x5",
                  "x2", "x3", "x4",
                  "x2", "x3", "x5",
                  "x2", "x4", "x5",
                  "x3", "x4", "x5"), ncol = 3, nrow = 10, byrow = T)

# function
extract_lpd <- function(x){
  
  g <- x$pointwise[, "elpd_loo"]
  return(g)
  
}

# fit, estimate, and analyze
alz_lm(rep = rep, c = c, ni = ni, nj = nj, icc = icc, df1 = df1, df2 = df2, cov1 = cov1, cov3 = cov3)
alz_s2(rep = rep, c = c, ni = ni, nj = nj, icc = icc, df1 = df1, df2 = df2, cov1 = cov1, cov3 = cov3)
alz_all(rep = rep, c = c, ni = ni, nj = nj, icc = icc, df1 = df1, df2 = df2, cov1 = cov1, cov3 = cov3)










