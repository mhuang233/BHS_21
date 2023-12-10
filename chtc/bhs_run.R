### ::: Run ::: ###
{
  rm(list = ls())
  library(loo)
  library(rstan)
  # library(gtools) for combinations
  library(rstanarm)
  library(LaplacesDemon)
  library(tidyverse)
  
  source("bhs_fun.R")
  options(scipen = 999)
  options(mc.cores = parallel::detectCores())
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
sd = as.numeric(arguments[[5]])
sigma = as.numeric(arguments[[6]])

# ##################### CHTC - Ends ######################### - DON'T RUN ENDS #

set.seed(c)

# dgf
{
  ni <- ni
  nj <- nj
  gamma00 <- 400
  gamma01 <- -20
  gamma02 <- 6
  gamma03 <- -2
  gamma04 <- 8
  u_0 <- sd
  w_0 <- sd
  sigma <- sigma # when sigma = 0, sigma = 5, sigma = 10, sigma = 30
  
  x <- matrix(99, nrow = ni*nj, ncol = 4)
  y <- matrix(99, nrow = ni*nj, ncol = 1) 
  r <- matrix(99, nrow = ni*nj, ncol = 1) 
  U <- matrix(99, nrow = ni*nj, ncol = 3)
  W <- matrix(99, nrow = ni*nj, ncol = 4)
  sim <- matrix(99, nrow = ni*nj, ncol = 9) 
}


df <- dgf(ni, nj, gamma00, gamma01, gamma02, gamma03, gamma04, u_0, w_0, sigma)

#nr <- nrow(df)
#df1 <- df[1:(nr/2), ]
#df2 <- df[(nr/2+1):nr, ]

# function
#extract_lpd <- function(x){
#  g <- loo::loo(x, r_eff = relative_eff(exp(x), chain_id = 1:nrow(x)))$pointwise
#  return(g)
#}

extract_lpd <- function(x){
  
  g <- x$pointwise[, "elpd_loo"]
  return(g)
  
}

model_strings <- c(
  "y ~ x1 + (1|i)",
  "y ~ x2 + (1|i)",
  "y ~ x3 + (1|i)",
  "y ~ x4 + (1|i)",
  "y ~ x1 + x2 + (1|i)",
  "y ~ x1 + x3 + (1|i)",
  "y ~ x1 + x4 + (1|i)",
  "y ~ x2 + x3 + (1|i)",
  "y ~ x2 + x4 + (1|i)",
  "y ~ x3 + x4 + (1|i)",
  "y ~ x1 + x2 + x3 + (1|i)",
  "y ~ x1 + x2 + x4 + (1|i)",
  "y ~ x1 + x3 + x4 + (1|i)",
  "y ~ x2 + x3 + x4 + (1|i)",
  "y ~ x1 + x2 + x3 + x4 + (1|i)"
)

var_names <- c(
  "x1",
  "x2",
  "x3",
  "x4",
  "x1_x2",
  "x1_x3",
  "x1_x4",
  "x2_x3",
  "x2_x4",
  "x3_x4",
  "x1_x2_x3",
  "x1_x2_x4",
  "x1_x3_x4",
  "x2_x3_x4",
  "x1_x2_x3_x4"
) 

# run
alz(rep, c, ni, nj, icc, df, model_strings, var_names)



