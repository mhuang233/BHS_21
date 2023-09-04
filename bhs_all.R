
### ::: All for simulation ::: ###

# #- DON'T RUN EXCEPT CHTC #################### CHTC - Starts ##################

args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
}

arguments = commandArgs(trailingOnly=TRUE)
n_iter = as.numeric(arguments[[1]])
n_rep = as.numeric(arguments[[2]])
ni = as.numeric(arguments[[3]])
nj = as.numeric(arguments[[4]])
seed = as.numeric(arguments[[5]])

# ##################### CHTC - Ends ########################### - DON'T RUN ENDS
source("fun_bhs.R")

### Library set up ###
{
  
  library(tidyverse)
  library(rstanarm)
  library(rstan)
  library(loo)
  # library(stats)
  # library(bayesplot)
  library(LaplacesDemon)
  library(parallel)
  options(scipen = 999)
  options(mc.cores = detectCores())
  
}


### Parameter specification ###
{
  n_iter = n_iter
  n_rep = n_rep
  ni = ni # of cluster
  nj = nj # of students in each school
  N = ni*nj
  gamma00 = 400
  gamma01 = 14 # binary for public
  gamma02 = 10 # excs
  gamma03 = 9 # homepos
  gamma04 = -20 # pisadiff
  gamma05 = 6 # gfofail
  gamma06 = -2 # belong
  gamma07 = 7 # compete
  gamma08 = 1 # adaptivity
  gamma09 = 20 # metsum
  gamma10 = -14 # metgoal
  gamma11 = -1 # swbp
  gamma12 = 7 # workmast
  gamma13 = 8 # teachint
  gamma14 = -4 # icres
  u_0 <- .4
  u_1 <- .2
  u_2 <- .6
}




# Output
out <- rep(n_rep = n_rep, ni, nj, gamma00, gamma01, gamma02, gamma03, gamma04, gamma05, gamma06,
           gamma07, gamma08, gamma09, gamma10, gamma11, gamma12, gamma13, 
           gamma14, u_0, u_1, u_2)

#path0 <- getwd()

#path1 <- paste0(path, "/ratio")

#dir.create(file.path(path1), recursive = TRUE)

saveRDS(out, file = paste0("bhs_", n_iter, "_", n_rep, "_", N, ".rds"), compress = T)

# Copy files and move to the staging folder

if (N = 4624) {
  file.copy(from = paste0("bhs_", n_iter, "_", n_rep, "_", N, ".rds"),
            to = paste0("/staging/mhuang233/bhs_", n_iter, "_", n_rep, "_", N, ".rds"))
  
  file.remove(paste0("bhs_", n_iter, "_", n_rep, "_", N, ".rds"))
}


# Create txt file for argument
# bhs <- data.frame(n_iter = c(rep(1:5, 3), 1:10), n_rep = c(rep(20, 15), rep(10, 10)),
#                   ni = c(rep(10, 5), rep(20, 5), rep(40, 5), rep(150, 10)),
#                   nj = c(rep(10, 5), rep(20, 5), rep(20, 5), rep(35, 10)),
#                   gb = c(rep(15, 5), rep(25, 5), rep(35, 5), rep(50, 10)))

# write.table(bhs, file = "bhs.txt", row.names = F, col.names = F)

# pars <- sapply(out, "[[", "$pars")
#bias <- as.data.frame(sapply(out, "[[", "bias"))
#rel_bias <- as.data.frame(sapply(out, "[[", "rel_bias"))
#kld <- as.data.frame(sapply(out, "[[", "kld"))


# LPD and weights
#lpd_bs <- as.data.frame(sapply(out, "[[", "lpd_bs"))
#lpd_bhs <- as.data.frame(sapply(out, "[[", "lpd_bhs"))
#w_bs_r <- as.data.frame(sapply(out, "[[", "w_bs_r"), row.names = c("m1","m2","m3","m4"))
#w_bhs_r <- as.data.frame(sapply(out, "[[", "w_bhs_m"), row.names = c("m1","m2","m3","m4"))


# parameters
#est_m <- as.data.frame(sapply(out, "[[", "est_m"))
#est_sd <- as.data.frame(sapply(out, "[[", "est_sd"))
#ci_l <- as.data.frame(sapply(out, "[[", "ci_l"))
#ci_u <- as.data.frame(sapply(out, "[[", "ci_u"))
#eff <- as.data.frame(sapply(out, "[[", "eff"))
#rhat <- as.data.frame(sapply(out, "[[", "rhat"))


#names(lpd_bs) <- paste0("V",1:n_rep)
#names(lpd_bhs) <- paste0("V",1:n_rep)
#names(bias) <- paste0("V",1:n_rep)
#names(rel_bias) <- paste0("V",1:n_rep)

#names(est_m) <- paste0("V",1:n_rep)
#names(est_sd) <- paste0("V",1:n_rep)
#names(ci_l) <- paste0("V",1:n_rep)
#names(ci_u) <- paste0("V",1:n_rep)
#names(eff) <- paste0("V",1:n_rep)
#names(rhat) <- paste0("V",1:n_rep)


#saveRDS(bias, file = paste0("bias_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(rel_bias, file = paste0("rel_bias_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(lpd_bs, file = paste0("lpd_bs_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(lpd_bhs, file = paste0("lpd_bhs_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(w_bs_r , file = paste0("w_bs_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(w_bhs_r, file = paste0("w_bhs_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(est_m, file = paste0("est_m_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(est_sd, file = paste0("est_sd_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(ci_l, file = paste0("cil_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(ci_u, file = paste0("ciu_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(eff, file = paste0("eff_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(rhat, file = paste0("rhat_", ni, "_", nj, ".rds"), compress = T)
#saveRDS(kld, file = paste0("kld_", ni, "_", nj, ".rds"), compress = T)


