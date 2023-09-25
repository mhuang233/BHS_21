
### Combine and summarize the results for CHTC ###

### :::::: Summarize the results :::::: ###
{
  library(data.table)
  library(rstanarm)
  library(rstan)
  library(parallel)
  library(xtable)
  library(ggplot2)
  library(ggthemes)
  library(kableExtra)
  library(LaplacesDemon)
  library(abind)
  library(dplyr)
  library(stringr)
  library(tidyverse)
}

### load the results
est <- c("rhat", "eff", "kld", "lpd", "rel_bias", "wr")

alz <- function(N, n_iter, est){
  a <- list.files(pattern = paste0("_", N, ".rds"), full.names=T)
  b <- lapply(a, readRDS)
  
  if (est == "wr"){
    
    for(k in 1:n_iter){
      c <- as.data.frame(lapply(b[[k]], "[[", est))
      assign(paste0("bhs_", N, "_", k, "_", est), c)
    }
    
    d <- mget(ls(pattern = paste0("_", est)))%>%
      abind(along = 3)%>%
      apply(., 1:2, mean)
    

  }else{
    
    for(k in 1:n_iter){
      c <- as.data.frame(sapply(b[[k]], "[[", est))%>%
        rowMeans()
      assign(paste0("bhs_", N, "_", k, "_", est), c)
    }
    
    d <- do.call(rbind, mget(ls(pattern = paste0("_", est))))%>%
      colMeans()
    
  }
  
  return(d)
  # save(d, file = paste0("out_", N, "_", est, ".RData"))
}

### For different n
N <- c(500, 4500, 4900)
n_iter <- 499

for (k in 1:length(est)){
  
  for (h in 1:length(N)){
      
      f <- alz(N[h], n_iter, est[k])
      assign(paste0("out_", N[h], "_", est[k]), f, envir = .GlobalEnv)

    
  }
  
  e <- do.call(cbind, mget(ls(pattern = paste0("_", est[k]))))
  assign(paste0("out_", est[k]), e, envir = .GlobalEnv)
  
}

# Weights
colnames(out_wr) <- c("bs_100", "pbma_100", "pbmabb_100", "bhs_100",
                      "bs_400", "pbma_400", "pbmabb_400", "bhs_400",
                      "bs_4500", "pbma_4500", "pbmabb_4500", "bhs_4500",
                      "bs_4624", "pbma_4624", "pbmabb_4624", "bhs_4624",
                      "bs_800", "pbma_800", "pbmabb_800", "bhs_800")

out_wts <- out_wr[, c(1, 2, 3, 4,
                      5, 6, 7, 8,
                      17, 18, 19, 20,
                      9, 10, 11, 12, 
                      13, 14, 15, 16)]
xtable(out_wts)

# Predictive performance
out_kld <- out_kld[, c(1, 2, 5, 3, 4)]
rownames(out_kld) <- c("bs", "pbma", "pbmabb", "bhs")

out_lpd <- out_lpd[, c(1, 2, 5, 3, 4)]
rownames(out_lpd) <- c("bs", "pbma", "pbmabb", "bhs")

xtable(out_kld)

xtable(out_lpd)

# Graph

kld <- data.frame(group = rep(c("BS", "BHS"), each = 4),
                  size = rep(c(100, 400, 800, 4500), 2),
                  value = c(out_kld[1, ], out_kld[2, ]))

lpd <- data.frame(group = rep(c("BS", "BHS"), each = 4),
                  size = rep(c(100, 400, 800, 4500), 2),
                  value = c(out_lpd[1, ], out_lpd[2, ]))

# Plots
# kld
k <- ggplot(kld, aes(x=size, y=value, group=group)) +
  geom_line(aes(linetype=group))+
  geom_point(aes(shape=group), size=3) +
  xlab("Sample Size") + ylab("KLDs")

k + theme_bw()

# lpd
l <- ggplot(lpd, aes(x=size, y=value, group=group)) +
  geom_line(aes(linetype=group))+
  geom_point(aes(shape=group), size=3) +
  xlab("Sample Size") + ylab("LPDs")

l + theme_bw()

# save pic as 570 and 350

# set seed

# seed <- sample(1:1000, 500)
ratio <- data.frame(n_iter = rep(1:500, 3),
                    n_rep = rep(1, 4500),
                    ni = rep(c(10, 30, 70), each = 500),
                    nj = rep(c(50, 150, 70), each = 500),
                    seed = rep(seed, 3), 
                    gb = rep(c(5, 20, 20), each = 500))

ratio <- data.frame(n_iter = 1:500,
                    n_rep = rep(1, 500),
                    ni = rep(50, 500),
                    nj = rep(10, 500),
                    seed = seed, 
                    gb = rep(8, 500))

write.table(ratio, "ratio.txt", row.names = F, col.names = F)




