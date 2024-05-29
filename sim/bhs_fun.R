
### ::: SIMULATION FUNCTIONS ::: ###
# icc <- var_b / (var_b + var_w)
# choose icc from .10, .15, .20 ###

# .10 = .10/(.10 + .90)
# .15 = .15/(.15 + .85)
# .20 = .20/(.20 + .80)

### Data Generation Functions ###
dgf <- function(ni, nj, gamma00, gamma01, gamma02, gamma03, gamma04, gamma05,
                u_0, u_1, u_2, w_0, w_1, w_2, sigma){
  
  # For random effects
  corU <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = T)
  sdU <- matrix(c(u_0, 0, 0, 0, u_1, 0, 0, 0, u_2), nrow = 3, ncol = 3, byrow = T)
  covU <- sdU%*%corU%*%sdU
  I <- covU
  
  # For fixed effects
  corW <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = T)
  sdW <- matrix(c(w_0, 0, 0, 0, w_1, 0, 0, 0, w_2), nrow = 3, ncol = 3, byrow = T)
  covW <- sdW%*%corW%*%sdW
  J <- covW
  
  # For outer loop or the macro level
  for (i in 1:ni){
    
    U[i,1] <- I[1,1]*rnorm(1, 0, .5)
    U[i,2] <- I[1,2]*U[i,1] + I[2,2]*rnorm(1)
    U[i,3] <- I[1,3]*U[i,1] + I[1,2]*U[i,2] + I[3,3]*rnorm(1)
    
    x[i, 4] <- rbinom(1, 1, .5) 
    x[i, 5] <- rnorm(1, .15, 1) 
    
    # For inner loop or the micro level
    for (j in 1:nj){
      
      x[j,1] <- J[1,1]*rnorm(1, 0, .5)
      x[j,2] <- J[1,2]*x[j,1] + J[2,2]*rnorm(1, .12, 1) 
      x[j,3] <- J[1,3]*x[j,1] + J[1,2]*x[j,2] + J[3,3]*rnorm(1, -.3, 1)
  
    }
    
    ind <- 1
    
    for (i in 1:ni){
      for (j in 1:nj) {
        
        # For response and error
        r[ind, 1] <- sigma*rnorm(1, 0, 1)
        y[ind, 1] <- gamma00 + gamma01*x[j,1] + gamma02*x[j,2] + gamma03*x[j,3] + gamma04*x[j,4] + gamma05*x[j,5] +
           U[i,1] + U[i,2]*x[i,4] + U[i,3]*x[i,5] + r[ind, 1]
        
        # Pull out the generated data
        tmp <- c(i, j, y[ind, 1], x[j, 1], x[j, 2], x[j, 3], x[j, 4], x[j, 5], U[i,1], U[i, 2], U[i, 3], r[ind, 1])
        
        sim[ind, ] <- tmp
        
        ind <- ind + 1
      }
    }
  }
  
  colnames(sim) <- c("i", "j", "y", "x1", "x2", "x3", "x4", "x5", "u_0", "u_1", "u_2", "r")
  
  sim <- as.data.frame(sim)
  
  return(sim)
  
}


# Fit, Est and Alz
alz_all <- function(rep, c, ni, nj, icc, df1, df2, cov1, cov3){
  
  # ::: stack all the models ::: #
  
  for (k in 1:15){
    
    if (k <= 5){
      as <- stan_lm(y ~ df1[, paste0(cov1[k])], data = df1, 
                    prior = R2(location = .5, what = 'mean'), 
                    iter = iter, chains = 4,
                    adapt_delta=.999,thin=10)
      
      bs <- loo(log_lik(as), r_eff = NA)
      
      assign(paste0("model_all_", cov1[k]), as)
      assign(paste0("loo_all_",cov1[k]), bs)
      
      # check eff and rhat
      es <- as$stan_summary[, c("n_eff", "Rhat")] %>% colMeans()
      assign(paste0("converge_all_", cov1[k]), es)
      
    }else{
      covx <- paste0(cov3[(k-5), ], collapse = " + ")
      covz <- paste0(cov3[(k-5), ], collapse = "_")
      
      model_string_s2 <- paste0("y ~ ", covx, " + (1 + x4 + x5|i)")
      as <- stan_lmer(as.formula(model_string_s2), data = df1, 
                prior_intercept = student_t(3, 400, 10), 
                prior_covariance = decov(scale = 0.50), 
                iter = iter, chains = 4,
                adapt_delta=.999, thin=10)
      
      bs <- loo(log_lik(as), r_eff = NA)
      
      assign(paste0("model_all_", covz), as)
      assign(paste0("loo_all_",covz), bs)
      
      # check eff and rhat
      es <- as$stan_summary[, c("n_eff", "Rhat")] %>% colMeans()
      assign(paste0("converge_all_", covz), es)
    }

  }
  
  # all
  cs <- ls(pattern = "model_all_x", all.names = T)
  m_all <- do.call(list, mget(cs))
  
  # all
  ds <- ls(pattern = "loo_all_x", all.names = T)
  loos_all <- do.call(list, mget(ds))
  
  time_all_bs <- system.time(w_bs_all <- loo::loo_model_weights(loos_all, method = "stacking"))
  time_all_pbma <- system.time(w_pbma_all <- loo::loo_model_weights(loos_all, method = "pseudobma", BB = FALSE))
  time_all_pbmabb <- system.time(w_pbmabb_all <- loo::loo_model_weights(loos_all, method = "pseudobma"))
  
  ### for all
  n_draws <- n_draws # based on # iteration
  ypred_bs_all <- matrix(NA, nrow = n_draws, ncol = nobs(m_all[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_bs_all), size = 1, prob = w_bs_all)
    ypred_bs_all[d, ] <- posterior_predict(m_all[[k]], draws = 1)
  }
  
  y_bs_all <- colMeans(ypred_bs_all)
  
  d1_all <- density(y_bs_all, kernel = c("gaussian"))$y
  d0_all <- density(df1$y, kernel = c("gaussian"))$y
  
  kld1_all <- KLD(d1_all, d0_all)$sum.KLD.py.px
  
  # pbma
  ypred_pbma_all <- matrix(NA, nrow = n_draws, ncol = nobs(m_all[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_pbma_all), size = 1, prob = w_pbma_all)
    ypred_pbma_all[d, ] <- posterior_predict(m_all[[k]], draws = 1)
  }
  
  y_pbma_all <- colMeans(ypred_pbma_all)
  
  d2_all <- density(y_pbma_all, kernel = c("gaussian"))$y
  
  kld2_all <- KLD(d2_all, d0_all)$sum.KLD.py.px
  
  # pbmabb
  ypred_pbmabb_all <- matrix(NA, nrow = n_draws, ncol = nobs(m_all[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_pbmabb_all), size = 1, prob = w_pbmabb_all)
    ypred_pbmabb_all[d, ] <- posterior_predict(m_all[[k]], draws = 1)
  }
  
  y_pbmabb_all <- colMeans(ypred_pbmabb_all)
  
  d3_all <- density(y_pbmabb_all, kernel = c("gaussian"))$y
  
  kld3_all <- KLD(d3_all, d0_all)$sum.KLD.py.px
  
  # bhs
  d_discrete = 1
  X =  df2[,4:8] # 5 in total 
   
  lpd_point_all <- do.call(cbind, lapply(loos_all, extract_lpd))
  
  dt_bhs_all <- list(X = X, N = nrow(X), d = ncol(X), d_discrete = d_discrete,
                    lpd_point = lpd_point_all, K = ncol(lpd_point_all), tau_mu = 1,
                    tau_sigma = 1, tau_discrete = .5, tau_con = 1)
  
  time_all_bhs <- system.time(
    fit_bhs_all <- stan("bhs_con.stan", data = dt_bhs_all, iter= iter, chains = 4))

  # Obtain the weights and the softmax function
  wts_bhs_all <- rstan::extract(fit_bhs_all, pars = 'w')$w
  w_bhs_all <- apply(wts_bhs_all, c(2,3), mean)
  w_bhs_m_all <- as.matrix(apply(wts_bhs_all, 3, mean))
  
  # Obtain the KLD
  ypred_bhs_all <- matrix(NA, nrow = n_draws, ncol = nobs(m_all[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:15, size = 1, prob = w_bhs_m_all)
    ypred_bhs_all[d, ] <- posterior_predict(m_all[[k]], draws = 1)
  }
  
  y_bhs_all <- colMeans(ypred_bhs_all)
  
  # KLD
  d4_all <- density(y_bhs_all, kernel = c("gaussian"))$y
  kld4_all <- KLD(d4_all, d0_all)$sum.KLD.py.px
  
  ws_all <- data.frame(as.matrix(w_bs_all), as.matrix(w_pbma_all), as.matrix(w_pbmabb_all), w_bhs_m_all)
  klds_all <- rbind(kld1_all, kld2_all, kld3_all, kld4_all)
  
  time_all <- data.frame(as.numeric(time_all_bs), as.numeric(time_all_pbma), 
                        as.numeric(time_all_pbmabb), as.numeric(time_all_bhs))
  
  cnames <- c("bs","pbma", "pbmabb", "bhs")  
  rownames(klds_all) <- cnames
  colnames(ws_all) <- cnames
  
  save(ws_all, klds_all, file = paste0(rep, "_", c, "_", ni, "_", nj, "_", icc, "_all_out.RData"))
  
  assign(paste0("rep_", rep, "_", c, "_seed", ni, "_", nj, "_", icc, "_all_ws"), ws_all)
  assign(paste0("rep_", rep, "_", c, "_seed", ni, "_", nj, "_", icc, "_all_kld"), klds_all)
  assign(paste0("rep_", rep, "_", c, "_seed", ni, "_", nj, "_", icc, "_all_time"), time_all)
  
  save(list = c(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_all_ws"), 
                paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_all_kld"),
                paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_all_time")),
       file = paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_all_out.RData"))
  
  file.copy(from = paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_all_out.RData"),
            to = paste0("/staging/mhuang233/rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_all_out.RData"))
  
  file.remove(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_all_out.RData"))

}

### create arguments for CHTC ###
# ni <- rep(c(10, 70, 10, 30, 50, 150), each = 3) # clusters
# nj <- rep(c(10, 70, 50, 150, 10, 30), each = 3) # students
# icc <- rep(c(.1, .15, .2), 6)
# c <- sample(1:1000, 100)

# arg <- data.frame(reps = rep(1:100, each = 18),
#                  cs = rep(c, each = 18),
#                  nis = rep(ni, 100),
#                  njs = rep(nj, 100),
#                 iccs = rep(icc, 100))

# write.table(arg, file = "args.txt", row.names = F, col.names = F)



