
### ::: SIMULATION FUNCTIONS ::: ###
# icc <- var_b / (var_b + var_w)
# choose icc from .10, .15, .20 ###

# .10 = .10/(.10 + .90)
# .15 = .15/(.15 + .85)
# .20 = .20/(.20 + .80)

### Data Generation Functions ###
dgf <- function(ni, nj, gamma00, gamma01, gamma02, gamma03, gamma04, gamma05,
                u_0, u_1, u_2, w_0, w_1, w_2, sigma){
  
  # For random effects and ICC
  corU <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = T)
  sdU <- matrix(c(u_0, 0, 0, 0, u_1, 0, 0, 0, u_2), nrow = 3, ncol = 3, byrow = T)
  covU <- sdU%*%corU%*%sdU
  I <- covU
  
  # For fixed effects and ICC
  corW <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = T)
  sdW <- matrix(c(w_0, 0, 0, 0, w_1, 0, 0, 0, w_2), nrow = 3, ncol = 3, byrow = T)
  covW <- sdW%*%corW%*%sdW
  J <- covW
  
  # For outer loop or the macro level
  for (i in 1:ni){
    
    U[i,1] <- I[1,1]*rnorm(1, 0, .5)
    U[i,2] <- I[1,2]*U[i,1] + I[2,2]*rnorm(1)
    U[i,3] <- I[1,3]*U[i,1] + I[1,2]*U[i,2] + I[3,3]*rnorm(1)
    
    x[i, 4] <- rbinom(1, 1, .5) # public
    x[i, 5] <- rnorm(1, .15, 1) # ict
    
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
        
        # Pull out the data
        
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
  
  # ::: stack all models ::: #
  
  for (k in 1:15){
    
    if (k <= 5){
      as <- stan_lm(y ~ df1[, paste0(cov1[k])], data = df1, 
                    prior = R2(location = .5, what = 'mean'), 
                    iter = 10000, chains = 4,
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
                iter = 10000, chains = 4,
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
  n_draws <- 4000
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
  # Build the model
  d_discrete = 1
  X =  df2[,4:8] # 5 in total 
   
  lpd_point_all <- do.call(cbind, lapply(loos_all, extract_lpd))
  
  dt_bhs_all <- list(X = X, N = nrow(X), d = ncol(X), d_discrete = d_discrete,
                    lpd_point = lpd_point_all, K = ncol(lpd_point_all), tau_mu = 1,
                    tau_sigma = 1, tau_discrete = .5, tau_con = 1)
  
  time_all_bhs <- system.time(
    fit_bhs_all <- stan("bhs_stan.stan", data = dt_bhs_all, chains = 4))

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
  
  # lpd_bhs <- lpd_point*w_bhs_r
  
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
  
  # n_eff and rhat
  #f_all <- ls(pattern = "converge_all_", all.names = T)
  #converge_all <- do.call(rbind, mget(f_all))
  
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


# ::: Linear models ::: #

alz_lm <- function(rep, c, ni, nj, icc, df1, df2, cov1, cov3){

  for (k in 1:nrow(cov1)){
    
    a_lm <- stan_lm(y ~ df1[, paste0(cov1[k])], data = df1, 
                    prior = R2(location = .5, what = 'mean'),
                    iter = 10000, chains = 4,
                    adapt_delta=.999,thin=10)
    
    b_lm <- loo(log_lik(a_lm), r_eff = NA)
    
    assign(paste0("lms_", cov1[k]), a_lm)
    assign(paste0("loo_lm_", cov1[k]), b_lm)
    
    # check eff and rhat
    e_lm <- a_lm$stan_summary[, c("n_eff", "Rhat")] %>% colMeans()
    assign(paste0("converge_lm_", cov1[k]), e_lm)
  }
  
  # lms
  c_lm <- ls(pattern = "lms_x", all.names = T)
  lms <- do.call(list, mget(c_lm))
  
  # loo
  d_lm <- ls(pattern = "loo_lm_x", all.names = T)
  loo_lms <- do.call(list, mget(d_lm))
  
  time_lm_bs <- system.time(w_bs_lm <- loo::loo_model_weights(loo_lms, method = "stacking"))
  time_lm_pbma <- system.time(w_pbma_lm <- loo::loo_model_weights(loo_lms, method = "pseudobma", BB = FALSE))
  time_lm_pbmabb <- system.time(w_pbmabb_lm <- loo::loo_model_weights(loo_lms, method = "pseudobma"))
  
  # For KLD
  n_draws <- 4000
  
  # bs
  ypred_bs_lm <- matrix(NA, nrow = n_draws, ncol = nobs(lms[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_bs_lm), size = 1, prob = w_bs_lm)
    ypred_bs_lm[d, ] <- posterior_predict(lms[[k]], draws = 1)
  }
  
  y_bs_lm <- colMeans(ypred_bs_lm)
  
  d1_lm <- density(y_bs_lm, kernel = c("gaussian"))$y
  d0_lm <- density(df1$y, kernel = c("gaussian"))$y
  
  kld1_lm <- KLD(d1_lm, d0_lm)$sum.KLD.py.px
  
  # pbma
  ypred_pbma_lm <- matrix(NA, nrow = n_draws, ncol = nobs(lms[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_pbma_lm), size = 1, prob = w_pbma_lm)
    ypred_pbma_lm[d, ] <- posterior_predict(lms[[k]], draws = 1)
  }
  
  y_pbma_lm <- colMeans(ypred_pbma_lm)
  
  d2_lm <- density(y_pbma_lm, kernel = c("gaussian"))$y
  
  kld2_lm <- KLD(d2_lm, d0_lm)$sum.KLD.py.px
  
  # pbmabb
  ypred_pbmabb_lm <- matrix(NA, nrow = n_draws, ncol = nobs(lms[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_pbmabb_lm), size = 1, prob = w_pbmabb_lm)
    ypred_pbmabb_lm[d, ] <- posterior_predict(lms[[k]], draws = 1)
  }
  
  y_pbmabb_lm <- colMeans(ypred_pbmabb_lm)
  
  d3_lm <- density(y_pbmabb_lm, kernel = c("gaussian"))$y
  
  kld3_lm <- KLD(d3_lm, d0_lm)$sum.KLD.py.px
  
  # bhs
  # Build the model
  d_discrete = 1
  X = df2[, 4:8]
  lpd_point_lms <- cbind(loo_lms[[1]]$pointwise[, "elpd_loo"],
                         loo_lms[[2]]$pointwise[, "elpd_loo"],
                         loo_lms[[3]]$pointwise[, "elpd_loo"],
                         loo_lms[[4]]$pointwise[, "elpd_loo"],
                         loo_lms[[5]]$pointwise[, "elpd_loo"])
  
  
  dt_bhs_lm <- list(X = X, N = nrow(X), d = ncol(X), d_discrete = d_discrete,
                    lpd_point = lpd_point_lms, K = ncol(lpd_point_lms), tau_mu = 1,
                    tau_sigma = 1, tau_discrete = .5, tau_con = 1)
  
  time_lm_bhs <- system.time(fit_bhs_lm <- stan("bhs_stan.stan", data = dt_bhs_lm, chains = 4, iter = 10000))
  
  
  # Obtain the weights and the softmax function
  wts_bhs_lm <- rstan::extract(fit_bhs_lm, pars = 'w')$w
  w_bhs_lm <- apply(wts_bhs_lm, c(2,3), mean)
  w_bhs_m_lm <- as.matrix(apply(wts_bhs_lm, 3, mean))
  
  # Obtain the KLD
  ypred_bhs_lm <- matrix(NA, nrow = n_draws, ncol = nobs(lms[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:nrow(cov1), size = 1, prob = w_bhs_m_lm)
    ypred_bhs_lm[d, ] <- posterior_predict(lms[[k]], draws = 1)
  }
  
  y_bhs_lm <- colMeans(ypred_bhs_lm)
  
  # lpd_bhs <- lpd_point*w_bhs_r
  
  # KLD
  d4_lm <- density(y_bhs_lm, kernel = c("gaussian"))$y
  kld4_lm <- KLD(d4_lm, d0_lm)$sum.KLD.py.px
  
  ws_lm <- data.frame(as.matrix(w_bs_lm), as.matrix(w_pbma_lm), as.matrix(w_pbmabb_lm), w_bhs_m_lm)
  cnames <- c("bs","pbma", "pbmabb", "bhs")
  
  klds_lm <- rbind(kld1_lm, kld2_lm, kld3_lm, kld4_lm)
  rownames(klds_lm) <- cnames
  colnames(ws_lm) <- cnames
  
  # n_eff and rhat
  f_lm <- ls(pattern = "converge_lm_", all.names = T)
  converge_lm <- do.call(rbind, mget(f_lm))
  
  # time
  time_lm <- data.frame(as.numeric(time_lm_bs), as.numeric(time_lm_pbma), 
                                  as.numeric(time_lm_pbmabb), as.numeric(time_lm_bhs))
  
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_time"), time_lm)
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_ws"), ws_lm)
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_kld"), klds_lm)
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_converge"), converge_lm)
  
  save(list = c(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_ws"), 
                paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_kld"),
                paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_converge"),
                paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_time")),
       file = paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_out.RData"))
  
  file.copy(from = paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_out.RData"),
            to = paste0("/staging/mhuang233/rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_out.RData"))
  
  file.remove(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_lm_out.RData"))
  
  }





alz_s2 <- function(rep, c, ni, nj, icc, df1, df2, cov1, cov3){
  # ::: 2 random slopes model ::: #
  
  for (k in 1:nrow(cov3)){
    
    covx <- paste0(cov3[k, ], collapse = " + ")
    covz <- paste0(cov3[k, ], collapse = "_")
    
    model_string_s2 <- paste0("y ~ ", covx, " + (1 + x4 + x5|i)")
    
    a_s2 <- stan_lmer(as.formula(model_string_s2), data = df1, 
                      prior_intercept = student_t(3, 400, 10), 
                      prior_covariance = decov(scale = 0.50), 
                      iter = 10000, chains = 4,
                      adapt_delta=.999, thin=10)
    
    b_s2 <- loo(log_lik(a_s2), r_eff = NA)
    
    assign(paste0("s2s_", covz), a_s2)
    assign(paste0("loo_s2_", covz), b_s2)
    
    # check eff and rhat
    e_s2 <- a_s2$stan_summary[, c("n_eff", "Rhat")] %>% colMeans()
    assign(paste0("converge_s2_", covz), e_s2)
  }
  
  # model_s2s
  c_s2 <- ls(pattern = "s2s_x", all.names = T)
  s2s <- do.call(list, mget(c_s2))
  
  # loo
  d_s2 <- ls(pattern = "loo_s2_x", all.names = T)
  loo_s2s <- do.call(list, mget(d_s2))
  
  time_s2_bs <- system.time(w_bs_s2 <- loo::loo_model_weights(loo_s2s, method = "stacking"))
  time_s2_pbma <- system.time(w_pbma_s2 <- loo::loo_model_weights(loo_s2s, method = "pseudobma", BB = FALSE))
  time_s2_pbmabb <- system.time(w_pbmabb_s2 <- loo::loo_model_weights(loo_s2s, method = "pseudobma"))
  
  # For KLD
  n_draws <- 4000
  
  # bs
  ypred_bs_s2 <- matrix(NA, nrow = n_draws, ncol = nobs(s2s[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_bs_s2), size = 1, prob = w_bs_s2)
    ypred_bs_s2[d, ] <- posterior_predict(s2s[[k]], draws = 1)
  }
  
  y_bs_s2 <- colMeans(ypred_bs_s2)
  
  d1_s2 <- density(y_bs_s2, kernel = c("gaussian"))$y
  d0_s2 <- density(df1$y, kernel = c("gaussian"))$y
  
  kld1_s2 <- KLD(d1_s2, d0_s2)$sum.KLD.py.px
  
  # pbma
  ypred_pbma_s2 <- matrix(NA, nrow = n_draws, ncol = nobs(s2s[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_pbma_s2), size = 1, prob = w_pbma_s2)
    ypred_pbma_s2[d, ] <- posterior_predict(s2s[[k]], draws = 1)
  }
  
  y_pbma_s2 <- colMeans(ypred_pbma_s2)
  
  d2_s2 <- density(y_pbma_s2, kernel = c("gaussian"))$y
  
  kld2_s2 <- KLD(d2_s2, d0_s2)$sum.KLD.py.px
  
  # pbmabb
  ypred_pbmabb_s2 <- matrix(NA, nrow = n_draws, ncol = nobs(s2s[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_pbmabb_s2), size = 1, prob = w_pbmabb_s2)
    ypred_pbmabb_s2[d, ] <- posterior_predict(s2s[[k]], draws = 1)
  }
  
  y_pbmabb_s2 <- colMeans(ypred_pbmabb_s2)
  
  d3_s2 <- density(y_pbmabb_s2, kernel = c("gaussian"))$y
  
  kld3_s2 <- KLD(d3_s2, d0_s2)$sum.KLD.py.px
  
  # bhs
  # Build the model
  d_discrete = 1
  X =  df2[,4:8] # 5 in total 
  
  extract_lpd <- function(x){
    
    g <- x$pointwise[, "elpd_loo"]
    return(g)
    
  }
  
  lpd_point_s2 <- do.call(cbind, lapply(loo_s2s, extract_lpd))
  
  dt_bhs_s2 <- list(X = X, N = nrow(X), d = ncol(X), d_discrete = d_discrete,
                    lpd_point = lpd_point_s2, K = ncol(lpd_point_s2), tau_mu = 1,
                    tau_sigma = 1, tau_discrete = .5, tau_con = 1)
  
  time_s2_bhs <- system.time(fit_bhs_s2 <- stan("bhs_stan.stan", data = dt_bhs_s2))# 
  
  
  # Obtain the weights and the softmax function
  wts_bhs_s2 <- rstan::extract(fit_bhs_s2, pars = 'w')$w
  w_bhs_s2 <- apply(wts_bhs_s2, c(2,3), mean)
  w_bhs_m_s2 <- as.matrix(apply(wts_bhs_s2, 3, mean))
  
  # Obtain the KLD
  ypred_bhs_s2 <- matrix(NA, nrow = n_draws, ncol = nobs(s2s[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:nrow(cov3), size = 1, prob = w_bhs_m_s2)
    ypred_bhs_s2[d, ] <- posterior_predict(s2s[[k]], draws = 1)
  }
  
  y_bhs_s2 <- colMeans(ypred_bhs_s2)
  
  # lpd_bhs <- lpd_point*w_bhs_r
  
  # KLD
  d4_s2 <- density(y_bhs_s2, kernel = c("gaussian"))$y
  kld4_s2 <- KLD(d4_s2, d0_s2)$sum.KLD.py.px
  
  ws_s2 <- data.frame(as.matrix(w_bs_s2), as.matrix(w_pbma_s2), as.matrix(w_pbmabb_s2), w_bhs_m_s2)
  klds_s2 <- rbind(kld1_s2, kld2_s2, kld3_s2, kld4_s2)
  cnames <- c("bs","pbma", "pbmabb", "bhs")
  
  rownames(klds_s2) <- cnames
  colnames(ws_s2) <- cnames
  
  # n_eff and rhat
  f_s2 <- ls(pattern = "converge_s2_", all.names = T)
  converge_s2 <- do.call(rbind, mget(f_s2))
  
  # time
  time_s2 <- data.frame(as.numeric(time_s2_bs), as.numeric(time_s2_pbma), 
                        as.numeric(time_s2_pbmabb), as.numeric(time_s2_bhs))
  
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_time"), time_s2)
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_ws"), ws_s2)
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_kld"), klds_s2)
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_converge"), converge_s2)
  
  save(list = c(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_ws"), 
                paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_kld"),
                paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_converge"),
                paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_time")),
                file = paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_out.RData"))
  
  file.copy(from = paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_out.RData"),
            to = paste0("/staging/mhuang233/rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_out.RData"))
  
  file.remove(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_", icc, "_s2_out.RData"))
  
}


### create arguments for CHTC ###
# ni <- rep(c(10, 70, 10, 30, 50, 150), each = 3) # clusters
# nj <- rep(c(10, 70, 50, 150, 10, 30), each = 3) # students
# icc <- rep(c(.1, .2, .3), 6)
# c <- sample(1:1000, 100)

# arg <- data.frame(reps = rep(1:100, each = 18),
#                  cs = rep(c, each = 18),
#                  nis = rep(ni, 100),
#                  njs = rep(nj, 100),
#                 iccs = rep(icc, 100))

# write.table(arg, file = "args.txt", row.names = F, col.names = F)



