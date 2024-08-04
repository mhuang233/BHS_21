
### ::: SIMULATION FUNCTIONS ::: ###

### === data generating functions === ###
means <- c(1.5, 0.07, -0.04, 0.14, -0.07, 
           -0.02, 0.28, 0.06, -0.09, 0.14, 
           0.31, -0.15, 0.17, 0.02, 0.24, 
           0.25, 0.18, -0.26)

sds <- c(0.5, 1.02, 1.16, 1.12, 1.07, 
         1, 0.99, 1.02, 1.01, 1.08,
         1.01, 1.02, 0.99, 0.94, 1.03,
         1.04, 0.95, 0.98)

# predictors in order:
# 1-4: Female, ESCS, HOMEPOS, ICTRES
# 5-8: JOYREAD, PISADIFF, SCREADCOMP, SCREADDIFF
# 9-15: METASUM, GFOFAIL, MASTGOAL, SWBP, WORKMAST, ADAPTIVITY, COMPETE
# 16-18:PERFEED, TEACHINT, BELONG

dgf <- function(ni, nj, gamma00, gamma01, gamma02, gamma03, gamma04, u_0, w_0, sigma){
  
  # for random effects and sd
  corU <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = T)
  sdU <- matrix(c(u_0, 0, 0, 0, u_0, 0, 0, 0, u_0), nrow = 3, ncol = 3, byrow = T)
  covU <- sdU%*%corU%*%sdU
  I <- covU
  
  # for fixed effects and sd
  corW <- diag(4)
  sdW <- diag(4)*w_0
  covW <- sdW%*%corW%*%sdW
  J <- covW
  
  # for outer loop or the macro level --> school levels
  for (i in 1:ni){
    
    U[i,1] <- I[1,1]*rnorm(1, 0, .5)
    U[i,2] <- I[1,2]*U[i,1] + I[2,2]*rnorm(1)
    U[i,3] <- I[1,3]*U[i,1] + I[1,2]*U[i,2] + I[3,3]*rnorm(1)
    
    # x[i, 4] <- rbinom(1, 1, .5) # public
    # x[i, 5] <- rnorm(1, .15, 1) # sd
    
    # For inner loop or the micro level --> individual levels
    for (j in 1:nj){
      
      x[j,1] <- J[1,1]*rnorm(1, means[5], sds[5])
      x[j,2] <- J[1,2]*x[j,1] + J[2,2]*rnorm(1, means[6], sds[6]) 
      x[j,3] <- J[1,3]*x[j,1] + J[2,3]*x[j,2] + J[3,3]*rnorm(1, means[7], means[7])
      x[j,4] <- J[1,4]*x[j,1] + J[2,4]*x[j,2] + J[3,4]*x[j,3] + J[4,4]*rnorm(1, means[8], sds[8])
      
    }
    
    ind <- 1
    
    for (i in 1:ni){
      for (j in 1:nj) {
        
        # for response and error
        
        r[ind, 1] <- sigma*rnorm(1, 0, 1)
        y[ind, 1] <- gamma00 + gamma01*x[j,1] + gamma02*x[j,2] + gamma03*x[j,3] + gamma04*x[j,4] + U[i,1] + r[ind, 1]
        
        # pull out the data
        
        tmp <- c(i, j, y[ind, 1], x[j, 1], x[j, 2], x[j, 3], x[j, 4], U[i,1], r[ind, 1])
        
        sim[ind, ] <- tmp
        
        ind <- ind + 1
      }
    }
  }
  
  colnames(sim) <- c("i", "j", "y", "x1", "x2", "x3", "x4", "u_0", "r")
  
  sim <- as.data.frame(sim)
  
  return(sim)
  
}

### === fit the model and obtain the results === ###
# icc <- var_b / (var_b + var_w)
# choose icc from .10, .20, .30 from Heges' range form .1 to .25.

alz <- function(rep, c, ni, nj, sd, df, model_strings, var_names){
  
  for (k in 1:15){
    
    a <- stan_lmer(as.formula(model_strings[k]), data = df, 
                   prior_intercept = student_t(3, 400, 10), 
                   prior_covariance = decov(scale = 0.50), 
                   iter = 10000, adapt_delta=.999, thin=10)
    
    b <- log_lik(a, merge_chains = FALSE) #loo(log_lik(a), r_eff = NA)
    
    assign(paste0("bms_", var_names[k]), a)
    assign(paste0("loglik_", var_names[k]), b)
    
    # check eff and rhat
    e <- a$stan_summary[, c("n_eff", "Rhat")] %>% colMeans()
    assign(paste0("converge_", var_names[k]), e)
   
  }
  
  f <- ls(pattern = "bms_x", all.names = T)
  bms_all <- do.call(list, mget(f))
  
  # all
  d <- ls(pattern = "loglik_x", all.names = T)
  loglik_all <- do.call(list, mget(d))
  print(bms_all)
  
  time_loo <- system.time(loo_bms <- lapply(loglik_all, loo, cores = 4))
  
  time_bs <- time_loo + system.time(
    w_bs <- loo::loo_model_weights(loo_bms, method = "stacking"))
  
  time_pbma <- time_loo + system.time(
    w_pbma <- loo::loo_model_weights(loo_bms, method = "pseudobma", BB = FALSE))
  
  time_pbmabb <- time_loo + system.time(
    w_pbmabb <- loo::loo_model_weights(loo_bms, method = "pseudobma"))
  
  ### for all
  n_draws <- nrow(as.matrix(a))
  
  ypred_bs <- matrix(NA, nrow = n_draws, ncol = nobs(bms_all[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_bs), size = 1, prob = w_bs)
    ypred_bs[d, ] <- posterior_predict(bms_all[[k]], draws = 1)
  }
  
  y_bs <- colMeans(ypred_bs);print(ypred_bs);print(y_bs)
  
  d1 <- density(y_bs, kernel = c("gaussian"))$y
  d0 <- density(df$y, kernel = c("gaussian"))$y
  
  kld1 <- KLD(d1, d0)$sum.KLD.py.px
  
  # pbma
  ypred_pbma <- matrix(NA, nrow = n_draws, ncol = nobs(bms_all[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_pbma), size = 1, prob = w_pbma)
    ypred_pbma[d, ] <- posterior_predict(bms_all[[k]], draws = 1)
  }
  
  y_pbma <- colMeans(ypred_pbma);print(ypred_pbma);print(ypred_pbma)
  
  d2 <- density(y_pbma, kernel = c("gaussian"))$y
  
  kld2 <- KLD(d2, d0)$sum.KLD.py.px
  
  # pbmabb
  ypred_pbmabb <- matrix(NA, nrow = n_draws, ncol = nobs(bms_all[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:length(w_pbmabb), size = 1, prob = w_pbmabb)
    ypred_pbmabb[d, ] <- posterior_predict(bms_all[[k]], draws = 1)
  }
  
  y_pbmabb <- colMeans(ypred_pbmabb);print(ypred_pbmabb);print(ypred_pbmabb)
  
  d3 <- density(y_pbmabb, kernel = c("gaussian"))$y
  
  kld3 <- KLD(d3, d0)$sum.KLD.py.px
 
  # bhs
  # build the model
  X <- df %>% select(x1, x2, x3, x4)
  N <- nrow(X)
  d <- ncol(X)
  lpd_point <- do.call(cbind, lapply(loo_bms, extract_lpd))
  K <- ncol(lpd_point)

  dt_bhs <- list(N = N, d = d, K = K, X = X, 
                 lpd_point = lpd_point, tau_mu = 1, tau_con = 1)

  time_bhs <- system.time(
    fit_bhs <- stan("bhs_con.stan", data = dt_bhs, iter = 10000, save_warmup = FALSE))
  
  # Obtain the weights and the softmax function
  wts_bhs <- rstan::extract(fit_bhs, pars = 'w')$w
  w_bhs <- apply(wts_bhs, c(2,3), mean)
  w_bhs_m <- as.matrix(apply(wts_bhs, 3, mean))
  
  # Obtain the KLD
  ypred_bhs <- matrix(NA, nrow = n_draws, ncol = nobs(bms_all[[1]]))
  for (d in 1:n_draws) {
    k <- sample(1:15, size = 1, prob = w_bhs_m)
    ypred_bhs[d, ] <- posterior_predict(bms_all[[k]], draws = 1)
  }
  
  y_bhs <- colMeans(ypred_bhs);print(ypred_bhs);print(y_bhs)
  
  # KLD
  d4 <- density(y_bhs, kernel = c("gaussian"))$y
  kld4 <- KLD(d4, d0)$sum.KLD.py.px
  
  ws <- data.frame(as.matrix(w_bs), as.matrix(w_pbma), as.matrix(w_pbmabb), w_bhs_m)
  klds <- rbind(kld1, kld2, kld3, kld4) %>% as.data.frame()
  
  time_bs;time_pbma;time_pbmabb;time_bhs
  time <- rbind(t(data.matrix(time_bs)), t(data.matrix(time_pbma)), 
                     t(data.matrix(time_pbmabb)), t(data.matrix(time_bhs))) %>% as.data.frame()
  
  rnames <- c("bs","pbma", "pbmabb", "bhs")  
  rownames(klds) <- rnames
  colnames(ws) <- rnames
  rownames(time) <- rnames
  
  print(klds);print(ws);print(time)

  
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_ws"), ws)
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_kld"), klds)
  assign(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_time"), time)
  
  sw <- paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_ws")
  dlk <- paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_kld")
  emit <- paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_time")

  save(list = c(sw,dlk,emit), file = paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_out.RData"))
  
  file.copy(from = paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_out.RData"),
            to = paste0("/staging/mhuang233/rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_out.RData"))
  
  
  file.remove(paste0("rep_", rep, "_seed", c, "_", ni, "_", nj, "_sd", sd, "_sigma", sigma, "_out.RData"))
  
}

