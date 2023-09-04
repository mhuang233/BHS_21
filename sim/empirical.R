### ::: Empirical study with 4 models ::: ###
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

# set up
df0 <- read.csv("pisa2018.BayesBook.csv")

df <- df0 %>%
  mutate(StdID = 1:nrow(df0)) %>% 
  dplyr::select(SchoolID, StdID, PV1READ, Female, ESCS, HOMEPOS, ICTRES,
                JOYREAD, PISADIFF, SCREADCOMP, SCREADDIFF,
                METASUM, GFOFAIL, MASTGOAL, SWBP, WORKMAST, ADAPTIVITY, COMPETE,
                PERFEED, TEACHINT, BELONG)

sch <- table(df$SchoolID);sch # check the number of students in each school.

df1 <- df[1:2419,] # training
df2 <- df[2420:4838,] # testing

extract_lpd <- function(x){
  
  g <- x$pointwise[, "elpd_loo"]
  return(g)
  
}

# Fit the models
bms <- list()

bms[[1]] <- stan_lmer(
  PV1READ ~ Female + ESCS + HOMEPOS + ICTRES + (1 + ICTRES|SchoolID), data = df1, 
  prior_intercept = student_t(3, 470, 100),
  iter = 5000, chains = 4,
  adapt_delta=.999,thin=10)


bms[[2]] <- stan_lmer(
  PV1READ ~ JOYREAD + PISADIFF + SCREADCOMP + SCREADDIFF + (1 |SchoolID),
  data = df1, prior_intercept = student_t(3, 470, 100),iter = 5000, chains = 4,
  adapt_delta=.999,thin=10)


bms[[3]] <- stan_lmer(
  PV1READ ~ METASUM + GFOFAIL + MASTGOAL + SWBP + WORKMAST + ADAPTIVITY + COMPETE + (1 |SchoolID),
  data = df1, prior_intercept = student_t(3, 470, 100),iter = 5000, chains = 4,
  adapt_delta=.999,thin=10)


bms[[4]] <- stan_lmer(
  PV1READ ~  PERFEED + BELONG + TEACHINT + (1 + TEACHINT|SchoolID),
  data = df1, prior_intercept = student_t(3, 470, 100),iter = 5000, chains = 4,
  adapt_delta=.999,thin=10)


# compute loo
loo_bms <- lapply(bms, loo, cores = 4)

# compute stacking weights
w_bs <- loo_model_weights(loo_bms, method = "stacking")
w_pbma <- loo_model_weights(loo_bms, method = "pseudobma", BB = F)
w_pbmabb <- loo_model_weights(loo_bms, method = "pseudobma", BB = T)

# compute the log density point
lpd_points <- do.call(cbind, lapply(loo_bms, extract_lpd))

# BHS
d_discrete <- 1

X <- df2 %>% 
  select(-c(SchoolID, StdID, PV1READ))

dt_bhs <- list(X = X, N = nrow(X), d = ncol(X), d_discrete = d_discrete,
                  lpd_point = lpd_points, K = ncol(lpd_points), tau_mu = 1,
                  tau_sigma = 1, tau_discrete = .5, tau_con = 1)

fit_bhs <- stan("bhs_discon.stan", data = dt_bhs, iter = 5000)

wts_bhs <- rstan::extract(fit_bhs, pars = 'w')$w
w_bhs <- apply(wts_bhs, c(2,3), mean)
w_bhs_m <- as.matrix(apply(wts_bhs, 3, mean))

### ::: KLDs ::: ###
n_draws <- nrow(as.matrix(bms[[1]])) 

# bs
ypred_bs <- matrix(NA, nrow = n_draws, ncol = nobs(bms[[1]]))
for (d in 1:n_draws) {
  k <- sample(1:length(w_bs), size = 1, prob = w_bs)
  ypred_bs[d, ] <- posterior_predict(bms[[k]], draws = 1)
}

y_bs <- colMeans(ypred_bs)

d1 <- density(y_bs, kernel = c("gaussian"))$y
d0 <- density(df2$PV1READ, kernel = c("gaussian"))$y

kld1 <- KLD(d1, d0)$sum.KLD.py.px

# pbma
ypred_pbma <- matrix(NA, nrow = n_draws, ncol = nobs(bms[[1]]))
for (d in 1:n_draws) {
  k <- sample(1:length(w_pbma), size = 1, prob = w_pbma)
  ypred_pbma[d, ] <- posterior_predict(bms[[k]], draws = 1)
}

y_pbma <- colMeans(ypred_pbma)

d2 <- density(y_pbma, kernel = c("gaussian"))$y

kld2 <- KLD(d2, d0)$sum.KLD.py.px

# pbmabb
ypred_pbmabb <- matrix(NA, nrow = n_draws, ncol = nobs(bms[[1]]))
for (d in 1:n_draws) {
  k <- sample(1:length(w_pbmabb), size = 1, prob = w_pbmabb)
  ypred_pbmabb[d, ] <- posterior_predict(bms[[k]], draws = 1)
}

y_pbmabb <- colMeans(ypred_pbmabb)

d3 <- density(y_pbmabb, kernel = c("gaussian"))$y

kld3 <- KLD(d3, d0)$sum.KLD.py.px

# bhs
ypred_bhs <- matrix(NA, nrow = n_draws, ncol = nobs(bms[[1]]))
for (d in 1:n_draws) {
  k <- sample(1:4, size = 1, prob = w_bhs_m)
  ypred_bhs[d, ] <- posterior_predict(bms[[k]], draws = 1)
} 

y_bhs <- colMeans(ypred_bhs)

d4 <- density(y_bhs, kernel = c("gaussian"))$y
kld4 <- KLD(d4, d0)$sum.KLD.py.px


# results
d4 <- density(y_bhs, kernel = c("gaussian"))$y
kld4 <- KLD(d4, d0)$sum.KLD.py.px

ws <- data.frame(as.matrix(w_bs), as.matrix(w_pbma), as.matrix(w_pbmabb), w_bhs_m) %>% 
  round(3)
klds <- rbind(kld1, kld2, kld3, kld4)
cnames <- c("bs","pbma", "pbmabb", "bhs")

rownames(klds) <- cnames
colnames(ws) <- cnames

# check the descriptive statistics for simulation
library(psych)

means <- describe(X) %>% select(mean)
sds <- describe(X) %>% select(sd)




