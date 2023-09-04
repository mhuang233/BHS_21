
### ::: Real data analysis ::: ###

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


df0 <- read.csv("pisa2018.BayesBook.csv")

df <- df0 %>%
  dplyr::select(SchoolID, CNTSTUID, Female, ESCS, METASUM, PERFEED, HOMEPOS, 
                ADAPTIVITY, TEACHINT, ICTRES, ATTLNACT, COMPETE, 
                WORKMAST, GFOFAIL, SWBP, MASTGOAL, BELONG, SCREADCOMP, 
                PISADIFF, Public, PV1MATH)

### Subset
sch <- table(df$SchoolID)
dt0  <- subset(df, SchoolID %in% names(sch[sch > 20]))

# Check
library(tidyverse)

dt0 %>%
  group_by(SchoolID)%>%
  summarise(n=n()) # 149 in total

# Randomly select 20 students in each group
dt <- dt0 %>% group_by(SchoolID) %>% slice_sample(n = 20)

# For the reduced dataset
df1 <- dt[1:200, ]
df2 <- dt[201:400, ]

# For the full dataset
nr <- nrow(dt)
nr_half <- nr/2
df1 <- dt[1:nr_half, ]
df2 <- dt[(nr_half+1):nr, ]


# Set up
{
  bs_r <- list()
  bhs_r <- list()
  pars <- list()
  loo_bs_r <- list()
  # df1 <- df[1:(nrow(df)/2), ] # split into half
  # df2 <- df[(nrow(df)/2+1):nrow(df), ]
  
  # Softmax function
  softmax <- function (x) {
    x_tilde <- c(x, 0)
    return (exp(x_tilde) / sum(exp(x_tilde)))
  }
}



bs_r[[1]] <- stan_lmer(
  PV1MATH ~ Public + ESCS + HOMEPOS + (1 |SchoolID), data = df1, 
  prior_intercept = student_t(3, 400, 10),
  prior_covariance = decov(scale = 0.50),
  iter = 5000, chains = 4, QR=TRUE,
  adapt_delta=.999,thin=10)

bs_r[[2]] <- stan_lmer(
  PV1MATH ~ PISADIFF + GFOFAIL + BELONG + TEACHINT + (1 + TEACHINT |SchoolID), data = df1, 
  prior_intercept = student_t(3, 400, 10), 
  prior_covariance = decov(scale = 0.50),
  iter = 5000, chains = 4, QR=TRUE,
  adapt_delta=.999,thin=10)

bs_r[[3]] <- stan_lmer(
  PV1MATH ~  COMPETE + ADAPTIVITY + METASUM + (1 |SchoolID), data = df1, 
  prior_intercept = student_t(3, 400, 10),
  prior_covariance = decov(scale = 0.50),
  iter = 5000, chains = 4, QR=TRUE,
  adapt_delta=.999,thin=10)

bs_r[[4]] <- stan_lmer(
  PV1MATH ~ MASTGOAL + SWBP + WORKMAST + ICTRES + (1 + ICTRES |SchoolID), data = df1, 
  prior_intercept = student_t(3, 400, 10),
  prior_covariance = decov(scale = 0.50),
  iter = 5000, chains = 4, QR=TRUE,
  adapt_delta=.999,thin=10)

#bs_r[[5]] <- stan_lmer(
#  y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + (1 |i) + (0 + x13 |i) + (0 + x14 |i), data = df1, 
#  prior_intercept = student_t(3, 400, 10),
#  prior_covariance = decov(scale = 0.50),
#  iter = 10000, chains = 4, QR=TRUE,
#  adapt_delta=.999,thin=10)

# SAVE the parameter, rhat and ess

r1 <- c("(Intercept)", "x1", "x2", "x3", "sigma", "Sigma[i:(Intercept),(Intercept)]")
r2 <- c("(Intercept)", "x4", "x5", "x6", "x13", "sigma", "Sigma[i:(Intercept),(Intercept)]", "Sigma[i:x13,x13]")
r3 <- c("(Intercept)", "x7", "x8", "x9", "sigma", "Sigma[i:(Intercept),(Intercept)]")
r4 <- c("(Intercept)", "x10", "x11", "x12", "x14", "sigma", "Sigma[i:(Intercept),(Intercept)]", "Sigma[i:x14,x14]")
cs <- c("n_eff", "Rhat")

r1 <- c(1, 2, 3, 4, 11, 12)
r2 <- c(1, 2, 3, 4, 5, 18, 19, 21)
r3 <- c(1, 2, 3, 4, 11, 12)
r4 <- c(1, 2, 3, 4, 5, 18, 19, 21)
cs <- c(11, 12)


par1 <- as.data.frame(bs_r[[1]]$stan_summary[r1, cs])
par2 <- as.data.frame(bs_r[[2]]$stan_summary[r2, cs])
par3 <- as.data.frame(bs_r[[3]]$stan_summary[r3, cs])
par4 <- as.data.frame(bs_r[[4]]$stan_summary[r4, cs])


rownames(par1) <- c("1Intercept", "x1", "x2", "x3", "1_sigma", "1_tau0")
rownames(par2) <- c("2Intercept", "x4", "x5", "x6", "x13", "2_sigma", "2_tau0", "2_tau13")
rownames(par3) <- c("3Intercept", "x7", "x8", "x9" ,"3_sigma", "3_tau0")
rownames(par4) <- c("4Intercept", "x10", "x11", "x12", "x14","4_sigma", "4_tau0", "4_tau14")

pars <- do.call("rbind", list(par1, par2, par3, par4));print(pars)


# true_val <- data.frame(t = c(400, 14, 10,  9, 10, .4, 
#                             400, -20, 6, -2, 8, 10, .4, .2,
#                             400, 7, 1, 20, 10, .4,
#                             400, -14, -1, 7, -4, 10, .4, .6))


# bias <- pars[,1] - true_val

# rownames(bias) <- c("1Intercept", "x1", "x2", "x3", "1_sigma", "1_tau0",
#                    "2Intercept", "x4", "x5", "x6", "x13", "2_sigma", "2_tau0", "2_tau13",
#                    "3Intercept", "x7", "x8", "x9" ,"3_sigma", "3_tau0",
#                    "4Intercept", "x10", "x11", "x12", "x14","4_sigma", "4_tau0", "4_tau14")


#rel_bias = bias/true_val
#bias_perc <- as.data.frame(cbind(bias,rel_bias))
# colnames(bias_perc) <- c("bias", "rel_bias")

bias <- data.frame(m1 = bs_r[[1]]$linear.predictors - bs_r[[1]]$y,
                   m2 = bs_r[[2]]$linear.predictors - bs_r[[2]]$y,
                   m3 = bs_r[[3]]$linear.predictors - bs_r[[3]]$y,
                   m4 = bs_r[[4]]$linear.predictors - bs_r[[4]]$y)

rel_bias <- colMeans(bias/bs_r[[1]]$y)

### LOO and weights
loo_bs_r[[1]] <- loo(log_lik(bs_r[[1]]))
loo_bs_r[[2]] <- loo(log_lik(bs_r[[2]]))
loo_bs_r[[3]] <- loo(log_lik(bs_r[[3]]))
loo_bs_r[[4]] <- loo(log_lik(bs_r[[4]]))

w_bs <- loo_model_weights(loo_bs_r, method = "stacking")
w_pbma <- loo_model_weights(loo_bs_r, method = "pseudobma", BB=FALSE)
w_pbmabb <- loo_model_weights(loo_bs_r, method = "pseudobma")


# Obtain the LPD
lpd_point <- as.matrix(cbind(loo_bs_r[[1]]$pointwise[, "elpd_loo"],
                             loo_bs_r[[2]]$pointwise[, "elpd_loo"],
                             loo_bs_r[[3]]$pointwise[, "elpd_loo"],
                             loo_bs_r[[4]]$pointwise[, "elpd_loo"]))


# For KLD
# bs
n_draws <- nrow(as.matrix(bs_r[[1]]));print(n_draws)
ypred_bs <- matrix(NA, nrow = n_draws, ncol = nobs(bs_r[[1]]))
for (d in 1:n_draws) {
  k <- sample(1:length(w_bs), size = 1, prob = w_bs)
  ypred_bs[d, ] <- posterior_predict(bs_r[[k]], draws = 1)
}

y_bs <- colMeans(ypred_bs)

d1 <- density(y_bs, kernel = c("gaussian"))$y
d0 <- density(df1$PV1MATH, kernel = c("gaussian"))$y

kld1 <- KLD(d1, d0)$sum.KLD.py.px

# pbma

#n_draws <- nrow(as.matrix(bs_r[[1]]));print(n_draws)
ypred_bma <- matrix(NA, nrow = n_draws, ncol = nobs(bs_r[[1]]))
for (d in 1:n_draws) {
  k <- sample(1:length(w_pbma), size = 1, prob = w_pbma)
  ypred_bma[d, ] <- posterior_predict(bs_r[[k]], draws = 1)
}

y_bma <- colMeans(ypred_bma)

d1 <- density(y_bma, kernel = c("gaussian"))$y
d0 <- density(df1$PV1MATH, kernel = c("gaussian"))$y

kld2 <- KLD(d1, d0)$sum.KLD.py.px


# pbmabb
#n_draws <- nrow(as.matrix(bs_r[[1]]));print(n_draws)
ypred_bmabb <- matrix(NA, nrow = n_draws, ncol = nobs(bs_r[[1]]))
for (d in 1:n_draws) {
  k <- sample(1:length(w_pbmabb), size = 1, prob = w_pbmabb)
  ypred_bmabb[d, ] <- posterior_predict(bs_r[[k]], draws = 1)
}

y_bmabb <- colMeans(ypred_bmabb)

d1 <- density(y_bmabb, kernel = c("gaussian"))$y
d0 <- density(df1$PV1MATH, kernel = c("gaussian"))$y

kld3 <- KLD(d1, d0)$sum.KLD.py.px


### ::: BHS ::: ###

# Build the model
d_discrete = 1
X =  df2[,c("ESCS", "HOMEPOS", "TEACHINT", "ICTRES", "PISADIFF",
            "METASUM", "GFOFAIL", "MASTGOAL", "SWBP", "WORKMAST", 
            "ADAPTIVITY", "Public", "COMPETE", "BELONG")] # 14 in total 

stan_bsr <- list(X = X, N = nrow(X), d = ncol(X), d_discrete = d_discrete,
                 lpd_point = lpd_point, K = ncol(lpd_point), tau_mu = 1,
                 tau_sigma = 1, tau_discrete = .5, tau_con = 1)

fit_bhs_r <- stan("bhs_stan.stan", data = stan_bsr, chains = 4, iter = 5000)


# Obtain the weights and the softmax function
wts_bhs_r <- rstan::extract(fit_bhs_r, pars = 'w')$w
w_bhs_r <- apply(wts_bhs_r, c(2,3), mean)
w_bhs_m <- as.matrix(apply(wts_bhs_r, 3, mean))

# wsoft_r <- softmax(w_bhs_r[-length(w_bhs_r)])


# Obtain the KLD
ypred_bhs_r <- matrix(NA, nrow = n_draws, ncol = nobs(bs_r[[1]]))
for (d in 1:n_draws) {
  k <- sample(1:4, size = 1, prob = w_bhs_m)
  ypred_bhs_r[d, ] <- posterior_predict(bs_r[[k]], draws = 1)
}

y_bhs_r <- colMeans(ypred_bhs_r)

# lpd_bhs <- lpd_point*w_bhs_r

# KLD
d2 <- density(y_bhs_r, kernel = c("gaussian"))$y
d0 <- density(df1$PV1MATH, kernel = c("gaussian"))$y
kld4 <- KLD(d2, d0)$sum.KLD.py.px


# summarize the weights and lpd

wr <- data.frame(as.matrix(w_bs), as.matrix(w_pbma), as.matrix(w_pbmabb), w_bhs_m)
colnames(wr) <- c("bs","pbma", "pbmabb", "bhs")

lpd_bhs <- lpd_point*w_bhs_r
lpd_bhs_m <- as.matrix(apply(lpd_bhs, 2, mean))
lpd_m <- as.matrix(apply(lpd_point, 2, mean))

bhs_log <- rstan::extract(fit_bhs_r, pars = 'log_lik')$log_lik
lpd_bs_m <- lpd_m*w_bs

lpd <- t(data.frame(bs = colMeans(lpd_m*wr[,1]), 
                    pbma = colMeans(lpd_m*wr[,2]),
                    pbmabb = colMeans(lpd_m*wr[,3]),
                    bhs = colMeans(lpd_bhs_m)))


#kld sum
kld <- rbind(kld1, kld2, kld3, kld4)


# est_m <- pars[,1]
# est_sd <- pars[,2]
# ci_l <- pars[, 3]
# ci_u <- pars[, 4]
eff <- pars$n_eff
rhat <- pars$Rhat




KLD(df1$PV1MATH, y_bhs_r)$sum.KLD.py.px

KLD(df1$PV1MATH, y_bmabb)$sum.KLD.py.px

KLD(df1$PV1MATH, y_bma)$sum.KLD.py.px

KLD(df1$PV1MATH, y_bs)$sum.KLD.py.px


