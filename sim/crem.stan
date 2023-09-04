/*A simple example of an crossed hierarchical model
*based on the Penicillin data from the lme4 package
*/

data {
  int<lower=0> N;//number of observations
  int<lower=0> n_sample;//number of samples
  int<lower=0> n_plate;//number of plates
  int<lower=1,upper=n_sample> sample_id[N];//vector of sample indeces
  int<lower=1,upper=n_plate> plate_id[N];//vector of plate indeces
  vector[N] y;
}
parameters {
  vector[n_sample] gamma;//vector of sample deviation from the average 
  vector[n_plate] delta;//vector of plate deviation from the average
  real<lower=0> mu;//average diameter value
  real<lower=0> sigma_gamma;//standard deviation of the gamma coeffs
  real<lower=0> sigma_delta;//standard deviation of the delta coeffs
  real<lower=0> sigma_y;//standard deviation of the observations
}
transformed parameters {
  vector[N] y_hat;

  for (i in 1:N)
    y_hat[i] = mu + gamma[sample_id[i]] + delta[plate_id[i]];
}
model {
  //prior on the scale coefficient
//weakly informative priors, see section 6.9 in STAN user guide
  sigma_gamma ~ cauchy(0,2.5);
  sigma_delta ~ cauchy(0,2.5);
  sigma_y ~ gamma(2,0.1);
  //get sample and plate level deviation
  gamma ~ normal(0, sigma_gamma);
  delta ~ normal(0, sigma_delta);
  //likelihood
  y ~ normal(y_hat, sigma_y);
}
generated quantities {
//sample predicted values from the model for posterior predictive checks
  real y_rep[N];
  for(n in 1:N)
    y_rep[n] = normal_rng(y_hat[n],sigma_y);
}