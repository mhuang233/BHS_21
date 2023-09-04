
// BHS stan for multilevel structure.

// w = \mu_k + sum(\alpha_{mkg}*f)

// \mu_k | \mu_0,\tau_mu \sim N(\mu_0,\tau_mu)
// \alpha_mkg | \sigma_gk \sim N(alpha_gk, \sigma_gk)
// \alpha_gk | \tau_alpha_gk \sim N(0, \tau_alpha_gk)
// \sigma_gk | \tau_sigma_gk \sim N(0, \tau_sigma_gk)

data {
  int<lower=0> n;//number of observations
  int<lower=0> k;//number of models
  int<lower=0> l;//number of schools
  int<lower=1,upper=k> school_id[n];//vector of sample indeces
  int<lower=1,upper=l> model_id[n];//vector of plate indeces
  
  matrix[n, k] exp_lpd_point;
  real<lower = 0> tau_mu;
  real<lower = 0> tau_k;
  real<lower = 0> tau_l;
}

parameters{
  
  vector[k-1] mu_k;
  real<lower=0> mu_0;
  
  vector[k-1] gamma;//vector of sample deviation from the average 
  vector[l] delta;//vector of plate deviation from the average
  
  vector[k-1] alpha_aux;
  vector[k-1] sigma_aux;
  
}

transformed parameters{
  
  vector[k-1] alpha;
  simplex[k] w[n];
  matrix[n, k] f;

  
  
}
