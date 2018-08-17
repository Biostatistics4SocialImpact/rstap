
// General Linear STAP for a binomial outcome
functions {
#include /functions/common_functions.stan
#include /functions/binomial_likelihoods.stan
}
data {
  // declares N, K, Z, Q, zbar 
#include /data/NKZ.stan
  int<lower=0> y[N];         // outcome: number of successes
  int<lower=0> trials[N];    // number of trials
  // declares has_intercept, link, prior_dist, prior_dist_for_intercept  
#include /data/data_glm.stan
  // declares has_weights, weights, has_offset, offset 
#include /data/weights_offset.stan
  int<lower=5,upper=5> family;
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_scale_{mean, scale, df}_for_aux
#include /data/hyperparameters.stan
}
transformed data {
  real aux = not_a_number();
}
parameters {
    real<upper=(link == 4 ? 0.0 : positive_infinity())> gamma[has_intercept];
  // declares z_beta,z_delta,theta_s,theta_t,X, X_tilde
#include /parameters/parameters_glm.stan
}
transformed parameters {
  // defines beta, delta 
#include /tparameters/tparameters_glm.stan
}
model {
#include /model/make_eta.stan
  if (has_intercept == 1 ) {
    if (link != 4) eta = eta + gamma[1];
    else eta = gamma[1] + eta - max(eta);
  }
  else{
#include /model/eta_no_intercept.stan
  }
  // Log-likelihood 
  if (has_weights == 0) {  // unweighted log-likelihoods
    real dummy;  // irrelevant but useful for testing
    dummy = ll_binom_lp(y, trials, eta, link);
  }
  else if (has_weights == 1) 
    target += dot_product(weights, pw_binom(y, trials, eta, link));
#include /model/priors_glm.stan
}
generated quantities {
  real alpha[has_intercept];
  real mean_PPD = 0;
  if (has_intercept == 1) {
    alpha[1] = gamma[1] - dot_product(zbar, delta);
  }
  {
    vector[N] pi;
#include /model/make_eta.stan
    if (has_intercept == 1) {
      if (link != 4) eta = eta + gamma[1];
      else {
        real shift = max(eta);
        eta = gamma[1] + eta - shift;
        alpha[1] = alpha[1] - shift;
      }
    }
    else {
#include /model/eta_no_intercept.stan
    }
    
    pi = linkinv_binom(eta, link);
    for (n in 1:N) mean_PPD = mean_PPD + binomial_rng(trials[n], pi[n]);
    mean_PPD = mean_PPD / N;
  }
}
