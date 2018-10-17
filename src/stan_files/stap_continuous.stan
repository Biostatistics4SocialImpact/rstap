// GLM for a Gaussian, Gamma or Inverse Gamma outcome
functions {
#include /functions/common_functions.stan
#include /functions/continuous_likelihoods.stan
}
data{
    // declares N, K, Z, zbar
#include /data/NKZ.stan
    int<lower=0> len_y; //length of y
    real lb_y;
    real<lower=lb_y> ub_y;
    vector<lower=lb_y,upper=ub_y>[len_y] y; //continuous outcome
    int<lower=1, upper = 3> family; // 1 gaussian, 2 gamma,  3 inv-gaussian
    // declares has_intercept, link, prior_dist, prior_dist_for_intercept
#include /data/data_glm.stan
    //declares has_weights, weights, has_offset, offset
#include /data/weights_offset.stan
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_{mean, scale, df}_for_aux
#include /data/hyperparameters.stan
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_} regularization
#include /data/glmer_stuff.stan
  // declares num_not_zero, w, v, u
#include /data/glmer_stuff2.stan
}
transformed data {
  vector[family == 3 ? len_y : 0] sqrt_y;
  vector[family == 3 ? len_y : 0] log_y;
  real sum_log_y = family == 1 ? not_a_number() : sum(log(y));
  int<lower=1> V[special_case ? t: 0, len_y] = make_V(len_y, special_case ? t: 0, v);
#include /tdata/tdata_glm.stan
  is_continuous = 1;
  if (family == 3) {
    sqrt_y = sqrt(y);
    log_y = log(y);
  }
}
parameters{
    real<lower=make_lower(family,link),upper=make_upper(family,link)> gamma[has_intercept];
  // declares z_beta,z_
#include /parameters/parameters_glm.stan
  real<lower=0> aux_unscaled; // interpretation depends on family!
}
transformed parameters {
  // aux has to be defined first in the hs case
  real aux = prior_dist_for_aux == 0 ? aux_unscaled : (prior_dist_for_aux <= 2 ? 
             prior_scale_for_aux * aux_unscaled + prior_mean_for_aux :
             prior_scale_for_aux * aux_unscaled);

  // defines beta, delta, X, X_tilde
#include /tparameters/tparameters_glm.stan
  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* tau * aux;
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p, 
                             aux, tau, scale, zeta, rho, z_T);
      b = make_b(z_b, theta_L, p, l);
    }
  }
}
model {
#include /model/make_eta.stan
  if( t > 0){
  #include /model/eta_add_Wb.stan
  }
  if (has_intercept == 1) {
    if ((family == 1 || link == 2) || (family == 4 && link != 5)) eta = eta + gamma[1];
    else if (family == 4 && link == 5) eta = eta - max(eta) + gamma[1];
    else eta = eta - min(eta) + gamma[1];
  }
  else{
#include /model/eta_no_intercept.stan
  }
  if (prior_dist_for_aux > 0 && prior_scale_for_aux > 0) {
    real log_half = -0.693147180559945286;
    if (prior_dist_for_aux == 1)
      target += normal_lpdf(aux_unscaled | 0, 1) - log_half;
    else if (prior_dist_for_aux == 2)
      target += student_t_lpdf(aux_unscaled | prior_df_for_aux, 0, 1) - log_half;
    else 
     target += exponential_lpdf(aux_unscaled | 1);
  }
  if(has_weights == 0){ // unweighted log-likelihoods
    if (family == 1) {
      if (link == 1) 
        target += normal_lpdf(y | eta, aux);
      else if (link == 2) 
        target += normal_lpdf(y | exp(eta), aux);
      else 
        target += normal_lpdf(y | inv(eta), aux);
    }
    else if (family == 2) {
      target += GammaReg(y, eta, aux, link, sum_log_y);
    }
    else if (family == 3) {
      target += inv_gaussian(y, linkinv_inv_gaussian(eta, link), 
                             aux, sum_log_y, sqrt_y);
    }
    }
    else {
        vector[N] summands;
        if (family == 1) summands = pw_gauss(y, eta, aux, link);
        else if (family == 2) summands = pw_gamma(y, eta, aux, link);
        else if (family == 3) summands = pw_inv_gaussian(y, eta, aux, link, log_y, sqrt_y);
        target += dot_product(weights, summands);
    }

#include /model/priors_glm.stan
  if (t > 0) decov_lp (z_b, z_T, rho, zeta, tau, regularization, del, shape, t, p);
}
generated quantities {
  real alpha[has_intercept];
  vector[Q] adj_beta;
  real mean_PPD = 0;
  if (has_intercept == 1) {
    alpha[1] = gamma[1] - dot_product(zbar, delta) - dot_product(colmeans(X),beta ./ colsds(X) );
  }
  {
    vector[N] eta_z;
#include /model/make_eta.stan
  if( t > 0){
#include /model/eta_add_Wb.stan
  }
    if (has_intercept == 1) {
      if (make_lower(family,link) == negative_infinity() &&
          make_upper(family,link) == positive_infinity()) eta = eta + gamma[1];
      else if (family == 4 && link == 5) {
        real max_eta = max(eta);
        alpha[1] = alpha[1] - max_eta;
        eta = eta - max_eta + gamma[1];
      }
      else {
        real min_eta = min(eta);
        alpha[1] = alpha[1] - min_eta;
        eta = eta - min_eta + gamma[1];
      }
    }
    else {
#include /model/eta_no_intercept.stan
    }
    {
    if (family == 1) {
      if (link > 1) eta = linkinv_gauss(eta, link);
      for (n in 1:len_y) mean_PPD = mean_PPD + normal_rng(eta[n], aux);
    }
    else if (family == 2) {
      if (link > 1) eta = linkinv_gamma(eta, link);
      for (n in 1:len_y) mean_PPD = mean_PPD + gamma_rng(aux, aux / eta[n]);
    }
    else if (family == 3) {
      if (link > 1) eta = linkinv_inv_gaussian(eta, link);
      for (n in 1:len_y) mean_PPD = mean_PPD + inv_gaussian_rng(eta[n], aux);
    }
    mean_PPD = mean_PPD / len_y;
    }
    adj_beta = beta ./ colsds(X); 
  }
}
