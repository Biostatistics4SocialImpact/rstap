// General Linear STAP  for a Bernoulli outcome
functions {
#include /functions/common_functions.stan
#include /functions/bernoulli_likelihoods.stan
}
data {
  // dimensions
  int<lower=0> K;        // number of predictors
  int<lower=1> N[2];     // number of observations where y = 0 and y = 1 respectively
  int<lower=1> NN;      // number of total observations
  vector[K] zbar;        // vector of column-means of rbind(Z0, X1)
  matrix[N[1],K] Z0;   // centered (by zbar) predictor matrix | y = 0
  matrix[N[2],K] Z1;   // centered (by zbar) predictor matrix | y = 1
  int y_0[N[1]]; // indices for y == 0
  int y_1[N[2]]; // indices for y == 1

  //STAP data
  int<lower=0> Q; // number of STAP,SAP, and TAPs 
  int<lower=0,upper=1> log_ar[Q]; // 0 = no log; 1 = log
  int<lower=0,upper=2> stap_code[Q]; // 0 = sap ; 1 = tap ; 2 = stap
  int<lower=0> Q_t; // number of taps
  int<lower=0> Q_s; // number of saps
  int<lower=0,upper=Q-Q_t-Q_s> Q_st; // number of staps
  int<lower=0> num_s_wei; // number of weibull s(t)aps
  int<lower=0> num_t_wei; // number of weibull (s)taps
  int<lower=0> M; // Max number of BEF's within inclusion distance
  vector<lower=0>[(Q_s+Q_st) > 0 ? M : 0] dists_crs[Q_s + Q_st]; //distance crs matrix 
  vector<lower=0>[(Q_t+Q_st) > 0 ? M : 0] times_crs[Q_t + Q_st]; // time crs matrix
  int<lower=0> weight_mat[Q,2]; // spatial-temporal weight functions 0 = no component; 1 = erf; 2 = cerf; 3 = exp; 4 = cexp;
  // Meta data
  int u_s[(Q_s + Q_st) > 0 ? NN: 0 , (Q_s + Q_st) > 0 ? ((Q_s + Q_st) * 2) : 0]; //index holder array for distances
  int u_t[(Q_t + Q_st) > 0 ? NN: 0 , (Q_t + Q_st) > 0 ? ((Q_s + Q_st) * 2) : 0]; // index holder array for times
  real<lower=0> max_distance; // Max distance for spatial scales
  real<lower=0> max_time; // Max Time
  
  // declares prior_PD, has_intercept, link, prior_dist, prior_dist_for_intercept
#include /data/data_glm.stan
  
  int<lower=5,upper=5> family;
  
  // weights
  int<lower=0,upper=1> has_weights;  // 0 = No, 1 = Yes
  vector[has_weights ? N[1] : 0] weights0;
  vector[has_weights ? N[2] : 0] weights1;
  
  // offset
  int<lower=0,upper=1> has_offset;  // 0 = No, 1 = Yes
  vector[has_offset ? N[1] : 0] offset0;
  vector[has_offset ? N[2] : 0] offset1;
  
  // declares prior_{mean, scale, df}, prior_{mean, scale, df}_for_intercept, prior_{mean, scale, df}_for_aux
#include /data/hyperparameters.stan
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, {len_}concentration, {len_}regularization
#include /data/glmer_stuff.stan

  // more glmer stuff
  int<lower=0> num_non_zero[2];     // number of non-zero elements in the Z matrices
  vector[num_non_zero[1]] w0;       // non-zero elements in the implicit Z0 matrix
  vector[num_non_zero[2]] w1;       // non-zero elements in the implicit Z1 matrix
  int<lower=0, upper = q - 1> v0[num_non_zero[1]]; // column indices for w0
  int<lower=0, upper = q - 1> v1[num_non_zero[2]]; // column indices for w1
  // where the non-zeros start in each row of Z0
  int<lower=0, upper = rows(w0) + 1> u0[t > 0 ? N[1] + 1 : 0];  
  // where the non-zeros start in each row of Z1
  int<lower=0, upper = rows(w1) + 1> u1[t > 0 ? N[2] + 1 : 0];  
  int<lower=0, upper=1> special_case;     // whether we only have to deal with (1|group)
}
transformed data {
  real aux = not_a_number();
  int<lower=1> V0[special_case ? t : 0,N[1]] = make_V(N[1], special_case ? t : 0, v0);
  int<lower=1> V1[special_case ? t : 0,N[2]] = make_V(N[2], special_case ? t : 0, v1);
  // defines hs, len_z_T, len_var_group, delta, pos
#include /tdata/tdata_glm.stan
}
parameters {
  real<upper=(link == 4 ? 0.0 : positive_infinity())> gamma[has_intercept];
  // declares z_beta, global, local, z_b, z_T, rho, zeta, tau
#include /parameters/parameters_glm.stan
}
transformed parameters {
  // defines beta, b, theta_L
#include /tparameters/tparameters_bernoulli.stan
  if (t > 0) {
    if (special_case) {
      int start = 1;
      theta_L = scale .* tau;
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p, 
                             1.0, tau, scale, zeta, rho, z_T);
      b = make_b(z_b, theta_L, p, l);
    }
  }
}
model {
  // defines eta0, eta1
#include /model/make_eta_bern.stan
  if (has_intercept == 1) {
    if (link != 4) {
      eta0 = gamma[1] + eta0;
      eta1 = gamma[1] + eta1;
    }
    else {
      real shift = fmax(max(eta0), max(eta1));
      eta0 = gamma[1] + eta0 - shift;
      eta1 = gamma[1] + eta1 - shift;
    }
  }
  if (has_weights == 0) {
    real dummy = ll_bern_lp(eta0, eta1, link, N);
  }
  else  {  // weighted log-likelihoods
    target += dot_product(weights0, pw_bern(0, eta0, link));
    target += dot_product(weights1, pw_bern(1, eta1, link));
  }
  
#include /model/priors_glm.stan
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, del, shape, t, p);
}
generated quantities {
  real alpha[has_intercept];
  vector[Q] adj_beta;
  real mean_PPD = 0;
  if (has_intercept == 1) {
     alpha[1] = gamma[1] - dot_product(zbar, delta) - dot_product(colmeans(X),beta ./ colsds(X));
  }
  {
    vector[N[1]] pi0;
    vector[N[2]] pi1;
    // defines eta0, eta1
#include /model/make_eta_bern.stan
    if (has_intercept == 1) {
      if (link != 4) {
        eta0 = gamma[1] + eta0;
        eta1 = gamma[1] + eta1;
      }      
      else {
        real shift;
        shift = fmax(max(eta0), max(eta1));
        eta0 = gamma[1] + eta0 - shift;
        eta1 = gamma[1] + eta1 - shift;
        alpha[1] = alpha[1] - shift;
      }
    }
    
  pi0 = linkinv_bern(eta0, link);
  pi1 = linkinv_bern(eta1, link);
  for (n in 1:N[1]) mean_PPD = mean_PPD + bernoulli_rng(pi0[n]);
  for (n in 1:N[2]) mean_PPD = mean_PPD + bernoulli_rng(pi1[n]);
    
    mean_PPD = mean_PPD / NN;
    adj_beta = beta ./colsds(X);
  }
}
