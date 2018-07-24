  // Log-priors for coefficients
       if (prior_dist == 1) target += normal_lpdf(z_delta | 0, 1);
  else if (prior_dist == 2) target += normal_lpdf(z_delta | 0, 1); // Student t via Cornish-Fisher expansion
  else if (prior_dist == 3) { // hs
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_delta | 0, 1);
  }
  else if (prior_dist == 4) { // hs+
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_delta | 0, 1);
  }
  else if (prior_dist == 5) { // laplace
    target += normal_lpdf(z_delta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
  }
  else if (prior_dist == 6) { // lasso
    target += normal_lpdf(z_delta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
    target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
  }
  else if (prior_dist == 7) { // product_normal
    target += normal_lpdf(z_delta | 0, 1);
  }
  /* else prior_dist is 0 and nothing is added */
  
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1)  // normal
      target += normal_lpdf(gamma | prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2)  // student_t
      target += student_t_lpdf(gamma | prior_df_for_intercept, prior_mean_for_intercept, 
                               prior_scale_for_intercept);
    /* else prior_dist is 0 and nothing is added */
  }

  // Log-priors for stap_coefficients 
       if (prior_dist_for_stap == 1) target += normal_lpdf(z_beta | 0, 1);
  else if (prior_dist_for_stap == 2) target += normal_lpdf(z_beta | 0, 1); // Student t via Cornish-Fisher expansion
  else if (prior_dist_for_stap == 3) { // hs
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
  }
  else if (prior_dist_for_stap == 4) { // hs+
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
  }
  else if (prior_dist_for_stap == 5) { // laplace
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix_stap[1] | 1);
  }
  else if (prior_dist_for_stap == 6) { // lasso
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix_stap[1] | 1);
    target += chi_square_lpdf(one_over_lambda_stap[1] | prior_df[1]);
  }
  else if (prior_dist_for_stap == 7) { // product_normal
    target += normal_lpdf(z_beta | 0, 1);
  }


  // Log-priors for theta-scale
  if(prior_dist_for_theta == 1)
    target += normal_lpdf(theta|prior_mean_for_theta, prior_scale_for_theta);
  else if(prior_dist_for_theta == 8){
    target += lognormal_lpdf(theta|prior_mean_for_theta, prior_scale_for_theta);
  }
  else if(prior_dist_for_theta == 9){
    target += beta_lpdf(theta|prior_mean_for_theta, prior_df_for_theta);
  }
