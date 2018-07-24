
  vector[K] delta;
  vector[Q] beta;
  matrix[N,Q] X;
  matrix[N,Q] X_tilde;

  //construction of X, X_tilde
  for(n in 1:N){
    for(q in 1:Q){
        if(u_array[n,q,1]>u_array[n,q,2])
            X[n,q] = 0;
        else
            X[n,q] = sum( erfc(dists_crs[q][u_array[n,q,1]:u_array[n,q,2]] * inv(theta[q])) );
    }
  }

  X_tilde = centerscale(X);
  if(prior_dist == 0) delta = z_delta;
  else if (prior_dist == 1) delta = z_delta .* prior_scale + prior_mean;
  else if (prior_dist == 2) for (k in 1:K) {
    delta[k] = CFt(z_delta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
  }
  else if (prior_dist == 5) // laplace
    delta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_delta;
  else if (prior_dist == 6) // lasso
    delta = prior_mean + one_over_lambda[1] * prior_scale .* sqrt(2 * mix[1]) .* z_delta;
  else if (prior_dist == 7) { // product_normal
    int z_pos = 1;
    for (k in 1:K) {
      delta[k] = z_delta[z_pos];
      z_pos = z_pos + 1;
      for (n in 2:num_normals[k]) {
        delta[k] = delta[k] * z_delta[z_pos];
        z_pos = z_pos + 1;
      }
      delta[k] = delta[k] * prior_scale[k] ^ num_normals[k] + prior_mean[k];
    }
  }


  if(prior_dist_for_stap == 0) beta = z_beta;
  else if (prior_dist_for_stap == 1) beta = z_beta .* prior_scale_for_stap + prior_mean_for_stap;
  else if (prior_dist_for_stap == 2) for (q in 1:Q) {
    beta[q] = CFt(z_beta[q], prior_df_for_stap[q]) * prior_scale_for_stap[q] + prior_mean_for_stap[q];
  }
  else if (prior_dist_for_stap == 5) // laplace
    beta = prior_mean_for_stap + prior_scale_for_stap .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist_for_stap == 6) // lasso
    beta = prior_mean_for_stap + one_over_lambda[1] * prior_scale_for_stap .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist_for_stap == 7) { // product_normal
    int z_pos = 1;
    for (q in 1:Q) {
      beta[q] = z_beta[z_pos];
      z_pos = z_pos + 1;
      for (n in 2:num_normals_for_stap[q]) {
        beta[q] = beta[q] * z_delta[z_pos];
        z_pos = z_pos + 1;
      }
      beta[q] = delta[q] * prior_scale_for_stap[q] ^ num_normals_for_stap[q] + prior_mean_for_stap[q];
    }
  }
