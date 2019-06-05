  vector[K] delta;
  vector[Q] beta;
  matrix[NN,Q] X;
  matrix[NN,Q] X_tilde;
  vector[q] b;
  vector[len_theta_L] theta_L;

  //construction of X, X_tilde
  {
      int cnt_s = 1;
      int cnt_t = 1;
      int cnt_shape_s = 1;
      int cnt_shape_t = 1;
      for(q_ix in 1:Q){
          for(n in 1:NN){
            if(stap_code[q_ix] == 0)
                X[n,q_ix] = assign_exposure(log_ar[q_ix], weight_mat[q_ix,1], u_s, dists_crs[cnt_s], theta_s[cnt_s], shape_s, cnt_shape_s, q_ix, n);
            else if(stap_code[q_ix] == 1)
                X[n,q_ix] = assign_exposure(log_ar[q_ix], weight_mat[q_ix,2], u_t, times_crs[cnt_t], theta_t[cnt_t], shape_t, cnt_shape_t, q_ix, n);
            else
                X[n,q_ix] = assign_st_exposure(log_ar[q_ix], weight_mat[q_ix,1], weight_mat[q_ix,2], u_s, u_t, dists_crs[cnt_s], times_crs[cnt_t], theta_s[cnt_s], theta_t[cnt_t], shape_s, shape_t, cnt_shape_s, cnt_shape_t, q_ix, n);
        }
            if(stap_code[q_ix] == 0 || stap_code[q_ix] == 2){
                cnt_s = cnt_s + 1;
                cnt_shape_s += weight_mat[q_ix,1] == 5 ? 1 : 0 ;
            }
            if(stap_code[q_ix] == 1 || stap_code[q_ix] == 2){
                cnt_t = cnt_t + 1;
                cnt_shape_t += weight_mat[q_ix,2] == 6 ? 1 : 0 ; 
            }
      }
  }

  X_tilde = centerscale(X);

  // addition of priors to 
  if(prior_dist == 0) delta = z_delta;
  else if (prior_dist == 1) delta = z_delta .* prior_scale + prior_mean;
  else if (prior_dist == 2) {
        for (k in 1:K) {
            delta[k] = CFt(z_delta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
        }
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
  else if (prior_dist_for_stap == 2) for (q_ix in 1:Q) {
    beta[q_ix] = CFt(z_beta[q_ix], prior_df_for_stap[q_ix]) * prior_scale_for_stap[q_ix] + prior_mean_for_stap[q_ix];
  }
  else if (prior_dist_for_stap == 5) // laplace
    beta = prior_mean_for_stap + prior_scale_for_stap .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist_for_stap == 6) // lasso
    beta = prior_mean_for_stap + one_over_lambda[1] * prior_scale_for_stap .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist_for_stap == 7) { // product_normal
    int z_pos = 1;
    for (q_ix in 1:Q) {
      beta[q_ix] = z_beta[z_pos];
      z_pos = z_pos + 1;
      for (n in 2:num_normals_for_stap[q_ix]) {
        beta[q_ix] = beta[q_ix] * z_delta[z_pos];
        z_pos = z_pos + 1;
      }
      beta[q_ix] = delta[q_ix] * prior_scale_for_stap[q_ix] ^ num_normals_for_stap[q_ix] + prior_mean_for_stap[q_ix];
    }
  }
