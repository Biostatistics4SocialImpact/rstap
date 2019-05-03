  vector[N] eta;  // linear predictor
  eta = Z * delta + X_tilde * beta;
  if(num_bar>0)
    eta = eta + X_bar[,bar_arr] * beta_bar;
  //if (has_offset == 1) eta = eta + offset;
