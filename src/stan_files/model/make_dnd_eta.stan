  vector[N] eta;  // linear predictor
  eta = Z * delta + X_tilde * beta;
  if (has_offset == 1) eta = eta + offset;
