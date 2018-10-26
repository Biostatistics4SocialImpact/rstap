
  vector[N[1]] eta0;
  vector[N[2]] eta1;
  eta0 = rep_vector(0.0, N[1]);
  eta1 = rep_vector(0.0, N[2]);

  if (has_intercept == 0) {
    real tmp;
    tmp = dot_product(zbar, delta);
    eta0 = eta0 + tmp;
    eta1 = eta1 + tmp;
  }else{
    eta0 =  eta0 + Z0 * delta + X_tilde[y_0,] * beta;
    eta1 = eta1 + Z1 * delta + X_tilde[y_1,] * beta;
  }
  if (has_offset == 1) {
    eta0 = eta0 + offset0;
    eta1 = eta1 + offset1;
  }
  if (special_case) for (i in 1:t) {
    eta0 = eta0 + b[V0[i]];
    eta1 = eta1 + b[V1[i]];
  }
  else if (t > 0) {
    eta0 = eta0 + csr_matrix_times_vector(N[1], q, w0, v0, u0, b);
    eta1 = eta1 + csr_matrix_times_vector(N[2], q, w1, v1, u1, b);
  }
