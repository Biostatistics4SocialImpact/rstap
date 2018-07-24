
  // correction to eta if model has no intercept (because X is centered)

  eta = eta + dot_product(zbar,delta) + dot_product(colmeans(X), beta); 
