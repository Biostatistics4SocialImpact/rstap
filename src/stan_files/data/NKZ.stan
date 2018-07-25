
  // dimensions
  int<lower=0> N;  // number of observations
  int<lower=0> K;  // number of fixed predictors
  int<lower=0> Q; // number of STAP's
  int<lower=0> M; // Max number of BEF's within inclusion distance
  
  // data
  vector[K] zbar;               // predictor means
  matrix[N,K] Z;       // centered predictor matrix 
  vector<lower=0>[M] dists_crs[Q]; //distance matrix
  int u_array[N,Q*2]; //index holder array
  real<lower=0> max_distance; // Max distance
