
  // dimensions
  int<lower=0> N;  // number of observations
  int<lower=0> K;  // number of fixed predictors
  int<lower=0> Q; // number of STAP,SAP, and TAPs
  int<lower=0,upper=1> log_ar[Q]; // 0 = no log; 1 = log
  int<lower=0,upper=2> stap_code[Q]; // 0 = sap ; 1 = tap ; 2 = stap
  int<lower=0> Q_t; // number of taps
  int<lower=0> Q_s; // number of saps
  int<lower=0,upper=Q-Q_t-Q_s> Q_st; // number of staps
  int<lower=0> M; // Max number of BEF's within inclusion distance
  
  // data
  vector[K] zbar;               // predictor means
  matrix[N,K] Z;       // centered predictor matrix 
  vector<lower=0>[M] dists_crs[Q_s+Q_st]; //distance crs matrix 
  vector<lower=0>[M] times_crs[Q_t + Q_st]; // time crs matrix
  int<lower=0> w[Q,2]; // spatial-temporal weight functions 0 = no component; 1 = erf; 2 = cerf; 3 = exp; 4 = cexp;
  // Meta data
  int u_s[N,(Q_s+Q_st)*2]; //index holder array for distances
  int u_t[N,(Q_s+Q_st)*2]; // index holder array for times
  real<lower=0> max_distance; // Max distance for spatial scales
  real<lower=0> max_time; // Max Time
