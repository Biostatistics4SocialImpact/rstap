
  vector[prior_dist == 7 ? sum(num_normals) : K] z_delta;
  vector[prior_dist_for_theta == 7 ? sum(num_normals) : Q] z_beta;
  vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
  vector<lower=0>[Q] mix_stap[prior_dist_for_stap == 5 || prior_dist_for_stap == 6];
  real<lower=0> one_over_lambda[prior_dist == 6];
  real<lower=0> one_over_lambda_stap[prior_dist_for_stap == 6];
  vector<lower=0,upper=max_distance>[Q] theta;

