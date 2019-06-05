  vector[prior_dist == 7 ? sum(num_normals) : K] z_delta;
  vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
  vector<lower=0>[prior_dist_for_stap == 7 ? sum(num_normals) : Q] z_beta;
  vector<lower=0>[Q] mix_stap[prior_dist_for_stap == 5 || prior_dist_for_stap == 6];
  real<lower=0> one_over_lambda[prior_dist == 6];
  real<lower=0> one_over_lambda_stap[prior_dist_for_stap == 6];
  real<lower=0,upper=max_distance> theta_s[Q_s+Q_st]; //scale for spatial weight function
  real<lower=0> shape_s[num_s_wei]; // shape for spatial weight function
  real<lower=0> shape_t[num_t_wei]; // shape for temporal weight function
  real<lower=0,upper=max_time> theta_t[Q_t+Q_st];//scale for temporal weight function
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0, upper =1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;
