
  // hyperparameter values are set to 0 if there is no prior
  vector<lower=0>[K] prior_scale;
  vector<lower=0>[Q] prior_scale_for_stap;
  vector<lower=0>[Q] prior_scale_for_theta;
  real<lower=0> prior_scale_for_intercept;
  real<lower=0> prior_scale_for_aux;
  vector[K] prior_mean;
  vector[Q] prior_mean_for_stap;
  vector[Q] prior_mean_for_theta;
  real prior_mean_for_intercept;
  real<lower=0> prior_mean_for_aux;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
  real<lower=0> prior_df_for_aux;
  vector<lower=0>[Q] prior_df_for_theta;
  vector<lower=0>[Q] prior_df_for_stap;

  int<lower=2> num_normals[prior_dist == 7 ? K : 0];
  int<lower=2> num_normals_for_stap[prior_dist_for_stap == 7 ? Q:0];
