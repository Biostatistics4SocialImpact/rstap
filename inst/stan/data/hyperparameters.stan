
  // hyperparameter values are set to 0 if there is no prior
  vector<lower=0>[K] prior_scale;
  vector<lower=0>[Q] prior_scale_for_stap;
  vector<lower=0>[(Q_s + Q_st) >0 ? (Q_s + Q_st) : 0] prior_scale_for_theta_s;
  vector<lower=0>[(Q_t + Q_st) >0 ? (Q_t + Q_st) : 0] prior_scale_for_theta_t;
  real<lower=0> prior_scale_for_intercept;
  real<lower=0> prior_scale_for_aux;
  vector[K] prior_mean;
  vector[Q] prior_mean_for_stap;
  vector[(Q_s + Q_st) > 0 ? (Q_s + Q_st) : 0] prior_mean_for_theta_s;
  vector[(Q_t + Q_st) > 0 ? (Q_t + Q_st) : 0] prior_mean_for_theta_t;
  real prior_mean_for_intercept;
  real<lower=0> prior_mean_for_aux;
  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
  real<lower=0> prior_df_for_aux;
  vector<lower=0>[(Q_s + Q_st) >0 ? (Q_s + Q_st) : 0] prior_df_for_theta_s;
  vector<lower=0>[(Q_t + Q_st) >0 ? (Q_t + Q_st) : 0] prior_df_for_theta_t;
  vector<lower=0>[Q] prior_df_for_stap;
  
  int<lower=2> num_normals[prior_dist == 7 ? K : 0];
  int<lower=2> num_normals_for_stap[prior_dist_for_stap == 7 ? Q : 0];
