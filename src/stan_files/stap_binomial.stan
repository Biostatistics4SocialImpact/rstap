// STAP  model for a binomial model with one of five link functions 
functions{
    #include /functions/common_functions.stan
    #include /functions/binomial_likelihoods.stan
}
data {
    int<lower=0> N; // number of subjects
    int<lower=0> q; // number of ef_types
    int<lower=0> p; // number of subj specific covariates
    int<lower=0> link; //
    int<lower=0> trials[N];
    int<lower=0> y[N];
    matrix[N,p] Z;
    #include "prior_data.stan" // prior_data
    int<lower=1> M; // Maximum number of EFs within boundary distance
    int u[N,q,2]; // index holder array  
    vector<lower=0>[M] spatial_mat[q];
    real d_constraint; // maximum distance
}
parameters {
    real beta_naught;
    vector[p] beta_one;
    vector[q] beta_two;
   real<lower=0, upper = 1> thetas[q];
}
transformed parameters{
    matrix [N,q] X;
    matrix [N,q] X_tilde;
    for(n_ix in 1:N){
        for(q_ix in 1:q){
            if(u[n_ix,q_ix,1]>u[n_ix,q_ix,2])
                X[n_ix,q_ix] = 0;
            else
                X[n_ix,q_ix] = sum( erfc(spatial_mat[q_ix][u[n_ix,q_ix,1]:u[n_ix,q_ix,2]] * inv(d_constraint * thetas[q_ix])));
        }
    } 
    X_tilde = centerscale(X);

}
model{ 
    #include "priors_stap_lm.stan"
    target += pw_binom(y, trials, beta_naught + Z * beta_one + X_tilde * beta_two, link);
}
