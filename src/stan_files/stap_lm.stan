// STKAP model for a Gaussian outcome with no link function
data {
    int<lower=0> N; // number of subjects
    int<lower=0> q; // number of ef_types
    int<lower=0> p; // number of subj specific covariates
    int<lower=0> qt; // number of total efs
    #include "prior_data.stan" // prior_data
    int<lower=1> M; // Maximum number of EFs within boundary distance
    int u[q,N,2]; // placement array  
    vector<lower=0>[M] spatial_mat[q];
    matrix[N,p] Z;
    vector[N] y;
}
parameters {
    vector[p+1] beta_one;
    vector[q] beta_two;
    real<lower=0> sigma;
    real<lower=0> thetas[q];
}
model { 
    #include "priors_stkap_lm.stan"
    {
        matrix[N,q] X_tilde;
        for(n_ix in 1:N){
            for(q_ix in 1:q)
                X_tilde[n_ix,q_ix] = sum(erfc(spatial_mat[q_ix][u[q_ix,n_ix,1]:u[q_ix,n_ix,2]] * inv(thetas[q_ix]) ) );
        }
        y ~ normal( Z * beta_one +  X_tilde * beta_two,sigma);
    }
}
