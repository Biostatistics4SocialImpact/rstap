// STAP model for a Gaussian outcome with no link function
data {
    int<lower=0> N; // number of subjects
    int<lower=0> q; // number of ef_types
    int<lower=0> p; // number of subj specific covariates
    #include "prior_data.stan" // prior_data
    int<lower=1> M; // Maximum number of EFs within boundary distance
    int u[q,N,2]; // index holder array  
    vector<lower=0>[M] spatial_mat[q];
    real d_constraint; // maximum distance
    matrix[N,p] Z;
    vector[N] y;
}
parameters {
    real beta_naught;
    vector[p] beta_one;
    vector[q] beta_two;
    real<lower=0> sigma;
    real<lower=0, upper = d_constraint> thetas[q];
}
model { 
    #include "priors_stap_lm.stan"
    {
        matrix[N,q] X_tilde;
        for(n_ix in 1:N){
            for(q_ix in 1:q)
                X_tilde[n_ix,q_ix] = sum( erfc(spatial_mat[q_ix][u[q_ix,n_ix,1]:u[q_ix,n_ix,2]] * inv(thetas[q_ix]) ) ); 
        }
        y ~ normal( beta_naught + Z * beta_one +  X_tilde * beta_two,sigma);
    }
}
