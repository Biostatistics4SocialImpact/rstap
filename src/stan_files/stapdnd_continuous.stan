// GLM for a Gaussian outcome
functions {
#include /functions/common_functions.stan
#include /functions/continuous_likelihoods.stan
}
data{
#include /data/dnd_data.stan
    vector[N] y;


}
parameters{
    real<lower=0> theta_s;
    real gamma;
    vector[Q] beta;
    vector[num_bar] beta_bar;
    vector[K] delta;
    real<lower=0> aux_unscaled;
}
transformed parameters {

    real aux = aux_unscaled;
    matrix[N,1] X;
    matrix[N,1] X_bar;
    matrix[N,1] X_delta;
    matrix[N,1] X_tilde;
    for(q_ix in 1:Q){
        for(n_ix in 1:N){
            if(u_s[n_ix,q_ix] > u_s[n_ix,q_ix])
                X[n_ix,q_ix] = 0;
            else
                X[n_ix,q_ix] = sum(exp(- dists_crs[1,u_s[n_ix,(q_ix*2) -1]:u_s[n_ix,(q_ix * 2 )]] * inv(theta_s)));
        }
    }

    X_bar = subj_mat' * ((subj_mat * X) .* subj_n);
    X_delta = X - X_bar;
    X_tilde[,1] = (X_delta[,1] - mean(X_delta[,1])) / sd(X_delta[,1]);
}
model {

#include /model/make_dnd_eta.stan
    y ~ normal(gamma + eta,aux);
    aux ~ cauchy(0,5);
    beta ~ normal(0,3);
    delta ~ normal(0,3);
    theta_s ~ lognormal(0,1);
    gamma ~ normal(25,5);
}
generated quantities{

    real alpha;
    real mean_PPD;
    vector[1] adj_beta;
    adj_beta[1] = beta[1] / sd(X_delta[,1]);
    alpha = gamma - dot_product(zbar,delta) - mean(X_tilde[,1]) * adj_beta[1];
    mean_PPD = 0;
    {
#include /model/make_dnd_eta.stan
    eta = eta + gamma;
    for(n_ix in 1:N)
        mean_PPD = mean_PPD + normal_rng(eta[n_ix],aux);
    }
   mean_PPD = mean_PPD / N;

}
