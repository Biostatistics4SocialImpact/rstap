// GLM for a Gaussian outcome
functions {
#include /functions/common_functions.stan
#include /functions/continuous_likelihoods.stan
}
data{

    #include /data/dnd_data.stan
    int N;
    int M;
    int n;
    int K;
    vector[N] y;
    vector[K] zbar;
    matrix[N,K] Z;
    matrix[1,M] dists_crs;
    int u_s[N,2];
    matrix[n,N] subj_mat;
    matrix[n,1] subj_n;
}
parameters{
    real<lower=0> theta_s;
    real gamma;
    vector[1] beta;
    vector[K] delta;
    real<lower=0> aux;
}
transformed parameters {

    matrix[N,1] X;
    matrix[N,1] X_bar;
    matrix[N,1] X_delta;
    matrix[N,1] X_tilde;
    for(n_ix in 1:N){
        if(u_s[n_ix,1] > u_s[n_ix,2])
            X[n_ix,1] = 0;
        else
            X[n_ix,1] = sum(exp(- dists_crs[1,u_s[n_ix,1]:u_s[n_ix,2]] * inv(theta_s)));
    }

    X_bar = subj_mat' * ((subj_mat * X) .* subj_n);
    X_delta = X - X_bar;
    X_tilde[,1] = (X_delta[,1] - mean(X_delta[,1])) / sd(X_delta[,1]);
}
model {

    vector[N] eta;
    eta = gamma +  Z*delta + X_tilde*beta;
    y ~ normal(eta,aux);
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
    vector[N] eta;
    eta = alpha +  Z*delta + X_tilde*beta;

    for(n_ix in 1:N)
        mean_PPD = mean_PPD + normal_rng(eta[n_ix],aux);
    }
   mean_PPD = mean_PPD / N;

}
