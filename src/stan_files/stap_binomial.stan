// STKAP model for a Gaussian outcome with no link function
functions{
  /**
    * Scale/center Matrix M
    * @param M matrix
    * @return A matrix M' s.t. M*1 = 0
    */
    matrix centerscale(matrix M){
        matrix[cols(M), rows(M)] out;
        for(q in 1:cols(M))  
            out[q] = ((M[,q] - mean(M[,q]))/sd(M[,q]))';
        return(out');
   }
   /** 
    * Get column means,sds from Matrix M
    * @param M matrix
    * @return a vector m
    */ 
    vector get_meansds(matrix M){
        vector[2*cols(M)] m;
        for(q in 1:2*cols(M)){
            if(q<=cols(M))
                m[q] = mean(M[,q]);
            else
                m[q] = sd(M[,(q-cols(M))]);
        }
        return(m);
    }
  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_binom(vector eta, int link) {
    if (link == 1)      return(inv_logit(eta)); // logit
    else if (link == 2) return(Phi(eta)); // probit
    else if (link == 3) return(atan(eta) / pi() + 0.5);  // cauchit
    else if (link == 4) return(exp(eta)); // log
    else if (link == 5) return(inv_cloglog(eta)); // cloglog
    else reject("Invalid link");
    return eta; // never reached
  }
  
  /**
  * Increment with the unweighted log-likelihood
  * @param y An integer array indicating the number of successes
  * @param trials An integer array indicating the number of trials
  * @param eta A vector of linear predictors
  * @param link An integer indicating the link function
  * @return lp__
  */
  real ll_binom_lp(int[] y, int[] trials, vector eta, int link) {
    if (link == 1) target += binomial_logit_lpmf(y | trials, eta);
    else if (link <  4) target += binomial_lpmf( y | trials, linkinv_binom(eta, link));
    else if (link == 4) {  // log
      for (n in 1:num_elements(y)) {
        target += y[n] * eta[n];
        target += (trials[n] - y[n]) * log1m_exp(eta[n]);
        target += lchoose(trials[n], y[n]);
      }
    }
    else if (link == 5) {  // cloglog
      for (n in 1:num_elements(y)) {
        real neg_exp_eta = -exp(eta[n]);
        target += y[n] * log1m_exp(neg_exp_eta);
        target += (trials[n] - y[n]) * neg_exp_eta;
        target += lchoose(trials[n], y[n]);
      }
    }
    else reject("Invalid link");
    return target();
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_binom(int[] y, int[] trials, vector eta, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link == 1) {  // logit
      for (n in 1:N) 
        ll[n] = binomial_logit_lpmf(y[n] | trials[n], eta[n]);
    }
    else if (link <= 5) {  // link = probit, cauchit, log, or cloglog
      vector[N] pi = linkinv_binom(eta, link); // may be unstable
      for (n in 1:N) ll[n] = binomial_lpmf(y[n] | trials[n], pi[n]) ;
    }
    else reject("Invalid link");
    return ll;
  }
}
data {
    int<lower=0> N; // number of subjects
    int<lower=0> q; // number of bef covariates
    int<lower=0> p; // number of subj specific covariates
    int<lower=0> link; //
    int<lower=0> trials[N];
    int<lower=0> y[N];
    matrix[N,p] Z;
    int<lower=0,upper=1> log_transform[q]; // log transform indicator for BEF covariates
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
            else if(log_transform[q])
                X[n_ix,q_ix] = log(1+sum( erfc(spatial_mat[q_ix][u[n_ix,q_ix,1]:u[n_ix,q_ix,2]] * inv(d_constraint * thetas[q_ix]))));
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
