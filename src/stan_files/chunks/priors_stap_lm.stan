// increment log_probability with appropriate priors


// intercept priors
    if(prior_intercept == 1) target += normal_lpdf(beta_naught|prior_intercept_loc,prior_intercept_scale);
    else if(prior_intercept == 2) target += student_t_lpdf(beta_naught|prior_intercept_df, prior_intercept_loc,prior_intercept_scale);
    else if(prior_intercept == 3) target += cauchy_lpdf(beta_naught|prior_intercept_loc,prior_intercept_scale);
    else if(prior_intercept == 4) target += lognormal_lpdf(beta_naught|prior_intercept_loc,prior_intercept_scale);
    else if(prior_intercept >= 5 && prior_intercept <=8){ 
        print("this prior not yet implimented");
    }
    else if (prior_intercept == 8) target += exponential_lpdf(beta_naught | prior_intercept_loc);

// beta_one priors
    if(prior_beta_one == 1) target += normal_lpdf(beta_one|prior_beta_one_loc,prior_beta_one_scale);
    else if(prior_beta_one == 2) target += student_t_lpdf(beta_one| prior_beta_one_df,prior_beta_one_loc,prior_beta_one_scale);
    else if(prior_beta_one == 3) target += cauchy_lpdf(beta_one | prior_beta_one_loc, prior_beta_one_scale);
    else if(prior_intercept == 4) target += lognormal_lpdf(beta_one|prior_beta_one_loc,prior_beta_one_scale);
    else if(prior_beta_one >= 5 && prior_beta_one <=8){ 
        print("This prior is not yet implimented");
    }
    else if(prior_beta_one == 8) target += exponential_lpdf(beta_one | prior_beta_one_loc);

// beta_two priors 

    if(prior_beta_two == 1) target += normal_lpdf(beta_two| prior_beta_two_loc, prior_beta_two_scale);
    else if(prior_beta_two == 2) target += student_t_lpdf(beta_two| prior_beta_two_df,prior_beta_two_loc,prior_beta_two_scale);
    else if(prior_beta_two == 3) target += cauchy_lpdf(beta_two | prior_beta_two_loc, prior_beta_two_scale);
    else if(prior_beta_two == 4) target += lognormal_lpdf(beta_two | prior_beta_two_loc, prior_beta_two_scale);
    else if(prior_beta_two >= 5 && prior_beta_two <=8){ 
        print("This prior is not yet implimented");
    }
    else if(prior_beta_two == 8) target += exponential_lpdf(beta_two | prior_beta_two_loc);

// sigma_priors

    if(prior_sigma == 1) target += normal_lpdf(sigma | prior_sigma_loc, prior_sigma_scale);
    else if(prior_sigma == 2) target += student_t_lpdf(sigma| prior_sigma_df,prior_sigma_loc,prior_sigma_scale);
    else if(prior_sigma == 3) target += cauchy_lpdf(sigma | prior_sigma_loc, prior_sigma_scale);
    else if(prior_sigma == 4) target += lognormal_lpdf(sigma | prior_sigma_loc, prior_sigma_scale);
    else if(prior_sigma >= 5 && prior_sigma <=8){ 
        print("This prior is not yet implimented");
    }
    else if(prior_sigma == 8) target += exponential_lpdf(sigma | prior_sigma_loc);

// theta_priors

    for(q_ix in 1:q){
        if(prior_theta[q_ix] == 1) target += normal_lpdf(thetas[q_ix] | prior_theta_loc[q_ix], prior_theta_scale[q_ix]);
        else if(prior_theta[q_ix]  == 2) target += student_t_lpdf(thetas[q_ix]| prior_theta_df[q_ix],prior_theta_loc[q_ix], prior_theta_scale[q_ix]);
        else if(prior_theta[q_ix] == 3) target += cauchy_lpdf(thetas[q_ix] | prior_theta_loc[q_ix], prior_theta_scale[q_ix]);
        else if(prior_theta[q_ix] == 4) target += lognormal_lpdf(thetas[q_ix] | prior_theta_loc[q_ix], prior_theta_scale[q_ix]);
        else if(prior_theta[q_ix] == 5) target += beta_lpdf(thetas[q_ix] | prior_theta_loc[q_ix], prior_theta_scale[q_ix]);
        else if(prior_theta[q_ix] >= 6 && prior_theta <=8){ 
            print("This prior is not yet implimented");
        }
        else if(prior_theta[q_ix] == 8) target += exponential_lpdf(thetas[q_ix] | prior_theta_loc[q_ix]);
    }
