#' Fitting Linear STKAP models
#'
#'@param y vector of continuous outcomes
#'@param Z n*p design matrix of subject specific covariates
#'@param dists_csr distance crs array
#'@param u crs index array
#'@param max_distance inclusion distance
#'@param prior_intercept
#'@param prior_beta_one
#'@param prior_BEFS
#'@param prior_theta
#'@param prior_sigma
#'@export stap_lm
stap_lm <- function(y, Z, dists_csr, u, max_distance = 3,
                    prior = normal(),
                    prior_intercept = normal(location = 25, scale = 5),
                    prior_stap = normal(),
                    prior_theta = list(theta_one = normal(location = 1.5, scale = .5)), 
                    prior_sigma = cauchy(location = 0, scale = 5),
                    ...){

    
    if(max_distance<max(dists_csr))
        stop("max_distance must be the maximum possible distance amongst all distances in dists_csr")

    beta_naught_p <- assign_dist(prior_intercept)
    beta_one_p <- assign_dist(prior_beta_one)
    beta_two_p <- assign_dist(prior_beta_two)
    theta_p  <- lapply(prior_theta,assign_dist)
    sigma_p <- assign_dist(prior_sigma)


    standata <- list(N = length(y),
                     q = nrow(dists_csr),
                     d_constraint = max_distance,
                     p = ncol(Z),
                     M = ncol(dists_csr),
                     spatial_mat = dists_csr,
                     u = u,
                     y = y,
                     Z = Z,
                     prior_intercept = beta_naught_p$prior_dist,
                     prior_intercept_loc = beta_naught_p$prior_mean,
                     prior_intercept_scale = beta_naught_p$prior_scale,
                     prior_intercept_df = beta_naught_p$prior_df,
                     prior_beta_one = beta_one_p$prior_dist,
                     prior_beta_one_df = beta_one_p$prior_df,
                     prior_beta_one_loc = beta_one_p$prior_mean,
                     prior_beta_one_scale = beta_one_p$prior_scale,
                     prior_beta_two = beta_two_p$prior_dist,
                     prior_beta_two_df = beta_two_p$prior_df,
                     prior_beta_two_loc = beta_two_p$prior_mean,
                     prior_beta_two_scale = beta_two_p$prior_scale,
                     prior_sigma = sigma_p$prior_dist,
                     prior_sigma_df = sigma_p$prior_df,
                     prior_sigma_loc = sigma_p$prior_mean,
                     prior_sigma_scale = sigma_p$prior_scale,
                     prior_theta = array(sapply(theta_p, "[[","prior_dist",simplify='array')),
                     prior_theta_df = array(sapply(theta_p, "[[","prior_df",simplify='array')),
                     prior_theta_loc = array(sapply(theta_p, "[[","prior_mean", simplify='array')),
                     prior_theta_scale = array(sapply(theta_p, "[[","prior_scale",simplify='array')))


    stanfit <- stanmodels$stap_lm


    sampling_args <- set_sampling_args(object = stanfit, data = standata,
                                       chains = chains,
                                       iter = iter,
                                       cores = cores,
                                       user_dots = control,
                                       warmup = warmup)

    out <- do.call(rstan::sampling,sampling_args)
}
