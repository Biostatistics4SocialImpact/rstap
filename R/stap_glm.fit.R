#' Fitting Generalized Linear STKAP models
#'
#'@param y n length vector or n x 2 matrix of outcomes
#'@param dists_crs q x M matrix of distances between outcome observations and environmental features where q is the number of spatial covariates, and M is the maximum number of environmental features amongst all q features
#'@param max_distance the upper bound of distance for which all
#'@param Z n x p design matrix of subject specific covariates
#'@export stap_glm
stap_glm.fit <- function(y, Z, dists_crs, u, max_distance = 3L,
                         weights = rep(1,NROW(y)),
                         offset = rep(0, NROW(y)),
                         family = stats::gaussian(),
                         ...,
                         prior_intercept = normal(),
                         prior_beta_one = normal(),
                         prior_beta_two = normal(),
                         prior_theta = list(theta_one = normal()),
                         prior_aux = cauchy(location = 0L, scale = 5L),
                         adapt_delta = NULL){

    family <- validate_family(family)
    supported_families <- c("binomial","gaussian","poisson")
    fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
    if(!length(fam))
        stop("'family' must be one of ", paste(supported_families, collapse = ', '))
    if(max_distance < max(dists_crs))
        stop("max_distance must be the maximum possible distance amongst all distances in dists_crs")
    if(max(u) != ncol(dists_crs)){
      ix <- which(u==max(u), arr.ind = T)
      if(length(dim(u))==2)
        if(any(u[ix[,1],1]<=u[ix[,1],2]))
          stop("Maximum index value of u must be equivalent to the number of columns in the distance matrix")
      if(any(u[ix[,1],ix[,2],1]<= u[ix[,1],ix[,2],2]))
        stop("Maximum index value of u must be equivalent to the number of columns in the distance matrix")
    }
    supported_links <- supported_glm_links(supported_families[fam])
    link <- which(supported_links == family$link)
    if(!length(link))
        stop("'link' must be one of", paste( supported_links, collapse = ', '))
    if(length(prior_theta)!=nrow(dists_crs))
        stop("Insufficient number of priors set for spatial scale parameter")
    
    if (binom_y_prop(y, family, weights))
        stop("To specify 'y' as proportion of successes and 'weights' as ",
             "number of trials please use stan_glm rather than calling ",
             "stan_glm.fit directly.")

    y <- validate_glm_outcome_support(y,family)
    if(is.binomial(family$family) && NCOL(y) == 2L){
        trials <- as.integer(y[,1L] + y[,2L])
        y <- as.integer(y[,1L])
    }
    
    # useless assignments to pass R CMD check
    has_intercept <- 
        prior_df <- prior_df_for_intercept <- prior_df_for_aux <- prior_df_for_smooth <-
        prior_dist <- prior_dist_for_intercept <- prior_dist_for_aux <- prior_mean <- 
        prior_mean_for_intercept <- prior_mean_for_aux <- prior_scale <-
        prior_scale_for_intercept <- prior_scale_for_aux <- prior_autoscale <- 
        prior_autoscale_for_intercept <- prior_autoscale_for_aux <- 
        prior_autoscale_for_smooth <- global_prior_scale <- global_prior_df <-
    
    
    beta_naught_p <- assign_dist(prior_intercept)
    beta_one_p <- assign_dist(prior_beta_one)
    beta_two_p <- assign_dist(prior_beta_two)
    theta_p  <- lapply(prior_theta, assign_dist)
    

    standata <- list(q = nrow(dists_crs),
                     d_constraint = max_distance,
                     p = ncol(Z),
                     M = ncol(dists_crs),
                     spatial_mat = dists_crs,
                     u = u,
                     Z = Z,
                     link = link,
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
                     prior_theta = array(sapply(theta_p, "[[","prior_dist",simplify='array')),
                     prior_theta_df = array(sapply(theta_p, "[[","prior_df",simplify='array')),
                     prior_theta_loc = array(sapply(theta_p, "[[","prior_mean", simplify='array')),
                     prior_theta_scale = array(sapply(theta_p, "[[","prior_scale",simplify='array')))

    if(family$family == 'binomial' && ncol(y) == 2L){
        trials <- y[, 1L]
        y <- y[,2]
        standata$y <- y
        standata$trials <- trials
        standata$N <- length(y)
        stanfit <- stanmodels$stap_binomial
    } else if (ncol(y) == 1L || is.null(ncol(y))){
        ## still to be implemented
    }

    sampling_args <- set_sampling_args(object = stanfit, data = standata,
                                       chains = chains,
                                       iter = iter,
                                       cores = cores,
                                       user_dots = control,
                                       warmup = warmup)

    out <- do.call(rstan::sampling,sampling_args)
}

# internal ---------------------------------------------------------------------------------------------------------

# @param family_name: string naming the family
# @return character vector of supported link functions for the family
supported_glm_links <- function(family_name){
    switch(
      family_name,
      binomial = c("logit","probit","cauchit", "log","cloglog"),
      gaussian = c("identity", "log", "inverse"),
      poisson = c("log", "identity", "sqrt"),
      stop("unsupported family")
    )
}

# Verify that outcome values match support implied by family object
#
# @param y outcome variable
# @param family family object
# @return y (possibly slightly modified) unless an error is thrown
#
validate_glm_outcome_support <- function(y, family){
  .is_count <- function(x) {
    all(x >= 0) && all(abs(x - round(x)) < .Machine$double.eps^0.5)
  }

  fam <- family$family

  if (fam!='binomial') {
    # make sure y has ok dimensions (matrix only allowed for binomial models)
    if (length(dim(y)) > 1) {
      if (NCOL(y) == 1) {
        y <- y[, 1]
      } else {
        stop("Except for binomial models the outcome variable ",
             "should not have multiple columns.",
             call. = FALSE)
      }
    }

    # check that values match support for non-binomial models
    if (fam!='gaussian') {
      return(y)
    } else if (fam=='poisson' && !.is_count(y)) {
      stop("All outcome values must be counts for Poisson models",
           call. = FALSE)
    }
  } else { # binomial models
    if (NCOL(y) == 1L) {
      if (is.numeric(y) || is.logical(y))
        y <-  as.integer(y)
      if (is.factor(y))
        y <- fac2bin(y)
      if (!all(y %in% c(0L, 1L)))
        stop("All outcome values must be 0 or 1 for Bernoulli models.",
             call. = FALSE)
    } else if (isTRUE(NCOL(y) == 2L)) {
      if (!.is_count(y))
        stop("All outcome values must be counts for binomial models.",
             call. = FALSE)
    } else {
      stop("For binomial models the outcome should be a vector or ",
           "a matrix with 2 columns.",
           call. = FALSE)
    }
  }
  return(y)
}

