#' Fitting Generalized Linear STKAP models
#'
#'@param y vector of continuous outcomes
#'@param X vector of environmental features integer coding
#'@param Z n*p design matrix of subject specific covariates
#'@export stap_glm
stap_glm <- function(y, Z, dists_csr, u, max_distance = 3, family = gaussian(),
                         prior_intercept = normal(location = 25, scale = 5),
                         prior_beta_one = normal(location = 0, scale = 5),
                         prior_beta_two = normal(location = 0, scale = 5),
                         prior_theta = list(theta_one = normal(location = 1.5, scale = .5)),
                         prior_sigma = cauchy(location = 0, scale = 5),
                         chains = 1,
                         iter = 500,
                         warmup = floor(iter/2),
                         thin = 1,
                         cores = getOption("mc.cores", 1L),
                         seed = Sys.Date(),
                         control = list(...),
                         ...){
 
    family <- validate_family(family)
    supported_families <- c("binomial","gaussian","poisson")
    fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
    if(!length(fam)){
        stop("'family' must be one of ", paste(supported_families, collapse = ', '))

    supported_links <- supported_glm_links(supported_families[fam])
    link <- which(supported_links == family$link)
    if(!length(link))
        stop("'link' must be one of", paste( supported_links, collapse = ', '))

    beta_naught_p <- assign_dist(prior_intercept)
    beta_one_p <- assign_dist(prior_beta_one)
    beta_two_p <- assign_dist(prior_beta_two)
    theta_p  <- lapply(prior_theta, assign_dist)
    sigma_p <- assign_dist(prior_sigma)

    y <- validate_glm_outcome_support(y,family)
    if(family$family == 'binomial')

    standata <- list(N = length(y),
                     q = nrow(dists_csr),
                     d_constraint = max_distance,
                     p = ncol(Z),
                     M = ncol(dists_csr),
                     spatial_mat = dists_csr,
                     u = u,
                     y = y,
                     Z = Z,
                     family = fam,
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
                     prior_sigma = sigma_p$prior_dist,
                     prior_sigma_df = sigma_p$prior_df,
                     prior_sigma_loc = sigma_p$prior_mean,
                     prior_sigma_scale = sigma_p$prior_scale,
                     prior_theta = array(sapply(theta_p, "[[","prior_dist",simplify='array')),
                     prior_theta_df = array(sapply(theta_p, "[[","prior_df",simplify='array')),
                     prior_theta_loc = array(sapply(theta_p, "[[","prior_mean", simplify='array')),
                     prior_theta_scale = array(sapply(theta_p, "[[","prior_scale",simplify='array')))


    stanfit <- stanmodels$stap_glm


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
supported_glm_links(family_name){
    switch(family_name,
           binomial = c("logit","log"),
           gaussian = c("identity","log"),
           poisson = c("log"),
           stop("unsupported family")
           )
}

# Verify that outcome values match support implied by family object
#
# @param y outcome variable
# @param family family object
# @return y (possibly slightly modified) unless an error is thrown
#
validate_glm_outcome_support <- function(y, family) {
  .is_count <- function(x) {
    all(x >= 0) && all(abs(x - round(x)) < .Machine$double.eps^0.5)
  }
  
  fam <- family$family
  
  if (!(fam=='binomial')) {
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
    if (is.gaussian(fam)) {
      return(y)
    } else if (is.poisson(fam) && !.is_count(y)) {
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

