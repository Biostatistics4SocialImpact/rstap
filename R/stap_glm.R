#' Fitting Generalized Linear STAP models
#'
#'@param y n length vector or n x 2 matrix of outcomes
#'@param dists_crs q x M matrix of distances between outcome observations and
#' environmental features where q is the number of spatial covariates, 
#' and M is the maximum number of environmental features amongst all q features
#' @param u crs array
#' @param family Same as \code{\link[stats]{glm}} for gaussian, binomial, and poisson
#' @param subject_data
#' @param distance_data
#' @param id_key name of column to join on between subject_data and distance_data
#' @param max_distance the inclusion distance; upper bound for all elements of dists_crs
#' @param weights 
#' @details The \code{stap_glm} function is similar in syntax to 
#' \code{\link[rstanarm]{stan_glm}} except instead of performing full bayesian
#' inference for a generalized linear model stap_glm incorporates spatial 
#' as detailed in in --need to add citation --
#'@export stap_glm
<<<<<<< HEAD
stap_glm <- function(formula,
                     family = stats::gaussian(),
                     subject_data,
                     distance_data,
                     id_key = NULL,
                     max_distance,
                     weights,
                     subset,
                     na.action = NULL,
                     offset = NULL,
                     model = TRUE,
                     z = FALSE,
                     y = TRUE,
                     x = FALSE,
                     contrasts = NULL,
                     ...,
                     prior = normal(),
                     prior_intercept = normal(),
                     prior_stap = normal(),
                     prior_theta = normal(location = max_distance/2, scale = 10),
                     prior_aux = cauchy(location = 0L, scale = 5L),
                     adapt_delta = NULL){
    crs_data <- extract_stap_components(formula,distance_data,
                                        subject_data, id_key, 
                                        max_distance)
    formula <- get_stapless_formula(formula)
    family <- validate_family(family)
    validate_glm_formula(formula)
    subject_data <- validate_data(subject_data, if_missing = environment(formula))
    call <- match.call(expand.dots = TRUE)
    mf <-  match.call(expand.dots = FALSE)
    mf$formula <- formula
    m <- match(c("formula","subset", "weights", "na.action", "offset"),
               table = names(mf), nomatch=0L)
    mf <- mf[c(1L,m)]
    mf$data <- subject_data
    mf$drop.unused.levels <- TRUE
    
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame()) 
    mf <- check_constant_vars(mf)
    mt <- attr(mf, "terms")
    Y <- array1D_check(model.response(mf, type = "any"))
    if(is.empty.model(mt))
        stop("No intercept or predictors specified.", call. = FALSE)
    Z <- model.matrix(mt, mf, contrasts)
    weights <- validate_weights(as.vector(model.weights(mf)))
    offset <- validate_offset(as.vector(model.offset(mf)), y = Y)
    if(binom_y_prop(Y,family, weights)) {
        y1 <- as.integer(as.vector(Y) * weights)
        Y <- cbind(y1, y0 = weights - y1)
        weights <- double(0)
=======
stap_glm <- function(y, Z, dists_crs, u, max_distance = 3L, family = gaussian(),
                         prior_intercept = normal(location = 25L, scale = 5L),
                         prior_beta_one = normal(location = 0L, scale = 5L),
                         prior_beta_two = normal(location = 0L, scale = 5L),
                         prior_theta = list(theta_one = normal(location = 1.5, scale = .5)),
                         prior_sigma = cauchy(location = 0L, scale = 5L),
                         log_transform = rep(1,nrow(dists_crs)),
                         chains = 1L,
                         iter = 500L,
                         warmup = floor(iter/2L),
                         thin = 1L,
                         cores = getOption("mc.cores", 1L),
                         seed = Sys.Date(),
                         control = list(...),
                         ...){

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

    beta_naught_p <- assign_dist(prior_intercept)
    beta_one_p <- assign_dist(prior_beta_one)
    beta_two_p <- assign_dist(prior_beta_two)
    theta_p  <- lapply(prior_theta, assign_dist)
    if(fam == 2)
        sigma_p <- assign_dist(prior_sigma)
    y <- validate_glm_outcome_support(y,family)

    standata <- list(q = nrow(dists_crs),
                     d_constraint = max_distance,
                     p = ncol(Z),
                     M = ncol(dists_crs),
                     spatial_mat = dists_crs,
                     u = u,
                     Z = Z,
                     link = link,
                     log_transform = log_transform,
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
>>>>>>> 4c0736a2f508d6bbf01f7230411ee5a007076eb4
    }
    stapfit <- stap_glm.fit(z = Z, y = Y, weights = weights,
                            dists_crs = crs_data$d_mat, u = crs_data$u,
                            max_distance = max_distance,
                            offset = offset, family = family,
                            prior = prior,
                            prior_intercept = prior_intercept,
                            prior_stap = prior_stap,
                            prior_aux = prior_aux,
                            prior_theta = prior_theta,
                            adapt_delta = adapt_delta,
                            ...)

    sel <- apply(Z, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
    Z <- Z[ , !sel, drop = FALSE]
    fit <- nlist(stapfit, family, formula, subject_data,
                 distance_data,
                 dists_crs = crs_data$d_mat,
                 u = crs_data$u,
                 offset, weights, z = Z, y = Y,
                 model = mf,  terms = mt, call,
                 na.action = attr(mf, "na.action"),
                 contrasts = attr(Z, "contrasts"),
                 stan_function = "stap_glm")
    out <- stapreg(fit)
    out$xlevels <- .getXlevels(mt, mf)
    if (!x)
        out$x <- NULL
    if(!z)
        out$z <- NULL
    if (!y)
        out$y <- NULL
    if (!model)
        out$model <- NULL
    return(out)
}

