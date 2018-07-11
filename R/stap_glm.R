#' Fitting Generalized Linear STKAP models
#'
#'@param y n length vector or n x 2 matrix of outcomes
#'@param dists_crs q x M matrix of distances between outcome observations and
#' environmental features where q is the number of spatial covariates, 
#' and M is the maximum number of environmental features amongst all q features
#' @param u crs array
#' @D_M the inclusion distance; upper bound for all elements of dists_crs
#' @family Same as \code{\link[stats]{glm}} for gaussian, binomial, and poisson
#' 

#' @details The \code{stap_glm} function is similar in syntax to 
#' \code{\link[rstanarm]{stan_glm}} except instead of performing full bayesian
#' inference for a generalized linear model stap_glm incorporates spatial 
#' information akin to the description in --need to add citation
#'@export stap_glm
stap_glm <- function(formula,
                     family = stats::gaussian(),
                     subject_data,
                     distance_data,
                     id_key = NULL,
                     D_M,
                     weights,
                     subset,
                     na.action = NULL,
                     offset = NULL,
                     model = TRUE,
                     x = FALSE,
                     y = TRUE,
                     contrasts = NULL,
                     ...,
                     prior = normal(),
                     prior_intercept = normal(),
                     prior_stap = normal(),
                     prior_theta = list(theta_one = normal(location = 1.5, scale = .5)),
                     prior_aux = cauchy(location = 0L, scale = 5L),
                     adapt_delta = NULL){

    crs_data <- extract_stap_components(formula,distance_data,
                                        subject_data, id_key, 
                                        max_distance = D_M)
    formula <- crs_data$formula
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
    }
    debug(stap_glm.fit) 
    stapfit <- stap_glm.fit(z = Z, y = Y, weights = weights,
                            dists_crs = crs_data$d_mat, u = crs_data$u,
                            max_distance = D_M,
                            offset = offset, family = family,
                            prior = prior,
                            prior_intercept = prior_intercept,
                            prior_stap = prior_stap,
                            prior_aux = prior_aux,
                            prior_theta = prior_theta,
                            adapt_delta = adapt_delta,
                            ...) ## time to debug stap_glm.fit

    sel <- apply(X, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
    X <- X[ , !sel, drop = FALSE]
    fit <- nlist(stanfit, family, formula, data, offset, weights,
                 x = X, y = Y, model = mf,  terms = mt, call,
                 na.action = attr(mf, "na.action"),
                 contrasts = attr(X, "contrasts"),
                 stan_function = "stap_glm")
    # out <- stapreg(fit)
    # out$xlevels <- .getXlevels(mt, mf)
    # if (!x) 
    #     out$x <- NULL
    # if (!y) 
    #     out$y <- NULL
    # if (!model) 
    #     out$model <- NULL
    # 
    # return(out)
}

