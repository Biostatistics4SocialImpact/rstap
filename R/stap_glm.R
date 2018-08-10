#' Fitting Generalized Linear STAP models
#'
#' @param formula 
#' @param family Same as \code{\link[stats]{glm}} for gaussian, binomial, and poisson
#' @param subject_data
#' @param distance_data
#' @param time_data
#' @param id_key name of column to join on between subject_data and bef_data
#' @param max_distance the inclusion distance; upper bound for all elements of dists_crs
#' @param weights 
#' @details The \code{stap_glm} function is similar in syntax to 
#' \code{\link[rstanarm]{stan_glm}} except instead of performing full bayesian
#' inference for a generalized linear model stap_glm incorporates spatial 
#' as detailed in in --need to add citation --
#'@export stap_glm
stap_glm <- function(formula,
                     family = gaussian(),
                     subject_data,
                     distance_data,
                     time_data,
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
                     prior_theta = normal(location = 30, scale = 10),
                     prior_aux = cauchy(location = 0L, scale = 5L),
                     adapt_delta = NULL){
    stap_data <- extract_stap_data(formula)
    Q_t <- lapply(stap_data,function(x) x$stap_type =='temporal')
    Q_s <- lapply(stap_data,function(x) x$stap_type =='spatial')
    Q_st <- lapply(stap_data,function(x) x$stap_type =='spatial-temporal')
    crs_data <- extract_crs_data(stap_data,
                                 distance_data,
                                 time_data,
                                 id_key)
    original_formula <- formula
    stapless_formula <- get_stapless_formula(formula)
    family <- validate_family(family)
    validate_glm_formula(formula)
    subject_data <- validate_data(subject_data, if_missing = environment(stapless_formula))
    call <- match.call(expand.dots = TRUE)
    mf <-  match.call(expand.dots = FALSE)
    mf$formula <- stapless_formula
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
    stapfit <- stap_glm.fit(z = Z, y = Y, weights = weights,
                            dists_crs = stap_data$d_mat, u_s = stap_data$u_s,
                            times_crs = stap_data$t_mat, u_t = stap_data$u_t,
                            weight_functions = stap_data$w,
                            stap_code = stap_data$stap_code,
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
    fit <- nlist(stapfit, family,
                 formula = original_formula,
                 subject_data,
                 bef_data,
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

