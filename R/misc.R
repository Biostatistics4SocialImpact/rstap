#' Set arguments for sampling
#'
#' Prepare a list of arguments to use with \code{rstan::sampling} via
#' \code{do.call}.
#'
#' @param object The stanfit object to use for sampling.
#' @param user_dots The contents of \code{...} from the user's call to
#'   the \code{stan_*} modeling function.
#' @param user_adapt_delta The value for \code{adapt_delta} specified by the
#'   user.
#' @param prior Prior distribution list (can be NULL).
#' @param ... Other arguments to \code{\link[rstan]{sampling}} not coming from
#'   \code{user_dots} (e.g. \code{data}, \code{pars}, \code{init}, etc.)
#' @return A list of arguments to use for the \code{args} argument for
#'   \code{do.call(sampling, args)}.
set_sampling_args <- function(object, prior, user_dots = list(),
                              user_adapt_delta = NULL, user_max_tree_depth = 10, ...) {
    args <- list(object = object, ...)
    unms <- names(user_dots)
    for (j in seq_along(user_dots)) {
        args[[unms[j]]] <- user_dots[[j]]
    }
    defaults <- list(adapt_delta = .85, max_treedepth = 18)
    args$control$adapt_delta <- defaults$adapt_delta
    args$control$max_treedepth <- defaults$max_treedepth

    args$save_warmup <- FALSE
    return(args)
}


#' create a named list using specified names or, if names are omitted using the
#' names of the objects in the list
#
# @param ... Objects to include in the list.
# @return A named list.
nlist <- function(...){
    m <- match.call()
    out <- list(...)
    no_names <- is.null(names(out))
    has_name <- if(no_names) FALSE else nzchar(names(out))
    if (all(has_name))
        return(out)
    nms <- as.character(m)[-1L]
    if(no_names) {
        names(out) <- nms
    } else{
        names(out)[!has_name] <- nms[!has_name]
   }
   return(out)
}



#' assign appropriate numeric coding for specified distribution
#'
#' @param prior the named list returned by the specific prior function
#' this code is taken from the rstanarm package
#' @return list with all possible parameter values appropriately populated
assign_dist <- function(prior){
    if(!is.list(prior))
        stop(sQuote(deparse(substitute(prior))), "should be a named list")

    prior_dist_name <- prior$dist
    prior_scale <- prior$scale
    prior_mean <- prior$location
    prior_df <- prior$df
    prior_mean [is.na(prior_mean)] <- 0
    prior_df[is.na(prior_df)] <- 1
    if(prior_dist_name == "normal") prior_dist <- 1L
    else if(prior_dist_name == "t") prior_dist <- 2L
    else if(prior_dist_name == "cauchy") prior_dist <- 3L
    else if(prior_dist_name == "lognormal") prior_dist <- 4L
    else if(prior_dist_name == "beta") prior_dist <- 5L
    else if(prior_dist_name == "product_normal") prior_dist <- 6L
    else if(prior_dist_name == NA){
        prior_scale <- 0L
        prior_dist <- 99L ## won't be used
    }
    else prior_dist <- 8L

    nlist(prior_dist, prior_mean,
          prior_scale, prior_df,
         prior_dist_name)
}



# Convert 2-level factor to 0/1
fac2bin <- function(y) {
  if (!is.factor(y)) 
    stop("Bug found: non-factor as input to fac2bin.", 
         call. = FALSE)
  if (!identical(nlevels(y), 2L)) 
    stop("Bug found: factor with nlevels != 2 as input to fac2bin.", 
         call. = FALSE)
  as.integer(y != levels(y)[1L])
}


# Check for stap_glmer syntax in formula for non-glmer models
#
# @param f the model \code{formula}.
# @return Nothin is returned but an error might be thrown
validate_glm_formula <- function(f) {
    if (any(grepl("\\|", f)))
        stop("Using '|' in model formula not allowed. ",
             "Maybe you meant to use 'stap_(g)lmer'?", call. = FALSE)
}

# Validate data argument
#
# Make sure that, if specified, data is a data frame. If data is not missing
# then dimension reduction is also performed on variables (i.e., a one column
# matrix inside a data frame is converted to a vector).
#
# @param data User's data argument
# @param if_missing Object to return if data is missing/null
# @return If no error is thrown, data itself is returned if not missing/null,
#   otherwise if_missing is returned.
#
drop_redundant_dims <- function(data) {
    drop_dim <- sapply(data, function(v) is.matrix(v) && NCOL(v) == 1)
    data[, drop_dim] <- lapply(data[, drop_dim, drop=FALSE], drop)
    return(data)
}
validate_data <- function(data, if_missing = NULL) {
    if (missing(data) || is.null(data)) {
        warn_data_arg_missing()
        return(if_missing)
    }
    if (!is.data.frame(data)) {
        stop("'data' must be a data frame.", call. = FALSE)
    }
    
    drop_redundant_dims(data)
}

# Throw a warning if 'data' argument to modeling function is missing
warn_data_arg_missing <- function() {
    warning(
        "Omitting the 'data' argument is not allowed in rstap",
        "This is because some post-estimation functions (in particular 'update', 'loo', 'kfold') ", 
        "are not guaranteed to work properly unless 'data' is specified as a data frame.",
        call. = FALSE
    )
}

# Check if any variables in a model frame are constants
# (the exception is that a constant variable of all 1's is allowed)
# 
# @param mf A model frame or model matrix
# @return If no constant variables are found mf is returned, otherwise an error
#   is thrown.
check_constant_vars <- function(mf) {
    # don't check if columns are constant for binomial
    mf1 <- if (NCOL(mf[, 1]) == 2) mf[, -1, drop=FALSE] else mf
    
    lu1 <- function(x) !all(x == 1) && length(unique(x)) == 1
    nocheck <- c("(weights)", "(offset)", "(Intercept)")
    sel <- !colnames(mf1) %in% nocheck
    is_constant <- apply(mf1[, sel, drop=FALSE], 2, lu1)
    if (any(is_constant)) {
        stop("Constant variable(s) found: ", 
             paste(names(is_constant)[is_constant], collapse = ", "), 
             call. = FALSE)
    }
    return(mf)
}

# Check weights argument
# 
# @param w The \code{weights} argument specified by user or the result of
#   calling \code{model.weights} on a model frame.
# @return If no error is thrown then \code{w} is returned.
validate_weights <- function(w) {
    if (missing(w) || is.null(w)) {
        w <- double(0)
    } else {
        if (!is.numeric(w)) 
            stop("'weights' must be a numeric vector.", 
                 call. = FALSE)
        if (any(w < 0)) 
            stop("Negative weights are not allowed.", 
                 call. = FALSE)
    }
    return(w)
}

# Check offset argument
#
# @param o The \code{offset} argument specified by user or the result of calling
#   \code{model.offset} on a model frame.
# @param y The result of calling \code{model.response} on a model frame.
# @return If no error is thrown then \code{o} is returned.
validate_offset <- function(o, y) {
    if (is.null(o)) {
        o <- double(0)
    } else {
        if (length(o) != NROW(y))
            stop(gettextf("Number of offsets is %d but should be %d (number of observations)",
                          length(o), NROW(y)), domain = NA, call. = FALSE)
    }
    return(o)
}


#' Check family argument
#
# @param f the \code{family} argument specified by user (or default)
#'@return If no error is thrown than either \code{f} itself is returned
#' (if already a family) or the family object created from \code{f} is 
#' returned if \code{f} a string or function. Code adapted from \pkg{rstanarm}.
validate_family <-  function(f) {
    if(is.character(f))
        f <-  get(f, mode = 'function', envir = parent.frame(2))
    if(is.function(f))
        f  <- f()
    if(!is(f,'family'))
        stop("'family' must be a family.", call. = F)
    
    return(f)
}

# Maybe broadcast 
#
# @param x A vector or scalar.
# @param n Number of replications to possibly make. 
# @return If \code{x} has no length the \code{0} replicated \code{n} times is
#   returned. If \code{x} has length 1, the \code{x} replicated \code{n} times
#   is returned. Otherwise \code{x} itself is returned.
maybe_broadcast <- function(x, n) {
  if (!length(x)) {
    rep(0, times = n)
  } else if (length(x) == 1L) {
    rep(x, times = n)
  } else {
    x
  }
}
# Check and set scale parameters for priors
#
# @param scale Value of scale parameter (can be NULL).
# @param default Default value to use if \code{scale} is NULL.
# @param link String naming the link function or NULL.
# @return If a probit link is being used, \code{scale} (or \code{default} if
#   \code{scale} is NULL) is scaled by \code{dnorm(0) / dlogis(0)}. Otherwise
#   either \code{scale} or \code{default} is returned.
set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link) || is.null(link))
  if (is.null(scale)) 
    scale <- default
  if (isTRUE(link == "probit"))
    scale <- scale * dnorm(0) / dlogis(0)
  
  return(scale)
}

# If y is a 1D array keep any names but convert to vector (used in stan_glm)
#
# @param y Result of calling model.response
array1D_check <- function(y) {
  if (length(dim(y)) == 1L) {
    nms <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nms)) 
      names(y) <- nms
  }
  return(y)
}
# Check for a binomial model with Y given as proportion of successes and weights 
# given as total number of trials
# 
binom_y_prop <- function(y, family, weights) {
  if (!is.binomial(family$family)) 
    return(FALSE)

  yprop <- NCOL(y) == 1L && 
    is.numeric(y) && 
    any(y > 0 & y < 1) && 
    !any(y < 0 | y > 1)
  if (!yprop)
    return(FALSE)
  
  wtrials <- !identical(weights, double(0)) && 
    all(weights > 0) && 
    all(abs(weights - round(weights)) < .Machine$double.eps^0.5)
  isTRUE(wtrials)
}
# Test for a given family
#
# @param x A character vector (probably x = family(fit)$family)
is.binomial <- function(x) x == "binomial"
is.gaussian <- function(x) x == "gaussian"
is.gamma <- function(x) x == "Gamma"
is.ig <- function(x) x == "inverse.gaussian"
is.nb <- function(x) x == "neg_binomial_2"
is.poisson <- function(x) x == "poisson"
is.beta <- function(x) x == "beta" || x == "Beta regression"
