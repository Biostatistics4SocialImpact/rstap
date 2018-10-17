# Part of the rstap package for estimating model parameters
# Copyright (C)  2018 Trustees of University of Michigan
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# Set arguments for sampling 
#
# Prepare a list of arguments to use with \code{rstan::sampling} via
# \code{do.call}.
#
# @param object The stanfit object to use for sampling.
# @param user_dots The contents of \code{...} from the user's call to
#   the \code{stan_*} modeling function.
# @param user_adapt_delta The value for \code{adapt_delta} specified by the
#   user.
# @param prior Prior distribution list (can be NULL).
# @param ... Other arguments to \code{\link[rstan]{sampling}} not coming from
#   \code{user_dots} (e.g. \code{data}, \code{pars}, \code{init}, etc.)
# @return A list of arguments to use for the \code{args} argument for 
#   \code{do.call(sampling, args)}.
set_sampling_args <- function(object, prior, user_dots = list(), 
                              user_adapt_delta = NULL, ...) {
  args <- list(object = object, ...)
  unms <- names(user_dots)
  for (j in seq_along(user_dots)) {
    args[[unms[j]]] <- user_dots[[j]]
  }
  defaults <- default_stan_control(prior = prior,
                                   adapt_delta = user_adapt_delta)

  if (!"control" %in% unms) {
    # no user-specified 'control' argument
    args$control <- defaults
  } else {
    # user specifies a 'control' argument
    if (!is.null(user_adapt_delta)) {
      # if user specified adapt_delta argument to stan_* then
      # set control$adapt_delta to user-specified value
      args$control$adapt_delta <- user_adapt_delta
    } else {
      # use default adapt_delta for the user's chosen prior
      args$control$adapt_delta <- defaults$adapt_delta
    }
    if (is.null(args$control$max_treedepth)) {
      # if user's 'control' has no max_treedepth set it to rstanarm default
      args$control$max_treedepth <- defaults$max_treedepth
    }
  }
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
    if(any(grepl("stap",f)))
        stop("error - stap component found in amended formula - please
             report bug")
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
    if (!is.data.frame(data) || missing(data) || is.null(data)) {
        stop("'subject_data' must be a supplied data frame.", call. = FALSE)
    }

    drop_redundant_dims(data)
}


#' Validate newdata argument for posterior_predict, log_lik, etc.
#'
#' Doesn't check if the correct variables are included (that's done in pp_data),
#' just that newdata is either NULL or a data frame with no missing values. Also
#' drops any unused dimensions in variables (e.g. a one column matrix inside a
#' data frame is converted to a vector).
#' 
#' @param x User's 'newdata' argument
#' @return Either NULL or a data frame
#'
validate_newdata <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.data.frame(x)) {
    stop("If 'newdata' is specified it must be a data frame.", call. = FALSE)
  }
  if (any(is.na(x))) {
    stop("NAs are not allowed in 'newdata'.", call. = FALSE)
  }
  
  x <- as.data.frame(x)
  drop_redundant_dims(x)
}
#' Validate prediction data for posterior_predict, log_lik, etc.
#' 
#' Doesn't check if the correct variables are included
#' just that the new subject, distance and/or time data are 
#' data frames if they are included and returns a NULL value 
#' In the return list if they are not. Ensures that if a new
#' subject data is submitted, that a new distance or
#' time dataset are also included.
#' @param newsdata newsubject data
#' @param newddata
#' @param newtdata
#' @return a list with a dataframe or NULL for each of the possible new datasets
validate_predictiondata <- function(newsdata, newddata, newtdata) {

    if(is.null(newsdata) & is.null(newddata) & is.null(newtdata))
        return(NULL)

    new_subjs <- !is.null(newsdata)
    new_d <- !is.null(newddata)
    new_t <- !is.null(newtdata)
    if((new_subjs & !new_d) & (new_subjs & !new_t))
        stop("If new subject data is specified it must be submitted with either new distance data, new time data or both", .call = F)
    newsdata <- check_data_frame(newsdata,new_subjs)
    newddata <- check_data_frame(newddata,new_d)
    newtdata <- check_data_frame(newtdata,new_t)
    return(nlist(newsubjdata = newsdata, newdistdata = newddata, newtimedata = newtdata))

}

check_data_frame <- function(x,indicator){
    if(indicator){
        if(!is.data.frame(x))
            stop("If new data is specified it must be a data frame. ", call. = F)
    }
    if(any(is.na(x)))
        stop("NAs are not allowed in any specified 'newdata'.", call. = F)
    if(!indicator)
        return(NULL)
    x <- as.data.frame(x)
    drop_redundant_dims(x)
    return(x)
}

#' Validate distance_data
#'
#' Make sure that data is a data frame.
#'
#' @param distance_data User's distance_data argument
#' @return If no error is thrown, the column index
#' for the distance data is returned. If no distance_data is supplied NULL type returned.
validate_distancedata <- function(distance_data, max_distance ) {
    if(missing(distance_data)|| is.null(distance_data))
        return(NULL)
    if(!is.data.frame(distance_data) || any(is.na(distance_data)))
        stop("if distance_data is supplied it must be supplied as a dataframe with no NA values")
    num_dbl <- sum(sapply(1:ncol(distance_data),
                      function(x) all(is.double(as.matrix(distance_data[,x])))))
    if(num_dbl!=1)
        stop("distance_data should be a data frame with only one numeric column - see `?stap_glm`")
    dcol_ix <- sum(sapply(1:ncol(distance_data), function(x) all(is.double(as.matrix(distance_data[,x])))*x))
    if(sum(distance_data[,dcol_ix]<=max_distance)==0) 
        stop("exclusion distance results in no BEFs included in the model")
    return(dcol_ix)
}

#' Validate time_data
#'
#' Make sure that time_data is a data frame, return time column index.
#'
#' @param time_data User's time_data argument
#' @return If no error is thrown, the index corresponding to the column
#' holding the time data is returned. If no time_data is supplied NULL type returned.
validate_timedata <- function(time_data){
    if(missing(time_data) || is.null(time_data))
        return(NULL)
    else if(!is.data.frame(time_data))
        stop("time_data dataframe must be supplied to function")

    num_dbl <- sum(sapply(1:ncol(time_data),
                      function(x) all(is.double(as.matrix(time_data[,x])))))
    if(num_dbl!=1)
        stop("time_data should be a data frame with only one numeric column - see `?stap_glm'")
    tcol_ix <- sum(sapply(1:ncol(time_data), function(x) all(is.double(as.matrix(time_data[,x])))*x))
    return(tcol_ix)
}

#' get_stapless_formula
#'
#' Get formula for typical covariates
#'
#' @param f formula from stap_glm
#' @return formula without ~ stap() components
#'
get_stapless_formula <- function(f){
    
    with_bars <- f
    f <- lme4::nobars(f)
    stap_ics <- which(all.names(f)%in% c("stap","stap_log"))
    sap_ics <- which(all.names(f) %in% c("sap","sap_log"))
    tap_ics <- which(all.names(f) %in% c("tap","tap_log"))
    if(!length(stap_ics) & !length(sap_ics) & !length(tap_ics))
        stop("No covariates designated as 'stap','sap',or 'tap'  in formula")
    stap_nms <- all.names(f)[stap_ics + 1]
    sap_nms <- all.names(f)[sap_ics + 1]
    tap_nms <- all.names(f)[tap_ics + 1]
    not_needed <- c(stap_nms,sap_nms,tap_nms,"cexp","exp","erf","cerf","wei","cwei") 
    formula_components <- all.vars(f)[!(all.vars(f) %in% not_needed)]
    bar_components <- sapply(lme4::findbars(with_bars),paste_bars)
    formula_components <- c(formula_components,bar_components)
    if(!attr(terms(f),"intercept"))
        formula_components <- c(formula_components,"0")
    if(grepl("cbind",all.names(f))[2]){
        new_f1 <- paste0("cbind(",formula_components[1],", ",formula_components[2], ")", " ~ ")
        ix <- 3
    }
    else{
        new_f1 <- paste0(formula_components[1],' ~ ')
        ix <- 2
    }
    new_f2 <- paste(formula_components[ix:length(formula_components)],collapse = "+")
    new_f <- paste0(new_f1,new_f2)
    return(as.formula(new_f, env = environment(f)))
}

# Throw a warning if 'data' argument to modeling function is missing
stop_data_arg_missing <- function() {
    stop(
        "Omitting the 'data' argument is not allowed in rstap",
        "This is because some pre-estimation merging of bef_data and",
        "post-estimation functions (in particular 'update', 'loo', 'kfold') ",
        "are not guaranteed to work properly unless 'data' is specified as a data frame.",
        call. = FALSE
    )
}

#' Check if any variables in a model frame are constants
#' (the exception is that a constant variable of all 1's is allowed)
#'
#' @param mf A model frame or model matrix
#' @return If no constant variables are found mf is returned, otherwise an error
#' is thrown.
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

# Return names of the last dimension in a matrix/array (e.g. colnames if matrix)
#
# @param x A matrix or array
last_dimnames <- function(x) {
  ndim <- length(dim(x))
  dimnames(x)[[ndim]]
}

#' Check weights argument
#'
#' @param w The \code{weights} argument specified by user or the result of
#'   calling \code{model.weights} on a model frame.
#' @return If no error is thrown then \code{w} is returned.
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
#' @param f the \code{family} argument specified by user (or default)
#' @return If no error is thrown than either \code{f} itself is returned
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

# Grep for "b" parameters (ranef)
#
# @param x Character vector (often rownames(fit$stan_summary))
# @param ... Passed to grep
b_names <- function(x, ...) {
  grep("^b\\[", x, ...)
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
# Check that a stanfit object (or list returned by rstan::optimizing) is valid
#
check_stanfit <- function(x) {
  if (is.list(x)) {
    if (!all(c("par", "value") %in% names(x)))
      stop("Invalid object produced please report bug")
  }
  else {
    stopifnot(is(x, "stanfit"))
    if (x@mode != 0)
      stop("Invalid stapfit object produced please report bug")
  }
  return(TRUE)
}

# Default control arguments for sampling
#
# Called by set_sampling_args to set the default 'control' argument for
# \code{rstan::sampling} if none specified by user. This allows the value of
# \code{adapt_delta} to depend on the prior.
#
# @param prior Prior distribution list (can be NULL).
# @param adapt_delta User's \code{adapt_delta} argument.
# @param max_treedepth Default for \code{max_treedepth}.
# @return A list with \code{adapt_delta} and \code{max_treedepth}.
default_stan_control <- function(prior, adapt_delta = NULL,
                                 max_treedepth = 15L) {
  if (!length(prior)) {
    if (is.null(adapt_delta)) adapt_delta <- 0.95
  } else if (is.null(adapt_delta)) {
    adapt_delta <- switch(prior$dist,
                          "R2" = 0.99,
                          "hs" = 0.99,
                          "hs_plus" = 0.99,
                          "lasso" = 0.99,
                          "product_normal" = 0.99,
                          0.95) # default
  }
  nlist(adapt_delta, max_treedepth)
}

# Test if an object is a stapreg object
#
# @param x The object to test.
is.stapreg <- function(x) inherits(x, "stapreg")

# Throw error if object isn't a stapreg object
#
# @param x The object to test.
validate_stapreg_object <- function(x, call. = FALSE) {
  if (!is.stapreg(x))
    stop("Object is not a stapreg object.", call. = call.)
}

# If a is NULL (and Inf, respectively) return b, otherwise just return a
# @param a,b Objects
`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}
`%ORifINF%` <- function(a, b) {
  if (a == Inf) b else a
}

#' Extract X, Y or Z from a stapreg object
#' 
#' @keywords internal
#' @export
#' @templateVar stapregArg object
#' @template args-stapreg-object
#' @param ... Other arguments passed to methods. For a \code{stanmvreg} object
#'   this can be an integer \code{m} specifying the submodel.
#' @return For \code{get_x} and \code{get_z}, a matrix. For \code{get_y}, either
#'   a vector or a matrix, depending on how the response variable was specified.
get_y <- function(object, ...) UseMethod("get_y")

#' @rdname get_y
#' @export
get_z <- function(object, ...) UseMethod("get_z")

#' @rdname get_y
#' @export
get_x <- function(object, ...) UseMethod("get_x")

#' @rdname get_y
#' @export
get_w <- function(object, ...) UseMethod("get_w")

#' @rdname get_y
#' @export
get_z.default <- function(object, ...){
    object[["z"]] %ORifNULL% model.matrix(object)
}
#' @export
get_y.default <- function(object, ...) {
  object[["y"]] %ORifNULL% model.response(model.frame(object))
}
get_x.default <- function(object, ...)
    object[["x"]]


#' @export
get_x.lmerMod <- function(object, ...) {
  object[["x"]]
}

#' @export
get_z.lmerMod <- function(object, ...){
    object$glmod$X %ORifNULL% stop("X not found")
}

#' @export
get_w.lmerMod <- function(object, ...) {
  Zt <- object$glmod$reTrms$Zt %ORifNULL% stop("Z not found")
  t(Zt)
}

# Get inverse link function
#
# @param x A stapreg object, family object, or string. 
#   this can be an integer \code{m} specifying the submodel.
# @return The inverse link function associated with x.
linkinv <- function(x, ...) UseMethod("linkinv")
linkinv.stapreg <- function(x, ...) {
  family(x)$linkinv
}
linkinv.family <- function(x, ...) {
  x$linkinv
}

# Wrapper for rstan::summary
# @param stanfit A stanfit object created using rstan::sampling or rstan::vb
# @return A matrix of summary stats
make_stap_summary <- function(stanfit){
    levs <- c(.5, .8, .95)
    qq <- ( 1 - levs ) / 2
    probs <- sort(c(.5, qq, 1 - qq))
    rstan::summary(stanfit, probs = probs, digits = 10)$summary
}


check_reTrms <- function(reTrms) {
  stopifnot(is.list(reTrms))
  nms <- names(reTrms$cnms)
  dupes <- duplicated(nms)
  for (i in which(dupes)) {
    original <- reTrms$cnms[[nms[i]]]
    dupe <- reTrms$cnms[[i]]
    overlap <- dupe %in% original
    if (any(overlap))
      stop("rstap  does not permit formulas with duplicate group-specific terms.\n", 
           "In this case ", nms[i], " is used as a grouping factor multiple times and\n",
           dupe[overlap], " is included multiple times.\n", 
           "Consider using || or -1 in your formulas to prevent this from happening.")
  }
  return(invisible(NULL))
}

# Issue warning if high rhat values
#
# @param rhats Vector of rhat values.
# @param threshold Threshold value. If any rhat values are above threshold a
#   warning is issued.
check_rhats <- function(rhats, threshold = 1.1, check_lp = FALSE) {
  if (!check_lp)
    rhats <- rhats[!names(rhats) %in% c("lp__", "log-posterior")]

  if (any(rhats > threshold, na.rm = TRUE))
    warning("Markov chains did not converge! Do not analyze results!",
            call. = FALSE, noBreaks. = TRUE)
}


# Get the correct column name to use for selecting the median
#
# @param algorithm String naming the estimation algorithm (probably
#   \code{fit$algorithm}).
# @return Either \code{"50%"} or \code{"Median"} depending on \code{algorithm}.
select_median <- function(algorithm) {
  switch(algorithm,
         sampling = "50%",
         meanfield = "50%",
         fullrank = "50%",
         optimizing = "Median",
         stop("Bug found (incorrect algorithm name passed to select_median)",
              call. = FALSE))
}

# Methods for creating linear predictor
#
# Make linear predictor vector from x and point estimates for delta and beta
# or linear predictor matrix from x and full posterior sample of delta and beta.
#
# @param delta_beta A vector or matrix of parameter estimates
# @param x Predictor matrix.
# @param offset Optional offset vector.
# @return A vector or matrix.
linear_predictor <- function(delta_beta, x, offset = NULL) {
  UseMethod("linear_predictor")
}
linear_predictor.default <- function(delta_beta, x, offset = NULL) {
  eta <- as.vector(if (NCOL(x) == 1L) x * delta_beta else x %*% delta_beta)
  if (length(offset))
    eta <- eta + offset

  return(eta)
}
linear_predictor.matrix <- function(delta_beta, x, offset = NULL) {
  if (NCOL(delta_beta) == 1L)
    delta_beta <- as.matrix(delta_beta)
  eta <- delta_beta %*% t(x)
  if (length(offset))
    eta <- sweep(eta, 2L, offset, `+`)

  return(eta)
}

# Regex parameter selection
#
# @param x stapreg object
# @param regex_pars Character vector of patterns
grep_for_pars <- function(x, regex_pars) {
  validate_stapreg_object(x)
  stopifnot(is.character(regex_pars))
  out <- unlist(lapply(seq_along(regex_pars), function(j) {
    grep(regex_pars[j], rownames(x$stap_summary), value = TRUE) 
  }))
  if (!length(out))
    stop("No matches for 'regex_pars'.", call. = FALSE)
  
  return(out)
}

# Combine pars and regex_pars
#
# @param x stapreg object
# @param pars Character vector of parameter names
# @param regex_pars Character vector of patterns
collect_pars <- function(x, pars = NULL, regex_pars = NULL) {
  if (is.null(pars) && is.null(regex_pars))
    return(NULL)
  if (!is.null(pars))
    pars[pars == "varying"] <- "b"
  if (!is.null(regex_pars))
    pars <- c(pars, grep_for_pars(x, regex_pars))
  unique(pars)
}
# Test if stapreg object used stan_(g)lmer
#
# @param x A stapreg object.
is.mer <- function(x) {
  stopifnot(is.stapreg(x))
  check1 <- inherits(x, "lmerMod")
  check2 <- !is.null(x$glmod)
  if (check1 && !check2) {
    stop("Bug found. 'x' has class 'lmerMod' but no 'glmod' component.")
  } else if (!check1 && check2) {
    stop("Bug found. 'x' has 'glmod' component but not class 'lmerMod'.")
  }
  isTRUE(check1 && check2)
}

#' @importFrom lme4 glmerControl
make_glmerControl <- function(...) {
  glmerControl(check.nlev.gtreq.5 = "ignore",
               check.nlev.gtr.1 = "stop",
               check.nobs.vs.rankZ = "ignore",
               check.nobs.vs.nlev = "ignore",
               check.nobs.vs.nRE = "ignore", ...)  
}

# Check if a fitted model (stapreg object) has weights
# 
# @param x stapreg object
# @return Logical. Only TRUE if x$weights has positive length and the elements
#   of x$weights are not all the same.
#
model_has_weights <- function(x) {
  wts <- x[["weights"]]
  if (!length(wts)) {
    FALSE
  } else if (all(wts == wts[1])) {
    FALSE
  } else {
    TRUE
  }
}

# Test if stapreg object used stan_nlmer
#
# @param x A stanreg object.
is.nlmer <- function(x) {
  is.mer(x) && inherits(x, "nlmerMod")
}

# Get the posterior sample size
#
# @param x A stapreg object
# @return the posterior sample size (or size of sample from approximate posterior)
posterior_sample_size <- function(x) {
 validate_stapreg_object(x)
 pss <- x$stapfit@sim$n_save
  sum(pss - x$stapfit@sim$warmup2)
}

paste_scale <- function(names)
    paste0(names,"_scale")

get_weight_function <- function(weight_code){
    switch(weight_code,function(x,y) { pracma::erf(x/y)} ,
           function(x,y){ pracma::erfc(x/y)},
           function(x,y){ exp(-x/y)}, 
           function(x,y){1- exp(-x/y)},
           function(x,y){exp(- (x/y)^y)},
           function(x,y){1 - exp(-(x/y)^y)})
}

paste_bars <- function(bar_element){
    l <- length(all.names(bar_element))
    all_names <- all.names(bar_element)
    if(l==2)
        paste0("(1",all_names[1],all_names[2],")")
    else
        paste0("(",paste(all_names[2:(l-1)],collapse = "+"),all_names[1],all_names[l], ")")
}
