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
    if (missing(data) || is.null(data)) {
        warn_data_arg_missing()
        return(if_missing)
    }
    if (!is.data.frame(data)) {
        stop("'subject_data' must be a supplied data frame.", call. = FALSE)
    }
    
    drop_redundant_dims(data)
}
# extract_stap_components
#
# extract stap components from formula and create crs matrices
# 
# @param formula that designates model expression including stap covariates 
# @param distance_data
# @return If no error is thrown a list with the crs data matrix, index matrix u
#           and corresponding covariate names is returned
#
extract_stap_components <- function(formula, distance_data, subject_data,
                                    id_key, max_distance){
    dcol_ix <- validate_distancedata(distance_data,max_distance)
    new_formula <- get_stapless_formula(formula)
    stap_covs <- all.names(formula)[which(all.names(formula)=='stap')+1]
    stap_col_ics <- apply(distance_data, 1, function(x) which(x %in% stap_covs))
    if(!all(stap_col_ics))
        stop("The stap_covariates must all be in (only) one column
             of the distance dataframe as a character or factor variable.
             see `?stap_glm`")
    stap_col <- colnames(distance_data)[stap_col_ics[1]]
    dcol <- colnames(distance_data)[dcol_ix]
    ##ensure subjects that have zero exposure are included
    ddata <- lapply(stap_covs, function(x) distance_data[which((distance_data[,stap_col]==x &
                                                                   distance_data[,dcol]<= max_distance)),])
    if(any(lapply(ddata,nrow)==0)){
        missing <- stap_covs[which(sapply(ddata,nrow)==0)]
        stap_covs <- stap_covs[which(sapply(ddata,nrow)!=0)]
        print(paste("The following stap_covariates are not present in distance_data:",
              paste(missing,collapse = ', ')))
        print("These will be omitted from the analysis")
        ddata <- lapply(ddata,function(x) if(nrow(x)!=0) x)
        ddata[sapply(ddata,is.null)] <- NULL
    }
    M <- max(sapply(ddata, nrow))
    mddata <- lapply(ddata,function(y) merge(subject_data[,id_key], y, by = eval(id_key),
                                            all.x = T) )
    d_mat <- lapply(mddata,function(x) x[!is.na(x[,dcol]),dcol])
    d_mat <- matrix(Reduce(rbind,lapply(d_mat,function(x) if(length(x)!=M) c(x,rep(0,M)) else x)),
                    nrow = length(stap_covs), ncol = M)
    freq <- lapply(mddata, function(x) xtabs(~ get(id_key) + get(stap_col),
                                         data = x, addNA = TRUE)[,1])
    u <- lapply(freq,function(x) cbind(
        replace(dplyr::lag(cumsum(x)),
                is.na(dplyr::lag(cumsum(x))),0)+1,
                cumsum(x)))
    u <- array(Reduce(function(x,y) abind::abind(x,y,along = 2), u), 
               dim = c(nrow(subject_data), length(stap_covs), 2) )
    
    return(list(d_mat = d_mat, u = u, formula = new_formula))
}
# Validate distance_data
#
# Make sure that data is a data frame. 
#
# @param distance_data User's distance_data argument
# @return If no error is thrown, distance_column is returned
#
validate_distancedata <- function(distance_data, max_distance ) {
    if(missing(distance_data) || is.null(distance_data) || 
       !is.data.frame(distance_data)) 
        stop("distance_data dataframe must be supplied to function")
    dcol_ix <- sum(sapply(1:ncol(distance_data), function(x) all(is.double(as.matrix(distance_data[,x])))*x))
    if(dcol_ix==0)
        stop("distance_data should be a data frame with only one numeric column - see `?stap_glm`")
    if(sum(distance_data[,dcol_ix]<=max_distance)==0) 
        stop("exclusion distance results in no BEFs included in the model")
    return(dcol_ix)
}


# get_stapless_formula
#
# Get formula for typical covariates
#
# @param f formula from stap_glm
# @return formula without ~ stap() components
#
get_stapless_formula <- function(f){
    stap_ics <- which(all.names(f)=='stap')
    if(!length(stap_ics))
        stop("No Stap Covariates designated in formula")
    stap_nms <- all.names(f)[stap_ics+1]
    formula_components <- all.vars(f)[!(all.vars(f)%in%stap_nms)]
    new_f <- paste(paste0(formula_components[1],' ~ '),formula_components[2:length(formula_components)],collapse = '+')
    return(as.formula(new_f))
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
check_stapfit <- function(x) {
  if (is.list(x)) {
    if (!all(c("par", "value") %in% names(x)))
      stop("Invalid object produced please report bug")
  }
  else {
    stopifnot(is(x, "stapfit"))
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

# Throw error if object isn't a stanreg object
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

# Return the appropriate stub for variable names
#
# @param object A stanmvreg object
get_stub <- function(object) {
  if (is.jm(object)) "Long" else if (is.mvmer(object)) "y" else NULL  
} 
