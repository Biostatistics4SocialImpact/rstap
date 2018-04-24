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
    else if(prior_dist_name == "laplace") prior_dist <- 4L
    else if(prior_dist_name == "lasso") prior_dist <- 5L
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
