#' Summarize the priors used for an rstap model
#' 
#' The \code{prior_summary} method provides a summary of the prior distributions
#' used for the parameters in a given model. In some cases the user-specified
#' prior does not correspond exactly to the prior used internally by
#' \pkg{rstap} (see the sections below). Especially in these cases, but also
#' in general, it can be much more useful to visualize the priors. 
#' @aliases prior_summary
#' @export
#'
#' @templateVar stapregArg object
#' @template args-stapreg-object
#' @param digits Number of digits to use for rounding.
#' @param ... Currently ignored by the method for stapreg objects.
#' 
#' @section Intercept (after predictors centered): 
#'   For \pkg{rstap} modeling functions that accept a \code{prior_intercept} 
#'   argument, the specified prior for the intercept term applies to the 
#'   intercept after \pkg{rstap} internally centers the predictors so they 
#'   each have mean zero. The estimate of the intercept returned to the user 
#'   correspond to the intercept with the predictors as specified by the user 
#'   (unmodified by \pkg{rstap}), but when \emph{specifying} the prior the 
#'   intercept can be interpreted as the expected outcome when the predictors are
#'   set to their means.    
#'
#' @section Adjusted scales: For some models you may see "\code{adjusted scale}"
#'   in the printed output and adjusted scales included in the object returned 
#'   by \code{prior_summary}. These adjusted scale values are the prior scales 
#'   actually used by \pkg{rstap} and are computed by adjusting the prior 
#'   scales specified by the user to account for the scales of the predictors 
#'   (as described in the documentation for the \code{\link[=priors]{autoscale}}
#'   argument). To disable internal prior scale adjustments set the 
#'   \code{autoscale} argument to \code{FALSE} when setting a prior using one of
#'   the distributions that accepts an \code{autoscale} argument. For example,
#'   \code{normal(0, 5, autoscale=FALSE)} instead of just \code{normal(0, 5)}.
#'   
#' @return A list of class "prior_summary.stapreg", which has its own print
#'   method.
#'   
#' @seealso The \link[=priors]{priors help page} and the \emph{Prior
#'   Distributions} vignette from the \pkg{rstanarm} package.
#' 
prior_summary.stapreg <- function(object, digits = 2,...) {
  x <- object[["prior.info"]]
  if (is.null(x)) {
    message("Priors not found in stapreg object.")
    return(invisible(NULL))
  }  
  structure(x, class = "prior_summary.stapreg",
            model_name = deparse(substitute(object)), 
            stan_function = object$stan_function,
            print_digits = digits)
  
}

#' @export
#' @method print prior_summary.stapreg
print.prior_summary.stapreg <- function(x, digits, ...) {
  if (missing(digits))
    digits <- attr(x, "print_digits") %ORifNULL% 2
  .dig <- digits
  .fr2 <- function(y, .digits = .dig, ...) format(y, digits = .digits, ...)
  .fr3 <- function(y, .nsmall = .dig) .fr2(y, nsmall = .nsmall)
  formatters <- list(.fr2, .fr3)
  model_name <- attr(x, "model_name")
  stan_function <- attr(x, "stan_function")
  
  msg <- paste0("Priors for model '", model_name, "'")
  cat(msg, "\n------")
  
  if (!is.null(x[["prior_intercept"]]))
      .print_scalar_prior(
        x[["prior_intercept"]], 
        txt = paste0("Intercept", " (after predictors centered)"), 
        formatters
      )
  if (!is.null(x[["prior"]]))
      .print_vector_prior(
        x[["prior"]], 
        txt = paste0("\nCoefficients")  , 
        formatters = formatters
      )
  if(!is.null(x[["prior_stap"]]))
      .print_vector_prior(
       x[["prior_stap"]],
       txt = "Stap Coefficients",
       formatters = formatters
       )
  if(!is.null(x[["prior_theta"]]))
      .print_vector_prior(
       x[["prior_theta"]],
       txt = "Stap Scales",
       formatters = formatters
       )

  if (!is.null(x[["prior_aux"]])) {
      aux_name <- x[["prior_aux"]][["aux_name"]]
      aux_dist <- x[["prior_aux"]][["dist"]]
      if (aux_dist %in% c("normal", "student_t", "cauchy"))
          x[["prior_aux"]][["dist"]] <- paste0("half-", aux_dist)
          .print_scalar_prior(
            x[["prior_aux"]], 
            txt = paste0("\nAuxiliary (", aux_name, ")"), 
            formatters
          )
    }    
  
  # unique to stap_(g)lmer
  if (!is.null(x[["prior_covariance"]]))
    .print_covariance_prior(x[["prior_covariance"]], txt = "\nCovariance", formatters)
  

  cat("\n------\n")
  cat("See help('prior_summary.stapreg') for more details\n")
  invisible(x)
}


# internal ----------------------------------------------------------------


# 
# @param x numeric vector
# @param formatter a formatting function to apply (see .fr2, .fr3 above)
# @param N the maximum number of values to include before replacing the rest
#   with '...'
.format_pars <- function(x, formatter, N = 3) {
  K <- length(x)
  if (K < 2)
    return(x)
  paste0(
    "[", 
    paste(c(formatter(x[1:min(N, K)]), if (N < K) "..."), 
          collapse = ","), 
    "]"
  )
}

# Print priors for intercept/coefs (called internally by print.prior_summary.stapreg)
#
# @param p named list of prior stuff
# @param txt header to be printed
# @param formatters a list of two formatter functions like .fr2, .fr3 (defined
#   in prior_summary.stapreg). The first is used for format all numbers except
#   for adjusted scales, for which the second function is used. This is kind of
#   hacky and should be replaced at some point.
# 
.print_scalar_prior <- function(p, txt = "Intercept", formatters = list()) {
  stopifnot(length(formatters) == 2)
  .f1 <- formatters[[1]]
  .f2 <- formatters[[2]]
  cat(paste0("\n", txt, "\n ~"),
      if (is.na(p$dist)) {
        "flat"
      } else if (p$dist == "exponential") {
        paste0(p$dist,"(rate = ", .f1(p$rate), ")")
      } else { # normal, student_t, cauchy
        if (is.null(p$df)) {
          paste0(p$dist,"(location = ", .f1(p$location), 
                 ", scale = ", .f1(p$scale),")")
        } else {
          paste0(p$dist, "(df = ", .f1(p$df), 
                 ", location = ", .f1(p$location), 
                 ", scale = ", .f1(p$scale), ")")
        }
      }
  )
  if (!is.null(p$adjusted_scale))
    cat("\n     **adjusted scale =", .f2(p$adjusted_scale), 
        if (p$dist == "exponential") ("(adjusted rate = 1/adjusted scale)"))
}
.print_vector_prior <- function(p, txt = "Coefficients", formatters = list()) {
  stopifnot(length(formatters) == 2)
  .f1 <- formatters[[1]]
  .f2 <- formatters[[2]]
  
  if (!(p$dist %in% c("R2", NA))) {
    if (p$dist %in% c("normal", "student_t", "cauchy", "laplace", "lasso", "product_normal", "lognormal")) {
      p$location <- .format_pars(p$location, .f1)
      p$scale <- .format_pars(p$scale, .f1)
      if (!is.null(p$df))
        p$df <- .format_pars(p$df, .f1)
      if (!is.null(p$adjusted_scale))
        p$adjusted_scale <- .format_pars(p$adjusted_scale, .f2)
    } else if (p$dist %in% c("hs_plus")) {
      p$df1 <- .format_pars(p$df, .f1)
      p$df2 <- .format_pars(p$scale, .f1)
    } else if (p$dist %in% c("hs")) {
      p$df <- .format_pars(p$df, .f1)
    } else if (p$dist %in% c("product_normal"))
      p$df <- .format_pars(p$df, .f1)
  }
  cat(paste0("\n", txt, "\n ~"),
      if (is.na(p$dist)) {
        "flat"
      } else if (p$dist %in% c("normal", "student_t", "cauchy", 
                               "laplace", "lasso", "product_normal","lognormal")) {
        if (is.null(p$df)) {
          paste0(p$dist, "(location = ", .f1(p$location), 
                 ", scale = ", .f1(p$scale), ")")
        } else {
          paste0(p$dist, "(df = ", .f1(p$df), 
                 ", location = ", .f1(p$location), 
                 ", scale = ", .f1(p$scale),")")
        }
      } else if (p$dist %in% c("hs_plus")) {
        paste0("hs_plus(df1 = ", .f1(p$df1), ", df2 = ", .f1(p$df2), ")")
      } else if (p$dist %in% c("hs")) {
        paste0("hs(df = ", .f1(p$df), ")")
      } else if (p$dist %in% c("R2")) {
        paste0("R2(location = ", .f1(p$location), ", what = '", p$what, "')")
      })
  
  if (!is.null(p$adjusted_scale)) cat("\n     **adjusted scale =", .f2(p$adjusted_scale))
}
.print_covariance_prior <- function(p, txt = "Covariance", formatters = list()) {
  .f1 <- formatters[[1]]
  p$regularization <- .format_pars(p$regularization, .f1)
  p$concentration <- .format_pars(p$concentration, .f1)
  p$shape <- .format_pars(p$shape, .f1)
  p$scale <- .format_pars(p$scale, .f1)
  cat(paste0("\n", txt, "\n ~"),
      paste0(p$dist, "(",  
             "reg. = ", .f1(p$regularization),
             ", conc. = ", .f1(p$concentration), ", shape = ", .f1(p$shape),
             ", scale = ", .f1(p$scale), ")")
  )
}

