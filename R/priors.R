# Part of the rstap package for estimating model parameters
# Copyright (C) 2018 Trustees of the University of Michigan
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

#' Prior distributions and options
#' 
#' The functions described on this page are used to specify the prior-related #' arguments of the various modeling functions in the \pkg{rstap} package (to
#' view the priors used for an existing model see \code{\link{prior_summary}}). 
#' The default priors used in the various \pkg{rstap} modeling functions are
#' intended to be \emph{weakly informative} in that they provide moderate
#' regularlization and help stabilize computation. For many applications the
#' defaults will perform well, but prudent use of more informative priors is
#' encouraged. All of the priors here are informed by the priors in \pkg{rstanarm}, though it should be noted that the heirarchical shape priors are not included.
#' 
#' @name priors
#' @param location Prior location. In most cases, this is the prior mean, but
#'   for \code{cauchy} (which is equivalent to \code{student_t} with
#'   \code{df=1}), the mean does not exist and \code{location} is the prior
#'   median. The default value is \eqn{0}. 
#' @param scale Prior scale. The default depends on the family (see
#'   \strong{Details}).
#' @param df Prior degrees of freedom. The default is \eqn{1} for 
#'   \code{student_t}, in which case it is equivalent to \code{cauchy}.
#'   For the \code{product_normal}   prior, the degrees of freedom 
#'   parameter must be an integer (vector) that is  at least \eqn{2} (the default).
#' @param autoscale A logical scalar, defaulting to \code{TRUE}. If \code{TRUE} 
#'   then the scales of the priors on the intercept and regression coefficients 
#'   may be additionally modified internally by \pkg{rstanarm} in the following 
#'   cases. First, for Gaussian models only, the prior scales for the intercept, 
#'   coefficients, and the auxiliary parameter \code{sigma} (error standard 
#'   deviation) are multiplied by \code{sd(y)}. Additionally --- not only for 
#'   Gaussian models --- if the \code{QR} argument to the model fitting function
#'   (e.g. \code{stap_glm}) is \code{FALSE} then: for a predictor with only one 
#'   value nothing is changed; for a predictor \code{x} with exactly two unique 
#'   values, we take the user-specified (or default) scale(s) for the selected 
#'   priors and divide by the range of \code{x}; for a predictor \code{x} with 
#'   more than two unique values, we divide the prior scale(s) by \code{sd(x)}.
#'   
#' @details The details depend on the family of the prior being used:
#' \subsection{Student t family}{
#'   Family members:
#'   \itemize{
#'   \item \code{normal(location, scale)}
#'   \item \code{student_t(df, location, scale)}
#'   \item \code{cauchy(location, scale)}
#'   }
#'   Each of these functions also takes an argument \code{autoscale} which is relevant 
#'   if used for any of the non-stap related parameters. It is not used otherwise.
#'   
#'   For the prior distribution for the intercept, \code{location}, 
#'   \code{scale}, and \code{df} should be scalars. As the 
#'   degrees of freedom approaches infinity, the Student t distribution 
#'   approaches the normal distribution and if the degrees of freedom are one,
#'   then the Student t distribution is the Cauchy distribution.
#'   
#'   If \code{scale} is not specified it will default to \eqn{10} for the
#'   intercept and \eqn{2.5} for the other coefficients, unless the probit link
#'   function is used, in which case these defaults are scaled by a factor of 
#'   \code{dnorm(0)/dlogis(0)}, which is roughly \eqn{1.6}.
#'   
#'   If the \code{autoscale} argument is \code{TRUE} (the default), then the 
#'   scales will be further adjusted as described above in the documentation of 
#'   the \code{autoscale} argument in the \strong{Arguments} section.
#' }
#' \subsection{Laplace family}{
#'   Family members:
#'   \itemize{
#'   \item \code{laplace(location, scale)}
#'   \item \code{lasso(df, location, scale)}
#'   }
#'   Each of these functions also takes an argument \code{autoscale}.
#'   
#'   The Laplace distribution is also known as the double-exponential 
#'   distribution. It is a symmetric distribution with a sharp peak at its mean 
#'   / median / mode and fairly long tails. This distribution can be motivated 
#'   as a scale mixture of normal distributions and the remarks above about the 
#'   normal distribution apply here as well.
#'   
#'   The lasso approach to supervised learning can be expressed as finding the
#'   posterior mode when the likelihood is Gaussian and the priors on the 
#'   coefficients have independent Laplace distributions. It is commonplace in
#'   supervised learning to choose the tuning parameter by cross-validation,
#'   whereas a more Bayesian approach would be to place a prior on \dQuote{it},
#'   or rather its reciprocal in our case (i.e. \emph{smaller} values correspond
#'   to more shrinkage toward the prior location vector). We use a chi-square
#'   prior with degrees of freedom equal to that specified in the call to
#'   \code{lasso} or, by default, 1. The expectation of a chi-square random
#'   variable is equal to this degrees of freedom and the mode is equal to the
#'   degrees of freedom minus 2, if this difference is positive.
#'   
#'   It is also common in supervised learning to standardize the predictors 
#'   before training the model. We do not recommend doing so. Instead, it is
#'   better to specify \code{autoscale = TRUE} (the default value), which 
#'   will adjust the scales of the priors according to the dispersion in the
#'   variables. See the documentation of the \code{autoscale} argument above 
#'   and also the \code{\link{prior_summary}} page for more information.
#' }
#' \subsection{Product-normal family}{
#'   Family members:
#'   \itemize{
#'   \item \code{product_normal(df, location, scale)}
#'   }
#'   The product-normal distribution is the product of at least two independent 
#'   normal variates each with mean zero, shifted by the \code{location}
#'   parameter. It can be shown that the density of a product-normal variate is
#'   symmetric and infinite at \code{location}, so this prior resembles a
#'   \dQuote{spike-and-slab} prior for sufficiently large values of the
#'   \code{scale} parameter. For better or for worse, this prior may be
#'   appropriate when it is strongly believed (by someone) that a regression
#'   coefficient \dQuote{is} equal to the \code{location}, parameter even though
#'   no true Bayesian would specify such a prior.
#'   
#'   Each element of \code{df} must be an integer of at least \eqn{2} because
#'   these \dQuote{degrees of freedom} are interpreted as the number of normal
#'   variates being multiplied and then shifted by \code{location} to yield the
#'   regression coefficient. Higher degrees of freedom produce a sharper
#'   spike at \code{location}.
#'   
#'   Each element of \code{scale} must be a non-negative real number that is
#'   interpreted as the standard deviation of the normal variates being
#'   multiplied and then shifted by \code{location} to yield the regression
#'   coefficient. In other words, the elements of \code{scale} may differ, but
#'   the k-th standard deviation is presumed to hold for all the normal deviates
#'   that are multiplied together and shifted by the k-th element of
#'   \code{location} to yield the k-th regression coefficient. The elements of 
#'   \code{scale} are not the prior standard deviations of the regression
#'   coefficients. The prior variance of the regression coefficients is equal to
#'   the scale raised to the power of \eqn{2} times the corresponding element of
#'   \code{df}. Thus, larger values of \code{scale} put more prior volume on
#'   values of the regression coefficient that are far from zero.
#' }
#' @return A named list to be used internally by the \pkg{rstap} model
#'   fitting functions.
#' @seealso The various vignettes for the \pkg{rstanarm} and \pkg{rstap} packages also discuss 
#'   and demonstrate the use of some of the supported prior distributions.
#' 
#' #' 
#' @references
#' Gelman, A., Jakulin, A., Pittau, M. G., and Su, Y. (2008). A weakly
#' informative default prior distribution for logistic and other regression
#' models. \emph{Annals of Applied Statistics}. 2(4), 1360--1383.
#' 
#' # Can assign priors to names
#' N05 <- normal(0, 5)
#' 
#' 
NULL


#' @rdname priors
#' @export
normal <- function( location = 0, scale = NULL, autoscale = TRUE) {
    validate_parameter_value(scale)
    nlist(dist = "normal", df = NA, location, scale, autoscale)
}

#' @rdname priors
#' @export
student_t <- function(df = 1, location = 0, scale = NULL, autoscale = TRUE) {
    validate_parameter_value(scale)
    validate_parameter_value(df)
    nlist(dist = "t", df, location, scale, autoscale)
}


#' @rdname priors
#' @export
cauchy <- function(location = 0, scale = NULL, autoscale = TRUE) { 
    student_t(df = 1, location = location, scale = scale, autoscale)
}

#' @rdname priors
#' @export
laplace <- function(location = 0, scale = NULL, autoscale = TRUE) {
  nlist(dist = "laplace", df = NA, location, scale, autoscale)
}

#' @rdname priors
#' @export
lasso <- function(df = 1, location = 0, scale = NULL, autoscale = TRUE) {
  nlist(dist = "lasso", df, location, scale, autoscale)
}

#' @rdname priors
#' @export
product_normal <- function(df = 2, location = 0, scale = 1) {
  validate_parameter_value(df)
  stopifnot(all(df >= 2), all(df == as.integer(df)))
  validate_parameter_value(scale)
  nlist(dist = "product_normal", df, location, scale)
}


#' @rdname priors
#' @export
#' @param rate Prior rate for the exponential distribution. Defaults to 
#'      \code{1}. For the exponential distribution the rate parameter is the 
#'      \emph{reciprocal} of the mean.
exponential <- function(rate = 1, autoscale = TRUE) {
    stopifnot(length(rate) == 1)
    validate_parameter_value(rate)
    nlist(dist = "exponential",
          df = NA, location = NA, scale = 1/rate,
          autoscale)
}

#' @rdname priors
#' @export
log_normal <- function(location = 0, scale = 1){
    nlist(dist = "lognormal", location, scale, df = NA)
}    

#' @rdname priors
#' @export
#' @param regularization Exponent for an LKJ prior on the correlation matrix in
#'   the \code{decov} or \code{lkj} prior. The default is \eqn{1}, implying a 
#'   joint uniform prior.
#' @param concentration Concentration parameter for a symmetric Dirichlet 
#'   distribution. The default is \eqn{1}, implying a joint uniform prior.
#' @param shape Shape parameter for a gamma prior on the scale parameter in the
#'   \code{decov} prior. If \code{shape} and \code{scale} are both \eqn{1} (the
#'   default) then the gamma prior simplifies to the unit-exponential
#'   distribution.
decov <- function(regularization = 1, concentration = 1, 
                  shape = 1, scale = 1) {
  validate_parameter_value(regularization)
  validate_parameter_value(concentration)
  validate_parameter_value(shape)
  validate_parameter_value(scale)
  nlist(dist = "decov", regularization, concentration, shape, scale)
}

#' @rdname priors
#' @export
lkj <- function(regularization = 1, scale = 10, df = 1, autoscale = TRUE) {
  validate_parameter_value(regularization)
  validate_parameter_value(scale)
  validate_parameter_value(df)
  nlist(dist = "lkj", regularization, scale, df, autoscale)
}

# internal ------------------------------------------------------------------------

# Check for positive scale or df parameter (NULL ok)
#
# @param x the value to check
# @return either an error is thrown or \code{TRUE} is returned invisibly
validate_parameter_value <- function(x) {
    nm <- deparse(substitute(x))
    if (!is.null(x)) {
        if (!is.numeric(x))
            stop (nm, "should be NULL or numeric", .call = FALSE)
        if (any(x<=0))
            stop (nm, "should be positive", .call = FALSE)
    }
    invisible(TRUE)
}
