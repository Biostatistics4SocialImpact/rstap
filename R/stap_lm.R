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

#' Bayesian spatial-temporal aggregated predictor(STAP) models via Stan
#'
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Bayesian inference for linear modeling with regularizing priors on the model
#' parameters that are driven by prior beliefs about \eqn{R^2}, the proportion
#' of variance in the outcome attributable to the predictors. See
#' \code{\link{priors}} for an explanation of this critical point.
#' \code{\link{stan_glm}} with \code{family="gaussian"} also estimates a linear
#' model with normally-distributed errors and allows for various other priors on
#'
#' @export
#' @templateVar fun stan_lm
#' @templateVar fitfun stan_lm.fit or stan_lm.wfit
#' @templateVar pkg stats
#' @templateVar pkgfun lm
#' @templateVar rareargs model,offset,weights
#' @templateVar rareargs2 na.action,singular.ok,contrasts
#' @template return-stapreg-object
#' @template return-stapfit-object
#' @template args-formula-data-subset
#' @template args-same-as-rarely
#' @template args-same-as-rarely-2
#' @template args-x-y
#' @template args-dots
#' @template args-adapt_delta
#' 
#' @param w Same as in \code{\link[stats]{lm.wfit}} but rarely specified.
#' @param prior Must be a call to \code{\link{R2}} with its 
#'   \code{location} argument specified or \code{NULL}, which would
#'   indicate a standard uniform prior for the \eqn{R^2}.
#' @param prior_intercept Either \code{NULL} (the default) or a call to
#'   \code{\link{normal}}. If a \code{\link{normal}} prior is specified
#'   without a \code{scale}, then the standard deviation is taken to be
#'   the marginal standard deviation of the outcome divided by the square
#'   root of the sample size, which is legitimate because the marginal
#'   standard deviation of the outcome is a primitive parameter being
#'   estimated.
#'   
#'   \strong{Note:} The prior distribution for the intercept is set so it
#'   applies to the value \emph{when all predictors are centered}. If you prefer
#'   to specify a prior on the intercept without the predictors being
#'   auto-centered, then you have to omit the intercept from the
#'   \code{\link[stats]{formula}} and include a column of ones as a predictor,
#'   in which case some element of \code{prior} specifies the prior on it,
#'   rather than \code{prior_intercept}. Regardless of how
#'   \code{prior_intercept} is specified, the reported \emph{estimates} of the
#'   intercept always correspond to a parameterization without centered
#'   predictors (i.e., same as in \code{glm}).
#'
#'
#' @details The \code{stan_lm} function is similar in syntax to the 
#'   \code{\link[stats]{lm}} function but rather than choosing the parameters to
#'   minimize the sum of squared residuals, samples from the posterior 
#'   distribution are drawn using MCMC. The \code{stan_lm} function 
#'   has a formula-based interface and would usually be called by users but
#'   the \code{stap_lm.fit} and \code{stap_lm.wfit} functions might be called 
#'   by other functions that parse the data themselves and are analgoous to
#'   \code{\link[stats]{lm.fit}} and \code{\link[stats]{lm.wfit}} respectively.
#'      
#'   In addition to estimating \code{sigma} --- the standard deviation of the
#'   normally-distributed errors --- this model estimates a positive parameter
#'   called \code{log-fit_ratio}. If it is positive, the marginal posterior 
#'   variance of the outcome will exceed the sample variance of the outcome
#'   by a multiplicative factor equal to the square of \code{fit_ratio}.
#'   Conversely if \code{log-fit_ratio} is negative, then the model underfits.
#'   Given the regularizing nature of the priors, a slight underfit is good.
#'   
#'   Finally, the posterior predictive distribution is generated with the
#'   predictors fixed at their sample means. This quantity is useful for
#'   checking convergence because it is reasonably normally distributed
#'   and a function of all the parameters in the model.
#'   
#'   
#' @references 
#' Lewandowski, D., Kurowicka D., and Joe, H. (2009). Generating random
#' correlation matrices based on vines and extended onion method. 
#' \emph{Journal of Multivariate Analysis}. \strong{100}(9), 1989--2001.
#' 
#' @seealso 
#' The vignettes for \code{stap_lm}, which has more
#' thorough descriptions and examples.
#' \url{https://biostatistics4socialimpact.github.io}
#' 
#' Also see \code{\link{stap_glm}}, which --- if \code{family =
#' gaussian(link="identity")} --- also estimates a stap linear model with
#' normally-distributed errors but specifies different priors.
#'   
#'   
stap_lm <- function(formula, subject_data, bef_data,
                    id_key = NULL, max_distance,
                    subset, weights, na.action,
                    model = TRUE,
                    x = FALSE, y = FALSE, z = FALSE,
                    singular.ok = TRUE, contrasts = NULL, offset, ...,
                    prior = R2(stop("'location' must be specified'")),
                    prior_intercept = NULL,
                    adapt_delta = NULL){

    crs_data <- extract_stap_components(formula, distance_data,
                                        subject_data, id_key,
                                        max_distance)
}
