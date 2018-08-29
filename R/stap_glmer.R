# Part of the rstap package for estimating model parameters
# Copyright (c) 2018
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

#' Bayesian spatial-temporal generalized linear models with group-specific terms via Stan
#' 
#' Bayesian inference for stap-glms with group-specific coefficients that have 
#' unknown covariance matrices with flexible priors.
#' 
#' @export
#' @templateVar armRef (Ch. 11-15)
#' @templateVar fun stap_glmer, stan_lmer 
#' @templateVar pkg lme4
#' @templateVar pkgfun glmer
#' @template return-stanreg-object
#' @template see-also
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_aux
#' @template args-prior_covariance
#' @template args-prior_PD
#' @template args-algorithm
#' @template args-adapt_delta
#' @template args-QR
#' @template reference-gelman-hill
#' @template reference-muth
#' 
#' @param formula,data Same as for \code{\link[lme4]{glmer}}. 
#' Note that the subject_data argument must be provided in addition to either distance_data and/or time_data.
#' @param family Same as for \code{\link[lme4]{glmer}} except limited to gaussian, binomial and poisson 
#' @param subset,weights,offset Same as \code{\link[stats]{glm}}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but rarely 
#'   specified.
#' @param ... For \code/stap_glmer}, further arguments passed to 
#'   \code{\link[rstan]{sampling}} (e.g. \code{iter}, \code{chains}, 
#'   \code{cores}, etc.) or to \code{\link[rstan]{vb}} (if \code{algorithm} is 
#'   \code{"meanfield"} or \code{"fullrank"}). For \code{stan_lmer} and 
#'   \code{stap_glmer.nb}, \code{...} should also contain all relevant arguments
#'   to pass to \code{stap_glmer} (except \code{family}).
#'
#' @details The \code{stap_glmer} function is similar in syntax to 
#'   \code{\link[lme4]{glmer}} but rather than performing (restricted) maximum 
#'   likelihood estimation of generalized linear models, Bayesian estimation is 
#'   performed via MCMC. The Bayesian model adds priors on the 
#'   regression coefficients (in the same way as \code{\link{stan_glm}}) and
#'   priors on the terms of a decomposition of the covariance matrices of the
#'   group-specific parameters. See \code{\link{priors}} for more information
#'   about the priors.
#'   
#'   The \code{stap_lmer} function is equivalent to \code{stap_glmer} with 
#'   \code{family = gaussian(link = "identity")}. 
#'   
#'   The \code{stap_glmer.nb} function, which takes the extra argument 
#'   \code{link}, is a wrapper for \code{stap_glmer} with \code{family = 
#'   \link{neg_binomial_2}(link)}.
#'   
#'   
#' @seealso The vignette for \code{stap_glmer} ... still to be written
#'    
#' @importFrom lme4 glFormula
#' @importFrom Matrix Matrix t
stap_glmer <- 
  function(formula,
           family = gaussian,
           subject_data = NULL,
           distance_data = NULL,
           time_data = NULL,
           id_key = NULL,
           max_distance,
           subset,
           weights,
           na.action = getOption("na.action", "na.omit"),
           offset,
           contrasts = NULL,
           ...,
           prior = normal(),
           prior_intercept = normal(),
           prior_stap = normal(),
           prior_theta = log_normal(location = 1L, scale = 1L),
           prior_aux = exponential(),
           prior_covariance = decov(),
           adapt_delta = NULL) {
  
  stap_data <- extract_stap_data(formula)
  crs_data <- extract_crs_data(stap_data,
                               subject_data,
                               distance_data,
                               time_data,
                               id_key,
                               max_distance)
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  data <- validate_data(subject_data , if_missing = environment(formula))
  family <- validate_family(family)
  mc[[1]] <- quote(lme4::glFormula)
  mc$control <- make_glmerControl()
  mc$data <- subject_data
  mc$prior <- mc$prior_intercept <- mc$prior_covariance <- mc$prior_aux <-
    mc$prior_PD <- mc$algorithm <- mc$scale <- mc$concentration <- mc$shape <-
    mc$adapt_delta <- mc$... <- mc$QR <- mc$sparse <- NULL
  glmod <- eval(mc, parent.frame())
  X <- glmod$X
  y <- glmod$fr[, as.character(glmod$formula[2L])]
  if (is.matrix(y) && ncol(y) == 1L)
    y <- as.vector(y)

  offset <- model.offset(glmod$fr) %ORifNULL% double(0)
  weights <- validate_weights(as.vector(model.weights(glmod$fr)))
  if (binom_y_prop(y, family, weights)) {
    y1 <- as.integer(as.vector(y) * weights)
    y <- cbind(y1, y0 = weights - y1)
    weights <- double(0)
  }
  if (is.null(prior)) 
    prior <- list()
  if (is.null(prior_intercept)) 
    prior_intercept <- list()
  if (is.null(prior_aux)) 
    prior_aux <- list()
  if (is.null(prior_covariance))
    stop("'prior_covariance' can't be NULL.", call. = FALSE)
  group <- glmod$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_aux = prior_aux, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = group, QR = QR,  ...)
  sel <- apply(X, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  X <- X[ , !sel, drop = FALSE]
  Z <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = cbind(X, Z), y = y, data, call, terms = NULL, model = NULL,
               na.action = attr(glmod$fr, "na.action"), contrasts, algorithm, glmod, 
               stan_function = "stap_glmer")
  out <- stanreg(fit)
  class(out) <- c(class(out), "lmerMod")
  
  return(out)
}


#' @rdname stap_glmer
#' @export
stan_lmer <- 
  function(formula,
           data = NULL,
           subset,
           weights,
           na.action = getOption("na.action", "na.omit"),
           offset,
           contrasts = NULL,
           ...,
           prior = normal(),
           prior_intercept = normal(),
           prior_aux = exponential(),
           prior_covariance = decov(),
           prior_PD = FALSE,
           algorithm = c("sampling", "meanfield", "fullrank"),
           adapt_delta = NULL,
           QR = FALSE) {
  if ("family" %in% names(list(...))) {
    stop(
      "'family' should not be specified. ", 
      "To specify a family use stap_glmer instead of stan_lmer."
    )
  }
  mc <- call <- match.call(expand.dots = TRUE)
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc[[1L]] <- quote(stap_glmer)
  mc$REML <- NULL
  mc$family <- "gaussian"
  out <- eval(mc, parent.frame())
  out$call <- call
  out$stan_function <- "stan_lmer"
  return(out)
}


#' @rdname stap_glmer
#' @export
#' @param link For \code{stap_glmer.nb} only, the link function to use. See 
#'   \code{\link{neg_binomial_2}}.
#' 
stap_glmer.nb <- 
  function(formula,
           data = NULL,
           subset,
           weights,
           na.action = getOption("na.action", "na.omit"),
           offset,
           contrasts = NULL,
           link = "log",
           ...,
           prior = normal(),
           prior_intercept = normal(),
           prior_aux = exponential(),
           prior_covariance = decov(),
           prior_PD = FALSE,
           algorithm = c("sampling", "meanfield", "fullrank"),
           adapt_delta = NULL,
           QR = FALSE) {
    
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call(expand.dots = TRUE)
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc[[1]] <- quote(stap_glmer)
  mc$REML <- mc$link <- NULL
  mc$family <- neg_binomial_2(link = link)
  out <- eval(mc, parent.frame())
  out$call <- call
  out$stan_function <- "stap_glmer.nb"
  return(out)
}
