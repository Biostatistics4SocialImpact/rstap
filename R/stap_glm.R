# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3 # of the License, or (at your option) any later version.
# # This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Bayesian generalized spatial-temporal aggregated predictor(STAP) models via Stan
#'
#' Generalized linear modeling with spatial temporal aggregated predictors using
#' prior distributions for the coefficients, intercept, spatial-temporal scales, and auxiliary parameters.
#'
#' @export
#' @templateVar pkg stats
#' @templateVar pkgfun glm
#' @templateVar sameargs offset,weights
#' @templateVar rareargs contrasts
#' @templateVar fun stap_glm
#' @templateVar fitfun stan_glm.fit
#' @template return-stapreg-object
#' @template return-stapfit-object
#' @template see-also
#' @template args-same-as
#' @template args-same-as-rarely
#' @template args-dots
#' @template reference-gelman-hill
#' @template reference-muth
#'
#' @param formula Similar as for \code{\link[stats]{glm}} with the addition of \code{stap},\code{sap} \code{tap} terms as needed
#' @param family Same as \code{\link[stats]{glm}} for gaussian, binomial, and poisson families.
#' @param benvo built environment \code{\link[rbenvo]{benvo}} object from the \code{rbenvo} package containing the relevant subject-BEF data
#' @param prior_list a call to \code{\link{stap_prior}}
#' @details The \code{stap_glm} function is similar in syntax to
#' \code{\link[rstanarm]{stan_glm}} except instead of performing full bayesian
#' inference for a generalized linear model stap_glm incorporates spatial-temporal covariates
#' @seealso The various vignettes for \code{stap_glm} at
#'   \url{https://biostatistics4socialimpact.github.io/rstap/articles} and the \href{http://arxiv.org/abs/1812.10208}{preprint} article.  
#'
#'@export stap_glm
#'@examples
#'
#' fit_glm <- stap_glm(formula = y ~ sex + sap(Fast_Food),
#'					   benvo = FFR_example,
#'                      family = gaussian(link = 'identity'),
#'						prior = stap_priors()
#'                      chains = 1, iter = 300, # for speed of example only
#'                      refresh = -1, verbose = FALSE) 
#'
stap_glm <- function(formula,
                     family = gaussian(),
					 benvo,
					 prior_list = stap_priors(),
                     ...
                     ){

    spec <- stap_specification(formula,benvo)
    family <- validate_family(family)
    validate_glm_formula(spec$stapless_formula)
	mf <- rbenvo::subject_design(benvo,spec$stapless_formula)

    stapfit <- stap_glm.fit(y = mf$y, z = mf$X, 
                            dists_crs = spec$d_mat,
                            u_s = spec$u_s,
                            times_crs = spec$t_mat,
                            u_t = spec$u_t,
							stapmat = spec$stapinfo,
							family = family,
                            prior_list = prior_list,
							group = NULL,
                            ...)

    sel <- apply(Z, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
    Z <- Z[ , !sel, drop = FALSE]
    fit <- nlist(stapfit, family,
                 formula = original_formula,
				 spec = spec,
                 subject_data,
                 distance_data,
                 time_data,
                 dists_crs = crs_data$d_mat,
                 times_crs = crs_data$t_mat,
                 u_s = crs_data$u_s,
                 u_t = crs_data$u_t,
                 model = mf, call,
                 stan_function = "stap_glm")
    out <- stapreg(fit)
    return(out)
}

#' @rdname stap_glm
#' @export
stap_lm <- 
  function(formula,
		   benvo,
		   prior = stap_prior(),
           ...){

  if ("family" %in% names(list(...))) {
    stop(
      "'family' should not be specified. ", 
      "To specify a family use stap_glmer instead of stap_lmer."
    )
  }
  mc <- call <- match.call(expand.dots = TRUE)
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc[[1L]] <- quote(stap_glm)
  mc$family <- "gaussian"
  out <- eval(mc, parent.frame())
  out$call <- call
  out$stan_function <- "stap_lm"
  return(out)
}
