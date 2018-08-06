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

#' Bayesian generalized spatial-temporal aggregated predictor(STAP) models via Stan
#'
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Generalized linear modeling with optional prior distributions for the
#' coefficients, intercept, and auxiliary parameters.
#'
#' @export
#' @templateVar armRef (Ch. 3-6)
#' @templateVar pkg stats
#' @templateVar pkgfun glm
#' @templateVar sameargs offset,weights 
#' @templateVar rareargs na.action,contrasts
#' @templateVar fun stap_glm 
#' @templateVar fitfun stan_glm.fit
#' @template return-stapreg-object
#' @template return-stapfit-object
#' @template see-also
#' @template args-formula-data-subset
#' @template args-same-as
#' @template args-same-as-rarely
#' @template args-dots
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_aux
#' @template args-adapt_delta
#' @template reference-gelman-hill
#' @template reference-muth
#' 
#' @param family Same as \code{\link[stats]{glm}}, except negative binomial GLMs
#'   are also possible using the \code{\link{neg_binomial_2}} family object.
#' @param y In \code{stap_glm}, logical scalar indicating whether to
#'   return the response vector. In \code{stan_glm.fit}, a response vector.
#' @param x In \code{stap_glm}, logical scalar indicating whether to
#'   return the design matrix. In \code{stan_glm.fit}, a design matrix.

#' @details The \code{stan_glm} function is similar in syntax to 
#'   \code{\link[stats]{glm}} but rather than performing maximum likelihood 
#'   estimation of generalized linear models, full Bayesian estimation is 
#'   performed (if \code{algorithm} is \code{"sampling"}) via MCMC. The Bayesian
#'   model adds priors (independent by default) on the coefficients of the GLM.
#'   The \code{stan_glm} function calls the workhorse \code{stan_glm.fit}
#'   function, but it is also possible to call the latter directly.
#'   
#'   
#' @seealso The various vignettes for \code{stap_glm} at
#'   \url{https:biostatistics4socialimpact.github.io/rstap/articles}
#' 
stap_glm <- function(formula,
                     family = gaussian(),
                     subject_data,
                     distance_data,
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
                     prior_theta = normal(location = max_distance/2, scale = 10),
                     prior_aux = cauchy(location = 0L, scale = 5L),
                     adapt_delta = NULL){
    crs_data <- extract_stap_components(formula,distance_data,
                                        subject_data, id_key, 
                                        max_distance)
    original_formula <- formula
    formula <- get_stapless_formula(formula)
    family <- validate_family(family)
    validate_glm_formula(formula)
    subject_data <- validate_data(subject_data, if_missing = environment(formula))
    call <- match.call(expand.dots = TRUE)
    mf <-  match.call(expand.dots = FALSE)
    mf$formula <- formula
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
                            dists_crs = crs_data$d_mat, u = crs_data$u,
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
                 distance_data,
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

