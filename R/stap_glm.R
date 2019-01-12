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
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_aux
#' @template args-adapt_delta
#' @template reference-gelman-hill
#' @template reference-muth
#'
#' @param formula Same as for \code{\link[stats]{glm}}. Note that in-formula transformations will not be passed ot the final design matrix. Covariates that have "scale" in their name are not advised as this text is parsed for in the final model fit.
#' @param family Same as \code{\link[stats]{glm}} for gaussian, binomial, and poisson families.
#' @param subject_data a data.frame that contains data specific to the subject or subjects on whom the outcome is measured. Must contain one column that has the subject_ID  on which to join the distance and time_data
#' @param distance_data a (minimum) three column data.frame that contains (1) an id_key (2) The sap/tap/stap features and (3) the distances between subject with a given id and the built environment feature in column (2), the distance column must be the only column of type "double" and the sap/tap/stap features must be specified in the dataframe exactly as they are in the formula.
#' @param time_data same as distance_data except with time that the subject has been exposed to the built environment feature, instead of distance 
#' @param subject_ID  name of column(s) to join on between subject_data and bef_data
#' @param max_distance the inclusion distance; upper bound for all elements of dists_crs
#' @param max_time inclusion time; upper bound for all elements of times_crs
#' @param model  logical denoting whether or not to return the fixed covariates model frame object in the fitted object
#' @param y In \code{stap_glm}, logical scalar indicating whether to return the response vector. In \code{stan_glm.fit}, a response vector.
#' @param prior_stap prior for spatial-temporal aggregated predictors. Note that prior is set on the standardized latent covariates.
#' @param prior_theta prior for the spatial-temporal aggregated predictors' scale. Can either be a single prior or a prior nested within a list of lists.
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
#'                    subject_data = homog_subject_data[1:100,], # for speed of example only
#'                      distance_data = homog_distance_data,
#'                      family = gaussian(link = 'identity'),
#'                      subject_ID = 'subj_id',
#'                      prior = normal(location = 0, scale = 5, autoscale = FALSE),
#'                      prior_intercept = normal(location = 25, scale = 5, autoscale = FALSE),
#'                      prior_stap = normal(location = 0, scale = 3, autoscale = FALSE),
#'                      prior_theta = log_normal(location = 1, scale = 1),
#'                      prior_aux = cauchy(location = 0,scale = 5),
#'                      max_distance = max(homog_distance_data$Distance),
#'                      chains = 1, iter = 300, # for speed of example only
#'                      refresh = -1, verbose = FALSE) 
#'
stap_glm <- function(formula,
                     family = gaussian(),
                     subject_data = NULL,
                     distance_data = NULL,
                     time_data = NULL,
                     subject_ID = NULL,
                     max_distance = NULL,
                     max_time = NULL,
                     weights,
                     offset = NULL,
                     model = TRUE,
                     y = TRUE,
                     contrasts = NULL,
                     ...,
                     prior = normal(),
                     prior_intercept = normal(),
                     prior_stap = normal(),
                     prior_theta = log_normal(location = 1L, scale = 1L),
                     prior_aux = exponential(), 
                     adapt_delta = NULL){

    stap_data <- extract_stap_data(formula)
    crs_data <- extract_crs_data(stap_data,
                                 subject_data,
                                 distance_data,
                                 time_data,
                                 id_key = subject_ID,
                                 max_distance,
                                 max_time)
    original_formula <- formula
    stapless_formula <- get_stapless_formula(formula)
    family <- validate_family(family)
    validate_glm_formula(stapless_formula)
    subject_data <- validate_data(subject_data, if_missing = environment(stapless_formula))
    call <- match.call(expand.dots = TRUE)
    mf <-  match.call(expand.dots = FALSE)
    mf$formula <- stapless_formula
    m <- match(c("formula", "weights", "offset"),
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
    stapfit <- stap_glm.fit(y = Y, z = Z, 
                            dists_crs = crs_data$d_mat,
                            u_s = crs_data$u_s,
                            times_crs = crs_data$t_mat,
                            u_t = crs_data$u_t,
                            stap_data = stap_data,
                            weights = weights,
                            max_distance = crs_data$max_distance,
                            max_time = crs_data$max_time,
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
                 stap_data = stap_data,
                 subject_data,
                 distance_data,
                 time_data,
                 dists_crs = crs_data$d_mat,
                 times_crs = crs_data$t_mat,
                 u_s = crs_data$u_s,
                 u_t = crs_data$u_t,
                 max_distance = max_distance,
                 offset, weights, z = Z, y = Y,
                 model = mf,  terms = mt, call,
                 contrasts = attr(Z, "contrasts"),
                 stan_function = "stap_glm")
    out <- stapreg(fit)
    out$zlevels <- .getXlevels(mt, mf)
    if (!y)
        out$y <- NULL
    if (!model)
        out$model <- NULL
    return(out)
}

#' @rdname stap_glm
#' @export
stap_lm <- 
  function(formula,
           family = gaussian(),
           subject_data = NULL,
           distance_data = NULL,
           time_data = NULL,
           subject_ID = NULL,
           max_distance = NULL,
           max_time = NULL,
           weights,
           offset = NULL,
           model = TRUE,
           y = TRUE,
           contrasts = NULL,
           ...,
           prior = normal(),
           prior_intercept = normal(),
           prior_stap = normal(),
           prior_theta = log_normal(location = 1L, scale = 1L),
           prior_aux = exponential(),
           adapt_delta = NULL){

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
