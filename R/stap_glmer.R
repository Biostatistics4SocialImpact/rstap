# Part of the rstap package for estimating model parameters
# Copyright (c) 2018
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# # This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the # GNU General Public License for more details.
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
#' @templateVar fun stap_glmer, stap_lmer 
#' @templateVar pkg lme4
#' @templateVar pkgfun glmer
#' @template return-stapreg-object
#' @template see-also
#' @template args-prior_intercept
#' @template args-priors
#' @template args-prior_aux
#' @template args-prior_covariance
#' @template args-adapt_delta
#' @template reference-gelman-hill
#' @template reference-muth
#' 
#' @param formula Same as for \code{\link[lme4]{glmer}}. Note that in-formula
#'  transformations will not be passed to the final design matrix.Covariates that have "scale" in their name are not advised as this text is parsed for in the final model fit.
#' @param subject_data a data.frame that contains data specific to the subject or subjects on whom the outcome is measured. Must contain one column that has the subject_ID  on which to join the distance and time_data
#' @param distance_data a (minimum) three column data.frame that contains (1) an id_key (2) The sap/tap/stap features and (3) the distances between subject with a given id and the built environment feature in column (2), the distance column must be the only column of type "double" and the sap/tap/stap features must be specified in the dataframe exactly as they are in the formula.
#' @param time_data same as distance_data except with time that the subject has been exposed to the built environment feature, instead of distance 
#' @param subject_ID  name of column to join on between subject_data and bef_data
#' @param group_ID  name of column to join on between \code{subject_data} and bef_data that uniquely identifies the groups 
#'@param max_distance the upper bound on any and all distances included in the model 
#'@param max_time the upper bound on any and all times included in the model
#'@param family Same as for \code{\link[lme4]{glmer}} except limited to gaussian, binomial and poisson 
#'@param prior_theta,prior_stap priors for the spatial scale and spatial effect parameters, respectively
#'@param weights,offset Same as \code{\link[stats]{glm}}.
#'@param contrasts Same as \code{\link[stats]{glm}}, but rarely 
#'   specified.
#' @param ... For \code{stap_glmer}, further arguments passed to 
#'   \code{\link[rstan]{sampling}} (e.g. \code{iter}, \code{chains}, 
#'   \code{cores}, etc.). For \code{stap_lmer} 
#'    \code{...} should also contain all relevant arguments
#'   to pass to \code{stap_glmer} (except \code{family}).
#'
#' @details The \code{stap_glmer} function is similar in syntax to 
#'   \code{\link[lme4]{glmer}} but rather than performing (restricted) maximum 
#'   likelihood estimation of generalized linear models, Bayesian estimation is 
#'   performed via MCMC. The Bayesian model adds priors on the 
#'   regression coefficients (in the same way as \code{\link{stap_glm}}) and
#'   priors on the terms of a decomposition of the covariance matrices of the
#'   group-specific parameters. See \code{\link{priors}} for more information
#'   about the priors.
#'   
#'   The \code{stap_lmer} function is equivalent to \code{stap_glmer} with 
#'   \code{family = gaussian(link = "identity")}. 
#'   
#' @seealso The Longituinal \href{https://biostatistics4socialimpact.github.io/rstap/articles/longitudinal-I.html}{Vignette} for \code{stap_glmer}.
#'    
#' @importFrom lme4 glFormula
#' @importFrom Matrix Matrix t
#'@examples
#'\dontrun{
#' ## subset to only include id, class name and distance variables
#' distdata <- homog_longitudinal_bef_data[,c("subj_ID","measure_ID","class","dist")]
#' timedata <- homog_longitudinal_bef_data[,c("subj_ID","measure_ID","class","time")]
#' ## distance or time column must be numeric
#' timedata$time <- as.numeric(timedata$time) 
#' fit <- stap_glmer(y_bern ~ centered_income +  sex + centered_age + stap(Coffee_Shop) + (1|subj_ID),
#'                   family = binomial(link='logit'),
#'                   subject_data = homog_longitudinal_subject_data,
#'                   distance_data = distdata,
#'                   time_data = timedata,
#'                   subject_ID = 'subj_ID',
#'                   group_ID = 'measure_ID',
#'                   prior_intercept = normal(location = 25, scale = 4, autoscale = F),
#'                   prior = normal(location = 0, scale = 4, autoscale=F),
#'                   prior_stap = normal(location = 0, scale = 4),
#'                   prior_theta = list(Coffee_Shop = list(spatial = log_normal(location = 1,
#'                                                                              scale = 1),
#'                                                          temporal = log_normal(location = 1,
#'                                                                                scale = 1))),
#'                   max_distance = 3, max_time = 50,
#'                   chains = 4, refresh = -1, verbose = FALSE, 
#'                   iter = 1E3, cores = 1)
#'}
stap_glmer <- 
  function(formula,
           family = gaussian(),
           subject_data = NULL,
           distance_data = NULL,
           time_data = NULL,
           subject_ID = NULL,
           group_ID = NULL,
           max_distance = NULL,
           max_time = NULL,
           weights,
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
  original_formula <- formula
  stapless_formula <- get_stapless_formula(formula)
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  mc$formula <- stapless_formula
  data <- validate_data(subject_data , if_missing = environment(formula))
  family <- validate_family(family)
  mc[[1]] <- quote(lme4::glFormula)
  mc$control <- make_glmerControl()
  mc$data <- subject_data
  mc$prior <- mc$prior_intercept <- mc$prior_covariance <- mc$prior_aux <-
    mc$prior_PD <- mc$algorithm <- mc$scale <- mc$concentration <- mc$shape <-
    mc$adapt_delta <- mc$... <- mc$subject_data <- mc$distance_data <- mc$time_data <- 
      mc$subject_ID <- mc$group_ID <- mc$max_distance <- 
      mc$max_time <- mc$prior_stap <- mc$prior_theta <- NULL
  glmod <- eval(mc, parent.frame())
  Z <- glmod$X
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
  stap_data <- extract_stap_data(formula)
  crs_data <- extract_crs_data(stap_data,
                               subject_data,
                               distance_data,
                               time_data,
                               id_key = c(subject_ID,group_ID),
                               max_distance,
                               max_time)
  stapfit <- stap_glm.fit(y = y,z = Z, 
                          dists_crs = crs_data$d_mat,
                          u_s = crs_data$u_s,
                          times_crs = crs_data$t_mat,
                          u_t = crs_data$u_t,
                          stap_data = stap_data,
                          max_distance = crs_data$max_distance,
                          max_time =crs_data$max_time,
                          weights = weights,
                          offset = offset, family = family,
                          prior = prior,
                          prior_intercept = prior_intercept,
                          prior_stap = prior_stap,
                          prior_aux = prior_aux, 
                          prior_theta = prior_theta,
                          adapt_delta = adapt_delta,
                          group = group,  ...)
  sel <- apply(Z, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  Z <- Z[ , !sel, drop = FALSE]
  W <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(W) <- b_names(names(stapfit), value = TRUE)
  
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
               max_time = crs_data$max_time,
               offset, weights,
               z = Z, w = W, 
               y = y, data,
               call, terms = NULL,
               model = NULL,
               contrasts, glmod, 
               stan_function = "stap_glmer")
  out <- stapreg(fit)
  class(out) <- c(class(out), "lmerMod")
  
  return(out)
}


#' @rdname stap_glmer
#' @export
stap_lmer <- 
  function(formula,
           subject_data = NULL,
           distance_data = NULL,
           time_data = NULL,
           subject_ID = NULL,
           group_ID = NULL,
           max_distance = NULL,
           max_time = NULL,
           weights,
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
  if ("family" %in% names(list(...))) {
    stop(
      "'family' should not be specified. ", 
      "To specify a family use stap_glmer instead of stap_lmer."
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
  out$stan_function <- "stap_lmer"
  return(out)
}
