# Part of the rstap package for estimating model parameters
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

#' Draw from posterior predictive distribution
#'
#' The posterior predictive distribution is the distribution of the outcome
#' implied by the model after using the observed data to update our beliefs
#' about the unknown parameters in the model. Simulating data from the posterior
#' predictive distribution using the observed predictors is useful for checking
#' the fit of the model. Drawing from the posterior predictive distribution at
#' interesting values of the predictors also lets us visualize how a
#' manipulation of a predictor affects (a function of) the outcome(s). With new
#' observations of predictor variables we can use the posterior predictive
#' distribution to generate predicted outcomes.
#'
#' @aliases posterior_predict
#' @export
#'
#' @templateVar stapregArg object
#' @template args-stapreg-object
#' @param newsubjdata Optionally, a data frame of the subject-specific data
#'   in which to look for variables with which to predict.
#'   If omitted, the original datasets are used. If \code{newsubjdata}
#'   is provided and any variables were transformed (e.g. rescaled) in the data
#'   used to fit the model, then these variables must also be transformed in
#'   \code{newsubjdata}. This only applies if variables were transformed before
#'   passing the data to one of the modeling functions and \emph{not} if
#'   transformations were specified inside the model formula. Also see the Note
#'   section below for a note about using the \code{newdata} argument with with
#'   binomial models.
#' @param newdistdata If newsubjdata is provided a data frame of the subject-distance
#'       must also be given for models with a spatial component
#' @param newtimedata If newsubjdata is provided, a data frame of the subject-time data
#'       must also be given for models with a temporal component
#' @param subject_ID  name of column to join on between subject_data and bef_data
#' @param group_ID name of column to join on between \code{subject_data} and bef_data that uniquely identifies the correlated groups (e.g. visits,schools). Currently only one group (e.g. a measurement ID) can be accounted for in a spatial temporal setting. 
#' @param draws An integer indicating the number of draws to return. The default
#'   and maximum number of draws is the size of the posterior sample.
#' @param re.form If \code{object} contains \code{\link[=stap_glmer]{group-level}}
#'   parameters, a formula indicating which group-level parameters to
#'   condition on when making predictions. \code{re.form} is specified in the
#'   same form as for \code{\link[lme4]{predict.merMod}}. The default,
#'   \code{NULL}, indicates that all estimated group-level parameters are
#'   conditioned on. To refrain from conditioning on any group-level parameters,
#'   specify \code{NA} or \code{~0}. The \code{newdata} argument may include new
#'   \emph{levels} of the grouping factors that were specified when the model
#'   was estimated, in which case the resulting posterior predictions
#'   marginalize over the relevant variables.
#' @param fun An optional function to apply to the results. \code{fun} is found
#'   by a call to \code{\link{match.fun}} and so can be specified as a function
#'   object, a string naming a function, etc.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param offset A vector of offsets. Only required if \code{newsubjdata} is
#'   specified and an \code{offset} argument was specified when fitting the
#'   model.
#' @param ... optional arguments to pass to pp_args
#'   
#' @return A \code{draws} by \code{nrow(newdata)} matrix of simulations from the
#'   posterior predictive distribution. Each row of the matrix is a vector of 
#'   predictions generated using a single draw of the model parameters from the 
#'   posterior distribution. The returned matrix will also have class
#'   \code{"ppd"} to indicate it contains draws from the posterior predictive
#'   distribution.
#'
#' @note For binomial models with a number of trials greater than one (i.e., not
#'   Bernoulli models), if \code{newsubjdata} is specified then it must include all
#'   variables needed for computing the number of binomial trials to use for the
#'   predictions. For example if the left-hand side of the model formula is
#'   \code{cbind(successes, failures)} then both \code{successes} and
#'   \code{failures} must be in \code{newdata}. The particular values of
#'   \code{successes} and \code{failures} in \code{newdata} do not matter so
#'   long as their sum is the desired number of trials. If the left-hand side of
#'   the model formula were \code{cbind(successes, trials - successes)} then
#'   both \code{trials} and \code{successes} would need to be in \code{newsubjdata},
#'   probably with \code{successes} set to \code{0} and \code{trials} specifying
#'   the number of trials. 
#'   
#' @seealso \code{\link{pp_check}} for graphical posterior predictive checks.
#'   Examples of posterior predictive checking can also be found in the
#'   \pkg{rstanarm} vignettes and demos.
#'
#' \code{\link{predictive_error}} and \code{\link{predictive_interval}}.
#'
posterior_predict.stapreg <- function(object, 
                                      newsubjdata = NULL, 
                                      newdistdata = NULL,
                                      newtimedata = NULL,
                                      draws = NULL,
                                      subject_ID = NULL,
                                      group_ID = NULL,
                                      re.form = NULL, 
                                      fun = NULL,
                                      seed = NULL,
                                      offset = NULL, ...) {
  if (!is.null(seed))
    set.seed(seed)
  if (!is.null(fun))
    fun <- match.fun(fun)

  dots <- list(...)
  stanmat <- NULL
  
  prediction_data <- validate_predictiondata(newsubjdata, newdistdata, newtimedata)
  if(!is.null(prediction_data$newsubjdata)){
      newsubjdata <- prediction_data$newsubjdata
      if(!is.null(prediction_data$newdistdata))
          newdistdata <- prediction_data$newdistdata
      if(!is.null(prediction_data$newtimedata))
          newtimedata <- prediction_data$newtimedata
   } else{
       newsubjdata <- NULL
       newdistdata <- NULL
       newtimedata <- NULL
  }

  id_key <- if(is.null(group_ID)) subject_ID else c(subject_ID,group_ID)

  pp_data_args <- c(list(object,
                         newsubjdata = newsubjdata,
                         newdistdata = newdistdata,
                         newtimedata = newtimedata,
                         re.form = re.form,
                         offset = offset,
                         id_key = id_key),
                    dots)

  
  dat <- do.call("pp_data", pp_data_args)

  ppargs <- pp_args(object, data = pp_eta(object, dat, draws))

  if (is.binomial(family(object)$family)) {
    ppargs$trials <- pp_binomial_trials(object, newsubjdata)
  }

  ppfun <- pp_fun(object)
  ytilde <- do.call(ppfun, ppargs)
  if (!is.null(newsubjdata) && nrow(newsubjdata) == 1L)
    ytilde <- t(ytilde)
  if(!is.null(fun))
    ytilde <- do.call(fun, list(ytilde))
  
  if (is.null(newsubjdata)) 
      colnames(ytilde) <- rownames(model.frame(object))
  else 
      colnames(ytilde) <- rownames(newsubjdata)  
  
  
  structure(ytilde, class = c("ppd", class(ytilde)))
}

  
# internal ----------------------------------------------------------------

# functions to draw from the various posterior predictive distributions
pp_fun <- function(object) {
  suffix <- family(object)$family

  get(paste0(".pp_", suffix), mode = "function")
}

.pp_gaussian <- function(mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
    rnorm(ncol(mu), mu[s,], sigma[s])
  }))
}
.pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s, ])
  }))
}
.pp_poisson <- function(mu) {
  t(sapply(1:nrow(mu), function(s) {
    rpois(ncol(mu), mu[s, ])
  }))
}
.pp_neg_binomial_2 <- function(mu, size) {
  t(sapply(1:nrow(mu), function(s) {
    rnbinom(ncol(mu), size = size[s], mu = mu[s, ])
  }))
}
.pp_Gamma <- function(mu, shape) {
  t(sapply(1:nrow(mu), function(s) {
    rgamma(ncol(mu), shape = shape[s], rate = shape[s] / mu[s, ])
  }))
}
.rinvGauss <- function(n, mu, lambda) {
  # draw from inverse gaussian distribution
  mu2 <- mu^2
  y <- rnorm(n)^2
  z <- runif(n)
  tmp <- (mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * y^2))
  x <- mu + tmp / (2 * lambda)
  ifelse(z <= (mu / (mu + x)), x, mu2 / x)
}
.pp_inverse.gaussian <- function(mu, lambda) {
  t(sapply(1:nrow(mu), function(s) {
    .rinvGauss(ncol(mu), mu = mu[s,], lambda = lambda[s])
  }))
}


# create list of arguments to pass to the function returned by pp_fun
#
# @param object stapreg 
# @param data output from pp_eta (named list with eta and stanmat)
# @return named list
pp_args <- function(object, data) {
  stanmat <- data$stanmat
  eta <- data$eta
  stopifnot(is.stapreg(object), is.matrix(stanmat))
  inverse_link <- linkinv(object)

  args <- list(mu = inverse_link(eta))
  famname <- family(object)$family
  if (is.gaussian(famname)) {
    args$sigma <- stanmat[, "sigma"]
  } else if (is.gamma(famname)) {
    args$shape <- stanmat[, "shape"]
  } else if (is.ig(famname)) {
    args$lambda <- stanmat[, "lambda"]
  } else if (is.nb(famname)) {
    args$size <- stanmat[, "reciprocal_dispersion"]
  } 
  args
}

# create eta and stanmat (matrix of posterior draws)
#
# @param object A stapreg object
# @param data Output from pp_data()
# @param draws Number of draws
# @return Linear predictor "eta" and matrix of posterior draws "stanmat". 
pp_eta <- function(object, data, draws = NULL) {
  z <- data$z
  x <- data$x 
  S <- posterior_sample_size(object) 
  if (is.null(draws))
    draws <- S
  if (draws > S) {
    err <- paste0("'draws' should be <= posterior sample size (",
                  S, ").")
    stop(err)
  }
  some_draws <- isTRUE(draws < S)
  if (some_draws)
    samp <- sample(S, draws)

  stanmat <- if (is.null(data$w)) 
      as.matrix.stapreg(object) else as.matrix(object$stapfit)
  delta_sel <- seq_len(ncol(z)) 
  delta <- cbind(stanmat[, delta_sel , drop = FALSE],1)
  stap_coefs_nms <- grep("_scale", coef_names(object$stap_data), 
                          invert = T, value = T)

  beta <- stanmat[,stap_coefs_nms,drop =F]

  if (some_draws){
    delta <- delta[samp, , drop = FALSE]
    beta <- beta[samp,,drop=F] 
    x <- x[samp,,,drop=F]
  }
  num_draws <- if(some_draws) draws else S
  if(object$stap_data$Q == 1)
      x_beta <- sapply(1:length(num_draws), function(i) x[i,,] * beta[i,])
  else
      x_beta <- sapply(1:length(num_draws), function(i) x[i,,] %*% t(beta[i,,drop=F]))
   

  eta <- linear_predictor(delta, cbind(z,x_beta), data$offset)
  if (!is.null(data$w)) {
    b_sel <- grepl("^b\\[", colnames(stanmat)) 
    b <- stanmat[, b_sel, drop = FALSE]
    if (some_draws)
      b <- b[samp, , drop = FALSE]
    if (is.null(data$Z_names)) {
      b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
    } else {
      b <- pp_b_ord(b, data$Z_names)
    }
    eta <- eta + as.matrix(b %*% data$Zt)
  }
  
  out <- nlist(eta, stanmat)
  return(out)
}

pp_b_ord <- function(b, Z_names) {
  b_ord <- function(x) {
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    m <- grep(paste0("b[", sub(" (.*):.*$", " \\1:_NEW_\\1", x), "]"),
              colnames(b), fixed = TRUE)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    x <- strsplit(x, split = ":", fixed = TRUE)[[1]]
    stem <- strsplit(x[[1]], split = " ", fixed = TRUE)[[1]]
    x <- paste(x[1], x[2], paste0("_NEW_", stem[2]), x[2], sep = ":")
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    x <- paste(paste(stem[1], stem[2]), paste0("_NEW_", stem[2]), sep = ":")
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    stop("no matches bug")
  }
  ord <- sapply(Z_names, FUN = b_ord)
  b[, ord, drop = FALSE]
}

# Number of trials for binomial models
pp_binomial_trials <- function(object, newsubjdata = NULL) {
  
  y <- get_y(object) 
  is_bernoulli <- NCOL(y) == 1L
  
  if (is_bernoulli) {
    trials <- if (is.null(newsubjdata)) 
      rep(1, NROW(y)) else rep(1, NROW(newsubjdata))
  } else {
    trials <- if (is.null(newsubjdata)) 
      rowSums(y) else rowSums(eval(formula(object)[[2L]], newsubjdata))
  }
  return(trials)
}
