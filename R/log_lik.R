# Part of the rstap  package for estimating model parameters
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

#' Pointwise log-likelihood matrix
#'
#' For models fit using MCMC, the \code{log_lik} method returns the
#' \eqn{S} by \eqn{N} pointwise log-likelihood matrix, where \eqn{S} is the size
#' of the posterior sample and \eqn{N} is the number of data points.
#'
#' @aliases log_lik
#' @export
#'
#' @templateVar stapregArg object
#' @template args-stapreg-object
#' @template args-dots-ignored
#' @param newdata An optional data frame of new data (e.g. holdout data) to use
#'   when evaluating the log-likelihood. See the description of \code{newdata}
#'   for \code{\link{posterior_predict}}.
#' @param offset A vector of offsets. Only required if \code{newdata} is
#'   specified and an \code{offset} was specified when fitting the model.
#'
#' @return A \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the posterior  
#'   sample and \eqn{N} is the number of data points. 
#'   
log_lik.stapreg <- function(object, newdata = NULL, offset = NULL, ...) {

  newdata <- validate_newdata(newdata)
  calling_fun <- as.character(sys.call(-1))[1]
  dots <- list(...)
  
  args <- ll_args.stapreg(object, newdata = newdata, offset = offset, 
                          reloo_or_kfold = calling_fun %in% c("kfold", "reloo"), 
                          ...)
  fun <- ll_fun(object, m = m)
  out <- vapply(
      seq_len(args$N),
      FUN = function(i) {
        as.vector(fun(
          data_i = args$data[i, , drop = FALSE],
          draws = args$draws
        ))
      },
      FUN.VALUE = numeric(length = args$S)
    )
  if (is.null(newdata)) colnames(out) <- rownames(model.frame(object, m = m))
  else colnames(out) <- rownames(newdata)
  return(out)
}




# internal ----------------------------------------------------------------

# get log likelihood function for a particular model
# @param x stapreg object
# @return a function
ll_fun <- function(x, m = NULL) {
  validate_stapreg_object(x)
  f <- family(x, m = m)
  if (is.nlmer(x)) 
    return(.ll_nlmer_i)
  fun <- paste0(".ll_", family(x, m = m)$family, "_i")
  get(fun, mode = "function")
}

# get arguments needed for ll_fun
# @param object stapreg object
# @param newdata same as posterior predict
# @param offset vector of offsets (only required if model has offset term and
#   newdata is specified)
# @param reloo_or_kfold logical. TRUE if ll_args is for reloo or kfold
# @param ... For models without group-specific terms (i.e., not stan_[g]lmer), 
#   if reloo_or_kfold is TRUE and 'newdata' is specified then ... is used to 
#   pass 'newx' and 'stanmat' from reloo or kfold (bypassing pp_data). This is a
#   workaround in case there are issues with newdata containing factors with
#   only a single level. 
# @return a named list with elements data, draws, S (posterior sample size) and
#   N = number of observations
ll_args <- function(object, ...) UseMethod("ll_args")
ll_args.stapreg <- function(object, newdata = NULL, offset = NULL, m = NULL, 
                            reloo_or_kfold = FALSE, ...) {
  validate_stapreg_object(object)
  f <- family(object, m = m)
  draws <- nlist(f)
  has_newdata <- !is.null(newdata)
  
  dots <- list(...)
  
  z_betareg <- NULL
  if (has_newdata && reloo_or_kfold && !is.mer(object)) {
    x <- dots$newx
    z_betareg <- dots$newz # NULL except for some stan_betareg models
    if (!is.null(z_betareg)) {
      z_betareg <- as.matrix(z_betareg)
    }
    stanmat <- dots$stanmat
    form <- as.formula(formula(object)) # in case formula is string
    y <- eval(form[[2L]], newdata)
  } else if (has_newdata) {
    ppdat <- pp_data(object, as.data.frame(newdata), offset = offset, m = m)
    pp_eta_dat <- pp_eta(object, ppdat, m = m)
    eta <- pp_eta_dat$eta
    stanmat <- pp_eta_dat$stanmat
    z_betareg <- ppdat$z_betareg
    x <- ppdat$x
    form <- as.formula(formula(object, m = m))
    y <- eval(form[[2L]], newdata)
  } else {
    stanmat <- as.matrix.stapreg(object)
    x <- get_x(object, m = m)
    y <- get_y(object, m = m)
  }
  
    fname <- f$family
    if (is.nlmer(object)) {
      draws <- list(mu = posterior_linpred(object, newdata = newdata),
                    sigma = stanmat[,"sigma"])
      data <- data.frame(y)
      data$offset <- if (has_newdata) offset else object$offset
      if (model_has_weights(object)) {
        data$weights <- object$weights
      }
      data$i_ <- seq_len(nrow(data))  # for nlmer need access to i inside .ll_nlmer_i
      return(nlist(data, draws, S = NROW(draws$mu), N = nrow(data)))
      
    } else if (!is.binomial(fname)) {
      data <- data.frame(y, x)
      if (!is.null(z_betareg)) {
        data <- cbind(data, z_betareg)
      }
    } else {
      if (NCOL(y) == 2L) {
        trials <- rowSums(y)
        y <- y[, 1L]
      } else if (is_clogit(object)) {
        if (has_newdata) strata <- eval(object$call$strata, newdata)
        else strata <- model.frame(object)[,"(weights)"]
        strata <- as.factor(strata)
        successes <- aggregate(y, by = list(strata), FUN = sum)$x
        formals(draws$f$linkinv)$g <- strata
        formals(draws$f$linkinv)$successes <- successes
        trials <- 1L
      } else {
        trials <- 1
        if (is.factor(y)) 
          y <- fac2bin(y)
        stopifnot(all(y %in% c(0, 1)))
      }
      data <- data.frame(y, trials, x)
    }
    nms <- NULL  
    beta_sel <- if (is.null(nms)) seq_len(ncol(x)) else nms$y[[m]]
    draws$beta <- stanmat[, beta_sel, drop = FALSE]
    m_stub <- get_m_stub(m, stub = get_stub(object))
    if (is.gaussian(fname)) 
      draws$sigma <- stanmat[, paste0(m_stub, "sigma")]
    if (is.gamma(fname)) 
      draws$shape <- stanmat[, paste0(m_stub, "shape")]
    if (is.ig(fname)) 
      draws$lambda <- stanmat[, paste0(m_stub, "lambda")]
    if (is.nb(fname)) 
      draws$size <- stanmat[, paste0(m_stub, "reciprocal_dispersion")]
  
  data$offset <- if (has_newdata) offset else object$offset
  if (model_has_weights(object))
    data$weights <- object$weights
    
  if (is.mer(object)) {
    b_sel <- if (is.null(nms)) b_names(colnames(stanmat)) else nms$y_b[[m]]
    b <- stanmat[, b_sel, drop = FALSE]
    if (has_newdata) {
      Z_names <- ppdat$Z_names
      if (is.null(Z_names)) {
        b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
      } else {
        b <- pp_b_ord(b, Z_names)
      }
      if (is.null(ppdat$Zt)) z <- matrix(NA, nrow = nrow(x), ncol = 0)
      else z <- t(ppdat$Zt)
    } else {
      z <- get_z(object, m = m)
    }
    data <- cbind(data, as.matrix(z)[1:NROW(x),, drop = FALSE])
    draws$beta <- cbind(draws$beta, b)
  }
  
    out <- nlist(data, draws, S = NROW(draws$beta), N = nrow(data)) 
  return(out)
}

# log-likelihood function helpers -----------------------------------------
.weighted <- function(val, w) {
  if (is.null(w)) {
    val
  } else {
    val * w
  } 
}

.xdata <- function(data) {
  sel <- c("y", "weights","offset", "trials","strata")
  data[, -which(colnames(data) %in% sel)]
}
.mu <- function(data, draws) {
  eta <- as.vector(linear_predictor(draws$beta, .xdata(data), data$offset))
  draws$f$linkinv(eta)
}

# log-likelihood functions ------------------------------------------------
.ll_gaussian_i <- function(data_i, draws) {
  val <- dnorm(data_i$y, mean = .mu(data_i, draws), sd = draws$sigma, log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_binomial_i <- function(data_i, draws) {
  val <- dbinom(data_i$y, size = data_i$trials, prob = .mu(data_i, draws), log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_poisson_i <- function(data_i, draws) {
  val <- dpois(data_i$y, lambda = .mu(data_i, draws), log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_neg_binomial_2_i <- function(data_i, draws) {
  val <- dnbinom(data_i$y, size = draws$size, mu = .mu(data_i, draws), log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_Gamma_i <- function(data_i, draws) {
  val <- dgamma(data_i$y, shape = draws$shape, 
                rate = draws$shape / .mu(data_i,draws), log = TRUE)
  .weighted(val, data_i$weights)
}
.ll_inverse.gaussian_i <- function(data_i, draws) {
  mu <- .mu(data_i, draws)
  val <- 0.5 * log(draws$lambda / (2 * pi)) - 
    1.5 * log(data_i$y) -
    0.5 * draws$lambda * (data_i$y - mu)^2 / 
    (data_i$y * mu^2)
  .weighted(val, data_i$weights)
}
