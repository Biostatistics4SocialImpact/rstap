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
#' @param newdata Optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the model matrix is used. If \code{newdata}
#'   is provided and any variables were transformed (e.g. rescaled) in the data
#'   used to fit the model, then these variables must also be transformed in
#'   \code{newdata}. This only applies if variables were transformed before
#'   passing the data to one of the modeling functions and \emph{not} if
#'   transformations were specified inside the model formula. Also see the Note
#'   section below for a note about using the \code{newdata} argument with with
#'   binomial models.
#' @param draws An integer indicating the number of draws to return. The default
#'   and maximum number of draws is the size of the posterior sample.
#' @param re.form If \code{object} contains \code{\link[=stan_glmer]{group-level}}
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
#' @param offset A vector of offsets. Only required if \code{newdata} is
#'   specified and an \code{offset} argument was specified when fitting the
#'   model.
#'   
#' @return A \code{draws} by \code{nrow(newdata)} matrix of simulations from the
#'   posterior predictive distribution. Each row of the matrix is a vector of 
#'   predictions generated using a single draw of the model parameters from the 
#'   posterior distribution. The returned matrix will also have class
#'   \code{"ppd"} to indicate it contains draws from the posterior predictive
#'   distribution.
#'
#' @note For binomial models with a number of trials greater than one (i.e., not
#'   Bernoulli models), if \code{newdata} is specified then it must include all
#'   variables needed for computing the number of binomial trials to use for the
#'   predictions. For example if the left-hand side of the model formula is
#'   \code{cbind(successes, failures)} then both \code{successes} and
#'   \code{failures} must be in \code{newdata}. The particular values of
#'   \code{successes} and \code{failures} in \code{newdata} do not matter so
#'   long as their sum is the desired number of trials. If the left-hand side of
#'   the model formula were \code{cbind(successes, trials - successes)} then
#'   both \code{trials} and \code{successes} would need to be in \code{newdata},
#'   probably with \code{successes} set to \code{0} and \code{trials} specifying
#'   the number of trials. See the Examples section below and the
#'   \emph{How to Use the rstaparm Package} for examples.
#' @note For models estimated with \code{\link{stan_clogit}}, the number of 
#'   successes per stratum is ostensibly fixed by the research design. Thus, when
#'   doing posterior prediction with new data, the \code{data.frame} passed to
#'   the \code{newdata} argument must contain an outcome variable and a stratifying
#'   factor, both with the same name as in the original \code{data.frame}. Then, the 
#'   posterior predictions will condition on this outcome in the new data.
#'   
#' @seealso \code{\link{pp_check}} for graphical posterior predictive checks.
#'   Examples of posterior predictive checking can also be found in the
#'   \pkg{rstanarm} vignettes and demos.
#'
#' \code{\link{predictive_error}} and \code{\link{predictive_interval}}.
#'
posterior_predict.stapreg <- function(object, newsubjdata = NULL,newdistancedata = NULL,
                                      newtimedata = NULL, draws = NULL,
                                      re.form = NULL, fun = NULL, seed = NULL,
                                      offset = NULL, ...) {
  if (!is.null(seed))
    set.seed(seed)
  if (!is.null(fun))
    fun <- match.fun(fun)

  dots <- list(...)
  m <- NULL
  stanmat <- NULL
  
  newdata <- validate_newdata(newdata)
  pp_data_args <- c(list(object,
                         newdata = newdata,
                         re.form = re.form,
                         offset = offset),
                    dots)
  dat <- do.call("pp_data", pp_data_args)
  ppargs <- pp_args(object, data = pp_eta(object, dat, draws), m = m)

  if (is.binomial(family(object, m = m)$family)) {
    ppargs$trials <- pp_binomial_trials(object, newdata, m = m)
  }

  ppfun <- pp_fun(object, m = m)
  ytilde <- do.call(ppfun, ppargs)
  if (!is.null(newdata) && nrow(newdata) == 1L)
    ytilde <- t(ytilde)
  if (!is.null(fun))
    ytilde <- do.call(fun, list(ytilde))
  
  if (is.null(newdata)) colnames(ytilde) <- rownames(model.frame(object, m = m))
  else colnames(ytilde) <- rownames(newdata)  
  
  # if function is called from posterior_traj then add mu as attribute
  fn <- tryCatch(sys.call(-3)[[1]], error = function(e) NULL)
  if (!is.null(fn) && grepl("posterior_traj", deparse(fn), fixed = TRUE))
    return(structure(ytilde, mu = ppargs$mu, class = c("ppd", class(ytilde))))
  
  structure(ytilde, class = c("ppd", class(ytilde)))
}

  
# internal ----------------------------------------------------------------

# functions to draw from the various posterior predictive distributions
pp_fun <- function(object, m = NULL) {
  suffix <- if (is_polr(object)) "polr" else 
            if (is_clogit(object)) "clogit" else 
            family(object, m = m)$family

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
.pp_clogit <- function(mu, strata) {
  t(sapply(1:nrow(mu), function(s) {
    unlist(by(mu[s,], INDICES = list(strata), FUN = rmultinom, n = 1, size = 1))
  }))
}
.pp_beta <- function(mu, phi) {
  t(sapply(1:nrow(mu), function(s) {
    rbeta(ncol(mu), mu[s,] * phi[s], (1 - mu[s, ]) * phi[s])
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
.pp_polr <- function(eta, zeta, linkinv, alpha = NULL) {
  n <- ncol(eta)
  q <- ncol(zeta)
  if (!is.null(alpha)) {
    pr <- linkinv(eta)^alpha
    if (NROW(eta) == 1) {
      pr <- matrix(pr, nrow = 1)
    }
    t(sapply(1:NROW(eta), FUN = function(s) {
      rbinom(NCOL(eta), size = 1, prob = pr[s, ])
    }))
  } else {
    t(sapply(1:NROW(eta), FUN = function(s) {
      tmp <- matrix(zeta[s, ], n, q, byrow = TRUE) - eta[s, ]
      cumpr <- matrix(linkinv(tmp), ncol = q)
      fitted <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
      apply(fitted, 1, function(p) which(rmultinom(1, 1, p) == 1))
    }))
  }
}


# create list of arguments to pass to the function returned by pp_fun
#
# @param object stapreg 
# @param data output from pp_eta (named list with eta and stanmat)
# @param m optional integer specifying the submodel for stanmvreg objects
# @return named list
pp_args <- function(object, data, m = NULL) {
  stanmat <- data$stanmat
  eta <- data$eta
  stopifnot(is.stapreg(object), is.matrix(stanmat))
  if (is.stanmvreg(object) && is.null(m)) STOP_arg_required_for_stanmvreg(m)
  inverse_link <- linkinv(object, m = m)
  if (is.nlmer(object)) inverse_link <- function(x) return(x)

  if (is_polr(object)) {
    zeta <- stanmat[, grep("|", colnames(stanmat), value = TRUE, fixed = TRUE)]
    args <- nlist(eta, zeta, linkinv = inverse_link)
    if ("alpha" %in% colnames(stanmat)) # scobit
      args$alpha <- stanmat[, "alpha"]
    return(args)
  }
  else if (is_clogit(object)) {
    return(list(mu = inverse_link(eta)))
  }

  args <- list(mu = inverse_link(eta))
  famname <- family(object, m = m)$family
  m_stub <- get_m_stub(m, stub = get_stub(object))
  if (is.gaussian(famname)) {
    args$sigma <- stanmat[, paste0(m_stub, "sigma")]
  } else if (is.gamma(famname)) {
    args$shape <- stanmat[, paste0(m_stub, "shape")]
  } else if (is.ig(famname)) {
    args$lambda <- stanmat[, paste0(m_stub, "lambda")]
  } else if (is.nb(famname)) {
    args$size <- stanmat[, paste0(m_stub, "reciprocal_dispersion")]
  } else if (is.beta(famname)) {
    args$phi <- data$phi
    if (is.null(args$phi)) {
      args$phi <- linkinv(object$family_phi)(data$phi_linpred)
    }
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
  x <- data$x
  S <- if (is.null(stanmat)) posterior_sample_size(object) else nrow(stanmat)
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
  if (is.null(stanmat)) {
    stanmat <- if (is.null(data$Zt)) 
      as.matrix.stapreg(object) else as.matrix(object$stanfit)
  }
  nms <-  NULL  
  beta_sel <- if (is.null(nms)) seq_len(ncol(x)) else nms$y[[m]]
  beta <- stanmat[, beta_sel, drop = FALSE]
  if (some_draws)
    beta <- beta[samp, , drop = FALSE]
  eta <- linear_predictor(beta, x, data$offset)
  if (!is.null(data$Zt)) {
    b_sel <- if (is.null(nms)) grepl("^b\\[", colnames(stanmat)) else nms$y_b[[m]]
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
  if (is.nlmer(object)) {
    if (is.null(data$arg1)) eta <- linkinv(object)(eta)
    else eta <- linkinv(object)(eta, data$arg1, data$arg2)
    eta <- t(eta)
  }
  
  out <- nlist(eta, stanmat)
  
  if (inherits(object, "betareg")) {
    z_vars <- colnames(stanmat)[grepl("(phi)", colnames(stanmat))]
    omega <- stanmat[, z_vars]
    if (length(z_vars) == 1 && z_vars == "(phi)") {
      out$phi <- stanmat[, "(phi)"] 
    } else {
      out$phi_linpred <- linear_predictor(as.matrix(omega), as.matrix(data$z_betareg), data$offset)
    }
  }
  
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
pp_binomial_trials <- function(object, newdata = NULL, m = NULL) {
  
  y <- get_y(object, m) 
  is_bernoulli <- NCOL(y) == 1L
  
  if (is_bernoulli) {
    trials <- if (is.null(newdata)) 
      rep(1, NROW(y)) else rep(1, NROW(newdata))
  } else {
    trials <- if (is.null(newdata)) 
      rowSums(y) else rowSums(eval(formula(object, m = m)[[2L]], newdata))
  }
  return(trials)
}
