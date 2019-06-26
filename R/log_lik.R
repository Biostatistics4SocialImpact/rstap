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
#' @param newsubjdata Optionally, a data frame of the subject-specific data
#'   in which to look for variables with which to predict.
#'   If omitted, the original datasets are used. If \code{newsubjdata}
#'   is provided and any variables were transformed (e.g. rescaled) in the data
#'   used to fit the model, then these variables must also be transformed in
#'   \code{newsubjdata}.  Also see the Note
#'   section below for a note about using the \code{newsubjdata} argument with with
#'   binomial models.
#' @param newdistdata If newsubjdata is provided a data frame of the subject-distance
#'       must also be given for models with a spatial component - can be the same as original distance_dataframe
#' @param newtimedata If newsubjdata is provided, a data frame of the subject-time data
#'       must also be given for models with a temporal component
#' @param offset A vector of offsets. Only required if \code{newsubjdata} is
#'   specified and an \code{offset} was specified when fitting the model.
#'
#' @return A \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the posterior  
#'   sample and \eqn{N} is the number of data points. 
#'   
log_lik.stapreg <- function(object, 
                            newsubjdata = NULL,
                            newdistdata = NULL,
                            newtimedata = NULL,
                            offset = NULL, ...) {

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
  calling_fun <- as.character(sys.call(-1))[1]
  dots <- list(...)
  
  args <- ll_args.stapreg(object,
                          newsubjdata = newsubjdata,
                          newdistdata = newdistdata,
                          newtimedata = newtimedata,
                          offset = offset, 
                          ...)
  fun <- ll_fun(object)
  out <- vapply(
      seq_len(args$N),
      FUN = function(i) {
        as.vector(fun(
          data_i = args$data[i, ,drop=F],
          draws = args$draws))
      },
      FUN.VALUE = numeric(length = args$S)
    )
  if (is.null(newsubjdata)) colnames(out) <- rownames(model.frame(object))
  else colnames(out) <- rownames(newsubjdata)
  return(out)
}




# internal ----------------------------------------------------------------

# get log likelihood function for a particular model
# @param x stapreg object
# @return a function
ll_fun <- function(x) {
  validate_stapreg_object(x)
  f <- family(x)
  fun <- paste0(".ll_", family(x)$family, "_i")
  get(fun, mode = "function")
}

# get arguments needed for ll_fun
# @param object stapreg object
# @param newsubjdata same as posterior predict
# @param newdistdata same as posterior predict
# @param newtimdata same as posterior predict
# @param offset vector of offsets (only required if model has offset term and
#   newsubjdata is specified)
# @param reloo_or_kfold logical. TRUE if ll_args is for reloo or kfold
# @return a named list with elements data, draws, S (posterior sample size) and
#   N = number of observations
ll_args <- function(object, ...) UseMethod("ll_args")
ll_args.stapreg <- function(object,
                            newsubjdata = NULL,
                            newdistdata = NULL,
                            newtimedata = NULL,
                            offset = NULL, 
                            ...) {

  validate_stapreg_object(object)
  f <- family(object)
  draws <- nlist(f)
  has_newdata <- !is.null(newsubjdata)
  
  dots <- list(...)
  
  if (has_newdata) {
    ppdat <- pp_data(object,
                     newsubjdata = newsubjdata,
                     newdistdata = newdistdata,
                     newtimedata = newtimedata,
                      offset = offset)
    pp_eta_dat <- pp_eta(object, ppdat)
    eta <- pp_eta_dat$eta
    stanmat <- pp_eta_dat$stanmat
    x <- ppdat$x
    z <- ppdat$z
    y <- f$linkinv(eta)
  } else {
    stanmat <- as.matrix.stapreg(object)
    z <- get_z(object)
    x <- get_x(object)
    y <- get_y(object)
    beta <- beta_names(object$stap_data)
    if(any_dnd(object$stap_data)){
      subj_mat <- as.matrix(Matrix::fac2sparse(as.factor(object$glmod$reTrms$flist[[1]])))
      subj_n_diag <- solve(diag(as.vector(table(object$glmod$reTrms$flist[[1]]) )))
      X_bar <- array(dim=dim(x))
      for(i in 1:object$stap_data$Q)
        X_bar[,,i] <- t(t(subj_mat) %*% (subj_n_diag %*% (subj_mat %*% t(x[,,i]) )))
      X_tilde <- x - X_bar
    }
    beta_samps <- stanmat[, beta, drop=F]
    stap_exp <- matrix(NA,dim(x)[2],dim(x)[1]) # num_obs x num_samps
    if(any_dnd(object$stap_data)){
      for(subj_ix in 1:dim(X_tilde)[2]){
        X_ <- as.matrix(X_tilde[,subj_ix,])
        if(any_bar(object$stap_data))
          X_ <- cbind(as.matrix(X_tilde[,subj_ix,]),as.matrix(X_bar[,subj_ix,]) )
        stap_exp[subj_ix,] <- rowSums(X_ * beta_samps)
      }
    }
    else{
      for(subj_ix in 1:dim(x)[2])
          stap_exp[subj_ix,] <- rowSums(as.matrix(x[,subj_ix,]) * beta_samps)
    }
    colnames(stap_exp) <- paste0("stap_exp_",1:ncol(stap_exp))
  }
  
    fname <- f$family
    if (!is.binomial(fname)) {
      data <- data.frame(y,z,stap_exp)
    } else {
      if (NCOL(y) == 2L) {
        trials <- rowSums(y)
        y <- y[, 1L]
      }  else {
        trials <- 1
        if (is.factor(y)) 
          y <- fac2bin(y)
        stopifnot(all(y %in% c(0, 1)))
      }
      data <- data.frame(y, trials, z, stap_exp)
    }
    delta <- colnames(z) 
    draws$delta <- stanmat[, delta, drop = FALSE]
    if (is.gaussian(fname)) 
      draws$sigma <- stanmat[,  "sigma"]
    if (is.gamma(fname)) 
      draws$shape <- stanmat[, "shape"]
    if (is.ig(fname)) 
      draws$lambda <- stanmat[, "lambda"]
    if (is.nb(fname)) 
      draws$size <- stanmat[, "reciprocal_dispersion"]
  
  data$offset <- if (has_newdata) offset else object$offset
  if (model_has_weights(object))
    data$weights <- object$weights
  
  if (is.mer(object)) {
    b_sel <-  b_names(colnames(stanmat)) 
    b <- stanmat[, b_sel, drop = FALSE]
    if (has_newdata) {
      Z_names <- ppdat$Z_names
      if (is.null(Z_names)) {
        b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
      } else {
        b <- pp_b_ord(b, Z_names)
      }
      if (is.null(ppdat$Zt)) z <- matrix(NA, nrow = nrow(x), ncol = 0)
      else w <- t(ppdat$Zt)
    } else {
      w <- get_w(object)
    }
    if(!is.binomial(fname))
      data <- data.frame(y,z,as.matrix(w),stap_exp)
    else
      data <- data.frame(y,trials,z,as.matrix(w),stap_exp)
    draws$delta <- cbind(draws$delta, b)
  }
  
    out <- nlist(data, draws, S = NROW(draws$delta), N = nrow(data)) 
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
  stap_sel <-  grep("stap_exp",colnames(data))
  data[, -c(which(colnames(data) %in% sel),stap_sel)]
}

.nxdata <- function(data) {
  sel <- c("y", "weights","offset", "trials","strata")
  data[, -which(colnames(data) %in% sel)]
}

.stap_samples <- function(data){
    stap_sel <- grep("stap_exp",colnames(data))
    data[,stap_sel]
}

.mu <- function(data, draws) {
  eta <- as.vector(linear_predictor(draws$delta,.xdata(data), data$offset))
  stap_exp <- as.vector(as.numeric(.stap_samples(data)))
  draws$f$linkinv(eta + stap_exp)
}

# log-likelihood functions ------------------------------------------------
.ll_gaussian_i <- function(data_i, draws,log_switch = T) {
  val <- dnorm(data_i$y, mean = .mu(data_i, draws), sd = draws$sigma, log = log_switch)
  .weighted(val, data_i$weights)
}
.ll_binomial_i <- function(data_i, draws, log_switch = T) {
  val <- dbinom(data_i$y, size = data_i$trials, prob = .mu(data_i, draws) , log = log_switch)
  .weighted(val, data_i$weights)
}
.ll_poisson_i <- function(data_i, draws, log_switch = T) {
  val <- dpois(data_i$y, lambda = .mu(data_i, draws) , log = log_switch)
  .weighted(val, data_i$weights)
}
.ll_neg_binomial_2_i <- function(data_i, draws, log_switch = T) {
  val <- dnbinom(data_i$y, size = draws$size, mu = .mu(data_i, draws), log = log_switch)
  .weighted(val, data_i$weights)
}
.ll_Gamma_i <- function(data_i, draws, log_switch = T) {
  val <- dgamma(data_i$y, shape = draws$shape, 
                rate = draws$shape / ( .mu(data_i,draws)), log = log_switch)
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
