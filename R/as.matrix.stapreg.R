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

#' Extract the posterior sample via matrix 
#' 
#' The posterior sample ---the post-warmup 
#' draws from the posterior distribution--- can be extracted 
#' from a fitted model object as a matrix, data frame, or array. The 
#' \code{as.matrix} and \code{as.data.frame} methods merge all chains together, 
#' whereas the \code{as.array} method keeps the chains separate.  
#' @name as.matrix.stapreg
#' @method as.matrix stapreg
#' @export
#' @templateVar stapregArg x
#' @template args-stapreg-object
#' @template args-pars
#' @template args-regex-pars
#' @param include_X logical to include latent exposure covariate in samples
#' @param ... Ignored.
#'   
#' @return A matrix, data.frame, or array, the dimensions of which depend on
#'   \code{pars} and \code{regex_pars}, as well as the model and estimation
#'   algorithm (see the Description section above).
#' 
#' @seealso \code{\link{stapreg-methods}}
#' 
#' @examples
#' \donttest{
#' if (!exists("example_model")) example(example_model)
#' # Extract posterior sample after MCMC
#' draws <- as.matrix(example_model)
#' print(dim(draws))
#'
#' # For example, we can see that the median of the draws for the intercept 
#' # is the same as the point estimate rstanarm uses
#' print(median(draws[, "(Intercept)"]))
#' print(example_model$coefficients[["(Intercept)"]])
#' 
#' # The as.array method keeps the chains separate
#' draws_array <- as.array(example_model)
#' print(dim(draws_array)) # iterations x chains x parameters
#'}
as.matrix.stapreg <- function(x, ..., pars = NULL, regex_pars = NULL,include_X = FALSE) {
  pars <- collect_pars(x, pars, regex_pars)
  user_pars <- !is.null(pars)
  
  mat <- as.matrix(x$stapfit)
  if (!user_pars)
      pars <- exclude_lp_and_ppd(colnames(mat))
  if(!include_X)
    pars <- pars[grep("X_theta_*",pars,invert=T)]
  if (user_pars)
    check_missing_pars(mat, pars)
  
  mat <- mat[, pars, drop = FALSE]
  if (!is.mer(x))
    return(mat)
  unpad_reTrms(mat)
}

#' @rdname as.matrix.stapreg
#' @method as.array stapreg
#' @export
as.array.stapreg <- function(x, ..., pars = NULL, regex_pars = NULL, include_X = FALSE) {
  pars <- collect_pars(x, pars, regex_pars)

  arr <- as.array(x$stapfit)
  if (identical(arr, numeric(0)))
    STOP_no_draws()
  
  if (!is.null(pars)) {
    check_missing_pars(arr, pars)
  } else {
    pars <- exclude_lp_and_ppd(last_dimnames(arr))
  }
  if(!include_X)
    pars <- pars[grep("X_theta_*",pars,invert=T)]

  arr <- arr[, , pars, drop = FALSE]
  
  if (!is.mer(x))
    return(arr)
  unpad_reTrms(arr)
}


#' @rdname as.matrix.stapreg
#' @method as.data.frame stapreg
#' @export
as.data.frame.stapreg <- function(x, ..., pars = NULL, regex_pars = NULL) {
  mat <- as.matrix.stapreg(x, pars = pars, regex_pars = regex_pars, ...)
  as.data.frame(mat)
}



# internal ----------------------------------------------------------------
STOP_no_draws <- function() stop("No draws found.", call. = FALSE)

check_missing_pars <- function(x, pars) {
  notfound <- which(!pars %in% last_dimnames(x))
  if (length(notfound)) 
    stop(
      "No parameter(s) ", 
      paste(pars[notfound], collapse = ", "), 
      call. = FALSE
    )
}

exclude_lp_and_ppd <- function(pars) {
  grep(
    pattern = "mean_PPD|log-posterior|NA", 
    x = pars, 
    invert = TRUE, 
    value = TRUE
  )
}

