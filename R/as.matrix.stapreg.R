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

#' Extract the posterior sample
#' 
#' For models fit using MCMC  the posterior sample ---the post-warmup 
#' draws from the posterior distribution--- can be extracted 
#' from a fitted model object as a matrix, data frame, or array. The 
#' \code{as.matrix} and \code{as.data.frame} methods merge all chains together, 
#' whereas the \code{as.array} method keeps the chains separate.  
#' @method as.matrix stapreg
#' @export
#' @templateVar stapregArg x
#' @template args-stapreg-object
#' @template args-pars
#' @template args-regex-pars
#' @param ... Ignored.
#'   
#' @return A matrix, data.frame, or array, the dimensions of which depend on
#'   \code{pars} and \code{regex_pars}, as well as the model and estimation
#'   algorithm (see the Description section above).
#' 
#' @seealso \code{\link{stapreg-methods}}
#' 
as.matrix.stapreg <- function(x, ..., pars = NULL, regex_pars = NULL) {
  pars <- collect_pars(x, pars, regex_pars)
  user_pars <- !is.null(pars)
  
  mat <- as.matrix(x$stapfit)
    if (!user_pars)
      pars <- exclude_lp_and_ppd(colnames(mat))
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
as.array.stapreg <- function(x, ..., pars = NULL, regex_pars = NULL) {
  pars <- collect_pars(x, pars, regex_pars)

  arr <- as.array(x$stapfit)
  if (identical(arr, numeric(0)))
    STOP_no_draws()
  
  if (!is.null(pars)) {
    check_missing_pars(arr, pars)
  } else {
    pars <- exclude_lp_and_ppd(last_dimnames(arr))
  }
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

