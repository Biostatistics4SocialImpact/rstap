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

#' In-sample or out-of-sample predictive errors
#' 
#' This is a convenience function for computing \eqn{y - y^{rep}}{y - yrep} 
#' (in-sample, for observed \eqn{y}) or \eqn{y - \tilde{y}}{y - ytilde} 
#' (out-of-sample, for new or held-out \eqn{y}). The method for stapreg objects 
#' calls \code{\link{posterior_predict}} internally, whereas the method for
#' objects with class \code{"ppd"} accepts the matrix returned by 
#' \code{posterior_predict} as input and can be used to avoid multiple calls to 
#' \code{posterior_predict}.
#' 
#' @aliases predictive_error
#' @export
#' 
#' @param object Either a fitted model object returned by one of the 
#'   \pkg{rstap} modeling functions (a \link[=stapreg-objects]{stapreg 
#'   object}) or, for the \code{"ppd"} method, a matrix of draws from the 
#'   posterior predictive distribution returned by 
#'   \code{\link{posterior_predict}}.
#' @param newdata,draws,seed,offset,re.form Optional arguments passed to 
#'   \code{\link{posterior_predict}}. For binomial models, please see the
#'   \strong{Note} section below if \code{newdata} will be specified.
#' @template args-dots-ignored
#' 
#' @return A \code{draws} by \code{nrow(newdata)} matrix. If \code{newdata} is 
#'   not specified then it will be \code{draws} by \code{nobs(object)}.
#'   
#' @note The \strong{Note} section in \code{\link{posterior_predict}} about 
#'   \code{newdata} for binomial models also applies for
#'   \code{predictive_error}, with one important difference. For
#'   \code{posterior_predict} if the left-hand side of the model formula is 
#'   \code{cbind(successes, failures)} then the particular values of 
#'   \code{successes} and \code{failures} in \code{newdata} don't matter, only 
#'   that they add to the desired number of trials. \strong{This is not the case
#'   for} \code{predictive_error}. For \code{predictive_error} the particular
#'   value of \code{successes} matters because it is used as \eqn{y} when
#'   computing the error.
#' 
#' @seealso \code{\link[=posterior_predict.stapreg]{posterior_predict}} to draw
#'   from the posterior predictive distribution without computing predictive
#'   errors.
#'   
predictive_error.stapreg <-
  function(object,
           newsubjdata = NULL,
           newdistdata = NULL,
           draws = NULL,
           re.form = NULL,
           seed = NULL,
           offset = NULL,
           ...) {
    if ("y" %in% names(list(...)))
      stop("Argument 'y' should not be specified if 'object' is a stapreg object.")
    
    y <- if (is.null(newdata))
      get_y(object) else eval(formula(object)[[2L]], newdata)
    
    fam <- family(object)$family
    if (is.binomial(fam) && NCOL(y) == 2)
      y <- y[, 1]
    
    ytilde <- posterior_predict(
      object,
      newdata = newdata,
      draws = draws,
      offset = offset,
      seed = seed,
      re.form = re.form
    )
    predictive_error.ppd(ytilde, y = y)
  }

#' @rdname predictive_error.stapreg
#' @export
#' @param y For the \code{"ppd"} method only, a vector of \eqn{y} values the 
#'   same length as the number of columns in the matrix used as \code{object}. 
#'   The method for stapreg objects takes \code{y} directly from the fitted 
#'   model object.
#'   
predictive_error.ppd <- function(object, y, ...) {
  ytilde <- unclass(object)
  rstantools::predictive_error(ytilde, y = y)
}
