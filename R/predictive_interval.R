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

#' Predictive intervals
#' 
#' The \code{predictive_interval} function computes Bayesian predictive intervals. 
#' The method for stapreg objects calls \code{\link{posterior_predict}}
#' internally, whereas the method for objects of class \code{"ppd"} accepts the
#' matrix returned by \code{posterior_predict} as input and can be used to avoid
#' multiple calls to \code{posterior_predict}.
#' 
#' @export
#' @aliases predictive_interval
#' 
#' @param object Either a fitted model object returned by one of the 
#'   \pkg{rstap} modeling functions (a \link[=stapreg-objects]{stapreg 
#'   object}) or, for the \code{"ppd"} method, a matrix of draws from the 
#'   posterior predictive distribution returned by 
#'   \code{\link{posterior_predict}}.
#' @template args-dots-ignored
#' @param draws,fun,offset,re.form,seed Passed to 
#'   \code{\link[=posterior_predict]{posterior_predict}}.
#' @param prob A number \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired
#'   probability mass to include in the intervals. The default is to report
#'   \eqn{90}\% intervals (\code{prob=0.9}) rather than the traditionally used
#'   \eqn{95}\% .
#' @param newsubjdata Optionally, a data frame of the subject-specific data
#'   in which to look for variables with which to predict.
#'   If omitted, the original datasets are used. If \code{newsubjdata}
#'   is provided and any variables were transformed (e.g. rescaled) in the data
#'   used to fit the model, then these variables must also be transformed in
#'   \code{newsubjdata}. This only applies if variables were transformed before
#'   passing the data to one of the modeling functions and \emph{not} if
#'   transformations were specified inside the model formula. Also see the Note
#'   section below for a note about using the \code{newsubjdata} argument with with
#'   binomial models.
#' @param subject_ID same as \code{\link{stap_glm}}
#' @param group_ID same as \code{\link{stap_glmer}}
#' @param newdistdata If newsubjdata is provided a data frame of the subject-distance
#'       must also be given for models with a spatial component
#' @param newtimedata If newsubjdata is provided, a data frame of the subject-time data
#' @return A matrix with two columns and as many rows as are in \code{newsubjdata}. 
#'   If \code{newsubjdata} is not provided then the matrix will have as many rows as
#'   the data used to fit the model. For a given value of \code{prob}, \eqn{p},
#'   the columns correspond to the lower and upper \eqn{100p}\% central interval
#'   limits and have the names \eqn{100\alpha/2}\% and \eqn{100(1 -
#'   \alpha/2)}\%, where \eqn{\alpha = 1-p}. For example, if \code{prob=0.9} is
#'   specified (a \eqn{90}\% interval), then the column names will be
#'   \code{"5\%"} and \code{"95\%"}, respectively.
#'   
#' @seealso \code{\link{predictive_error}}, \code{\link{posterior_predict}}, 
#'   \code{\link{posterior_interval}}
#' @examples
#' if (!exists("example_model")) example(example_model)
#' 
#' predictive_interval(example_model)
#'newdata <- data.frame(subj_ID = 1, measure_ID = 1, centered_income = -1, sex = 0, centered_age = -1) 
#' # newdata
#' predictive_interval(example_model, newsubjdata = newdata, newdistdata = distdata, newtimedata = timedata, subject_ID = "subj_ID", group_ID = "measure_ID")
#' 
predictive_interval.stapreg <-
  function(object,
           prob = 0.9,
           newsubjdata = NULL,
           newdistdata = NULL,
           newtimedata = NULL,
           draws = NULL,
           subject_ID = NULL,
           group_ID = NULL,
           re.form = NULL,
           fun = NULL,
           seed = NULL,
           offset = NULL,
           ...) {
    
    ytilde <- posterior_predict(
      object,
      newsubjdata = newsubjdata,
      newdistdata = newdistdata,
      newtimedata = newtimedata,
      subject_ID = subject_ID,
      group_ID = group_ID,
      draws = draws,
      seed = seed,
      re.form = re.form,
      offset = offset,
      fun = fun
    )
    predictive_interval.ppd(ytilde, prob = prob)
  }

#' @rdname predictive_interval.stapreg
#' @export
predictive_interval.ppd <- function(object, prob = 0.9, ...) {
  ytilde <- unclass(object)
  rstantools::predictive_interval(ytilde, prob = prob)
}
