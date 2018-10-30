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

#' Fitted model objects 
#' 
#' The \pkg{rstap} model-fitting functions return an object of class 
#' \code{'stapreg'}, which is a list containing at a minimum the components listed 
#' below. Each \code{stapreg} object will also have additional classes (e.g. 'glm')
#' and several additional components depending on the model and estimation 
#' algorithm. \cr
#' \cr
#' 
#' @name stapreg-objects 
#'   
#' @section Elements for \code{stapreg} objects:   
#' \describe{
#'   \item{\code{coefficients}}{
#'   Point estimates, as described in \code{\link{print.stapreg}}.
#'   }
#'   \item{\code{ses}}{
#'   Standard errors based on \code{\link[stats]{mad}}, as described in
#'   \code{\link{print.stapreg}}.
#'   }
#'   \item{\code{residuals}}{
#'   Residuals of type \code{'response'}.
#'   }
#'   \item{\code{fitted.values}}{
#'   Fitted mean values. For GLMs the linear predictors are transformed by the
#'   inverse link function.
#'   }
#'   \item{\code{linear.predictors}}{
#'   Linear fit on the link scale. For linear models this is the same as
#'   \code{fitted.values}.
#'   }
#'   \item{\code{covmat}}{
#'   Variance-covariance matrix for the coefficients based on draws from the
#'   posterior distribution, the variational approximation, or the asymptotic 
#'   sampling distribution, depending on the estimation algorithm.
#'   }
#'   \item{\code{model,x,y,z}}{
#'   If requested, the the model frame, model matrix and response variable used, 
#'   respectively. Note that z corresponds to the fixed covariates, z to the spatial aggregated covariates, and y the response.
#'   }
#'   \item{\code{family}}{
#'   The \code{\link[stats]{family}} object used.
#'   }
#'   \item{\code{call}}{
#'   The matched call.
#'   }
#'   \item{\code{formula}}{
#'   The model \code{\link[stats]{formula}}.
#'   }
#'   \item{\code{data,offset,weights}}{
#'   The \code{data}, \code{offset}, and \code{weights} arguments.
#'   }
#'   \item{\code{prior.info}}{
#'   A list with information about the prior distributions used.
#'   }
#'   \item{\code{stapfit,stan_summary}}{
#'   The object of \code{\link[rstan]{stanfit-class}} returned by RStan and a
#'   matrix of various summary statistics from the stapfit object.
#'   }
#'   \item{\code{rstan_version}}{
#'   The version of the \pkg{rstan} package that was used to fit the model.
#'   }
#' }
