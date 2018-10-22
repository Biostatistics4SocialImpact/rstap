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

#' Target average acceptance probability
#' 
#' Details about the \code{adapt_delta} argument to \pkg{rstap}'s modeling
#' functions - also found in the \pkg{rstanarm} documentation.
#' 
#' @name adapt_delta
#' @template reference-stan-manual
#'   
#' @details For the No-U-Turn Sampler (NUTS), the variant of Hamiltonian Monte 
#'   Carlo used used by \pkg{rstap}, \code{adapt_delta} is the target average
#'   proposal acceptance probability for adaptation. 
#'   
#'   The default value of \code{adapt_delta} is 0.95
#'   
#'   In general you should not need to change \code{adapt_delta} unless you see
#'   a warning message about divergent transitions, in which case you can
#'   increase \code{adapt_delta} from the default to a value \emph{closer} to 1
#'   (e.g. from 0.95 to 0.99, or from 0.99 to 0.999, etc). The step size used by
#'   the numerical integrator is a function of \code{adapt_delta} in that
#'   increasing \code{adapt_delta} will result in a smaller step size and fewer
#'   divergences. Increasing \code{adapt_delta} will typically result in a
#'   slower sampler, but it will always lead to a more robust sampler.
NULL
