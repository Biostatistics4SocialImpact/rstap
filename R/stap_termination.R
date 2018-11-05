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

#' Spatial-Temporal Exposure Termination-Maximization Estimates
#'
#' 
#'
#' @aliases stap_termination
#' @export
#'
#' @templateVar stapregArg object
#' @template args-stapreg-object
#' @template args-dots-ignored
#' @template args-pars
#' @template args-regex-pars
#' @param prob A number \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired
#'   probability mass to include in the intervals. The default is to report
#'   \eqn{90}\% intervals (\code{prob=0.9}) rather than the traditionally used
#'   \eqn{95}\% (see Details).
#' @param exposure_limit A number indicating the desired
#'          amount of exposure for which the function will return an estimate of
#'           distance/time. 
#'           Note that the exposure_limit corresponds to spatial exposure and 1-temporal exposure.
#' @return A matrix with two columns and as many rows as model parameters (or
#'   the subset of parameters specified by \code{pars} and/or
#'   \code{regex_pars}). For a given value of \code{prob}, \eqn{p}, the columns
#'   correspond to the lower and upper \eqn{100p}\% interval limits and have the
#'   names \eqn{100\alpha/2}\% and \eqn{100(1 - \alpha/2)}\%, where \eqn{\alpha
#'   = 1-p}. For example, if \code{prob=0.9} is specified (a \eqn{90}\%
#'   interval), then the column names will be \code{"5\%"} and \code{"95\%"},
#'   respectively.
stap_termination <- function(object,
                             prob = .9,
                             exposure_limit= 0.05,
                             pars = NULL,
                             regex_pars = NULL,
                             max_value = NULL,
                             ...
                             )
    UseMethod("stap_termination")

#' @rdname stap_termination
#' @export 
stap_termination.stapreg <- function(object,
                                     prob = .9,
                                     exposure_limit = 0.05,
                                     pars = NULL,
                                     regex_pars = NULL,
                                     max_value = NULL,
                                     ...){

    stap_data <- object$stap_data
    scls <- theta_names(stap_data)
    interval <- posterior_interval(object,prob = prob)
    scls <- intersect(rownames(interval),scls)
    interval <- interval[scls,,drop=F]
    low <- interval[,1, drop=T]
    up <- interval[,2, drop = T]
    med <- apply(as.matrix(object)[,scls,drop=F],2,median)
    lower <- median <- upper <- rep(0,length(med))
    scl_ix <- 1 
    max_dist <- if(!is.null(max_value)) max_value else object$max_distance
    max_time <- if(!is.null(max_value)) max_value else object$max_time
    for(ix in 1:stap_data$Q){
        if(stap_data$stap_code[ix] %in% c(0,2)){
            f <- get_weight_function(stap_data$weight_mats[scl_ix,1])
            lower[scl_ix] <- .find_root(function(x){ f(x,low[scl_ix]) - exposure_limit}, c(0,max_distance))
            median[scl_ix] <- .find_root(function(x){ f(x,med[scl_ix]) - exposure_limit}, c(0,max_distance))
            upper[scl_ix] <- .find_root(function(x){ f(x,up[scl_ix]) - exposure_limit},  c(0,max_distance))
            if(stap_data$stap_code[ix] == 2)
                scl_ix <- scl_ix + 1
        }
        if(stap_data$stap_code[ix] %in% c(1,2)){
            f <- get_weight_function(stap_data$weight_mats[ix,2])
            lower[scl_ix] <- .find_root(function(x){ f(x,low[scl_ix]) - (1-exposure_limit)}, c(0,max_time))
            median[scl_ix] <- .find_root(function(x){ f(x,med[scl_ix]) - (1-exposure_limit)}, c(0,max_time))
            upper[scl_ix] <- .find_root(function(x){ f(x,up[scl_ix]) - (1-exposure_limit)},  c(0,max_time)) 
        }
        scl_ix <- scl_ix + 1
    }
    out <- cbind(lower,median,upper)
    rownames(out) <- scls
    sel <- intersect(scls,rownames(posterior_interval(object,pars=pars,regex_pars=regex_pars)))
    out <- out[sel,,drop=F]
    return(out)
}

.find_root <- function(f, interval)
    uniroot(f, interval)$root

