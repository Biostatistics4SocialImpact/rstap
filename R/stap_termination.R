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
#' @param prob A number \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired
#'   probability mass to include in the intervals. The default is to report
#'   \eqn{90}\% intervals (\code{prob=0.9}) rather than the traditionally used
#'   \eqn{95}\% (see Details).
#' @param exposure_limit A number indicating the desired
#'          amount of exposure for which the function will return an estimate of
#'           distance/time. 
#'           Note that the exposure_limit corresponds to spatial exposure and 1-temporal exposure.
#' @param max_value by defuault the max_distance and/or time from the model's original input
#'        will be used to calculate the upper bound of possible terminating 
#'        distances/times - the max_value can be used to specify a new value for this value.
#' @return A matrix with three columns and as many rows as model parameters (or
#'   the subset of parameters specified by \code{pars} and/or
#'   \code{regex_pars}). For a given value of \code{prob}, \eqn{p}, the columns
#'   correspond to the lower and upper \eqn{100p}\% interval limits and have the
#'   names \eqn{100\alpha/2}\% and \eqn{100(1 - \alpha/2)}\%, where \eqn{\alpha
#'   = 1-p}. For example, if \code{prob=0.9} is specified (a \eqn{90}\%
#'   interval), then the column names will be \code{"5\%"} and \code{"95\%"},
#'   respectively.
#'@examples
#'\dontrun{
#' fit_glm <- stap_glm(formula = y ~ sex + sap(Fast_Food),
#'                    subject_data = homog_subject_data,
#'                      distance_data = homog_distance_data,
#'                      family = gaussian(link = 'identity'),
#'                      subject_ID = 'subj_id',
#'                      prior = normal(location = 0, scale = 5, autoscale = F),
#'                      prior_intercept = normal(location = 25, scale = 5, autoscale = F),
#'                      prior_stap = normal(location = 0, scale = 3, autoscale = F),
#'                      prior_theta = log_normal(location = 1, scale = 1),
#'                      prior_aux = cauchy(location = 0,scale = 5),
#'                      max_distance = max(homog_distance_data$Distance),
#'                      chains = CHAINS, iter = ITER,
#'                      refresh = -1,verbose = F)
#'terminal_points <- stap_termination(fit_glm, prob = .9, exposure_limit = 0.01)
#'}
#'@examples
#' \dontrun{
#' fit_glm <- stap_glm(formula = y ~ sex + sap(Fast_Food),
#'                    subject_data = homog_subject_data,
#'                      distance_data = homog_distance_data,
#'                      family = gaussian(link = 'identity'),
#'                      subject_ID = 'subj_id',
#'                      prior = normal(location = 0, scale = 5, autoscale = F),
#'                      prior_intercept = normal(location = 25, scale = 5, autoscale = F),
#'                      prior_stap = normal(location = 0, scale = 3, autoscale = F),
#'                      prior_theta = log_normal(location = 1, scale = 1),
#'                      prior_aux = cauchy(location = 0,scale = 5),
#'                      max_distance = max(homog_distance_data$Distance),
#'                      chains = CHAINS, iter = ITER,
#'                      refresh = -1,verbose = F)
#'terminal_vals <- stap_termination(fit_glm, prob = .9, exposure_limit = 0.01)
#'}
#'
stap_termination <- function(object,
                             prob = .9,
                             exposure_limit= 0.05,
                             pars = NULL,
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
                                     max_value = NULL,
                                     ...){

    if(prob>=1 || prob<=0)
        stop("prob must be a number between 0 and 1")
    stap_data <- object$stap_data
    scls <- theta_names(stap_data)
    shps <- shape_names(stap_data)
    mat <- as.matrix(object)
    scls <- intersect(colnames(mat),scls)
    shps <- intersect(colnames(mat),shps)
    scls <- mat[,scls,drop=F]
    shps <- mat[,shps,drop=F]
    lower <- median <- upper <- rep(0,stap_data$Q)
    scl_ix <- 1 
    shp_ix <- 1
    max_distance <- if(!is.null(max_value)) max_value else object$max_distance
    max_time <- if(!is.null(max_value)) max_value else object$max_time
    for(ix in 1:stap_data$Q){
        shape_s <- (stap_data$weight_mats[ix,1]==5)
        shape_t <- (stap_data$weight_mats[ix,2]==6)
        if(stap_data$stap_code[ix] %in% c(0,2)){
            f <- get_weight_function(stap_data$weight_mats[scl_ix,1])
            fvec <- purrr::map2_dbl(scls[,scl_ix],shps[,shp_ix],function(x,y) .find_root(function(z){ f(z,x,y) - exposure_limit},c(0,max_distance)))
            lower[scl_ix] <- quantile(fvec, (.5 - prob/2) )
            median[scl_ix] <- median(fvec)
            upper[scl_ix] <- quantile(fvec,(.5 + prob/2))
            if(shape_s)
                shp_ix <- shp_ix + 1
            if(stap_data$stap_code[ix] == 2)
                scl_ix <- scl_ix + 1
        }
        if(stap_data$stap_code[ix] %in% c(1,2) ){
            f <- get_weight_function(stap_data$weight_mats[scl_ix,1])
            fvec <- purrr::map2_dbl(scls[,scl_ix],shps[,shp_ix],function(x,y) .find_root(function(z){ f(z,x,y) - exposure_limit},c(0,max_time)))
            lower[scl_ix] <- quantile(fvec, (.5 - prob/2) )
            median[scl_ix] <- median(fvec)
            upper[scl_ix] <- quantile(fvec,(.5 + prob/2))
            if(shape_t)
                shp_ix <- shp_ix + 1
        }
        scl_ix <- scl_ix + 1
    }
    out <- matrix(apply(cbind(lower,median,upper),1,sort),ncol=3)
    rownames(out) <- if(any_bar(stap_data)) grep("_bar",beta_names(stap_data),invert = TRUE, value = T) else beta_names(stap_data)
    colnames(out) <- c(paste0((.5-prob/2),"%"),"median",paste0((.5+prob/2),"%"))
    if(!is.null(pars))
        return(out[pars,,drop=F])
    else
        return(out)
}

.find_root <- function(f, interval)
    tryCatch(uniroot(f, interval)$root,warning = function(w){
                 print(paste("Error in root solving function, likely need to increase max value",w))}))

