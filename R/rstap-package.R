#' The 'rstap' package.
#' 
#' @description rstap is a package that implements spatial-temporal aggregated predictor functions in R. This allows for the modeling of features impact on measured subjects that can be related either through space or time.
#' 
#' @docType package
#' @name rstap-package
#' @aliases rstap
#' @useDynLib rstap, .registration = TRUE
#'
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling 
#' @importFrom utils capture.output
#' @importFrom pracma erfc
#' @importFrom pracma erf
#' @importFrom pracma erfcinv
#' @importFrom pracma erfinv
#' @importFrom dplyr lag
#' @import stats
#' @import bayesplot
#' @import rstantools
#' @export log_lik posterior_linpred posterior_predict posterior_interval
#' @export predictive_interval predictive_error prior_summary
#' 
#' @references 
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.3. http://mc-stan.org
#' 
NULL
