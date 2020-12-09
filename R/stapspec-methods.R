

#' Stap Specification Methods
#' 
#' @name stapspec-methods
#'
#' @param x stapspec object
#' @details methods for stapspec objects. Function names are typically self explanatory
#'


#' @rdname stapspec-methods
#' @export 
get_Q <- function(x) UseMethod("get_Q")


#' @describeIn get_Q retrieve total number of STAPS
#' @export 
get_Q.stapspec <- function(x) return(nrow(x$stapinfo))


