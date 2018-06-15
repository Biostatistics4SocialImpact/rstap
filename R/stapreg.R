
# Create a stapreg object
#
# @param object A list provided by one of the \code{stap_*} modeling functions.
# @return A stanreg object
#
stapreg <- function(object){

    stanfit <- object$stanfit
    family <- object$family
    y <- object$y
    x <- object$x
    nvars <- ncol(x)
    nobs <- NROW(y)
    ynames <- if(is.matrix(y)) rownames(y) else names(y)

}
    
