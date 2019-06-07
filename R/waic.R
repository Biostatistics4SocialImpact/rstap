#' WAIC
#'
#' @rdname waic
#' @param x a stapreg object
#'
#'
waic.stapreg <- function(x) {
    stop("waic is currently in development. This should be functioning in a future release.")

    args <- ll_args(x)
    out <- .waic(ll_fun(x), data = args$data,
                 draws = args$draws, is.mer(x))
    return(out)
}
.waic <- function(ll_fun, data, draws, is_mer){
    
    
    if(is_mer){
        stop("waic for lmer objects is still in development")
    }else{
        LPPD <- sum(sapply(1:nrow(data),function(i) log(mean(ll_fun(data[i,,drop=F], draws, log_switch = F)))))
        P_2 <- sum(sapply(1:nrow(data),function(i) var( ll_fun(data[i,,drop=F], draws, log_switch = T))))
    }   
    

   return(-2 * (LPPD - P_2)) 
}
