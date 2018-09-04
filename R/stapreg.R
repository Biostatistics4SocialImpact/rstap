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

#' Create a stapreg object
#'
#' @param object A list provided by one of the \code{stap_*} modeling functions.
#' @return A stanreg object
#'
stapreg <- function(object){

    stapfit <- object$stapfit
    stap_data <- object$stap_data
    family <- object$family
    y <- object$y
    Z <- object$z
    nvars <- ncol(Z) + stap_data$Q_s*2 + stap_data $Q_t*2 +stap_data$Q_st*3
    nobs <- NROW(y)
    ynames <- if(is.matrix(y)) rownames(y) else names(y)
    stap_summary <- make_stap_summary(stapfit)
    coefs <- stap_summary[1:nvars, "50%"]
    stanmat <- as.matrix(stapfit)[,names(coefs), drop = F ]
    x <- .calculate_stap_X(object$dists_crs, object$times_crs,
                          object$u_s, object$u_t,
                          stanmat[,coef_names(stap_data),drop=F],
                          stap_data)
    X_tilde <- array(NA, dim(x))
    for(n_ix in 1:dim(x)[1]) X_tilde[n_ix,,] <- as.matrix(scale(x[n_ix,,]))
    colnames(stanmat) <- c(colnames(Z),
                           coef_names(stap_data))
    ses <- apply(stanmat, 2L, mad)
    covmat <- cov(stanmat)
    check_rhats(stap_summary[,"Rhat"])
    delta_beta <- coefs[grep("_scale",names(coefs),invert= TRUE)]


   #linear predictor, fitted values 
    eta <- linear_predictor(delta_beta, cbind(Z,apply(X_tilde,c(2,3),median)),object$offset)
    mu <- family$linkinv(eta)


    if(NCOL(y) == 2L){
        # residuals of type 'response', (glm which does deviance residuals by default)
        residuals <- y[,1L]/ rowSums(y) - mu
    } else {
        ytmp <- if(is.factor(y)) fac2bin(y) else y
        residuals <- ytmp - mu
    }
    names(eta) <- names(mu) <- names(residuals) <- ynames
    

    out <- nlist(
        coefficients = coefs, 
        ses = ses,
        fitted.values = mu,
        linear.predictors = eta,
        residuals,
        covmat,
        y,
        x,
        z = Z,
        model = object$model,
        max_distance = object$max_distance,
        family,
        offset = if (any(object$offset != 0)) object$offset else NULL,
        weights = object$weights, 
        prior.weights = object$weights, 
        contrasts = object$contrasts, 
        na.action = object$na.action,
        formula = object$formula, 
        terms = object$terms,
        prior.info = attr(stapfit, "prior.info"),
        stap_summary,  
        stapfit = stapfit,
        stap_data = stap_data,
        rstan_version = utils::packageVersion("rstan"),
        call = object$call, 
        # sometimes 'call' is no good (e.g. if using do.call(stap_glm, args)) so
        # also include the name of the modeling function (for use when printing,
        # etc.)
        stan_function = object$stan_function
      )

    structure(out, class = c("stapreg", "glm","lm"))

}

.calculate_stap_X <- function(dists_crs, times_crs, u_s, u_t, scales, stap_data){

    n <- nrow(u_s)
    q <- stap_data$Q
    X <- array(NA,dim=c(nrow(scales), n, q))
    stap_code <- stap_data$stap_code
    

    cnt_s <- 1
    cnt_t <- 1
    scl_ix <- 1
    for(q_ix in 1:q){
        for(n_ix in 1:n){
            if(stap_code[q_ix] == 0)
                X[,n_ix,q_ix] <- assign_weight(u_s, dists_crs[cnt_s,], scales[,scl_ix], stap_data$log_switch[q_ix], 
                                               stap_data$weight_mats[q_ix,1],n_ix,cnt_s)
            else if(stap_code[q_ix] == 1)
                X[,n_ix,q_ix] <- assign_weight(u_t, times_crs[cnt_t,], scales[,scl_ix], stap_data$log_switch[q_ix], 
                                               stap_data$weight_mats[q_ix,2],n_ix,cnt_t)
            else{
                X[,n_ix,q_ix] <- assign_weight(u_s, dists_crs[cnt_s,], scales[,scl_ix], stap_data$log_switch[q_ix],
                                               stap_data$weight_mats[q_ix,1],n_ix,cnt_s)
                scl_ix <- scl_ix + 1
                X[,n_ix,q_ix] <- X[,n_ix,q_ix] * assign_weight(u_t, times_crs[cnt_t,], scales[,scl_ix], stap_data$log_switch[q_ix],
                                               stap_data$weight_mats[q_ix,2],n_ix,cnt_t)
            }
        }
        scl_ix <- scl_ix + 1
        if(stap_code[q_ix] == 0 || stap_code[q_ix] == 2)
            cnt_s <- cnt_s + 1
        else if(stap_code[q_ix] == 1 || stap_code[q_ix] == 2)
            cnt_t <- cnt_t + 1
   }
   return(X) 
}


assign_weight <- function(u, crs_data, scales, log_code, weight_code,n,q){

    w <- get_weight_function(weight_code)
    if(u[n,(q*2)-1]>u[n,(q *2)])
        return(0)
    else
        out <- sapply(scales, function(z) w(crs_data[u[n,(q*2)-1]:u[n,(q*2)]],z)) 
    if(log_code)
        return(log(sum(out)))
    else
        return(sum(out))
}
