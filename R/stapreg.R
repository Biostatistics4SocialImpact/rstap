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
    family <- object$family
    y <- object$y
    Z <- object$z
    nvars <- ncol(Z) + 2*nrow(object$dists_crs)
    n_stap_vars <- nrow(object$dists_crs)
    n_fixef_vars <- ncol(Z)
    stap_data <- object$stap_data
    nobs <- NROW(y)
    ynames <- if(is.matrix(y)) rownames(y) else names(y)
    stap_summary <- make_stap_summary(stapfit)
    coefs <- stap_summary[1:nvars, "50%"]
    sc_coefs <- coefs[grep("_scale",names(coefs))]
    stanmat <- as.matrix(stapfit)[,names(coefs), drop = F ]
    X <- calculate_stap_X(object$dists_crs,object$u_s, ## Need to adjust this to accomodate tap and stap covariates too
                                stanmat[,names(sc_coefs),drop = F])
    X_tilde <- array(NA,dim=dim(X))
    for(n_ix in 1:dim(X)[1]) X_tilde[n_ix,,] <- as.matrix(scale(X[n_ix,,]))
    colnames(stanmat) <- c(colnames(Z),
                           rownames(object$dists_crs),
                           paste0(rownames(object$dists_crs),'_spatial_scale'))
    ses <- apply(stanmat, 2L, mad)
    covmat <- cov(stanmat)
    check_rhats(stap_summary[,"Rhat"])
    delta_beta <- coefs[grep("_scale",names(coefs),invert = TRUE)]
    

    # linear predictor, fitted values
    eta <- linear_predictor(delta_beta, cbind(Z,apply(X_tilde,c(2,3),median)), object$offset)
    mu <- family$linkinv(eta)
    
    if (NCOL(y) == 2L) {
        #residuals of type 'response', (glm which does 'deviance' residuals by default)
        residuals <- y[,1L] / rowSums(y) - mu
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
        x = X,
        x_tilde = X_tilde,
        z = Z,
        n_stap_vars = n_stap_vars,
        n_fixef_vars = n_fixef_vars,
        model = object$model, 
        data = object$data, 
        family,
        offset = if (any(object$offset != 0)) object$offset else NULL,
        weights = object$weights, 
        prior.weights = object$weights, 
        contrasts = object$contrasts, 
        na.action = object$na.action,
        formula = object$formula, 
        terms = object$terms,
        prior.info = attr(stapfit, "prior.info"),
        algorithm = object$algorithm,
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
    

# stapreg-internal --------------------------------------------------------

calculate_stap_X <- function(dists_crs, u, scales){
    ags <- expand.grid(n_ix = 1:nrow(u), q_ix = nrow(dists_crs))
    X <- array(NA,dim=c(nrow(scales),nrow(u),nrow(dists_crs)))
    for(n_ix in 1:nrow(u)){
        for(q_ix in 1:nrow(dists_crs)){
            if(u[n_ix,(q_ix*2)-1]>u[n_ix,(q_ix*2)])
                X[,n_ix,q_ix] <- 0
            else
                X[,n_ix,q_ix] <- sapply(scales[,q_ix], 
                                        function(z) sum(pracma::erfc(dists_crs[u[n_ix,(q_ix*2)-1]:u[n_ix,(q_ix*2)]]/z)))
        }
    }
   return(X) 
}
