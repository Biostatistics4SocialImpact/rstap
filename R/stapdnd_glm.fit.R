# Part of the rstap2 package for estimating model parameters
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# # This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the # GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Fitting Generalized Linear STAP Difference in Difference models
#'
#'@template args-adapt_delta
#'@param y n length vector or n x 2 matrix of outcomes
#'@param z n x p design matrix of subject specific covariates
#'@param dists_crs (q_s+q_st) x M matrix of distances between outcome 
#'observations and built environment features with a hypothesized spatial scale
#'@param u_s n x (q *2) matrix of compressed row storage array indices for dists_crs
#'@param times_crs (q_t+q_st) x M matrix of times where the outcome observations
#' were exposed to the built environment features with a hypothesized temporal scale
#'@param u_t n x (q*2) matrix of compressed row storage array  indices for times_crs
#'@param weight_functions a Q x 2 matrix with integers coding the appropriate weight function for each STAP
#'@param stap_data object of class "stap_data" that contains information on all the spatial-temporal predictors in the model
#'@param max_distance the upper bound on any and all distances included in the model 
#'@param max_time the upper bound on any and all times included in the model
#'@param subject_mat n X N matrix that contains column indices denoting subject specific observations
#'@param subj_n n x 1 length vector that contains the number of observations for each subject
#'@param weights weights to be added to the likelihood observation for a given subject
#'@param offset offset term to be added to the outcome for a given subject
#'@param family distributional family - only binomial gaussian or poisson currently allowed
#'@param prior,prior_intercept,prior_stap,prior_theta,prior_aux see \code{stap_glm} for more information
#'@param group list of of group terms from \code{lme4::glmod}
#'@param ... optional arguments passed to the sampler - e.g. iter,warmup, etc.
#'@export stapdnd_glm
stapdnd_glm.fit <- function(y, z, w, 
                            dists_crs, u_s,
                            times_crs, u_t,
                            subj_matrix,
                            subj_n,
                            weight_functions,
                            stap_data,
                            max_distance = max(dists_crs),
                            max_time = max(times_crs),
                            weights = rep(1, NROW(y)),
                            offset = rep(0, NROW(y)),
                            family = stats::gaussian(),
                             ...,
                             prior = normal(),
                             prior_intercept = normal(),
                             prior_stap = normal(),
                             group = list(), 
                             prior_theta = list(theta_one = normal()),
                             prior_aux = cauchy(location = 0L, scale = 5L),
                             adapt_delta = NULL){

    family <- validate_family(family)
    supported_families <- c("gaussian")
    fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
    if(!length(fam))
        stop("'family' must be one of ", paste(supported_families, collapse = ', '))
    if(max_distance < max(dists_crs))
        stop("max_distance must be the maximum possible distance amongst all distances in dists_crs")

    supported_links <- supported_glm_links(supported_families[fam])
    link <- which(supported_links == family$link)
    if(!length(link))
        # stop("'link' must be one of", paste( supported_links, collapse = ', '))
    if (binom_y_prop(y, family, weights))
        stop("To specify 'y' as proportion of successes and 'weights' as ",
             "number of trials please use stan_glm rather than calling ",
             "stan_glm.fit directly.")

    y <- validate_glm_outcome_support(y,family)
    if(is.binomial(family$family) && NCOL(y) == 2L){
        trials <- as.integer(y[,1L] + y[,2L])
        y <- as.integer(y[,1L])
    }

    # useless assignments to pass R CMD check
    has_intercept <-
        prior_df <- prior_df_for_intercept <- prior_df_for_aux  <-
        prior_dist <- prior_dist_for_intercept <- prior_dist_for_aux <- prior_mean <-
        prior_mean_for_intercept <- prior_mean_for_aux <- prior_scale <-
        prior_scale_for_intercept <- prior_scale_for_aux <- prior_autoscale <-
        prior_autoscale_for_intercept <- prior_autoscale_for_aux <- ztemp <-
        zbar <- prior_dist_for_stap <- prior_mean_for_stap <- 
        prior_scale_for_stap <- prior_df_for_stap <- prior_dist_for_theta <- 
        prior_scale_for_theta <- prior_df_for_theta <- prior_mean_for_theta <- NULL

    z_stuff <- center_z(z)
    
    for (i in names(z_stuff)) # ztemp, zbar, has_intercept
        assign(i, z_stuff[[i]])
    nvars <- ncol(ztemp)
    
    ok_dists <- nlist("normal", student_t = "t", "cauchy",
                      "laplace", "lasso", "product_normal")
    ok_intercept_dists <- ok_dists[1:3]
    ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
    
    prior_stuff <- handle_glm_prior(
        prior,
        nvars,
        link = family$link,
        default_scale = 2.5,
        ok_dists = ok_dists
    )
    
    for (i in names(prior_stuff))
        assign(i, prior_stuff[[i]])
    
    prior_intercept_stuff <- handle_glm_prior(
        prior_intercept,
        nvars = 1,
        default_scale = 10,
        link = family$link,
        ok_dists = ok_intercept_dists
    )
    
    names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), "_for_intercept")
    for (i in names(prior_intercept_stuff))
        assign(i, prior_intercept_stuff[[i]])
    
    prior_stap_stuff <- handle_glm_prior(
        prior_stap,
        nvars = stap_data$Q,
        link = family$link,
        default_scale = 2.5,
        ok_dists = ok_dists
    )
    names(prior_stap_stuff) <- paste0(names(prior_stap_stuff), "_for_stap")
    for(i in names(prior_stap_stuff))
        assign(i, prior_stap_stuff[[i]])
    
    prior_aux_stuff <-
        handle_glm_prior(
            prior_aux,
            nvars = 1,
            default_scale = 1,
            link = NULL, # don't need to adjust scale based on logit vs probit
            ok_dists = ok_aux_dists
        )
    # prior_{dist, mean, scale, df, dist_name, autoscale}_for_aux
    names(prior_aux_stuff) <- paste0(names(prior_aux_stuff), "_for_aux")
    if (is.null(prior_aux)) {
        prior_aux_stuff$prior_scale_for_aux <- Inf
    }
    for (i in names(prior_aux_stuff))
        assign(i, prior_aux_stuff[[i]])
    
    #prior_{dist, mean, scale, df, dist_name, autoscale}_for_theta
    if(is.null(prior_theta$dist))
        prior_theta  <- handle_theta_stap_prior(prior_theta,
                                                ok_dists = nlist("normal","lognormal"),
                                                stap_code = stap_data$stap_code,
                                                coef_names = grep("_scale",coef_names(stap_data),value = T, invert = T)
        )
    else{
        prior_theta_stuff <-
            handle_glm_prior(
                prior_theta,
                nvars = nrow(dists_crs),
                default_scale = 1,
                link = NULL,
                ok_dists = nlist("normal","lognormal")
            )
        names(prior_theta_stuff) <- paste0(names(prior_theta_stuff),"_for_theta")
        for(i in names(prior_theta_stuff))
            assign(i, prior_theta_stuff[[i]])
    }
    
    
    famname <- supported_families[fam]
    is_bernoulli <- is.binomial(famname) && all(y %in% 0:1)
    is_nb <- is.nb(famname)
    is_gaussian <- is.gaussian(famname)
    is_gamma <- is.gamma(famname)
    is_ig <- is.ig(famname)
    is_continuous <- is_gaussian || is_gamma || is_ig 
    
    # require intercept for certain family and link combinations
    if (!has_intercept) {
        linkname <- supported_links[link]
        needs_intercept <- !is_gaussian && linkname == "identity" ||
            is_gamma && linkname == "inverse" ||
            is.binomial(famname) && linkname == "log"
        if (needs_intercept)
            stop("To use this combination of family and link ",
                 "the model must have an intercept.")
    }
    
    if (is_gaussian) {
        ss <- sd(y)
        if (prior_dist > 0L && prior_autoscale)
            prior_scale <- ss * prior_scale
        if (prior_dist_for_intercept > 0L && prior_autoscale_for_intercept)
            prior_scale_for_intercept <-  ss * prior_scale_for_intercept
        if (prior_dist_for_aux > 0L && prior_autoscale_for_aux)
            prior_scale_for_aux <- ss * prior_scale_for_aux
    }
    if (prior_dist > 0L && prior_autoscale) {
        min_prior_scale <- 1e-12
        prior_scale <- pmax(min_prior_scale, prior_scale /
                                apply(ztemp, 2L, FUN = function(x) {
                                    num.categories <- length(unique(x))
                                    z.scale <- 1
                                    if (num.categories == 2) {
                                        z.scale <- diff(range(x))
                                    } else if (num.categories > 2) {
                                        z.scale <- sd(x)
                                    }
                                    return(z.scale)
                                }))
    }
    prior_scale <-
        as.array(pmin(.Machine$double.xmax, prior_scale))
    prior_scale_for_intercept <-
        min(.Machine$double.xmax, prior_scale_for_intercept)
    
    if (length(weights) > 0 && all(weights == 1)) weights <- double()
    if (length(offset)  > 0 && all(offset  == 0)) offset  <- double()
    if(all(is.na(u_s))){
        u_s <- array(double(),dim=c(0,0))
        dists_crs <- array(double(),dim=c(0,0))
        max_distance <- 0 
    }
    if(all(is.na(u_t))){
        u_t <- array(double(),dim=c(0,0))
        times_crs <- array(double(),dim=c(0,0))
        max_time <- 0 
    }
    
    # create entries in the data block of the .stan file
    standata <- nlist(
        subj_matrix = subj_matrix,
        subj_n = subj_n,
        N = nrow(ztemp),
        K = ncol(ztemp),
        Q = stap_data$Q, 
        Q_s = stap_data$Q_s, 
        Q_t = stap_data$Q_t,
        Q_st = stap_data$Q_st,
        weight_mat = stap_data$weight_mats,
        log_ar = stap_data$log_switch, 
        stap_code = stap_data$stap_code,
        M = ncol(dists_crs),
        zbar = as.array(zbar),
        family = stan_family_number(famname),
        link,
        max_distance = max_distance / get_space_constraint(stap_data,quantile= 0.975),
        max_time = max_time / get_time_constraint(stap_data,0.975),
        u_s = u_s,
        u_t = u_t,
        dists_crs = dists_crs,
        times_crs = times_crs,
        has_weights = length(weights) > 0,
        has_offset = length(offset) > 0,
        has_intercept,
        prior_dist,
        prior_mean,
        prior_scale,
        prior_df,
        prior_dist_for_intercept,
        prior_scale_for_intercept = c(prior_scale_for_intercept),
        prior_mean_for_intercept = c(prior_mean_for_intercept),
        prior_df_for_intercept = c(prior_df_for_intercept),
        prior_dist_for_intercept,
        prior_dist_for_stap, prior_mean_for_stap,
        prior_scale_for_stap = array(prior_scale_for_stap),
        prior_df_for_stap,
        prior_dist_for_aux = prior_dist_for_aux,
        num_normals = if(prior_dist == 7) as.integer(prior_df) else integer(0),
        num_normals_for_stap = if(prior_dist_for_stap == 7) as.integer(prior_df_for_stap) else integer(0)
        # mean,df,scale for aux added below depending on family
    )
}
