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
#'@param w n x k design matrix of random effects
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
#'@param subject_matrix n X N matrix that contains column indices denoting subject specific observations
#'@param subject_n n x 1 length vector that contains the number of observations for each subject
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
                            subject_matrix,
                            subject_n,
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
        subj_mat = subject_matrix,
        subj_n = subject_n,
        n = dim(subject_n)[1],
        bar_arr = if(any_bar(stap_data)) array(which(stap_data$bar_code == 1), 
                                               dim = num_bar(stap_data)) else integer(),
        num_dnd = num_dnd(stap_data),
        num_bar = num_bar(stap_data),
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
        num_s_wei = num_s_wei(stap_data),
        num_t_wei = num_t_wei(stap_data),
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
    
    if(is.null(prior_theta$dist)){
        if(stap_data$Q_st>0){
            standata$prior_dist_for_theta_s <- array(prior_theta$theta_s_dist)
            standata$prior_scale_for_theta_s <- array(prior_theta$theta_s_scale)
            standata$prior_df_for_theta_s <- array(prior_theta$theta_s_df)
            standata$prior_mean_for_theta_s <- array(prior_theta$theta_s_mean)
            standata$prior_dist_for_theta_t <- array(prior_theta$theta_t_dist)
            standata$prior_scale_for_theta_t <- array(prior_theta$theta_t_scale)
            standata$prior_mean_for_theta_t <- array(prior_theta$theta_t_mean)
            standata$prior_df_for_theta_t <- array(prior_theta$theta_t_df)
        } else if(stap_data$Q_s>0){
            standata$prior_dist_for_theta_s <- array(prior_theta$theta_s_dist)
            standata$prior_scale_for_theta_s <- array(prior_theta$theta_s_scale)
            standata$prior_df_for_theta_s <- array(prior_theta$theta_s_df)
            standata$prior_mean_for_theta_s <- array(prior_theta$theta_s_mean)
            #null entries
            standata$prior_dist_for_theta_t <- double()
            standata$prior_scale_for_theta_t <- double()
            standata$prior_df_for_theta_t <- double()
            standata$prior_mean_for_theta_t <- double()
        } else if(stap_data$Q_t>0){
            standata$prior_dist_for_theta_t <- array(prior_theta$theta_t_dist)
            standata$prior_scale_for_theta_t <- array(prior_theta$theta_t_scale)
            standata$prior_mean_for_theta_t <- array(prior_theta$theta_t_mean)
            standata$prior_df_for_theta_t <- array(prior_theta$theta_t_df)
            #null entries
            standata$prior_mean_for_theta_s <- double()
            standata$prior_dist_for_theta_s <- double()
            standata$prior_scale_for_theta_s <- double()
            standata$prior_df_for_theta_s <- double()
        }
    } else{
        if(stap_data$Q_st>0){
            standata$prior_dist_for_theta_s <- array(rep(prior_dist_for_theta,stap_data$Q_st+stap_data$Q_s))
            standata$prior_scale_for_theta_s <- array(prior_scale_for_theta)
            standata$prior_df_for_theta_s <-  array(prior_df_for_theta)
            standata$prior_mean_for_theta_s <-  array(prior_mean_for_theta)
            standata$prior_dist_for_theta_t <- array(rep(prior_dist_for_theta,stap_data$Q_st+stap_data$Q_s))
            standata$prior_scale_for_theta_t <- array(prior_scale_for_theta)
            standata$prior_mean_for_theta_t <- array(prior_mean_for_theta)
            standata$prior_df_for_theta_t <- array(prior_df_for_theta)
        } else if(stap_data$Q_t>0){
            standata$prior_dist_for_theta_t <- array(as.integer(prior_dist_for_theta))
            standata$prior_scale_for_theta_t <- array(prior_scale_for_theta)
            standata$prior_mean_for_theta_t <- array(prior_mean_for_theta)
            standata$prior_df_for_theta_t <- array(prior_df_for_theta)
            #null entries
            standata$prior_mean_for_theta_s <- double()
            standata$prior_dist_for_theta_s <- double()
            standata$prior_scale_for_theta_s <- double()
            standata$prior_df_for_theta_s <- double()
        } else if(stap_data$Q_s>0){
            standata$prior_dist_for_theta_s <- array(as.integer(prior_dist_for_theta))
            standata$prior_scale_for_theta_s <- array(prior_scale_for_theta)
            standata$prior_df_for_theta_s <-  array(prior_df_for_theta)
            standata$prior_mean_for_theta_s <-  array(prior_mean_for_theta)
            #null entries
            standata$prior_dist_for_theta_t <- double()
            standata$prior_scale_for_theta_t <- double()
            standata$prior_df_for_theta_t <- double()
            standata$prior_mean_for_theta_t <- double()
        }
    }
    
    
    #make a copy of user specification before modifying 'group' (used for keeping
    # track of priors)
    user_covariance <- if (!length(group)) NULL else group[["decov"]]
    
    if (length(group)) {
        check_reTrms(group)
        decov <- group$decov
        W <- t(group$Zt)
        group <-
            pad_reTrms(Ztlist = group$Ztlist,
                       cnms = group$cnms,
                       flist = group$flist)
        W <- group$Z
        p <- sapply(group$cnms, FUN = length)
        l <- sapply(attr(group$flist, "assign"), function(i)
            nlevels(group$flist[[i]]))
        t <- length(l)
        b_nms <- make_b_nms(group)
        g_nms <- unlist(lapply(1:t, FUN = function(i) {
            paste(group$cnms[[i]], names(group$cnms)[i], sep = "|")
        }))
        standata$t <- t
        standata$p <- as.array(p)
        standata$l <- as.array(l)
        standata$q <- ncol(W)
        standata$len_theta_L <- sum(choose(p, 2), p)
        if (is_bernoulli) {
            parts0 <- rstan::extract_sparse_parts(W[y == 0, , drop = FALSE])
            parts1 <- rstan::extract_sparse_parts(W[y == 1, , drop = FALSE])
            standata$num_non_zero <- c(length(parts0$w), length(parts1$w))
            standata$w0 <- parts0$w
            standata$w1 <- parts1$w
            standata$v0 <- parts0$v 
            standata$v1 <- parts1$v  
            standata$u0 <- parts0$u  
            standata$u1 <- parts1$u 
        } else {
            parts <- rstan::extract_sparse_parts(W)
            standata$num_non_zero <- length(parts$w)
            standata$w <- parts$w
            standata$v <- parts$v
            standata$u <- parts$u
        }
        standata$shape <- as.array(maybe_broadcast(decov$shape, t))
        standata$scale <- as.array(maybe_broadcast(decov$scale, t))
        standata$len_concentration <- sum(p[p > 1])
        standata$concentration <-
            as.array(maybe_broadcast(decov$concentration, sum(p[p > 1])))
        standata$len_regularization <- sum(p > 1)
        standata$regularization <-
            as.array(maybe_broadcast(decov$regularization, sum(p > 1)))
        standata$special_case <- all(sapply(group$cnms, FUN = function(x) {
            length(x) == 1 && x == "(Intercept)" }))
    } else { # not multilevel
        standata$t <- 0L
        standata$p <- integer(0)
        standata$l <- integer(0)
        standata$q <- 0L
        standata$len_theta_L <- 0L
        if (is_bernoulli) {
            standata$num_non_zero <- rep(0L, 2)
            standata$w0 <- standata$w1 <- double(0)
            standata$v0 <- standata$v1 <- integer(0)
            standata$u0 <- standata$u1 <- integer(0)
        } else {
            standata$num_non_zero <- 0L
            standata$w <- double(0)
            standata$v <- integer(0)
            standata$u <- integer(0)
        }
        standata$special_case <- 0L
        standata$shape <- standata$scale <- standata$concentration <-
            standata$regularization <- rep(0, 0)
        standata$len_concentration <- 0L
        standata$len_regularization <- 0L
    }
    
    if (!is_bernoulli) {
        standata$Z <- array(ztemp, dim = dim(ztemp))
        standata$y <- y
        standata$weights <- weights
        standata$offset <- offset
    }
    
    if (is_continuous) {
        standata$ub_y <- Inf
        standata$lb_y <- if (is_gaussian) -Inf else 0
        standata$prior_scale_for_aux <- prior_scale_for_aux %ORifINF% 0
        standata$prior_df_for_aux <- c(prior_df_for_aux)
        standata$prior_mean_for_aux <- c(prior_mean_for_aux)
        standata$len_y <- length(y)
        stanfit <- stanmodels$stapdnd_continuous
    } else if (is.binomial(famname)) {
        standata$prior_scale_for_aux <-
            if (!length(group) || prior_scale_for_aux == Inf)
                0 else prior_scale_for_aux
        standata$prior_mean_for_aux <- 0
        standata$prior_df_for_aux <- 0
        if (is_bernoulli) {
            y0 <- y == 0
            y1 <- y == 1
            standata$y_0 <- which(y0)
            standata$y_1 <- which(y1)
            standata$N <- c(sum(y0), sum(y1))
            standata$NN <- nrow(ztemp)
            standata$Z0 <- ztemp[y0, , drop = FALSE]
            standata$Z1 <- ztemp[y1, , drop = FALSE]
            standata$w_W0 = double(0)
            standata$v_W0 = integer(0)
            standata$u_W0 = integer(0)
            standata$w_W1 = double(0)
            standata$v_W1 = integer(0)
            standata$u_W1 = integer(0)
            if (length(weights)) {
                # nocov start
                # this code is unused because weights are interpreted as number of
                # trials for binomial glms
                standata$weights0 <- weights[y0]
                standata$weights1 <- weights[y1]
                # nocov end
            } else {
                standata$weights0 <- double(0)
                standata$weights1 <- double(0)
            }
            if (length(offset)) {
                # nocov start
                standata$offset0 <- offset[y0]
                standata$offset1 <- offset[y1]
                # nocov end
            } else {
                standata$offset0 <- double(0)
                standata$offset1 <- double(0)
            }
            # stanfit <- stanmodels$stapdnd_bernoulli
        } else {
            standata$trials <- trials
            # stanfit <- stanmodels$stapdnd_binomial
        }
    } else if (is.poisson(famname)) {
        standata$prior_scale_for_aux <- prior_scale_for_aux %ORifINF% 0
        standata$prior_mean_for_aux <- 0
        standata$prior_df_for_aux <- 0
        # stanfit <- stanmodels$stapdnd_count
    } else if (is_nb) {
        standata$prior_scale_for_aux <- prior_scale_for_aux %ORifINF% 0
        standata$prior_df_for_aux <- c(prior_df_for_aux)
        standata$prior_mean_for_aux <- c(prior_mean_for_aux)
        stanfit <- stanmodels$count
    } else if (is_gamma) {
        # nothing
    } else { # nocov start
        # family already checked above
        stop(paste(famname, "is not supported."))
    } # nocov end
    
    prior_info <- summarize_glm_prior(
        user_prior = prior_stuff,
        user_prior_intercept = prior_intercept_stuff,
        user_prior_stap = prior_stap_stuff,
        user_prior_theta = if(!is.null(prior_theta$dist)) prior_theta_stuff else prior_theta,
        user_prior_aux = prior_aux_stuff,
        has_intercept = has_intercept,
        has_predictors = nvars > 0,
        adjusted_prior_scale = prior_scale,
        adjusted_prior_intercept_scale = prior_scale_for_intercept,
        adjusted_prior_aux_scale = prior_scale_for_aux,
        family = family
    )
    
    pars <- c(if (has_intercept) "alpha",
              "delta",
              "adj_beta",
              if(any_bar(stap_data)) "beta_bar",
              if(stap_data$Q_s + stap_data$Q_st>0) "theta_s",
              if(stap_data$Q_t + stap_data$Q_st >0) "theta_t",
              if(num_s_wei(stap_data)>0) "shape_s",
              if(num_t_wei(stap_data)>0) "shape_t",
              if(length(group)) "b",
              if(standata$len_theta_L) "theta_L",
              if (is_continuous | is_nb) "aux",
              "X",
              "mean_PPD"
              )
    
    sampling_args <- set_sampling_args(
        object = stanfit,
        prior = prior,
        user_dots = list(...),
        user_adapt_delta = adapt_delta,
        data = standata,
        pars = pars,
        show_messages = FALSE)
    
    stapfit <- do.call(sampling, sampling_args)
    check <- try(check_stanfit(stapfit))
    if (!isTRUE(check)) return(standata)
    if(standata$len_theta_L){
        thetas_ref <- rstan::extract(stapfit, pars = "theta_L", inc_warmup = FALSE,
                                     permuted = FALSE)
        cnms <- group$cnms
        nc <- sapply(cnms, FUN = length)
        nms <- names(cnms)
        Sigma <- apply(thetas_ref, 1:2, FUN = function(theta) {
            Sigma <- lme4::mkVarCorr(sc = 1, cnms, nc, theta, nms)
            unlist(sapply(Sigma, simplify = FALSE,
                          FUN = function(x) x[lower.tri(x,TRUE)]))
        })
        l <- length(dim(Sigma))
        end <- utils::tail(dim(Sigma), 1L)
        shift <- grep("^theta_L", names(stapfit@sim$samples[[1]]))[1] - 1L
        if(l==3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
            stapfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain]
        }
        else for (chain in 1:end) {
            stapfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
        }
        Sigma_nms <- lapply(cnms, FUN = function(grp) {
            nm <- outer(grp, grp, FUN = paste, sep = ",")
            nm[lower.tri (nm , diag = TRUE)]
        })
        for(j in seq_along(Sigma_nms)){
            Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
        }
        Sigma_nms <- unlist(Sigma_nms)
    }
    new_names <- c(if (has_intercept) "(Intercept)",
                   colnames(ztemp),
                   coef_names(stap_data),
                   if(length(group) && length(group$flist)) c (paste0("b[", b_nms, "]")),
                   if(standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                   if (is_gaussian) "sigma",
                   if (is_gamma) "shape",
                   if (is_ig) "lambda",
                   if (is_nb) "reciprocal_dispersion",
                   paste0("X_theta_",1:(nrow(ztemp))),
                   "mean_PPD",
                   "log-posterior")
    stapfit@sim$fnames_oi <- new_names
    return(structure(stapfit, prior.info = prior_info))
}
