#' Fitting Generalized Linear STKAP models
#'
#'@param y n length vector or n x 2 matrix of outcomes
#'@param Z n x p design matrix of subject specific covariates
#'@param dists_crs q x M matrix of distances between outcome observations and
#' environmental features where q is the number of spatial covariates,
#'  and M is the maximum number of environmental features amongst all q features
#'@param max_distance the upper bound of distance for which all
#'@export stap_glm
stap_glm.fit <- function(y, z, dists_crs, u_s,
                         times_crs, u_t,
                         weight_functions,
                         stap_data,
                         max_distance = 3L,
                         weights = rep(1,NROW(y)),
                         offset = rep(0, NROW(y)),
                         family = stats::gaussian(),
                         ...,
                         prior = normal(),
                         prior_intercept = normal(),
                         prior_stap = normal(),
                         group = list(), ## group terms not yet implemented
                         prior_theta = list(theta_one = normal()),
                         prior_aux = cauchy(location = 0L, scale = 5L),
                         adapt_delta = NULL){

    family <- validate_family(family)
    supported_families <- c("binomial","gaussian","poisson")
    fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
    if(!length(fam))
        stop("'family' must be one of ", paste(supported_families, collapse = ', '))
    if(max_distance < max(dists_crs))
        stop("max_distance must be the maximum possible distance amongst all distances in dists_crs")
    ## need to insert "check u" function here to make sure dimensions are correct for the crs data and corresponding index
    ## conditional on stap condition
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
        prior_autoscale_for_intercept <- prior_autoscale_for_aux <- NULL

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
        nvars = length(stap_data),
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
    prior_theta_stuff <-
        handle_glm_prior(
            prior_theta,
            nvars = nrow(dists_crs),
            default_scale = 1,
            link = NULL,
            ok_dists = nlist("normal","lognormal","beta")
        )
    names(prior_theta_stuff) <- paste0(names(prior_theta_stuff),"_for_theta")
    for(i in names(prior_theta_stuff))
        assign(i, prior_theta_stuff[[i]])


    famname <- supported_families[fam]
    is_bernoulli <- is.binomial(famname) && all(y %in% 0:1)
    is_nb <- is.nb(famname)
    is_gaussian <- is.gaussian(famname)
    is_gamma <- is.gamma(famname)
    is_ig <- is.ig(famname)
    is_beta <- is.beta(famname)
    is_continuous <- is_gaussian || is_gamma || is_ig || is_beta

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
    if((is.na(u_t))){
        u_t <- array(double(),dim=c(0,0))
        times_crs <- array(double(),dim=c(0,0))
        max_time <- 0 
    }

    # create entries in the data block of the .stan file
    standata <- nlist(
        N = nrow(ztemp),
        K = ncol(ztemp),
        Q = length(stap_data),
        Q_s = sum(sapply(stap_data, function(x) x$stap_type == 'spatial')),
        Q_t = sum(sapply(stap_data, function(x) x$stap_type == 'temporal')),
        Q_st = sum(sapply(stap_data, function(x) x$stap_type == 'spatial-temporal')),
        w = t(sapply(stap_data, function(x) x$weight_code)),
        log_ar = array(sapply(stap_data, function(x) x$log_switch), dim = length(stap_data)),
        stap_code = array(sapply(stap_data, function(x) x$stap_code), dim = length(stap_data)),
        M = ncol(dists_crs),
        zbar = as.array(zbar),
        family = stan_family_number(famname),
        link,
        max_distance = max_distance / pracma::erfcinv(0.975),
        max_time = max_time / pracma::erfcinv(0.975),
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
        prior_dist_for_stap, prior_mean_for_stap,
        prior_scale_for_stap = array(prior_scale_for_stap),
        prior_df_for_stap,
        prior_dist_for_theta,
        prior_scale_for_theta = array(prior_scale_for_theta),
        prior_mean_for_theta = array(prior_mean_for_theta),
        prior_df_for_theta = array(prior_df_for_theta),
        prior_dist_for_aux = prior_dist_for_aux,
        num_normals = if(prior_dist == 7) as.integer(prior_df) else integer(0),
        num_normals_for_stap = if(prior_dist_for_stap == 7) as.integer(prior_df_for_stap) else integer(0)
        # mean,df,scale for aux added below depending on family
    )
    #make a copy of user specification before modifying 'group' (used for keeping
    # track of priors)
    user_covariance <- if (!length(group)) NULL else group[["decov"]]

    if (length(group) && length(group$flist)) {
        if (length(group$strata)) {
            # standata$clogit <- TRUE
            # standata$J <- nlevels(group$strata)
            # standata$strata <- c(as.integer(group$strata)[y == 1],
            #                      as.integer(group$strata)[y == 0])
        }
        # check_reTrms(group)
        # decov <- group$decov
        if (is.null(group$SSfun)) {
            # standata$SSfun <- 0L
            # standata$input <- double()
            # standata$Dose <- double()
        } else {
            # standata$SSfun <- group$SSfun
            # standata$input <- group$input
            # if (group$SSfun == 5) standata$Dose <- group$Dose
            # else standata$Dose <- double()
        }
        # Z <- t(group$Zt)
        # group <-
        #     pad_reTrms(Ztlist = group$Ztlist,
        #                cnms = group$cnms,
        #                flist = group$flist)
        # Z <- group$Z
        # p <- sapply(group$cnms, FUN = length)
        # l <- sapply(attr(group$flist, "assign"), function(i)
        #     nlevels(group$flist[[i]]))
        # t <- length(l)
        # b_nms <- make_b_nms(group)
        # g_nms <- unlist(lapply(1:t, FUN = function(i) {
        #     paste(group$cnms[[i]], names(group$cnms)[i], sep = "|")
        # }))
        # standata$t <- t
        # standata$p <- as.array(p)
        # standata$l <- as.array(l)
        # standata$q <- ncol(Z)
        # standata$len_theta_L <- sum(choose(p, 2), p)
        if (is_bernoulli) {
            # parts0 <- extract_sparse_parts(Z[y == 0, , drop = FALSE])
            # parts1 <- extract_sparse_parts(Z[y == 1, , drop = FALSE])
            # standata$num_non_zero <- c(length(parts0$w), length(parts1$w))
            # standata$w0 <- parts0$w
            # standata$w1 <- parts1$w
            # standata$v0 <- parts0$v - 1L
            # standata$v1 <- parts1$v - 1L
            # standata$u0 <- parts0$u - 1L
            # standata$u1 <- parts1$u - 1L
        } else {
            # parts <- extract_sparse_parts(Z)
            # standata$num_non_zero <- length(parts$w)
            # standata$w <- parts$w
            # standata$v <- parts$v - 1L
            # standata$u <- parts$u - 1L
        }
        # standata$shape <- as.array(maybe_broadcast(decov$shape, t))
        # standata$scale <- as.array(maybe_broadcast(decov$scale, t))
        # standata$len_concentration <- sum(p[p > 1])
        # standata$concentration <-
        #     as.array(maybe_broadcast(decov$concentration, sum(p[p > 1])))
        # standata$len_regularization <- sum(p > 1)
        # standata$regularization <-
        #     as.array(maybe_broadcast(decov$regularization, sum(p > 1)))
        # standata$special_case <- all(sapply(group$cnms, FUN = function(x) {
        #     length(x) == 1 && x == "(Intercept)" }))
    } else { # not multilevel
        if (length(group)) {
            # standata$clogit <- TRUE
            # standata$J <- nlevels(group$strata)
            # standata$strata <- c(as.integer(group$strata)[y == 1],
            #                      as.integer(group$strata)[y == 0])
        }
        ## to be added later when group terms implemented
        # standata$t <- 0L
        # standata$p <- integer(0)
        # standata$l <- integer(0)
        # standata$q <- 0L
        # standata$len_theta_L <- 0L
        if (is_bernoulli) {
            # standata$num_non_zero <- rep(0L, 2)
            # standata$w0 <- standata$w1 <- double(0)
            # standata$v0 <- standata$v1 <- integer(0)
            # standata$u0 <- standata$u1 <- integer(0)
        } else {
            # standata$num_non_zero <- 0L
            # standata$w <- double(0)
            # standata$v <- integer(0)
            # standata$u <- integer(0)
        }
        # standata$special_case <- 0L
        # standata$shape <- standata$scale <- standata$concentration <-
        #     standata$regularization <- rep(0, 0)
        # standata$len_concentration <- 0L
        # standata$len_regularization <- 0L
        # standata$SSfun <- 0L
        # standata$input <- double()
        # standata$Dose <- double()
    }


    if (!is_bernoulli) {
        standata$Z <- array(ztemp, dim = dim(ztemp))
        # standata$nnz_Z <- 0L
        # standata$w_Z <- double(0)
        # standata$v_Z <- integer(0)
        # standata$u_Z <- integer(0)
        standata$y <- y
        standata$weights <- weights
        standata$offset <- offset
        # standata$K_smooth <- ncol(S)
        # standata$S <- S
        # standata$smooth_map <- smooth_map
    }

    if (is_continuous) {
        standata$ub_y <- Inf
        standata$lb_y <- if (is_gaussian) -Inf else 0
        standata$prior_scale_for_aux <- prior_scale_for_aux %ORifINF% 0
        standata$prior_df_for_aux <- c(prior_df_for_aux)
        standata$prior_mean_for_aux <- c(prior_mean_for_aux)
        standata$len_y <- length(y)
        stanfit <- stanmodels$continuous
    } else if (is.binomial(famname)) {
        standata$prior_scale_for_aux <-
            if (!length(group) || prior_scale_for_aux == Inf)
                0 else prior_scale_for_aux
        standata$prior_mean_for_aux <- 0
        standata$prior_df_for_aux <- 0
        if (is_bernoulli) {
            y0 <- y == 0
            y1 <- y == 1
            standata$N <- c(sum(y0), sum(y1))
            standata$X0 <- array(ztemp[y0, , drop = FALSE], dim = c(1, sum(y0), ncol(ztemp)))
            standata$X1 <- array(ztemp[y1, , drop = FALSE], dim = c(1, sum(y1), ncol(ztemp)))
            standata$nnz_X0 = 0L
            standata$w_X0 = double(0)
            standata$v_X0 = integer(0)
            standata$u_X0 = integer(0)
            standata$nnz_X1 = 0L
            standata$w_X1 = double(0)
            standata$v_X1 = integer(0)
            standata$u_X1 = integer(0)
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
            stanfit <- stanmodels$bernoulli
        } else {
            standata$trials <- trials
            stanfit <- stanmodels$binomial
        }
    } else if (is.poisson(famname)) {
        standata$prior_scale_for_aux <- prior_scale_for_aux %ORifINF% 0
        standata$prior_mean_for_aux <- 0
        standata$prior_df_for_aux <- 0
        stanfit <- stanmodels$count
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
              "beta",
              "theta_s",
              "theta_t",
              if (is_continuous | is_nb) "aux",
              "mean_PPD")

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
    new_names <- c(if (has_intercept) "(Intercept)",
                   colnames(ztemp),
                   rownames(dists_crs),
                   paste0(rownames(dists_crs),'_spatial_scale'),
                   if (is_gaussian) "sigma",
                   if (is_gamma) "shape",
                   if (is_ig) "lambda",
                   if (is_nb) "reciprocal_dispersion",
                   "mean_PPD",
                   "log-posterior")
    stapfit@sim$fnames_oi <- new_names
    return(structure(stapfit, prior.info = prior_info))

}

# internal ---------------------------------------------------------------------------------------------------------

# @param family_name: string naming the family
# @return character vector of supported link functions for the family
supported_glm_links <- function(family_name){
    switch(
      family_name,
      binomial = c("logit","probit","cauchit", "log","cloglog"),
      gaussian = c("identity", "log", "inverse"),
      poisson = c("log", "identity", "sqrt"),
      stop("unsupported family")
    )
}
# Family number to pass to Stan
# @param famname string naming the family
# @return an integer family code
stan_family_number <- function(famname) {
    switch(
        famname,
        "gaussian" = 1L,
        "Gamma" = 2L,
        "inverse.gaussian" = 3L,
        "beta" = 4L, ## need to get rid of this eventually
        "Beta regression" = 4L,
        "binomial" = 5L,
        "poisson" = 6L,
        "neg_binomial_2" = 7L,
        stop("Family not valid.")
    )
}

# Add extra level _NEW_ to each group
#
# @param Ztlist ranef indicator matrices
# @param cnms group$cnms
# @param flist group$flist
pad_reTrms <- function(Ztlist, cnms, flist) {
    stopifnot(is.list(Ztlist))
    l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
    p <- sapply(cnms, FUN = length)
    n <- ncol(Ztlist[[1]])
    for (i in attr(flist, "assign")) {
        levels(flist[[i]]) <- c(gsub(" ", "_", levels(flist[[i]])),
                                paste0("_NEW_", names(flist)[i]))
    }
    for (i in 1:length(p)) {
        Ztlist[[i]] <- rbind(Ztlist[[i]], Matrix(0, nrow = p[i], ncol = n, sparse = TRUE))
    }
    Z <- t(do.call(rbind, args = Ztlist))
    return(nlist(Z, cnms, flist))
}

# Drop the extra reTrms from a matrix x
#
# @param x A matrix or array (e.g. the posterior sample or matrix of summary
#   stats)
# @param columns Do the columns (TRUE) or rows (FALSE) correspond to the
#   variables?
unpad_reTrms <- function(x, ...) UseMethod("unpad_reTrms")
unpad_reTrms.default <- function(x, ...) {
    if (is.matrix(x) || is.array(x))
        return(unpad_reTrms.array(x, ...))
    keep <- !grepl("_NEW_", names(x), fixed = TRUE)
    x[keep]
}

unpad_reTrms.array <- function(x, columns = TRUE, ...) {
    ndim <- length(dim(x))
    if (ndim > 3)
        stop("'x' should be a matrix or 3-D array")

    nms <- if (columns)
        last_dimnames(x) else rownames(x)
    keep <- !grepl("_NEW_", nms, fixed = TRUE)
    if (length(dim(x)) == 2) {
        x_keep <- if (columns)
            x[, keep, drop = FALSE] else x[keep, , drop = FALSE]
    } else {
        x_keep <- if (columns)
            x[, , keep, drop = FALSE] else x[keep, , , drop = FALSE]
    }
    return(x_keep)
}

make_b_nms <- function(group, m = NULL, stub = "Long") {
    group_nms <- names(group$cnms)
    b_nms <- character()
    m_stub <- if (!is.null(m)) get_m_stub(m, stub = stub) else NULL
    for (i in seq_along(group$cnms)) {
        nm <- group_nms[i]
        nms_i <- paste(group$cnms[[i]], nm)
        levels(group$flist[[nm]]) <- gsub(" ", "_", levels(group$flist[[nm]]))
        if (length(nms_i) == 1) {
            b_nms <- c(b_nms, paste0(m_stub, nms_i, ":", levels(group$flist[[nm]])))
        } else {
            b_nms <- c(b_nms, c(t(sapply(paste0(m_stub, nms_i), paste0, ":",
                                         levels(group$flist[[nm]])))))
        }
    }
    return(b_nms)
}


# Create "prior.info" attribute needed for prior_summary()
#
# @param user_* The user's prior, prior_intercept, prior_covariance, and
#   prior_aux specifications. For prior and prior_intercept these should be
#   passed in after broadcasting the df/location/scale arguments if necessary.
# @param has_intercept T/F, does model have an intercept?
# @param has_predictors T/F, does model have predictors?
# @param adjusted_prior_*_scale adjusted scales computed if using autoscaled priors
# @param family Family object.
# @return A named list with components 'prior', 'prior_intercept', and possibly
#   'prior_covariance' and 'prior_aux' each of which itself is a list
#   containing the needed values for prior_summary.
summarize_glm_prior <-
    function(user_prior,
             user_prior_intercept,
             user_prior_aux,
             has_intercept,
             has_predictors,
             adjusted_prior_scale,
             adjusted_prior_intercept_scale,
             adjusted_prior_aux_scale,
             family) {
        rescaled_coef <-
            user_prior$prior_autoscale &&
            has_predictors &&
            !is.na(user_prior$prior_dist_name) &&
            !all(user_prior$prior_scale == adjusted_prior_scale)
        rescaled_int <-
            user_prior_intercept$prior_autoscale_for_intercept &&
            has_intercept &&
            !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
            (user_prior_intercept$prior_scale_for_intercept != adjusted_prior_intercept_scale)
        rescaled_aux <- user_prior_aux$prior_autoscale_for_aux &&
            !is.na(user_prior_aux$prior_dist_name_for_aux) &&
            (user_prior_aux$prior_scale_for_aux != adjusted_prior_aux_scale)

        if (has_predictors && user_prior$prior_dist_name %in% "t") {
            if (all(user_prior$prior_df == 1)) {
                user_prior$prior_dist_name <- "cauchy"
            } else {
                user_prior$prior_dist_name <- "student_t"
            }
        }
        if (has_intercept &&
            user_prior_intercept$prior_dist_name_for_intercept %in% "t") {
            if (all(user_prior_intercept$prior_df_for_intercept == 1)) {
                user_prior_intercept$prior_dist_name_for_intercept <- "cauchy"
            } else {
                user_prior_intercept$prior_dist_name_for_intercept <- "student_t"
            }
        }
        if (user_prior_aux$prior_dist_name_for_aux %in% "t") {
            if (all(user_prior_aux$prior_df_for_aux == 1)) {
                user_prior_aux$prior_dist_name_for_aux <- "cauchy"
            } else {
                user_prior_aux$prior_dist_name_for_aux <- "student_t"
            }
        }
        prior_list <- list(
            prior =
                if (!has_predictors) NULL else with(user_prior, list(
                    dist = prior_dist_name,
                    location = prior_mean,
                    scale = prior_scale,
                    adjusted_scale = if (rescaled_coef)
                        adjusted_prior_scale else NULL,
                    df = if (prior_dist_name %in% c
                             ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
                        prior_df else NULL
                )),
            prior_intercept =
                if (!has_intercept) NULL else with(user_prior_intercept, list(
                    dist = prior_dist_name_for_intercept,
                    location = prior_mean_for_intercept,
                    scale = prior_scale_for_intercept,
                    adjusted_scale = if (rescaled_int)
                        adjusted_prior_intercept_scale else NULL,
                    df = if (prior_dist_name_for_intercept %in% "student_t")
                        prior_df_for_intercept else NULL
                ))
        )

        aux_name <- .rename_aux(family)
        prior_list$prior_aux <- if (is.na(aux_name))
            NULL else with(user_prior_aux, list(
                dist = prior_dist_name_for_aux,
                location = if (!is.na(prior_dist_name_for_aux) &&
                               prior_dist_name_for_aux != "exponential")
                    prior_mean_for_aux else NULL,
                scale = if (!is.na(prior_dist_name_for_aux) &&
                            prior_dist_name_for_aux != "exponential")
                    prior_scale_for_aux else NULL,
                adjusted_scale = if (rescaled_aux)
                    adjusted_prior_aux_scale else NULL,
                df = if (!is.na(prior_dist_name_for_aux) &&
                         prior_dist_name_for_aux %in% "student_t")
                    prior_df_for_aux else NULL,
                rate = if (!is.na(prior_dist_name_for_aux) &&
                           prior_dist_name_for_aux %in% "exponential")
                    1 / prior_scale_for_aux else NULL,
                aux_name = aux_name
            ))

        return(prior_list)
    }

# rename aux parameter based on family
.rename_aux <- function(family) {
    fam <- family$family
    if (is.gaussian(fam)) "sigma" else
        if (is.gamma(fam)) "shape" else
            if (is.ig(fam)) "lambda" else
                if (is.nb(fam)) "reciprocal_dispersion" else NA
}
# Verify that outcome values match support implied by family object
#
# @param y outcome variable
# @param family family object
# @return y (possibly slightly modified) unless an error is thrown
#
validate_glm_outcome_support <- function(y, family){
  .is_count <- function(x) {
    all(x >= 0) && all(abs(x - round(x)) < .Machine$double.eps^0.5)
  }

  fam <- family$family

  if (fam!='binomial') {
    # make sure y has ok dimensions (matrix only allowed for binomial models)
    if (length(dim(y)) > 1) {
      if (NCOL(y) == 1) {
        y <- y[, 1]
      } else {
        stop("Except for binomial models the outcome variable ",
             "should not have multiple columns.",
             call. = FALSE)
      }
    }

    # check that values match support for non-binomial models
    if (fam!='gaussian') {
      return(y)
    } else if (fam=='poisson' && !.is_count(y)) {
      stop("All outcome values must be counts for Poisson models",
           call. = FALSE)
    }
  } else { # binomial models
    if (NCOL(y) == 1L) {
      if (is.numeric(y) || is.logical(y))
        y <-  as.integer(y)
      if (is.factor(y))
        y <- fac2bin(y)
      if (!all(y %in% c(0L, 1L)))
        stop("All outcome values must be 0 or 1 for Bernoulli models.",
             call. = FALSE)
    } else if (isTRUE(NCOL(y) == 2L)) {
      if (!.is_count(y))
        stop("All outcome values must be counts for binomial models.",
             call. = FALSE)
    } else {
      stop("For binomial models the outcome should be a vector or ",
           "a matrix with 2 columns.",
           call. = FALSE)
    }
  }
  return(y)
}

