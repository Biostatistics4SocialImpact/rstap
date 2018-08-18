# Part of the rstap package for estimating model parameters
# Copyright (C) 2013, 2014, 2015, 2016, 2017 
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

# Center a matrix x and return extra stuff
#
# @param z A design matrix
center_z <- function(z) {
  z <- as.matrix(z)
  has_intercept <- if (ncol(z) == 0) 
    FALSE else grepl("(Intercept", colnames(z)[1L], fixed = TRUE)
  
  ztemp <- if (has_intercept) z[, -1L, drop=FALSE] else z
  if (has_intercept) {
    zbar <- colMeans(ztemp)
    ztemp <- sweep(ztemp, 2, zbar, FUN = "-")
  }
  else zbar <- rep(0, ncol(ztemp))
  
  sel <- apply(ztemp, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  if (any(sel)) {
    # drop any column of x with < 2 unique values (empty interaction levels)
    # exception is column of 1s isn't dropped 
    warning("Dropped empty interaction levels: ",
            paste(colnames(ztemp)[sel], collapse = ", "))
    ztemp <- ztemp[, !sel, drop = FALSE]
    zbar <- zbar[!sel]
  }
  
  return(nlist(ztemp, zbar, has_intercept))
}

# Deal with priors
#
# @param prior A list
# @param nvars An integer indicating the number of variables
# @param default_scale Default value to use to scale if not specified by user
# @param link String naming the link function.
# @param ok_dists A list of admissible distributions.
handle_glm_prior <- function(prior, nvars, default_scale, link,
                             ok_dists = nlist("normal", student_t = "t", 
                                              "cauchy", "hs", "hs_plus", 
                                              "laplace", "lasso", "product_normal")) {
  if (!length(prior))
    return(list(prior_dist = 0L, prior_mean = as.array(rep(0, nvars)),
                prior_scale = as.array(rep(1, nvars)),
                prior_df = as.array(rep(1, nvars)), prior_dist_name = NA,
                global_prior_scale = 0, global_prior_df = 0,
                slab_df = 0, slab_scale = 0,
                prior_autoscale = FALSE))

  if (!is.list(prior)) 
    stop(sQuote(deparse(substitute(prior))), " should be a named list")
  
  prior_dist_name <- prior$dist
  prior_scale <- prior$scale
  prior_mean <- prior$location
  prior_df <- prior$df
  prior_mean[is.na(prior_mean)] <- 0
  prior_df[is.na(prior_df)] <- 1
  global_prior_scale <- 0
  global_prior_df <- 0
  slab_df <- 0
  slab_scale <- 0
  if (!prior_dist_name %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist_name %in% 
             c("normal", "t", "cauchy", "laplace", "lasso", 
               "product_normal", 'lognormal','beta')) {
    if (prior_dist_name == "normal") prior_dist <- 1L
    else if (prior_dist_name == "t") prior_dist <- 2L
    else if (prior_dist_name == "laplace") prior_dist <- 5L
    else if (prior_dist_name == "lasso") prior_dist <- 6L
    else if (prior_dist_name == "product_normal") prior_dist <- 7L
    else if(prior_dist_name == 'lognormal') prior_dist <- 8L
    else if(prior_dist_name =='beta') prior_dist <- 9L
    prior_scale <- set_prior_scale(prior_scale, default = default_scale, 
                                   link = link)
  } else if (prior_dist_name %in% c("hs", "hs_plus")) {
    prior_dist <- ifelse(prior_dist_name == "hs", 3L, 4L)
    global_prior_scale <- prior$global_scale
    global_prior_df <- prior$global_df
    slab_df <- prior$slab_df
    slab_scale <- prior$slab_scale
  } else if (prior_dist_name %in% "exponential") {
    prior_dist <- 3L # only used for scale parameters so 3 not a conflict with 3 for hs
  } 
  
  prior_df <- maybe_broadcast(prior_df, nvars)
  prior_df <- as.array(pmin(.Machine$double.xmax, prior_df))
  prior_mean <- maybe_broadcast(prior_mean, nvars)
  prior_mean <- as.array(prior_mean)
  prior_scale <- maybe_broadcast(prior_scale, nvars)

  nlist(prior_dist, 
        prior_mean, 
        prior_scale, 
        prior_df, 
        prior_dist_name, 
        global_prior_scale,
        global_prior_df,
        slab_df,
        slab_scale,
        prior_autoscale = isTRUE(prior$autoscale))
}

#' extract stap data from formula and create stapcov object 
#' 
#' @param formula that designates model expression including stap covariates 
#' @return stap_list - a list for each stap containing the covariate name, stap type
#'        weight function and log indicator 
#'
extract_stap_data <- function(formula){


   
    all_names <- all.names(formula)
    staps <- c("sap","tap","stap","sap_log","tap_log","stap_log")
    stap_covs <- all_names[which(all_names%in% staps) +1]
    stap_code <- get_stap_code(all_names,stap_covs)
    weight_code <- get_weight_code(all_names,stap_covs,stap_code)
    log_code <- sapply(all_names[which(all_names %in% staps)], function(x) grepl("_log",x))*1
    out <- lapply(1:length(stap_covs),function(x) list(covariate = stap_covs[x],
                                                    stap_type = get_stap_name(stap_code[x]),
                                                    stap_code = stap_code[x],
                                                    weight_function = get_weight_name(weight_code[x,]),
                                                    weight_code = weight_code[x,],
                                                    log_switch = log_code[x]))
    stap_data(out)
}

get_weight_name <- function(code){
    list("spatial" = weight_switch(code[1]),
         "temporal" = weight_switch(code[2]))
}

weight_switch <- function(num){
    switch(num+1,"none", "erf","cerf","exp",'cexp')
}

get_stap_name <- function(code)
    switch(code+1,"spatial","temporal","spatial-temporal")


get_weight_code <- function(all_names, stap_covs, stap_code){
    w <- matrix(0,nrow = length(stap_covs),ncol=2)
    w_codes <- list("erf"=1,"cerf"=2,"exp"=3,"cexp"=4)
    for(ix in 1:length(stap_covs)){
        temp <- all_names[which(all_names == stap_covs[ix])+1]
        if(stap_code[ix] %in% c(0,2)){
            if(temp %in% c("cerf","cexp"))
                w[ix,1] <- w_codes[[temp]]
            else
                w[ix,1] <- 2
        }else if(stap_code[ix] == 1){
            if(temp %in% c("erf","exp"))
                w[ix,2] <- w_codes[[temp]]
            else
                w[ix,2] <- 1
        }
        if(stap_code[ix] == 2){
            temp <- all_names[which(all_names == stap_covs[ix])+2]
            if(temp %in% c("erf","exp"))
               w[ix,2] <- w_codes[[temp]]
            else
                w[ix,2] <- 1
        }
    }
    return(w)
}

#' Get stap coding from formula
#'
#' @param  all_names character vector from calling all.names(formula)
#' @return vector of length equal to number of staps + saps + taps
#' with the appropriate coding for each appropriate predictor
get_stap_code <- function(all_names,stap_covs){
    staps <- list("sap"=0,"tap"=1,"stap"=2,
                  "sap_log" = 0, "tap_log" = 1, "stap_log" = 2)
    sapply(stap_covs,function(x) as.vector(staps[[all_names[which(all_names == x)-1]]]))
}

#' extract crs data
#'
#' @param stap_data_
#' @param distance_data
#' @param time_data
#' @param id_key
#' @param max_distance
extract_crs_data <- function(stap_data, distance_data, time_data, id_key, max_distance){

    dcol_ix <- validate_distancedata(distance_data,max_distance)
    tcol_ix <- validate_timedata(time_data)
    if(is.null(dcol_ix) & is.null(tcol_ix))
        stop("Neither distance_data, nor time_data submitted to function",",at least one is neccessary for rstap functions")


    if(stap_data$t_only){
        stap_covs <- stap_data$covariates
        t_col_ics <- apply(time_data, 1, function(x) which( x %in% stap_covs))
        if(!all(t_col_ics)) stop("Stap covariates must all be in (only) one column
                                 of the distance dataframe as a character or factor variable.
                                 See '?stap_glm'")
        stap_col <- colnames(time_data)[t_col_ics[1]]
        tcol <- colnames(time_data)[tcol_ix]
        tdata <- lapply(stap_covs, function(x) time_data[which(time_data[,stap_col] == x), ])
        if(any(lapply(tdata,nrow)==0)){ 
           missing <- stap_covs[which(sapply(tdata,nrow)==0)]
           stap_covs <- stap_covs[which(sapply(tdata,nrow)!=0)]
           warning(paste("The following stap covariates are not present in time_data:",
                       paste(missing, collapse = ", ")),
                   "These will be omitted from the analysis")
           tdata <- lapply(tdata,function(x) if(nrow(x)!=0) x)
           tdata[sapply(tdata,is.null)] <- NULL
       }
        M <- max(sapply(tdata,nrow))
        mddata <- lapply(tdata,function(y) merge(subject_data[,id_key, drop = F], y, by = eval(id_key),
                                                 all.x = T))
        t_mat <- lapply(mddata,function(x) x[!is.na(x[,tcol]),tcol])
        t_mat <- matrix(Reduce(rbind,lapply(t_mat,function(x) if(length(x)!=M) c(x,rep(0,M-length(x))) else x)),
                        nrow = length(stap_covs), ncol = M)
        rownames(t_mat) <- stap_covs
        freq <- lapply(mddata, function(x) xtabs(~ get(id_key) + get(stap_col),
                                                 data = x, addNA = TRUE)[,1])
        u_t <- lapply(freq,function(x) cbind(
                                             replace(dplyr::lag(cumsum(x)),
                                             is.na(dplyr::lag(cumsum(x))),0)+1,
                                             cumsum(x)))
        u_t <- abind::abind(u_t)
        dimnames(u_t) <- NULL
        return(list(d_mat = NA, t_mat = t_mat, u_t = u_t, u_s = NA))
     }else if(stap_data$d_only){
         stap_covs <- stap_data$covariates
        d_col_ics <- apply(distance_data, 1, function(x) which(x %in% stap_covs))
        if(!all(d_col_ics)) stop("Stap - of any kind - covariates must all be in (only) one column
                                 of the distance dataframe as a character or factor variable.
                                 See '?stap_glm'")
        stap_col <- colnames(distance_data)[d_col_ics[1]]
        dcol <- colnames(distance_data)[dcol_ix]

        ##ensure subjects that have zero exposure are included
        ddata <- lapply(stap_covs, function(x) distance_data[which((distance_data[,stap_col]==x &
                                                                       distance_data[,dcol]<= max_distance)),])
        if(any(lapply(ddata,nrow)==0)){
            missing <- stap_covs[which(sapply(ddata,nrow)==0)]
            stap_covs <- stap_covs[which(sapply(ddata,nrow)!=0)]
            print(paste("The following stap_covariates are not present in distance_data:",
                  paste(missing,collapse = ', ')))
            print("These will be omitted from the analysis")
            ddata <- lapply(ddata,function(x) if(nrow(x)!=0) x)
            ddata[sapply(ddata,is.null)] <- NULL
        }
        M <- max(sapply(ddata, nrow))
        mddata <- lapply(ddata,function(y) merge(subject_data[,id_key, drop = F], y, by = eval(id_key),
                                                all.x = T) )
        d_mat <- lapply(mddata,function(x) x[!is.na(x[,dcol]),dcol])
        d_mat <- matrix(Reduce(rbind,lapply(d_mat,function(x) if(length(x)!=M) c(x,rep(0,M-length(x))) else x)),
                        nrow = length(stap_covs), ncol = M)
        rownames(d_mat) <- stap_covs
        freq <- lapply(mddata, function(x) xtabs(~ get(id_key) + get(stap_col),
                                             data = x, addNA = TRUE)[,1])
        u <- lapply(freq,function(x) cbind(
            replace(dplyr::lag(cumsum(x)),
                    is.na(dplyr::lag(cumsum(x))),0)+1,
                    cumsum(x)))
        u_s <- abind::abind(u)
        dimnames(u_s) <- NULL
        return(list(d_mat = d_mat, t_mat = NA, u_t = NA, u_s = u_s))
    } else{
        sap_covs <- sap_covs(stap_data) 
        tap_covs <- tap_covs(stap_data)
        stap_covs_only <- stap_covs(stap_data)
        d_col_ics <- apply(distance_data, 1, function(x) which(x %in% sap_covs))
        t_col_ics <- apply(time_data, 1, function(x) which(x %in% tap_covs))
        if(!all(d_col_ics) && !all(t_col_ics) && !all(dst_col_ics) && !all(tst_col_ics))
            stop("Stap covariates - of any kind - must all be in (only) one column
                 of the distance dataframe as a character or factor variable. See '?stap_glm'")
        stap_dcol <- colnames(distance_data)[d_col_ics[1]]
        stap_tcol <- colnames(time_data)[t_col_ics[1]]
        dcol <- colnames(distance_data)[dcol_ix]
        tcol <- colnames(time_data)[tcol_ix]

        ##ensure subjects that have zero exposure are included
        ddata <- lapply(setdiff(sap_covs,stap_covs_only), function(x) distance_data[which((distance_data[,stap_col] == x &
                                                                                           distance_data[,dcol] <= max_distance)),])
        if(any(lapply(ddata,nrow)==0)){
            missing <- stap_covs[which(sapply(ddata,nrow)==0)]
            stap_covs <- stap_covs[which(sapply(ddata,nrow)!=0)]
            print(paste("The following stap_covariates are not present in distance_data:",
                  paste(missing,collapse = ', ')))
            print("These will be omitted from the analysis")
            ddata <- lapply(ddata,function(x) if(nrow(x)!=0) x)
            ddata[sapply(ddata,is.null)] <- NULL
        }

        tdata <- lapply(setdiff(tap_covs,stap_covs_only), function(x) time_data[which(time_data[,stap_col] == x),])

        if(any(lapply(tdata,nrow)==0)){
           missing <- stap_covs[which(sapply(tdata,nrow)==0)]
           stap_covs <- stap_covs[which(sapply(tdata,nrow)!=0)]
           warning(paste("The following stap covariates are not present in time_data:",
                       paste(missing, collapse = ", ")),
                   "These will be omitted from the analysis")
           tdata <- lapply(tdata,function(x) if(nrow(x)!=0) x)
           tdata[sapply(tdata,is.null)] <- NULL
       }
        M <-  max(sapply(tdata,nrow))
        if(M != max(sapply(ddata,nrow)))
            stop("Something wrong")

        mtdata <- lapply(tdata, function(x) merge(subject_data[,id_key, drop=F], y, by = eval(id_key),
                                                  all.x = T) )
        t_mat <- lapply(mtdata, function(x) x[!is.na(x[,tcol]),tcol, drop = F])
        t_mat <- matrix(Reduce(rbind,lapply(t_mat, function(x) if(length(x)!=M) c(x,rep(0,M-length(x))) else x)),
                        nrow= length(tap_covs), ncol = M)
        rownames(t_mat) <- tap_covs
        freq <- lapply(mtdata, function(x) xtabs(~ get(id_key) + get(stap_col), 
                                                 data = x, addNA = TRUE)[,1])
        u_t <- lapply(freq, function(x) cbind(
                                            replace(dplyr::lag(cumsum(x)),
                                                    is.na(dplyr::lag(cumsum(x))),0) +1,
                                            cumsum(x)))
        u_t <- abind::abind(u_t)
        dimnames(u_t) <- NULL
        

        mddata <- lapply(ddata, function(y) merge(subject_data[,id_key, drop = F], y, by = eval(id_key),
                                                 all.x = T))

        d_mat <- lapply(mdata, function(x) x[!is.na(x[,dcol]),dcol, drop = F])
        d_mat <- matrix(Reduce(rbind,lapply(d_mat, function(x) if (length(x)!=M) c(x,rep(0,M-length(x))) else x)),
                        nrow = length(sap_covs), ncol = M)
        rownames(d_mat) <- sap_covs
        freq <- lapply(mddata, function(x) xtabs(~get(id_key) + get(stap_col),
                                                 data= x, addNA = TRUE)[,1])
        u_s <- lapply(freq, function(x) cbind(
                                            replace(dplyr::lag(cumsum(x)),
                                            is.na(dplyr::lag(cumsum(x))),0) + 1,
                                            cumsum(x)))
        u_s <- abind::abind(u_s)
        dimnames(u_s) <- NULL
        return(list(d_mat = d_mat, t_mat = t_mat, u_s = u_s, u_t = u_t))
    }
}


