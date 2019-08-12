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
               "product_normal", 'lognormal')) {
    if (prior_dist_name == "normal") prior_dist <- 1L
    else if (prior_dist_name == "t") prior_dist <- 2L
    else if (prior_dist_name == "laplace") prior_dist <- 5L
    else if (prior_dist_name == "lasso") prior_dist <- 6L
    else if (prior_dist_name == "product_normal") prior_dist <- 7L
    else if(prior_dist_name == 'lognormal') prior_dist <- 8L
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


# handles stap_theta priors
#
# @param prior A list
# @param nvars An integer indicating the number of variables
# @param link String naming the link function.
# @param ok_dists A list of admissible distributions.
handle_theta_stap_prior <- function(prior,ok_dists,stap_code,default_scale,coef_names){


        if(!length(prior))
            stop("We highly reccomend against using a flat prior",
                 "on the spatial-temporal scales, as it will make
                 sampling take a very long time")

        if (!is.list(prior) ) 
            stop(sQuote(deparse(substitute(prior))), " should be a list of lists")


        check_theta_priors(prior,stap_code,coef_names)

        theta_s_dist <- rep(0,sum(stap_code==0) + sum(stap_code==2))
        theta_t_dist <- rep(0,sum(stap_code==1) + sum(stap_code==2))
        theta_s_scale <- rep(0,sum(stap_code==0) + sum(stap_code==2))
        theta_t_scale <- rep(0,sum(stap_code==1) + sum(stap_code==2))
        theta_s_mean <- rep(0,sum(stap_code==0) + sum(stap_code==2))
        theta_t_mean <- rep(0,sum(stap_code==1) + sum(stap_code==2))
        theta_s_df <- rep(0,sum(stap_code==0) + sum(stap_code==2))
        theta_t_df <- rep(0,sum(stap_code==1) + sum(stap_code==2))
        prior_theta <- list()


        for(i in 1:length(prior)){
            if(stap_code[i] %in% c(0,2)){
                prior_dist_name <- prior[[i]]$spatial$dist
                if(!prior_dist_name %in% unlist(ok_dists)){
                    stop("The prior distribution should be one of",
                         paste(names(ok_dists), collapse = ", "))
                    
              } else if (prior_dist_name %in% 
                         c("normal", "t", "cauchy", "laplace", "lasso", 
                           "product_normal", 'lognormal')) {
                if (prior_dist_name == "normal") prior_dist <- 1L
                else if (prior_dist_name == "t") prior_dist <- 2L
                else if (prior_dist_name == "laplace") prior_dist <- 5L
                else if (prior_dist_name == "lasso") prior_dist <- 6L
                else if (prior_dist_name == "product_normal") prior_dist <- 7L
                else if(prior_dist_name == 'lognormal') prior_dist <- 8L
                theta_s_dist[i] <- prior_dist
              }
                theta_s_scale[i] <- if(is.null(prior[[i]]$spatial$scale)) default_scale else prior[[i]]$spatial$scale
                theta_s_mean[i] <- if(is.na(prior[[i]]$spatial$location)) 0 else prior[[i]]$spatial$location
                theta_s_df[i] <- if(is.na(prior[[i]]$spatial$df)) 0 else prior[[i]]$spatial$df
               } 
            if(stap_code[i] %in% c(1,2)){
                prior_dist_name <- prior[[i]]$temporal$dist
                if(!prior_dist_name %in% unlist(ok_dists)){
                    stop("The prior distribution should be one of",
                         paste(names(ok_dists), collapse = ", "))
                    
                } else if (prior_dist_name %in% 
                           c("normal", "t", "cauchy", "laplace", "lasso", 
                             "product_normal", 'lognormal')) {
                    if (prior_dist_name == "normal") prior_dist <- 1L
                    else if (prior_dist_name == "t") prior_dist <- 2L
                    else if (prior_dist_name == "laplace") prior_dist <- 5L
                    else if (prior_dist_name == "lasso") prior_dist <- 6L
                    else if (prior_dist_name == "product_normal") prior_dist <- 7L
                    else if(prior_dist_name == 'lognormal') prior_dist <- 8L
                    theta_t_dist[i] <- prior_dist
                }
                theta_t_scale[i] <- if(is.null(prior[[i]]$temporal$scale)) default_scale else prior[[i]]$temporal$scale
                theta_t_mean[i] <- if(is.na(prior[[i]]$temporal$location)) 0 else prior[[i]]$temporal$location
                theta_t_df[i] <- if(is.na(prior[[i]]$temporal$df)) 0 else prior[[i]]$temporal$df
                }
            }
        if(any(stap_code %in% c(0,2))){
            prior_theta$theta_s_dist <- theta_s_dist
            prior_theta$theta_s_scale <- theta_s_scale
            prior_theta$theta_s_mean <- theta_s_mean
            prior_theta$theta_s_df <- theta_s_df
         }
        if(any(stap_code %in% c(1,2))){
            prior_theta$theta_t_dist <- theta_t_dist
            prior_theta$theta_t_scale <- theta_t_scale
            prior_theta$theta_t_mean <- theta_t_mean
            prior_theta$theta_t_df <- theta_t_df
        }
        return(prior_theta)
}

# checks theta prior if long list
# 
# @param prior_list list of stap priors
# @param stap_code 
# @param coef_names
check_theta_priors <- function(prior_list,stap_code,coef_names){

    if(any( !(names(prior_list) %in% coef_names))  | length(stap_code) != length(prior_list)  )
        stop("if assigning any individual priors for spatial-temporal scales - ALL scales must be assigned an
             appropriately named prior")

    chck <- sapply(1:length(stap_code), function(x) switch(stap_code[x]+1, names(prior_list[[x]]) == "spatial",
                                                  names(prior_list[[x]]) == "temporal",
                                                  all(names(prior_list[[x]])==c("spatial","temporal"))))
    if(!all(chck==T))
        stop("if assigning any individual priors for spatial-temporal scales -
             ALL scales must be assigned a prior and named appropriately")


}

# extract stap data from formula and create stapcov object 
# 
# @param formula that designates model expression including stap covariates 
#
extract_stap_data <- function(formula){

    all_names <- all.names(formula)
    staps <- c("sap","tap","stap")
    staps <- c(staps,paste0(staps,"_log"),paste0(staps,"_dnd"),paste0(staps,"_bar"),
               paste0(staps,"_dnd_bar"),paste0(staps,"_bar_dnd"))
    stap_covs <- all_names[which(all_names%in% staps) +1]
    if(length(stap_covs)==0)
        stop("No stap covariates specified")
    stap_code <- get_stap_code(all_names,stap_covs)
    dnd_code <- sapply(all_names[which(all_names %in% staps)], function(x) grepl("_dnd",x))*1
    bar_code <- sapply(all_names[which(all_names %in% staps)], function(x) grepl("_bar",x))*1
    weight_code <- get_weight_code(all_names,stap_covs,stap_code)
    log_code <- sapply(all_names[which(all_names %in% staps)], function(x) grepl("_log",x))*1
    out <- lapply(1:length(stap_covs),function(x){ list(covariate = stap_covs[x],
                                                    stap_type = get_stap_name(stap_code[x]),
                                                    stap_code = stap_code[x],
                                                    dnd_code = dnd_code[x],
                                                    bar_code = bar_code[x],
                                                    weight_function = get_weight_name(weight_code[x,]),
                                                    weight_code = weight_code[x,],
                                                    log_switch = log_code[x])} )
    return(stap_data(out))
}

get_weight_name <- function(code){
    list("spatial" = weight_switch(code[1]),
         "temporal" = weight_switch(code[2]))
}

weight_switch <- function(num){
    switch(num+1,"none", "erf","cerf","exp",'cexp',"wei","cwei")
}

get_stap_name <- function(code)
    switch(code+1,"spatial","temporal","spatial-temporal")


get_weight_code <- function(all_names, stap_covs, stap_code){

    w <- matrix(0,nrow = length(stap_covs),ncol=2)
    w_codes <- list("erf"=1,"cerf"=2,"exp"=3,"cexp"=4,"wei"=5,"cwei"=6)
    for(ix in 1:length(stap_covs)){
        temp <- all_names[which(all_names == stap_covs[ix])+1]
        if(stap_code[ix] %in% c(0,2)){
            if(temp %in% c("erf","cexp","cwei"))
                stop("erf,cwei and the complementary  exponential are reserved for temporal accumulation only")
            if(temp %in% c("cerf","exp","wei"))
                w[ix,1] <- w_codes[[temp]]
            else
                w[ix,1] <- 2
        }else if(stap_code[ix] == 1){
            if(temp %in% c("cerf","exp","wei"))
                stop("cerf, wei and exponential functions  are reserved for spatial decay only")
            if(temp %in% c("erf","cexp","cwei"))
                w[ix,2] <- w_codes[[temp]]
            else
                w[ix,2] <- 1
        }
        if(stap_code[ix] == 2){
            temp <- all_names[which(all_names == stap_covs[ix])+2]
            if(temp %in% c("erf","cexp","cwei"))
               w[ix,2] <- w_codes[[temp]]
            else
                w[ix,2] <- 1
        }
    }
    return(w)
}

# Get stap coding from formula
#
# @param  all_names character vector from calling all.names(formula)
# @param names of the stap covariates
# @return vector of length equal to number of staps + saps + taps
# with the appropriate coding for each appropriate predictor
get_stap_code <- function(all_names,stap_covs){

    staps <- c("sap"=0,"tap"=1,"stap"=2,
                "sap_log" = 0, "tap_log" = 1, "stap_log" = 2,
                "sap_dnd" = 0, "tap_dnd" = 1, "stap_dnd" = 2,
                "sap_bar" = 0, "tap_bar" = 1, "stap_bar" = 2,
                "sap_dnd_bar" = 0, "tap_dnd_bar"=1,"stap_dnd_bar" =2,
               "sap_bar_dnd" = 0, "tap_bar_dnd"=1,"stap_dnd_bar" =2)
    sapply(unique(stap_covs),function(x) as.vector(staps[all_names[which(all_names == x)-1]]))
}

# extract crs data
#
# @param stap_data the stap data object extracted from \code{extract_stap_data}
# @param subject_data the subject_data data.frame
# @param distance_data the distance data.frame (optional)
# @param time_data the time data.frame (optional)
# @param id_key string of the id column(s) name to join on across subject, distance and time data. 
# @param max_distance  the maximum distance in distance_data
# @param max_time  the maximum distance in time_data 
# @return a list of the crs data for the spatial and/or temporal data as appropriate 
extract_crs_data <- function(stap_data, subject_data, distance_data, 
                             time_data, id_key, max_distance, max_time){
    
    dcol_ix <- validate_distancedata(distance_data,max_distance)
    tcol_ix <- validate_timedata(time_data)
    if(is.null(dcol_ix) & is.null(tcol_ix))
        stop("Neither distance_data, nor time_data submitted to function",",at least one is neccessary for rstap functions")
    if(is.null(max_distance) & !is.null(distance_data)) max_distance <- max(distance_data[,dcol_ix])
    if(is.null(max_time) & !is.null(time_data)) max_time <- max(time_data[,tcol_ix])

    if(stap_data$t_only){
        stap_covs <- stap_data$covariates 
        t_col_ics <- unlist(apply(time_data, 1, function(x) which( x %in% stap_covs))) ## get column index for 
        .check_bef_data(t_col_ics,F)
        stap_var <- colnames(time_data)[t_col_ics[1]]
        time_var <- colnames(time_data)[tcol_ix]
        tdata <- purrr::map(stap_covs,function(x) dplyr::filter(time_data,!!dplyr::sym(stap_var) == x,
                                                                !!dplyr::sym(time_var) <= max_time))        
        mtdata <- purrr::map(tdata, function(y) dplyr::left_join(subject_data,y,by=id_key))
        M <- max(purrr::map_dbl(mtdata, nrow))
        t_mat <- .get_crs_mat(mtdata,time_var,M,stap_data$Q, stap_covs) 
        u_t <- .get_crs_u(mtdata,subject_data,id_key,stap_var,stap_covs)

        return(list(d_mat = NA, t_mat = t_mat, u_s = NA, u_t = u_t,
                    max_distance = max_distance,
                    max_time = max_time))
     }else if(stap_data$d_only){
        stap_covs <- stap_data$covariates
        d_col_ics <- unlist(apply(distance_data, 1, function(x) which(x %in% stap_covs)))
        .check_bef_data(d_col_ics)
        stap_var <- colnames(distance_data)[d_col_ics[1]]
        dist_var <- colnames(distance_data)[dcol_ix]
        ddata <- purrr::map(stap_covs,function(x) dplyr::filter(distance_data,!!dplyr::sym(stap_var) == x,!!dplyr::sym(dist_var) <=max_distance))
        mddata <- purrr::map(ddata, function(y) dplyr::left_join(subject_data,y,by=id_key))
        M <- max(purrr::map_dbl(mddata, nrow))
        d_mat <- .get_crs_mat(mddata, dist_var, M,stap_data$Q, stap_covs)
        u_s <- .get_crs_u(mddata, subject_data, id_key, stap_var, stap_covs)

        return(list(d_mat = d_mat, t_mat = NA,  u_s = u_s, u_t = NA,
                    max_distance = max_distance,
                    max_time = max_time))
    } 
    else{
        sap_covs <- sap_covs(stap_data) 
        tap_covs <- tap_covs(stap_data)
        stap_covs <- stap_covs(stap_data)
        sap_stap <- union(sap_covs,stap_covs)
        tap_stap <- union(tap_covs,stap_covs)
        d_col_ics <- unlist(apply(distance_data, 1,
                                  function(x) which(x %in% sap_stap)))
        t_col_ics <- unlist(apply(time_data, 1,
                                  function(x) which(x %in% tap_stap)))
        if(!all(d_col_ics) && !all(t_col_ics) && !all(d_col_ics) && !all(t_col_ics))
            stop("sap tap or stap covariates  must all be in (only) one column
                 of the distance dataframe as a character or factor variable. See '?stap_glm'",.call=F)
        .check_bef_data(d_col_ics)
        .check_bef_data(t_col_ics,F)
        stap_dvar <- colnames(distance_data)[d_col_ics[1]]
        stap_tvar <- colnames(time_data)[t_col_ics[1]]
        dist_var <- colnames(distance_data)[dcol_ix]
        time_var <- colnames(time_data)[tcol_ix]
        ddata <- purrr::map(sap_stap,function(x) dplyr::filter(distance_data,!!dplyr::sym(stap_dvar) == x,!!dplyr::sym(dist_var) <=max_distance))
        mddata <- purrr::map(ddata, function(y) dplyr::left_join(subject_data,y,by = id_key))
        tdata <- purrr::map(tap_stap,function(x) dplyr::filter(time_data,!!dplyr::sym(stap_tvar) == x,
                                                                !!dplyr::sym(time_var) <= max_time))
        mtdata <- purrr::map(tdata, function(y) dplyr::left_join(subject_data,y,by=id_key))
        M <- max(purrr::map_dbl(mddata, nrow))
        if(M != max(sapply(mtdata,nrow)))
            stop("Something wrong please report bug")
        t_mat <-  .get_crs_mat(mtdata, time_var, M, stap_data$Q_t + stap_data$Q_st, tap_stap)
        u_t <- .get_crs_u(mtdata, subject_data, id_key, stap_tvar, tap_stap)
        d_mat <- .get_crs_mat(mddata, dist_var, M,stap_data$Q, stap_covs)
        u_s  <- .get_crs_u(mddata,subject_data,id_key,stap_dvar,sap_stap)

        return(list(d_mat = d_mat, t_mat = t_mat, u_s = u_s, u_t = u_t, 
                    max_distance = max_distance, max_time = max_time))
    }
}


.check_bef_data <- function(col_ics,distance=T){
    if(length(col_ics)==0)
        stop("No rows in ", if(distance) "distance" else "time", " data found with designated stap covariate", .call =F)
    else if(!all(col_ics)) 
        stop("Stap covariates must all be in (only) one column of the time/distance dataframe as a character or factor variable.
                                  See '?stap_glm'")
}
        
# handle missing stap
.handle_missing_stap <- function(data,stap_covs,type_of_data = "time"){
    if(!any(purrr::map_int(data,nrow)==0))
        return(data)
    else if(any(purrr::map_lgl(data, function(x) any(is.na(x)))))
        stop("No NA data values allowed in BEF time or distance data")
    else{
        missing <- stap_covs[which(sapply(data,nrow)==0)]
        stap_covs <- stap_covs[which(sapply(data,nrow)!=0)]
        warning(paste("The following stap covariates are not present in ", type_of_data, 
                      paste(missing, collapse = ", ")), 
                "These will be ommitted from the analysis")
        data <- purrr::map(data, function(x) if(nrow(x)!=0) x)
        data[sapply(data,is.null)] <- NULL
        return(data)
    }
}



# Get time or distance csr matrix from data
.get_crs_mat <- function(list_data, col_name, M, Q,labels){

    crs_mat <- lapply(list_data, function(x) x[!is.na(x[,col_name]),col_name,drop=T])
    crs_mat <- matrix(Reduce(rbind, lapply(crs_mat,function(x) if(length(x)!=M) c(x,rep(0,M-length(x))) else x)),
                      nrow = Q, ncol = M)
    rownames(crs_mat) <- labels
    return(crs_mat)
}

.get_crs_u <- function(list_data,subject_data,id_key,stap_col,stap_covs){

    args <- dplyr::syms(id_key)
    n <- NULL ## R CMD CHECK
    id_data <- subject_data %>% dplyr::distinct(!!! args)
    ldata <- purrr::map(list_data,function(df){
        u_s <- df %>% dplyr::group_by(!!! args ) %>% 
        dplyr::count(!! dplyr::sym(stap_col)) %>% 
            dplyr::right_join(id_data,by=id_key) %>% 
            dplyr::arrange(!!! args) %>% 
            dplyr::ungroup() %>% 
            dplyr::select(n) %>%
            dplyr::mutate(start = replace(dplyr::lag(cumsum(n)),
                                          is.na(dplyr::lag(cumsum(n))),0) +1,
                          end = cumsum(n)) %>%  
            dplyr::select(-n) %>% as.matrix()
        names(u_s) <- c("","")
        return(u_s)
    })
    u_crs <- Reduce(cbind,ldata)
    return(u_crs)
}

#' get_stapless_formula
#'
#' Get formula for typical covariates
#'
#' @param f formula from stap_glm
#' @return formula without ~ stap() components
#'
get_stapless_formula <- function(f){
    
    with_bars <- f
    f <- lme4::nobars(f)
    stap_ics <- which(all.names(f)%in% c("stap","stap_log","stap_dnd_bar",
                                         "stap_dnd","stap_bar_dnd"))
    sap_ics <- which(all.names(f) %in% c("sap","sap_log","sap_dnd_bar",
                                         "sap_dnd","sap_bar_dnd"))
    tap_ics <- which(all.names(f) %in% c("tap","tap_log","tap_dnd_bar",
                                         "tap_dnd","tap_bar_dnd"))
    if(!length(stap_ics) & !length(sap_ics) & !length(tap_ics))
        stop("No covariates designated as 'stap','sap',or 'tap'  in formula", .call = F)
    stap_nms <- all.names(f)[stap_ics + 1]
    sap_nms <- all.names(f)[sap_ics + 1]
    tap_nms <- all.names(f)[tap_ics + 1]
    not_needed <- c(stap_nms,sap_nms,tap_nms,"cexp","exp","erf","cerf","wei","cwei") 
    formula_components <- all.vars(f)[!(all.vars(f) %in% not_needed)]
    bar_components <- sapply(lme4::findbars(with_bars),paste_bars)
    formula_components <- c(formula_components,bar_components)
    if(any(grepl("scale",formula_components)))
        stop("Don't use variable names with the word `scale` in them - this will cause problems with rstap methods downstream", call.=F)
    if(!attr(terms(f),"intercept"))
        formula_components <- c(formula_components,"0")
    if(grepl("cbind",all.names(f))[2]){
        new_f1 <- paste0("cbind(",formula_components[1],", ",formula_components[2], ")", " ~ ")
        ix <- 3
    }
    else{
        new_f1 <- paste0(formula_components[1],' ~ ')
        ix <- 2
    }
    new_f2 <- paste(formula_components[ix:length(formula_components)],collapse = "+")
    new_f <- paste0(new_f1,new_f2)
    return(as.formula(new_f, env = environment(f)))
}

