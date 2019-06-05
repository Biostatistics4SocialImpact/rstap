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

#' Create a stap_data object 
#'
#' @param object  a named list of objects containing information about the staps for a given model 
#' @return an object of class "stap_data"
stap_data <- function(object) {
    
    stap_code <- array(sapply(object, function(x) x$stap_code), dim = length(object))
    Q_t <- sum(stap_code == 1) 
    Q_s <- sum(stap_code == 0)
    Q_st <- sum(stap_code == 2)
    Q <- length(object)
    
    covariates <- sapply(object, function(x) x$covariate) 
    weight_mats <-  t(sapply(object, function(x) x$weight_code))
    log_switch <- array(sapply(object, function(x) x$log_switch), dim = Q)
    dnd_code <- array(sapply(object, function(x) x$dnd_code), dim = Q)
    bar_code <- array(sapply(object, function(x) x$bar_code), dim = Q)
    t_only <- Q_t == Q
    d_only <- Q_s == Q
    any_t <- Q_t > 0
    any_s <- Q_s > 0
    any_st <- Q_st > 0 
    dt <- (t_only == d_only)
    out <- nlist(covariates,
                 stap_code,
                 weight_mats,
                 log_switch,
                 dnd_code,
                 bar_code,
                 Q_t,
                 Q_s,
                 Q_st,
                 Q,
                 t_only,
                 d_only,
                 any_s,
                 any_t,
                 any_st,
                 dt)
    out <- structure(out, class = c("stap_data"))
    check_dups(out)
    return(out)
}

coef_names <- function(x)
    UseMethod("coef_names")

beta_names <- function(x)
    UseMethod("beta_names")

theta_names <- function(x)
    UseMethod("theta_names")

sap_covs <- function(x)
    UseMethod("sap_covs")

tap_covs <- function(x)
    UseMethod("tap_covs")

stap_covs <- function(x)
    UseMethod("stap_covs")

any_sap <- function(x)
    UseMethod("any_sap")

any_tap <- function(x)
    UseMethod("any_tap")

any_stap <- function(x)
    UseMethod("any_stap")

check_dups <- function(x)
    UseMethod("check_dups")

any_dnd <- function(x)
    UseMethod("any_dnd")

any_bar <- function(x)
    UseMethod("any_bar")

num_norm <- function(x)
    UseMethod("num_norm")

num_dnd <- function(x)
    UseMethod("num_dnd")

num_dnd <- function(object)
    return(sum(object$dnd_code))

num_bar <- function(x)
    UseMethod("num_dnd")

num_s_wei <- function(x)
    UseMethod("num_s_wei")

num_t_wei <- function(x)
    UseMethod("num_t_wei")

any_sbar <- function(x)
    UseMethod("any_sbar")

any_tbar <- function(x)
    UseMethod("any_tbar")

bar_ics <- function(x)
    UseMethod("any_tbar")

any_norm <- function(x)
    UseMethod("any_norm")

num_bar <- function(object)
    return(sum(object$bar_code))

bar_ics <- function(object)
    return(which(object$bar_code==1))

dnd_ics <- function(object)
    return(which(object$dnd_code==1))

any_norm <- function(object)
    return(num_norm(object)>0)

num_norm.stap_data <- function(object){
    bar_array <- which(object$bar_code==1)
    num_bar <- sum(object$bar_code)
    dnd_array <- which(object$dnd_code==1)
    num_dnd <- num_dnd(object)
    return(object$Q - (num_bar + num_dnd - length(intersect(bar_array,dnd_array)) ))
    
}

num_s_wei.stap_data <- function(object){

    sum(object$weight_mats[,1] == 5)
}

num_t_wei.stap_data <- function(object){

    sum(object$weight_mats[,2] == 6)
}
    
any_dnd.stap_data <- function(object)
    return(any(object$dnd_code==1))

any_bar.stap_data <- function(object)
    return(any(object$bar_code==1))

sap_covs.stap_data <- function(object)
    object$covariates[which(object$stap_code==0)]

tap_covs.stap_data <- function(object)
    object$covariates[which(object$stap_code==1)]

stap_covs.stap_data <- function(object)
    object$covariates[which(object$stap_code==2)]

any_sap.stap_data <- function(object)
    object$any_s

any_tap.stap_data <- function(object)
    object$any_t

any_stap.stap_data <- function(object)
    object$any_st

any_sbar.stap_data <- function(object){
    if(length(intersect( which(object$stap_code == 0),
                         which(object$bar_code == 1))))
        return(TRUE)
    else
        return(FALSE)
}

any_tbar.stap_data <- function(object){
    if(length(intersect( which(object$stap_code == 1),
                         which(object$bar_code == 1))))
        return(TRUE)
    else
        return(FALSE)
}

coef_names.stap_data <- function(object){
    get_name <- function(a,d,b,x,y){
        space_shape <- a[1] >4
        time_shape <- a[2] > 4
        space_time_shape <- space_shape && time_shape
        dnd <- d > 0
        bar <- b >0
        if(dnd)
            name <- paste0(y,"_dnd")
        if(bar)
            name <- paste0(y,"_bar")
        else
            name <- y
        if(space_shape)
            name <- c(name,switch(x+1,
                              paste0(y,c("_spatial_scale","_spatial_shape")),
                              paste0(y,"_temporal_scale"),
                              c(paste0(y,c("_spatial_scale","_spatial_shape")),
                                paste0(y,"_temporal_scale"))))
        else if(time_shape)
            name <- c(name,switch(x+1,
                                  paste0(y,c("_spatial_scale")),
                                  paste0(y,c("_temporal_scale","_temporal_shape")),
                                  c(paste0(y,"_spatial_scale"),
                                    paste0(y,c("_temporal_scale","_temporal_shape")) ) ))
        else if(space_time_shape)
            name <- c(name,switch(x+1,
                                  paste0(y,c("_spatial_scale","_spatial_shape")),
                                  paste0(y,c("_temporal_scale","_temporal_shape")),
                                  c(paste0(y,"_spatial_scale","_spatial_shape"),
                                    paste0(y,c("_temporal_scale","_temporal_shape")) ) ))
        else
            name <- c(name,switch(x+1,
                                  paste0(y,c("_spatial_scale")),
                                  paste0(y,c("_spatial_scale")),
                                  paste0(y,c("_spatial_scale","temporal_scale"))))
        
        
        return(name)
            
    }
    out <- as.vector(sapply(1:object$Q,function(z) get_name(object$weight_mat[z,],
                                                            object$bar_code[z],
                                                            object$dnd_code[z],
                                                            object$stap_code[z],
                                                            object$covariates[z]
                                                            )))
    out <- unlist(out)
    scales <- out[grep("_scale",out)]
    shapes <- out[grep("_shape",out)]
    out <- out[grep("_scale|_shape",out,invert=T)]
    out <- c(out,scales,shapes)
    
    return(unlist(out))
}

beta_names.stap_data <- function(object){
    nms <- coef_names(object)
    bt_nms <- nms[grep("_scale",nms,invert=T,value=T)]
    bt_nms <- bt_nms[grep("_shape",bt_nms,invert=T,value=T)]
    return(bt_nms)
}

theta_names.stap_data <- function(object){
    nms <- coef_names(object)
    th_nms <- nms[grep("_scale",nms,value=T)]
    th_nms <- th_nms[grep("shape",th_nms,value=T)]
    return(th_nms)
}

check_dups.stap_data <- function(object){
    sap <- sap_covs(object)
    tap <- tap_covs(object)
    stap <- stap_covs(object)
    if(length(sap) != length(unique(sap)))
        stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
    else if(length(tap) != length(unique(tap)))
        stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
    else if(length(stap) != length(unique(stap)))
       stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
    else if(!all.equal(union(sap,stap),c(sap,stap)))
       stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
    else if(!all.equal(union(tap,stap),c(tap,stap)))
        stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
}

coef_names.default <- function(object)
    warning("coef_names is only used for stap_data classes")

sap_covs.default <- function(object)
    warning("sap_covs is only used for stap_data classes")

tap_covs.default <- function(object)
    warning("sap_covs is only used for stap_data classes")

any_tap.default <- function(object)
    warning("any_tap is only used for stap_data classes")
