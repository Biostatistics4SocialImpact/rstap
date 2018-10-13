# Part of the rstap  package for estimating model parameters
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

pp_data <-
  function(object,
           newsubjdata = NULL,
           newdistdata = NULL,
           newtimedata = NULL,
           re.form = NULL,
           offset = NULL,
           ...) {
    validate_stapreg_object(object)
    if (is.mer(object)) {
       out <- .pp_data_mer(object, 
                           newsubjdata = newsubjdata, 
                           newdistdata = newdistdata,
                           newtimedata = newtimedata, 
                           re.form = re.form, ...)
       if (!is.null(offset)) out$offset <- offset
      return(out)
    }
    .pp_data(object, newsubjdata = newsubjdata, newdistdata = newdistdata,
             newtimedata = newtimedata, offset = offset, ...)
  }

# for models without lme4 structure
.pp_data <- function(object, newsubjdata = NULL, newdistdata = NULL, newtimedata = NULL,
                     offset = NULL, id_key, ...) {
  if (is.null(newsubjdata)) {
    z <- get_z(object)
    x <- get_x(object) 
    if (is.null(offset)) {
      offset <- object$offset %ORifNULL% rep(0, nrow(x))
    }
    return(nlist(z, x, offset))
  }

  offset <- .pp_data_offset(object, newsubjdata, offset)
  Terms <- delete.response(terms(object))
  m <- model.frame(Terms, newsubjdata, xlev = object$zlevels)
  z <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

  crs_data <- extract_crs_data(object$stap_data,
                               subject_data = newsubjdata,
                               distance_data = newdistdata,
                               time_data = newtimedata,
                               id_key = id_key,
                               max_distance = object$max_distance,
                               max_time = object$max_time)
  stanmat <- as.matrix(object$stapfit)
  x <- .calculate_stap_X(dists_crs = crs_data$d_mat, 
                         times_crs = crs_data$t_mat,
                         u_s = crs_data$u_s, u_t = crs_data$u_t,
                         scales = stanmat[,theta_names(object$stap_data), drop = F],
                         stap_data = object$stap_data)
  
      
  
  return(nlist(z,x, offset))
}


# for models fit using stap_(g)lmer 
.pp_data_mer <- function(object, newsubjdata, newdistdata, newtimedata,
                         re.form, id_key, ...) {

  z <- .pp_data_mer_z(object, newsubjdata,  ...)
  w <- .pp_data_mer_w(object, newsubjdata, re.form, ...)
  if(is.null(newsubjdata))
      x <- get_x(object)
  else{
      crs_data <- extract_crs_data(object$stap_data,
                                   newsubjdata,
                                   newdistdata,
                                   newtimedata,
                                   id_key,
                                   object$max_distance,
                                   object$max_time)
      stanmat <- as.matrix(object$stapfit)
      x <- .calculate_stap_X(dists_crs = crs_data$d_mat, 
                             times_crs = crs_data$t_mat,
                             u_s = crs_data$u_s, u_t = crs_data$u_t,
                             scales = stanmat[,theta_names(object$stap_data), drop = F],
                             stap_data = object$stap_data)
  }
      
  offset <- model.offset(model.frame(object))
  if (!missing(newsubjdata) && (!is.null(offset) || !is.null(object$call$offset))) {
    offset <- try(eval(object$call$offset, newdata), silent = TRUE)
    if (!is.numeric(offset)) offset <- NULL
  }
  return(nlist(z, x, offset = offset, Wt = w$Wt, w_names = w$w_names))
}


# the functions below are heavily based on a combination of 
# lme4:::predict.merMod and lme4:::mkNewReTrms, although they do also have 
# substantial modifications
.pp_data_mer_z <- function(object, newsubjdata, ...) {
  z <- get_z(object)
  if (is.null(newsubjdata)) return(z)
  form <-  attr(object$glmod$fr, "formula")
  L <- length(form)
  form[[L]] <- lme4::nobars(form[[L]])
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  Terms <- terms(object)
  mf <- model.frame(object)
  ff <- formula(form)
  vars <- rownames(attr(terms.formula(ff), "factors"))
  mf <- mf[vars]
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  mfnew <- model.frame(delete.response(Terms), newsubjdata, xlev = orig_levs)
  z <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(z, "contrasts"))
  return(z)
}

.pp_data_mer_w <- function(object, newsubjdata, re.form = NULL,
                           allow.new.levels = TRUE, na.action = na.pass, 
                           ...) {
  NAcheck <- !is.null(re.form) && !is(re.form, "formula") && is.na(re.form)
  fmla0check <- (is(re.form, "formula") && 
                   length(re.form) == 2 && 
                   identical(re.form[[2]], 0))
  if (NAcheck || fmla0check) return(list())
  if (is.null(newsubjdata) && is.null(re.form)) {
    W <- get_w(object) 
    return(list(Wt = t(W)))
  } else if (is.null(newsubjdata)) {
    rfd <- mfnew <- model.frame(object)
  } else {
    terms_fixed <- delete.response(terms(object, fixed.only = TRUE))
    mfnew <- model.frame(terms_fixed, newsubjdata, na.action = na.action)
    newdata.NA <- newsubjdata
    if (!is.null(fixed.na.action <- attr(mfnew,"na.action"))) {
      newdata.NA <- newdata.NA[-fixed.na.action,]
    }
    tt <- delete.response(terms(object, random.only = TRUE))
    rfd <- model.frame(tt, newdata.NA, na.action = na.pass)
    if (!is.null(fixed.na.action))
      attr(rfd,"na.action") <- fixed.na.action
  }
  if (is.null(re.form)) 
    re.form <- justRE(formula(object))
  if (!inherits(re.form, "formula"))
    stop("'re.form' must be NULL, NA, or a formula.")
  if (length(fit.na.action <- attr(mfnew,"na.action")) > 0) {
    newdata <- newdata[-fit.na.action,]
  }
  ReTrms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), rfd)
  if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
    stop("NAs are not allowed in prediction data",
         " for grouping variables unless 'allow.new.levels' is TRUE.")
  ns.re <- names(re <- ranef(object))
  nRnms <- names(Rcnms <- ReTrms$cnms)
  if (!all(nRnms %in% ns.re))
    stop("Grouping factors specified in re.form that were not present in original model.")
  new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
  Wt <- ReTrms$Zt
  W_names <- make_b_nms(ReTrms)
  w <- nlist(Wt = ReTrms$Zt, W_names)
  return(w)
}


# handle offsets ----------------------------------------------------------
null_or_zero <- function(x) {
  isTRUE(is.null(x) || all(x == 0))
}

.pp_data_offset <- function(object, newdata = NULL, offset = NULL) {
  if (is.null(newdata)) {
    # get offset from model object (should be null if no offset)
    if (is.null(offset)) 
      offset <- object$offset %ORifNULL% model.offset(model.frame(object))
  } else {
    if (!is.null(offset))
      stopifnot(length(offset) == nrow(newdata))
    else {
      # if newdata specified but not offset then confirm that model wasn't fit
      # with an offset (warning, not error)
      if (!is.null(object$call$offset) || 
          !null_or_zero(object$offset) || 
          !null_or_zero(model.offset(model.frame(object)))) {
        warning(
          "'offset' argument is NULL but it looks like you estimated ", 
          "the model using an offset term.", 
          call. = FALSE
        )
      }
      offset <- rep(0, nrow(newdata))
    }
  }
  return(offset)
}
