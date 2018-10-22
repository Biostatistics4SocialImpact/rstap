# Part of the rstap package for estimating model parameters
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# # This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Methods for stapreg objects
#' 
#' The methods documented on this page are actually some of the least important 
#' methods defined for \link[=stapreg-objects]{stapreg} objects. The most 
#' important methods are documented separately, each with its own page. Links to
#' those pages are provided in the \strong{See Also} section, below.
#' 
#' @name stapreg-methods
#' @aliases VarCorr fixef sigma
#' 
#' @templateVar stapregArg object,x
#' @template args-stapreg-object
#' @param ... Ignored, except by the \code{update} method. See
#'   \code{\link{update}}.
#' 
#' @details The methods documented on this page are similar to the methods 
#'   defined for objects of class 'lm', 'glm', 'glmer', etc. However there are a
#'   few key differences:
#'   
#' \describe{
#' \item{\code{residuals}}{
#' Residuals are \emph{always} of type \code{"response"} (not \code{"deviance"}
#' residuals or any other type). }
#' \item{\code{coef}}{
#' Medians are used for point estimates. See the \emph{Point estimates} section
#' in \code{\link{print.stapreg}} for more details.
#' }
#' \item{\code{se}}{
#' The \code{se} function returns standard errors based on 
#' \code{\link{mad}}. See the \emph{Uncertainty estimates} section in
#' \code{\link{print.stapreg}} for more details.
#' }
#' \item{\code{confint}}{
#' For models fit using optimization, confidence intervals are returned via a 
#' call to \code{\link[stats]{confint.default}}. If \code{algorithm} is 
#' \code{"sampling"}, \code{"meanfield"}, or \code{"fullrank"}, the
#' \code{confint} will throw an error because the
#' \code{\link{posterior_interval}} function should be used to compute Bayesian 
#' uncertainty intervals.
#' }
#' }
#' 
#' @seealso 
#' \itemize{
#'  \item The \code{\link[=print.stapreg]{print}},
#'    \code{\link[=summary.stapreg]{summary}}, and \code{\link{prior_summary}} 
#'    methods for stapreg objects for information on the fitted model.
#'  \item The \code{\link[=plot.stapreg]{plot}} method to plot estimates and
#'    diagnostics.
#'  \item The \code{\link{pp_check}} method for graphical posterior predictive
#'    checking.
#'  \item The \code{\link{posterior_predict}} and \code{\link{predictive_error}}
#'    methods for predictions and predictive errors.
#'  \item The \code{\link{posterior_interval}} and \code{\link{predictive_interval}}
#'    methods for uncertainty intervals for model parameters and predictions.
#'  \item \code{\link{log_lik}} method for  computing the log-likelihood 
#'   of (possibly new) data.
#'  \item The \code{\link[=as.matrix.stapreg]{as.matrix}}, \code{as.data.frame}, 
#'    and \code{as.array} methods to access posterior draws.
#' }
#' 
NULL

#' @rdname stapreg-methods
#' @export
coef.stapreg <- function(object, ...) {
  if (is.mer(object)) 
    return(coef_mer(object, ...))
  object$coefficients
}

#' @rdname stapreg-methods
#' @export
#' @param parm For \code{confint}, an optional character vector of parameter
#'   names.
#' @param level For \code{confint}, a scalar between \eqn{0} and \eqn{1}
#'   indicating the confidence level to use.
#'
confint.stapreg <- function(object, parm, level = 0.95, ...) {
    stop("Please use posterior_interval() to obtain", 
         "Bayesian interval estimates.", 
         call. = FALSE)
}

#' @rdname stapreg-methods
#' @export
fitted.stapreg <- function(object, ...)  
    return(object$fitted.values)

#' @rdname stapreg-methods
#' @export 
nobs.stapreg <- function(object, ...) {
  nrow(model.frame(object))
}

#' Retrieves number of Spatial-temporal aggregated predictors 
#' @export
#' @keywords internal
#' @param object A fitted model object
#' @param ... Arguments to methods.
#' @return number of spatial temporal aggregated predictors
#' @export
#' @seealso \code{\link{nstap.stapreg}} 
nstap <- function(object)
    UseMethod("nstap")

#' Retrieves number of temporal aggregated predictors 
#' @export
#' @keywords internal
#' @param object A fitted model object
#' @return number of spatial aggregated predictors
#' @export
#' @seealso \code{\link{nstap.stapreg}} 
nsap <- function(object)
    UseMethod("nsap")

#' Retrieves number of temporal aggregated predictors 
#' @export
#' @keywords internal
#' @param object A fitted model object
#' @return number of temporal aggregated predictors
#' @export
#' @seealso \code{\link{nstap.stapreg}} 
ntap <- function(object)
    UseMethod("ntap")

#' @rdname stapreg-methods
#' @export
nstap.stapreg <- function(object)
    return(object$stap_data$Q_st)

#' @rdname stapreg-methods
#' @export
ntap.stapreg <- function(object)
    return(object$stap_data$Q_t)

#' @rdname stapreg-methods
#' @export
nsap.stapreg <- function(object)
    return(object$stap_data$Q_s)

#' @export 
#' @keywords internal
#' @name nfix
#' @param object a A fitted model object
#' @param ... Arguments to methods.
#' @return number of fixed effects -  neither staps nor random effects
#' @export
#' @seealso \code{\link{nfix.stapreg}}
nfix <- function(object, ...)
    UseMethod("nfix")

#' @rdname stapreg-methods
#' @export
nfix.stapreg <- function(object,...)
    return(ncol(get_z(object)))

#' @rdname stapreg-methods
#' @export 
residuals.stapreg <- function(object, ...) {
  object$residuals
}

#' Extract standard errors
#' 
#' Generic function for extracting standard errors from fitted models.
#' 
#' @export
#' @keywords internal
#' @param object A fitted model object.
#' @param ... Arguments to methods.
#' @return Standard errors of model parameters.
#' @seealso \code{\link{se.stapreg}}
#' 
se <- function(object, ...) UseMethod("se")

#' @rdname stapreg-methods
#' @export
se.stapreg <- function(object, ...) {
  object$ses
}


#' @rdname stapreg-methods
#' @export 
#' @param correlation For \code{vcov}, if \code{FALSE} (the default) the
#'   covariance matrix is returned. If \code{TRUE}, the correlation matrix is
#'   returned instead.
#'
vcov.stapreg <- function(object, correlation = FALSE, ...) {
  out <- object$covmat
  if (!correlation) return(out)
  cov2cor(out)
}


#' @rdname stapreg-methods
#' @export
#' @export fixef
#' @importFrom lme4 fixef
#' 
fixef.stapreg <- function(object, ...) {
  coefs <- object$coefficients
  coefs <- coefs[b_names(names(coefs), invert = TRUE)]
  return(coefs)
}

#' @rdname stapreg-methods
#' @export
#' @export ngrps
#' @importFrom lme4 ngrps
#' 
ngrps.stapreg <- function(object, ...) {
  vapply(.flist(object), nlevels, 1)  
}

#' @rdname stapreg-methods
#' @export
#' @export ranef
#' @importFrom lme4 ranef
#' 
ranef.stapreg <- function(object, ...) {
  all_names <- object$stapfit@sim$fnames_oi
  sel <- b_names(all_names)
  ans <- object$stap_summary[sel, '50%']
  # avoid returning the extra levels that were included
  ans <- ans[!grepl("_NEW_", names(ans), fixed = TRUE)]
  fl <- .flist(object)
  levs <- lapply(fl, levels)
  asgn <- attr(fl, "assign")
  cnms <- .cnms(object)
  fl <- fl
  asgn <- asgn
  levs <- levs
  cnms <- cnms
  nc <- vapply(cnms, length, 1L)
  nb <- nc * vapply(levs, length, 1L)
  nbseq <- rep.int(seq_along(nb), nb)
  ml <- split(ans, nbseq)
  for (i in seq_along(ml)) {
    ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE, 
                      dimnames = list(NULL, cnms[[i]]))
  }
  ans <- lapply(seq_along(fl), function(i) {
    data.frame(do.call(cbind, ml[i]), row.names = levs[[i]], 
               check.names = FALSE)
  })
  names(ans) <- names(fl)
  structure(ans, class = "ranef.mer")
}


#' @rdname stapreg-methods
#' @export
#' @rawNamespace if(getRversion()>='3.3.0') importFrom(stats, sigma) else
#'   importFrom(lme4,sigma)
#'
sigma.stapreg <- function(object, ...) {
  if (!("sigma" %in% rownames(object$stap_summary))) 
    return(1)
  
  object$stap_summary["sigma", '50%']
}

#' @rdname stapreg-methods
#' @param sigma Ignored (included for compatibility with
#'   \code{\link[nlme]{VarCorr}}).
#' @export
#' @export VarCorr
#' @importFrom nlme VarCorr
#' @importFrom stats cov2cor
VarCorr.stapreg <- function(x, sigma = 1, ...) {
  dots <- list(...) 
  mat <- if ("stanmat" %in% names(dots)) as.matrix(dots$stanmat) else as.matrix(x)
  cnms <- .cnms(x)
  useSc <- "sigma" %in% colnames(mat)
  if (useSc) sc <- mat[,"sigma"] else sc <- 1
  Sigma <- apply(mat[,grepl("^Sigma\\[", colnames(mat)), drop = FALSE],2,median)
  nc <- vapply(cnms, FUN = length, FUN.VALUE = 1L)
  nms <- names(cnms)
  ncseq <- seq_along(nc)
  if (length(Sigma) == sum(nc * nc)) { # stanfit contains all Sigma entries
    spt <- split(Sigma, rep.int(ncseq, nc * nc))
    ans <- lapply(ncseq, function(i) {
      Sigma <- matrix(0, nc[i], nc[i])
      Sigma[,] <- spt[[i]]
      rownames(Sigma) <- colnames(Sigma) <- cnms[[i]]
      stddev <- sqrt(diag(Sigma))
      corr <- cov2cor(Sigma)
      structure(Sigma, stddev = stddev, correlation = corr)
    })       
  } else { # stanfit contains lower tri Sigma entries
    spt <- split(Sigma, rep.int(ncseq, (nc * (nc + 1)) / 2))
    ans <- lapply(ncseq, function(i) {
      Sigma <- matrix(0, nc[i], nc[i])
      Sigma[lower.tri(Sigma, diag = TRUE)] <- spt[[i]]
      Sigma <- Sigma + t(Sigma)
      diag(Sigma) <- diag(Sigma) / 2
      rownames(Sigma) <- colnames(Sigma) <- cnms[[i]]
      stddev <- sqrt(diag(Sigma))
      corr <- cov2cor(Sigma)
      structure(Sigma, stddev = stddev, correlation = corr)
    })    
  }
  names(ans) <- nms
  structure(ans, sc = mean(sc), useSc = useSc, class = "VarCorr.merMod")
}

# Exported but doc kept internal ----------------------------------------------

#' family method for stapreg objects
#'
#' @keywords internal
#' @export
#' @param object,... See \code{\link[stats]{family}}.
family.stapreg <- function(object, ...) object$family

#' model.frame method for stapreg objects
#' 
#' @keywords internal
#' @export
#' @param formula,... See \code{\link[stats]{model.frame}}.
#' @param fixed.only See \code{\link[lme4]{model.frame.merMod}}.
#' 
model.frame.stapreg <- function(formula, fixed.only = FALSE, ...) {
  if (is.mer(formula)) {
    fr <- formula$glmod$fr
    if (fixed.only) {
      ff <- formula(formula, fixed.only = TRUE)
      vars <- rownames(attr(terms.formula(ff), "factors"))
      fr <- fr[vars]
    }
    return(fr)
  }
  
  NextMethod("model.frame")
}

#' model.matrix method for stapreg objects
#' 
#' @keywords internal
#' @export
#' @param object,... See \code{\link[stats]{model.matrix}}.
#' 
model.matrix.stapreg <- function(object, subject_data) {
  if (inherits(object, "gamm4")) return(object$jam$X)
  if (is.mer(object)) return(object$glmod$X)
  return(model.matrix(terms(object), data = subject_data))
}



#' formula method for stapreg objects
#' 
#' @keywords internal
#' @export
#' @param x A stapreg object.
#' @param ... Can contain \code{fixed.only} and \code{random.only} arguments 
#'   that both default to \code{FALSE}.
#' 
formula.stapreg <- function(x, ..., m = NULL) {
  if (is.mer(x)) return(formula_mer(x, ...))
  x$formula
}

#' terms method for stapreg objects
#' @export
#' @keywords internal
#' @param x,fixed.only,random.only,... See lme4:::terms.merMod.
#' 
terms.stapreg <- function(x, ..., fixed.only = TRUE, random.only = FALSE) {
  if (!is.mer(x)){
      tt <- NextMethod("terms")
      return(tt)
  }
      
  
  fr <- x$glmod$fr
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)
  
  Terms <- attr(fr, "terms")
  if (fixed.only) {
    Terms <- terms.formula(formula(x, fixed.only = TRUE))
    attr(Terms, "predvars") <- attr(terms(fr), "predvars.fixed")
  } 
  if (random.only) {
    Terms <- terms.formula(lme4::subbars(formula.stapreg(x, random.only = TRUE)))
    attr(Terms, "predvars") <- attr(terms(fr), "predvars.random")
  }
  
  return(Terms)
}
waic <- function(x)
    UseMethod("waic")


# internal ----------------------------------------------------------------
.glmer_check <- function(object) {
  if (!is.mer(object))
    stop("This method is for stan_glmer and stan_lmer models only.", 
         call. = FALSE)
}
.cnms <- function(object, ...) UseMethod(".cnms")
.cnms.stapreg <- function(object, ...) {
  .glmer_check(object)
  object$glmod$reTrms$cnms
}
.flist <- function(object, ...) UseMethod(".flist")
.flist.stapreg <- function(object, ...) {
  .glmer_check(object)
  as.list(object$glmod$reTrms$flist)
}

coef_mer <- function(object, ...) {
  if (length(list(...))) 
    warning("Arguments named \"", paste(names(list(...)), collapse = ", "), 
            "\" ignored.", call. = FALSE)
  fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
  ref <- ranef(object)
  refnames <- unlist(lapply(ref, colnames))
  missnames <- setdiff(refnames, names(fef))
  nmiss <- length(missnames)
  if (nmiss > 0) {
    fillvars <- setNames(data.frame(rbind(rep(0, nmiss))), missnames)
    fef <- cbind(fillvars, fef)
  }
  val <- lapply(ref, function(x) fef[rep.int(1L, nrow(x)), , drop = FALSE])
  for (i in seq(a = val)) {
    refi <- ref[[i]]
    row.names(val[[i]]) <- row.names(refi)
    nmsi <- colnames(refi)
    if (!all(nmsi %in% names(fef))) 
      stop("Unable to align random and fixed effects.", call. = FALSE)
    for (nm in nmsi) 
      val[[i]][[nm]] <- val[[i]][[nm]] + refi[, nm]
  }
  structure(val, class = "coef.mer")
}

justRE <- function(f, response = FALSE) {
  response <- if (response && length(f) == 3) f[[2]] else NULL
  reformulate(paste0("(", vapply(lme4::findbars(f), 
                                 function(x) paste(deparse(x, 500L), 
                                                   collapse = " "), 
                                 ""), ")"), 
              response = response)
}
formula_mer <- function (x, fixed.only = FALSE, random.only = FALSE, ...) {
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)
  
  fr <- x$glmod$fr
  if (is.null(form <- attr(fr, "formula"))) {
    if (!grepl("lmer$", deparse(getCall(x)[[1L]]))) 
      stop("Can't find formula stored in model frame or call.", call. = FALSE)
    form <- as.formula(formula(getCall(x), ...))
  }
  if (fixed.only) {
    form <- attr(fr, "formula")
    form[[length(form)]] <- lme4::nobars(form[[length(form)]])
  }
  if (random.only)
    form <- justRE(form, response = TRUE)
  
  return(form)
}

