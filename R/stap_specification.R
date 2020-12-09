
#' Creates STAP model specification
#' 
stap_specification <- function(formula,benvo){


	finfo <- .get_formula_info(formula)
	stapless_formula <- finfo$stapless_formula
	terms <- finfo$stapinfo$term
	comps <- finfo$stapinfo$component
	check_terms(benvo,terms,comps)
	M <- max(purrr::map_dbl(terms,function(x) nrow(benvo$sub_bef_data[[x]])))
	u_crs <- .get_ucrs(benvo,terms,comps) 
	d_mat <- .get_dmat(benvo,terms,comps,M)
	t_mat <- .get_tmat(benvo,terms,comps,M)

	return(structure(list(stapinfo=finfo$stapinfo,
	                       stapless_formula = stapless_formula,
	                      u_s  = u_crs$u_s,
	                      u_t = u_crs$u_t,
	                      d_mat = d_mat,
	                      t_mat = t_mat),class ="stapspec"))

}


.Q <- function(x) nrow(x)
.Q_s <- function(x) sum(x$component == "Distance")
.Q_st <- function(x) sum(x$component == "Distance-Time")
.Q_t <- function(x) sum(x$component == "Time")
.num_bar <- function(x) sum(x$bw == TRUE)
.num_s_wei <- function(x) sum(x$s_function == "wei")
.num_t_wei <- function(x) sum(x$t_function == "wei")
.any_dnd <- function(x) sum(x$bw)>0
.coef_names <- function(x) c(paste0(x$term,"_scale")


check_terms <- function(benvo,terms,comps){

	purrr::map_dbl(terms,function(x,y) rbenvo::term_check(benvo,x))
	stopifnot(all(purrr::map2_lgl(terms,comps, function(x,y) y %in% rbenvo::component_lookup(benvo,x))))

}


.get_ucrs <- function(benvo,terms,component){

	s_ics <- which(component %in% c("Distance","Distance-Time"))
	u_s <- Reduce(cbind,purrr::map2(terms[s_ics],component[s_ics],function(x,y) extract_ucrs(benvo,x,y) ))

	t_ics <- which(component %in% c("Time","Distance-Time"))
	u_t <- Reduce(cbind,purrr::map2(terms[t_ics],component[t_ics],function(x,y) extract_ucrs(benvo,x,y) ))

	return(list(u_s = u_s,
				u_t = u_t))
}

.get_dmat <- function(benvo,terms,component,M){

	dmats <- Reduce(rbind,purrr::map2(terms,component,function(x,y) extract_dmat(benvo,x,y,M)))
	if(!is.null(dmats))
		dmats <- as.matrix(dmats)
	
	return(dmats)
}

.get_tmat <- function(benvo,terms,component,M){

	tmats <- Reduce(rbind,purrr::map2(terms, component, function(x,y) extract_tmat(benvo,x,y,M)))
	if(!is.null(tmats))
		tmats <- as.matrix(tmats)
	
	return(tmats)
}

.get_formula_info <- function(formula){


    with_bars <- lme4::findbars(formula)
    f <- lme4::nobars(formula)
	stapinfo <- .get_bef_info(f)
	if(is.null(stapinfo))
		stop("No covariates designated as ",paste0(c("stap","sap","tap"),collapse=","),call. = F)
	new_f <- .get_new_f(f,with_bars,stapinfo$term)

	return(list(stapless_formula = new_f,
	            stapinfo = stapinfo))

}

.get_bef_info <- function(f){


	vs <- formula.tools::rhs.vars(f)
	vs <- vs[stringr::str_detect(vs,"stap|sap|tap")]
	out <- Reduce(rbind,lapply(vs,function(x) eval(parse(text = x) ) ))
	
	return(out)
}

.get_new_f <- function(f,with_bars,not_needed){

	formula_components <- all.vars(f)[!(all.vars(f) %in% not_needed)]

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

	if(length(with_bars)){
		mer_f <- paste0(lapply(with_bars,function(x) paste0("(",deparse(x),")")),collapse = " + ")
		new_f <- paste0(new_f," + ",mer_f)
	}
	return(new_f)

}

