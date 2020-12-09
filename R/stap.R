
#' Spatial Temporal Aggregated Predictor formula function
#' 
#' @param name string name
#' @param s_function one of "wei","exp","erf"
#' @param t_function one of "wei", "exp","erf" 
#'
#' @export
stap <- function(name = "",
				 s_function = "wei",
				 t_function = "wei",
				 bw = FALSE,
				 log = FALSE
				 ){

	.base_stap(name,
			   component = "Distance-Time",
			   s_function = s_function,
			   t_function = t_function,
			   bw = bw ,
			   log = log)
}

#' Spatial  Aggregated Predictor formula function
#' 
#' @param name string name
#' @param s_function one of "wei","exp","erf"
#' @param t_function one of "wei", "exp","erf" 
#'
#' @export
sap <- function(name = "",
				s_function = "wei",
				bw = FALSE,
				log = FALSE){

	.base_stap(name,
			   s_function = s_function,
			   bw = bw, 
			   log = log)
}

#' Temporal Aggregated Predictor formula function
#' 
#' @param name string name
#' @param s_function one of "wei","exp","erf"
#' @param t_function one of "wei", "exp","erf" 
#'
#' @export
tap <- function(name = "",
				t_function = "wei",
				bw = FALSE,
				log = FALSE){

	.base_stap(name,
			   component = "Time",
			   t_function = t_function)
}

.base_stap <- function(name,
					   component = "Distance",
					   s_function = NA,
					   t_function = NA,
					   bw = FALSE,
					   log = FALSE){
	check_f <- function(x){
		if(!is.na(x))
			stopifnot(x %in% c("wei","exp","erf"))
			
	}
	stopifnot(component %in% c("Distance","Time","Distance-Time"))
	check_f(s_function)
	check_f(t_function)

	return(data.frame(term = name, 
					  component = component,
					  component_code = .cstr_to_int(component),
					  s_function = s_function,
					  s_code = .fs_str_to_int(s_function),
					  t_function = t_function,
					  t_code = .ft_str_to_int(t_function),
					  bwi = bw,
					  log = log))
}

.cstr_to_int <- function(stri){
	return(dplyr::case_when(stri=="Distance" ~ 0,
							stri == "Time" ~ 1,
							stri == "Distance-Time" ~ 2))
}

.fs_str_to_int <- function(stri){
	return(dplyr::case_when(stri == "wei" ~ 5,
					 stri == "exp" ~ 3,
					 stri == "erf" ~ 2,
					 is.na(stri) ~ 0,
					 TRUE ~ 0))

}

.ft_str_to_int <- function(stri){
	return(dplyr::case_when(stri == "wei" ~ 6,
					 stri == "exp" ~ 4,
					 stri == "erf" ~ 1,
					 is.na(stri) ~ 0,
					 TRUE ~ 0))
}
