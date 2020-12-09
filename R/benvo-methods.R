#' Extract start, length entries for compressed row storage format
#'
#' @param benvo built environment object from \pkg{rbenvo}
#' @param term string of term in benvo
#' @param component string of benvo term component (e.g. Distance, Time, or Distance-Time)
#'
#' @export
extract_ucrs <- function(benvo,term,component) UseMethod ("extract_ucrs")


#' @describeIn extract_ucrs
#' @export
extract_ucrs.benvo <- function(benvo,term,component){


	id <- rbenvo::get_id(benvo)
	jdf <- rbenvo::joinvo(benvo,term,component)
	u <- jdf %>% 
		dplyr::group_by_at(id) %>% 
		dplyr::count() %>% 
		dplyr::ungroup() %>%
		dplyr::select(n) %>% 
		dplyr::mutate(start = replace(dplyr::lag(cumsum(n)),
									  is.na(dplyr::lag(cumsum(n))),0) +1,
					  end = cumsum(n)) %>% 
		dplyr::select(-n) %>% as.matrix()
	names(u) <- NULL
	return(u)
}


#' Extract distance vector for compressed row storage format
#'
#' @param benvo built environment object from \pkg{rbenvo}
#' @param term string of term in benvo
#' @param component string of benvo term component (e.g. Distance, Time, or Distance-Time)
#'
#' @export
extract_dmat <- function(benvo, term, component,M) UseMethod("extract_dmat")


#' @describeIn extract_dmat
#' @export
extract_dmat.benvo <- function(benvo,term,component,M){

	if(component == "Time")
		return(NULL)
	id <- rbenvo::get_id(benvo)

	out <- benvo$sub_bef_data[[term]] %>% 
		dplyr::arrange_at(id) %>% 
		dplyr::pull(Distance)

	if(length(out != M))
		out <- matrix(c(out,rep(0,(M-length(out)))),ncol=M,nrow=1)
	out <- matrix(out,ncol=M,nrow=1)

	return(out)
}


#' Extract time vector for compressed row storage format
#'
#' @param benvo built environment object from \pkg{rbenvo}
#' @param term string of term in benvo
#' @param component string of benvo term component (e.g. Distance, Time, or Distance-Time)
#'
#' @export
extract_tmat <- function(benvo,term,component, M) UseMethod("extract_tmat")

#' @describeIn
#' @export
extract_tmat.benvo <- function(benvo,term,component,M){


	if(component == "Distance")
		return(NULL)

	out <- benvo$sub_bef_data[[term]] %>% 
		dplyr::arrange_at(id) %>% 
		dplyr::pull(Time)

	if(length(out != M))
		out <- c(out,rep(0,(M-length(out))))

	return(out)
}

