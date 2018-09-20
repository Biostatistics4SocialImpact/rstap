#' Cross-Sectional Built Environment Data 
#' 
#' A dataset containing the distances between simulated subject
#' and built environment  positions within a 3 x 3 cubic square.
#' 
#' @format a dataframe with 30,250 rows and 3 variables:
#' \describe{
#'  \item{subj_id}{The subject associated with the distance}
#'  \item{BEF}{The built_environment class associated with the distance and subject}
#'  \item{dist}{The (euclidean) distance between the subject and built environment feature}
#' }
#' @source \url{https://biostatistics4socialimpact.github.io/rstap/articles/introduction.html}
"homog_distance_data"

#' Cross-Sectional Subject Data
#'
#' A dataset containing the subject specific information
#' corresponding to the built environment set-up in the homog_distance_data
#' 
#' @format a dataframe with 550 rows and 3 columns 
#' \describe{
#' \item{subj_id}{The subject unique identifier}
#' \item{y}{A simulated outcome - this example is meant to be akin to BMI}
#' \item{sex}{A simultaed binary covariate meant to represent something like sex}
#' }
#' @source \url{https:://biostatistics4socialimpact.github.io/rstap/articles/introduction.html}

#' Longitudinal Built Environment Data
#'
#' A dataset containing built environment data across time.
#' These data are simulated within a unit square and across arbitrary dates
#' with business closure and openings included in the dataset.
#' @format a dataframe with 17,875 rows and 10 columns
#' \describe{
#' \item{subj_ID}{The subject unique identifier}
#' \item{measure_ID}{The measurement unique identifier}
#' \item{bef_ID}{The Built Environment Unique identifier}
#' \item{measure_date}{The date at which the subject was measured}
#' \item{date_open}{The date at which the business opened}
#' \item{date_close}{The date at which the business may have closed; NA if the business is still open}
#' \item{date}{The date at which the subject first moved to the location associated with the distance and time with the built environment feature}
#' \item{class}{The kind of built environment feature. Only one is in the simulated dataset - "Coffee Shop"}
#' \item{dist}{The distance between the subject and BEF at the date to be associated with the measure ID}
#' \item{time}{The time for which the subject was "exposed" to the BEF at corresponding distance} 
#' }
#' @source \url{https://biostatistics4socialimpact.github.io/rstap/articles/longitudinal-I.html}
