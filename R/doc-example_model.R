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

#' Example model
#' 
#' A model for use in \pkg{rstap} examples. 
#' 
#' @name example_model
#' @format Calling \code{example("example_model")} will run the model in the 
#'   Examples section, below, and the resulting stapreg object will then be
#'   available in the global environment. The \code{chains} and \code{iter}
#'   arguments are specified to make this example small in size. In practice,
#'   we recommend that they be left unspecified in order to use the default
#'   values (4 and 2000 respectively) or increased if there are convergence
#'   problems. The \code{cores} argument is optional and on a multicore system,
#'   the user may well want to set that equal to the number of chains being
#'   executed.
#'   
#' @seealso The Longituinal \href{https://biostatistics4socialimpact.github.io/rstap/articles/longitudinal-I.html}{Vignette} for \code{stap_glmer}.
#'
#' @examples
#'  ## following lines make example run faster
#' distdata <- subset(homog_longitudinal_bef_data[,c("subj_ID","measure_ID","class","dist")],subj_ID<=25)
#' timedata <- subset(homog_longitudinal_bef_data[,c("subj_ID","measure_ID","class","time")],subj_ID<=25)
#' timedata$time <- as.numeric(timedata$time)
#' subjdata <- subset(homog_longitudinal_subject_data,subj_ID<=25)
#' example_model <- 
#'   stap_glmer(y_bern ~ centered_income +  sex + centered_age + stap(Coffee_Shop) + (1|subj_ID),
#'              family = gaussian(),
#'              subject_data = subjdata,
#'              distance_data = distdata,
#'              time_data = timedata,
#'              subject_ID = 'subj_ID',
#'              group_ID = 'measure_ID',
#'              prior_intercept = normal(location = 25, scale = 4, autoscale = F),
#'              prior = normal(location = 0, scale = 4, autoscale=F),
#'              prior_stap = normal(location = 0, scale = 4),
#'              prior_theta = list(Coffee_Shop = list(spatial = log_normal(location = 1,
#'                                                                              scale = 1),
#'                                                          temporal = log_normal(location = 1,
#'                                                                                scale = 1))),
#'              max_distance = 3, max_time = 50,
#'              # chains, cores, and iter set to make the example small and fast
#'              chains = 1, iter = 5E2, cores = 1)
#' example_model
NULL
