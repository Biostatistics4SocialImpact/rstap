.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
  packageStartupMessage("rstap: Spatial Temporal Aggregated Predictors \n  was built under R version 3.5.0")
}
