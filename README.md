[<img src = "https://avatars1.githubusercontent.com/u/28572271?s=400&u=4cfc3435602d8ad1cc847faa0000caa418713ce4&v=4" height = "42" width = "42"/>](https://biostatistics4socialimpact.github.io)
# rstap
<!---
[![Build Status](https://travis-ci.org/Biostatistics4SocialImpact/rstap.svg?branch=master)](https://travis-ci.org/Biostatistics4SocialImpact/rstap)
-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rstap?color=green)](http://cran.r-project.org/package=rstap)

## Spatial-Temporal Aggregated Predictor Models Implemented in R

This is an R package that fits spatial temporal aggregated predictor models using [Stan](http://mc-stan.org) (via the **rstan** package) for the back-end
estimation. The primary target audience is researchers interested in the effect of built environment features (BEFs) on human health, though other
applications are possible. See the package's [website](https://biostatistics4socialimpact.github.io/rstap) for an [introduction](https://biostatistics4socialimpact.github.io/rstap/articles/Introduction.html). Currently count, binomial, and continuous outcomes are supported.


### Installation

#### Development Version

To install the current development version from GitHub, first make sure that you can install the **rstan**
package and C++ toolchain by following these
[instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

Once **rstan** is successfully installed, you can install **rstap** from
GitHub using the **devtools** package by executing the following in R:

```r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("biostatistics4socialimpact/rstap")
```

Note that vignettes for this package are separately available from the 
[rstap website](https://biostatistics4socialimpact.github.io/rstap). 

If installation fails, please let us know by [filing an issue](https://github.com/biostatistics4socialimpact/rstap/issues).

