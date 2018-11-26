[<img src = "https://avatars1.githubusercontent.com/u/28572271?s=400&u=4cfc3435602d8ad1cc847faa0000caa418713ce4&v=4" height = "42" width = "42"/>](https://biostatistics4socialimpact.github.io)
# rstap
[![Build Status](https://travis-ci.org/Biostatistics4SocialImpact/rstap.svg?branch=master)](https://travis-ci.org/Biostatistics4SocialImpact/rstap)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Spatial-Temporal Aggregated Predictor Models Implemented in R

This is an R package that fits spatial temporal aggregated predictor models   
using [Stan](http://mc-stan.org) (via the **rstan** package) for the back-end
estimation. The primary target audience is researchers interested in the effect of built environment features (BEFs) on human health, though other
applications are possible. See the package's [website](https://biostatistics4socialimpact.github.io/rstap) for an [introduction](https://biostatistics4socialimpact.github.io/rstap/articles/Introduction.html).


### Installation

#### CRAN

To install the latest version uploaded to CRAN using the following code in R:

```r
install.packages("rstap", dependences = TRUE)
```

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
install_github("biostatistics4socialimpact/rstap", build_vignettes = FALSE)
```

You can switch `build_vignettes` to `TRUE` but it takes longer to install and the 
vignettes are already separately available from the 
[stap website](https://biostatistics4socialimpact.github.io/rstap). 

If installation fails, please let us know by [filing an issue](https://github.com/biostatistics4socialimpact/rstap/issues).

