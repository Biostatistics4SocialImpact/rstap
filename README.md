<img src = "docs/reference/figures/rstap_hex.png" height = "75" width = "75"/>
## `rstap`: Spatial-Temporal Aggregated Predictor Models Implemented in R
<!---
[![Build Status](https://travis-ci.org/Biostatistics4SocialImpact/rstap.svg?branch=master)](https://travis-ci.org/Biostatistics4SocialImpact/rstap)
-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rstap?color=green)](http://cran.r-project.org/package=rstap)

## About

This is an R package that fits spatial temporal aggregated predictor models using [Stan](http://mc-stan.org) (via the **rstan** package) for the back-end
estimation. The primary target audience is researchers interested in the effect of built environment features (BEFs) on human health, though other
applications are possible. See the package's [website](https://biostatistics4socialimpact.github.io/rstap) for an [introduction](https://biostatistics4socialimpact.github.io/rstap/articles/Introduction.html). Currently count, binomial, and continuous outcomes are supported.


## Installation

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

If installation fails, or you encounter other problems, please let us know by [filing an issue](https://github.com/biostatistics4socialimpact/rstap/issues).


## Contributing

Both examples and base code are welcome. Whether you're commiting a case study or a helping me flesh out further functionality of the package. Contact me via atpvyc at umich dot edu if interested.

## How to cite this package

Please use the citation associated with the arxiv [preprint](https://arxiv.org/abs/1812.10208).

## Acknowledgments 

This work was developed with support from NIH grant R01-HL131610 (PI: Sanchez).


