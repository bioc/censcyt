# censcyt

[![Build Status](https://travis-ci.org/retogerber/censcyt.svg?branch=main)](https://travis-ci.org/retogerber/censcyt)
[![codecov](https://codecov.io/gh/retogerber/censcyt/branch/main/graph/badge.svg)](https://codecov.io/gh/retogerber/censcyt)

## Summary

`censcyt` is an R package extending the [diffcyt](https://github.com/lmweber/diffcyt) pipeline for differential discovery in high-dimensional cytometry via high-resolution clustering by including differential abundance analysis in the presence of a covariate subject to right censoring.


## Vignettes

The main workflow can be found in the Bioconductor [package vignette of diffcyt](http://bioconductor.org/packages/release/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html). 
An example use of the `censcyt` methods for differential abundance analysis with a covariate subject to right censoring is 
available in the [package vignette](https://github.com/retogerber/censcyt/blob/main/vignettes/censored_covariate.Rmd)


## Development version

To install directly from GitHub run the following code:
```{r}
# First install 'devtools' package from CRAN
install.packages("devtools")

# Then install 'censcyt' package from GitHub
devtools::install_github("gerberreto/censcyt")
```

