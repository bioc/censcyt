# censcyt

[![Build Status](https://travis-ci.org/retogerber/censcyt.svg?branch=main)](https://travis-ci.org/retogerber/censcyt)
[![codecov](https://codecov.io/gh/retogerber/censcyt/branch/main/graph/badge.svg)](https://codecov.io/gh/retogerber/censcyt)

## Summary

`censcyt` is an R package extending the [diffcyt](https://github.com/lmweber/diffcyt) (differential discovery in high-dimensional cytometry via high-resolution clustering) pipeline. `censcyt` (**Cens**ored diff**cyt**) includes methods for differential abundance analysis in the presence of a covariate subject to right censoring. It uses the *reversed* association testing approach (like `diffcyt`) meaning the censored variable (e.g. survival time) is modeled as a predictor. Classical survival analysis methods on the other hand model the censored variable as the response. See also our [preprint on biorxiv](https://www.biorxiv.org/content/10.1101/2020.11.09.374447v1).


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

