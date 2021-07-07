# censcyt
[![Bioc devel](http://bioconductor.org/shields/build/devel/bioc/censcyt.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/censcyt)
[![Bioc release](http://bioconductor.org/shields/build/release/bioc/censcyt.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/censcyt)
[![R build status](https://github.com/retogerber/censcyt/workflows/R-CMD-check/badge.svg)](https://github.com/retogerber/censcyt/actions)
[![codecov](https://codecov.io/gh/retogerber/censcyt/branch/main/graph/badge.svg)](https://codecov.io/gh/retogerber/censcyt)

## Summary

`censcyt` is an R package extending the [diffcyt](https://github.com/lmweber/diffcyt) (differential discovery in high-dimensional cytometry via high-resolution clustering) pipeline. `censcyt` (**Cens**ored diff**cyt**) includes methods for differential abundance analysis in the presence of a covariate subject to right censoring. It uses the *reversed* association testing approach (like `diffcyt`) meaning the censored variable (e.g. survival time) is modeled as a predictor. Classical survival analysis methods on the other hand model the censored variable as the response. See also our [preprint on biorxiv](https://www.biorxiv.org/content/10.1101/2020.11.09.374447v1).


## Vignettes

The main workflow can be found in the Bioconductor [package vignette of diffcyt](http://bioconductor.org/packages/release/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html). 
An example use of the `censcyt` methods for differential abundance analysis with a covariate subject to right censoring is 
available in the [package vignette](http://bioconductor.org/packages/devel/bioc/vignettes/censcyt/inst/doc/censored_covariate.html) on bioconductor.


## Install

To install from bioconductor:
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("censcyt")
```


### Development version

To install directly from GitHub run the following code:
```{r}
# First install 'devtools' package from CRAN
install.packages("devtools")

# Then install 'censcyt' package from GitHub
devtools::install_github("retogerber/censcyt")
```

