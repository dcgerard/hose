
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hose: Higher-Order Spectral Estimators

[![Build
Status](https://travis-ci.org/dcgerard/hose.svg?branch=master)](https://travis-ci.org/dcgerard/hose)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/dcgerard/hose?branch=master&svg=true)](https://ci.appveyor.com/project/dcgerard/hose)
[![CRAN
status](https://www.r-pkg.org/badges/version/hose)](https://CRAN.R-project.org/package=hose)

# Summary

`hose` is a package designed for working with higher-order spectral
estimators, which were first introduced in Gerard and Hoff (2017). These
estimators are based on the higher-order singular value decomposition of
De Lathauwer, De Moor, and Vandewalle (2000) and are useful when your
data exhibit tensor-specific structure, such as having approximately low
multilinear rank. This code will allow you to:

  - Calculate Stein’s unbiased risk estimate (SURE) for all higher-order
    spectral estimators that are weakly differentiable and satisfy mild
    integrability conditions.
  - Calculate the mode-specific soft-thresholding estimator that
    minimizes the SURE using a coordinate descent algorithm.
  - Iterate through all the possible multilinear ranks of a mean tensor
    and choose the multilinear rank the minimizes the SURE.
  - Calculate a generalized SURE, motivated by generalized cross
    validation (Josse and Sardy 2016), for all higher-order spectral
    estimators.
  - Calculate the SURE of estimators that apply higher-order spectral
    shrinkage to sub-tensors of the overall data tensor.
  - Calculate the SURE for estimators that individually shrink elements
    of the core array of the HOSVD of the data tensor.
  - Non-parametrically estimate the variance to use in these SURE
    procedures.

The main functions are:

  - `get_c()`: Pre-format the data before applying mode-specific
    singular value shrinkage.
  - `tensor_var_est()`: Estimate the variance of the data from multiple
    options.
  - `soft_coord()`: Estimate the underlying low-rank mean tensor via
    soft-thresholding.

# Citation

If you find these methods useful, please cite

> **Gerard, David**, and Peter Hoff. 2017. “Adaptive Higher-Order
> Spectral Estimators.” *Electron. J. Statist.* 11 (2). The Institute of
> Mathematical Statistics; the Bernoulli Society: 3703–37.
> <https://doi.org/10.1214/17-EJS1330>.

Or, using BibTex:

``` tex
@ARTICLE{gerard2017adaptive,
    AUTHOR = {David Gerard and Peter Hoff},
     TITLE = {Adaptive higher-order spectral estimators},
   JOURNAL = {Electron. J. Statist.},
  FJOURNAL = {Electronic Journal of Statistics},
      YEAR = {2017},
    VOLUME = {11},
    NUMBER = {2},
     PAGES = {3703-3737},
      ISSN = {1935-7524},
       DOI = {10.1214/17-EJS1330},
      SICI = {1935-7524(2017)11:2<3703:AHOSE>2.0.CO;2-Q},
}
```

# Installation

You can install from CRAN in the usual way:

``` r
install.packages("hose")
```

Or, to install the latest (unstable) version, run the following code in
R:

``` r
install.packages(c("tensr", "softImpute", "RMTstat", "devtools"))
devtools::install_github("dcgerard/hose")
```

# Vignette

I’ve provided a vignette demonstrating the methods available in `hose`.
You can find it
[here](http://dcgerard.github.io/code/sure_example.html). Or you can
build the vignette on install with

``` r
install.packages("devtools")
devtools::install_github("dcgerard/hose", build_vignettes = TRUE)
```

and access the vignette by running the following code in R:

``` r
utils::vignette("sure_example", package = "hose")
```

# References

<div id="refs" class="references">

<div id="ref-lathauwer2000multilinear">

De Lathauwer, L., B. De Moor, and J. Vandewalle. 2000. “A Multilinear
Singular Value Decomposition.” *SIAM Journal on Matrix Analysis and
Applications* 21 (4): 1253–78.
<https://doi.org/10.1137/S0895479896305696>.

</div>

<div id="ref-gerard2017adaptive">

Gerard, David, and Peter Hoff. 2017. “Adaptive Higher-Order Spectral
Estimators.” *Electron. J. Statist.* 11 (2): 3703–37.
<https://doi.org/10.1214/17-EJS1330>.

</div>

<div id="ref-josse2016adaptive">

Josse, Julie, and Sylvain Sardy. 2016. “Adaptive Shrinkage of Singular
Values.” *Statistics and Computing* 26 (3): 715–24.
<https://doi.org/10.1007/s11222-015-9554-9>.

</div>

</div>
