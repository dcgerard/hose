<!-- README.md is generated from README.Rmd. Please edit that file -->
HOSE: Higher-Order Spectral Estimators
======================================

[![Build Status](https://travis-ci.org/dcgerard/hose.svg?branch=master)](https://travis-ci.org/dcgerard/hose) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dcgerard/hose?branch=master&svg=true)](https://ci.appveyor.com/project/dcgerard/hose) [![Coverage Status](https://img.shields.io/codecov/c/github/dcgerard/hose/master.svg)](https://codecov.io/github/dcgerard/hose?branch=master)

Summary
=======

`hose` is a package designed for working with higher-order spectral estimators. These estimators are based on the higher-order singular value decomposition of De Lathauwer et. al. (2000) and are useful when your data exhibit tensor-specific structure, such as having approximately low multilinear rank. This code will allow you to:

-   Calculate Stein's unbiased risk estimate (SURE) for all higher-order spectral estimators that are weakly differentiable and satisfy mild integrability conditions.
-   Calculate the mode-specific soft-thresholding estimator that minimizes the SURE using a coordinate descent algorithm.
-   Iterate through all the possible multilinear ranks of a mean tensor and choose the multilinear rank the minimizes the SURE.
-   Calculate a generalized SURE, motivated by generalized cross validation \[Sardy, 2012\], for all higher-order spectral estimators.
-   Calculate the SURE of estimators that apply higher-order spectral shrinkage to sub-tensors of the overall data tensor.
-   Calculate the SURE for estimators that individually shrink elements of the core array of the HOSVD of the data tensor.
-   Non-parametrically estimate the variance to use in these SURE procedures.

All details of these methods may be found in

Gerard, D., & Hoff, P. (2015). [Adaptive higher-order spectral estimators](http://arxiv.org/pdf/1505.02114v1.pdf). *arXiv preprint arXiv:1505.02114*.

A vignette demonstrating the methods available in `hose` may be found [here](http://home.uchicago.edu/~dcgerard/sure_example.html).

Installation
============

To install, run the following code in R:

``` r
install.packages(c("tensr", "softImpute", "RMTstat", "devtools"))
devtools::install_github("dcgerard/hose")
```

References
==========

**Gerard, D.**, & Hoff, P. (2015). [Adaptive higher-order spectral estimators](http://arxiv.org/pdf/1505.02114v1.pdf). *arXiv preprint arXiv:1505.02114*.

Lieven De Lathauwer, Bart De Moor, and Joos Vandewalle. [A multilinear singular value decomposition](http://epubs.siam.org/doi/abs/10.1137/S0895479896305696) . *SIAM J. Matrix Anal. Appl.*, 21(4):1253–1278 (electronic), 2000. ISSN 0895-4798. doi: 10.1137/S0895479896305696.

Sylvain Sardy. [Smooth blockwise iterative thresholding: a smooth fixed point estimator based on the likelihood’s block gradient](http://dx.doi.org/10.1080/01621459.2012.664527). *J. Amer. Statist. Assoc.*, 107(498):800–813, 2012. ISSN 0162-1459. doi: 10.1080/01621459.2012.664527.
