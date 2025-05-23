
<!-- README.md is generated from README.Rmd. Please edit that file -->

# propensity <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/malcolmbarrett/propensity/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/malcolmbarrett/propensity/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/malcolmbarrett/propensity/graph/badge.svg)](https://app.codecov.io/gh/malcolmbarrett/propensity)
<!-- badges: end -->

The goal of propensity is to calculate propensity scores and weights for
a wide variety of research questions.

propensity is under very early development.

## Installation

You can install the development version of propensity from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("r-causal/propensity")
```

## Example

propensity is under very early development. Currently, it supports
calculating propensity score weights for binary exposures:

``` r
library(propensity)
propensity_scores <- c(.1, .3, .4, .3)
x <- c(0, 0, 1, 0)

# ATE weights
wt_ate(propensity_scores, .exposure = x)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[4]>
#> [1] 1.111111 1.428571 2.500000 1.428571

# Stabilized ATE weights
wt_ate(propensity_scores, .exposure = x, stabilize = TRUE)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[4]>
#> [1] 0.2777778 0.3571429 0.6250000 0.3571429

# ATO weights
wt_ato(propensity_scores, .exposure = x)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ato}[4]>
#> [1] 0.1 0.3 0.6 0.3
```
