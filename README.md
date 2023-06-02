
<!-- README.md is generated from README.Rmd. Please edit that file -->

# propensity <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->
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

Currently, propensity supports calculating propensity score weights for
binary exposures:

``` r
library(propensity)

propensity_scores <- c(.1, .3, .4, .3)
x <- c(0, 0, 1, 0)

# Average treatment effect (ATE) weights
wt_ate(propensity_scores, .exposure = x)
#> [1] 1.111111 1.428571 2.500000 1.428571

# Stabilized ATE weights
wt_ate(propensity_scores, .exposure = x, stabilize = TRUE)
#> [1] 0.2777778 0.3571429 0.6250000 0.3571429

# Average treatment effect in the overlap (ATO) weights
wt_ato(propensity_scores, .exposure = x)
#> [1] 0.1 0.3 0.6 0.3
```
