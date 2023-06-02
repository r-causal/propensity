
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
binary exposures.

There are [5 common types of
weights](https://www.r-causal.org/chapters/chapter-10.html) that you
might wish to use in a propensity score setting:

1.  **Average treatment effect (ATE):** used when the target population
    is the entire population of interest. An example question this
    answers is “Should a marketing campaign be rolled out to all
    eligible people?”

2.  **Average treatment effect among the treated (ATT):** used when the
    target population is the exposed/treated population. An example
    question this answers is “Should we stop our marketing campaign to
    those currently receiving it?”

3.  **Average treatment effect among the unexposed (ATU):** used when
    the target population is the unexposed/untreated/control population.
    An example question this answers is “Should we send our marketing
    campaign to those not currently receiving it?”

4.  **Average treatment effect among the evenly matchable (ATM):** used
    when the target population is those deemed “evenly matchable” by
    some distance metric. An example question this answers is “Should we
    send our marketing campaign to those of with similar demographic
    characteristics?”

5.  **Average treatment effect among the overlap (ATO):** used when the
    target population is the same as the ATM setting, however the ATO
    weights are slightly attenuated with improved variance properties.

Each of these weights can optionally be “stabilized” to prevent extreme
weights that lead to wide confidence intervals.

The below example shows how we would generate ATE and ATO weights for 4
participants (1 exposed and 3 unexposed) with known probabilities of
being exposed (propensity scores).

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
