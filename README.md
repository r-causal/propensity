
<!-- README.md is generated from README.Rmd. Please edit that file -->

# propensity <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/r-causal/propensity/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-causal/propensity/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/r-causal/propensity/graph/badge.svg)](https://app.codecov.io/gh/r-causal/propensity)
<!-- badges: end -->

## Overview

propensity makes it easy to calculate propensity score weights and use
them to estimate causal effects. It supports:

- Six estimands for binary exposures (ATE, ATT, ATU, ATO, ATM, and
  entropy weights)
- Binary, categorical, and continuous exposures
- Trimming, truncation, and calibration for extreme propensity scores
- Inverse probability weighted estimation with standard errors that
  account for propensity score estimation

You can learn more in `vignette("propensity")`.

## Installation

You can install propensity from [CRAN](https://cran.r-project.org/)
with:

``` r
install.packages("propensity")
```

You can install the development version of propensity from
[GitHub](https://github.com/r-causal/propensity) with:

``` r
# install.packages("pak")
pak::pak("r-causal/propensity")
```

## Usage

``` r
library(propensity)

# Simulate data with a confounder, binary exposure, and binary outcome
n <- 200
x1 <- rnorm(n)
z <- rbinom(n, 1, plogis(0.5 * x1))
y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
dat <- data.frame(x1, z, y)

# Step 1: Fit a propensity score model
ps_mod <- glm(z ~ x1, data = dat, family = binomial())

# Step 2: Calculate ATE weights and fit a weighted outcome model
wts <- wt_ate(ps_mod)
outcome_mod <- glm(y ~ z, data = dat, family = binomial(), weights = wts)

# Step 3: Estimate causal effects with correct standard errors
ipw(ps_mod, outcome_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Propensity Score Model:
#>   Call: glm(formula = z ~ x1, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z, family = binomial(), data = dat, weights = wts) 
#> 
#> Estimates:
#>         estimate  std.err      z  ci.lower ci.upper conf.level p.value  
#> rd      0.142304 0.070204 2.0270 0.0047068  0.27990       0.95 0.04266 *
#> log(rr) 0.280314 0.142195 1.9713 0.0016172  0.55901       0.95 0.04869 *
#> log(or) 0.573392 0.286710 1.9999 0.0114518  1.13533       0.95 0.04551 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

By default, `ipw()` computes standard errors by M-estimation, stacking
the propensity score and outcome estimating equations so the uncertainty
of estimating the propensity scores is carried into the standard errors.
Set `se_method = "linearization"` for the influence-function method
(binary exposures with an exposure-only outcome model).

## Estimands

Each weight function targets a different population:

| Estimand    | Target population           | Function                       |
|-------------|-----------------------------|--------------------------------|
| **ATE**     | Entire population           | `wt_ate()`                     |
| **ATT**     | Treated units               | `wt_att()`                     |
| **ATU**     | Untreated units             | `wt_atu()` (alias: `wt_atc()`) |
| **ATO**     | Overlap population          | `wt_ato()`                     |
| **ATM**     | Matched population          | `wt_atm()`                     |
| **Entropy** | Entropy-balanced population | `wt_entropy()`                 |

ATO and ATM weights are bounded by construction, making them a good
alternative when ATE weights are highly variable.

## Learn more

- [Causal Inference in R](https://www.r-causal.org/) – A book on causal
  inference methods in R
- `vignette("propensity")` – Getting started with propensity score
  weighting
- [propensity package
  documentation](https://r-causal.github.io/propensity/) – Full
  reference and articles
