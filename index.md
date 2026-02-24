# propensity

## Overview

propensity makes it easy to calculate propensity score weights and use
them to estimate causal effects. It supports:

- Six estimands for binary exposures (ATE, ATT, ATU, ATO, ATM, and
  entropy weights)
- Binary, categorical, and continuous exposures
- Trimming, truncation, and calibration for extreme propensity scores
- Inverse probability weighted estimation with standard errors that
  account for propensity score estimation

You can learn more in
[`vignette("propensity")`](https://r-causal.github.io/propensity/articles/propensity.md).

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
#>         estimate std.err       z ci.lower ci.upper conf.level   p.value    
#> rd       0.14230 0.07038 2.02194   0.0044  0.28025       0.95 0.0431831 *  
#> log(rr)  0.28031 0.10770 2.60262   0.0692  0.49141       0.95 0.0092513 ** 
#> log(or)  0.57339 0.16200 3.53950   0.2559  0.89090       0.95 0.0004009 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

[`ipw()`](https://r-causal.github.io/propensity/reference/ipw.md) uses
linearization to account for uncertainty in the estimated propensity
scores when computing standard errors.

## Estimands

Each weight function targets a different population:

| Estimand    | Target population           | Function                                                                                                                                                 |
|-------------|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| **ATE**     | Entire population           | [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)                                                                                  |
| **ATT**     | Treated units               | [`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md)                                                                                  |
| **ATU**     | Untreated units             | [`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md) (alias: [`wt_atc()`](https://r-causal.github.io/propensity/reference/wt_ate.md)) |
| **ATO**     | Overlap population          | [`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md)                                                                                  |
| **ATM**     | Matched population          | [`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md)                                                                                  |
| **Entropy** | Entropy-balanced population | [`wt_entropy()`](https://r-causal.github.io/propensity/reference/wt_ate.md)                                                                              |

ATO and ATM weights are bounded by construction, making them a good
alternative when ATE weights are highly variable.

## Learn more

- [Causal Inference in R](https://www.r-causal.org/) – A book on causal
  inference methods in R
- [`vignette("propensity")`](https://r-causal.github.io/propensity/articles/propensity.md)
  – Getting started with propensity score weighting
- [propensity package
  documentation](https://r-causal.github.io/propensity/) – Full
  reference and articles
