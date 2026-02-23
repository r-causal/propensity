
<!-- README.md is generated from README.Rmd. Please edit that file -->

# propensity <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/r-causal/propensity/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-causal/propensity/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/r-causal/propensity/graph/badge.svg)](https://app.codecov.io/gh/r-causal/propensity)
<!-- badges: end -->

propensity provides a toolkit for propensity score analysis in causal
inference. It supports multiple causal estimands across binary,
categorical, and continuous exposures, handles extreme propensity scores
through trimming, truncation, and calibration, and estimates causal
effects with valid standard errors via inverse probability weighting.

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

## Quick Start

Estimate a causal effect in three steps:

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

`ipw()` performs inverse probability weighted estimation. It accounts
for uncertainty in the estimated propensity scores when computing
standard errors, producing valid confidence intervals and p-values.

## Multiple Estimands

propensity supports six causal estimands for binary exposures:

``` r
ps <- fitted(ps_mod)

wt_ate(ps, z)     # Average Treatment Effect
wt_att(ps, z)     # Average Treatment Effect on the Treated
wt_atu(ps, z)     # Average Treatment Effect on the Untreated
wt_ato(ps, z)     # Average Treatment Effect, Overlap population
wt_atm(ps, z)     # Average Treatment Effect, Matched population
wt_entropy(ps, z) # Entropy-weighted ATE
```

Choose your estimand based on your research question:

| Estimand | Target population | Function |
|----|----|----|
| **ATE** | Entire population | `wt_ate()` |
| **ATT** | Treated units | `wt_att()` |
| **ATU** | Untreated units | `wt_atu()` (alias: `wt_atc()`) |
| **ATO** | Overlap population – units with the most equipoise between groups | `wt_ato()` |
| **ATM** | Matched population – mimics 1:1 matching | `wt_atm()` |
| **Entropy** | Entropy-balanced compromise between ATE and overlap | `wt_entropy()` |

ATO and ATM produce naturally bounded weights, making them good
alternatives when ATE weights are highly variable.

## Flexible Input

Weight functions accept propensity scores in several forms:

``` r
# GLM object -- extracts fitted values and exposure automatically (recommended)
wt_ate(ps_mod)

# Numeric vector of propensity scores
wt_ate(ps, z)

# Data frame of class probabilities (e.g., from multinomial models)
ps_df <- data.frame(control = 1 - ps, treated = ps)
wt_ate(ps_df, z)
```

## Handling Extreme Weights

Propensity scores near 0 or 1 produce extreme weights that inflate
variance. propensity offers four complementary strategies:

1.  **Switch estimand** – `wt_ato()` and `wt_atm()` produce bounded
    weights by design.
2.  **Trim** – remove observations with extreme scores (sets them to
    `NA`).
3.  **Truncate** – winsorize extreme scores to a fixed range.
4.  **Calibrate** – adjust scores so they better reflect true treatment
    probabilities.

``` r
ps <- fitted(ps_mod)

# Diagnose: inspect the weight distribution
summary(wt_ate(ps, z))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.380   1.806   1.956   2.000   2.154   2.902

# 1. Switch estimand -- overlap weights are bounded by design
summary(wt_ato(ps, z))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2752  0.4463  0.4889  0.4912  0.5358  0.6554

# 2. Trim -- remove observations with extreme propensity scores
ps_trimmed <- ps_trim(ps, method = "adaptive")
summary(wt_ate(ps_trimmed, z))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.380   1.806   1.956   2.000   2.154   2.902

# 3. Truncate -- bound extreme propensity scores to a range
ps_truncated <- ps_trunc(ps, lower = 0.05, upper = 0.95)
summary(wt_ate(ps_truncated, z))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.380   1.806   1.956   2.000   2.154   2.902

# 4. Calibrate -- adjust scores to better reflect true probabilities
ps_calibrated <- ps_calibrate(ps, z)
summary(wt_ate(ps_calibrated, z))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.401   1.804   1.956   2.000   2.156   2.881
```

After trimming, refit the propensity score model on the retained subset
so the scores reflect the trimmed population:

``` r
ps_refitted <- ps_refit(ps_trimmed, ps_mod)
wts_refitted <- wt_ate(ps_refitted, z)
```

## Advanced Features

### Weight Stabilization

Stabilized weights can reduce variance. Stabilization is supported for
`wt_ate()` and `wt_cens()`:

``` r
summary(wt_ate(ps, z, stabilize = TRUE))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.7036  0.9017  0.9817  0.9999  1.0732  1.4536
```

### Continuous Exposures

For continuous treatments, weights use the ratio of the marginal to
conditional density. Stabilization is strongly recommended:

``` r
# Fit a model for the continuous exposure
continuous_exposure <- rnorm(n, mean = 0.5 * x1)
dat$a <- continuous_exposure
ps_continuous <- glm(a ~ x1, data = dat, family = gaussian())

# Density-ratio weights (stabilization strongly recommended)
wts_continuous <- wt_ate(ps_continuous, stabilize = TRUE)
```

### Categorical Exposures

For multi-level treatments, supply a matrix or data frame of class
probabilities with one column per treatment level:

``` r
# Multinomial propensity scores (one column per treatment level)
ps_matrix <- matrix(c(0.3, 0.5, 0.2), ncol = 3, nrow = n, byrow = TRUE)
categorical_exposure <- factor(sample(1:3, n, replace = TRUE))

wt_ate(ps_matrix, categorical_exposure, exposure_type = "categorical")

# For ATT with categorical exposures, specify the focal level
wt_att(ps_matrix, categorical_exposure, .focal_level = "2")
```

### Censoring Weights

`wt_cens()` calculates inverse probability of censoring weights for
survival or longitudinal analyses. These address informative censoring –
not treatment assignment:

``` r
# Model the probability of being uncensored
cens_mod <- glm(uncensored ~ x1 + x2, data = dat, family = binomial())

# Censoring weights (uses the same formula as wt_ate())
wts_cens <- wt_cens(cens_mod)

# Combine with treatment weights for a doubly-weighted analysis
wts_combined <- wt_ate(ps_mod) * wts_cens
```

### Calibration Methods

`ps_calibrate()` supports two calibration methods:

``` r
# Logistic calibration (default) -- fits a logistic regression of exposure on
# predicted propensity scores
ps_calibrate(ps, z, method = "logistic")

# Isotonic regression calibration -- fits a monotone step function; useful for
# non-smooth relationships with large samples
ps_calibrate(ps, z, method = "isoreg")
```

## Learn More

- [propensity package
  documentation](https://r-causal.github.io/propensity/) – Full
  reference and articles
- [Causal Inference in R](https://www.r-causal.org/) – A comprehensive
  guide to causal inference methods in R
