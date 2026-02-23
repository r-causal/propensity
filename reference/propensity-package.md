# propensity: A Toolkit for Calculating and Working with Propensity Scores

propensity provides tools for propensity score analysis in causal
inference. Calculate propensity score weights for a variety of causal
estimands, handle extreme propensity scores through trimming,
truncation, and calibration, and estimate causal effects with inverse
probability weighting. The package supports binary, categorical, and
continuous exposures.

## Weight functions

Calculate propensity score weights for different causal estimands:

- [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md):
  Average treatment effect (ATE) weights

- [`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md):
  Average treatment effect on the treated (ATT) weights

- [`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md):
  Average treatment effect on the untreated (ATU) weights
  ([`wt_atc()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  is an alias)

- [`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md):
  Average treatment effect for the evenly matchable (ATM) weights

- [`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md):
  Average treatment effect for the overlap population (ATO) weights

- [`wt_entropy()`](https://r-causal.github.io/propensity/reference/wt_ate.md):
  Entropy balancing weights

- [`wt_cens()`](https://r-causal.github.io/propensity/reference/wt_ate.md):
  Censoring weights

## Propensity score modifications

Handle extreme propensity scores before calculating weights:

- [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md):
  Trim observations with extreme propensity scores

- [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md):
  Truncate (winsorize) extreme propensity scores

- [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md):
  Calibrate propensity scores to improve balance

- [`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md):
  Re-estimate the propensity score model after trimming

## Estimation

- [`ipw()`](https://r-causal.github.io/propensity/reference/ipw.md):
  Inverse probability weighted estimator with variance estimation that
  accounts for propensity score estimation uncertainty

## PSW class

The [`psw()`](https://r-causal.github.io/propensity/reference/psw.md)
class represents propensity score weights with metadata about the
estimand and modifications applied:

- [`psw()`](https://r-causal.github.io/propensity/reference/psw.md),
  [`as_psw()`](https://r-causal.github.io/propensity/reference/psw.md),
  [`is_psw()`](https://r-causal.github.io/propensity/reference/psw.md):
  Create and test propensity score weights

- [`estimand()`](https://r-causal.github.io/propensity/reference/psw.md):
  Query the causal estimand

- [`is_stabilized()`](https://r-causal.github.io/propensity/reference/psw.md):
  Check if weights are stabilized

## See also

- [`vignette("propensity")`](https://r-causal.github.io/propensity/articles/propensity.md)
  for a getting started guide

- The [package website](https://r-causal.github.io/propensity/) for full
  documentation

## Author

**Maintainer**: Malcolm Barrett <malcolmbarrett@gmail.com>
([ORCID](https://orcid.org/0000-0003-0299-5825)) \[copyright holder\]
