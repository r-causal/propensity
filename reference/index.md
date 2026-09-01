# Package index

## Overview

A tour of the package, from weights to effect estimates.

- [`propensity`](https://r-causal.github.io/propensity/reference/propensity-package.md)
  [`propensity-package`](https://r-causal.github.io/propensity/reference/propensity-package.md)
  : propensity: A Toolkit for Calculating and Working with Propensity
  Scores

## Weights

Propensity score weights for each estimand, across binary, categorical,
and continuous exposures.

- [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  [`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  [`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  [`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  [`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  [`wt_entropy()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  [`wt_atc()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  [`wt_cens()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  : Calculate propensity score weights

## Weight objects

Construct propensity score weights and read what they record.

- [`new_psw()`](https://r-causal.github.io/propensity/reference/psw.md)
  [`psw()`](https://r-causal.github.io/propensity/reference/psw.md)
  [`is_psw()`](https://r-causal.github.io/propensity/reference/psw.md)
  [`is_stabilized()`](https://r-causal.github.io/propensity/reference/psw.md)
  [`stabilization_score()`](https://r-causal.github.io/propensity/reference/psw.md)
  [`as_psw()`](https://r-causal.github.io/propensity/reference/psw.md) :
  Propensity Score Weight Vectors
- [`exposure_type()`](https://r-causal.github.io/propensity/reference/exposure_type.md)
  [`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md)
  [`format(`*`<propensity_density_meta>`*`)`](https://r-causal.github.io/propensity/reference/exposure_type.md)
  [`print(`*`<propensity_density_meta>`*`)`](https://r-causal.github.io/propensity/reference/exposure_type.md)
  : What a set of weights records about the exposure
- [`numerator_model()`](https://r-causal.github.io/propensity/reference/numerator_model.md)
  : The model a set of weights was stabilized on

## Joint treatments

Weights for a sequence of two treatments, from a factorized pair of
treatment models.

- [`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
  [`is_joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
  : Record the two treatment models of a joint exposure
- [`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
  [`is_joint_wt()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
  [`joint_wt_meta()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
  : Product weights for a joint intervention on two treatments

## Continuous exposure densities

The density families a continuous exposure’s weights are built from.

- [`dens_normal()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
  [`dens_laplace()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
  [`dens_t()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
  [`dens_kernel()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
  [`dens_fn()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
  : Density specifications for continuous exposures

## Trimming

Restrict a set of propensity scores to a region of overlap, and refit
the model on what remains.

- [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
  : Trim Propensity Scores

- [`ps_trim_meta()`](https://r-causal.github.io/propensity/reference/ps_trim_meta.md)
  :

  Extract trimming metadata from a `ps_trim` object

- [`is_ps_trimmed()`](https://r-causal.github.io/propensity/reference/is_ps_trimmed.md)
  : Test whether propensity scores have been trimmed

- [`is_unit_trimmed()`](https://r-causal.github.io/propensity/reference/is_unit_trimmed.md)
  : Identify which units were trimmed

- [`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md)
  : Refit a Propensity Score Model on Retained Observations

- [`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md)
  : Check if propensity scores have been refit

## Truncation

Bound extreme propensity scores in place rather than dropping them.

- [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
  : Truncate (Winsorize) Propensity Scores

- [`ps_trunc_meta()`](https://r-causal.github.io/propensity/reference/ps_trunc_meta.md)
  :

  Extract truncation metadata from a `ps_trunc` object

- [`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md)
  : Test whether propensity scores have been truncated

- [`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md)
  : Identify which units were truncated

## Calibration and tilting

Recalibrate fitted propensity scores, and evaluate the tilting function
an estimand implies.

- [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
  : Calibrate propensity scores
- [`is_ps_calibrated()`](https://r-causal.github.io/propensity/reference/is_ps_calibrated.md)
  : Check if propensity scores are calibrated
- [`ps_tilt()`](https://r-causal.github.io/propensity/reference/ps_tilt.md)
  : Propensity score tilting functions

## Effect estimation

Inverse probability weighted estimation and the surfaces its results are
read through.

- [`ipw(`*`<joint_wt_models>`*`)`](https://r-causal.github.io/propensity/reference/ipw-methods.md)
  [`ipw(`*`<multinom>`*`)`](https://r-causal.github.io/propensity/reference/ipw-methods.md)
  [`ipw(`*`<lm>`*`)`](https://r-causal.github.io/propensity/reference/ipw-methods.md)
  [`ipw(`*`<glm>`*`)`](https://r-causal.github.io/propensity/reference/ipw-methods.md)
  : Inverse Probability Weighted Estimation
- [`print(`*`<ipw_diagnostic_se>`*`)`](https://r-causal.github.io/propensity/reference/print.ipw_diagnostic_se.md)
  [`as.data.frame(`*`<ipw_diagnostic_se>`*`)`](https://r-causal.github.io/propensity/reference/print.ipw_diagnostic_se.md)
  [`tidy(`*`<ipw_diagnostic_se>`*`)`](https://r-causal.github.io/propensity/reference/print.ipw_diagnostic_se.md)
  : Print a result whose standard errors are a diagnostic
- [`tidy(`*`<ipw>`*`)`](https://r-causal.github.io/propensity/reference/tidy.ipw.md)
  : Tidy an inverse probability weighted result
- [`glance(`*`<ipw>`*`)`](https://r-causal.github.io/propensity/reference/glance.ipw.md)
  : Glance at an inverse probability weighted result
- [`augment(`*`<ipw>`*`)`](https://r-causal.github.io/propensity/reference/augment.ipw.md)
  : Augment an inverse probability weighted result with per-observation
  columns
- [`tidy(`*`<ipw_pooled>`*`)`](https://r-causal.github.io/propensity/reference/tidy.ipw_pooled.md)
  : Tidy a pooled inverse probability weighted result
- [`glance(`*`<ipw_pooled>`*`)`](https://r-causal.github.io/propensity/reference/glance.ipw_pooled.md)
  : Glance at a pooled inverse probability weighted result

## Re-exports

- [`reexports`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`ipw`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`as_marginal`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`as_conditional`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`is_causal_wt`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`estimand`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`estimand<-`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`tidy`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`glance`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`augment`](https://r-causal.github.io/propensity/reference/reexports.md)
  [`pool_ipw`](https://r-causal.github.io/propensity/reference/reexports.md)
  : Objects exported from other packages
