# propensity 0.1.0.9000 (development version)

* `ipw()` now defaults to M-estimation sandwich standard errors, computed by
  stacking the propensity score, outcome, and estimand estimating equations
  with the deli package. Point estimates are unchanged for the `ate` estimand
  and for saturated (exposure-only) outcome models; the tilted estimands with a
  covariate-adjusted outcome model are corrected as described below. Set
  `se_method = "linearization"` to restore the previous influence-function
  method (binary exposures only).

* Corrected the marginal-mean standardization for the tilted estimands (`att`,
  `atu`, `atm`, `ato`, and `entropy`) on the M-estimation path. With a
  covariate-adjusted outcome model, the stacked marginal means previously
  averaged the counterfactual predictions over the full sample, reporting an
  ATE-type contrast for every estimand. The marginal means are now standardized
  to the estimand's tilted target population, so point estimates and standard
  errors change for those configurations. The `ate` estimand and saturated
  outcome models are unaffected.

* `se_method = "linearization"` now requires an exposure-only, offset-free
  outcome model and errors otherwise, directing you to
  `se_method = "mestimation"`. The linearization influence functions are those
  of the Hajek weighted-mean estimator, which match the reported g-computation
  estimates only for an outcome model of the exposure alone; the M-estimation
  path handles covariate-adjusted outcome models correctly.

* Corrected the linearization standard errors for `log(rr)` and `log(or)`,
  which were scaled by the reciprocal of the risk ratio and odds ratio,
  respectively (underestimated when the ratio exceeds one, overestimated
  otherwise). Risk-difference standard errors are unaffected.

* Corrected the linearization standard errors for probit and cloglog propensity
  score models. The influence function for the propensity score coefficients
  omitted the GLM score factor `1 / (p (1 - p) g'(p))`, which is 1 only for the
  canonical logit link, so probit and cloglog standard errors were mis-scaled.
  Logit standard errors are unchanged.

* `ipw()` now supports categorical exposures through a `nnet::multinom()`
  propensity score model and continuous exposures through an `lm()` or
  gaussian-family `glm()` propensity score model with a weighted marginal
  structural outcome model. A binary or categorical exposure supports every
  estimand; a continuous exposure supports the ATE.

* `ipw()` now supports the `atu` and `entropy` estimands for binary and
  categorical exposures, alongside `ate`, `att`, `atm`, and `ato`.

* `ipw()` gained guards that error, with guidance to refit, when the outcome
  model weights were trimmed, truncated, or calibrated; are missing; or do not
  match the propensity score model. It also rejects a propensity score model fit
  with case weights and an outcome model fit with an offset term.

* `ipw()` now rejects an outcome model whose family or link it cannot stack,
  on both standard error methods, with an informative error instead of silently
  treating it as a binomial or identity-linear model. A poisson, quasipoisson,
  Gamma, or inverse-gaussian outcome model, or a gaussian model with a
  non-identity link, errors and directs you to a binomial or quasibinomial model
  or a gaussian identity link.

* For a continuous exposure, `ipw()` now validates the propensity score and
  marginal structural model links at entry with honest errors. A non-identity
  gaussian propensity score link (previously surfaced as a misleading
  weights-mismatch error) and a marginal structural model link outside identity,
  logit, and log (previously an error only after fitting) now fail fast, naming
  the offending link.

* `ipw()` now converts a factor or logical outcome response to the model's 0/1
  coding on both standard error methods, following glm's convention (the first
  factor level is failure, every other level is success). Previously a factor
  response crashed the M-estimation solve or produced `NA` linearization
  standard errors with factor-arithmetic warnings.

* Factor and logical exposures now work on the linearization path, recoded to
  0/1 with the second factor level (or `TRUE`) as the exposed group, matching
  the M-estimation path. Previously a factor exposure crashed the linearization
  variance with factor-arithmetic errors.

* The binary M-estimation path now reconstructs the propensity score design from
  `.data` when the model frame is unavailable (for example, a model fit with
  `model = FALSE` whose fitting data is gone) and reports an informative error
  directing you to supply `.data`, matching the categorical and continuous
  paths. Previously it failed with a raw "object not found" error that `.data`
  did not resolve.

* `wt_ate()` and `wt_cens()` now record a user-supplied `stabilization_score`
  attribute on the returned weights, readable with the new `stabilization_score()`
  accessor.

* Fixed `broom::tidy(glm_fit, conf.int = TRUE)` failing on GLMs weighted by
  `psw` vectors. `confint.glm()` builds profile-likelihood intervals via
  `profile.glm()`, which refits through `glm.fit()`; the refit indexes
  `weights[good]` with a matrix subscript, which `[.vctrs_vctr` rejected.
  Added a `[.psw` method that falls back to base R linear indexing for
  matrix/array subscripts and delegates everything else to `[.vctrs_vctr`.

* Comparison operators on `psw` (`==`, `!=`, `<`, `>`, `<=`, `>=`) now
  short-circuit `vec_equal()` / `vec_compare()` and return a logical vector
  silently. Previously each comparison fired a `propensity_class_downgrade`
  warning via `vec_ptype2.psw.double()`, producing 100+ warnings during a
  single `tidy(glm, conf.int = TRUE)` call. Combine and cast paths still warn.

* Added a `NEWS.md` file to track changes to the package.
