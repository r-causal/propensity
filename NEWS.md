# propensity 0.1.0.9000 (development version)

* Arithmetic and subsetting on a `psw` vector now carry the record left by a
  modified propensity score (`ps_trim_meta`, `ps_trunc_meta`, and
  `ps_calib_meta`) whenever the result comes back at the length the record was
  written for. Every such operation previously discarded these records, so
  `is_unit_trimmed()` on the result of something as ordinary as
  `weights / sum(weights)` reported every unit as retained while
  `is_ps_trimmed()` still reported the weights as trimmed. An operation that
  shortens the weights drops the record, silently: these records also travel by
  routes vctrs does not control, so a warning on the one route it does control
  would be neither complete nor about anything the user wrote.

* `is_unit_trimmed()` and `is_refit()` on a `psw` vector now raise an error of
  class `propensity_missing_meta_error`, pointing at the `ps_trim` object the
  weights were built from, instead of answering from a trimming record they
  cannot use. `is_unit_trimmed()` answers by position, so it refuses both when
  weights marked as trimmed carry no record and when the record does not cover
  the vector it is given. The second case is reachable with no subsetting at all
  in the user's code: `model.frame()` drops the `NA`-weighted rows from a
  weights column in C and re-attaches the original variable's attributes to the
  shortened result, so the weights read back out of an outcome model fit on
  trimmed weights carry a record written for rows that are no longer there, and
  previously reported trimmed units at stale positions among rows that had all
  been retained. `is_refit()` reads a single flag rather than a position, so it
  answers from any record present and refuses only when the record is absent
  entirely. Weights that were never trimmed are unaffected.

* `is_unit_trimmed()` on an empty `psw` now returns an empty logical vector.
  Where the weights were empty but still carried a trimming record, indexing by
  the positions that record named grew the result instead, returning a vector as
  long as the original weights, padded with `NA`, for a vector of no
  observations.

* The attributes describing a categorical exposure (`n_categories`,
  `category_names`, and `focal_category`) now survive `psw` arithmetic and
  subsetting, including subsetting that shortens the weights. They name the
  exposure levels rather than the units, so they mean the same thing at any
  length. Previously any arithmetic or subsetting dropped them, which left
  `ipw()` unable to resolve the focal level for categorical `"att"` and `"atu"`
  weights that had passed through a data frame operation.

* Casting a double vector, an integer vector, a `ps_trim`, or a `ps_trunc` to a
  `psw` prototype now carries every metadata field of the prototype rather than
  its estimand alone. The stabilized, trimmed, truncated, and calibrated flags
  and the stabilization score were dropped, so such a cast returned weights that
  described themselves as neither stabilized nor modified whatever the prototype
  said, with nothing signaled. A per-observation stabilization score transfers
  only when its length matches the data being cast, since it describes the
  prototype's observations rather than the incoming data; at any other length it
  is dropped silently and the result stays marked as stabilized.

* `ipw(se_method = "linearization")` now rejects a `ps_link` naming a link other
  than the one the propensity score model was fit with, naming both. The score
  factor, the weight derivatives, and the correction matrix are all derived from
  `ps_link`, so a mismatched value scaled the estimation correction by the wrong
  derivative and changed the standard errors while leaving the estimates
  untouched, with nothing signaled. `ps_link` left at its `NULL` default, or set
  to the model's own link, is unaffected.

* `ipw()` now checks that `.data` rebuilds the design its propensity score model
  was fit to, reporting a disagreement as a data error naming both widths and
  the usual cause, a column whose type differs from the fitting data. Such a
  `.data` previously produced a raw error about a names attribute for a
  categorical exposure, and would have produced a non-conformable multiply
  elsewhere. A covariate recoded without changing the design's width is
  unaffected, since it carries the same numbers.

* `ipw()` now reports a `.data` missing any covariate either model was fit with
  as a missing column, alongside the exposure and outcome it already checked.
  Previously only the exposure and outcome were checked and a missing covariate
  surfaced as a raw object-not-found error from `model.matrix()`.

* `ipw()` now rejects an outcome model with a matrix response, such as
  `cbind(successes, failures)`, with an error naming the response shape, on
  every exposure type and both standard error methods. Such a model previously
  failed further downstream and differently depending on the call: without
  `.data` it reported an outcome twice the length of the exposure, and with
  `.data` it reported a missing column named `"cbind"`, which cannot exist. The
  matrix-response guards for the propensity score model and the outcome model
  now share the error class `propensity_ipw_response_error`;
  `propensity_ipw_exposure_error` is reserved for the guards concerned with what
  the exposure means rather than what shape a response has.

* `ipw()` now reports an outcome model whose fitting data can no longer be
  reached with a classed error naming `outcome_mod`, on every exposure type and
  both standard error methods, instead of the raw object-not-found error the
  refitting attempt produced. The message directs the user to refit
  `outcome_mod`, and says why `.data` is not the remedy it is for a propensity
  model in the same state: the weights are read from the outcome model's own
  frame, which `.data` cannot supply.

* `ipw(se_method = "mestimation")` now rejects a propensity score model that
  separates the exposure, with an error that says so and names how many
  observations are affected. A separating model has no finite maximum likelihood
  estimate, and rebuilding its propensity scores for the estimating equations
  gives probabilities of exactly zero or one, which leave the corresponding
  weights undefined. Previously this surfaced as an error about the outcome
  model's weights and the focal level for the `"ate"` estimand, and as an error
  from the solver for the others, neither of which pointed at the propensity
  model. `se_method = "linearization"` is unaffected: it takes its propensity
  scores from `predict()`, whose inverse link cannot return exactly zero or one.

* `ipw(se_method = "linearization")` now checks the outcome model's weights
  against the propensity score model, as the M-estimation path already did. The
  linearization path predicts the propensity scores and takes the weights from
  the outcome model without requiring the two to describe the same analysis, so
  a `.data` whose covariate values differ from the ones the weights were built
  on, or weights that lost their per-observation stabilization score in an
  operation that changed their length, changed the standard errors with nothing
  signaled. The point estimates were unaffected in both cases, which made the
  output look untouched. Both now raise the same error the M-estimation path
  raises.

* `ipw()` now reports a `.data` whose row count disagrees with the fitted models
  as a problem with `.data`, on every exposure type and both standard error
  methods. The linearization path already did. The M-estimation paths sized
  their designs to `.data` instead, so the disagreement surfaced later as a
  weight-consistency failure, whose message asks about how the weights were
  built and, for a binary or categorical exposure, about the focal level. None
  of that is the cause when the mistake is the wrong data frame. A `.data` of
  the right size but the wrong content is still reported as a weights
  inconsistency, since a row count cannot detect it.

* `ipw()` now rebuilds its design matrices under the contrasts the supplied
  models were fit with. The counterfactual and propensity score designs were
  previously rebuilt under whatever coding was in force at the time of the call,
  which for a model fit with non-default contrasts is a different design than
  the one its coefficients belong to. Because the two are multiplied together
  positionally, a categorical exposure carrying non-default contrasts, an
  ordered factor exposure, a two-level factor exposure with non-default
  contrasts, or a `contrasts` option changed between fitting and calling `ipw()`
  produced silently wrong estimates, sometimes reversed in sign. A propensity
  score model adjusted for a covariate with non-default contrasts instead
  produced a spurious warning and, when `.data` was supplied, a weight-mismatch
  error against weights that were correct. All of these are now fixed. Results
  from models fit under the default contrasts are unchanged.

* The `glm` methods of `wt_ate()`, `wt_att()`, `wt_atu()`, `wt_atm()`,
  `wt_ato()`, `wt_entropy()`, and `wt_cens()` now supply the numeric methods
  with the probability of the level resolved as focal. Fitted values give the
  probability of the response's second level, and they are subtracted from one
  whenever the focal level resolves to the first level instead. Weights built from a fitted model with a flipped focal level were
  previously inverted, weighting each unit by the other group's probability;
  they are now correct. Calls that leave the focal level at its default are
  unchanged.

* The first argument of `ipw()` is now `wt_mod` rather than `ps_mod`, on the
  generic and on every method. It still takes the fitted weighting model, but a
  call that supplied it by name, such as `ipw(ps_mod = fit, outcome_mod = out)`,
  now errors. `.data`, `estimand`, `ps_link`, and `conf_level` follow `...` in
  the method signatures and must be supplied by name. A value passed into one
  of those slots by position, such as the `ipw(ps_mod, outcome_mod, dat)` form
  that earlier releases accepted, falls into `...` and errors.

* `ipw()` now errors on anything left in `...`, which the methods do not use.
  This catches a misspelled argument name and an argument supplied by position
  where the signature requires a name; both were previously accepted in silence
  and the default used.

* The `ipw` result object renames its `ps_mod` component to `wt_mod` to match,
  so code reading `result$ps_mod` now gets `NULL` rather than the fitted
  weighting model.

* `ipw()` is an S3 generic. The generic now belongs to the causalgenerics
  package, which propensity imports from, re-exporting `ipw()` and registering
  its `glm`, `multinom`, `lm`, and default methods against it. `estimand()`,
  `estimand<-`, and `is_causal_wt()` come from causalgenerics as well, and the
  `causal_wts` methods that `psw` objects inherit are defined there rather than
  here. Calling any of them is unchanged, but a method registered against a
  propensity-owned `ipw` generic is no longer dispatched to; register it
  against `causalgenerics::ipw()` instead. The causalgenerics and deli packages
  are new hard dependencies.

* The minimum R version is now 4.3, raised from 4.2.0.

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
  `se_method = "mestimation"` for a covariate-adjusted outcome model. The
  linearization influence functions are those
  of the Hajek weighted-mean estimator, which match the reported g-computation
  estimates only for an outcome model of the exposure alone; the M-estimation
  path handles covariate-adjusted outcome models correctly.

* Both standard error methods now require an outcome model that can represent a
  baseline at every exposure level. `se_method = "linearization"` requires an
  intercept outright: without one the marginal mean under no exposure was
  pinned at the link's zero point rather than estimated, so the reported
  estimates silently stopped matching the Hajek means the influence functions
  describe. `se_method = "mestimation"` rejects a counterfactual design matrix
  that is identically zero, which is what a numeric no-intercept coding such as
  `y ~ z - 1` produces at the zero-coded level. The marginal mean there was
  fixed by the outcome link rather than estimated, and the reported effect could
  carry the opposite sign with less than half the standard error. A saturated
  factor coding such as `y ~ 0 + zf` is a reparameterization rather than a
  missing baseline, and still works on the M-estimation path.

* Corrected the linearization standard errors for `log(rr)` and `log(or)`,
  which were scaled by the reciprocal of the risk ratio and odds ratio,
  respectively (underestimated when the ratio exceeds one, overestimated
  otherwise). Risk-difference standard errors are unaffected.

* Corrected the linearization standard errors for probit and cloglog propensity
  score models. The influence function for the propensity score coefficients
  omitted the GLM score factor `1 / (p (1 - p) g'(p))`, which is 1 only for the
  canonical logit link, so probit and cloglog standard errors were mis-scaled.
  Logit standard errors are unchanged.

* Corrected the linearization standard errors for stabilized ATE weights. The
  weight derivative omitted the stabilizer, so the propensity score correction
  was divided by the group constant with nothing multiplying it back. Standard
  errors change for any analysis weighted by stabilized `wt_ate()` weights, by
  an amount that tracks how far the stabilizer sits from one: on one
  representative sample the risk-difference standard error fell by about 3
  percent under the default marginal stabilization, by about 3.6 percent under
  a `stabilization_score` taken from a numerator model, and by about 2.7
  percent under a scalar score of 0.5, while a score averaging close to one
  moved it by about 0.1 percent. The log risk ratio and log odds ratio standard
  errors move in the same direction by somewhat less. Unstabilized weights are
  unaffected.

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

* For a continuous exposure, `ipw()` now accepts only a plain `lm()` or a
  gaussian `glm()` propensity score model and errors on any subclass, such as a
  robust or an additive fit, naming the class it was given. The stacked system
  carries an ordinary least-squares score block, so a subclass was previously
  solved to the least-squares root (a 3 percent estimate drift on a robust fit
  in one measured sample) or had its design reconstructed as something other
  than the one it was fit with.

* `ipw()` now rejects an outcome model whose formula omits the exposure, for
  binary and categorical exposures and on both standard error methods. The
  counterfactual designs are built by setting the exposure to each level in
  turn, so such a model gave one identical design per level: with `.data` the
  result was a table of near-zero estimates with `NaN` standard errors, and
  without `.data` the call reported a misleading "supply `.data`" hint that led
  straight into it.

* The binary M-estimation path now validates the propensity score link at entry
  against logit, probit, and cloglog, whether the link comes from `ps_link` or
  from the model's family. A log or identity link was previously accepted
  silently, and a cauchit link failed late with an internal message.

* Binary `wt_att()` and `wt_atu()` now record the focal level they were built
  with whenever the caller names `.focal_level` or `.reference_level`, and
  `ipw()` errors when that level is not the one it treats as focal
  (the second sorted level of the exposure), naming both levels and directing
  you to relevel and refit. Previously the linearization path silently mirrored
  the estimand's derivative roles, and the M-estimation path reported an
  unrelated weights mismatch.

* For a categorical exposure, `ipw()` now resolves the exposure levels against
  the fitted `nnet::multinom()` model's own training levels rather than the
  ordering `.data` implies. An exposure column supplied as character, or as a
  factor with a different level order, now gives the same answer as the column
  the model was fit to, and a value the model never saw errors naming it.

* `ipw()` now converts a factor or logical outcome response to the model's 0/1
  coding on both standard error methods, following glm's convention (the first
  factor level is failure, every other level is success). Previously a factor
  response crashed the M-estimation solve or produced `NA` linearization
  standard errors with factor-arithmetic warnings.

* Factor exposures now work on the linearization path, recoded to 0/1 with the
  second factor level as the exposed group, matching the M-estimation path;
  logical exposures are recoded the same way. Previously a factor exposure
  crashed the linearization variance with factor-arithmetic errors.

* The binary M-estimation path now reconstructs the propensity score design from
  `.data` when the model frame is unavailable (for example, a model fit with
  `model = FALSE` whose fitting data is gone) and reports an informative error
  directing you to supply `.data`, matching the categorical and continuous
  paths. Previously it failed with a raw "object not found" error that `.data`
  did not resolve.

* The linearization path now reconstructs the propensity score design through
  the same extractor as the M-estimation path, so a propensity score model
  whose fitting data is gone rebuilds its design from `.data` and reports the
  same informative error when `.data` is missing. It also rejects a `.data`
  whose row count disagrees with the fitted models; such a `.data` was
  previously accepted silently and shrank the reported standard errors by
  roughly the square root of the row ratio.

* Requesting the `atu` or `entropy` estimand with `se_method = "linearization"`
  now errors with the documented message directing you to
  `se_method = "mestimation"`, matching the categorical and continuous paths.
  Previously it raised a bare internal argument-matching error.

* The error raised when `ipw()` cannot determine the estimand from the weights is
  now attributed to `ipw()` rather than the internal `check_estimand()` helper.

* A more-than-two-level exposure on the linearization path now raises the
  informative binary-only error instead of an internal
  "message must be a character vector" crash.

* A `.focal_level` longer than one now errors informatively, naming the
  argument, instead of raising a raw length-coercion error from an internal
  comparison.

* Supplying a non-`NULL` `ps_link` with a multinomial or continuous propensity
  score model now errors, instead of silently ignoring it; `ps_link` applies
  only to a binomial glm propensity model on the binary path. The continuous
  rejection covers both entries, the `lm()` method and the gaussian branch of
  the `glm()` method.

* A propensity score model with a matrix response (such as
  `cbind(successes, failures)`) now errors with one consistent message on both
  standard error methods, instead of diverging into a misleading case-weights or
  adjusted-outcome error; a binary exposure must be a single-column response.

* `print()` of an `ipw` object now formats the estimate and standard error as the
  coefficient pair and the z statistic as the test statistic; previously the
  columns were misassigned, so the z statistic printed at full precision and the
  lower confidence bound was truncated.

* The M-estimation solve no longer produces `NaN` at saturated or extreme
  propensity scores: the binary entropy tilt and the categorical propensity
  score reconstruction are both guarded at the endpoints of the unit interval
  and against overflow in the linear predictor.

* `wt_ate()` and `wt_cens()` now record a user-supplied `stabilization_score`
  attribute on the returned weights, readable with the new `stabilization_score()`
  accessor.

* A per-observation `stabilization_score` is now dropped, with a new
  `propensity_stabilization_score_warning`, when an operation changes a `psw`
  vector to a different, non-zero length, since the score cannot be re-indexed
  to the observations that remain. A scalar score, length-preserving operations
  such as arithmetic, and zero-length results such as `w[0]`, which leave no
  observations for a score to line up against, are unaffected.

* Fixed `broom::tidy(glm_fit, conf.int = TRUE)` failing on GLMs weighted by
  `psw` vectors. `confint.glm()` builds profile-likelihood intervals via
  `profile.glm()`, which refits through `glm.fit()`; the refit indexes
  `weights[good]` with a matrix subscript, which `[.vctrs_vctr` rejected. A
  matrix or array subscript on a `psw` vector now falls back to base R linear
  indexing; every other subscript still goes through `[.vctrs_vctr`.

* Comparison operators on `psw` (`==`, `!=`, `<`, `>`, `<=`, `>=`) now
  short-circuit `vec_equal()` / `vec_compare()` and return a logical vector
  silently. Previously each comparison fired a `propensity_class_downgrade`
  warning via `vec_ptype2.psw.double()`, producing 100+ warnings during a
  single `tidy(glm, conf.int = TRUE)` call. Combine and cast paths still warn.

* Added a `NEWS.md` file to track changes to the package.
