# propensity 0.1.0.9000 (development version)

* `wt_ate()` and `wt_cens()` gain a `.density` argument, which chooses the
  family of the conditional density a continuous exposure's weights are a ratio
  of. That density was a normal one and nothing else before, which is a strong
  claim to make about the residuals of a model for a dose. `.density` accepts
  the strings `"normal"` (the default), `"laplace"`, and `"kernel"`; a
  specification built by the new `dens_normal()`, `dens_laplace()`,
  `dens_t()`, `dens_kernel()`, or `dens_fn()`; or a bare function of one
  argument. A tail heavier than the normal's, from `dens_t()` or `"laplace"`,
  holds down the weight of a unit whose exposure the model fits poorly, and a
  kernel estimate assumes no shape at all, at the cost of a density that is not
  a smooth function of the model's parameters. `.density` sits after `...`, so
  it is supplied by name, and it describes a continuous exposure alone: a
  family other than the default is refused for a binary or categorical
  exposure, whose weights are not a ratio of densities, the way `.sigma`
  already was.

  Both densities in the ratio are now evaluated on a standardized residual,
  `(A - mu) / sigma` for the conditional density and the exposure standardized
  by its own mean and standard deviation for the marginal one, and each is
  divided by the spread that standardized it. That factor is the Jacobian of
  the change of variable, and it returns both densities to the exposure's own
  units, so every family is read on one scale and the normal family returns the
  weights the package returned before, to within a rounding error in the last
  binary digit. Whatever the density gives back is checked before it becomes a
  weight: one finite, non-negative value for each standardized residual, and
  not zero at every one of them. Anything else is an error of class
  `propensity_density_error` whose message reports the standardized residuals
  it failed at. The weights record the family they were built from, alongside
  the numerator that stabilized them and where the residual spread came from;
  read it back with `density_meta()`.

* `wt_ate()` and `wt_cens()` gain a `numerator` argument, which chooses how the
  marginal density that stabilizes a continuous exposure's weights is arrived
  at. `"marginal"`, the default and the behavior of every earlier version,
  reads the family `.density` names at the population mean and standard
  deviation of the exposure; those two moments are parameters of the weights,
  and `ipw()` estimates them alongside the rest of its parameter vector.
  `"integrated"` marginalizes the conditional density numerically instead,
  averaging it over the units at each of 50 points spanning the exposure and
  interpolating that average back to each observed exposure, which is what
  WeightIt has done since version 2.0.0 and is what makes the two packages'
  continuous weights agree. It estimates no parameters of its own.

  The default did not change, and it was studied before it was left alone: a
  simulation study run for this package compared the two over linear, skewed,
  heteroskedastic, heavy-tailed, and bimodal designs. The integrated numerator
  never improved bias, root mean squared error, or weighted covariate balance
  by more than Monte Carlo noise, including in the scenarios drawn to favor
  marginalizing, and in the heavy-tailed cells it was worse: it inflated the
  root mean squared error, and in some replicates the interpolated density came
  back negative. `"integrated"` is offered for agreement with WeightIt, and for
  a conditional density whose marginal has no closed form in the same family.
  Because it is read off an interpolation rather than a formula, it can dip
  below zero where the density on the grid comes close to it; the interpolated
  values are held to the same output check any density is held to, so that
  becomes an error of class `propensity_density_error` rather than a negative
  weight. An integrated numerator needs a conditional density to marginalize,
  so it is refused with an error of class `propensity_numerator_error` for a
  binary or categorical exposure, with `stabilize = FALSE`, with a
  `stabilization_score`, and with any `.sigma`. A unit whose exposure or fitted
  conditional mean is missing is not read by the average and has no weight, as
  it has none under any other numerator, and a model with nothing to condition
  on gives weights of exactly one rather than the grid's approximation to them.
  The weights record the numerator they were built from; read it back with
  `density_meta()`.

* `wt_ate()` and `wt_cens()` now read a continuous exposure's conditional means
  from the model that fit them, rather than from a `glm()` alone. Both gain an
  `lm` method, which a `MASS::rlm()` reaches by inheritance, and their `glm`
  method reads a `gaussian()` fit under any of its links as well as an
  `mgcv::gam()` fit with it. Each of these classes reports its conditional mean
  on the scale of the exposure, so a log, inverse, or square root link never has
  to be undone, and every one of them is spread by the same pooled residual root
  mean square. `rlm` reports a robust scale estimate of its own in `fit$s`,
  which resists the extreme residuals rather than pooling all of them, and is
  used only when it is asked for with `.sigma = fit$s`. A family whose
  spread changes with its fitted values, such as `poisson()` or
  `quasi(variance = "mu")`, describes a different density for every unit, which
  a single spread cannot stand in for, and is refused with an error of class
  `propensity_model_family_error`. `quasi(variance = "constant")` is the
  gaussian variance under another name and is accepted.

  The binary path is now held to the same rule from the other side, which
  changes what some calls return. Only `binomial()` and `quasibinomial()` fit
  the probability those weights divide by, so a gaussian model of a zero-one
  response, an `lm()` included, is refused with the same error class rather than
  having its fitted values read as propensity scores. A linear probability model
  is not held to the unit interval, so the refusal reads the model rather than
  the values it fitted, and comes before the check on the range of a propensity
  score. The estimands whose weights are not a ratio of densities take no model
  of a continuous exposure: `wt_att()`, `wt_atu()`, `wt_atm()`, `wt_ato()`, and
  `wt_entropy()` have no `lm` method.

* A continuous-exposure `ipw()` fit whose outcome model reads the exposure
  through more than one design column, such as `y ~ A + I(A^2)` or a basis like
  `poly(A, 2)` or `splines::ns(A, 3)`, now records the conditional reading and
  says so once. No coefficient of a curve is the effect of the exposure, since
  the dose response has a different slope at every dose, so there is no marginal
  reading of such a fit to present. `ipw()` refuses `effects = "marginal"` with
  an error of class `propensity_ipw_effects_error`, and the result declares the
  conditional reading as the only one it supports, so `as_marginal()`, `coef()`,
  `vcov()`, `confint()`, `tidy()`, and `as.data.frame()` refuse it with an error
  of class `causalgenerics_unsupported_reading_marginal`. Marginalizing the
  curve over the observed doses is a separate estimand that this package does
  not compute; use the marginaleffects package on the conditional result:
  `avg_slopes()` for slopes, `avg_comparisons()` for contrasts, and
  `avg_predictions()` for causal dose-response functions
  (<https://marginaleffects.com/chapters/interactions.html>). A fit whose
  exposure enters through one column is unchanged. A caller who names no
  reading includes a wrapper that forwards the whole `effects` default, which
  is a set of readings rather than one of them, so such a call is announced and
  given the conditional reading rather than refused for the reading the default
  resolves to first.

* A binary or categorical `ipw()` result now reports the counterfactual mean at
  each exposure level, one row per level under the effect label `"mean"`, ahead
  of the contrasts built from those means. A risk difference is now read beside
  the two risks it is a difference of, on either standard error route and within
  each stratum of a `.by` request. The rows are covered by the stored covariance
  and carried by `coef()`, `vcov()`, `confint()`, `tidy()`, `as.data.frame()`,
  and the pooled results, and `exponentiate = TRUE` leaves them alone, a
  counterfactual risk being no kind of ratio.

* The p-value of a reported row now comes from the upper tail of the normal
  distribution directly rather than from one minus the lower tail, which carries
  no precision past about 1e-16 and returned an exact zero for any test
  statistic beyond about 8. The counterfactual mean rows often make such a
  statistic, a fitted risk being large beside its standard error. The printed
  report is unchanged, since a p-value that small prints as `< 2.2e-16` either
  way, but the value the result stores now carries its magnitude.

* The `estimates` frame of a binary `ipw()` result gains the `contrast` column a
  categorical result already carried, naming the level each mean row belongs to
  and the pair of levels each contrast row compares, such as `"1 vs 0"`. Row
  labels on every surface of the result name that contrast, so code matching a
  label such as `"rd"` exactly needs `"rd 1 vs 0"`.

* The exposure-type machinery that resolves `exposure_type` and detects the type
  of an exposure now lives in causalgenerics, where the other packages in the
  ecosystem can reach it, and propensity calls it there rather than keeping its
  own copy. The refusal of a type a weight function does not answer says the
  same thing it did, but it now carries the classes
  `causalgenerics_unsupported_exposure_type` and `causalgenerics_error` in place
  of `propensity_wt_not_supported_error`. Code that catches that refusal by
  class needs the new one.

* `ipw()` gains a `.by` argument, which reports the effect of the exposure
  within each level of a modifier alongside the effect for the sample as a
  whole. The reported rows come in three blocks: the overall rows unchanged from
  a fit without `.by`, one block per level of the modifier, and one block per
  level against the reference level, so the difference between two strata is
  reported rather than left to be taken by hand. A stratum's effect is
  g-computation restricted to that stratum, which means it is the effect in that
  stratum's own covariate distribution, and every row is a parameter of the same
  stacked system, so the standard error of a difference between strata accounts
  for the two having been estimated together. Fitting each stratum on its own
  subset reports the same stratum effects with no covariance between them, which
  leaves that difference untestable. The stratum rows and their contrasts report
  the risk difference and the log risk ratio but no odds ratio, since an odds
  ratio is noncollapsible and a difference of two of them moves with the outcome
  distribution in each stratum whether or not the effect does. `.by` requires
  `se_method = "mestimation"` and a binary or categorical exposure.

* `ipw()` now reports a joint intervention on two discrete treatments. Cross
  them with `causalgenerics::joint_exposure()`, which records the crossing on the
  vector, and the result is written in the two treatments rather than as each
  cell against a reference cell: one counterfactual mean per cell, each
  treatment's simple effects within the levels of the other, and the interaction
  between them, on the risk difference scale as the additive interaction and on
  the log risk ratio scale as the multiplicative one. The simple effects include
  the comparisons that are not against the reference cell, which the
  vs-reference reporting cannot express at all. The interaction is reported once,
  under the first treatment's framing, because it is symmetric in the two. Every
  row is a parameter of one stacked system, so an interaction equals the
  corresponding double difference of the cell means exactly.

* `wt_joint()` and `joint_wt_models()` build the weight for a joint intervention
  from one treatment model per treatment, which is the sequential factorization
  the weight really has: the first treatment given the covariates, times the
  second given the first treatment and the covariates. Pass the container to
  `ipw()` as `wt_mod` with an outcome model reading both treatments and the
  reported surface is the declared crossing's, row for row. Prefer this route
  when the two treatments call for different adjustment sets, or when the
  dependence of the second on the first is what you want to model. The
  factorization is checked rather than assumed: `joint_wt_models()` refuses a
  second model that does not condition on the first treatment, since the product
  of two marginal weights is a different quantity that nothing downstream can
  tell apart, and `wt_joint()` requires a continuous component to be stabilized.
  The second treatment may be a dose, recorded with an `lm()` or an
  identity-link gaussian `glm()`, in which case there are no cells to report and
  the surface is the marginal structural model's own coefficients. A model
  written in bare treatment terms reports rows naming the treatment each one
  varies and where it is evaluated; every other treatment-reading model, a
  transformed term or a basis such as `poly()` or `splines::ns()` among them,
  reports one row per coefficient named after the coefficient it reports. Both
  treatment models are stacked alongside the outcome model, so the standard
  errors account for having estimated both.

* A new vignette, `vignette("effect-modification")`, separates effect
  modification from interaction as causal questions and works both on one
  dataset where they give different answers. It covers when a stabilizing
  numerator may condition on the modifier, why no odds ratio is reported for
  effect modification, and the reasoning behind the checks on joint weights.

* `tidy()` on a pooled result gains an `effects` argument, which reports the
  reading it names rather than the one the pooled result records: `"conditional"`
  returns the pooled coefficients of the outcome models, one row per coefficient,
  and `"marginal"` returns the pooled causal contrasts. `NULL`, the default,
  follows the reading the pooled result records. The pooling has already happened
  by the time the method sees the result, so naming a reading reports a table the
  pooling built rather than pooling a second time, and the pooled result is left
  as it is. `"fixed"` is not among the values it takes, since `mice::pool()` asks
  for that reading from the tidier of each imputation's own fit rather than from
  a pooled one.

* `pool_ipw()` now pools both readings of the results it is given, storing the
  one they record and keeping the other beside it, so `as_marginal()` and
  `as_conditional()`, both re-exported here, move a pooled result between the two
  readings after the pooling rather than only before it. A reading the pooling
  could not compute is refused rather than reported under the other one's name:
  the conditional reading needs the covariance the joint estimation of the
  weights and the outcome implies, and a set of `se_method = "linearization"`
  fits records none, so such a pool holds the marginal reading alone and asking
  it for the other errors under the causalgenerics classes
  `causalgenerics_pool_missing_surface_conditional` and
  `causalgenerics_pool_missing_surface`, carrying the explanation the pooling
  recorded. The Multiple imputation section of `ipw()` documents both routes to a
  reading, the one taken before the pooling and the one taken after it.

* `pool_ipw()` is now re-exported, so the whole multiple-imputation workflow can
  be written without reaching for a second namespace for its last step: fit the
  analysis once per imputed dataset inside `mice::with()`, then pool the results.
  `tidy()` and `glance()` methods report what it returns, the first as the pooled
  estimates in the columns broom conventions use, the second as one row
  describing the pooled fit. `ipw()` gains a Multiple imputation section
  documenting the workflow, why the propensity model belongs inside the
  per-imputation expression, and how the pooling compares with `mice::pool()`.

* `tidy()` now tolerates two of the arguments `mice::pool()` passes every tidier
  it calls, so an `ipw()` result can be fit once per imputed dataset and pooled by
  Rubin's rules. `parametric` is accepted and ignored, and `effects = "fixed"`
  asks for the reading the result records, which is what `NULL` asks for. The
  dots stay closed, so every other stray argument is still refused, and the
  accessors underneath still accept only `"marginal"` and `"conditional"`. A
  categorical result pools grouped by `term` and `contrast` together, so each
  effect measure of each contrast is combined with itself rather than with its
  neighbor.

* The column naming the contrast on a categorical `ipw()` result is now called
  `contrast` in the estimates table the result stores, where it was called
  `comparison`. `contrast` is the name every surface of a result already
  reported, so `tidy()`, `as.data.frame()`, and the printed table are unchanged
  and only the stored frame moves; code reading `result$estimates$comparison`
  needs updating. Results built by earlier versions stay readable, because
  causalgenerics accepts either name when it reads a stored frame.

* The marginal reading of `tidy()` is now the result's own `as.data.frame()`
  table read as a tibble rather than a second assembly of the same columns. The
  rows, the columns, their order, and their values are the same as before, and
  the covariance of the effects that frame attaches as an attribute is the one
  thing the tibble does not carry, a tidied table being its columns. The
  interval rebuilder and the confidence level validator this method kept for
  itself go with the second assembly.

* `tidy()` now refuses a bad `conf.int`, `conf.level`, or `exponentiate` under
  the causalgenerics condition naming the argument at fault, such as
  `causalgenerics_invalid_argument_conf.level`, rather than under a class of this
  package's own. These are the three arguments the method shares with
  `as.data.frame()`, which is where they are now validated, in both readings and
  before either is assembled, so one value cannot be well formed in one reading
  of a result and refused in the other. Code catching
  `propensity_conf_level_error` by name needs updating. The levels `conf.level`
  accepts are unchanged, while `conf.int` and `exponentiate` now take a single
  `TRUE` or `FALSE` and nothing else, where a bare `if` had accepted anything it
  could read as one.

* `ipw()` gains an `effects` argument on every method, which records the reading
  the result it builds presents: `"marginal"`, the default, for the
  population-averaged causal contrasts every method reported before, or
  `"conditional"` for the outcome model's coefficient surface. Both surfaces are
  computed whichever value is named, so the argument settles which one the result
  presents and nothing else. `as_marginal()` and `as_conditional()`, which move a
  result between the two readings afterwards, are re-exported here, so a result
  can be flipped without loading causalgenerics. A printed result names its
  reading beside the estimand and again over the table it decides.

* In the conditional reading, `coef()` reports the outcome model's coefficients,
  and `vcov()` and `confint()` report them against the block of the joint
  sandwich that every route stacking estimating equations attaches to the outcome
  model: `se_method = "mestimation"` for a binary exposure, and the categorical
  and continuous routes. These are the coefficients of the weighted outcome model
  with the uncertainty of estimating the weights carried into them, rather than
  the ones the model's own fitting routine reports. A linearization fit stacks no
  such system and so has no such block: its conditional reading errors from
  `vcov()` and `confint()`, and prints the coefficients under a note saying the
  standard errors are not reported.

* `tidy()` gains an `effects` argument, which reads a result in the reading it
  names rather than in the one the result records: `"conditional"` returns the
  outcome model's coefficient surface, one row per coefficient, in the columns
  and the order the marginal table of causal contrasts uses, so rows of the two
  stack. Its standard errors are the block of the joint sandwich the accessors
  report, which keeps the tidied table and the printed one the same numbers, and
  a linearization fit stacks no such system and so errors rather than reporting
  the covariance the outcome model computed for itself. `NULL`, the default,
  follows the reading the result records. `glance()` and `augment()` describe the
  fit and its observations rather than its estimates and report the same thing in
  either reading.

* In the conditional reading, `exponentiate` follows the link of the outcome
  model rather than the labels of the rows, a coefficient table having none to
  pick out: a `logit` or a `log` link puts every coefficient on a scale an
  exponential undoes, and the estimate and, when they were asked for, the
  interval bounds move to the natural scale while the standard error, the
  statistic, and the p-value stay describing the link scale. Every other link
  errors rather than exponentiating coefficients that are not on such a scale.

* An `ipw()` result now answers the standard model accessors: `coef()` and
  `confint()` for the reported effects, `vcov()` for their covariance, `nobs()`
  and `df.residual()` for the counts describing the fit, and `weights()` for the
  `psw` vector the outcome model was fit with. The methods belong to
  causalgenerics, which owns the shared result class; what propensity supplies is
  the covariance of the effects it reports, recorded on the estimates table.
  Coefficients are named for the effect measure, and for the effect measure and
  the contrast together where a categorical exposure reports one row per
  contrast. `glance()` reports the counts through the same accessors, so the
  two surfaces cannot disagree. This raises the causalgenerics requirement to
  0.1.0.9000.

* Under `se_method = "mestimation"` the `wt_mod` and `outcome_mod` an `ipw()`
  result holds now carry their block of the joint sandwich, so `vcov()` on either
  reports a covariance that accounts for the whole system having been estimated
  from the same data rather than the one the model's own fitting routine reports.
  The models passed to `ipw()` are left as they were fit, and everything else
  about the stored copies is unchanged, so predictions and coefficients read the
  same from either. Linearization stacks no such system, so its component models
  are stored exactly as they arrived.

* New `tidy()`, `glance()`, and `augment()` methods describe an `ipw()` result in
  the three shapes broom defines: one row per estimate, one row for the fit, and
  one row per observation. Each is registered against the generic the generics
  package owns, and all three generics are re-exported here, so the verbs are
  available from propensity itself without loading broom. The generics and tibble
  packages move from Suggests to Imports and are now hard dependencies.

* `tidy()` returns the estimates as a tibble under the column names broom
  conventions use, `term`, `estimate`, `std.error`, `statistic`, and `p.value`,
  with `contrast` naming the contrast on a categorical exposure. `conf.int`
  adds the `conf.low` and `conf.high` bounds, `conf.level` reports them at a
  level other than the one the result was fit at by rebuilding them from the
  estimate and its standard error, and `exponentiate` puts the `log(rr)` and
  `log(or)` rows on their natural scale exactly as the result's own
  `as.data.frame()` method does. Nothing is dropped and nothing is re-estimated.

* `glance()` returns one row describing the fit rather than its estimates: the
  estimand the weights target, the number of observations the standard errors
  were estimated from, and the residual degrees of freedom of the stacked
  M-estimation system, which are those observations less the parameters the
  system solves for. A result reporting several effect measures, or several
  contrasts of a categorical exposure, still returns exactly one row.
  `se_method = "linearization"` stacks nothing and records no parameter count, so
  the observations are the outcome model's and the degrees of freedom are `NA`.

* `augment()` returns the data the outcome model was fit on with the propensity
  score, the weights, the fitted values, and the residuals attached as
  dot-prefixed columns, one row per observation and no observation dropped. A
  categorical exposure carries a probability for every level and so gets one
  `.propensity_<level>` column per level in place of `.propensity`. Unlike
  broom's `augment()` methods, the model frame's `(weights)` column is not
  carried through: the weights appear once, as `.weights`, the `psw` vector the
  outcome model was fit with, so the class and the estimand they record travel
  with them. Pass the modeling data to `data` to carry those columns on a frame
  that also holds the covariates the outcome formula left out.

* `augment()` refuses a result whose propensity score model and outcome model
  were fit to different rows, with an error of class
  `propensity_augment_alignment_error`, rather than reporting one observation's
  propensity score beside another observation's fitted value. Two models of the
  same data most often part this way over missing values, each dropping the rows
  a variable of its own is missing on. What is compared is the number of
  observations each model produced an answer for and, when the outcome model's
  frame names the exposure, the exposure values of the two model frames position
  by position, reading a factor and the numbers its labels spell as one encoding
  of one exposure. Exposure values that disagree prove the two models hold
  different observations; exposure values that agree prove nothing, two different
  sets of rows being free to carry the same sequence of values, and two encodings
  that cannot be read onto each other, such as a factor of `"a"` and `"b"`
  against a recoding of it as 0 and 1, prove nothing either way. The check is
  one-sided in that direction deliberately: it refuses only what it can prove, so
  a fit whose models do describe the same observations is never refused over the
  labels the rows of either frame carry or the encoding its exposure is written
  in.

* `augment()` refuses a frame that already holds a column it would add, with an
  error of class `propensity_augment_column_error` naming every column that
  clashes. Such a frame previously returned one naming the same column twice, and
  a `.resid` column of the caller's was written over whenever the frame also held
  the outcome to difference. The names in the way are the ones the call would
  actually add, so `.resid` is among them only for a frame it would be added to,
  and the propensity columns of a categorical fit are the `.propensity_<level>`
  columns rather than `.propensity`.

* `augment()` refuses a result whose outcome model was fit without weights, under
  `propensity_ipw_weights_missing_error`, the class `ipw()` reports the same fact
  under. An `ipw()` outcome model is weighted by construction, so `.weights` has
  nothing to report for a result built around one that is not.

* `wt_entropy()` is documented as computing entropy weights, which tilt the
  propensity score by its own entropy in the sense of Zhou, Matsouaka, and
  Thomas (2020). The package overview called them entropy balancing weights, and
  the weight function help, the README, and the getting started vignette named an
  entropy-balanced population. Entropy balancing is a different method, solving
  for weights that satisfy exact covariate moment constraints rather than tilting
  a fitted propensity score. Behavior is unchanged.

* The weight functions document the open interval they require of a categorical
  propensity score matrix, and how it differs from the narrower rule `ipw()`
  applies to the scores it rebuilds from its own propensity score model. A
  separated `nnet::multinom()` fit reaches the endpoints readily, and neither
  `ps_trim()` nor `ps_trunc()` accepts such a matrix either, so the remedy is to
  bound the fitted probabilities away from 0 and 1 or to refit. `ipw()`,
  `ps_trim()`, and `ps_trunc()` carry a pointer to that text. Behavior is
  unchanged.

* `ps_trim()` and `ps_trunc()` document what `method = "cr"` does when the
  exposure groups' propensity score distributions do not overlap at all.
  `ps_trim()` trims every observed unit, a truthful record of an empty overlap
  region, while `ps_trunc()` errors with class `propensity_no_overlap_error`,
  since there is no range left to bound the scores to. Behavior is unchanged.

* A deprecated argument now reports the call that supplied it. lifecycle decides
  whether a deprecation belongs to the caller or to the package that raised it,
  and the warnings for `.treated`, `.untreated`, and `ps` were raised from
  helpers that named no calling environment, so every one of them arrived with a
  bullet asking the reader to report an issue against propensity for an argument
  they had written themselves.

* `wt_ate()` and `wt_cens()` now refuse a `.sigma` that holds a value at or
  below zero, a value without a bound, or nothing but missing values, with an
  error of class `propensity_sigma_error`. A standard deviation reached
  `dnorm()` as it arrived: a negative one came back as `NaN` weights under a
  base R warning about them, a zero came back as infinite weights, and a missing
  one came back missing.

* `ps_trim()`, `ps_trunc()`, `ps_calibrate()`, and `ps_tilt()` now name
  `.propensity` when they are called without propensity scores, with an error of
  class `propensity_missing_arg_error`. Dispatch reported base R's missing
  argument message, which names the formal rather than the two spellings the
  scores may arrive under.

* `ps_trim()` and `ps_trunc()` now name the threshold they were given and the
  `1/k` the width of a propensity score matrix imposes on it. The refusal read
  `Invalid trimming threshold (delta >= 1/k)`, which left the caller to work out
  both numbers, one of which was their own.

* `ps_trim()` and `ps_trunc()` now refuse `method = "cr"` when one exposure
  group has no observed propensity score, with an error of class
  `propensity_no_data_error`. The bounds are the lowest score among the focal
  units and the highest among the reference units, and over none of them `min()`
  and `max()` returned `Inf` and `-Inf` under a base R warning: `ps_trim()` then
  trimmed every unit, and `ps_trunc()` reported distributions that do not
  overlap, which describes a different problem.

* `ps_trunc()` now validates the percentiles it bounds a matrix of categorical
  propensity scores at, as it already did for a vector. Bounds that crossed
  pinned every score to the bound on the far side of the other, and a
  probability outside the unit interval was refused in base R's words about
  `quantile()`, an argument the caller never wrote.

* `ps_trim()` and `ps_trunc()` no longer record a bound the method never read.
  `method = "adaptive"` and `method = "cr"` announce that `lower` and `upper` are
  ignored and then work their cutoffs out from the scores, but the ignored bounds
  were kept in the metadata, so two trimmings that ran the same rule to the same
  cutoff described themselves differently and combining them warned and fell back
  to numeric.

* The percentile cutoffs `ps_trim()` and `ps_trunc()` record are no longer named
  for the probability they were read at. `quantile()` names its result `"5%"` or
  `"95%"`, which says nothing about the cutoff and reappeared wherever the cutoff
  was printed or compared.

* `ps_trim()` and `ps_trunc()` now refuse `.focal_level`, `.reference_level`,
  `.treated`, and `.untreated` on the categorical path, with an error of class
  `propensity_unsupported_arg_error` naming the argument. A categorical exposure
  is described by one column of scores per level, so no level is treated as
  focal; the arguments were declared and never read, and the data frame method
  dropped them on the way to the matrix method.

* A condition raised by `ps_trim()` or `ps_trunc()` now names the function the
  caller wrote, whatever shape the propensity scores arrived in. A data frame of
  scores reaches the matrix method on the categorical path, and the vector
  method on the binary one, by a call rather than by dispatch, so every error
  and warning either route raised was reported against `ps_trim.matrix()`,
  `ps_trunc.matrix()`, `ps_trim.default()`, or `ps_trunc.default()`, methods the
  caller has no reason to know of, while the same condition on a matrix or on a
  bare vector named the generic.

* Refusing to cast one `ps_calib` to another now names the disagreement, as the
  `ps_trim` and `ps_trunc` casts already do. A `ps_calib` is printed with its
  method but not with whether the fit was smoothed, so two types that differ
  only in that rendered identically and the refusal read as a type that cannot
  be converted to itself.

* `ps_calibrate(smooth = TRUE)` no longer requires mgcv when it is going to fall
  back to a straight line. Whether enough distinct propensity scores are present
  to place the knots of a spline is now settled first, so a calibration that
  fits a logistic regression is no longer refused over a package it takes no
  part in.

* Printing a `ps_trim` column inside a tibble is documented as slicing the
  column for display, which raises the record-drop warning. The warning is
  truthful and describes the vector built for the display; the column is
  unchanged.

* `ps_trim()`, `ps_trunc()`, `ps_calibrate()`, and `ps_tilt()` now take the
  propensity scores in `.propensity`, the name the weight functions already read
  them under. These four called the argument `ps`, so a call written against one
  family was refused by the other. The old name still works and is deprecated:
  supplying `ps` warns and reaches the same result, supplying both `ps` and
  `.propensity` is an error naming the two spellings, and `ps` will be removed in
  a future release. One consequence is inherent to the rename: a call that names
  `ps` and leaves a later argument positional, such as
  `ps_trim(ps = x, "adaptive")`, binds the positional argument to `.propensity`
  rather than to `method`, and is then refused for supplying the scores twice. A
  call that names `ps` must name the arguments after it as well. The error the
  weight functions raise for a propensity score outside (0, 1) now names
  `.propensity` too, which is what those functions have always called the
  argument it reports on.

* `ps_refit()` now refits a model whose formula transforms a term, such as
  `z ~ log(x)` or a spline basis, without being handed the data. The default
  `.data` came from the stored model frame, which holds each term already
  computed, so refitting looked for a variable named `x` among columns named
  `log(x)` and failed. The data the model names are now read back and cut down
  by row name to the rows the model analyzed, so a transformation is recomputed
  from the retained rows alone, which is what refitting on those rows means. A
  spline's knots therefore move to where the retained rows place them. A model
  fit without a data argument names none, and its variables are read out of the
  formula's environment instead.

* `ps_refit()` now refits a model whose `weights` or `offset` names a column of
  the data it was fit on. The stored model frame records both under fixed names
  rather than the ones the call reads, so refitting from it reported the column
  as a missing object. The recovered data are the frame the original call read,
  which carries that column, so the weights and offsets follow the retained
  rows. A `weights` or `offset` vector held outside the data still cannot follow
  them and raises an error about differing variable lengths.

* `ps_refit()` now names the likely cause when the data recovered from `model`
  hold fewer rows than the propensity scores describe. Scores predicted from a
  fit with `na.action = na.exclude` are padded back to the full length of the
  data, so they outnumber the rows the fit read. The refusal now says so and
  points at an `na.action` that drops those rows.

* `ps_refit()` no longer puts a `subset` from the original model call to work a
  second time. The retained rows are already a narrowing of the sample that
  `subset` chose, so re-applying it selected among rows it was never about and
  returned propensity scores that were quietly wrong. The `subset` is now
  dropped from the refit call, on both the default and the explicit `.data`
  route. A `subset` passed through `...` is an instruction of its own and is
  still honored. How `ps_refit()` recovers its data, and the limits on
  `weights` and `offset` read from outside the formula, are now documented.

* `ps_calibrate(method = "isoreg")` no longer runs a step described as
  preventing extrapolation beyond the observed range. Each of the two isotonic
  fits was raised to its own minimum over the exposure group that reads it,
  which every value that group reads already clears, so the step changed
  nothing. The calibrated scores are unchanged.

* `ps_calibrate()` now says so when `smooth = TRUE` cannot be honored. A spline
  needs enough distinct propensity scores to place its knots, so with fewer
  than 10 among the observations it reads, those with both a score and an
  exposure recorded, the fit falls back to logistic regression without one. The
  returned metadata recorded the fallback, but nothing at the point of the call
  said the spline that was asked for was never fit. The announcement respects
  `options(propensity.quiet)`, and the threshold is now documented.

* The propensity scores `ps_calibrate()` accepts are documented as the closed
  interval, including exactly 0 and exactly 1. Every other route refuses the
  endpoints, since a score there divides into no usable weight. Calibration
  takes them deliberately, because repairing scores a model pushed to the ends
  of the interval is part of what it is for. The logistic calibration curve
  maps them back inside the interval; isotonic calibration can return a score
  at an endpoint when its pooled block is pure, and such scores are rejected by
  the weight functions in turn.

* `ps_calibrate(smooth = FALSE)` now returns a calibrated score for every unit
  with an observed propensity score, including units whose exposure is missing.
  The fit reads only the units with a recorded exposure, as it always did, but
  the curve it produces was read back only over those same units, leaving a
  unit with a missing exposure uncalibrated. The smooth fit already predicted
  over every unit. Isotonic calibration still returns `NA` for a unit with a
  missing exposure, which has no group fit to be read from, and the behavior of
  each method is now documented.

* `ps_calibrate()` now refuses a matrix of propensity scores with an error of
  class `propensity_type_error` naming the shape as the problem. A one-column
  matrix was flattened and its dimensions dropped without a word, and a matrix
  with more than one column failed a length check that compared its cells
  against the observations in `.exposure` and reported the mismatch as a length
  problem, which said nothing about the shape. A one-dimensional array is still
  accepted, as the weight functions accept it.

* The unused `ps_calib_matrix` class is removed, along with its
  `is_ps_calibrated()` and `print()` methods. No calibration route ever
  constructed one.

* Comparing a `ps_calib` with a number, with an integer, or with another
  `ps_calib` is now silent and returns a plain logical vector. Every comparison
  settled its type through the numeric downgrade first, so asking which scores
  clear a threshold reported dropping metadata the answer never depended on.
  Sizes are still recycled through vctrs, so lengths with no common size raise
  an error rather than take base R's answer.

* `ps_calibrate()` now reads the values of trimmed or truncated propensity
  scores once, up front, so calibrating them no longer reports a class
  conversion the caller never asked for. The range check and each of the
  logistic, smooth, and isotonic fits compared or modeled the scores in their
  original class, which announced the downgrade once per comparison. The
  calibration produced is the one the same values give as a plain numeric
  vector.

* Combining a `ps_calib` with an integer vector now gives a numeric result
  holding the calibrated scores, and casting a `ps_calib` to integer now raises
  a lossy-cast error. The combination found no common type for the two and
  refused to run at all, and the cast refused without naming what an integer
  would have cost. This is the answer `ps_trim`, `ps_trunc`, and `psw` already
  give.

* Casting a numeric or integer vector to a `ps_calib` now carries the
  calibration of the target. It described the result as having been calibrated
  by a method named `"unknown"`, which no argument to `ps_calibrate()` accepts,
  without smoothing.

* `is_unit_truncated()` now answers per unit for a `psw` vector built from
  truncated propensity scores, reading the positions out of the truncation
  record as it already did for a `ps_trunc`. It returned the single flag
  `is_ps_truncated()` answers, which was one value for a query asked once per
  observation. Weights marked as truncated whose record no longer describes
  them raise an error of class `propensity_missing_meta_error` rather than
  report the wrong units, matching `is_unit_trimmed()`.

* Casting one `psw` to another now compares the whole description of the
  weights, the estimand, the stabilization status and score, and the trimmed,
  truncated, and calibrated flags, and refuses with an incompatible-cast error
  naming the field they disagree on. It returned the values unexamined. Most of
  what is compared goes unmentioned when a `psw` is printed, so the refusal
  names the disagreement rather than rendering two identical-looking types.
  Subassignment rests on
  that comparison: `w[1] <- value` casts the replacement to `w`'s type and then
  leaves base R to keep `w`'s attributes, so weights for one estimand could be
  written into weights for another under the target's description, and
  `rbind()` on data frames of weights did the same. Combining is unaffected,
  since it settles its type through `vec_ptype2()`, which still warns and
  returns numeric.

* Combining a `psw` with an integer vector now gives a numeric result holding
  the weights, rather than an integer vector of every weight rounded away. The
  combination reported the class it dropped and then silently changed the
  numbers it kept. This is the answer `ps_trim` and `ps_trunc` already give.
  Casting a `psw` to integer still raises a lossy-cast error when the weights
  have fractional parts, since there the caller named the type.

* Casting one `ps_trim` to another now compares the trimming parameters and the
  refit flag, and casting one `ps_trunc` to another compares the truncation
  parameters, which is the comparison the combine already made. Each refuses
  with an incompatible-cast error naming what disagrees, instead of returning
  the values described by the target. Neither class is printed with the
  parameters being compared, so the refusal names them rather than rendering
  two identical-looking types.

* The line a `ps_trim` is printed under now names the number of observations
  the trimmed count is out of, as in `ps_trim; trimmed 9 of 20`. It ended at
  `of ` with nothing after it, and an object whose record had been dropped was
  described with a trailing separator that introduced nothing.

* The line a `ps_trunc` is printed under now reports its bounds to three
  significant digits. A bound read off the scores is a score, and `"pctl"` and
  `"cr"` put every digit of it in the description.

* Printing a `ps_trim` or `ps_trunc` matrix of categorical propensity scores no
  longer dumps the modification record after the scores. The record is a set of
  index vectors as long as the data, and the header already summarizes it.
  Read it with `ps_trim_meta()` or `ps_trunc_meta()`.

* `ps_trim()` and `ps_trunc()` now announce which column they take from a data
  frame of propensity scores given for a binary exposure. The rule is unchanged,
  the second column of a two column data frame and the first column otherwise,
  and is now documented. `options(propensity.quiet = TRUE)` silences the
  announcement.

* Casting a numeric vector into a `ps_trim` or a `ps_trunc` is documented as a
  type operation rather than a trimming or a truncation: the result is described
  by the target and can hold scores outside the target's cutoffs or bounds,
  including 0 and 1, without being trimmed or pinned.

* Two `ps_trunc` objects are compatible for combining only if they agree on the
  percentiles their bounds were requested at, as well as on the bounds
  themselves. Two objects truncated at the same scores from different
  percentiles are described differently, and now warn and combine as numeric,
  as other metadata mismatches do.

* Combining `ps_trim` objects trimmed the same way, or `ps_trunc` objects
  truncated the same way, now keeps the description of that trimming or
  truncation on the result. The prototype the combination is built on carries
  the metadata across instead of trimming or truncating an empty vector again,
  which had no scores to read a cutoff off and no exposure to be handed:
  `method = "pctl"` came back with missing quantiles, `ps_trim(method = "pref")`
  and `ps_trunc(method = "cr")` raised an error about a missing `.exposure`
  argument the caller had not omitted, `ps_trunc(method = "pctl")` failed the
  cast against the prototype it had just built, and the refit flag set by
  `ps_refit()` was dropped. The result still carries no record of which units
  were modified, since concatenation appends one set of observations to another.

* Casting a numeric or integer vector to `ps_trim` or `ps_trunc` now returns it
  under the trimming or truncation of the target, and records the positions for
  the values arriving. Previously the result reported `method = "unknown"` and,
  for `ps_trunc`, missing bounds, whatever it had been cast to.

* Two `ps_trim` objects are now compatible for combining only if they agree on
  the cutoffs their trimming settled on, matching the rule `ps_trunc` already
  applied to its bounds. Two objects trimmed with the same method and the same
  percentiles can still have been trimmed at different scores, and the result
  would have reported one object's cutoffs over the other object's units.
  Incompatible objects warn and combine as numeric, as other metadata mismatches
  do.

* Combining a `ps_trim` or a `ps_trunc` with an integer vector now gives a
  numeric result holding the propensity scores, rather than an integer vector of
  zeros and ones. Casting either class to integer now raises a lossy-cast error
  instead of rounding every score away, as `psw` already did.

* `ps_trim()` and `ps_trunc()` now follow one policy for a propensity score
  that arrives missing: it joins neither the retained nor the modified
  positions, `is_unit_trimmed()` and `is_unit_truncated()` report `FALSE` for
  it, and the value propagates as `NA`. A matrix row with a missing cell passes
  through `ps_trim()` unchanged and comes back missing throughout from
  `ps_trunc()`, which has to divide by the row sum to renormalize. The
  `"adaptive"`, `"pctl"`, `"cr"`, and `"optimal"` cutoffs are now worked out
  from the complete scores or rows, so each is the cutoff the same call would
  produce with the missing observations dropped; `"pref"` continues to center
  its preference scores on the proportion exposed across the whole sample, which
  is a fact about the exposure rather than about the scores. Previously
  `ps_trim()` recorded an arrived-missing score as one it had trimmed under
  `method = "ps"` and `method = "pref"`, `ps_trunc()` recorded a matrix row it
  could not renormalize as one it had pinned to a bound whenever an observed
  cell of that row fell outside them, a matrix row with a missing cell was
  blanked out by `ps_trim()` and its observed scores lost with it, and the
  `"adaptive"`, `"pctl"`, `"cr"`, and `"optimal"` methods either raised a bare
  base error or trimmed the whole sample.

* `ps_trim(method = "pref")`, `ps_trim(method = "cr")`, and
  `ps_trunc(method = "cr")` now refuse an exposure with missing values, with an
  error of class `propensity_missing_value_error`. Their cutoffs come from the
  exposure groups, and a unit that belongs to neither left every cutoff missing,
  which trimmed the whole sample without a word or recorded bounds that were
  never applied.

* `ps_trim()` now requires `lower` below `upper` for `method = "pctl"` and
  `method = "pref"`, and `ps_trunc()` for `method = "pctl"`, with the error of
  class `propensity_range_error` that `method = "ps"` already raised. Bounds the
  wrong way around describe an empty interval, so every unit fell outside it.
  Percentile bounds outside `[0, 1]` are refused with the same class, naming the
  valid range rather than passing the bare `quantile()` error about `probs`, an
  argument neither function takes.

* `ps_trunc(method = "cr")` now refuses exposure groups whose propensity score
  distributions do not overlap, with an error of class
  `propensity_no_overlap_error`. The bounds crossed, every score was pinned to
  the bound on the far side of the other, and the result reported a common range
  the groups did not have.

* `ps_trim()` and `ps_trunc()` now accept a propensity score matrix or data
  frame with `method` left at its default. The generic offers every method and
  the matrix methods matched against the two each supports, so the unevaluated
  default was compared against a set it could not belong to and a call that
  named no method was refused for naming the wrong one.

* A method the categorical path does not define is now refused with an error of
  class `propensity_method_error` naming the methods that path does support,
  rather than with the bare argument-matching error of class `rlang_error`
  raised before.

* `ps_trim()` now refuses `method = "optimal"` on a vector of propensity
  scores, with an error of class `propensity_wt_not_supported_error` pointing
  at the matrix or data frame input the method is defined for. Optimal trimming
  is defined over the rows of a propensity score matrix; on a vector it fell
  through to common-range trimming and recorded itself as `"optimal"`, so the
  result misreported what had been done to it.

* `ps_trim()` now refuses a categorical trimming threshold at or above `1/k`
  for `k` exposure levels, with an error of class `propensity_range_error`,
  matching `ps_trunc()`. Such a threshold cannot be met by every column of a row
  that sums to one, and the untrimmed scores were returned instead with the
  rejected threshold recorded in the metadata, reporting a trimming that never
  happened.

* `ps_trunc()` now refuses `method = "cr"` with no `.exposure`, with an error of
  class `propensity_missing_arg_error`. The common range method reached the
  binary transform with nothing to transform, and the caller was told that the
  exposure could not be converted rather than that it was never supplied. The
  matching refusals in `ps_trim()` name `.exposure`, the argument the function
  takes, rather than `exposure`, which it does not.

* `ps_trim()` and `ps_trunc()` document why their categorical thresholds
  differ: `ps_trim()` defaults to 0.1, following common-support trimming
  practice, and `ps_trunc()` to 0.01, a gentle winsorization that keeps every
  unit.

* The weight functions now accept a matrix `.propensity` on a binary exposure,
  reading it exactly as they read the equivalent data frame: the column holding
  the probability of the resolved focal level, chosen by the rules documented
  under `.propensity_col`. A matrix survived every check unchanged, because
  comparison and coercion are elementwise, and then failed where the weights
  were given their class with a dimensionality error against an argument the
  caller never named. Matrices on categorical exposures are unaffected.

* `.sigma` is now validated: it must be numeric, must hold either a single
  value or one per observation, and applies only to continuous exposures.
  Anything else is refused with an error of class `propensity_sigma_error`.
  `.sigma` was read straight into a normal density, so a value that was not a
  number reached the density as though it were one and a length that neither
  matched nor divided the exposure recycled into weights nobody asked for.
  `.sigma` also sits in the third position, so an exposure type supplied
  without a name arrived there and was absorbed without a word.

* The data frame methods now announce the detected exposure type once. Each
  resolved the type to choose a column and then handed the unresolved argument
  to the numeric method, which resolved it again, so a call that made one
  decision announced it twice.

* `wt_cens()` now names itself, rather than `wt_ate()`, in the deprecation
  warning for `.treated` and `.untreated` on its numeric method. The numeric
  method delegated the deprecated arguments to the ATE machinery, which raised
  the warning against itself.

* `wt_cens()` now supports binary and continuous exposures only. Remaining
  uncensored is a two-level event, and a categorical exposure was answered with
  ATE weights carrying the `"uncensored"` estimand. Naming a categorical
  exposure, or supplying one that detection reads as categorical, is now an
  error of class `propensity_wt_not_supported_error`.

* `.focal_level` and `.reference_level` must now hold a single level, on every
  route that reads them. Each is compared against the exposure elementwise, so
  more than one level recycled across the observations and sorted the units
  into groups nobody named: two focal levels alternated down a binary exposure,
  and two reference levels left every unit coded as reference. The refusal has
  class `propensity_focal_level_error`.

* An exposure type a weight function does not support is now refused with an
  error of class `propensity_wt_not_supported_error` whether it was named or
  detected, and the message lists the types that function does support. Naming
  one was previously reported as an unrecognized value.

* Stabilized categorical weights now refuse an exposure with no observed
  values. The default stabilizer is the share of units at each level, so with
  none observed every share was 0 / 0 and the weights came back missing
  everywhere. The refusal has class
  `propensity_stabilize_categorical_error`. Without stabilization such an
  exposure still returns missing weights, which is the answer to that call.

* The refusal of a focal or reference level the exposure never takes now
  reports an exposure with no observed values as such, instead of listing the
  levels present and rendering an empty sentence.

* The categorical range refusal now names `.propensity`, the argument the
  scores arrived in, matching the binary refusal.

* An `.exposure` with dimensions is now refused on the binary path, with an
  error of class `propensity_binary_transform_error` that names the shape it
  was given. The length rule reads `length()`, which counts cells for a matrix
  and columns for a data frame, so a shape of either kind could match the
  propensity scores and survive the 0 and 1 coding unchanged, because
  comparison and coercion are elementwise, and then fail where the weights are
  given their class, reporting a dimensionality error against an argument the
  caller never named. `ps_trim()`, `ps_trunc()`, and `ps_calibrate()` took the
  same exposure without complaint and read a matrix in storage order.

* `ps_trim()` and `ps_trunc()` now name the function you called when they
  refuse a propensity score matrix. The refusal was attributed to the frame the
  matrix method was dispatched from, which is your own function when either is
  called from one and no call at all when either is called at the top level.
  The messages themselves are unchanged.

* `ps_trim()` and `ps_trunc()` document the focal level contract against the
  argument they take. The text was inherited from the weight functions and
  named `.propensity`, which neither function has, and it carried a requirement
  that belongs to `wt_att()` and `wt_atu()`. Behavior is unchanged.

* The getting started example in `README.Rmd` fits its weighted outcome model
  with `quasibinomial()`, matching the documentation examples. A weighted
  `binomial()` fit warns about non-integer successes, and the warning was
  suppressed rather than avoided.

* The weight functions now warn when the column of a `.propensity` data frame
  chosen for the focal level is one of several carrying that name, naming the
  count and the column read. A data frame is allowed repeated names, through
  `check.names = FALSE` or a later assignment, and the match took the first of
  them: the choice was made on nothing the caller had expressed, between
  columns that may hold different numbers, and the column it landed on was as
  likely to be the wrong one. The warning has class
  `propensity_df_duplicate_column_warning`, and the column read is unchanged.
  Setting `.propensity_col` answers the question the warning asks and is silent.

* The binary path now resolves the default focal level among the values the
  exposure takes rather than among the levels a factor declares. A factor
  answers for every level it declares, and subsetting one without
  `droplevels()` leaves declared levels behind, so a level no observation held
  could sit second and be taken as focal: every unit was then coded as
  reference, and the weights described an exposure nobody had. `ps_trim()`,
  `ps_trunc()`, and `ps_calibrate()` read the same coding and were wrong in
  the same way. A `.propensity` data frame is matched to columns by these same
  levels, so a frame named for every level the exposure holds now matches by
  name, and the warning raised when no column matches no longer carries a
  `droplevels()` hint, which was a workaround for this defect.

* `.focal_level` and `.reference_level` are now checked against the values the
  exposure takes on the binary path, with an error of class
  `propensity_focal_level_error` that names the levels present. A level the
  exposure never takes sorted every unit into one group in silence: a focal
  level nobody holds left every unit reference, and a reference level nobody
  holds left every unit focal, so a misspelled level returned weights rather
  than a refusal. The categorical path already refused this. The check covers
  factor, character, 0/1, and logical exposures.

* Exposure levels are now counted from the values an exposure takes, so a
  missing value is no longer counted as a level of its own. A two-level
  exposure carrying missing values counted three, which took it off the binary
  path entirely: at a sample size where the categorical heuristic fired it was
  refused for holding too few levels, and at a smaller one it fell through to
  the continuous density weights without a word. Such an exposure is now coded
  from the levels it takes, with the missing values carried through to the
  weights. The `glm` methods resolve a named focal level for it as they do for
  any other exposure, so naming the level the fit treats as the reference now
  reads the complement of the fitted values there as well.

* Continuous exposure weights are now the ratio of two fully normalized normal
  densities, and `.sigma` sets the spread of the conditional density one
  observation at a time. Both densities were evaluated as standard normal
  ordinates at z-scores, which drops the `1/sigma` factor that makes a normal
  density integrate to one, and the conditional spread was the pooled residual
  standard deviation in every case, so `.sigma` was documented and accepted but
  never read: no value of it changed a single weight. The weight values change,
  and stabilized weights now sit closer to a mean of one.

* The `glm` methods of `wt_ate()` and `wt_cens()` no longer extract
  `influence(model)$sigma` for a continuous exposure. The conditional density
  uses the pooled residual standard deviation of the exposure around the fitted
  values unless `.sigma` names observation-level ones. The extraction had no
  effect on the result while `.sigma` was ignored, and now that `.sigma` is read
  it would silently make every model route a leave-one-out calculation, which is
  a different estimator from the one the numeric route computes from the same
  fitted values. Pass `.sigma = influence(model)$sigma` to ask for those
  standard deviations.

* `ipw()` accepts continuous weights built with the pooled default alone. Its
  stacked estimating equations carry a single pooled residual variance
  alongside the propensity score coefficients, so weights spread by
  observation-level standard deviations are a different function of the data
  with no counterpart in that system and cannot be reproduced at any parameter
  value. The weight consistency check refuses them, and its message now names
  `.sigma` as a cause. Estimates and standard errors for weights built the
  supported way are unchanged, since the normalization rescales those weights
  by a constant.

* `ps_calibrate()` documents the focal level contract its `.exposure` coding
  actually follows, matching the weight functions: every binary coding honors
  `.focal_level`, an exposure with no level named defaults to its higher level,
  and naming the other level reverses the coding, so `ps` must then hold the
  probability of the named level. The documentation said only that coding was
  determined automatically. Behavior is unchanged.

* Errors and warnings from the weight functions now name the function you
  called on every route into them. A rejection from `wt_ate()` and friends
  previously reported an internal frame, or no call at all, whenever the
  propensity score arrived as a data frame, as a fitted model, or as trimmed,
  truncated, or calibrated scores, and whenever the exposure was categorical.
  Even on a plain numeric propensity score, an exposure that cannot be coded
  0 and 1 was refused in the name of the internal helper that holds the weight
  formula. A categorical rejection reported the function you called it from,
  which named your own code rather than the weight function. The messages
  themselves are unchanged.

* The weight functions now reject a `call` argument that is neither a call nor
  an environment with an error of class `propensity_call_arg_error`. The
  argument names the frame a condition is attributed to, and the weight
  generics pass their arguments on to their methods, so a value of any other
  type was accepted in silence and then reported the condition system's own
  internals when a later check failed.

* A `stabilization_score` with dimensions is now rejected. A matrix holding one
  value per observation passed the length rule and was silently flattened into
  the order its values happened to be stored in.

* `wt_ate()` now reports a `.propensity` and `.exposure` of different lengths
  before it checks the stabilization score. A score matching `.propensity` was
  reported as the wrong length for the weights, which described the score
  rather than the mismatch that made it wrong.

* New `ps_tilt()` returns the tilting function h of the propensity score that
  defines an estimand's target population, with methods for a numeric
  propensity score (binary exposure) and for a matrix or data frame of
  per-level probabilities (categorical exposure). A weight is h divided by the
  propensity score of the exposure level a unit received, and an h-weighted mean
  standardizes counterfactual predictions to the same population, which for the
  `"atm"`, `"ato"`, and `"entropy"` estimands is the only route to the target:
  those populations are not subsets of the rows and cannot be reached by
  filtering. The weight functions and the `ipw()` estimating equations now read
  their tilt from `ps_tilt()` rather than each carrying their own copy. Weights,
  estimates, and standard errors are unchanged.

* `ipw()` now requires the propensity score model's response to be the exposure
  column itself for a continuous exposure, as it already did for a binary one,
  and rejects anything else with an error of class
  `propensity_ipw_response_error`. `ipw()` names the exposure by reading that
  left-hand side, so a response written as `cbind(a1, a2)` or as a
  transformation such as `log(a)` gives back several names where one was
  expected. Neither was caught: without `.data` the call stopped on an internal
  assertion about `.exposure_name` having length 2, and with `.data` it asked
  for a column named `"log"`, which cannot exist. A propensity score model
  written the same way but fit with `glm(family = gaussian())` reached the
  binary guard instead and was told it had a matrix response it does not have.
  Compute the transformation into its own column and fit the propensity score
  model on that column.

* The propensity score response error now describes the shape it found.
  A `cbind(successes, failures)` response keeps the matrix wording, whether the
  shape is read from the model frame or, for a fit made with `model = FALSE`
  whose data are gone, from the formula; a left-hand side that is a call but not
  a matrix, `factor(z) ~ x` for instance, is now reported as an expression
  rather than a single column, naming the expression. It previously described
  every such model as having a matrix response, which sent the writer of
  `factor(z) ~ x` looking for a matrix they had not written.

* The `se_method = "linearization"` outcome model rejection now distinguishes an
  outcome model that carries the exposure and more from one that does not carry
  the exposure at all. An intercept-only outcome model, `y ~ 1`, is adjusted for
  nothing and was nevertheless reported as "adjusted for terms beyond" the
  exposure; it now reads that the model does not include the exposure. The error
  class and the redirect to `se_method = "mestimation"` are unchanged.

* `ipw()` now rejects an estimand that names no estimand at all, whether it
  arrives as the `estimand` argument or as the estimand the weights record, with
  an error of class `propensity_ipw_estimand_error` listing the estimands it
  accepts. Neither source was checked against that list, and each was reported
  by whatever noticed it first: an argument checked against weights carrying an
  estimand came back as a disagreement between two estimands when only one of
  them exists, and one checked against plain numeric weights was not compared at
  all and reached the weighted means, where base R reported that `x` and `w` had
  different lengths. Weights built by hand with `psw()`, which records the
  estimand it is given, reached the same report on every path. A valid estimand
  that disagrees with the weights still reports the disagreement, and a valid
  estimand that the exposure type does not support still reports that.

* `ipw()` and `psw()` now require the estimand to be a single string. Membership
  was tested with `%in%`, which reads its left side through `as.character()`, so
  a value of another type matched the name it prints as and carried on in its
  own type: `estimand = list("ate")`, which is what single-bracket indexing of an
  options list returns, reached the weighted means as a list and stopped there in
  base R's terms, and `estimand = factor("att")` reached the tilt, which selects
  a branch with `switch()` and reads a factor as its integer level code, so a
  fit could complete under an estimand nobody asked for with nothing but base R's
  note about the coercion to say so. `psw()` reports the type with an error of
  class `propensity_estimand_type_error`; `ipw()` reports it with
  `propensity_ipw_estimand_error`, whichever of the two sources it came from.

* `ipw()` now warns whenever a fit reports standard errors that collapsed to
  essentially zero, with a warning of class
  `propensity_ipw_degenerate_se_warning` naming the affected effects. This was
  reported only where the solver had also said it could not pin the solution
  down, so a fit whose numbers came back cleanly and were still a false
  certainty said nothing: an outcome that does not vary within an exposure arm
  returned its contrast with a standard error of 1e-17, a p-value of zero, and
  an interval of no width. What is reported is a standard error of exactly zero,
  or one so small beside the estimate that the test statistic the two make is
  larger than double precision can carry any information at; the units the
  outcome is measured in do not enter, so a healthy fit of an outcome recorded
  on a scale a thousand million times finer is not reported. Both standard error
  methods now report it, once per fit. The estimates are unchanged and stay in
  the output.

* Two errors about exposures that are not binary, one per standard error method,
  now carry the `propensity_ipw_exposure_error` class the other exposure guards
  use; the mismatched-length and estimand errors from `ipw()` now carry
  `propensity_length_error` and `propensity_ipw_estimand_error`. All four
  previously arrived as the bare `propensity_error` every condition in the
  package shares and could not be caught apart from one another.

* `ipw()` now rejects a propensity score model or an outcome model whose design
  has columns the fit could not tell apart, with an error of class
  `propensity_ipw_rank_error` naming the model and the columns it has no
  coefficient for. `lm()` and `glm()` drop such a column from the fit and record
  `NA` for its coefficient rather than failing, while the design keeps the
  column, and `ipw()` multiplies the two together position by position. The `NA`
  then propagated: from the propensity score model it reached the rebuilt
  propensity scores and stopped the call with base R's "missing value where
  TRUE/FALSE needed" under `se_method = "mestimation"` and with an exactly
  singular system from LAPACK under `se_method = "linearization"`, and for a
  continuous exposure it was reported as a disagreement about the estimand the
  weights were built for. From the outcome model it reached the M-estimator's
  starting values, which reported that `stacked_equations` returned non-finite
  values at `init`. A column duplicated exactly and a column duplicated to
  within the tolerance the fit pivots at are both dropped by the fit and both
  reported here; a merely correlated column that the fit kept a coefficient for
  is untouched. The `NA` is the whole signal, so the rejection covers the models
  that record one, which are those fit with `lm()` and `glm()`. A categorical
  propensity score model fit with `nnet::multinom()` optimizes rather than
  pivots and returns finite coefficients for a dependent column, so its design
  and its coefficient matrix still agree and it is not rejected; a duplicated
  column there surfaces through the solver warning below, which asks you to
  check both designs for columns that duplicate one another.

* `ipw()` now warns, with a warning of class `propensity_ipw_solver_warning`,
  when the estimating equations behind the estimates have no unique root at the
  values the solver returned, naming the effects whose standard errors came back
  degenerate and saying that those standard errors are not meaningful. A
  binomial outcome with no events in one exposure arm, for instance, reported a
  log odds ratio of 20.4 with a standard error of 6e-36, a p-value of zero, and
  an interval of no width, and the only signal was a warning from the solver
  about rank deficiency in the estimating equations, in vocabulary from neither
  `ipw()` nor the models you fit. That warning is now replaced by this one. The
  estimates themselves are unchanged and stay in the output: they are the
  g-computation point estimates, and it is the inference around them that is
  empty.

* `ipw()` now reports a fit whose variance could not be built at all, with an
  error of class `propensity_ipw_variance_error`. When the derivatives of the
  stacked estimating equations come back non-finite there is nothing to invert,
  and this was said twice in vocabulary from neither `ipw()` nor the models you
  fit: a warning that a bread matrix contained `NA` values, raised as the
  variance was abandoned, and then, as the standard errors were read, an error
  that the fit had no variance to compute inference from, which discarded the
  estimates on the way past. Both are now replaced by one report, once per fit.
  Where the fit carries a structure that explains it, the report names it: an
  exposure arm or level in which the outcome never varies, which leaves the
  outcome model with no finite fit inside it, is named along with what the
  outcome does there. A binomial outcome that is an event for every exposed
  observation, and a categorical exposure with no events in one level, are both
  reported this way. Where no arm is degenerate, as for a design the fits kept
  every coefficient for but that all but repeats itself, the report says that
  the equations have no finite derivative at the solution and where to look.
  This is an error rather than a warning, as it was before: the standard errors,
  the intervals, and the p-values all come from the variance, so there is
  nothing to return.

* `ipw()` now rejects an outcome model that reads the exposure through a call
  alongside anything else, such as `y ~ z + x1 + I(z * x1)`, when `.data` is not
  supplied, with an error of class `propensity_ipw_exposure_error` naming the
  term and directing you to supply `.data`. Without `.data` the counterfactual
  designs come from the outcome model's own model frame, and a term recorded
  under a call keeps its fitted values whatever exposure value is written beside
  it, so the interaction was dropped from the marginal means with nothing
  signaled and the estimates moved away from the ones the same model gives with
  `.data`, which are also the ones `y ~ z * x1` gives on either route. What
  decides this is the variable the model frame holds rather than the term it
  builds, so a call reading the exposure is rejected wherever it sits, including
  inside an interaction such as `y ~ z + x2 + I(z * x1):x2`, which was the same
  silent movement. An interaction written with `*` or `:` over plain columns,
  such as `y ~ z * x1`, holds those columns in the frame and forms the product
  from the one being set, and is unaffected on both routes.

* `ipw()` now rejects a `.data` column supplied as a logical vector where the
  model recorded a factor, an ordered factor, or a character vector, with the
  same `propensity_ipw_data_error` the other type mismatches raise, naming the
  column and both types. A logical column carries no levels for the rebuild to
  re-level against the fit, so `model.matrix()` codes it on its own terms,
  `FALSE` first. Against a two-level fitted factor the widths agreed and nothing
  noticed: a logical whose `TRUE` marked the fit's second level reproduced the
  fitted design by coincidence, and one whose `TRUE` marked its first level
  silently swapped the two levels' coefficients and moved the estimates. The
  rejection is keyed on the type rather than on the values, because a logical
  column carries nothing that says which level its `TRUE` was meant to name. A
  column fit as a factor may still be supplied as character or as an ordered
  factor.

* `ipw()` now rejects a `.data` column supplied as a factor, an ordered factor,
  or a character vector where the model recorded a logical vector, with the same
  `propensity_ipw_data_error` and the same reasoning read the other way. A fit
  that recorded a logical column recorded no levels for it either, so it took
  `FALSE` as the reference, while each of those three brings a reference of its
  own: the widths agree, and which level the fitted coefficient is paired with
  is then decided by whatever the column declares or sorts to. A factor
  declaring `TRUE` before `FALSE` reproduced the fitted design reversed and
  moved the estimates with nothing signaled, and a character column matched only
  because `"FALSE"` sorts before `"TRUE"`. Supply the column as the logical
  vector the models were fit on, or refit them on the factor.

* `ipw()` now reports missing values in the `.data` columns the models read,
  with an error of class `propensity_ipw_data_error` naming the columns and the
  number of incomplete rows. Every design is rebuilt with `model.frame()`, which
  drops a row missing a value in any column it reads, while the weights, the
  exposure, and the outcome values keep every row of `.data`. The two were then
  recycled against each other: base R warned twice that one object was not a
  multiple of the other, and the mismatch surfaced as weights that failed their
  consistency check, whose message is about how the weights were built rather
  than about the data that was passed. Models fit on data with missing values
  dropped those rows themselves, so calling `ipw()` without `.data` is
  unaffected.

* `ipw()` now reports a `.data` outcome column supplied as a factor that
  declares no levels, with an error of class `propensity_ipw_data_error`. Such a
  column holds nothing but missing values and has no first level for the check
  that compares it against the level the outcome model was fit to treat as the
  failure, which stopped as a subscript out of bounds.

* `ipw()` now rebuilds an outcome model that reads the exposure through a term
  that has to work its levels out from the values it sees, such as
  `y ~ factor(z) + x1` fit on a numerically coded exposure, when `.data` is
  supplied, and rejects it when it is not. `ipw()` estimates the marginal means
  by setting the exposure column to one value at a time and rebuilding the
  outcome design. With `.data` the term is re-evaluated at each value and put
  back on the levels the model was fit with, which gives the estimates a model
  fit on a column converted to a factor beforehand gives; it previously failed
  inside base R as "contrasts can be applied only to factors with 2 or more
  levels". Without `.data` the design comes from the outcome model's own model
  frame, where the term is held at the values it was fit on and the
  counterfactual value is ignored, so every level would be given the fitted
  design and every contrast between them would be zero; that route now raises an
  error of class `propensity_ipw_exposure_error` naming the term and directing
  you to supply `.data`, where it previously reported that the exposure was
  missing from the model frame. With `.data`, a term that recomputes at the
  value being set, such as `I(z^2)` or `cut(z, breaks)`, is unaffected and keeps
  working. Without `.data` nothing recomputes at all, so such a term is now
  rejected on that route as well.

* Every design `ipw()` rebuilds from `.data` now uses the levels the model
  recorded rather than the levels the supplied column declares. The propensity
  score design already did, but the counterfactual outcome designs did not, so a
  factor covariate re-leveled after fitting, or stored as a character vector
  whose fitted levels were not in alphabetical order, paired each level with
  another level's coefficient and moved the estimates with nothing signaled. A
  column that declares a level no observation carries is accepted and rebuilds
  the fitted design, since the fits drop unused levels themselves; a column
  holding a value the fit never saw now raises an error of class
  `propensity_ipw_data_error` naming the column and the value, where it
  previously failed inside base R as a report that a factor had new levels.

* `ipw()` now checks every `.data` column against the type its model was fit on,
  in both directions, and raises an error of class `propensity_ipw_data_error`
  naming the column and both types. Only a column fit as a factor and supplied
  as numeric was checked before. A numeric column supplied as a factor or as a
  logical failed as a non-conformable multiply, or, when the two codings took
  the same number of design columns, rebuilt a design of zeros and ones the
  model was never fit on and moved the estimates with nothing signaled. A column
  fit as a factor may still be supplied as character or as an ordered factor,
  which rebuild the fitted design.

* `ipw()` now rejects a character column supplied under the outcome model's
  response name, and a response supplied as a factor where the model was fit on
  a numeric or logical column, with an error of class
  `propensity_ipw_data_error`. A model can never be fit on a character response,
  so such a column is always a `.data` mismatch: it previously reached the value
  conversion, which coerced it with an unclassed "NAs introduced by coercion"
  warning and then failed inside the M-estimator's own vocabulary, or, on the
  linearization path, left the standard errors non-finite. A factor response
  whose first level is not the level the model was fit to treat as the failure
  is rejected with the same class: `ipw()` codes a factor outcome as an
  indicator for its non-first levels, so a response re-leveled the other way
  described the opposite outcome. It reversed the sign of every estimate on the
  M-estimation path, and left the point estimates alone and moved the standard
  errors on the linearization path.

* `ipw()` now rejects a binary exposure column supplied as character where the
  outcome model holds the exposure as a factor. The counterfactual rebuild sets
  the column to one value at a time, and a character column carries no levels of
  its own, so it held the single level it was set to and failed inside base R.

* `ipw()` now rejects a binary factor exposure supplied in `.data` on the fitted
  levels in a different order, with an error of class
  `propensity_ipw_data_error` naming the column. `ipw()` treats the second level
  of a binary exposure as the exposed group, so a re-leveled column contrasts
  the levels the other way round. This did error, at the check that recomputes
  the weights from the propensity score model, whose message is about how the
  weights were built rather than about the data that was passed. A level no observation carries is accepted wherever
  it sits, as it is for a covariate. A categorical exposure is unaffected: it is
  resolved against the propensity score model's own level order before anything
  reads it.

* `ipw()` now checks that the counterfactual outcome designs rebuilt from
  `.data` are as wide as the design the outcome model was fit to, the mirror of
  the check the propensity design already gets. A `.data` column whose levels
  differ from the fitting data previously reached the marginal means as a raw
  non-conformable multiply.

* `psw()` and the weight functions that take a `stabilization_score` now check
  it where it is recorded. A score must be numeric, positive, and finite, and
  must hold either a single value, which scales every weight, or one value per
  observation, which scales each weight by its own; anything else raises an
  error of class `propensity_stabilization_score_error`. A score of the wrong
  length was previously recycled into the weights without comment, so two values
  supplied for four observations stabilized half the weights on the wrong
  multiplier, and a zero, negative, or missing score produced weights that were
  zero, negative, or missing.

* A `stabilization_score` is now stored as a double whatever numeric type it
  arrives as. The metadata two sets of weights carry is compared for identity,
  so a score written `1L` and a score written `1` read as two different scores:
  combining the weights or operating on them warned about a conflict and dropped
  the score from the result. They now describe the same stabilization and carry
  through.

* `ps_trim()` and `ps_trunc()` objects no longer carry a record of which units
  were modified onto a result holding a different number of observations. An
  operation that changed the number of observations without going through `[`
  kept the original positions, so `vctrs::vec_slice()`, and with it
  `dplyr::filter()`, joins, and group-wise summaries, returned a short column
  whose record still described the rows before the filter: `is_unit_trimmed()`
  on a two-element result answered with a four-element vector naming rows that
  were no longer there. The record is now dropped when it cannot be re-indexed,
  with a warning of class `propensity_trim_record_warning` or
  `propensity_trunc_record_warning`, and `is_unit_trimmed()`,
  `is_unit_truncated()`, and `ps_refit()` raise an error of class
  `propensity_missing_meta_error` on an object whose record is absent or was
  written for a different number of observations, rather than answer from stale
  positions. The values, the class, and the method and its cutoffs or bounds are
  untouched by the drop, and a record that covers the object it is on is carried
  exactly as before, so arithmetic, renaming, and subassignment at the same
  length are unaffected.

  A record is checked against the object it is on by counting observations, so
  an operation that keeps the count and reorders them is neither re-indexed nor
  refused. `vctrs::vec_slice(x, 5:1)` and `dplyr::arrange()` return a column
  whose positions still name the order the observations arrived in, and
  `is_unit_trimmed()`, `is_unit_truncated()`, and `ps_refit()` answer from
  those positions: the first two name the wrong units and the third refits on
  the wrong rows. Subsetting with `[` is handed the subscript and re-indexes, so
  reorder with `[`, or order the propensity scores before trimming or
  truncating them.

* Combining `ps_trim()` or `ps_trunc()` objects with `c()` now drops the record
  of which units were modified. Concatenation appends one set of observations to
  another, so the positions either record names would describe units from the
  other input. The combined result previously worked its record out from the
  values, which reported a propensity score that arrived missing as one
  `ps_trim()` had removed, and a score that arrived equal to a bound as one
  `ps_trunc()` had winsorized.

* Subsetting a `ps_trim()` or `ps_trunc()` object with a subscript that names a
  position more than once now reports that unit at every place it holds.
  `x[c(1, 1, 2)]` recorded only the first copy of position 1 and left position 2
  in neither the trimmed nor the retained set.

* `unique()` and `rep()` on `ps_trim()` and `ps_trunc()` objects now return a
  record written for the result. Both previously carried the record from the
  object they were called on, so `is_unit_trimmed()` on the three unique values
  of a five-element vector answered with five elements naming positions the
  result does not hold. Each knows which position every value it returns came
  from, and now re-indexes the record onto it. `unique()` also rejects a
  non-default `incomparables` with an error of class
  `propensity_unsupported_arg_error` rather than compare the values it names
  along with the rest.

* Subsetting a `ps_trim()` or `ps_trunc()` matrix by row name, or by a logical
  subscript holding `NA`, now maps the record onto the rows the subset returns.
  Row names were compared against integer positions and matched nothing, so
  every row of the result was recorded as neither trimmed nor retained and
  `is_unit_trimmed()` reported all of them as retained. A logical subscript
  holding `NA` left that row out of the record while base R returned it as a row
  of `NA`, shifting every position after it, and an integer subscript holding
  `NA` failed with base R's message about a missing value where `TRUE`/`FALSE`
  is needed. A row taken by an `NA` subscript names no observation and so is
  reported as neither trimmed nor retained, which is what `[` on a vector
  already did.

* `vctrs::vec_cast()` to a `ps_trunc()` object now returns the documented
  metadata structure. The result carried the method and bounds but no
  `truncated_idx`, unlike the `ps_trim()` equivalents.

* An arithmetic operation between two `psw` objects now keeps the metadata the
  two operands agree on. Multiplying weights by weights, as in combining weights
  against confounding with weights against censoring, previously rebuilt the
  result from six fields alone and discarded the record of which units a
  modified propensity score trimmed, truncated, or calibrated, the attributes
  describing a categorical exposure, and the stabilization score. The result was
  visibly inconsistent: `w * 1` kept the trimming record and `w * psw(1)` did
  not, so `is_unit_trimmed()` answered for the first and refused for the second.
  Each of those attributes now carries when only one operand records it and when
  both record the same value, so weights marked as stabilized keep the score
  they were stabilized on rather than reading as stabilized by the default. A
  score is carried only when the result is stabilized, which takes both
  operands; a result that is not stabilized drops the score without comment,
  since the score describes weights it no longer stands for. Operands that
  record an attribute differently drop it, since neither value describes the
  result, and a single warning of class
  `propensity_metadata_conflict_warning` names every attribute dropped that way.
  The rule is applied per operation, so an attribute one operation drops for a
  disagreement can be carried again by a later operation whose operands agree.
  The six fields merge as before: two different estimands are pasted together,
  the result is stabilized only when both operands are, and it is marked as
  trimmed, truncated, or calibrated when either operand is.

* An arithmetic operation between a `psw` that names an estimand and one that
  names none now records the estimand the one operand names. The two labels were
  pasted together as though both named something, so
  `psw(w, estimand = "ate") * psw(v)` reported an estimand of `"ate, "`.

* Combining `psw` objects with `c()` now keeps the attributes describing a
  categorical exposure when the inputs agree on them. These name the exposure
  levels rather than the units, so they mean the same thing at the combined
  length; inputs that name different levels drop them with the same warning
  arithmetic raises, once for each attribute dropped, whatever order the inputs
  were given in. The records left by a modified propensity score are still
  dropped by `c()`, whether or not the inputs agree on them: concatenation
  appends one set of observations to another, so the positions a record names
  would describe units from the other input.

* An unrecognized `exposure_type` given to any of the weight functions is now
  attributed to the function the caller named rather than to the internal helper
  that validates the argument. `wt_ate(ps, z, exposure_type = "wrong")`
  previously raised its error under `match_exposure_type()`, an unexported name
  that appears in no documentation and that gives no hint of which call in a
  pipeline was at fault. The rejection now reports `wt_ate()`, `wt_att()`,
  `wt_cens()`, and the rest, on every route into them: a numeric vector, a data
  frame of propensity scores, a fitted model, and a trimmed, truncated, or
  calibrated propensity score alike. The message is unchanged, and the valid set
  each function reports still reflects the exposure types that function
  supports.

* `ipw(se_method = "mestimation")` now warns with a class of its own when a
  reported effect is undefined at the marginal means the solver reached. The
  means are free parameters of the estimating equations, so a fit whose exposure
  arm is all events, or all non-events, sends one past the range its log risk
  ratio or log odds ratio is defined on. Such a fit previously emitted base R's
  unclassed `NaNs produced`, which named neither the effect nor the contrast and
  which no handler could select on. The warning now carries class
  `propensity_ipw_contrast_warning`, names the effect and, for a categorical
  exposure, the contrast it belongs to, and reports once per effect per contrast
  however often the solver revisits those means. The estimates, the
  standard errors, and the convergence behavior of every fit are unchanged. The
  warning follows the solver rather than the reported estimates, so it is raised
  whenever the path taken to the root leaves the range, finite differences for
  the bread included; an effect defined at the root but within a finite
  difference of the boundary is reported too.

* The error `ipw()` raises when the weights it recomputes from `wt_mod` disagree
  with the ones `outcome_mod` was fit with now reports that comparison rather
  than declaring the weights inconsistent, and lists the causes as
  possibilities. Weights that are exactly right reach this error: a `.data` whose
  values differ from the data the models were fit to moves the recomputed
  weights on its own, and so does trimming, truncating, or normalizing the
  weights after `wt_mod` was fit. The message now names those two causes
  alongside a mismatched estimand or focal level, and the refit remedy is
  offered for the case where the weights are at fault instead of reading as the
  single next step. The error class, the tolerance the comparison uses, and the
  per-exposure-type focal hints are unchanged.

* `ipw()` now accepts an outcome model whose response is a transformation of a
  column, such as `log(y)` or `scale(y)`, when `.data` is supplied. The response
  is evaluated against `.data` rather than looked up by name, so it is computed
  the way the fit computed it, and the transformed and untransformed routes give
  the same estimates. Such a call previously asked for a column named after the
  function wrapping the response, which no correct `.data` could supply. The
  column report is fixed along with it: a `.data` missing the response now names
  the variable the response reads, on the binary, categorical, and continuous
  paths and on both standard error methods. A matrix response such as
  `cbind(successes, failures)` is still rejected for its shape.

* `ipw()` now rejects a `.data` column supplied as numeric where either model
  recorded it as a factor, with an error of class `propensity_ipw_data_error`
  naming the column, the levels the fit recorded for it, and what to supply
  instead. Designs rebuilt from `.data` are rebuilt under the coding the fit
  recorded, so such a column previously reached `stats::model.matrix()` with a
  contrast specification that could not be applied to it and failed inside base
  R, as an unclassed warning that the variable is not a factor followed by an
  error that contrasts apply only to factors. Neither mentioned `.data`, and
  neither said what to supply. Both the propensity score model and the outcome
  model are checked, on both standard error methods. The exposure is exempt,
  since the counterfactual designs set that column themselves, and a character
  column is still accepted: `stats::model.frame()` re-levels it against the
  recorded levels and rebuilds the design the model was fit to.

* `ipw(se_method = "linearization")` now rejects a propensity score model that
  separates the exposure, as the M-estimation path already did, with the same
  error of class `propensity_ipw_separation_error` naming how many observations
  are affected. Both paths now refuse a fit at the same threshold: the fitted
  linear predictors, put through the link's inverse, give a probability of
  exactly zero or one for at least one observation, whose weight is then
  undefined. The linearization path previously accepted such a fit and returned
  an estimate with a small standard error beside it, because its propensity
  scores come from `predict()`, whose inverse link is bounded away from zero and
  one; every weight was therefore finite and nothing downstream failed. Fits
  short of saturation are unaffected. A propensity score model fit with a link
  the estimating equations cannot invert, such as `cauchit`, is not checked for
  saturation at all; such a model is still rejected for its link, on both
  standard error methods and with the same error as before.

* The `ps_link` argument of `ipw()` is deprecated. Both standard error methods
  accept only the link the propensity score model was fit with, so the argument
  can only restate what `wt_mod` already supplies, and supplying it now warns.
  The warning is advisory and nothing else changes: the matching value is still
  accepted and gives the same result as omitting the argument, and a value that
  names a different link is still rejected. A multinomial or continuous
  propensity score model, which has no link to override, still rejects the
  argument outright rather than deprecating it.

* A data frame of propensity scores now resolves the column holding the
  probability of the focal level by name. When the frame has a column named for
  every level of a binary `.exposure`, the weight functions read the column
  named for the level that `.focal_level` or `.reference_level` resolves to,
  wherever that column sits. Selection previously fixed on the second column, so
  a caller who named the level whose column came first, or who arranged the
  columns the other way round, silently received weights built on the
  probability of the wrong level. A frame whose names do not cover the levels
  still falls back to the second column, now with a warning of class
  `propensity_df_column_warning` naming the column used when a level was named
  and could not be matched. Supplying `.propensity_col` overrides all of this
  and never warns, and calls that name no level are unchanged.

* `.focal_level` and `.reference_level` are now honored for 0/1 numeric and
  logical exposures, which previously ignored them. A caller who named `0` as
  the focal level of a 0/1 exposure, or `FALSE` for a logical one, silently
  received weights built with the opposite level as focal: for `wt_att()` and
  `wt_atu()` those are the weights for the other estimand, and for the
  remaining estimands they are the weights of a model fit the other way round.
  These codings now behave as a two-level factor already does, so `.propensity`
  must hold the probability of whichever level is named. The `glm` methods
  derive that probability themselves, subtracting the fitted values from one
  when the named level is not the response's success level. Calls that name no
  level are unchanged: the higher level, `1` or `TRUE`, is still focal. As a
  consequence, `ipw()` now rejects `"att"` or `"atu"` weights built on a 0/1 or
  logical exposure with the lower level named as focal, since its own binary
  path always treats the higher level as focal. Those weights previously
  recorded the higher level and passed the check, because they were in fact
  built for it.

* A 0/1 exposure now reaches the same weights whether it is stored as double or
  as integer. Integer storage was not recognized as a 0/1 coding, so it took the
  two-level fallback path, which announced the focal level it chose. The weights
  agreed, but the messages did not, and only the double form ignored the named
  levels.

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

  Whether a record covers the weights is decided by the number of observations
  the record was written for, which is the test a `ps_trim` applies to itself,
  so weights and the propensity scores they were built from answer alike. It
  was previously decided by adding the retained and the trimmed positions
  together, which falls short of the observations whenever a position belongs to
  neither: `x[c(1, NA, 2)]` gives three observations and two recorded positions,
  and `is_unit_trimmed()` answered on the `ps_trim` while refusing on weights
  built from it.

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
  when its length matches the data being cast, since it describes the
  prototype's observations rather than the incoming data; at any other nonzero
  length it is dropped silently and the result stays marked as stabilized. A
  cast to no observations keeps it, having brought no observations for it to
  contradict, and the result is itself a prototype.

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
  model. `se_method = "linearization"` refuses the same fits, described above.

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
