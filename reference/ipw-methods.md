# Inverse Probability Weighted Estimation

The `joint_wt_models` method estimates the effects of a joint
intervention on two treatments from the pair of fitted treatment models
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
records and a weighted outcome model that reads both treatments.
Standard errors are computed by M-estimation. Neither
`se_method = "linearization"` nor `se_method = "robust"` is available
here: neither path solves a stacked system, and the cell means and every
contrast built from them are parameters of one.

The only supported estimand is `"ate"`, which is what the product
weights
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
builds target.

Either treatment may be categorical rather than binary, fit with
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html), in
which case its block is the multinomial score that fit solves and the
surface reports the whole crossing its levels make. Such a component may
be stabilized like any other, its numerator estimated in a block of its
own; see **Joint exposures**.

The second treatment may be a dose, whatever the first one is, in which
case the surface is the marginal structural model's own coefficients
rather than the cells of a crossing, with one level contrast per
non-reference level of the first treatment and one dose slope at its
reference level. The dose model is read through the same registry the
single-treatment route reads, so an
[`stats::lm()`](https://rdrr.io/r/stats/lm.html), a gaussian
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) at an identity or a
log link, and a [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html)
fit with one of the psi functions MASS supplies are stacked, and the
dose component's weights may be any density ratio
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
builds. A dose model or a density the stacked system cannot
differentiate is refused rather than resampled, since this route has no
resampling method; see **Joint exposures**.

Either component may be stabilized, including on a numerator model the
caller fit and passed to
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)'s
`stabilize`. Each component's numerator is estimated in a block of its
own, so the standard errors account for it having been fitted; see
**Joint exposures** for what each block holds and for what a numerator
may condition on. A numerator the caller computed and passed as
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)'s
`stabilization_score` is read back off the product and held fixed,
contributing no block, as it is for a single treatment.

Every score this route stacks is unweighted, so a treatment model fit
with prior case weights is refused by the name it was recorded under in
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md).

This is the second of the two routes to a joint intervention. The other
declares the crossing with
[`causalgenerics::joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint_exposure.html)
and weights it with one multinomial model over the cells; see **Joint
exposures**. Both report the same surface, so the choice between them is
a modeling choice rather than a reporting one. Prefer this route when
the two treatments call for different adjustment sets, or when the
dependence of the second treatment on the first is what you want to
model directly: it weights through the sequential factorization
`f(A | L) f(E | A, L)`, so each treatment gets its own model and its own
covariates. Prefer the declared crossing when one model over the cells
is the natural specification.

The factorization is validated when the container is built rather than
here:
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
refuses a second model that does not condition on the first treatment,
and a pair that condition on each other. By the time a container reaches
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
those questions are settled.

The `multinom` method estimates causal effects for a categorical
exposure from a fitted
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
propensity score model and a weighted outcome model. Standard errors are
computed by M-estimation. Neither `se_method = "linearization"` nor
`se_method = "robust"` is available for a categorical exposure, and both
are refused with an error of class `propensity_method_error`: the
stacked system has a sandwich for every fit this exposure type accepts,
and the diagnostic sandwich describes two fitted cells that a
categorical result does not report.

For a K-level exposure, the counterfactual mean at each level is
reported first, under the effect label `"mean"`, and then the effects
for each non-reference level against the reference (first) factor level.
The `contrast` column names the level each mean row belongs to and the
pair of levels each effect row compares, as it does for a binary
exposure.

Two fits are refused with an error of class
`propensity_model_family_error`. A
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) fit to
two levels reports the single probability of a binary exposure rather
than one probability per level; fit such an exposure with
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) and
`family = binomial()` and pass it to the `glm` method instead. A fit to
a matrix response records no levels at all, since
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) reads a
matrix as counts; refit it with the exposure factor on the left-hand
side.

The `lm` method estimates the causal dose-response effect for a
continuous exposure from a fitted
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) (or gaussian-family
[`stats::glm()`](https://rdrr.io/r/stats/glm.html)) propensity score
model of the exposure and a weighted marginal structural outcome model.
A [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) fit with one of
the psi functions MASS supplies and an
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) of a gaussian
family reach this method by inheritance and are stacked at their own
scores, the additive fit with its smoothing parameters held at the
values it reports. A gaussian model is read through its link, which may
be identity or log; the other gaussian links, and every class the
stacked system has no score for, are refused, as are weights built with
a `"kernel"` density, which have no closed-form standard error. See
**Continuous propensity score models** in
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for what each refusal says, what to do about it, and what the additive
route conditions on. The only supported estimand is `"ate"`.

Standard errors are computed by M-estimation; neither the linearization
method nor the robust diagnostic is available for continuous exposures.
The package offers no resampling method, so a fit the stacked system
cannot write has no standard error here, and the refusal points at a
bootstrap of the whole pipeline written by hand. See **Continuous
propensity score models** in
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).

Every term of the marginal structural model that reads the exposure must
read the exposure and nothing else, and the reported effects are the
coefficients those terms contribute, labeled `"log(or)"` for a
logit-link outcome and `"log(rr)"` for a log-link one. A model with one
exposure coefficient keeps the eight-column contract with no contrast
column and labels its row `"slope"` at an identity link; a dose-response
curve such as `y ~ A + I(A^2)` reports one row per coefficient under
`"coef"`, gaining a contrast column that names them.

A curve is also the one shape of fit that has a single reading. An
exposure entering the outcome model through more than one design column,
whether written as `y ~ A + I(A^2)` or as a basis such as `poly(A, 2)`,
`splines::ns(A, 3)`, or `splines::bs(A, 3)`, has no coefficient that is
the effect of the exposure, because the dose response has a different
slope at every dose.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
records the conditional reading for such a fit and reports it once as a
message. Asking
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
itself for `effects = "marginal"` is refused with an error of class
`propensity_ipw_effects_error`; the result declares the conditional
reading as the only one it supports, so
[`causalgenerics::as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html),
[`stats::coef()`](https://rdrr.io/r/stats/coef.html),
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html),
[`stats::confint()`](https://rdrr.io/r/stats/confint.html),
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html),
and [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
refuse it from the result class instead, with an error of class
`causalgenerics_unsupported_reading_marginal`. Marginalizing the curve
over the observed doses is a separate estimand that this package does
not compute; use the marginaleffects package on the conditional result:
`avg_slopes()` for slopes, `avg_comparisons()` for contrasts, and
`avg_predictions()` for causal dose-response functions. A basis term
reading a covariate rather than the exposure is a covariate term however
many columns it expands to, so it contributes no row and leaves the
reading marginal.

A model frame records the term rather than the variables inside it, so
the frame of a fit whose exposure enters through
[`poly()`](https://rdrr.io/r/stats/poly.html) or a spline holds no
exposure column and `.data` must supply one. The propensity score design
is then rebuilt through the model's terms object, whose `predvars`
attribute records the basis the model was fit with, so the rebuilt
columns are the fitted ones.

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
is a bring-your-own-model (BYOM) inverse probability weighted estimator
for causal inference. You supply a fitted propensity score model and a
fitted weighted outcome model;
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
computes causal effect estimates with standard errors that correctly
account for the two-step estimation process.

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
supports binary, categorical, and continuous exposures. For a binary or
categorical exposure, a binary outcome returns the risk difference, log
risk ratio, and log odds ratio, and a continuous outcome returns the
difference in means. A continuous exposure is supplied through an
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or gaussian
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) propensity score
model with a weighted marginal structural outcome model whose link is
identity, logit, or log, and
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reports the single exposure coefficient of that model. A
[`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) fit with one of
the psi functions MASS supplies and an
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) fit of a gaussian
family are read as the equations their own coefficients solve; every
other subclass errors.

## Usage

``` r
# S3 method for class 'joint_wt_models'
ipw(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  .by = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization", "robust"),
  effects = c("marginal", "conditional")
)

# S3 method for class 'multinom'
ipw(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  .by = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization", "robust"),
  .focal_level = NULL,
  effects = c("marginal", "conditional")
)

# S3 method for class 'lm'
ipw(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  .by = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization", "robust"),
  effects = c("marginal", "conditional")
)

# S3 method for class 'glm'
ipw(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  .by = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization", "robust"),
  effects = c("marginal", "conditional")
)
```

## Arguments

- wt_mod:

  The weighting object: a fitted propensity score model that produced
  the weights, typically a logistic regression of class
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) with the exposure
  as the left-hand side of the formula.
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  is an S3 generic that dispatches on this object. The left-hand side
  must be the exposure column itself, for every exposure type: a matrix
  response such as `cbind(successes, failures)` and a transformed
  response such as `log(a)` or `factor(z)` both error, because
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  reads the exposure's name from that expression. Compute a transformed
  exposure into its own column and fit `wt_mod` on that column.

- outcome_mod:

  A fitted weighted outcome model of class
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html), with the outcome as
  the dependent variable and propensity score weights supplied via the
  `weights` argument. The weights should be created with a propensity
  weight function such as
  [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md).
  Supported outcome models are an
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html), a gaussian
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) with an identity
  link, and a binomial or quasibinomial
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) with a logit,
  probit, cloglog, log, or identity link; any other family (such as
  poisson or Gamma), an unsupported link (such as cauchit), or a
  non-identity gaussian link errors. A factor or logical outcome
  response is converted to `0`/`1` following glm's coding (the first
  factor level is the failure, every other level is the success).

- ...:

  These dots are for future extensions and must be empty. They separate
  the two models from the remaining arguments, which must therefore be
  supplied by name. Anything left in them, such as a `.data` argument
  given by position or a misspelled argument name, errors.

- .data:

  A data frame containing the exposure, outcome, and covariates. If
  `NULL` (the default),
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  extracts data from the model objects. Supply `.data` explicitly if the
  outcome model formula contains transformations that prevent extraction
  of the exposure variable from
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html), or
  if the propensity score model cannot reconstruct its design (for
  example, fit with `model = FALSE` with the fitting data no longer
  available);
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  then rebuilds the propensity design from `.data`. `.data` must be the
  data the models were fit to. It must have one row per observation the
  models were fit to, and each column must carry the type its model was
  fit on: a column fit as a factor may also be supplied as character or
  as an ordered factor, since every design rebuilt from `.data` is
  rebuilt under the levels the fit recorded, but a factor supplied as
  numeric, or a numeric supplied as a factor or a logical, errors. A row
  count or a column type that disagrees with the fitted models errors
  rather than rebuilding a design the models were not fit to.

  The levels a categorical column declares are its fit's to decide, so a
  column re-leveled after fitting is rebuilt in the fitted order, and
  one that declares a level no observation carries rebuilds the fitted
  design, since the fits drop unused levels themselves. A value the fit
  never saw has no coefficient to multiply and errors. The exposure
  column and the outcome model's response are read rather than rebuilt,
  and both are rejected when the order they declare contradicts the fit:
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  treats the second level of a binary exposure as the exposed group, and
  codes a factor response as an indicator for its non-first levels, so
  either one re-leveled the other way describes the opposite contrast.

  Supplying `.data` also decides how a term that transforms the exposure
  is handled. The counterfactual designs are built by setting the
  exposure column to each level in turn and rebuilding the outcome
  design from `.data`, so a term such as `as.numeric(a == "c")`,
  `I(z^2)`, or `factor(z)` is re-evaluated at the level being set rather
  than held at its observed value, and such a model is g-computation on
  the model as specified. Without `.data` there is nothing to
  re-evaluate: the designs come from the outcome model's own model
  frame, which holds each such term at the values it was fit on. A term
  that has to derive the exposure's levels from the values it sees, such
  as `factor(a)`, is therefore rejected on that route, with `.data`
  named as the remedy.

- .by:

  A single unquoted variable naming a modifier to report effects within,
  or `NULL` (the default) to report the effect for the sample as a whole
  and nothing else. The variable is looked up in `.data` when you supply
  one and in the outcome model's own model frame otherwise, and it must
  be a factor or a character vector with no missing values. The result
  then reports the effects it reports without `.by`, keyed by the group
  `"overall"`, followed by one set within each level of the modifier and
  one set contrasting each level against the first, and the `estimates`
  frame gains a `group` column naming the subgroup each row belongs to.
  See **Effect modification** for what those rows estimate and what they
  require.

- estimand:

  A character string specifying the causal estimand: one of `"ate"`,
  `"att"`, `"atu"`, `"atm"`, `"ato"`, or `"entropy"`. The available
  estimands depend on the exposure type: a binary or categorical
  exposure supports all six, while a continuous exposure supports only
  `"ate"`. For a categorical exposure, `"att"` and `"atu"` require a
  focal level (see `.focal_level`). See the **Estimands** section for
  the full support matrix. If `NULL`, the estimand is inferred from the
  weights in `outcome_mod`, which requires weights created with a
  propensity weight function such as
  [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md).

- ps_link:

  **\[deprecated\]** A character string specifying the link function
  used in the propensity score model: `"logit"`, `"probit"`, or
  `"cloglog"`. The link is always read from `wt_mod`, so `ps_link` can
  only restate what the fitted model already supplies and can simply be
  dropped. Supplying it warns.

  The argument applies only to a binomial
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) propensity score
  model on the binary path. A multinomial or continuous propensity score
  model has no link for `ps_link` to override, so a non-`NULL` value
  errors on both of those paths; leave it `NULL` there. `ps_link` cannot
  name a link other than the one `wt_mod` was fit with:
  `se_method = "linearization"` and `se_method = "robust"` error, and
  `se_method = "mestimation"` reports the resulting weights as
  inconsistent with the propensity score model.

- conf_level:

  Confidence level for intervals. Default is `0.95`.

- se_method:

  Method for standard error estimation. `"mestimation"` (the default)
  stacks the propensity score, outcome, and estimand estimating
  equations and returns the empirical sandwich variance.
  `"linearization"` uses the influence-function linearization of
  Kostouraki et al. (2024). Both account for the uncertainty of
  estimating the propensity scores.

  `"robust"` does not, and is a diagnostic rather than an inference. It
  is the linearization route with the correction for having estimated
  the propensity score dropped, so it reports the sandwich the weighted
  outcome model computes for itself: the delta-method reading of
  `sandwich::vcovHC(outcome_mod, type = "HC0")`, exactly. Ignoring that
  the weights were estimated generally understates the variance, and
  understates it most where the weights are least stable, so this is
  offered for reading beside one of the other two rather than in place
  of it. A result fit this way says so when printed and marks the tables
  it coerces to; see **Standard errors as a diagnostic** below.

  Which exposure types support which: `"mestimation"` supports every one
  of them; `"linearization"` and `"robust"` support a binary exposure
  alone. A binary or a categorical exposure has a sandwich for every fit
  it accepts. A few continuous fits have no sandwich, and the package
  computes no standard error for them: weights built with a `"kernel"`
  density are refused, and so are a
  [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) and an
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) fit whose score
  the stacked system cannot write. See **Continuous propensity score
  models** below.

  `"linearization"` and `"robust"` support only an outcome model of the
  exposure alone, fit with an intercept; a covariate-adjusted outcome
  model requires `"mestimation"`. For a binary or categorical exposure,
  every method requires an outcome model that can represent a baseline,
  so a numeric no-intercept coding such as `y ~ z - 1` errors on any of
  them. The standard errors of the conditional reading described under
  `effects` require `"mestimation"` as well: the linearization and
  robust routes store the outcome model unwrapped, so
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`stats::confint()`](https://rdrr.io/r/stats/confint.html), and
  [`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw.md)
  refuse `effects = "conditional"` there with a classed error rather
  than reporting the covariance the outcome model computed for itself.

- effects:

  The presentation mode the result records, either `"marginal"` (the
  default) or `"conditional"`. The marginal reading reports the
  population-averaged causal contrasts; the conditional reading reports
  the outcome model's coefficient surface. Both surfaces are computed
  whichever mode is named, so the argument settles which one the result
  presents and nothing else.
  [`causalgenerics::as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  and
  [`causalgenerics::as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  move a result between the two readings afterwards, and the accessors
  take an `effects` argument of their own for a single call.

  One kind of fit has a single reading rather than two. A continuous
  exposure entering `outcome_mod` through more than one design column,
  such as a [`poly()`](https://rdrr.io/r/stats/poly.html) or spline dose
  response, has no coefficient that is the effect of the exposure,
  because a curve has a different slope at every dose. Such a result
  records the conditional reading and says so once.
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  refuses `effects = "marginal"` itself with an error of class
  `propensity_ipw_effects_error`, and the result declares the
  conditional reading as the only one it supports, so the accessors
  refuse it from the result class with an error of class
  `causalgenerics_unsupported_reading_marginal`. See the continuous
  exposure section below.

  The covariance the conditional reading reports is the outcome block of
  the jointly estimated sandwich, which every route that stacks
  estimating equations attaches to the outcome model it stores:
  `se_method = "mestimation"` for a binary exposure, the categorical and
  joint routes, which run on M-estimation alone, and the continuous
  route under that method. A linearization or robust fit stacks no
  system and records no covariance of either kind, so
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`stats::confint()`](https://rdrr.io/r/stats/confint.html) error on
  its conditional reading, and printing it writes the coefficients under
  a note saying the standard errors are not reported rather than beside
  the ones the outcome model computed for itself.

- .focal_level:

  For the categorical (`multinom`) method with the `att` or `atu`
  estimand, the focal exposure level. If `NULL`, it is taken from the
  `focal_category` attribute of the weights in `outcome_mod`. An
  explicit value overrides the attribute.

## Value

Methods of
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
return an S3 object of class `ipw`. That result class is shared across
packages and its components, its
[`print()`](https://rdrr.io/r/base/print.html) method, and its
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
are documented at
[`causalgenerics::new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.html).
propensity's methods fill every component; three of them take
propensity-specific values:

- `estimand`:

  One of `"ate"`, `"att"`, `"atu"`, `"atm"`, `"ato"`, or `"entropy"`.

- `se_method`:

  One of `"mestimation"`, `"linearization"`, or `"robust"`.

- `fit`:

  The fitted M-estimator object when `se_method` is `"mestimation"`.
  Calling [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) or
  [`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
  on this object exposes the full stacked system of estimating
  equations, including the propensity score, outcome, and estimand
  parameters. It is `NULL` under `"linearization"` and `"robust"`,
  neither of which solves a system.

The result answers the standard model accessors, which causalgenerics
registers for the shared class:
[`stats::coef()`](https://rdrr.io/r/stats/coef.html) and
[`stats::confint()`](https://rdrr.io/r/stats/confint.html) for the
reported effects, [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html)
for their covariance,
[`stats::nobs()`](https://rdrr.io/r/stats/nobs.html) and
[`stats::df.residual()`](https://rdrr.io/r/stats/df.residual.html) for
the counts describing the fit, and
[`stats::weights()`](https://rdrr.io/r/stats/weights.html) for the
[`psw()`](https://r-causal.github.io/propensity/reference/psw.md) vector
the outcome model was fit with. Coefficients are named for the effect
measure and the contrast together, the contrast being the level a
counterfactual mean belongs to or the pair of levels an effect measure
compares, and for the subgroup as well where `.by` reports one row per
subgroup. Which surface
[`stats::coef()`](https://rdrr.io/r/stats/coef.html),
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`stats::confint()`](https://rdrr.io/r/stats/confint.html) report
follows the presentation mode the result records, described under
`effects` above.

Under `se_method = "mestimation"` the stored `wt_mod` and `outcome_mod`
are the models supplied, carrying their block of the joint sandwich.
Calling [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) on either
reports a covariance that accounts for the whole system having been
estimated from the same data, rather than the one the model's own
fitting routine reports. The models supplied are left as they were fit.
Linearization solves no such system, so its component models are stored
unchanged.

## Workflow

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
is designed around a three-step workflow:

1.  Fit a propensity score model (e.g., logistic regression of exposure
    on confounders).

2.  Calculate propensity score weights for your estimand (e.g.,
    [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md))
    and fit a weighted outcome model.

3.  Pass both models to
    [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
    to obtain causal effect estimates with correct standard errors.

You are responsible for specifying and fitting both models.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
handles the variance estimation.

## Estimands

The available estimands depend on the exposure type:

|             |        |                         |            |
|-------------|--------|-------------------------|------------|
| Estimand    | Binary | Categorical             | Continuous |
| `"ate"`     | yes    | yes                     | yes        |
| `"att"`     | yes    | yes (needs focal level) | no         |
| `"atu"`     | yes    | yes (needs focal level) | no         |
| `"atm"`     | yes    | yes                     | no         |
| `"ato"`     | yes    | yes                     | no         |
| `"entropy"` | yes    | yes                     | no         |

A continuous exposure supports only `"ate"`; any other request errors.
For a categorical exposure, `"att"` and `"atu"` require a focal level,
supplied through `.focal_level` or taken from the weights (see the
`multinom` method).

A binary exposure takes no focal level:
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
always treats the second sorted exposure level as the exposed group,
which is the second level of a factor, `TRUE` for a logical, and the
larger of two numeric values. Because `"att"` and `"atu"` weights are
mirror images of each other, weights built with
[`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
or
[`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
on the other level carry values that are correct for the reversed roles,
and
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
rejects them rather than correcting them as though the roles matched. To
target the other level, relevel the exposure so that it sorts second and
refit both models.

The estimand a weighted outcome model was fit for is not necessarily the
one a tool that averages its per-row predictions reports, so see
**Standardizing model predictions** in
[`ps_tilt()`](https://r-causal.github.io/propensity/reference/ps_tilt.md)
before taking such an average of a model fit for anything other than
`"ate"`.

## Effect measures

A binary or categorical exposure reports the counterfactual mean at each
exposure level before the contrasts built from those means. Those rows
carry the effect label `"mean"` and lead the table, one per level in
level order with the reference level first. Each is the outcome model's
prediction with every unit set to that level, averaged over the sample
under the estimand's tilt, so a binomial outcome reports the marginal
risk under each level and a continuous one reports the marginal mean.
Reading a risk difference beside the two risks it is a difference of is
what those rows are for: a difference of 0.05 is a different finding at
risks of 0.02 and 0.07 than at 0.50 and 0.55.

For a binary exposure, the reported measures depend on the outcome
model. A binary outcome
([`stats::glm()`](https://rdrr.io/r/stats/glm.html) with
`family = binomial()`) returns three measures:

- `rd`: Risk difference (marginal risk in exposed minus unexposed)

- `log(rr)`: Log risk ratio

- `log(or)`: Log odds ratio

A continuous outcome ([`stats::lm()`](https://rdrr.io/r/stats/lm.html)
or [`stats::glm()`](https://rdrr.io/r/stats/glm.html) with
`family = gaussian()`) returns only the difference in means (`diff`).

For a categorical exposure, the same measures are reported for each
non-reference level against the reference (first) factor level, so a
K-level exposure produces one block of measures per non-reference level.

The estimates table of either exposure carries a `contrast` column
naming what each row is about: the level a mean row belongs to, and the
pair of levels a contrast row compares, written `"1 vs 0"` or `"b vs a"`
in the levels the models were fit on. A binary exposure names its one
pair the same way a categorical one names each of its several, so the
two tables are read alike.

For a continuous exposure,
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reports the exposure coefficients of the weighted marginal structural
outcome model. A model with one exposure coefficient reports one row and
no `contrast` column, and a dose-response curve such as `y ~ A + I(A^2)`
reports one row per coefficient, with the `contrast` column naming the
coefficient each row belongs to.

The `effect` label names the scale: `log(or)` for a logit link and
`log(rr)` for a log link, since a coefficient there is a log odds ratio
or a log risk ratio per unit of whatever column it multiplies. At an
identity link the word depends on what the row claims. A model with one
exposure coefficient reports `slope`, because that coefficient is the
slope of the dose response everywhere. A model with several reports
`coef`, because a curve has a different slope at every dose and no one
of its coefficients is that slope; the row says which coefficient it is
and leaves the dose response to be built downstream from
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html).

Every term of the marginal structural model that reads the exposure must
read the exposure and nothing else. The requirement is on what a term
reads rather than on how it is written, so `A + I(A^2)` and `A + sin(A)`
are admitted and `A * x1` is not: a term reading a covariate as well
contributes a coefficient that depends on that covariate, and no row
could name the effect it stands for. Read the whole coefficient vector
from the returned `fit` object for a model this surface cannot report.

A basis matrix reads the exposure alone, so `poly(A, 2)`,
`splines::ns(A, 3)`, and `splines::bs(A, 3)` are admitted on the same
footing. One such term contributes one coefficient per basis column and
each of them is a row, reported under `coef` at an identity link like
any other multi-coefficient surface. The `contrast` column names the
coefficient rather than the term, a distinction a basis makes and
`A + I(A^2)` does not: the single term `poly(A, 2)` reports the
contrasts `poly(A, 2)1` and `poly(A, 2)2`. A basis reading a covariate
is a covariate term and contributes no row, and one reading both, such
as `poly(A, 2):x1`, is refused with every other mixed term.

A marginal structural model built on a basis of the exposure requires
`.data`. A model frame records the term rather than the variables inside
it, so the frame of such a fit carries the basis matrix and no column
holding the exposure, which is where
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reads the exposure from when `.data` is absent. Omitting it errors, and
the message names the exposure and directs to `.data`. Supplying it does
not refit the basis: the outcome design is read from the fit, and a
design rebuilt from `.data` is rebuilt through the fitted terms object,
whose `predvars` attribute records the basis that term was fit with.

Use
[`as.data.frame()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.html)
with `exponentiate = TRUE` to obtain risk ratios and odds ratios on
their natural scale.

## Effect modification

`.by` names a modifier and asks for the effect within each of its levels
alongside the effect for the sample as a whole. The reported rows come
in four blocks, in this order:

- the overall rows, unchanged from the fit without `.by`, the
  counterfactual means and then the contrasts, under the group label
  `"overall"`;

- the counterfactual mean at each exposure level within each level of
  the modifier, under a `"var = value"` label such as `"sex = female"`;

- one block of contrasts per level of the modifier, under the same
  labels;

- one block of contrasts per non-reference level against the reference
  level, which is the modifier's first level, under a
  `"var = value vs var = value"` label.

The last block reports no means of its own. Its rows compare the effects
in two strata, and a difference of two effects belongs to neither
stratum's population, so there is no level whose mean it would be.

A stratum's effect is g-computation restricted to that stratum: the
outcome model predicts each unit's outcome at each exposure level, and
the predictions are averaged over the units of the stratum, weighted by
the estimand's tilt. What a stratum row estimates is therefore the
effect in that stratum's own covariate distribution, and the last block
is a difference of two such effects. Those two populations differ
whenever the covariates are distributed differently across the strata,
so a difference between strata is a difference between two conditional
effects and not, on its own, evidence that the modifier changes the
effect. It is the estimand's population that the tilt fixes, and the
modifier's stratum that the indicator fixes, and the row reports the
effect where those two meet: an `att` fit reports the effect among the
exposed units of a stratum, not among all of its units. For a
categorical exposure the tilt an `att` fit standardizes over is the
focal level's own propensity score, so a stratum row there reports the
effect among the units of that stratum who look like the ones assigned
the focal level.

A categorical exposure already reports one contrast per non-reference
level, and `.by` crosses those contrasts with the strata: every block
holds each contrast on each reported measure, so a K-level exposure with
S strata reports the whole-sample block and then S + (S - 1) blocks of
(K - 1) contrasts, beside S blocks of K means. A row is named by three
things, the measure, the contrast, and the subgroup, and the `estimates`
frame carries `contrast` and `group` in that order after `effect`.
Within a contrast block the ordering is the one an ungrouped fit already
uses, contrast-major and measure-minor, so a contrast's measures sit
together; the blocks themselves run group-major.

The stratum rows and the stratum contrasts report the risk difference
and the log risk ratio for a binary outcome, and the difference in means
for a continuous one. They report no odds ratio. An odds ratio is
noncollapsible, so a difference of two of them moves with the outcome
distribution in each stratum whether or not the effect does, and
reporting it as effect modification would invite exactly the reading it
does not support. The overall rows keep all three measures.

Three requirements apply. The modifier must be a factor or a character
vector with no missing values: a missing value names no subgroup, and a
continuous variable names no fixed set of them, so cut it into groups
first. A level no unit carries is dropped rather than refused, since the
fits drop it too and an empty stratum has no mean to standardize. Every
one of the modifier's remaining levels must hold every exposure level,
since a stratum in which nobody took one of them identifies no contrast
against it and the outcome model would extrapolate one from the strata
that do; the remedy is a coarser modifier rather than a refit. Finally,
`.by` requires `se_method = "mestimation"` and a binary or categorical
exposure: the stratum means and their contrasts are parameters of the
stacked estimating equations, which is what puts them in the same
sandwich as everything else, and a continuous exposure reports its
marginal structural model's own coefficient rather than a set of
standardized means to modify. Each of those errors.

An outcome model with no term reading the exposure and the modifier
together warns rather than errors. The stratum effects are g-computation
on the model as specified, so such a model forces the same effect on
every stratum up to the covariate distribution each one has, and the
modification it reports is a property of the model rather than of the
data. The fit still returns.

Stabilizing the weights on the modifier is worth doing here. A
stabilized weight divides by the propensity score and multiplies by a
numerator, and the numerator may condition on any variable the outcome
model the effects are computed from already reads, without changing what
is estimated. The modifier is such a variable by construction, so a
numerator conditioning on it leaves the estimator consistent for the
same stratum effects and changes only how variable the weights are. That
is where it pays: the part of the exposure's variation the modifier
explains cancels between the numerator and the propensity score, and the
weights tighten by more than the default numerator tightens them.

The converse is the constraint, and it is the reason this is a statement
about the model rather than a free choice. A numerator conditioning on a
variable the outcome model does not read changes the estimand rather
than the precision, because the pseudo-population it builds is one in
which that variable still predicts the exposure, and confounding by it
remains (Robins et al. 2000; Cole and Hernán 2008; Hernán and Robins
2020, Chapter 12). Condition the numerator on the modifier and on
nothing the model does not hold. The same rule governs an
exposure-by-modifier interaction reported without `.by`, since what
matters is what the reported model reads and not which argument asked
for the rows.

Build the numerator with
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)'s
`stabilization_score`, evaluated at the exposure each unit actually
took; see **Stabilization** in
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for the recipe. The sandwich treats a supplied score as a known
constant, so the standard errors do not account for the numerator model
having been fitted. The default stabilizer is estimated in the stacked
system instead, which is one parameter wider as a result. There is a
third way: pass the fitted numerator model to `stabilize` itself, a
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) of a binary exposure
or an [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
of a categorical one, and its own equations join the stack, which is
what a supplied score has none of. That accounting has something to say
only where the reported model is not saturated in the variables the
numerator reads; a numerator on a modifier the reported model interacts
with fully divides out of every cell, so the estimates and the standard
errors are the ones the unstabilized weights give.

A stabilized fit reporting effects within the strata of a modifier, and
carrying the default numerator, warns with class
`propensity_ipw_by_stabilizer_warning` and returns: the default
numerator conditions on nothing, which is consistent for the same
stratum effects and leaves precision on the table. A numerator model the
caller supplied is read the same way, since a fitted model carries the
terms it read: one that never reads the modifier warns with the same
class, naming the terms it was built on, and one that reads the
modifier, on its own or through a transformation of it, is silent. A
`stabilization_score` says nothing either way, being a vector of numbers
that does not carry what it conditions on.

Every reported row is a parameter of one stacked system solved on one
sample, so the covariance couples them:
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) reports nonzero
entries between the overall rows, the stratum rows, and the stratum
contrasts, and the standard error of a stratum contrast accounts for the
two strata having been estimated together. Fitting each stratum on its
own subset would report the same stratum effects with no covariance
between them, leaving the difference between two strata untestable.

## Joint exposures

[`causalgenerics::joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint_exposure.html)
crosses two discrete treatments into one categorical exposure and
records the crossing on the vector. The result is a factor, so
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
fits it exactly as it fits any exposure over the same cells. What the
declaration changes is what gets reported: instead of each cell against
the reference cell, the result is written in the two treatments.

- One row per cell, under the effect label `"mean"`, carrying the
  counterfactual risk (or mean) had everyone been set to that cell. The
  cell is named in the `contrast` column and the group is `"overall"`.

- The simple effects: each treatment's effect within a fixed level of
  the other, with the `contrast` naming the treatment and the level it
  is set to (`"a: 1 vs 0"`) and the `group` naming the level the other
  treatment is held at (`"e = 1"`). These include the comparisons that
  are not against the reference cell, which the vs-reference reporting
  cannot express at all.

- The interaction: the difference between two of the first treatment's
  simple effects, keyed by the two levels of the second being compared
  (`"e = 1 vs e = 0"`). On the risk difference scale this is the
  additive interaction, and on the log risk ratio scale the
  multiplicative one, which is the log of the ratio of the two risk
  ratios.

A 2-by-2 crossing therefore reports fourteen rows for a binary outcome
and nine for a continuous one. The interaction is reported once, under
the first treatment's framing, rather than twice: interaction is
symmetric in the two treatments, so the difference between the first
treatment's two simple effects is the difference between the second
treatment's two, and reporting both would be one quantity under two
labels.

No odds ratio is reported anywhere on this surface. It is
noncollapsible, so neither a simple effect reported beside one nor a
difference of two of them says what it appears to.

Every row is a parameter of the same stacked system, so the covariance
couples them: the interaction covaries with each simple effect it is
built from, and with each cell mean underneath those. The cell means are
the marginal-mean block the categorical path already estimates, and the
simple effects and interaction rows are contrasts over that block, so an
interaction equals the corresponding double difference of the cell means
exactly rather than to within a solver residual.

Two restrictions apply. A declared crossing is reported for the `"ate"`
estimand alone: every cell mean on this surface standardizes to one
population, and a focal estimand standardizes each of them to a
population the simple effects and the interaction are not defined over,
so weights built for anything else error. And `.by` cannot be combined
with a declared crossing: effect modification of a joint intervention is
a three-way question, the interaction between two treatments within the
levels of a third variable, which this surface neither answers nor
sensibly projects. Either restriction is escaped the same way, by
dropping the declaration with `factor(x)` and reporting each cell
against the reference cell.

### Two treatment models instead of a declared crossing

A joint intervention can also be weighted through the sequential
factorization it really has, \\f(A \| L) f(E \| A, L)\\: one treatment
model per treatment, recorded with
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md),
and the product of their weights, built with
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md).
Pass the container as `wt_mod` and an outcome model that reads both
treatments as separate columns.

The reported surface is the one above, row for row and label for label:
the same cell means, simple effects, and interaction, under the same
`effect`, `contrast`, and `group` values. The two routes estimate the
same estimand through different models of the same propensity score, so
they agree as estimators rather than to the last bit, and they coincide
when the two parameterizations are saturated in the same covariates.

Either treatment may have more than two levels, fit with
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) and
weighted with categorical
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
weights. Such a component is stacked under the multinomial score its own
fit solves, which reads no order into the pair, so a categorical
treatment sits in either position. The surface grows with the crossing
rather than changing shape: a 2-by-3 crossing reports twenty-four rows
for a binary outcome and fifteen for a continuous one, being the six
cell means, each treatment's simple effects within each level of the
other, and the first treatment's interaction against each non-reference
level of the second. A categorical component needs no stabilization, and
a caller who asks for it anyway is answered: the marginal proportions of
its levels, a `stabilization_score`, and a fitted
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
numerator each get the block described under **Stabilizing numerators,
one per component**. A categorical treatment may also be paired with a
dose, whose vocabulary names one level contrast per non-reference level;
see the dose section below.

A [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) fit
to two levels is a logistic regression solved by a different optimizer,
and it is stacked as one.
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
records such a component as `"categorical"`, after the model's class,
but the block written for it is the binomial score, the weight factor is
the binary one
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
already builds from such a fit, and the pair reports the same
fourteen-row two-by-two surface a pair of binomial
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) fits reports.

Prefer the two-model route when the two treatments call for different
adjustment sets, or when the dependence of the second treatment on the
first is what you want to model directly, since each treatment then gets
its own model and its own covariates. Prefer the declared crossing when
one model over the cells is the natural specification. The factorization
itself is validated when the container is built rather than here:
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
refuses a second model that does not condition on the first treatment,
and a pair that condition on each other.

Both treatment models are stacked alongside the outcome model, so the
weights entering the outcome score are rebuilt from both propensity
score blocks on every evaluation and the reported standard errors
account for having estimated both. The same two restrictions apply: the
`"ate"` estimand only, and no `.by`. Standard errors come from
M-estimation alone.

### A joint intervention with a dose

The second treatment may be a dose rather than a second discrete
treatment, recorded with any of the classes the single-treatment route
reads a conditional density from, and weighted with a stabilized density
ratio, which
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
requires of a continuous component. A dose has no cells, so there is no
crossing to report counterfactual risks over: the surface is the
marginal structural model's own coefficients, which for an identity-link
model are exactly the weighted fit's coefficients.

The first treatment may be binary or categorical. What the position
settles is the dose: it is the second factor of the factorization, the
treatment whose model conditions on the other, and nothing here carries
the density of a first factor that is one, so a dose in the first slot
and a pair of doses are both refused by type with an error of class
`propensity_ipw_exposure_error`.

Which of two surfaces a fit reports is decided by the marginal
structural model alone. A model written in bare treatment terms and fit
with an intercept reports the vocabulary surface below, whose rows name
the treatment each one varies and where it is evaluated. Every other
treatment-reading model reports the coefficient surface, whose rows are
named after the coefficients they report.

#### The vocabulary surface

For `y ~ a * e` with `a` binary and `e` a dose, three coefficients carry
causal content and each is a row:

- the binary treatment's effect at a dose of zero, with the `contrast`
  naming the treatment and its levels (`"a: 1 vs 0"`) and the `group`
  naming where it is evaluated (`"e = 0"`);

- the dose's slope at the binary treatment's reference level
  (`"e: per unit"`, `"a = 0"`);

- the interaction, the change in the binary treatment's effect per unit
  of the dose, keyed by the one-unit step of the other treatment
  (`"a: 1 vs 0"`, `"e + 1 vs e"`).

The `effect` column names the scale, keeping the vocabulary each
treatment type already uses: a level contrast is a `diff` and a per-unit
change in the dose is a `slope` under an identity link, and both are
`log(or)` under a logit link and `log(rr)` under a log link. An additive
model, `y ~ a + e`, forces one effect for the binary treatment at every
dose and one slope at either of its levels, so it reports two rows under
the group `"overall"` and no interaction.

A categorical first treatment reports those same three readings with the
level contrasts run out over the levels it has: one contrast row per
non-reference level at a dose of zero, one slope row at the reference
level, and one interaction row per non-reference level. For `y ~ z * e`
with `z` taking `"lo"`, `"mid"`, and `"hi"`, that is five rows, in the
outcome design's column order:

- `"z: mid vs lo"` at `"e = 0"` and `"z: hi vs lo"` at `"e = 0"`;

- `"e: per unit"` at `"z = lo"`;

- `"z: mid vs lo"` and `"z: hi vs lo"`, each at `"e + 1 vs e"`.

A binary treatment is the one non-reference level of that rule rather
than a case of its own, so the three rows above are this surface at two
levels. An additive model reports the level contrasts and the slope
under the group `"overall"` and no interaction, K rows for K levels.
Every row is one coefficient of the weighted fit, so the surface grows
by a row per level rather than by a table of contrasts, and per-arm dose
slopes, which are not coefficients of this model, are built downstream
from [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html). No guard here reads how
thin an arm is: a level the dose model rarely sees has its conditional
density read far from any data, and the product weight carries that as a
heavy tail. Check the weights' effective sample size and their spread
within each level, as for any joint weight, before reading the rows.

Those three readings hold of a model that is linear in each treatment
and carries an intercept, and of no other, which is why such a model
reports this surface and no other model does. A model fit without one,
written `- 1` or `+ 0`, reports the coefficient surface: dropping the
intercept expands a factor treatment to an indicator for every level
rather than for the non-reference levels alone, and where the columns
survive it, as a 0/1 numeric treatment's do, the forced zero moves what
the coefficients mean instead, so no row of this vocabulary would be
true of the fit either way.

The columns are read as well as the formula: a factor treatment under a
coding other than treatment contrasts leaves the terms bare while
rescaling or recentering the column, so its coefficients are no longer
the effects these rows name, and such a fit is refused rather than
reported in this vocabulary. Refit the outcome model with the treatment
as an unordered factor under treatment contrasts. A treatment with two
levels also has a numeric coding whose bare term contributes that
column, 0 for the reference level and 1 for the other; a treatment with
more than two levels has none, so there it is the unordered factor
alone. The refusal belongs to this surface rather than to the fit. A
model the vocabulary has no reading for reports the coefficient surface
instead, and a model written without an intercept is one such, so an
ordered factor there contributes an `a.L:e` row named after the column
it multiplies, claiming nothing about a coding.

An interaction of the treatment with the dose written where the dose has
no term of its own, as `y ~ a + a:e`, is refused here too, as is a
crossing carrying no treatment term either, `y ~ a:e`. R codes a factor
in such an interaction with one indicator per level rather than with the
contrasts it carries, so the design gains a reference-level column no
row names and each column after it holds the level before it. Write the
crossing as `y ~ a * e`, which gives the dose a term of its own and
reports on this vocabulary.

The dose column is read the same way, against the dose the treatment
models were fit to. A dose transformed in the outcome model's own
formula, as `I(e / 10)` or `scale(e)`, leaves that column alone and is
no longer a bare term, so the fit reports the coefficient surface. A
dose rescaled or recentered in the data after the treatment models were
fit is refused on this surface instead: the models have to agree on what
the variable is, and the weights the outcome model carries were built
for the dose those models hold. Working on a rescaled dose is supported
either by rescaling it before every model is fit or by writing the
transformation into the outcome formula. Neither refusal reaches the
coefficient surface, which reports whatever columns the model holds and
names each row after its own column: a bare-term model on a rescaled
dose fit without an intercept, `y ~ a * e - 1`, reports the coefficients
and standard errors that reparameterization leaves.

#### The coefficient surface

A marginal structural model carrying a transformed or a basis treatment
term, such as `y ~ a * sin(e)`, `y ~ a + e + I(e^2) + a:e`, or
`y ~ a * splines::ns(e, 3)`, reports one row per treatment-reading
coefficient. The `contrast` column carries the coefficient name exactly
as the fit writes it, and there is no `group` column at all: the surface
makes no claim about where a row is evaluated, because for a curve there
is no one place. Build an interpretable dose response downstream from
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html), whose covariance between
the rows is what that takes. Marginalizing the curve over the observed
doses is a separate estimand that this package does not compute, here as
on the single-treatment route: use the marginaleffects package on the
result, `avg_slopes()` for slopes, `avg_comparisons()` for contrasts,
and `avg_predictions()` for causal dose-response functions.

Unlike the single-treatment route, this one announces nothing and
refuses nothing for such a fit. There the surface a basis exposure
leaves is the only reading the result has, so the reading is declared;
here both readings are computed whatever the marginal structural model
is, and `effects` is answered either way.

The `effect` label follows the outcome link, `coef` at an identity link
and the same `log(or)` and `log(rr)` as everywhere else. `coef` is a
word of its own rather than `diff` or `slope`, both of which carry an
evaluated-at claim on the vocabulary surface.

#### Covariates and mixed terms

A covariate entering on its own is admitted on either surface and
contributes no row, however many columns it expands to: it adjusts the
marginal structural model and the surface has nothing to say about its
coefficient. A term reading a treatment and a covariate together is
refused on both, since its coefficient is a change in an effect per unit
of that covariate and no row could name the effect it stands for.
Nothing here is standardized over the covariates, which is what the
declared crossing does instead.

Unlike the single-treatment route, this one needs no `.data` for a basis
marginal structural model. Without one, each treatment and each
treatment design is read off the model that owns it, so the outcome
model frame is never asked for a column it does not hold.

A `.data` the caller does supply is the frame every design here is
rebuilt from, on the terms the single-treatment routes rebuild theirs:
both treatment columns, both treatment designs, each component's
numerator design, and the counterfactual designs are built over the rows
every model read and under the coding each fit recorded. A column one of
those fits reads that is missing from `.data` errors with class
`propensity_columns_exist_error`, and one supplied as another type
errors with class `propensity_ipw_data_error`, naming the fits that
recorded it, a treatment model by the treatment it fits. A set of fits
made under `na.action = na.exclude` is reconciled with the frame they
were given by restricting `.data` to the rows they analyzed.

The dose model contributes its coefficients and the conditional variance
its density ratio divides by, all as parameters of the same stacked
system. Each component's stabilizing numerator is estimated in a block
of its own; see **Stabilizing numerators, one per component**.

#### Stabilizing numerators, one per component

Either component may be stabilized, and each component's numerator is
estimated in its own block of the stacked system, so the standard errors
account for that numerator having been fitted rather than reading it as
a constant. What the block holds follows what the component's weights
record. A binary component stabilized by the default numerator
contributes the one marginal proportion that numerator is, a categorical
one the \\k - 1\\ free proportions its levels leave, and a dose
stabilized by the default numerator the exposure's two marginal moments.
A component stabilized on a fitted model passed to
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)'s
`stabilize` contributes that model instead: one parameter per
coefficient for a binary component, one run of them per non-reference
level for a categorical one, and for a dose one per coefficient plus one
for the spread its density is read at. A dose stabilized with
`numerator = "integrated"` is built from the dose block and the data
alone and contributes nothing, and an unstabilized discrete component
contributes nothing either.

Those parameters are named `stab_<component>_<parameter>`, where
`<component>` is the name the component was recorded under in
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md):
a joint system carries two components, and a name saying only which role
a parameter plays would not say whose. A categorical component composes
that prefix with the multinomial layout its own block takes,
`<level>:<term>` for a fitted numerator and the level alone for a
marginal one, so a treatment `z` with levels `lo`, `mid`, and `hi`
reports `stab_z_mid` and `stab_z_hi` for the marginal proportions and
`stab_z_mid:(Intercept)`, `stab_z_mid:v`, `stab_z_hi:(Intercept)`, and
`stab_z_hi:v` for a numerator of `v`.

A numerator model is read through the guards its own single-treatment
route reads it through. A binary component's is read as an unpenalized
binomial fit at a logit, probit, or cloglog link; a categorical
component's as an
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) fit
that declares the treatment's levels in the order the component's own
model declares them; and a dose's goes through the registry under
**Continuous propensity score models**. A multinomial numerator fit
under another level order, a model of a response other than the
treatment it stabilizes, or one fit to a different set of observations
errors with class `propensity_ipw_numerator_error`; one fit with case
weights errors with class `propensity_ipw_ps_weights_error`. Every
numerator design here is one of the designs a `.data` the caller
supplies rebuilds, on the terms the single-treatment routes rebuild
theirs, so a numerator that keeps no model frame and whose fitting call
can no longer be evaluated errors with class `propensity_ipw_data_error`
without one, and is rebuilt from `.data` with one.

Both components' numerators arrive in the same argument, one `stabilize`
per component, so every refusal about one names the component it was
built for beside that argument, `stabilize` for `z`. `wt_mod` here is
the container the two treatment models arrived in rather than a model
anyone could refit, so a refusal reading a numerator against its
treatment model names that component too.

What a numerator may condition on is a statement about the reported
model rather than a free choice, and it is the same statement here as
for a single treatment. A numerator may condition on any variable the
marginal structural model already reads, such as a modifier `V`: \\P(Z
\| V)\\ for a discrete treatment and \\f(A \| V)\\ for a dose leave the
estimator consistent for the same effects and only tighten the weights.
A numerator conditioning on a variable the reported model does not read
changes what is estimated rather than the precision, because the
pseudo-population it builds is one in which that variable still predicts
the treatment, and the estimates are then read with `V` or not at all.
The accounting a stacked numerator adds has something to say only where
the reported model is not saturated in the variables the numerator
reads: a numerator that divides out of every cell of the reported model
moves neither the estimates nor their standard errors.

#### The dose model and the ratio its weights are

The dose model goes through the registry described under **Continuous
propensity score models**, so the classes and links this route stacks
are the ones the single-treatment route stacks: an
[`stats::lm()`](https://rdrr.io/r/stats/lm.html), a gaussian
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) fit with an identity
or a log link, a [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html)
fit with one of the psi functions MASS supplies, and an
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) of a gaussian
family read at its own penalized score. A gaussian dose model fit with
another link errors with class `propensity_ipw_link_error`, and a robust
or an additive fit whose score the system cannot write raises the same
conditions it raises for a single dose. The limitations of the additive
route are the ones listed under **Continuous propensity score models**,
and they apply here unchanged.

Every score this route stacks is the unweighted one, so a component
model fit with prior case weights is not at the root of the block
written for it and errors with class `propensity_ipw_ps_weights_error`,
naming the component rather than the container the two models arrived
in. Refit that component without `weights`.

The dose component's weights may be any ratio
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
builds: a heavier-tailed `.density`, `numerator = "integrated"`, or a
single `.sigma` the caller fixed.
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
records each component's ratio on the product, and the stacked system
rebuilds the dose factor from that record, so what is differentiated is
the weight the outcome model was fit with. An integrated numerator is
built from the dose block and the data alone, so it adds no
stabilization parameters; it is the conditional density averaged over
the units the dose model was fit to and read at points spanning their
dose, so it is rebuilt over those rows the way a single dose's is, and a
dose model that keeps no model frame and whose fitting call can no
longer be evaluated errors with class `propensity_ipw_data_error` under
it, `.data` or no `.data`. A fixed scalar `.sigma` is a known constant,
so the dose block is its coefficients alone and none of that number's
uncertainty is carried. A spread supplied for each observation has no
counterpart in the system and errors with class
`propensity_ipw_sigma_error`.

A component built with a `stabilization_score` is rebuilt from the score
itself.
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
records each component's score on the product, and the sandwich treats
it as a known constant here exactly as it does for a single treatment,
so it adds no stabilization parameters and none of its uncertainty is
carried. A product built by a version of this package that recorded no
score, one assembled by hand, or one subset after it was built, which
drops a score held per observation, says a component was stabilized
without saying what by: the system stands that component's own marginal
numerator in, and the weight-consistency check reports the difference,
naming the components whose numerator was stood in for. Rebuild such
weights with
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
to have the score read back. A record left whole over fewer rows is
refused instead: a model frame that drops incomplete rows re-attaches
the record as it was, leaving a component's score at the length the
product was built at, and a score that describes other observations than
the ones being weighted errors with class
`propensity_ipw_stabilization_score_error`.

Dose weights built with a `"kernel"` density are refused with class
`propensity_ipw_se_method_unavailable_error`. The bandwidth is not a
function of the parameters that the sandwich could differentiate, and
the package has no second method to fall back on, so the refusal points
at building the dose weights from a density the stacked system can
write, or at bootstrapping the whole fit by hand.

## Variance estimation

By default (`se_method = "mestimation"`), standard errors are computed
by M-estimation. The propensity score model, the weighted outcome model,
the marginal means, and the effect contrasts are stacked as a single
system of estimating equations, and the empirical sandwich variance of
that joint system is used. The stacked marginal means are standardized
to the estimand's tilted target population, so a non-`ate` estimand with
a covariate-adjusted outcome model reports the contrast for that target
population rather than the full-sample average. Stacking the models
accounts for the uncertainty introduced by estimating the propensity
scores, avoiding the underestimated standard errors that arise from
treating estimated weights as fixed. See Stefanski and Boos (2002) for
the M-estimation framework.

M-estimation standard errors are available for all exposure types:
binary (from a [`stats::glm()`](https://rdrr.io/r/stats/glm.html)
propensity score model), categorical (from a
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) model),
and continuous (from an
[`stats::lm()`](https://rdrr.io/r/stats/lm.html), a gaussian
[`stats::glm()`](https://rdrr.io/r/stats/glm.html), a
[`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html), or an
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) model; see
**Continuous propensity score models** below). The `atm` weight
`pmin(e, 1 - e)` is not differentiable at a propensity score of `0.5`;
deli's central finite difference straddles the kink there and averages
the one-sided slopes. Its effect on the variance is negligible unless
many observations sit at exactly `0.5`.

A binary exposure's propensity score model is read as an unpenalized
logistic fit, whatever the standard error method: the stacked system
solves that score, and both the linearization and the robust diagnostic
differentiate it. An
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) fit of a binomial
family is not at the root of that score but of a penalized one, so it is
refused with class `propensity_ipw_se_method_unavailable_error` rather
than described by a standard error for a model that was not fit. What is
refused is the standard error and not the weights:
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
reads the same fit and builds them from its fitted probabilities. Write
the shape each smooth term carries out as columns of a
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) formula, or bootstrap
the whole fit by hand. An additive dose model is stacked at its
penalized score instead, as described under **Continuous propensity
score models**.

M-estimation standard errors for a categorical exposure allocate memory
roughly linearly in the number of observations, on the order of 70 to 90
kilobytes per observation for a single fit, so a sample of 10,000
observations allocates on the order of hundreds of megabytes. The cost
comes from the stacked estimating-equation machinery rather than any
single term, and for very large samples you can expect long,
garbage-collection-heavy fits. This is expected behavior.

### Stabilizing numerators for a discrete exposure

A numerator estimated by a model the caller passed to
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)'s
`stabilize` is estimated here as well, whichever type the exposure is. A
binary exposure's numerator model contributes one parameter per
coefficient, named `stab_<term>`, and a categorical exposure's
contributes one per coefficient of its
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) fit,
named `stab_<level>:<term>`, so that the level a coefficient belongs to
is read off the name. The multinomial numerator's own score is stacked
alongside the propensity score model's, written as the same multinomial
estimating equation, so the standard errors account for the numerator
having been fitted rather than reading it as a constant. The default
stabilizer contributes what it estimates instead: the single marginal
proportion `stab_pi` for a binary exposure, and the K - 1 free marginal
proportions `stab_<level>` for a categorical one. A
`stabilization_score` contributes nothing at all, being a number the
caller fixed. What a dose's numerator contributes is described under
**Continuous propensity score models**.

Reading an estimated numerator as known would be conservative for these
weights rather than wrong, since stabilization is offered for the `ate`
estimand alone (Lunceford and Davidian 2004; Kostouraki et al. 2024),
and Shu et al. (2020) stack the numerator and the denominator of
stabilized weights for that reason. Stacking the multinomial score
therefore tightens the inference rather than rescuing it. Neither result
isolates the numerator's own contribution, and Kostouraki et al.'s
tutorial is written for a binary exposure, so what they give is the
direction of the correction rather than a bound on its size.

A numerator model's design is one of the designs
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
rebuilds when `.data` is supplied: it is rebuilt over the rows every
model read and under the coding the fit recorded, so a column only the
numerator reads is asked of `.data` and is held to the type the fit
recorded for it. A `.data` missing such a column errors with class
`propensity_columns_exist_error`, and one supplying it as another type
with class `propensity_ipw_data_error`. Without `.data` the design is
read off the fit itself, and a fit that keeps no model frame rebuilds
one by re-evaluating its fitting call. A
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) fit with
`model = FALSE` is such a fit, so one made inside a function whose frame
is gone errors with class `propensity_ipw_data_error`, and the remedy is
to supply `.data`.

A multinomial numerator model is held to the propensity score model's
levels in the propensity score model's own order. Everything its block
contributes is positional: the coefficients are stacked level-major, and
the softmax the numerator is rebuilt from reads the first level as the
reference the others are contrasted with, so a fit made under another
order is a different parameterization rather than a permutation of this
one. Such a fit is refused with an error of class
`propensity_ipw_numerator_error`, which is also what a model of a
response other than the exposure, or one fit to a different number of
observations, errors with. One fit with case weights errors with class
`propensity_ipw_ps_weights_error`, for the reason a weighted propensity
score model does.
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) keeps
no model frame, so the design its coefficients were fit over is
recovered by re-evaluating the fitting call; a fit whose call can no
longer be evaluated, such as one made inside a function whose frame is
gone, errors with class `propensity_ipw_data_error`, and the remedy is
to supply `.data`, which the numerator's design is rebuilt from.

A categorical exposure reports M-estimation standard errors alone, so a
multinomial numerator model meets neither of the other two methods. The
linearization's refusal of a fitted numerator model, below, and the
robust diagnostic's reading of every weight as a known constant are both
matters for a binary exposure.

### Continuous propensity score models

A continuous exposure's weights are a ratio of densities centered on the
conditional mean the propensity score model fits, so the stacked system
carries that model's own score. Four classes supply one: an
[`stats::lm()`](https://rdrr.io/r/stats/lm.html), a gaussian
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) fit with an identity
or a log link, a [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html)
fit with one of the psi functions MASS supplies, and a gaussian
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) fit with one of
the same two links. A model fit with a family that is neither gaussian
nor binomial is neither a dose model nor a binary propensity score
model, and errors with class `propensity_model_family_error` naming the
family it has. A gaussian model fit with any other link is refused by
name, because the coefficients its iteration stops at are not a tight
enough root to seed the stacked solve from; refit it with one of the two
links or as an [`stats::lm()`](https://rdrr.io/r/stats/lm.html). Any
other class errors, since solving it here would report a propensity
score model that was never fit. Each class is read as the class it is
rather than by what it inherits from, so a subclass of a supported fit
errors as well.

A [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) is stacked as
the psi score its coefficients are the root of, clipped where the fit
itself clipped. Each of the three psi functions MASS supplies has a loss
of its own in the stacked system:
[`MASS::psi.huber()`](https://rdrr.io/pkg/MASS/man/rlm.html),
[`MASS::psi.bisquare()`](https://rdrr.io/pkg/MASS/man/rlm.html), and
[`MASS::psi.hampel()`](https://rdrr.io/pkg/MASS/man/rlm.html). `rlm()`
clips the residual divided by its scale estimate at the psi's own
constants, and the stacked score clips the raw residual, so each
constant is multiplied by the scale the fit settled on. Those constants
are read off the psi the fit carries, which is where `rlm()` records a
caller's choice of them, so a retuned psi is stacked at the values it
was retuned to. That scale enters the location score as a known
constant. It solves an equation of its own that the stacked system does
not write, so none of its sampling variability is propagated, and the
standard errors are those of a psi score read at a scale treated as
fixed. The spread of the conditional density is a separate quantity, and
is the pooled residual spread there as it is for every other class.

[`MASS::psi.bisquare()`](https://rdrr.io/pkg/MASS/man/rlm.html) and
[`MASS::psi.hampel()`](https://rdrr.io/pkg/MASS/man/rlm.html) redescend,
and that is worth knowing about the standard errors they get. A
redescending psi's estimating equation has more than one root, so there
is no single set of coefficients it identifies and the solve reports
whichever root it is seeded at. The seed here is the fit's own
coefficients, so what comes back is the root the user has, and the
sandwich is the covariance of that root read locally. It carries no
claim about the other roots, and it does not describe the sampling
behavior of a procedure that might have converged to one of them
instead. Use a bootstrap of the whole fit if that behavior is what you
need.

Weights built with `.sigma = fit$s`, the robust scale, are read as the
fixed spread they are, but that combination is usually a poor choice:
the robust scale resists the very residuals the density is meant to
spread over, so on contaminated data the conditional density is far
narrower than the residuals are, and the ratio can collapse onto an
effective sample size of about one. Leave `.sigma` unset unless you have
a reason to hold the spread at a particular value.

Three robust fits are refused. One fit with a psi function of the
caller's own has no loss here to write its score as, and one fit with
`method = "MM"` reaches its coefficients through a high-breakdown start
that decides which root of a redescending psi the fit finishes at and
supplies the scale it clips at, neither of which is an equation this
system writes. A sandwich read locally at those coefficients would
therefore not describe how an MM fit behaves across samples. Both error
with class `propensity_ipw_robust_psi_error`. One whose iteration
stopped short of its own tolerance is the root of nothing, and errors
with class `propensity_ipw_convergence_error`. The first points to a
refit with one of the three psi functions or to a bootstrap of the whole
fit written by hand, the second to the default `method = "M"`, and the
last to a larger `maxit` or a looser `acc`, which is what a fit needs to
reach its own tolerance.

An [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) is stacked as
the penalized least squares it is, with its smoothing parameters held at
the values the fit reports. Its design is the smooth basis it was built
on rather than the columns its formula names, and its coefficients are
the root of the least-squares score of that basis less the penalty those
smoothing parameters define. An
[`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html) is the same fit
computed for larger data and is read the same way. The additive part of
an [`mgcv::gamm()`](https://rdrr.io/pkg/mgcv/man/gamm.html) fit is not:
it carries the `gam` class alone, without the `glm` and `lm` classes
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reads a model through, and errors with class `propensity_method_error`.
Refit it with [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) or
[`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html). A penalty
attached to a parametric term, such as one from `paraPen`, belongs to no
smooth term and is not one this route can place, so such a fit errors
with class `propensity_ipw_se_method_unavailable_error`. A fit made
under a smoothing floor from `min.sp` errors with the same class,
because mgcv adds that floor to the penalty and leaves it out of the
smoothing parameters the fitted object reports, so nothing in the fit
records the penalty it was made under; refit without `min.sp`, which
[`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html) ignores in any
case.

Four limitations come with the additive route, and each is a property of
the method rather than of a particular fit.

1.  The variance conditions on the smoothing parameters rather than
    estimating them, so none of the uncertainty of having chosen them is
    carried. A simulation put that at roughly three percent of the
    standard error at `n = 500`, roughly seven percent at `n = 200`, and
    less at larger samples. A fit handed its smoothing parameters and a
    fit that chose them report the same standard error here, which is
    the same limitation stated as an equality.

2.  What the propensity score block reports is the frequentist variance
    of a penalized fit rather than the Bayesian covariance mgcv reports
    by default, so [`vcov()`](https://rdrr.io/r/stats/vcov.html) on the
    stacked fit reads like `vcov(fit, freq = TRUE)` and will not match
    `vcov(fit)`.

3.  The envelope is the one every other continuous fit has: a gaussian
    family, an identity or a log link, and no prior weights. A fit
    outside it is refused where every other model of the wrong family,
    link, or weighting is.

4.  Coverage degrades as the weights grow heavier tails, and the
    interval is optimistic once the effective sample size falls much
    below 60 percent of the observations. That holds whatever the dose
    model is: the same simulation found it for a correctly specified
    [`stats::glm()`](https://rdrr.io/r/stats/glm.html) on the same data.

Every density
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
builds a ratio in is stacked as it was built: the normal, Laplace, and
Student t families, and a density you wrote yourself, under either the
marginal or the integrated numerator.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reads the family, the numerator, and the spread off the record the
weights carry and rebuilds the same ratio at each value of the parameter
vector, so the weights it differentiates are the weights the outcome
model was fit with. A marginal numerator contributes the exposure's two
marginal moments to the parameter vector; an integrated one is built
from the propensity score block and the data alone and contributes none.
That numerator is the conditional density averaged over the units the
propensity score model was fit to, read at points spanning their
exposure, so it is rebuilt over those rows even where `.data` leaves the
rest of the system standing on fewer of them. Rebuilding it needs the
design those rows enter through, which `.data` cannot stand in for: a
propensity score model that keeps no model frame, made inside a function
whose frame is gone, errors with class `propensity_ipw_data_error` under
`numerator = "integrated"`, and the remedy is to refit it where its data
is available. Weights carrying no such record, including any written
with [`psw()`](https://r-causal.github.io/propensity/reference/psw.md)
by hand, are read as the normal ratio the package built before the
record existed.

A numerator estimated by a model the caller passed to
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)'s
`stabilize` contributes that model: one parameter per coefficient, named
with a `stab_` prefix so the terms it shares with the propensity score
model stay apart from it, and one for the spread its density is read at,
`sigma2_n`. Its score and the moment its spread is the root of are
stacked alongside the propensity score model's, so the standard errors
account for the numerator having been fitted, which is what separates
this from a `stabilization_score`: a score the caller computed is
carried as a known constant and contributes no parameter at all. The
numerator model is read through the registry above, so a class, family,
or link whose score this system cannot write is refused there in the
terms that registry refuses a propensity score model in, naming
`stabilize`. A model of a response other than the exposure, or one fit
to a different set of observations, errors with class
`propensity_ipw_numerator_error`. One fit with case weights errors with
class `propensity_ipw_ps_weights_error`, for the reason a weighted
propensity score model does: the block written for it is its unweighted
score, and its coefficients are not at the root of that.

A dose numerator model's design is one of the designs
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
rebuilds when `.data` is supplied, on the terms a discrete exposure's
is: it is rebuilt over the rows every model read and under the coding
the fit recorded, so a column only the numerator reads is asked of
`.data` and is held to the type the fit recorded for it. A `.data`
missing such a column errors with class
`propensity_columns_exist_error`, and one supplying it as another type
with class `propensity_ipw_data_error`. Without `.data` the design is
read off the fit itself, and a fit that keeps no model frame rebuilds
one by re-evaluating its fitting call. A
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) fit with
`model = FALSE` is such a fit, so one made inside a function whose frame
is gone errors with class `propensity_ipw_data_error`, and the remedy is
to supply `.data`.

One choice is refused for the standard error rather than for the model.
A `"kernel"` density chooses its bandwidth from the residuals of the
propensity score model, which is not a function of the parameters that
the sandwich could differentiate, so weights built from it have no
closed-form standard error and raise an error of class
`propensity_ipw_se_method_unavailable_error`. The package offers no
resampling method of its own, so the remedy the refusal names is a
bootstrap of the whole pipeline written by hand: resample the rows,
refit the propensity score model, rebuild the weights with
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
and refit the outcome model on each resample, then read the spread of
the coefficient across the resamples. Doing it that way rather than
inside
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
keeps what each replicate re-estimates, such as a kernel bandwidth,
under the writer's control. It is also how to carry a smoothing choice
that the additive route conditions on: a bootstrap that refits the
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) each time
re-chooses the smoothing parameters along with everything else.

### Standard errors by linearization

Setting `se_method = "linearization"` instead uses the
influence-function linearization of Kostouraki et al. (2024). It is
available only for a binary exposure and only for the `ate`, `att`,
`ato`, and `atm` estimands; a categorical or continuous exposure, or the
`atu` or `entropy` estimand, requires `se_method = "mestimation"`. The
linearization also supports only an outcome model whose formula contains
the exposure alone. Its point estimates come from g-computation on the
fitted outcome model, while its influence functions are those of the
Hajek weighted-mean estimator; the two agree only for an exposure-only
outcome model. A covariate-adjusted outcome model (any term beyond the
exposure, including covariates, interactions, or transformed terms)
errors on the linearization path and requires
`se_method = "mestimation"`, which stacks adjusted outcome models
correctly. The conditional reading is restricted the same way:
linearization stacks no system, so the outcome model it stores carries
no block of one, and the standard errors that reading reports are
refused with a classed error rather than replaced by the covariance the
outcome model computed for itself, which treats the estimated weights as
fixed.

Weights stabilized on a fitted numerator model are refused here as well,
with an error of class `propensity_ipw_se_method_unavailable_error`. The
influence functions read the stabilizer as a known constant, which is
exact for the default stabilizer and for a `stabilization_score` the
caller fixed and is not exact for a model fit to the same data, so the
method is refused rather than a standard error reported for a numerator
nobody fit. `se_method = "mestimation"` estimates the numerator model
alongside everything else. `se_method = "robust"` reads every weight as
a known constant by design, that being what the diagnostic is, so it is
not refused.

The linearization outcome model must also carry an intercept, which is a
stricter requirement than the M-estimation path imposes. Under a numeric
coding such as `y ~ z - 1` the linear predictor under no exposure is
fixed at the link's zero point rather than estimated, so the
g-computation marginal means no longer match the Hajek means the
influence functions describe. The rejection is broader than that
mechanism: a saturated factor coding such as `y ~ 0 + zf` does estimate
both cell means, and the M-estimation path accepts it, but the
linearization influence functions are derived for the intercept
parameterization, so every no-intercept outcome model errors here. See
**Model requirements** for the baseline contract both methods impose.

### Standard errors as a diagnostic

Setting `se_method = "robust"` reports the sandwich the weighted outcome
model computes for itself, with the weights entering as known constants.
It is the linearization route with the correction for having estimated
the propensity score dropped, so the point estimates are the
linearization ones and only the standard errors move. What it reports
for the two counterfactual means and the contrasts built from them is
the delta-method reading of
`sandwich::vcovHC(outcome_mod, type = "HC0")`, exactly: the influence
variance is read as a plain mean of squares rather than as a sample
variance, so what comes back is the HC0 sandwich itself rather than a
finite-sample rescaling of it.

This is offered as a diagnostic and not as an inference. It ignores that
the weights were estimated, and the correction it drops is generally
what a propensity score model contributes to the variance, so the
standard errors it reports generally understate it, and understate it
most where the weights are least stable. Use `"mestimation"` or
`"linearization"` to report a standard error. Reach for the diagnostic
to see what the propensity correction is worth on a given fit, by
reading it beside one of them.

A result fit this way carries the class `"ipw_diagnostic_se"` ahead of
`"ipw"`, and printing it writes one line naming the method, so the
numbers are not read as accounting for the propensity score model.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw.md)
carry the same mark as an `ipw_se_diagnostic` attribute on the table.
Every requirement of the linearization route applies here: a binary
exposure, the `ate`, `att`, `ato`, and `atm` estimands, an outcome model
of the exposure alone fit with an intercept, no `.by`, and no
conditional reading. The categorical, continuous, and joint treatment
routes refuse it as they refuse linearization.

The printed mark does not survive
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html):
a pooled object is a summary of results rather than a result, and what
says its standard errors are the diagnostic ones is the method it
records, `se_method = "robust"`.

## Multiple imputation

With missing data, fit the whole analysis once per imputed dataset and
pool the results. Everything the analysis needs is rebuilt inside a
single expression, so each imputation gets its own propensity model, its
own weights, and its own outcome model:

    imp <- mice::mice(dat, m = 20, print = FALSE)
    fits <- with(imp, {
      ps <- glm(z ~ x1 + x2, family = binomial())
      w <- wt_ate(ps)
      om <- glm(y ~ z, family = quasibinomial(), weights = w)
      ipw(ps, om)
    })
    pool_ipw(fits)

[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
is the recommended verb. It reads the results themselves rather than a
tidied table of them, so the effect labels survive, a categorical result
keeps the contrast each row reports, and the complete-data degrees of
freedom fall back to the outcome models whenever a result records none
of its own. That fallback is written for the condition rather than for
one route: a fit under `se_method = "linearization"` or
`se_method = "robust"` records none, and so does a result another
package built on the same class from estimating equations that report no
residual degrees of freedom. It also takes the smallest degrees of
freedom across the pooled results, which is the conservative choice when
they differ.
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw_pooled.md)
and
[`glance()`](https://r-causal.github.io/propensity/reference/glance.ipw_pooled.md)
report what it returns.

[`mice::pool()`](https://amices.org/mice/reference/pool.html) also
works, for every exposure type, and is the right choice when a result
has to travel through the same pipeline as other analyses. Four things
are worth knowing about that route. It reaches the result through
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw.md),
which pools a categorical result grouped by `term` and `contrast`
together, a grouping mice has supported since 3.15.0.
[`mice::pool()`](https://amices.org/mice/reference/pool.html) groups
correctly but its [`summary()`](https://rdrr.io/r/base/summary.html)
prints only `term`, so a categorical pooled table shows each effect
measure repeated with no contrast label beside it, and the labels have
to be read off `pooled$pooled` instead. A result that records no
complete-data degrees of freedom of its own, a fit under
`se_method = "linearization"` or `se_method = "robust"` among them,
leaves the pooled degrees of freedom missing unless `dfcom` is passed
explicitly, as in
`mice::pool(fits, dfcom = df.residual(fits$analyses[[1]]$outcome_mod))`.
Every result of that kind has to be told, where
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
reads the outcome models for all of them unasked. That remedy is what
the package's requirement of mice 3.18.0 or later is for: mice 3.17.0
introduced a regression in the `dfcom` argument of `pool()`, and 3.18.0
is the version that repairs it. And the `exponentiate` argument of
[`summary()`](https://rdrr.io/r/base/summary.html) on a pooled `mipo` is
not the one these methods take: a `mipo` records no scale for the rows
it holds, so it exponentiates every one of them, returning a risk
difference as its exponential and still labeling the row `rd`.
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw_pooled.md)
and the pooled result's own frame exponentiate the rows reported on the
log scale alone and relabel those.

To pool the conditional reading, record it inside the same expression
with
[`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html),
and the pooled result reports the outcome models' coefficients rather
than the causal contrasts.

Recording the reading beforehand is one route to it rather than the only
one.
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
pools both readings of the results it is given, storing the one they
record and keeping the other beside it, so
[`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
and
[`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
move a pooled result between them after the pooling and report a table
the pooling already built rather than pooling a second time.
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw_pooled.md),
[`stats::coef()`](https://rdrr.io/r/stats/coef.html),
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html),
[`stats::confint()`](https://rdrr.io/r/stats/confint.html), and the
pooled result's own frame each take an `effects` argument naming a
reading for the one call, which leaves the pooled result as it is.
[`mice::pool()`](https://amices.org/mice/reference/pool.html) offers
neither: it reaches a result through
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw.md)
and returns a `mipo` holding the one table it asked for, so on that
route the reading still has to be recorded before the pooling.

A reading the pooling could not compute is refused rather than reported
under the other one's name. The conditional reading needs the covariance
the joint estimation of the weights and the outcome implies, and a fit
under `se_method = "linearization"` or `se_method = "robust"` stacks no
such system and records no such block, so a set of those fits pools the
marginal reading alone. The pooled result records why the other reading
is missing, and asking for it, whether by flipping the result or by
naming it for one call, errors with that recorded explanation under the
classes `causalgenerics_pool_missing_surface_conditional` and
`causalgenerics_pool_missing_surface`.

The analysis belongs inside [`with()`](https://rdrr.io/r/base/with.html)
rather than outside it because the propensity score model has to be
estimated once per imputation. Weights built from propensity scores
averaged across imputations are identical in every imputation, which
leaves no between-imputation variance for the weights to contribute and
no per-imputation estimation for the corrected standard errors to
account for; both components of the variance are lost. Fitting within
each imputation keeps the uncertainty of having estimated the weights
inside each result, and the pooling adds the uncertainty the imputation
itself contributed. Leyrat et al. (2019) compare the approaches and
recommend this one.

## Model requirements

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
cannot yet account for propensity scores that were trimmed
([`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)),
truncated
([`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)),
or calibrated
([`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md))
before weighting. The stacked estimating equations rebuild the weights
from the propensity score model on every evaluation, so a weight that is
no longer a deterministic function of that model breaks the sandwich
variance. Supplying weights built from a modified score errors on either
standard error method; refit the weights from the unmodified propensity
score model. An outcome model fit without weights also errors on either
method. The outcome model must not carry an offset term on either
method, since neither the stacked outcome score nor the linearization
influence functions thread an offset; supplying one errors.

For a binary or categorical exposure, the outcome model formula must
contain the exposure. The counterfactual designs are built by setting
the exposure to each level in turn, so a model fit on covariates alone
gives one identical design per level; it errors on either standard error
method.

For those two exposure types, the outcome model must also be able to
represent a baseline at every level, for example through an intercept or
through a saturated factor coding of the exposure. A numeric
no-intercept coding such as `y ~ z - 1` errors on both standard error
methods, for different reasons. On the M-estimation path the
counterfactual design at the zero-coded level is identically zero, so
the marginal mean there is fixed by the outcome link (`0.5` under a
logit or probit link, `0` for a linear model) rather than estimated from
the data. On the linearization and robust paths the intercept is
required outright, because without it the g-computation means stop
matching the Hajek means the influence functions describe. A saturated
factor coding such as `y ~ 0 + zf` is a reparameterization whose designs
are the level indicators, so the M-estimation path accepts it and
reproduces the with-intercept fit; the linearization and robust paths
still require the intercept. A no-intercept model that keeps a
covariate, such as `y ~ z + x1 - 1`, still estimates the marginal mean
from that covariate and runs on the M-estimation path. A continuous
exposure has no counterfactual designs, since
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reports the marginal structural model's own exposure coefficient, and is
unaffected.

Three further requirements apply to the weights and the propensity score
model. The outcome model weights must match the values implied by the
propensity score model; a mismatch errors, on both standard error
methods.

For a continuous exposure that requirement bears on the spread of the
conditional density.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
stacks a single pooled residual variance alongside the propensity score
coefficients, so the weights it rebuilds are the ones
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
produces with its pooled default. A single `.sigma` is taken instead as
a known constant: the weights are rebuilt at the number that was
supplied, and the stacked system carries none of that number's
uncertainty, which is what fixing a spread says. Weights built with an
observation-level `.sigma`, such as `influence(model)$sigma`, are a
different function of the data with no counterpart in the stacked
system, and are refused before anything is solved, with an error of
class `propensity_ipw_sigma_error`. Rebuild the weights with the pooled
default or with one number to use
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).

It bears on a `stabilization_score` in the same way. A score written per
observation is one value per unit, so it does not survive the rows being
restricted: subsetting the weights drops it, and a model frame that
drops incomplete rows leaves it at the length the weights were built at.
Either way the record still names a score as the numerator and the score
in hand no longer describes the observations being weighted, so
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
refuses before anything is solved, with an error of class
`propensity_ipw_stabilization_score_error`. Rebuild the weights on the
rows being analyzed, or stabilize on a single `stabilization_score`,
which scales every weight and survives any restriction.

Weights stabilized on a fitted numerator model carry a further
requirement, on `outcome_mod` rather than on the weights. A numerator
fit on a covariate divides every weight by a quantity that varies with
it, so the pseudo-population the weights build is one in which that
covariate still predicts the exposure. Confounding by it remains there,
and what the weighted fit identifies is the effect within its levels
rather than the marginal one, which only a marginal structural model
carrying it reports (Robins et al. 2000; Cole and Hernán 2008).
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
compares the two: the variables on the right-hand side of the numerator
model's formula against the variables `outcome_mod`'s formula reads
anywhere, and any the second does not read are reported once for the
call, in a warning of class `propensity_ipw_stabilizer_coverage_warning`
naming each of them. The comparison is of variables rather than of
terms, so a model reading `factor(v)`, `poly(v, 2)`, or
`splines::ns(v, 3)` reads `v`; a marginal structural model conditional
on a transformation of a covariate is conditional on the covariate. The
joint route asks the same question of each component's numerator and
reports the answers together. The estimates are the ones the fit gives
either way, so this is a warning and not a refusal: add the variables to
`outcome_mod`, or refit the model given to `stabilize` without them. The
default stabilizer, an intercept-only numerator model, and a
`stabilization_score` all condition on nothing this can name and say
nothing.

The propensity score model must not separate the exposure. A model whose
covariates predict the exposure without error has no finite maximum
likelihood estimate, and the propensity scores its coefficients imply
reach exactly zero or one, leaving the corresponding weights undefined.
For a binary exposure, whose propensity score model is a binomial
[`glm()`](https://rdrr.io/r/stats/glm.html), both standard error methods
reject such a fit at the same threshold: the fitted linear predictors,
put through the link's inverse, saturate for at least one observation.
Nothing short of saturation is rejected. A fitted model's
[`predict()`](https://rdrr.io/r/stats/predict.html) cannot show this on
its own, because the inverse link `glm` uses is bounded away from zero
and one, so a separated fit otherwise yields finite weights and a
standard error that looks ordinary. For a categorical exposure, which
runs on the M-estimation path alone, the threshold is narrower: only a
probability of exactly zero at the level a unit was actually assigned is
rejected, since that alone leaves the unit's own weight undefined, and a
softmax column for a level the unit was not assigned may reach zero
without the fit being refused. A continuous exposure has no saturating
inverse link and is not checked.

The weight functions apply a stricter rule to the propensity scores they
are handed: a categorical matrix holding an exact 0 or 1 anywhere is
refused, whichever level the cell belongs to. A separated
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) fit is
therefore stopped where its weights are built rather than here.
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
repairs such a matrix, its categorical matrix method reading the closed
interval so that a cell at an endpoint is bounded rather than refused,
while
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
holds the open interval as the weight functions do. Weights built from a
truncated score are their own refusal here, since
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
cannot yet account for a modified propensity score, so a separated fit
is repaired for the weight functions rather than for this one. See
**Propensity scores at 0 and 1** in
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for the remedy.

The propensity score model must also be fit without case weights, since
the stacked propensity score equations are unweighted and a weighted fit
would not sit at the score root; that requirement applies to
`se_method = "mestimation"` alone, because neither the linearization
path nor the robust one restacks the propensity score model.
Linearization still corrects for the uncertainty of estimating the
propensity scores; it does so through the influence functions rather
than through a stacked score. The robust diagnostic drops that
correction, which is what makes it a diagnostic.

## References

Stefanski LA, Boos DD. The calculus of M-estimation. *The American
Statistician*. 2002;56(1):29–38.
[doi:10.1198/000313002753631330](https://doi.org/10.1198/000313002753631330)

Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the
inverse probability-of-treatment weighting estimator: A tutorial for
different types of propensity score weights. *Statistics in Medicine*.
2024;43(13):2672–2694.
[doi:10.1002/sim.10078](https://doi.org/10.1002/sim.10078)

Lunceford JK, Davidian M. Stratification and weighting via the
propensity score in estimation of causal treatment effects: a
comparative study. *Statistics in Medicine*. 2004;23(19):2937–2960.
[doi:10.1002/sim.1903](https://doi.org/10.1002/sim.1903)

Shu D, Young JG, Toh S, Wang R. Variance estimation in inverse
probability weighted Cox models. *Biometrics*. 2020.
[doi:10.1111/biom.13332](https://doi.org/10.1111/biom.13332)

Robins JM, Hernán MA, Brumback B. Marginal structural models and causal
inference in epidemiology. *Epidemiology*. 2000;11(5):550–560.
[doi:10.1097/00001648-200009000-00011](https://doi.org/10.1097/00001648-200009000-00011)

Cole SR, Hernán MA. Constructing inverse probability weights for
marginal structural models. *American Journal of Epidemiology*.
2008;168(6):656–664.
[doi:10.1093/aje/kwn164](https://doi.org/10.1093/aje/kwn164)

Hernán MA, Robins JM. *Causal Inference: What If*. Boca Raton: Chapman &
Hall/CRC; 2020.

Leyrat C, Seaman SR, White IR, et al. Propensity score analysis with
partially observed covariates: How should multiple imputation be used?
*Statistical Methods in Medical Research*. 2019;28(1):3–19.
[doi:10.1177/0962280217713032](https://doi.org/10.1177/0962280217713032)

## See also

[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_entropy()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for calculating propensity score weights.

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
for handling extreme propensity scores before weighting.

## Examples

``` r
# Simulate data with a confounder, binary exposure, and binary outcome
set.seed(123)
n <- 200
x1 <- rnorm(n)
z <- rbinom(n, 1, plogis(0.5 * x1))
y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
dat <- data.frame(x1, z, y)

# Step 1: Fit a propensity score model
ps_mod <- glm(z ~ x1, data = dat, family = binomial())

# Step 2: Calculate ATE weights and fit a weighted outcome model
wts <- wt_ate(ps_mod)
#> ℹ Using exposure variable "z" from the propensity score model
#> ℹ Treating `.exposure` as binary
outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

# Step 3: Estimate causal effects with correct standard errors
result <- ipw(ps_mod, outcome_mod)
result
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x1, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
#> 
#> Marginal estimates:
#>                estimate  std.err       z  ci.lower ci.upper conf.level p.value
#> mean 0         0.439827 0.050624  8.6881 0.3406056  0.53905       0.95 < 2e-16
#> mean 1         0.582131 0.048746 11.9421 0.4865905  0.67767       0.95 < 2e-16
#> rd 1 vs 0      0.142304 0.070204  2.0270 0.0047068  0.27990       0.95 0.04266
#> log(rr) 1 vs 0 0.280314 0.142195  1.9713 0.0016172  0.55901       0.95 0.04869
#> log(or) 1 vs 0 0.573392 0.286710  1.9999 0.0114518  1.13533       0.95 0.04551
#>                   
#> mean 0         ***
#> mean 1         ***
#> rd 1 vs 0      *  
#> log(rr) 1 vs 0 *  
#> log(or) 1 vs 0 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Exponentiate log-RR and log-OR to get RR and OR
as.data.frame(result, exponentiate = TRUE)
#>   term contrast  estimate  std.error statistic      p.value
#> 1 mean        0 0.4398270 0.05062412  8.688092 3.685785e-18
#> 2 mean        1 0.5821312 0.04874616 11.942094 7.140277e-33
#> 3   rd   1 vs 0 0.1423042 0.07020402  2.027009 4.266153e-02
#> 4   rr   1 vs 0 1.3235458 0.14219501  1.971337 4.868533e-02
#> 5   or   1 vs 0 1.7742759 0.28670966  1.999906 4.551042e-02

# Continuous outcome example
y_cont <- 2 + 0.8 * z + 0.3 * x1 + rnorm(n)
dat$y_cont <- y_cont
outcome_cont <- lm(y_cont ~ z, data = dat, weights = wts)
ipw(ps_mod, outcome_cont)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x1, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: lm(formula = y_cont ~ z, data = dat, weights = wts) 
#> 
#> Marginal estimates:
#>             estimate  std.err       z ci.lower ci.upper conf.level   p.value
#> mean 0      1.980510 0.105395 18.7913  1.77394   2.1871       0.95 < 2.2e-16
#> mean 1      2.881082 0.090094 31.9788  2.70450   3.0577       0.95 < 2.2e-16
#> diff 1 vs 0 0.900572 0.136610  6.5923  0.63282   1.1683       0.95 4.331e-11
#>                
#> mean 0      ***
#> mean 1      ***
#> diff 1 vs 0 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Effect modification: the effect within each level of a factor, and the
# difference between them. The outcome model carries the exposure-by-modifier
# term that lets the effect differ across the strata.
dat$grp <- factor(rep(c("a", "b"), length.out = n))
ps_grp <- glm(z ~ x1 + grp, data = dat, family = binomial())
wts_grp <- wt_ate(ps_grp)
#> ℹ Using exposure variable "z" from the propensity score model
#> ℹ Treating `.exposure` as binary
outcome_grp <- glm(
  y ~ z * grp + x1,
  data = dat,
  family = quasibinomial(),
  weights = wts_grp
)
ipw(ps_grp, outcome_grp, .by = grp)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x1 + grp, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z * grp + x1, family = quasibinomial(), data = dat, 
#>     weights = wts_grp) 
#> 
#> Marginal estimates:
#>                                    estimate   std.err       z    ci.lower
#> mean 0 overall                     0.431743  0.050756  8.5062  0.33226201
#> mean 1 overall                     0.572814  0.048951 11.7019  0.47687248
#> rd 1 vs 0 overall                  0.141072  0.070930  1.9889  0.00205201
#> log(rr) 1 vs 0 overall             0.282732  0.146153  1.9345 -0.00372240
#> log(or) 1 vs 0 overall             0.568087  0.289473  1.9625  0.00073072
#> mean 0 grp = a                     0.388489  0.075881  5.1197  0.23976443
#> mean 1 grp = a                     0.643072  0.063420 10.1399  0.51877097
#> mean 0 grp = b                     0.474996  0.066425  7.1509  0.34480583
#> mean 1 grp = b                     0.502556  0.074269  6.7667  0.35699105
#> rd 1 vs 0 grp = a                  0.254583  0.099646  2.5549  0.05928135
#> log(rr) 1 vs 0 grp = a             0.503992  0.220170  2.2891  0.07246776
#> rd 1 vs 0 grp = b                  0.027560  0.098203  0.2806 -0.16491347
#> log(rr) 1 vs 0 grp = b             0.056401  0.200510  0.2813 -0.33659204
#> rd 1 vs 0 grp = b vs grp = a      -0.227023  0.138856 -1.6350 -0.49917523
#> log(rr) 1 vs 0 grp = b vs grp = a -0.447591  0.295023 -1.5171 -1.02582580
#>                                   ci.upper conf.level   p.value    
#> mean 0 overall                    0.531223       0.95 < 2.2e-16 ***
#> mean 1 overall                    0.668756       0.95 < 2.2e-16 ***
#> rd 1 vs 0 overall                 0.280091       0.95   0.04671 *  
#> log(rr) 1 vs 0 overall            0.569186       0.95   0.05305 .  
#> log(or) 1 vs 0 overall            1.135444       0.95   0.04971 *  
#> mean 0 grp = a                    0.537214       0.95 3.060e-07 ***
#> mean 1 grp = a                    0.767374       0.95 < 2.2e-16 ***
#> mean 0 grp = b                    0.605186       0.95 8.622e-13 ***
#> mean 1 grp = b                    0.648121       0.95 1.318e-11 ***
#> rd 1 vs 0 grp = a                 0.449885       0.95   0.01062 *  
#> log(rr) 1 vs 0 grp = a            0.935517       0.95   0.02207 *  
#> rd 1 vs 0 grp = b                 0.220033       0.95   0.77898    
#> log(rr) 1 vs 0 grp = b            0.449393       0.95   0.77849    
#> rd 1 vs 0 grp = b vs grp = a      0.045129       0.95   0.10206    
#> log(rr) 1 vs 0 grp = b vs grp = a 0.130643       0.95   0.12923    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# The same fit with the weights stabilized on the modifier. The numerator is
# the probability of the exposure each unit took given `grp`, which is the
# shape the default stabilizer has and conditions only on a variable the
# outcome model reads, so it tightens the weights without moving the estimand.
num_grp <- glm(z ~ grp, data = dat, family = binomial())
p_grp <- fitted(num_grp)
wts_stab <- wt_ate(
  ps_grp,
  stabilize = TRUE,
  stabilization_score = ifelse(dat$z == 1, p_grp, 1 - p_grp)
)
#> ℹ Using exposure variable "z" from the propensity score model
#> ℹ Treating `.exposure` as binary
sd(wts_grp)
#> [1] 0.402489
sd(wts_stab)
#> [1] 0.133368
outcome_stab <- glm(
  y ~ z * grp + x1,
  data = dat,
  family = quasibinomial(),
  weights = wts_stab
)
ipw(ps_grp, outcome_stab, .by = grp)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x1 + grp, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z * grp + x1, family = quasibinomial(), data = dat, 
#>     weights = wts_stab) 
#> 
#> Marginal estimates:
#>                                    estimate   std.err       z    ci.lower
#> mean 0 overall                     0.431772  0.050787  8.5016  0.33223140
#> mean 1 overall                     0.572883  0.048948 11.7038  0.47694619
#> rd 1 vs 0 overall                  0.141111  0.070936  1.9893  0.00207874
#> log(rr) 1 vs 0 overall             0.282784  0.146168  1.9347 -0.00369977
#> log(or) 1 vs 0 overall             0.568249  0.289503  1.9628  0.00083438
#> mean 0 grp = a                     0.388684  0.076110  5.1069  0.23951204
#> mean 1 grp = a                     0.642865  0.063444 10.1329  0.51851776
#> mean 0 grp = b                     0.474860  0.066369  7.1548  0.34477918
#> mean 1 grp = b                     0.502901  0.074211  6.7766  0.35744987
#> rd 1 vs 0 grp = a                  0.254181  0.099876  2.5450  0.05842836
#> log(rr) 1 vs 0 grp = a             0.503168  0.220709  2.2798  0.07058732
#> rd 1 vs 0 grp = b                  0.028041  0.098089  0.2859 -0.16420945
#> log(rr) 1 vs 0 grp = b             0.057374  0.200232  0.2865 -0.33507324
#> rd 1 vs 0 grp = b vs grp = a      -0.226140  0.139012 -1.6268 -0.49859853
#> log(rr) 1 vs 0 grp = b vs grp = a -0.445795  0.295444 -1.5089 -1.02485404
#>                                   ci.upper conf.level   p.value    
#> mean 0 overall                    0.531313       0.95 < 2.2e-16 ***
#> mean 1 overall                    0.668820       0.95 < 2.2e-16 ***
#> rd 1 vs 0 overall                 0.280143       0.95   0.04667 *  
#> log(rr) 1 vs 0 overall            0.569268       0.95   0.05303 .  
#> log(or) 1 vs 0 overall            1.135664       0.95   0.04966 *  
#> mean 0 grp = a                    0.537856       0.95 3.275e-07 ***
#> mean 1 grp = a                    0.767212       0.95 < 2.2e-16 ***
#> mean 0 grp = b                    0.604941       0.95 8.377e-13 ***
#> mean 1 grp = b                    0.648353       0.95 1.230e-11 ***
#> rd 1 vs 0 grp = a                 0.449934       0.95   0.01093 *  
#> log(rr) 1 vs 0 grp = a            0.935749       0.95   0.02262 *  
#> rd 1 vs 0 grp = b                 0.220292       0.95   0.77497    
#> log(rr) 1 vs 0 grp = b            0.449820       0.95   0.77447    
#> rd 1 vs 0 grp = b vs grp = a      0.046319       0.95   0.10379    
#> log(rr) 1 vs 0 grp = b vs grp = a 0.133265       0.95   0.13132    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Continuous exposure: an lm propensity model of the dose on covariates,
# weights, which are stabilized by default for a continuous exposure, and a
# weighted marginal structural outcome model
a <- 0.5 + 0.8 * x1 + rnorm(n)
y_dose <- 1 + 0.6 * a + 0.3 * x1 + rnorm(n)
dat$a <- a
dat$y_dose <- y_dose
ps_cont <- lm(a ~ x1, data = dat)
wts_cont <- wt_ate(ps_cont)
#> ℹ Using exposure variable "a" from the propensity score model
#> ℹ Treating `.exposure` as continuous
msm <- lm(y_dose ~ a, data = dat, weights = wts_cont)
ipw(ps_cont, msm)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: lm(formula = a ~ x1, data = dat) 
#> 
#> Outcome Model:
#>   Call: lm(formula = y_dose ~ a, data = dat, weights = wts_cont) 
#> 
#> Marginal estimates:
#>       estimate  std.err      z ci.lower ci.upper conf.level   p.value    
#> slope 0.371450 0.084413 4.4004    0.206   0.5369       0.95 1.081e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# A dose-response curve written as a basis matrix reports one row per basis
# coefficient. The model frame of such a fit holds the basis rather than the
# dose, so `.data` supplies the exposure.
msm_curve <- lm(y_dose ~ poly(a, 2), data = dat, weights = wts_cont)
ipw(ps_cont, msm_curve, .data = dat)
#> ℹ `ipw()` reports only the conditional reading because the exposure enters `outcome_mod` through more than one term, such as a spline or polynomial.
#> ℹ With a nonlinear dose-response, no single coefficient is the effect of the exposure, so there is no marginal effect to report.
#> ℹ Use the marginaleffects package to marginalize over the dose: `avg_slopes()` for slopes, `avg_comparisons()` for contrasts, and `avg_predictions()` for causal dose-response functions. See <https://marginaleffects.com/chapters/interactions.html>.
#> ℹ Set `effects = "conditional"` to silence this message.
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: conditional (outcome model) 
#> 
#> Weight Estimator:
#>   Call: lm(formula = a ~ x1, data = dat) 
#> 
#> Outcome Model:
#>   Call: lm(formula = y_dose ~ poly(a, 2), data = dat, weights = wts_cont) 
#> 
#> Conditional estimates (outcome model):
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept)  1.21010    0.11507 10.5160 < 2.2e-16 ***
#> poly(a, 2)1  6.35379    1.47843  4.2977 1.726e-05 ***
#> poly(a, 2)2 -2.10894    1.25245 -1.6839   0.09221 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Categorical exposure: a multinomial propensity model and per-level contrasts
set.seed(2)
n <- 300
x1 <- rnorm(n)
a <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
y <- rbinom(n, 1, plogis(-0.5 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.3 * x1))
dat_cat <- data.frame(x1, a, y)

ps_cat <- nnet::multinom(a ~ x1, data = dat_cat, trace = FALSE)
ps_mat <- predict(ps_cat, type = "probs")
wts_cat <- wt_ate(ps_mat, dat_cat$a, exposure_type = "categorical")
outcome_cat <- glm(y ~ a, data = dat_cat, family = quasibinomial(), weights = wts_cat)
ipw(ps_cat, outcome_cat)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: nnet::multinom(formula = a ~ x1, data = dat_cat, trace = FALSE) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ a, family = quasibinomial(), data = dat_cat, 
#>     weights = wts_cat) 
#> 
#> Marginal estimates:
#>                estimate  std.err       z ci.lower ci.upper conf.level   p.value
#> mean a         0.360894 0.047526  7.5936 0.267745  0.45404       0.95 3.111e-14
#> mean b         0.531874 0.048166 11.0425 0.437470  0.62628       0.95 < 2.2e-16
#> mean c         0.570088 0.051072 11.1625 0.469989  0.67019       0.95 < 2.2e-16
#> rd b vs a      0.170980 0.067926  2.5171 0.037847  0.30411       0.95  0.011831
#> log(rr) b vs a 0.387821 0.160396  2.4179 0.073451  0.70219       0.95  0.015610
#> log(or) b vs a 0.699153 0.283716  2.4643 0.143079  1.25523       0.95  0.013729
#> rd c vs a      0.209194 0.070015  2.9878 0.071966  0.34642       0.95  0.002810
#> log(rr) c vs a 0.457206 0.159807  2.8610 0.143990  0.77042       0.95  0.004223
#> log(or) c vs a 0.853695 0.294111  2.9026 0.277247  1.43014       0.95  0.003700
#>                   
#> mean a         ***
#> mean b         ***
#> mean c         ***
#> rd b vs a      *  
#> log(rr) b vs a *  
#> log(or) b vs a *  
#> rd c vs a      ** 
#> log(rr) c vs a ** 
#> log(or) c vs a ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
