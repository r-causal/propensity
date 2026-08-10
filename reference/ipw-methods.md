# Inverse Probability Weighted Estimation

The `joint_wt_models` method estimates the effects of a joint
intervention on two treatments from the pair of fitted treatment models
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
records and a weighted outcome model that reads both treatments.
Standard errors are computed by M-estimation; the linearization method
is not available. The only supported estimand is `"ate"`, which is what
the product weights
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
builds target.

The second treatment may be a dose, in which case the surface is the
marginal structural model's own coefficients rather than the cells of a
crossing; see **Joint exposures**.

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
computed by M-estimation; the linearization method is not available for
categorical exposures. For a K-level exposure, effects are reported for
each non-reference level against the reference (first) factor level, and
the estimates table gains a `contrast` column naming the pair of levels
each row compares.

The `lm` method estimates the causal dose-response effect for a
continuous exposure from a fitted
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) (or gaussian-family
[`stats::glm()`](https://rdrr.io/r/stats/glm.html)) propensity score
model of the exposure and a weighted marginal structural outcome model.
The only supported estimand is `"ate"`. Standard errors are computed by
M-estimation; the linearization method is not available for continuous
exposures. Every term of the marginal structural model that reads the
exposure must read the exposure and nothing else, and the reported
effects are the coefficients those terms contribute, labeled `"log(or)"`
for a logit-link outcome and `"log(rr)"` for a log-link one. A model
with one exposure coefficient keeps the eight-column contract with no
contrast column and labels its row `"slope"` at an identity link; a
dose-response curve such as `y ~ A + I(A^2)` reports one row per
coefficient under `"coef"`, gaining a contrast column that names them.

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
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or identity-link
gaussian [`stats::glm()`](https://rdrr.io/r/stats/glm.html) propensity
score model with a weighted marginal structural outcome model whose link
is identity, logit, or log, and
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reports the single exposure coefficient of that model. A subclass of
either propensity score model, such as a robust or an additive fit,
errors.

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
  se_method = c("mestimation", "linearization"),
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
  se_method = c("mestimation", "linearization"),
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
  se_method = c("mestimation", "linearization"),
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
  se_method = c("mestimation", "linearization"),
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
  `se_method = "linearization"` errors, and `se_method = "mestimation"`
  reports the resulting weights as inconsistent with the propensity
  score model.

- conf_level:

  Confidence level for intervals. Default is `0.95`.

- se_method:

  Method for standard error estimation. `"mestimation"` (the default)
  stacks the propensity score, outcome, and estimand estimating
  equations and returns the empirical sandwich variance.
  `"linearization"` uses the influence-function linearization of
  Kostouraki et al. (2024). Both account for the uncertainty of
  estimating the propensity scores. `"linearization"` supports only an
  outcome model of the exposure alone, fit with an intercept; a
  covariate-adjusted outcome model requires `"mestimation"`. For a
  binary or categorical exposure, both methods require an outcome model
  that can represent a baseline, so a numeric no-intercept coding such
  as `y ~ z - 1` errors on either. The standard errors of the
  conditional reading described under `effects` require `"mestimation"`
  as well: the linearization route stores its outcome model unwrapped,
  so [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html),
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

  The covariance the conditional reading reports is the outcome block of
  the jointly estimated sandwich, which every route that stacks
  estimating equations attaches to the outcome model it stores:
  `se_method = "mestimation"` for a binary exposure, and the categorical
  and continuous routes, which run on M-estimation alone. A
  linearization fit stacks no such system and records no such block, so
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

  Either `"mestimation"` or `"linearization"`.

- `fit`:

  The fitted M-estimator object when `se_method` is `"mestimation"`,
  otherwise `NULL`. Calling
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) or
  [`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
  on this object exposes the full stacked system of estimating
  equations, including the propensity score, outcome, and estimand
  parameters.

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
measure, and for the effect measure and the contrast together where a
categorical exposure reports one row per contrast. Which surface
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
non-reference level against the reference (first) factor level. The
estimates table gains a `contrast` column naming the pair of levels each
row compares (for example `"b vs a"`), so a K-level exposure produces
one block of measures per non-reference level.

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
in three blocks, in this order:

- the overall rows, unchanged from the fit without `.by`, under the
  group label `"overall"`;

- one block per level of the modifier, under a `"var = value"` label
  such as `"sex = female"`;

- one block per non-reference level against the reference level, which
  is the modifier's first level, under a `"var = value vs var = value"`
  label.

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
(K - 1) contrasts. A row is named by three things there, the measure,
the contrast, and the subgroup, and the `estimates` frame carries
`contrast` and `group` in that order after `effect`. Within a block the
ordering is the one an ungrouped categorical fit already uses,
contrast-major and measure-minor, so a contrast's measures sit together;
the blocks themselves run group-major.

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
which that variable still predicts the exposure. Condition the numerator
on the modifier and on nothing the model does not hold. The same rule
governs an exposure-by-modifier interaction reported without `.by`,
since what matters is what the reported model reads and not which
argument asked for the rows.

Build the numerator with
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)'s
`stabilization_score`, evaluated at the exposure each unit actually
took; see **Stabilization** in
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for the recipe. The sandwich treats a supplied score as a known
constant, so the standard errors do not account for the numerator model
having been fitted. The default stabilizer is estimated in the stacked
system instead, which is one parameter wider as a result.

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
treatment, recorded with an
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or an identity-link
gaussian [`stats::glm()`](https://rdrr.io/r/stats/glm.html) and weighted
with a stabilized density ratio, which
[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
requires of a continuous component. A dose has no cells, so there is no
crossing to report counterfactual risks over: the surface is the
marginal structural model's own coefficients, which for an identity-link
model are exactly the weighted fit's coefficients.

Which of two surfaces a fit reports is decided by the marginal
structural model alone. A model written in bare treatment terms reports
the vocabulary surface below, whose rows name the treatment each one
varies and where it is evaluated. Every other treatment-reading model
reports the coefficient surface, whose rows are named after the
coefficients they report.

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

Those three readings hold of a model that is linear in each treatment
and of no other, which is why such a model reports this surface and no
other model does. The columns are read as well as the formula: a factor
treatment under a coding other than treatment contrasts leaves the terms
bare while rescaling or recentering the column, so its coefficients are
no longer the effects these rows name, and such a fit is refused rather
than reported on either surface. Refit the outcome model with the
treatment as a 0/1 numeric or as an unordered factor under treatment
contrasts.

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
the rows is what that takes.

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
marginal structural model. Each treatment is read off its own treatment
model, so the outcome model frame is never asked for a column it does
not hold.

The dose model contributes its coefficients and the conditional variance
its density ratio divides by, and the stabilizing numerator contributes
the dose's two marginal moments, all as parameters of the same stacked
system.

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
and continuous (from an [`stats::lm()`](https://rdrr.io/r/stats/lm.html)
or identity-link gaussian
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) model). The `atm`
weight `pmin(e, 1 - e)` is not differentiable at a propensity score of
`0.5`; deli's central finite difference straddles the kink there and
averages the one-sided slopes. Its effect on the variance is negligible
unless many observations sit at exactly `0.5`.

M-estimation standard errors for a categorical exposure allocate memory
roughly linearly in the number of observations, on the order of 70 to 90
kilobytes per observation for a single fit, so a sample of 10,000
observations allocates on the order of hundreds of megabytes. The cost
comes from the stacked estimating-equation machinery rather than any
single term, and for very large samples you can expect long,
garbage-collection-heavy fits. This is expected behavior.

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
one route: a fit under `se_method = "linearization"` records none, and
so does a result another package built on the same class from estimating
equations that report no residual degrees of freedom. It also takes the
smallest degrees of freedom across the pooled results, which is the
conservative choice when they differ.
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
`se_method = "linearization"` among them, leaves the pooled degrees of
freedom missing unless `dfcom` is passed explicitly, as in
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
under `se_method = "linearization"` stacks no such system and records no
such block, so a set of those fits pools the marginal reading alone. The
pooled result records why the other reading is missing, and asking for
it, whether by flipping the result or by naming it for one call, errors
with that recorded explanation under the classes
`causalgenerics_pool_missing_surface_conditional` and
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
the data. On the linearization path the intercept is required outright,
because without it the g-computation means stop matching the Hajek means
the influence functions describe. A saturated factor coding such as
`y ~ 0 + zf` is a reparameterization whose designs are the level
indicators, so the M-estimation path accepts it and reproduces the
with-intercept fit; the linearization path still requires the intercept.
A no-intercept model that keeps a covariate, such as `y ~ z + x1 - 1`,
still estimates the marginal mean from that covariate and runs on the
M-estimation path. A continuous exposure has no counterfactual designs,
since
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reports the marginal structural model's own exposure coefficient, and is
unaffected.

Three further requirements apply to the weights and the propensity score
model. The outcome model weights must match the values implied by the
propensity score model; a mismatch errors, on both standard error
methods.

For a continuous exposure that requirement fixes the spread of the
conditional density.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
stacks a single pooled residual variance alongside the propensity score
coefficients, so the weights it rebuilds are the ones
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
produces with its pooled default. Weights built with an
observation-level `.sigma`, such as `influence(model)$sigma`, are a
different function of the data with no counterpart in the stacked
system, and the consistency check refuses them. Refit the weights
without `.sigma` to use
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).

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
therefore stopped where its weights are built rather than here, and
neither
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
nor
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
offers a way past that, since both validate a categorical matrix under
the same open interval. See **Propensity scores at 0 and 1** in
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for the remedy.

The propensity score model must also be fit without case weights, since
the stacked propensity score equations are unweighted and a weighted fit
would not sit at the score root; that requirement applies to
`se_method = "mestimation"` alone, because the linearization path does
not restack the propensity score model. It still corrects for the
uncertainty of estimating the propensity scores; it does so through the
influence functions rather than through a stacked score.

## References

Stefanski LA, Boos DD. The calculus of M-estimation. *The American
Statistician*. 2002;56(1):29–38.
[doi:10.1198/000313002753631330](https://doi.org/10.1198/000313002753631330)

Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the
inverse probability-of-treatment weighting estimator: A tutorial for
different types of propensity score weights. *Statistics in Medicine*.
2024;43(13):2672–2694.
[doi:10.1002/sim.10078](https://doi.org/10.1002/sim.10078)

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
#> ℹ Using exposure variable "z" from GLM model
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
#>         estimate  std.err      z  ci.lower ci.upper conf.level p.value  
#> rd      0.142304 0.070204 2.0270 0.0047068  0.27990       0.95 0.04266 *
#> log(rr) 0.280314 0.142195 1.9713 0.0016172  0.55901       0.95 0.04869 *
#> log(or) 0.573392 0.286710 1.9999 0.0114518  1.13533       0.95 0.04551 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Exponentiate log-RR and log-OR to get RR and OR
as.data.frame(result, exponentiate = TRUE)
#>   term  estimate  std.error statistic    p.value
#> 1   rd 0.1423042 0.07020402  2.027009 0.04266153
#> 2   rr 1.3235458 0.14219501  1.971337 0.04868533
#> 3   or 1.7742759 0.28670966  1.999906 0.04551042

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
#>      estimate std.err      z ci.lower ci.upper conf.level   p.value    
#> diff  0.90057 0.13661 6.5923  0.63282   1.1683       0.95 4.331e-11 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Effect modification: the effect within each level of a factor, and the
# difference between them. The outcome model carries the exposure-by-modifier
# term that lets the effect differ across the strata.
dat$grp <- factor(rep(c("a", "b"), length.out = n))
ps_grp <- glm(z ~ x1 + grp, data = dat, family = binomial())
wts_grp <- wt_ate(ps_grp)
#> ℹ Using exposure variable "z" from GLM model
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
#>                             estimate   std.err       z    ci.lower ci.upper
#> rd overall                  0.141072  0.070930  1.9889  0.00205201 0.280091
#> log(rr) overall             0.282732  0.146153  1.9345 -0.00372240 0.569186
#> log(or) overall             0.568087  0.289473  1.9625  0.00073072 1.135444
#> rd grp = a                  0.254583  0.099646  2.5549  0.05928135 0.449885
#> log(rr) grp = a             0.503992  0.220170  2.2891  0.07246776 0.935517
#> rd grp = b                  0.027560  0.098203  0.2806 -0.16491347 0.220033
#> log(rr) grp = b             0.056401  0.200510  0.2813 -0.33659204 0.449393
#> rd grp = b vs grp = a      -0.227023  0.138856 -1.6350 -0.49917523 0.045129
#> log(rr) grp = b vs grp = a -0.447591  0.295023 -1.5171 -1.02582580 0.130643
#>                            conf.level p.value  
#> rd overall                       0.95 0.04671 *
#> log(rr) overall                  0.95 0.05305 .
#> log(or) overall                  0.95 0.04971 *
#> rd grp = a                       0.95 0.01062 *
#> log(rr) grp = a                  0.95 0.02207 *
#> rd grp = b                       0.95 0.77898  
#> log(rr) grp = b                  0.95 0.77849  
#> rd grp = b vs grp = a            0.95 0.10206  
#> log(rr) grp = b vs grp = a       0.95 0.12923  
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
#> ℹ Using exposure variable "z" from GLM model
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
#>                             estimate   std.err       z    ci.lower ci.upper
#> rd overall                  0.141111  0.070936  1.9893  0.00207874 0.280143
#> log(rr) overall             0.282784  0.146168  1.9347 -0.00369977 0.569268
#> log(or) overall             0.568249  0.289503  1.9628  0.00083438 1.135664
#> rd grp = a                  0.254181  0.099876  2.5450  0.05842836 0.449934
#> log(rr) grp = a             0.503168  0.220709  2.2798  0.07058732 0.935749
#> rd grp = b                  0.028041  0.098089  0.2859 -0.16420945 0.220292
#> log(rr) grp = b             0.057374  0.200232  0.2865 -0.33507324 0.449820
#> rd grp = b vs grp = a      -0.226140  0.139012 -1.6268 -0.49859853 0.046319
#> log(rr) grp = b vs grp = a -0.445795  0.295444 -1.5089 -1.02485404 0.133265
#>                            conf.level p.value  
#> rd overall                       0.95 0.04667 *
#> log(rr) overall                  0.95 0.05303 .
#> log(or) overall                  0.95 0.04966 *
#> rd grp = a                       0.95 0.01093 *
#> log(rr) grp = a                  0.95 0.02262 *
#> rd grp = b                       0.95 0.77497  
#> log(rr) grp = b                  0.95 0.77447  
#> rd grp = b vs grp = a            0.95 0.10379  
#> log(rr) grp = b vs grp = a       0.95 0.13132  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Continuous exposure: an lm propensity model of the dose on covariates,
# stabilized weights, and a weighted marginal structural outcome model
a <- 0.5 + 0.8 * x1 + rnorm(n)
y_dose <- 1 + 0.6 * a + 0.3 * x1 + rnorm(n)
dat$a <- a
dat$y_dose <- y_dose
ps_cont <- lm(a ~ x1, data = dat)
wts_cont <- wt_ate(
  fitted(ps_cont),
  a,
  exposure_type = "continuous",
  stabilize = TRUE
)
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
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: lm(formula = a ~ x1, data = dat) 
#> 
#> Outcome Model:
#>   Call: lm(formula = y_dose ~ poly(a, 2), data = dat, weights = wts_cont) 
#> 
#> Marginal estimates:
#>                  estimate std.err       z ci.lower ci.upper conf.level
#> coef poly(a, 2)1   6.3538  1.4784  4.2977   3.4561  9.25145       0.95
#> coef poly(a, 2)2  -2.1089  1.2525 -1.6838  -4.5637  0.34582       0.95
#>                    p.value    
#> coef poly(a, 2)1 1.726e-05 ***
#> coef poly(a, 2)2   0.09221 .  
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
#>                estimate  std.err      z ci.lower ci.upper conf.level  p.value
#> rd b vs a      0.170980 0.067926 2.5171 0.037847  0.30411       0.95 0.011831
#> log(rr) b vs a 0.387821 0.160396 2.4179 0.073451  0.70219       0.95 0.015610
#> log(or) b vs a 0.699153 0.283716 2.4643 0.143079  1.25523       0.95 0.013729
#> rd c vs a      0.209194 0.070015 2.9878 0.071966  0.34642       0.95 0.002810
#> log(rr) c vs a 0.457206 0.159807 2.8610 0.143990  0.77042       0.95 0.004223
#> log(or) c vs a 0.853695 0.294111 2.9026 0.277247  1.43014       0.95 0.003700
#>                  
#> rd b vs a      * 
#> log(rr) b vs a * 
#> log(or) b vs a * 
#> rd c vs a      **
#> log(rr) c vs a **
#> log(or) c vs a **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
