# Inverse Probability Weighted Estimation

The `multinom` method estimates causal effects for a categorical
exposure from a fitted
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
propensity score model and a weighted outcome model. Standard errors are
computed by M-estimation; the linearization method is not available for
categorical exposures. For a K-level exposure, effects are reported for
each non-reference level against the reference (first) factor level, and
the estimates table gains a `comparison` column identifying the
contrast.

The `lm` method estimates the causal dose-response effect for a
continuous exposure from a fitted
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) (or gaussian-family
[`stats::glm()`](https://rdrr.io/r/stats/glm.html)) propensity score
model of the exposure and a weighted marginal structural outcome model.
The only supported estimand is `"ate"`. Standard errors are computed by
M-estimation; the linearization method is not available for continuous
exposures. The marginal structural model must contain exactly one
exposure term, and the reported effect is that single coefficient:
`"slope"` for an identity-link outcome, `"log(or)"` for a logit-link
outcome, and `"log(rr)"` for a log-link outcome. The estimates table
keeps the eight-column contract with no comparison column.

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
# S3 method for class 'multinom'
ipw(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization"),
  .focal_level = NULL
)

# S3 method for class 'lm'
ipw(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization")
)

# S3 method for class 'glm'
ipw(
  wt_mod,
  outcome_mod,
  ...,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95,
  se_method = c("mestimation", "linearization")
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
  as `y ~ z - 1` errors on either.

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
estimates table gains a `comparison` column identifying each contrast
(for example `"b vs a"`), so a K-level exposure produces one block of
measures per non-reference level.

For a continuous exposure,
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
reports the single exposure coefficient of the weighted marginal
structural outcome model. Its label follows the outcome link: `slope`
for an identity link, `log(or)` for a logit link, and `log(rr)` for a
log link. The reported coefficient is the one attached to the exposure
term of the marginal structural model. When that term is a
transformation of the exposure (for example an `I(A^2)`-only model), the
coefficient of the transformed term is reported under the same
link-based label.

Use
[`as.data.frame()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.html)
with `exponentiate = TRUE` to obtain risk ratios and odds ratios on
their natural scale.

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
correctly.

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
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x1, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
#> 
#> Estimates:
#>         estimate  std.err      z  ci.lower ci.upper conf.level p.value  
#> rd      0.142304 0.070204 2.0270 0.0047068  0.27990       0.95 0.04266 *
#> log(rr) 0.280314 0.142195 1.9713 0.0016172  0.55901       0.95 0.04869 *
#> log(or) 0.573392 0.286710 1.9999 0.0114518  1.13533       0.95 0.04551 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Exponentiate log-RR and log-OR to get RR and OR
as.data.frame(result, exponentiate = TRUE)
#>   effect  estimate    std.err        z    ci.lower  ci.upper conf.level
#> 1     rd 0.1423042 0.07020402 2.027009 0.004706801 0.2799015       0.95
#> 2     rr 1.3235458 0.14219501 1.971337 1.001618514 1.7489427       0.95
#> 3     or 1.7742759 0.28670966 1.999906 1.011517578 3.1122097       0.95
#>      p.value
#> 1 0.04266153
#> 2 0.04868533
#> 3 0.04551042

# Continuous outcome example
y_cont <- 2 + 0.8 * z + 0.3 * x1 + rnorm(n)
dat$y_cont <- y_cont
outcome_cont <- lm(y_cont ~ z, data = dat, weights = wts)
ipw(ps_mod, outcome_cont)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x1, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: lm(formula = y_cont ~ z, data = dat, weights = wts) 
#> 
#> Estimates:
#>      estimate std.err      z ci.lower ci.upper conf.level   p.value    
#> diff  0.90057 0.13661 6.5923  0.63282   1.1683       0.95 4.331e-11 ***
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
#> 
#> Weight Estimator:
#>   Call: lm(formula = a ~ x1, data = dat) 
#> 
#> Outcome Model:
#>   Call: lm(formula = y_dose ~ a, data = dat, weights = wts_cont) 
#> 
#> Estimates:
#>       estimate  std.err      z ci.lower ci.upper conf.level   p.value    
#> slope 0.371450 0.084413 4.4004    0.206   0.5369       0.95 1.081e-05 ***
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
#> 
#> Weight Estimator:
#>   Call: nnet::multinom(formula = a ~ x1, data = dat_cat, trace = FALSE) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ a, family = quasibinomial(), data = dat_cat, 
#>     weights = wts_cat) 
#> 
#> Estimates:
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
