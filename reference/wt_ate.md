# Calculate propensity score weights

Compute inverse probability weights for causal inference under different
estimands. Each function targets a different population:

- `wt_ate()`: **Average Treatment Effect** – the full population.

- `wt_att()`: **Average Treatment Effect on the Treated** – the treated
  (focal) group.

- `wt_atu()`: **Average Treatment Effect on the Untreated** – the
  untreated (reference) group. `wt_atc()` is an alias.

- `wt_atm()`: **Average Treatment Effect for the Evenly Matchable** –
  units with the most overlap.

- `wt_ato()`: **Average Treatment Effect for the Overlap Population** –
  weights proportional to overlap.

- `wt_entropy()`: **Entropy-weighted Average Treatment Effect** – the
  entropy-tilted population.

- `wt_cens()`: **Inverse probability of censoring weights** – uses the
  same formula as `wt_ate()` but labels the estimand `"uncensored"`. Use
  these to adjust for censoring in survival analysis, not for treatment
  weighting. Remaining uncensored is a two-level event, so `wt_cens()`
  supports binary and continuous exposures only.

`.propensity` accepts a numeric vector of predicted probabilities, a
`data.frame` or matrix of per-level probabilities, a fitted propensity
score model, or a modified propensity score created by
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md),
or
[`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md).

All functions return a
[`psw`](https://r-causal.github.io/propensity/reference/psw.md) object –
a numeric vector that tracks the estimand, stabilization status, and any
trimming or truncation applied.

## Usage

``` r
wt_ate(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = NULL,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_ate(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = NULL,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_att(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_att(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_atu(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_atu(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_atm(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_atm(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_ato(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_ato(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_entropy(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_entropy(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_atc(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

wt_cens(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = NULL,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_cens(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = NULL,
  stabilization_score = NULL,
  ...,
  .density = "normal",
  numerator = c("marginal", "integrated"),
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)
```

## Arguments

- .propensity:

  Propensity scores in one of several forms:

  - A **numeric vector** of predicted probabilities (binary/continuous).

  - A **data frame** or **matrix** with one column per exposure level.
    Both shapes are read the same way, for categorical exposures and for
    binary ones alike, including the choice of column described under
    `.propensity_col`. That argument is itself a formal of the data
    frame methods only, so selecting a column by name means passing a
    data frame. A frame of probabilities predicted from a fitted
    tidymodels classification model has its columns named
    `.pred_<level>`, and the level a column holds is read from the part
    of the name after that prefix. A frame holding a `.pred_class`
    column, which the same models return when no prediction type is
    named, carries predicted levels rather than probabilities and is
    refused with an error of class `propensity_df_class_column_error`.

  - A fitted **propensity score model** – fitted values are extracted
    automatically. For a binary exposure that is a
    [`glm()`](https://rdrr.io/r/stats/glm.html) fit with
    [`binomial()`](https://rdrr.io/r/stats/family.html) or
    [`quasibinomial()`](https://rdrr.io/r/stats/family.html). For a
    continuous exposure it is a model of the conditional mean:
    [`lm()`](https://rdrr.io/r/stats/lm.html), a
    [`glm()`](https://rdrr.io/r/stats/glm.html) fit with
    [`gaussian()`](https://rdrr.io/r/stats/family.html) under any of its
    links, an [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) fit
    with [`gaussian()`](https://rdrr.io/r/stats/family.html), or a
    [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html). For a
    categorical exposure it is an
    [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html),
    whose fitted values hold the probability of every level; a
    multinomial fit of two levels reports a single probability and is
    read as a model of a binary exposure. Handing one exposure type a
    model of another, or a family whose spread changes with its fitted
    values, is an error of class `propensity_model_family_error`. The
    one exception is the multinomial fit, which offers no continuous
    type at all, so naming one for it is an error of class
    `causalgenerics_unsupported_exposure_type`. See **Exposure types**
    and **Continuous exposures** in Details.

  - A modified propensity score created by
    [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
    [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
    [`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md),
    or
    [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md).

  For a binary exposure, `.propensity` is the probability of the focal
  level, so a numeric vector must be supplied on the scale of whichever
  level `.focal_level` or `.reference_level` resolves to. The model
  methods derive it from the fitted values, which give the probability
  of the response's second level, and subtract them from one when the
  resolved focal level is the first level instead. A data frame or
  matrix reduces to a single column, read on the same scale; see
  `.propensity_col` for how that column is chosen. A matrix reduces by
  those same rules, but `.propensity_col` belongs to the data frame
  methods alone, so a matrix whose column you want to name has to be
  converted with
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) first.

  For a continuous exposure, `.propensity` is the conditional mean of
  the exposure. A data frame of a single column is read as that mean; a
  data frame of several columns is refused with an error of class
  `propensity_df_ambiguous_column_error` unless `.propensity_col` names
  the column to read, since a continuous exposure has no levels the
  columns could be matched against.

- .exposure:

  The exposure (treatment) variable. For binary exposures, a numeric 0/1
  vector, logical, or two-level factor. For categorical exposures, a
  factor or character vector. For continuous exposures, a numeric
  vector. Optional when `.propensity` is a fitted model, in which case
  it is the model's response. Missing values are not counted as a level
  of their own, and are carried through to the weights as missing.

- .sigma:

  The residual spread of the conditional density for continuous
  exposures: a single standard deviation applied to every unit, or one
  for each observation (e.g., `influence(model)$sigma`). Optional: with
  none supplied, including when `.propensity` is a fitted model, the
  conditional density uses the pooled residual spread of `.exposure`
  around `.propensity`.

  The two shapes differ downstream. A single spread is a constant the
  weights can be rebuilt from, so it is recorded in
  [`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md)
  as `sigma_value` and
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  reads it as a spread held fixed, propagating none of its uncertainty.
  A spread supplied for each observation is a function of the data that
  nothing rebuilds, so the record holds where it came from and nothing
  more, and
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  refuses such weights; see **Continuous exposures** in Details.

  Must be numeric, and applies only to continuous exposures. `.sigma`
  sits in the third position, which is where a value meant for
  `exposure_type` arrives when it is supplied without a name, so
  anything else is refused with an error of class
  `propensity_sigma_error`. A spread you supply and a density built with
  `sigma_method = "mle"`, which estimates one from the residuals, are
  two instructions about the same quantity and are refused together with
  an error of class `propensity_density_error`.

- exposure_type:

  Type of exposure: `"auto"` (default), `"binary"`, `"categorical"`, or
  `"continuous"`. `"auto"` detects the type from `.exposure`. Not every
  weight function answers every type; see **Exposure types** in Details
  for which does which. Naming a type a function does not support, or
  supplying an exposure detection reads as one, is an error of class
  `causalgenerics_unsupported_exposure_type`.

- .focal_level:

  The value of `.exposure` representing the focal (treated) group. Every
  binary coding honors it: 0/1 numeric, logical, two-level factor, and
  two-level character exposures are all coded with the named level as
  focal, and a level the exposure never takes is an error. With no level
  named, a binary exposure defaults to its higher level, which is `1`
  for a 0/1 exposure, `TRUE` for a logical one, and the second of the
  two levels a factor or character exposure takes. Levels a factor
  declares but never takes are not candidates. Naming any other level
  reverses the coding, so `.propensity` must then hold the probability
  of the named level. Required for `wt_att()` and `wt_atu()` with
  categorical exposures.

- .reference_level:

  The value of `.exposure` representing the reference (control) group.
  For a binary exposure, naming it makes the exposure's other level
  focal, with the same consequence for `.propensity`, and a level the
  exposure never takes is an error. Automatically detected if not
  supplied.

- stabilize:

  Whether to multiply the weights by an estimate of the marginal
  probability of the exposure a unit took (binary or categorical) or of
  its density (continuous), and what that estimate is. It takes one of
  three forms:

  - A logical. `TRUE` and `FALSE` ask for stabilization or its absence
    outright, and an unstabilized continuous exposure reports that its
    weights are not the recommended ones. `NULL`, the default, reads the
    answer from the exposure type: a continuous exposure is stabilized,
    and a binary or categorical exposure is not.

  - A fitted model of the exposure, which stabilizes the weights on what
    that model estimates rather than on the marginal probability or
    density of the exposure.

    For a binary exposure the model is a
    [`binomial()`](https://rdrr.io/r/stats/family.html) fit and the
    numerator is its fitted probability of the level each unit took, so
    the weights are \\P(A = a_i \mid V_i) / f(a_i \mid X_i)\\ for the
    variables \\V\\ the numerator model reads. For a continuous exposure
    the numerator is the family `.density` names, read at the model's
    fitted mean and the root mean square of its residuals, so the
    weights are \\f(A \mid V) / f(A \mid X)\\. For a categorical
    exposure the model is an
    [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) fit
    and the numerator is its fitted probability of the level each unit
    took, read from the column named for that level, so the weights are
    \\P(Z = z_i \mid V_i) / f(z_i \mid X_i)\\. The multinomial fit has
    to report a probability for every level the exposure has, matched by
    name rather than by position.

    The model is recorded on the result, where
    [`numerator_model()`](https://r-causal.github.io/propensity/reference/numerator_model.md)
    reads it back, and
    [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
    estimates it alongside everything else so that the standard errors
    account for it having been fitted, whichever type the exposure is: a
    categorical exposure's multinomial score is stacked there the way a
    binary exposure's binomial one is. Conditioning the numerator on
    \\V\\ changes what is estimated unless the model the estimates are
    read from also reads \\V\\; see **Stabilization** in Details.

    Any [`lm()`](https://rdrr.io/r/stats/lm.html), or anything built on
    one, is read this way, as is an
    [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html).
    The model is held to the exposure: a fit whose spread changes with
    its fitted values, or a model of a conditional mean where a
    probability is needed, is refused with
    `propensity_model_family_error`, which is also what a model of a
    dose handed to a binary exposure, a
    [`binomial()`](https://rdrr.io/r/stats/family.html) model handed to
    a dose, a multinomial fit handed either of them, and a fit reporting
    some other set of levels than the exposure's are refused with. A
    model supplied together with `stabilization_score` or
    `numerator = "integrated"` is `propensity_numerator_error`. So is a
    model fit with case `weights`, which estimates the numerator in a
    reweighted sample rather than in the one the weights are being built
    for. A model with a fitted value for some other set of observations
    is `propensity_length_error`. That check counts the fitted values
    rather than reading which rows they came from, so a model fit to
    another dataset of the same length passes it: fit the numerator
    model on the data the weights are being built for.

  Anything else is refused with an error of class
  `propensity_stabilize_error`. Stabilization is only supported by
  `wt_ate()` and `wt_cens()`, since the other estimands' weights already
  carry a tilting function where the numerator would go and stabilizing
  them would move the population they target rather than their variance.
  See **Stabilization** in Details.

- stabilization_score:

  Optional stabilization multiplier to use instead of the default
  described under **Stabilization**: the marginal mean of `.exposure`,
  or its marginal normal density for a continuous exposure. Either a
  single value applied to every weight or a numeric vector holding one
  value per observation, which is multiplied into the weights
  observation by observation. Every value must be positive and finite,
  and any other length or value is refused with an error of class
  `propensity_stabilization_score_error`. A per-observation score is
  recorded on the result and is dropped, with a warning, by any
  operation that changes the length of the weight vector, since it
  cannot be re-indexed for the new length. Ignored when stabilization is
  off, whether because `stabilize` is `FALSE` or because the exposure
  type resolves the default that way.

- ...:

  These dots are for future extensions and must be empty.

- .density:

  The family of the conditional density for continuous exposures,
  described under **Continuous exposures** in Details. One of the
  strings `"normal"` (the default), `"laplace"`, or `"kernel"`; a
  specification built by
  [`dens_normal()`](https://r-causal.github.io/propensity/reference/dens_normal.md),
  [`dens_laplace()`](https://r-causal.github.io/propensity/reference/dens_normal.md),
  [`dens_t()`](https://r-causal.github.io/propensity/reference/dens_normal.md),
  [`dens_kernel()`](https://r-causal.github.io/propensity/reference/dens_normal.md),
  or
  [`dens_fn()`](https://r-causal.github.io/propensity/reference/dens_normal.md);
  or a bare function of the standardized residual, which is wrapped with
  [`dens_fn()`](https://r-causal.github.io/propensity/reference/dens_normal.md).

  A function is called on the whole vector of standardized residuals at
  once rather than on one residual at a time, and is called a second
  time on the standardized exposure when the marginal numerator
  stabilizes the weights. Each call must return one finite, non-negative
  value for each value it was given, not every one of them zero.
  Anything else is refused with an error of class
  `propensity_density_error`, which is also what a family other than the
  default gets for a binary or categorical exposure, whose weights are
  not a ratio of densities. The default is accepted for every exposure
  type and ignored outside a continuous one.

  `.density` sits after `...` and so can only be supplied by name. The
  family the weights were built from is recorded on the result and read
  back with
  [`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md).

- numerator:

  How the marginal density that stabilizes a continuous exposure's
  weights is obtained, described under **Stabilization** in Details.
  Either `"marginal"`, the default, which reads the family `.density`
  names at the mean and standard deviation of `.exposure` over the rows
  the propensity model kept, or `"integrated"`, which averages the
  conditional density over the units on a grid spanning `.exposure` and
  interpolates the result back to each observed exposure.

  `"integrated"` describes a continuous exposure alone, and needs a
  conditional density to marginalize: it is refused with an error of
  class `propensity_numerator_error` for a binary or categorical
  exposure, with `stabilize = FALSE`, with a `stabilization_score`, and
  with any `.sigma`. The default is accepted for every exposure type and
  ignored outside a continuous one.

  `numerator` sits after `...` and so can only be supplied by name. The
  numerator the weights were built from is recorded on the result and
  read back with
  [`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md).

- .treated:

  **\[deprecated\]** Use `.focal_level` instead.

- .untreated:

  **\[deprecated\]** Use `.reference_level` instead.

- .propensity_col:

  Column to use when `.propensity` is a data frame with a binary or
  continuous exposure. Accepts a column name (quoted or unquoted) or
  numeric index. For a binary exposure the selected column is read as
  the probability of the resolved focal level; for a continuous one it
  is read as the conditional mean of the exposure.

  With no column named and a binary exposure, the exposure's levels
  drive the choice. When the data frame has a column named for every
  level of `.exposure`, the column named for the resolved focal level is
  used, wherever it sits in the frame. Columns named `.pred_<level>`,
  the shape a fitted tidymodels classification model predicts
  probabilities in, are matched by the level after the prefix. Otherwise
  the second column is used, or the only column when the frame has just
  one. Falling back to position after `.focal_level` or
  `.reference_level` was supplied warns and reports the column used,
  since the named level could not be matched to a column; the fallback
  is silent when no level was named and when the frame has a single
  column.

  With no column named and a continuous exposure, a frame of a single
  column is read as the conditional mean, and a frame of several columns
  is refused with an error of class
  `propensity_df_ambiguous_column_error`. There are no levels to read
  the names against, so nothing distinguishes one column of conditional
  means from another and there is no position worth trusting.

  Ignored for categorical exposures, where all columns are used.

## Value

A [`psw`](https://r-causal.github.io/propensity/reference/psw.md) vector
(a double vector with class `psw`) carrying these attributes:

- `estimand`: character, e.g. `"ate"`, `"att"`, `"uncensored"`.

- `stabilized`: logical, whether stabilization was applied.

- `stabilization_score`: double, the user-supplied stabilization score
  (one value, or one per observation), or `NULL` if none was supplied.

- `trimmed`: logical, whether the propensity scores were trimmed.

- `truncated`: logical, whether the propensity scores were truncated.

- `calibrated`: logical, whether the propensity scores were calibrated.

- `exposure_type`: character, the type of exposure the weights were
  built for, read back with
  [`exposure_type()`](https://r-causal.github.io/propensity/reference/exposure_type.md).

- `density_meta`: for a continuous exposure, whose weights are a ratio
  of densities, the record of that ratio: the density family, the
  numerator that stabilized the weights, and where the residual spread
  came from. Read back with
  [`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md),
  and `NULL` for a binary or categorical exposure.

## Details

### Exposure types

All weight functions support binary exposures. `wt_ate()` and
`wt_cens()` also support continuous exposures. All except `wt_cens()`
support categorical exposures; naming a categorical exposure to
`wt_cens()`, or handing it one detection reads as categorical, is an
error of class `causalgenerics_unsupported_exposure_type`.

- **Binary**: `.exposure` is a two-level vector (e.g., 0/1, logical, or
  a two-level factor). `.propensity` is a numeric vector of P(treatment
  \| X), or a matrix or data frame holding it in one of its columns.

- **Categorical**: `.exposure` is a factor or character vector with 3+
  levels. `.propensity` must be a matrix or data frame with one column
  per level, where rows sum to 1. Every weight function but `wt_cens()`
  also takes a fitted
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) and
  reads those columns off its fitted values, matching them to the levels
  of `.exposure` by name. A multinomial fit of only two levels reports a
  single probability and is read as a model of a binary exposure. Those
  per-level probabilities are the generalized propensity score of a
  multi-valued exposure (Imbens, 2000), and weighting by them is the
  multi-arm inverse probability estimator (Feng et al., 2012).

- **Continuous**: `.exposure` is a numeric vector. `.propensity` is a
  vector of conditional means (fitted values). Weights are a ratio of
  densities, whose family is chosen by `.density`, and are stabilized
  unless `stabilize = FALSE` asks otherwise.

- **Auto** (default): Detects the exposure type from `.exposure`.

### Stabilization

Stabilization multiplies the base weight by an estimate of P(A) (binary)
or f_A(A) (continuous), or by an estimate of either conditional on
variables a model of the exposure reads, reducing variance. When no
`stabilization_score` is supplied, that estimate is the marginal mean of
`.exposure` for a binary or categorical exposure, and for a continuous
exposure the marginal density of the family `.density` names, evaluated
at the population mean and standard deviation of `.exposure`.
Stabilization is supported for ATE and censoring weights (`wt_ate()` and
`wt_cens()`) alone.

Whether it is applied is read from the exposure type unless `stabilize`
says so outright. A continuous exposure is stabilized, since the
unstabilized ratio carries the exposure's own units and is not a form of
the weights the package recommends; a binary or categorical exposure is
left unstabilized, and `stabilize = TRUE` asks for it there.

For a continuous exposure, `numerator` chooses how that marginal density
is arrived at. `"marginal"`, the default, reads the family `.density`
names at the mean and standard deviation of `.exposure` over the rows
the propensity model kept, which are the rows the weights describe: a
model that dropped rows leaves `.propensity` missing at them, and those
rows carry no weight. Those two moments are parameters of the weights,
and
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
estimates them alongside the rest of its parameter vector, so the
standard errors account for the numerator having been estimated.
`"integrated"` marginalizes the conditional density numerically instead:
it averages \\f\_{A\|X}(t \mid X_i)\\ over the units at each of 50
points spanning `.exposure`, then interpolates that average back to each
observed exposure with a cubic spline. It estimates no parameters of its
own, and it is what WeightIt has done since version 2.0.0.

The two agree whenever the family is normal and the fitted conditional
means are themselves normal, since an average of normal densities
centered on normally distributed means is again a normal density. They
separate as those means grow skewed or bimodal, or as the family's tails
grow heavier than the normal's.

The marginal numerator is the default because it is the published
convention (Robins et al., 2000; Naimi et al., 2014) and because it held
up: in a simulation study run for this package, over linear, skewed,
heteroskedastic, heavy-tailed, and bimodal designs, marginalizing never
improved bias, root mean squared error, or weighted covariate balance by
more than Monte Carlo noise, including in the scenarios drawn to favor
it, and in the heavy-tailed cells it was worse: it inflated the root
mean squared error, and in some replicates the interpolated density came
back negative. The integrated numerator is here for agreement with
WeightIt, and for a conditional density whose marginal has no closed
form in the same family. Being read off an interpolation rather than a
formula, it can dip below zero where the density on the grid comes close
to it, which is refused with an error of class
`propensity_density_error` rather than turned into a negative weight.

The default numerator conditions on nothing, and `stabilization_score`
replaces it with one of your own, which leaves `numerator` nothing to
choose: the two cannot be supplied together. The case that pays is a
numerator conditioning on a variable the model the estimates are read
from also reads, such as an effect modifier. Fit that numerator and
evaluate it at the exposure each unit actually took, which is the shape
the default already has:

    num <- glm(A ~ V, data = dat, family = binomial())
    p <- fitted(num)
    wt_ate(ps_mod, stabilize = TRUE,
           stabilization_score = ifelse(dat$A == 1, p, 1 - p))

Passing `fitted(num)` on its own is P(A = 1 \| V) for every unit, which
is not the probability of the exposure an untreated unit took, so the
untreated weights would carry the wrong numerator. Whether a numerator
may condition on `V` at all is a question about the reported model
rather than about the weights; see **Effect modification** in
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).

There is a third way of writing the same thing: pass the fitted
numerator model to `stabilize` itself, and what it estimates becomes the
numerator.

    num <- glm(A ~ V, data = dat, family = binomial())
    wt_ate(ps_mod, stabilize = num)

    num <- lm(A ~ V, data = dat)
    wt_ate(dose_mod, stabilize = num)

    num <- nnet::multinom(Z ~ V, data = dat, trace = FALSE)
    wt_ate(ps_mat, Z, exposure_type = "categorical", stabilize = num)

For a binary exposure the weights are then \\P(A = a_i \mid V_i) / f(a_i
\mid X_i)\\: the model's fitted probability read at the level each unit
actually took, which is \\p_i\\ for a unit with \\A_i = 1\\ and \\1 -
p_i\\ for a unit with \\A_i = 0\\, over the same denominator the
unstabilized weights divide by. Writing it as the model rather than as
the numbers is what keeps the untreated units from carrying \\P(A = 1
\mid V)\\, which is the mistake the `stabilization_score` form has to be
written around by hand.

For a continuous exposure they are \\f(A \mid V) / f(A \mid X)\\: the
family `.density` names, read at the numerator model's fitted mean and
the root mean square of its residuals, over the same family read at the
propensity score model's.

For a categorical exposure they are \\P(Z = z_i \mid V_i) / f(z_i \mid
X_i)\\: the multinomial model's fitted probability of the level each
unit took, gathered from the column named for that level, over the same
denominator again. The fit has to report a probability for every level
the exposure takes, and its columns are read by name rather than by
position, so a fit that declares those levels in another order builds
the same weights here.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
asks for the propensity score model's own level order in addition, since
the block it stacks reads the coefficients positionally. Fitting a
multinomial model of the exposure on the stabilization terms is what the
ipw package (van der Wal & Geskus, 2011) does for a multi-valued
exposure as well.

The same contract governs all three, and it is a statement about the
model the estimates will be read from rather than about the weights. A
numerator conditioning on `V` builds a pseudo-population in which `V`
still predicts the exposure, so confounding by `V` remains there and
what is estimated is the effect conditional on `V` rather than the
marginal one. The outcome model, or the marginal structural model, must
therefore include `V` (Robins et al., 2000; Cole & Hernán, 2008; Hernán
& Robins, 2020, Chapter 12). The limiting case says it from the other
end: with `V` the whole set of confounders the numerator is the
denominator, the weights collapse toward 1, and the weighting has
nothing left to do that the model is not already doing. The case that
pays is a numerator conditioning on an effect modifier the reported
model reads anyway.

Handing the model over rather than the numbers it evaluates to is what
lets
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
estimate it: the model's own estimating equations join the stacked
system, so the standard errors account for the numerator having been
fitted, where a `stabilization_score` is carried as a known constant.
That accounting has something to say only where the reported model is
not saturated in the variables the numerator reads. A numerator of a
binary or categorical exposure on an effect modifier is constant within
each cell of the modifier and the exposure, so a model saturated in
those cells fits the same coefficients with the numerator and without
it, and the standard errors do not move either.

What a numerator model is held to is its shape rather than its
provenance: its class and family, the levels it declares for a
categorical exposure, and one fitted value per observation. A fit to a
different dataset of the same length passes every one of those and
multiplies in a numerator belonging to other rows, without a word, so
fit the numerator model on the data the weights are being built for.

Only `wt_ate()` and `wt_cens()` stabilize at all, which follows from
what the other estimands are rather than from anything missing here. An
`att`, `atu`, `atm`, `ato`, or entropy weight already carries a tilting
function where a numerator would go, so an estimate of \\P(A)\\ or \\P(A
\mid V)\\ put there would move the population the weights target rather
than how variable they are. Asking for it is an error.

[`numerator_model()`](https://r-causal.github.io/propensity/reference/numerator_model.md)
reads the model back off the weights.

### Handling extreme weights

Extreme weights signal positivity violations, poor model fit, or limited
overlap. You can address them by:

- Choosing an overlap-focused estimand (`wt_ato()`, `wt_atm()`,
  `wt_entropy()`), which down-weight units in regions of poor overlap.

- Trimming
  ([`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md))
  or truncating
  ([`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md))
  propensity scores before computing weights.

- Calibrating weights with
  [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md).

- Stabilizing ATE and censoring weights: a continuous exposure is
  stabilized by default, and a binary or categorical one is with
  `stabilize = TRUE`.

- Choosing a heavier-tailed `.density` for a continuous exposure, which
  holds down the weight of a unit whose exposure the propensity score
  model fits poorly.

See the [halfmoon](https://CRAN.R-project.org/package=halfmoon) package
for weight diagnostics and visualization.

### Propensity scores at 0 and 1

Propensity scores must lie strictly inside (0, 1), since a score at
either endpoint leaves the weight undefined. A score of exactly 0 or 1
is refused with an error of class `propensity_range_error`. For a
categorical exposure that rule reads every cell of the matrix rather
than only the column holding each unit's observed level, which is the
same open interval the binary path enforces on the single score and its
complement.

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
applies a narrower rule to the propensity scores it rebuilds from its
own propensity score model. For a categorical exposure it refuses only a
probability of exactly zero at the level a unit was actually assigned,
since that alone leaves the unit's own weight undefined, and a column
for a level the unit was not assigned may reach zero without the fit
being refused. The difference is deliberate: the weight functions can
refuse a matrix that holds nothing
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
would have objected to. The weights still have to be built before
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
can be reached, so the refusal here is the one to resolve.

[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) reaches
the endpoints readily. Under separation the softmax puts the probability
at a unit's assigned level at exactly 1 in double precision, and the
columns for the other levels can underflow to exactly 0.
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
is the repair for such a matrix: its categorical matrix method reads the
closed interval \\\[0, 1\]\\, so it accepts a cell of exactly 0 or 1,
pins every score below its threshold up to it, and renormalizes each row
to sum to 1. What it returns lies strictly inside (0, 1) and weights
like any other truncated matrix. That relaxation belongs to the
categorical matrix method alone:
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
sets extreme scores to missing and gains nothing from an endpoint, so it
holds the open interval, and so does the binary vector path of
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md).
If the bounds a truncation computes for itself leave a score at an
endpoint anyway, which `method = "pctl"` can do when the lower quantile
of the scores is itself zero, the call is refused and names the bound it
computed; supply `lower` and `upper` explicitly instead. Refitting the
propensity score model so that it does not separate remains the
alternative.

### Weight formulas

#### Binary exposures

For binary treatments (\\A \in \\0, 1\\\\), with propensity score \\e(X)
= P(A=1 \mid X)\\:

- **ATE**: \\w = \frac{A}{e(X)} + \frac{1-A}{1-e(X)}\\

- **ATT**: \\w = A + \frac{(1-A) \cdot e(X)}{1-e(X)}\\

- **ATU**: \\w = \frac{A \cdot (1-e(X))}{e(X)} + (1-A)\\

- **ATM**: \\w = \frac{\min(e(X), 1-e(X))}{A \cdot e(X) + (1-A) \cdot
  (1-e(X))}\\

- **ATO**: \\w = A \cdot (1-e(X)) + (1-A) \cdot e(X)\\

- **Entropy**: \\w = \frac{h(e(X))}{A \cdot e(X) + (1-A) \cdot
  (1-e(X))}\\, where \\h(e) = -\[e \log(e) + (1-e) \log(1-e)\]\\

The entropy weight tilts the propensity score by \\h(e)\\, the entropy
of the score itself, in the sense of Zhou, Matsouaka, and Thomas (2020).
It is not entropy balancing, which solves for weights that satisfy exact
covariate moment constraints rather than tilting a fitted propensity
score.

#### Continuous exposures

Weights use the density ratio \\w = f_A(A) / f\_{A\|X}(A \mid X)\\,
where \\f_A\\ is the marginal density of the exposure and \\f\_{A\|X}\\
is the conditional density the propensity model estimates. Only
`wt_ate()` and `wt_cens()` support continuous exposures.

Both densities are read on a standardized residual. The conditional
density is evaluated at \\z_i = (A_i - \hat{A}\_i) / \sigma\\, where
\\\hat{A}\_i\\ is the fitted conditional mean in `.propensity` and
\\\sigma\\ is the residual spread; the marginal density is evaluated at
\\z^A_i = (A_i - \bar{A}) / s_A\\, where \\\bar{A}\\ and \\s_A\\ are the
mean and standard deviation of `.exposure` over the rows the propensity
model kept. Both moments are read over those rows whether the spread was
pooled from the model's own residuals or supplied through `.sigma`, so
the two halves of the ratio describe one set of units. Each density is
then divided by the spread that standardized it, the Jacobian of that
change of variable, which returns both to the exposure's own units so
that each integrates to one:

\$\$w_i = \frac{g(z^A_i) / s_A}{g(z_i) / \sigma}\$\$

`.density` chooses the family \\g\\:

|  |  |
|----|----|
| `.density` | \\g(z)\\ |
| `"normal"`, [`dens_normal()`](https://r-causal.github.io/propensity/reference/dens_normal.md) | the standard normal density, the default |
| `"laplace"`, [`dens_laplace()`](https://r-causal.github.io/propensity/reference/dens_normal.md) | \\\exp(-\lvert z \rvert) / 2\\ |
| [`dens_t()`](https://r-causal.github.io/propensity/reference/dens_normal.md) | Student's t with `df` degrees of freedom |
| `"kernel"`, [`dens_kernel()`](https://r-causal.github.io/propensity/reference/dens_normal.md) | a kernel density estimate of the standardized residuals |
| [`dens_fn()`](https://r-causal.github.io/propensity/reference/dens_normal.md), or a bare function | a density you write yourself |

The numerator and the denominator take the same family, so the choice
describes the shape of the residuals rather than the estimand. A heavier
tail than the normal's, such as
[`dens_t()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
or `"laplace"`, holds down the weight of a unit whose exposure the model
fits poorly, and `"kernel"` assumes no shape at all, at the cost of a
density that is not a smooth function of the model's parameters. A
`stabilization_score` replaces the numerator outright and
`stabilize = FALSE` leaves it out, so under either the family describes
the denominator alone. Under `numerator = "integrated"` it describes
both again, the numerator there being the same conditional density
averaged over the units rather than the same family read at the
exposure's own moments; see **Stabilization**.

A unit whose exposure or fitted conditional mean is missing has no
standardized residual, and so no weight: the density is asked only about
the units that have one, and a kernel is fit on those. A propensity
score
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
set aside is the ordinary way this arises.

The conditional means in `.propensity` can be read from the model that
fit them. `wt_ate()` and `wt_cens()` take a fitted
[`lm()`](https://rdrr.io/r/stats/lm.html), a
[`glm()`](https://rdrr.io/r/stats/glm.html) fit with
[`gaussian()`](https://rdrr.io/r/stats/family.html) under any of its
links, an [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) fit
with [`gaussian()`](https://rdrr.io/r/stats/family.html), and a
[`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html), and read the
exposure from the model's response when `.exposure` is not supplied.
Each of these reports its conditional mean on the scale of the exposure,
so a log, inverse, or square root link never has to be undone. Any other
subclass of [`lm()`](https://rdrr.io/r/stats/lm.html) is read as an
`lm`, because the fitted conditional means are all the weights ask of a
model.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
names the classes it supports instead of accepting a subclass on those
terms: it stacks each model's own estimating equations, which it can
only write for a score it knows.

A family whose spread changes with its fitted values, such as
[`poisson()`](https://rdrr.io/r/stats/family.html) or
`quasi(variance = "mu")`, describes a different density for every unit,
which a single spread cannot stand in for, and is refused with an error
of class `propensity_model_family_error`. `quasi(variance = "constant")`
is the gaussian variance under another name and is accepted. A model of
a binary exposure is held to the same rule from the other side: only
[`binomial()`](https://rdrr.io/r/stats/family.html) and
[`quasibinomial()`](https://rdrr.io/r/stats/family.html) fit the
probability those weights divide by, so a model of a conditional mean is
refused there whatever its fitted values happen to fall on, a linear
probability model included.

\\\sigma\\ is the pooled residual spread \\\sqrt{\mathrm{mean}((A -
\hat{A})^2)}\\ unless `.sigma` supplies a single standard deviation, or
one for each unit, which the model methods do not do on their own.
[`dens_t()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
offers a second estimator, the maximum likelihood scale of the t itself,
through `sigma_method = "mle"`; see **The spread of a t density** in its
documentation. The marginal density is spread by \\s_A\\ whatever the
conditional spread was arrived at, the spread of the exposure itself
being no business of the conditional model's.

Every model class is spread that same way,
[`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) included. `rlm`
reports a robust scale estimate of its own in `fit$s`, which resists the
extreme residuals rather than pooling all of them; pass `.sigma = fit$s`
to be spread by it instead. That is rarely what you want on the data a
robust fit is used for: the scale resists the very residuals the
conditional density is meant to spread over, so the ratio can
concentrate on a handful of units.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
stacks a robust fit's own score with that scale held fixed as a known
constant, and reads `.sigma = fit$s` as a fixed spread in the same way,
propagating the uncertainty of neither.

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
models the conditional density with a single conditional variance,
estimated jointly with the rest of the parameter vector, by the moment
the pooled spread is the root of or, for a density built with
`sigma_method = "mle"`, by the score of the t. A single `.sigma` is
taken instead as a known constant: the weights are rebuilt at the number
that was supplied, and the stacked system carries none of that number's
uncertainty. An observation-level `.sigma` is a different function of
the data with no counterpart in that estimating equation, and is refused
before anything is solved, with an error of class
`propensity_ipw_sigma_error`. Build weights with the pooled default, or
with one number, when the outcome model is headed for
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).

A conditional density of exactly zero at a unit's own exposure is the
continuous counterpart of a propensity score of exactly 0 or 1: the
weight it leaves is infinite, and every estimate built from it is
undefined. It is refused with an error of class
`propensity_density_error` naming the observations it happened at. A
light-tailed family at an outlying residual and a kernel read past the
range it was fit on both reach it, so the remedy is a family with
heavier tails, such as
[`dens_t()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
or
[`dens_kernel()`](https://r-causal.github.io/propensity/reference/dens_normal.md).
A zero in the numerator is not refused: a marginal density of zero is a
weight of zero, which is a weight like any other.

Stabilized weights from the normal family have a finite second moment
only while the marginal variance of the exposure stays below twice the
variance of the conditional density, \\\mathrm{Var}(A) \< 2\sigma^2\\.
Past that boundary the weights have no finite variance, and estimates
built from them are erratic however large the sample is. Reaching it is
reported with a message of class `propensity_density_variance_message`;
the weights are still returned, and the report is informational because
a well-fitting model reaches the boundary in ordinary use. Explaining
more of the exposure lowers \\\sigma^2\\ while the marginal variance is
fixed by the data, so a better fit tightens the boundary rather than
escaping it; a heavier-tailed family sits at a different one, and an
unstabilized weight has no boundary at all. The boundary is a property
of the normal family read against the exposure's own marginal density,
so nothing is reported for an unstabilized weight, a supplied
`stabilization_score`, `numerator = "integrated"`, another family, or an
observation-level `.sigma`.

#### Categorical exposures

For \\K\\-level treatments, weights take the tilting-function form \\w_i
= h(\mathbf{e}\_i) / e\_{i,Z_i}\\, where \\e\_{i,Z_i}\\ is the
propensity for unit \\i\\'s observed level and \\h(\cdot)\\ depends on
the estimand:

- **ATE**: \\h(\mathbf{e}) = 1\\

- **ATT**: \\h(\mathbf{e}) = e\_{\text{focal}}\\

- **ATU**: \\h(\mathbf{e}) = 1 - e\_{\text{focal}}\\

- **ATM**: \\h(\mathbf{e}) = \min(e_1, \ldots, e_K)\\

- **ATO**: \\h(\mathbf{e}) = \bigl(\sum_k 1/e_k\bigr)^{-1}\\

- **Entropy**: \\h(\mathbf{e}) = -\sum_k e_k \log(e_k)\\

## References

Barrett, M., D'Agostino McGowan, L., & Gerke, T. *Causal Inference in
R*. <https://www.r-causal.org/>

Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the
propensity score in observational studies for causal effects.
*Biometrika*, 70(1), 41–55.

Li, L., & Greene, T. (2013). A weighting analogue to pair matching in
propensity score analysis. *The International Journal of Biostatistics*,
9(2), 215–234. (ATM weights)

Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates
via propensity score weighting. *Journal of the American Statistical
Association*, 113(521), 390–400. (ATO weights)

Zhou, Y., Matsouaka, R. A., & Thomas, L. (2020). Propensity score
weighting under limited overlap and model misspecification. *Statistical
Methods in Medical Research*, 29(12), 3721–3756. (Entropy weights)

Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous
treatments. In *Applied Bayesian Modeling and Causal Inference from
Incomplete-Data Perspectives* (pp. 73–84).

Robins, J. M., Hernán, M. A., & Brumback, B. (2000). Marginal structural
models and causal inference in epidemiology. *Epidemiology*, 11(5),
550–560. (Stabilized weights and the conditional numerator)

Cole, S. R., & Hernán, M. A. (2008). Constructing inverse probability
weights for marginal structural models. *American Journal of
Epidemiology*, 168(6), 656–664. (Numerator covariates in the weighted
model)

Hernán, M. A., & Robins, J. M. (2020). *Causal Inference: What If*.
Chapman & Hall/CRC. (Chapter 12, stabilized weights and the conditional
numerator)

Imbens, G. W. (2000). The role of the propensity score in estimating
dose-response functions. *Biometrika*, 87(3), 706–710. (Generalized
propensity score of a multi-valued exposure)

Feng, P., Zhou, X.-H., Zou, Q.-M., Fan, M.-Y., & Li, X.-S. (2012).
Generalized propensity score for estimating the average treatment effect
of multiple treatments. *Statistics in Medicine*, 31(7), 681–697.

van der Wal, W. M., & Geskus, R. B. (2011). ipw: An R package for
inverse probability weighting. *Journal of Statistical Software*,
43(13), 1–23. (Multinomial numerator model for a multi-valued exposure)

Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014).
Constructing inverse probability weights for continuous exposures: a
comparison of methods. *Epidemiology*, 25(2), 292–299.

Austin, P. C., & Stuart, E. A. (2015). Moving towards best practice when
using inverse probability of treatment weighting (IPTW). *Statistics in
Medicine*, 34(28), 3661–3679.

Greifer, N. *WeightIt: Weighting for Covariate Balance in Observational
Studies*. R package. <https://CRAN.R-project.org/package=WeightIt> (the
integrated numerator, which WeightIt has computed for continuous
treatments since version 2.0.0)

## See also

- [`psw()`](https://r-causal.github.io/propensity/reference/psw.md) for
  the returned weight vector class.

- [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
  [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
  [`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md),
  and
  [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
  for modifying propensity scores before weighting.

- [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  for inverse-probability-weighted estimation of causal effects.

- [`ps_tilt()`](https://r-causal.github.io/propensity/reference/ps_tilt.md)
  for the tilting function each weight divides by the propensity score
  of the received exposure level, which also standardizes a
  g-computation estimate to the same target population.

## Examples

``` r
# -- Binary exposure, numeric propensity scores ----------------------
set.seed(123)
ps <- runif(100, 0.1, 0.9)
trt <- rbinom(100, 1, ps)

wt_ate(ps, trt)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[100]>
#>   [1] 1.492675 1.368655 1.745754 5.165661 1.173194 7.328950 2.094172 1.228599
#>   [9] 1.847923 1.870179 7.433103 1.861044 1.557495 2.262990 1.223002 1.219720
#>  [17] 1.422212 7.482363 1.568225 1.157940 1.232086 1.528485 1.632905 1.116800
#>  [25] 1.601115 3.001420 1.868276 1.738182 1.495501 1.278267 1.148871 5.612908
#>  [33] 2.878230 3.793252 1.135965 2.073670 3.410265 3.661309 2.820518 1.399190
#>  [41] 1.272653 1.759439 1.757406 1.653101 4.505402 1.267499 1.401399 1.896705
#>  [49] 1.455134 1.271840 1.158299 1.830697 1.352924 1.246136 1.822296 1.360961
#>  [57] 1.253173 1.423191 1.225436 1.665474 1.582048 1.213405 2.455942 1.469523
#>  [65] 1.330297 1.847790 1.336806 1.333490 1.359668 1.824369 1.421302 1.657339
#>  [73] 3.013373 1.111729 2.082235 1.381397 1.677439 1.694293 2.621656 1.232906
#>  [81] 3.391031 1.576182 2.303524 1.368819 1.222930 1.811313 1.126170 1.227836
#>  [89] 5.240410 4.165936 1.257160 1.606473 2.667996 1.598960 2.806635 1.333605
#>  [97] 1.377723 1.211939 1.899058 2.037508
wt_att(ps, trt)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = att}[100]>
#>   [1] 0.4926755 1.0000000 0.7457538 4.1656608 1.0000000 1.0000000 1.0941724
#>   [8] 1.0000000 1.0000000 0.8701789 6.4331026 0.8610445 1.0000000 1.2629898
#>  [15] 0.2230018 1.0000000 0.4222125 1.0000000 0.5682254 1.0000000 1.0000000
#>  [22] 1.0000000 1.0000000 1.0000000 1.0000000 2.0014200 1.0000000 1.0000000
#>  [29] 0.4955011 0.2782671 1.0000000 4.6129081 1.8782298 2.7932516 0.1359647
#>  [36] 1.0000000 2.4102647 1.0000000 1.0000000 0.3991897 0.2726533 0.7594392
#>  [43] 0.7574058 0.6531012 1.0000000 0.2674992 0.4013989 0.8967053 0.4551341
#>  [50] 1.0000000 0.1582988 0.8306973 1.0000000 0.2461361 1.0000000 0.3609610
#>  [57] 0.2531726 1.0000000 1.0000000 0.6654737 1.0000000 0.2134045 1.0000000
#>  [64] 0.4695226 1.0000000 0.8477904 1.0000000 1.0000000 1.0000000 0.8243692
#>  [71] 1.0000000 1.0000000 2.0133726 0.1117285 1.0000000 0.3813969 0.6774393
#>  [78] 1.0000000 1.0000000 0.2329063 1.0000000 1.0000000 1.0000000 1.0000000
#>  [85] 0.2229300 0.8113126 1.0000000 1.0000000 4.2404103 1.0000000 0.2571604
#>  [92] 1.0000000 1.0000000 1.0000000 1.0000000 0.3336052 1.0000000 0.2119390
#>  [99] 0.8990583 1.0375079
wt_atu(ps, trt)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = atu}[100]>
#>   [1] 1.0000000 0.3686554 1.0000000 1.0000000 0.1731942 6.3289497 1.0000000
#>   [8] 0.2285990 0.8479233 1.0000000 1.0000000 1.0000000 0.5574953 1.0000000
#>  [15] 1.0000000 0.2197205 1.0000000 6.4823626 1.0000000 0.1579396 0.2320863
#>  [22] 0.5284847 0.6329051 0.1167996 0.6011153 1.0000000 0.8682760 0.7381824
#>  [29] 1.0000000 1.0000000 0.1488715 1.0000000 1.0000000 1.0000000 1.0000000
#>  [36] 1.0736701 1.0000000 2.6613092 1.8205180 1.0000000 1.0000000 1.0000000
#>  [43] 1.0000000 1.0000000 3.5054016 1.0000000 1.0000000 1.0000000 1.0000000
#>  [50] 0.2718404 1.0000000 1.0000000 0.3529239 1.0000000 0.8222956 1.0000000
#>  [57] 1.0000000 0.4231912 0.2254357 1.0000000 0.5820478 1.0000000 1.4559422
#>  [64] 1.0000000 0.3302967 1.0000000 0.3368064 0.3334905 0.3596676 1.0000000
#>  [71] 0.4213022 0.6573389 1.0000000 1.0000000 1.0822347 1.0000000 1.0000000
#>  [78] 0.6942927 1.6216558 1.0000000 2.3910308 0.5761821 1.3035242 0.3688192
#>  [85] 1.0000000 1.0000000 0.1261698 0.2278362 1.0000000 3.1659355 1.0000000
#>  [92] 0.6064733 1.6679958 0.5989600 1.8066347 1.0000000 0.3777228 1.0000000
#>  [99] 1.0000000 1.0000000
wt_atm(ps, trt)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = atm}[100]>
#>   [1] 0.4926755 0.3686554 0.7457538 1.0000000 0.1731942 1.0000000 1.0000000
#>   [8] 0.2285990 0.8479233 0.8701789 1.0000000 0.8610445 0.5574953 1.0000000
#>  [15] 0.2230018 0.2197205 0.4222125 1.0000000 0.5682254 0.1579396 0.2320863
#>  [22] 0.5284847 0.6329051 0.1167996 0.6011153 1.0000000 0.8682760 0.7381824
#>  [29] 0.4955011 0.2782671 0.1488715 1.0000000 1.0000000 1.0000000 0.1359647
#>  [36] 1.0000000 1.0000000 1.0000000 1.0000000 0.3991897 0.2726533 0.7594392
#>  [43] 0.7574058 0.6531012 1.0000000 0.2674992 0.4013989 0.8967053 0.4551341
#>  [50] 0.2718404 0.1582988 0.8306973 0.3529239 0.2461361 0.8222956 0.3609610
#>  [57] 0.2531726 0.4231912 0.2254357 0.6654737 0.5820478 0.2134045 1.0000000
#>  [64] 0.4695226 0.3302967 0.8477904 0.3368064 0.3334905 0.3596676 0.8243692
#>  [71] 0.4213022 0.6573389 1.0000000 0.1117285 1.0000000 0.3813969 0.6774393
#>  [78] 0.6942927 1.0000000 0.2329063 1.0000000 0.5761821 1.0000000 0.3688192
#>  [85] 0.2229300 0.8113126 0.1261698 0.2278362 1.0000000 1.0000000 0.2571604
#>  [92] 0.6064733 1.0000000 0.5989600 1.0000000 0.3336052 0.3777228 0.2119390
#>  [99] 0.8990583 1.0000000
wt_ato(ps, trt)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ato}[100]>
#>   [1] 0.3300620 0.2693559 0.4271815 0.8064139 0.1476262 0.8635548 0.5224844
#>   [8] 0.1860648 0.4588520 0.4652918 0.8654667 0.4626673 0.3579435 0.5581067
#>  [15] 0.1823397 0.1801400 0.2968702 0.8663524 0.3623366 0.1363971 0.1883685
#>  [22] 0.3457573 0.3875945 0.1045842 0.3754354 0.6668244 0.4647472 0.4246864
#>  [29] 0.3313278 0.2176909 0.1295806 0.8218392 0.6525642 0.7363739 0.1196909
#>  [36] 0.5177632 0.7067676 0.7268737 0.6454552 0.2853006 0.2142400 0.4316371
#>  [43] 0.4309795 0.3950764 0.7780442 0.2110449 0.2864273 0.4727700 0.3127781
#>  [50] 0.2137378 0.1366649 0.4537601 0.2608601 0.1975194 0.4512416 0.2652251
#>  [57] 0.2020253 0.2973537 0.1839637 0.3995702 0.3679078 0.1758725 0.5928243
#>  [64] 0.3195069 0.2482880 0.4588131 0.2519485 0.2500884 0.2645261 0.4518654
#>  [71] 0.2964199 0.3966231 0.6681459 0.1004998 0.5197467 0.2760951 0.4038532
#>  [78] 0.4097832 0.6185617 0.1889083 0.7051044 0.3655555 0.5658826 0.2694433
#>  [85] 0.1822917 0.4479142 0.1120344 0.1855591 0.8091752 0.7599579 0.2045566
#>  [92] 0.3775185 0.6251868 0.3745935 0.6437014 0.2501529 0.2741646 0.1748760
#>  [99] 0.4734232 0.5092044
wt_entropy(ps, trt)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = entropy}[100]>
#>   [1] 0.9466884 0.7974021 1.1914845 2.5383087 0.4910630 2.9202760 1.4494516
#>   [8] 0.5903003 1.2746181 1.2917997 2.9354392 1.2847853 1.0158386 1.5532690
#>  [15] 0.5808315 0.5752272 0.8649741 2.9425314 1.0267968 0.4612870 0.5961433
#>  [22] 0.9855373 1.0902253 0.3741728 1.0595945 1.9101181 1.2903427 1.1850225
#>  [29] 0.9498151 0.6697735 0.4429918 2.6301699 1.8588898 2.1880067 0.4161138
#>  [36] 1.4360497 2.0632803 2.1467853 1.8339424 0.8365617 0.6611694 1.2030532
#>  [43] 1.2013433 1.1091726 2.3850345 0.6531902 0.8393281 1.3118818 0.9040828
#>  [50] 0.6599161 0.4620024 1.2611029 0.7765170 0.6192602 1.2544407 0.7872501
#>  [57] 0.6305931 0.8661618 0.5849629 1.1205924 1.0407229 0.5643262 1.6597604
#>  [64] 0.9206519 0.7455600 1.2745145 0.7545814 0.7499980 0.7855318 1.2560894
#>  [71] 0.8638679 1.1130997 1.9149498 0.3626234 1.4416708 0.8139569 1.1315051
#>  [78] 1.1466626 1.7427820 0.5975110 2.0565825 1.0348389 1.5766261 0.7976169
#>  [85] 0.5807093 1.2456605 0.3950015 0.5890166 2.5542633 2.2959657 0.6369461
#>  [92] 1.0648287 1.7647932 1.0574806 1.8278454 0.7501570 0.8092153 0.5617751
#>  [99] 1.3136430 1.4119476

# Stabilized ATE weights (reduces variance)
wt_ate(ps, trt, stabilize = TRUE)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate; stabilized}[100]>
#>   [1] 0.7612645 0.6706411 0.8903344 2.6344870 0.5748651 3.5911853 1.0680279
#>   [8] 0.6020135 0.9054824 0.9537912 3.7908823 0.9491327 0.7631727 1.1541248
#>  [15] 0.6237309 0.5976630 0.7253284 3.6663577 0.7997950 0.5673904 0.6037223
#>  [22] 0.7489575 0.8001235 0.5472318 0.7845465 1.5307242 0.9154552 0.8517094
#>  [29] 0.7627055 0.6519162 0.5629470 2.8625831 1.4678972 1.9345583 0.5793420
#>  [36] 1.0160984 1.7392350 1.7940415 1.3820538 0.7135867 0.6490532 0.8973140
#>  [43] 0.8962770 0.8430816 2.2076468 0.6464246 0.7147134 0.9673197 0.7421184
#>  [50] 0.6232018 0.5907324 0.9336556 0.6629327 0.6355294 0.8929248 0.6940901
#>  [57] 0.6391180 0.6973637 0.6004635 0.8493916 0.7752034 0.6188363 1.2034117
#>  [64] 0.7494566 0.6518454 0.9423731 0.6550351 0.6534103 0.6662371 0.9304283
#>  [71] 0.6964381 0.8120960 1.5368200 0.5669815 1.0202950 0.7045124 0.8554940
#>  [78] 0.8302034 1.2846113 0.6287822 1.6616051 0.7723292 1.1287269 0.6707214
#>  [85] 0.6236943 0.9237694 0.5518232 0.6016397 2.6726093 2.0413084 0.6411518
#>  [92] 0.7871719 1.3073180 0.7834904 1.3752510 0.6801387 0.6750841 0.6180889
#>  [99] 0.9685198 1.0391291

# Inspect the result
w <- wt_ate(ps, trt)
#> ℹ Treating `.exposure` as binary
estimand(w)
#> [1] "ate"
summary(w)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.112   1.317   1.591   2.044   2.047   7.482 

# -- Overlap-focused estimands handle extreme PS better --------------
ps_extreme <- c(0.01, 0.02, 0.98, 0.99, rep(0.5, 4))
trt_extreme <- c(0, 0, 1, 1, 0, 1, 0, 1)

max(wt_ate(ps_extreme, trt_extreme))
#> ℹ Treating `.exposure` as binary
#> [1] 2
max(wt_ato(ps_extreme, trt_extreme))
#> ℹ Treating `.exposure` as binary
#> [1] 0.5

# -- From a fitted GLM -----------------------------------------------
x1 <- rnorm(100)
x2 <- rnorm(100)
trt2 <- rbinom(100, 1, plogis(0.5 * x1 + 0.3 * x2))
ps_model <- glm(trt2 ~ x1 + x2, family = binomial)

# Exposure is extracted from the model automatically
wt_ate(ps_model)
#> ℹ Using exposure variable "trt2" from the propensity score model
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[100]>
#>   [1] 1.772299 2.816364 1.919886 2.059945 1.590514 2.019637 1.592626 3.675866
#>   [9] 1.731064 1.526937 2.194737 1.585679 1.604308 2.035997 1.836701 3.038426
#>  [17] 1.926720 2.527851 2.956648 1.456686 1.951071 2.308830 2.122857 1.859348
#>  [25] 1.333779 1.789153 2.089527 2.008161 1.809376 1.845397 1.236823 2.436497
#>  [33] 1.817852 2.352419 1.300776 1.637110 3.144659 2.924006 4.687625 1.406029
#>  [41] 2.396310 2.163061 2.510037 1.356974 1.430312 2.334509 3.102993 2.202745
#>  [49] 1.310198 1.607962 2.624988 2.535210 1.793952 1.647484 2.569296 2.061068
#>  [57] 3.046296 1.936586 1.517183 1.624578 2.826262 2.825176 1.554039 1.098856
#>  [65] 1.637046 1.507218 2.581387 1.554926 2.385541 1.698365 1.811373 1.670399
#>  [73] 1.760472 1.365580 1.911171 1.539045 2.208360 1.724933 1.617692 1.874702
#>  [81] 1.406166 3.813823 1.982021 2.553269 1.982884 1.957207 1.440986 1.560568
#>  [89] 1.603633 1.885138 1.584638 2.218637 1.630055 1.585865 1.872610 1.274698
#>  [97] 1.393716 1.394134 1.800220 1.778728

# -- Continuous exposure from a fitted model -------------------------
dose <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(100)
dose_model <- lm(dose ~ x1 + x2)

# The exposure is read from the model, and a continuous exposure is
# stabilized without being asked
w_dose <- wt_ate(dose_model)
#> ℹ Using exposure variable "dose" from the propensity score model
#> ℹ Treating `.exposure` as continuous

# A tail heavier than the normal's holds down the weight of a unit whose
# dose the model fits poorly
wt_ate(dose_model, .density = dens_t(df = 4))
#> ℹ Using exposure variable "dose" from the propensity score model
#> ℹ Treating `.exposure` as continuous
#> <psw{estimand = ate; stabilized}[100]>
#>   [1] 1.8646342 0.9733784 0.8752609 0.8688423 0.7479327 0.7798293 0.9016771
#>   [8] 0.7597835 0.6310464 0.6911070 0.7999152 0.7082642 0.8020596 0.7990152
#>  [15] 0.5044498 0.7911563 0.7468465 0.9357483 0.9238356 0.8522952 0.9048531
#>  [22] 0.5633915 0.6523243 1.1495582 0.3374106 0.7535244 1.3004884 0.7506499
#>  [29] 0.7409988 1.7593253 0.8070058 0.7424372 0.8403367 1.0590076 1.0785524
#>  [36] 0.3059582 1.0582182 1.3140401 0.8689788 0.9671610 1.1106228 1.4415967
#>  [43] 1.6841604 0.7055659 0.5493660 2.4236372 0.5417968 0.4050666 2.2863800
#>  [50] 0.5851142 0.7353432 1.6511412 1.0393921 0.7953566 1.1868726 1.0917254
#>  [57] 0.8642301 1.1409825 1.0153774 0.7270439 0.4827268 0.9199718 1.2137686
#>  [64] 0.4873167 0.5196829 1.0329573 0.7868673 0.6905860 0.6315537 0.9681556
#>  [71] 1.6544586 0.9066800 0.8659805 0.1618948 0.5878388 0.7318754 1.2574093
#>  [78] 0.8112021 0.9306485 0.9048382 0.8503056 0.6369130 0.9173769 0.7515608
#>  [85] 0.8473410 0.9312568 0.5513774 1.4042470 2.4415982 1.3000206 0.8634111
#>  [92] 1.0266033 1.0473099 0.8128075 0.3222650 0.4487284 0.9522881 1.0692692
#>  [99] 0.9326978 0.3828823
#> density:   t(df = 4)
#> numerator: marginal
#> sigma:     pooled

# The stabilizing numerator can marginalize the conditional density over the
# units instead of reading the same family at the exposure's own moments
w_int <- wt_ate(
  dose_model,
  .density = dens_t(df = 4),
  numerator = "integrated"
)
#> ℹ Using exposure variable "dose" from the propensity score model
#> ℹ Treating `.exposure` as continuous

# It can also be the conditional density a second model estimates, which is
# worth having when the model the estimates are read from also reads the
# variables that model conditions on
num_model <- lm(dose ~ x1)
w_cond <- wt_ate(dose_model, stabilize = num_model)
#> ℹ Using exposure variable "dose" from the propensity score model
#> ℹ Treating `.exposure` as continuous

# What each set of weights records about the ratio it is
density_meta(w_dose)
#> density:   normal
#> numerator: marginal
#> sigma:     pooled
density_meta(w_int)
#> density:   t(df = 4)
#> numerator: integrated
#> sigma:     pooled

# -- Data frame input ------------------------------------------------
ps_df <- data.frame(
  control = c(0.9, 0.7, 0.3, 0.1),
  treated = c(0.1, 0.3, 0.7, 0.9)
)
exposure <- c(0, 0, 1, 1)
wt_ate(ps_df, exposure)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[4]>
#> [1] 1.111111 1.428571 1.428571 1.111111
wt_ate(ps_df, exposure, .propensity_col = "treated")
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[4]>
#> [1] 1.111111 1.428571 1.428571 1.111111

# -- Censoring weights -----------------------------------------------
cens_ps <- runif(50, 0.6, 0.95)
cens_ind <- rbinom(50, 1, cens_ps)
wt_cens(cens_ps, cens_ind)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = uncensored}[50]>
#>  [1] 1.083078 1.265968 9.835512 1.243399 1.199165 4.524070 1.153452 1.091492
#>  [9] 1.127185 1.599997 1.662887 1.617418 1.107057 1.247371 3.446433 1.068590
#> [17] 1.239284 1.272276 3.764746 1.404803 1.132586 1.626884 3.667167 1.111347
#> [25] 1.184439 1.191086 3.594604 1.259226 1.543126 1.143215 1.098416 3.043524
#> [33] 7.675660 1.236511 1.069308 2.903103 1.275362 1.104237 1.105709 1.643950
#> [41] 1.062066 1.296044 1.358309 1.340223 1.580968 2.912338 1.347898 1.389593
#> [49] 1.341669 1.415592
estimand(wt_cens(cens_ps, cens_ind))  # "uncensored"
#> ℹ Treating `.exposure` as binary
#> [1] "uncensored"

# -- Categorical exposure from a multinomial fit ---------------------
set.seed(5)
x3 <- rnorm(100)
dose_level <- factor(sample(c("low", "mid", "high"), 100, replace = TRUE))
cat_model <- nnet::multinom(dose_level ~ x3, trace = FALSE)

# The fit reports one probability per level, and the exposure and its level
# order are read from the model
wt_ate(cat_model)
#> ℹ Using exposure variable "dose_level" from the propensity score model
#> ℹ Treating `.exposure` as categorical
#> <psw{estimand = ate}[100]>
#>   [1] 2.627904 2.829846 3.375768 2.703487 3.632186 3.233032 2.657388 2.644156
#>   [9] 3.113899 2.709502 3.539632 3.249417 2.609520 3.081511 2.610172 2.685315
#>  [17] 2.647216 3.656999 2.984049 2.675082 2.832620 2.823523 2.711684 3.445671
#>  [25] 2.850682 2.672211 3.575537 3.590855 3.224908 3.195973 2.725495 2.787071
#>  [33] 2.564519 2.812539 3.587106 3.489123 3.173263 3.038661 3.069665 3.304136
#>  [41] 2.694864 2.630910 3.060826 2.885123 2.658668 2.748277 2.624315 2.658366
#>  [49] 3.214893 2.691322 2.712655 3.357535 2.792821 2.647658 2.687614 2.621382
#>  [57] 2.865372 3.070283 3.316645 2.718015 3.155149 2.776314 3.189703 3.409303
#>  [65] 3.209603 2.668202 3.630985 3.117957 3.380608 3.112345 2.679760 3.291034
#>  [73] 2.728332 2.700167 2.734436 3.302117 3.493073 2.707990 2.996421 3.185624
#>  [81] 3.409654 3.072290 3.493442 3.048336 3.440244 3.224773 2.532538 2.657285
#>  [89] 3.061122 3.174848 2.783302 3.165692 2.909956 2.780810 2.813005 2.731683
#>  [97] 3.129433 3.179659 3.297815 3.317412

# The same weights from the fitted matrix, which needs the exposure and its
# type supplied
ps_mat <- predict(cat_model, type = "probs")
wt_ate(ps_mat, dose_level, exposure_type = "categorical")
#> <psw{estimand = ate}[100]>
#>   [1] 2.627904 2.829846 3.375768 2.703487 3.632186 3.233032 2.657388 2.644156
#>   [9] 3.113899 2.709502 3.539632 3.249417 2.609520 3.081511 2.610172 2.685315
#>  [17] 2.647216 3.656999 2.984049 2.675082 2.832620 2.823523 2.711684 3.445671
#>  [25] 2.850682 2.672211 3.575537 3.590855 3.224908 3.195973 2.725495 2.787071
#>  [33] 2.564519 2.812539 3.587106 3.489123 3.173263 3.038661 3.069665 3.304136
#>  [41] 2.694864 2.630910 3.060826 2.885123 2.658668 2.748277 2.624315 2.658366
#>  [49] 3.214893 2.691322 2.712655 3.357535 2.792821 2.647658 2.687614 2.621382
#>  [57] 2.865372 3.070283 3.316645 2.718015 3.155149 2.776314 3.189703 3.409303
#>  [65] 3.209603 2.668202 3.630985 3.117957 3.380608 3.112345 2.679760 3.291034
#>  [73] 2.728332 2.700167 2.734436 3.302117 3.493073 2.707990 2.996421 3.185624
#>  [81] 3.409654 3.072290 3.493442 3.048336 3.440244 3.224773 2.532538 2.657285
#>  [89] 3.061122 3.174848 2.783302 3.165692 2.909956 2.780810 2.813005 2.731683
#>  [97] 3.129433 3.179659 3.297815 3.317412

# A multinomial model of the exposure on a modifier stabilizes the weights on
# the probability of the level each unit took given that modifier
v <- rbinom(100, 1, 0.5)
cat_num_model <- nnet::multinom(dose_level ~ v, trace = FALSE)
w_cat <- wt_ate(cat_model, stabilize = cat_num_model)
#> ℹ Using exposure variable "dose_level" from the propensity score model
#> ℹ Treating `.exposure` as categorical
numerator_model(w_cat)
#> Call:
#> nnet::multinom(formula = dose_level ~ v, trace = FALSE)
#> 
#> Coefficients:
#>     (Intercept)          v
#> low   0.3022743 -0.4357960
#> mid  -0.6359876  0.9079251
#> 
#> Residual Deviance: 211.8462 
#> AIC: 219.8462 
```
