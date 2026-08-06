# Propensity score tilting functions

Every estimand this package supports targets a population whose
covariate distribution is the study population reweighted by a tilting
function \\h\\ of the propensity score. `ps_tilt()` evaluates \\h\\ at
each observation's propensity score.

The tilt is the numerator of every weight: a weight is \\h\\ divided by
the propensity score of the exposure level the unit actually received.
It is also what standardizes a g-computation estimate to a target
population, which for `"atm"`, `"ato"`, and `"entropy"` is the only
route to the estimand, since those populations are not subsets of the
rows and cannot be reached by filtering.

## Usage

``` r
ps_tilt(
  .propensity,
  estimand,
  ...,
  .focal_level = NULL,
  ps = lifecycle::deprecated()
)
```

## Arguments

- .propensity:

  Propensity scores. A numeric vector of \\P(Z = \text{focal} \mid X)\\
  for a binary exposure, or a matrix or data frame with one column per
  exposure level, named for that level, for a categorical exposure.

- estimand:

  One of `"ate"`, `"att"`, `"atu"`, `"atm"`, `"ato"`, or `"entropy"`.

- ...:

  These dots are for future extensions and must be empty.

- .focal_level:

  The exposure level the `"att"` and `"atu"` tilts target, matched
  against the column names of `.propensity`. A column named for the
  level exactly wins; failing that, a `.pred_` prefix is stripped, so
  `.pred_a` matches the level `a`. Required for those two estimands with
  a categorical `.propensity`, and accepted nowhere else: a numeric
  `.propensity` is already the probability of the focal level, and the
  remaining tilts treat every level alike.

- ps:

  **\[deprecated\]** Use `.propensity` instead. A call that names `ps`
  must name the arguments after it as well, since a positional argument
  binds to `.propensity`.

## Value

A plain double vector, unnamed, with one element per observation: the
length of `.propensity` for the numeric method, and the number of rows
of `.propensity` for the matrix and data frame methods.

## Tilting functions

For a binary exposure with propensity score \\e = P(Z = \text{focal}
\mid X)\\:

|  |  |  |
|----|----|----|
| estimand | \\h(e)\\ | target population |
| `"ate"` | \\1\\ | everyone |
| `"att"` | \\e\\ | the focal group |
| `"atu"` | \\1 - e\\ | the reference group |
| `"atm"` | \\\min(e, 1 - e)\\ | the evenly matchable |
| `"ato"` | \\e(1 - e)\\ | the overlap population |
| `"entropy"` | \\-e \log(e) - (1 - e) \log(1 - e)\\ | the entropy-tilted population |

For a categorical exposure with propensity score vector \\(e_1, \ldots,
e_K)\\ and focal level \\f\\: `"ate"` is \\1\\, `"att"` is \\e_f\\,
`"atu"` is \\1 - e_f\\, `"atm"` is \\\min_k e_k\\, `"ato"` is \\(\sum_k
1 / e_k)^{-1}\\, and `"entropy"` is \\-\sum_k e_k \log(e_k)\\.

Stabilization and censoring weights are not tilts. Stabilization
multiplies a weight by a marginal quantity that does not depend on the
covariates, and
[`wt_cens()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
reuses the `"ate"` formula, so neither has an entry here.

## Standardizing model predictions

A g-computation estimate averages an outcome model's per-row
counterfactual predictions, and where those predictions vary from row to
row it is the weights of that average that standardize the estimate to a
target population. The tilt is those weights, which is what the second
example below supplies by hand. An average taken with equal weight
standardizes a covariate-adjusted outcome model to the whole sample, so
such a model targets the ATE however its own fitting weights were built:
fit for a non-`"ate"` estimand and averaged that way, it reports the
full-sample contrast rather than the one it was weighted for. Tools that
average per-row model predictions, such as the `avg_comparisons()`
function in the marginaleffects package, average with equal weight
unless they are told otherwise, and the remedy is to hand them the tilt
as the weights of the average: `ps_tilt(ps, "att")` for a binary
exposure and `ps_tilt(ps, "att", .focal_level = "b")` for a categorical
one.

The requirement is easy to miss because an outcome model saturated in
the exposure hides it. Such a model predicts one value per exposure
level, so every average of its predictions returns the same contrast,
and there the estimand is settled by the weights the model was fit with
rather than by the weights of the average: a model fit with
[`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
weights reports the ATT whether the average is taken with equal weight
or with the tilt. The two averages come apart at the first covariate the
outcome model adjusts for, where the predictions vary from row to row
and the weights the average is taken under are what decides the
population the estimate describes.

## Propensity score range

Every propensity score in `.propensity` must lie strictly inside \\(0,
1)\\, the same requirement
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and the rest of the weight family impose before they tilt. The bound
holds for a matrix or data frame `.propensity` entry by entry, and each
row must sum to one on top of it. A fitted model that separates the
exposure can return a probability of exactly zero or one; those scores
have no weight to divide and are rejected here rather than tilted.

A missing propensity score gives a missing tilt under every estimand,
`"ate"` included, so an observation whose propensity score is unknown
never counts toward a tilted mean. A numeric `.propensity` propagates
`NA` position by position, and a matrix `.propensity` gives `NA` for any
row holding one: a probability vector with a missing entry is not one
this can tilt on, whichever level the tilt reads.

## Modified propensity scores

`ps_tilt()` takes plain propensity scores. A score modified by
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
or
[`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
carries a class of its own and has no method here; pass the scores
underneath it, with `as.numeric(x)` for a binary exposure or
`as.matrix(x)` for a categorical one. Units
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
set to `NA` stay `NA` through that extraction and take an `NA` tilt.

## See also

[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and the rest of the weight family, which divide the tilt by the
propensity score of the received exposure level.

## Examples

``` r
set.seed(1)
n <- 500
x <- rnorm(n)
z <- rbinom(n, 1, plogis(0.6 * x))
y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + x))
sim <- data.frame(x = x, z = z, y = y)
ps <- unname(
  predict(glm(z ~ x, data = sim, family = binomial()), type = "response")
)

# a weight is the tilt over the propensity of the received exposure level
received <- z * ps + (1 - z) * (1 - ps)
all.equal(
  as.numeric(wt_ato(ps, z, exposure_type = "binary")),
  ps_tilt(ps, "ato") / received
)
#> [1] TRUE

# tilted g-computation: standardize counterfactual predictions to the
# overlap population with an h-weighted mean, no weights in the outcome model
fit <- glm(y ~ z * x, data = sim, family = binomial())
m1 <- predict(fit, transform(sim, z = 1), type = "response")
m0 <- predict(fit, transform(sim, z = 0), type = "response")

h <- ps_tilt(ps, "ato")
weighted.mean(m1, h) - weighted.mean(m0, h)
#> [1] 0.1231355

# the same predictions standardized to everyone give the ATE instead
mean(m1) - mean(m0)
#> [1] 0.1210093
```
