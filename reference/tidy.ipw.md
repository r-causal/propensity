# Tidy an inverse probability weighted result

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
estimates of an
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
result as a tibble using the column names broom conventions use. For a
binary or categorical exposure there is one row per exposure level,
under the term `"mean"`, and then one row per effect measure per
contrast of levels; for a continuous exposure there is one row per
reported coefficient. The rows arrive in the order the result stores
them and nothing is dropped.

Those columns are the ones the result's own coercion surface reports, so
either reading is
[`as.data.frame()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.html)
read as a tibble rather than a second assembly of the same table. The
values are the ones the result already holds, nothing being
re-estimated: the confidence interval is the only thing rebuilt, and
only when the requested `conf.level` differs from the level the result
was fit at. The one thing the frame carries that the tibble does not is
the covariance it attaches as an attribute, a tidied table being its
columns.

`conf.int`, `conf.level`, and `exponentiate` are the arguments this
method shares with that surface, so the surface validates them, in both
readings and before either is assembled. A value it refuses is reported
under the causalgenerics condition naming the argument at fault, which
is what keeps one argument of one method from being well formed in one
reading of a result and refused in the other.

A result reports its effects in one of two readings, and
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
one the result records unless `effects` names the other for the call.
The marginal reading is the table of causal contrasts described above.
The conditional reading is the outcome model's coefficient surface: one
row per coefficient, with the standard errors of the block of the joint
estimation that carries the uncertainty of having estimated the weights
from the same data. The two readings return the same columns in the same
order, with one exception: the `contrast` column belongs to the marginal
reading. A binary and a categorical result both carry it there, naming
the exposure level each mean row belongs to and the pair of levels each
effect measure compares, and the conditional reading of either returns
the same table one column narrower. A coefficient names neither a level
nor a pair of them.

## Usage

``` r
# S3 method for class 'ipw'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  ...,
  effects = NULL,
  parametric = NULL
)
```

## Arguments

- x:

  An `ipw` object, as returned by
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).

- conf.int:

  Logical. Should the confidence interval bounds be returned in the
  `conf.low` and `conf.high` columns? Defaults to `FALSE`.

- conf.level:

  The confidence level of the interval, a single number between 0 and 1.
  Defaults to `0.95`, which is this method's own default rather than the
  level `x` was fit at. When the two differ, the bounds are recomputed
  as the normal approximation
  `estimate +/- qnorm(1 - (1 - conf.level) / 2) * std.error`, which is
  the interval
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  itself reports. Not used when `conf.int` is `FALSE`, though it must
  still be a valid level.

- exponentiate:

  Logical. Should the estimate and its bounds be exponentiated on the
  rows reported on the log scale? Defaults to `FALSE`. This behaves
  exactly as it does in
  [`as.data.frame()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.html):
  the `log(rr)` and `log(or)` rows have their estimate and bounds
  exponentiated and their `term` relabeled to `rr` and `or`, while every
  other row is left alone. The standard error, the test statistic, and
  the p-value describe the log scale estimate and stay there.

  The conditional reading has no rows labeled as ratios to pick out, so
  the link of the outcome model settles the question there: a `logit`
  link puts every coefficient on the log odds scale and a `log` link
  puts every coefficient on the log risk scale, and both are scales an
  exponential undoes. The estimate and, when they were asked for, the
  interval bounds are exponentiated, no term is relabeled, and the
  standard error, the test statistic, and the p-value describe the link
  scale and stay there. Every other link errors rather than
  exponentiating coefficients that are not on such a scale.

- ...:

  These dots are for future extensions and must be empty.

- effects:

  The reading to report, either `"marginal"` or `"conditional"`. `NULL`,
  the default, reports the reading the result records; any other value
  overrides it for the one call and leaves the result as it is. The
  marginal reading reports the population-averaged causal contrasts; the
  conditional reading reports the outcome model's coefficient surface.
  [`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  and
  [`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  move a result between the two readings.

  `"fixed"` is also accepted, and asks for the reading the result
  records, exactly as `NULL` does. It is the name mixed models give the
  effects that are not random, and
  [`mice::pool()`](https://amices.org/mice/reference/pool.html) asks
  every tidier it calls for it. A result reports one reading at a time
  and neither of them is random, so the name resolves to whichever one
  the result records. The accessors do not accept it: it is spelled
  here, where multiple imputation reaches the result, rather than in the
  surface underneath.

  The covariance the conditional reading reports is the outcome block of
  the jointly estimated sandwich, which every route that stacks
  estimating equations attaches to the outcome model it stores:
  `se_method = "mestimation"` for a binary exposure, the categorical and
  joint routes, which run on M-estimation alone, and the continuous
  route under that method. A linearization or robust fit stacks no
  system and records no covariance of either kind, so its conditional
  reading errors rather than reporting the covariance the outcome model
  computed for itself, which treats the estimated weights as fixed.

- parametric:

  Accepted and ignored.
  [`mice::pool()`](https://amices.org/mice/reference/pool.html) passes
  `parametric = TRUE` to every tidier it calls, to ask the models it was
  written against for their parametric coefficient table rather than for
  a smooth term. An `ipw` result reports one table, and the inference it
  reports is already parametric, so there is nothing for the argument to
  select and nothing is read from it. It is here so that a result can be
  pooled across multiply imputed datasets.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per estimate and the columns:

- `term`:

  The effect measure, such as `"mean"` for a counterfactual mean, or
  `"rd"`, `"log(rr)"`, `"log(or)"`, `"diff"`, `"slope"`, or `"coef"` for
  the measures built from them.

- `contrast`:

  What the row is about: the exposure level a `"mean"` row belongs to,
  such as `"0"` or `"b"`; the pair of levels an effect measure compares,
  such as `"1 vs 0"` or `"b vs a"`; or the coefficient name on a surface
  whose rows are named after their coefficients. Present only where a
  row reports one.

- `estimate`:

  The estimated effect.

- `std.error`:

  The standard error of the estimate.

- `statistic`:

  The z statistic, the estimate over its standard error.

- `p.value`:

  The two-sided p-value of that statistic.

- `conf.low`, `conf.high`:

  The interval bounds. Present only when `conf.int` is `TRUE`.

The conditional reading returns those columns in that order, holding one
row per coefficient of the outcome model: `term` is the coefficient's
name, `estimate` is the coefficient, `std.error` is the square root of
the diagonal of the corrected covariance, `statistic` is the estimate
over that standard error, and `p.value` is the two-sided normal p-value
of the statistic. There is no `contrast` column, whatever the exposure,
because a coefficient names neither an exposure level nor a pair of
them, and the bounds are built at the level `conf.level` asks for, the
stored ones belonging to the effects the marginal reading reports.

## See also

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for the estimator,
[`glance()`](https://r-causal.github.io/propensity/reference/glance.ipw.md)
for a one-row summary of the fit,
[`augment()`](https://r-causal.github.io/propensity/reference/augment.ipw.md)
for its per-observation columns,
[`as.data.frame()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.html)
for the result's own columns, and
[`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
and
[`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
for the reading a result records. For results fitted to multiply imputed
data, see the Multiple imputation section of
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html),
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
for the pooling, and
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw_pooled.md)
for what it returns.

## Examples

``` r
set.seed(123)
n <- 200
x1 <- rnorm(n)
z <- rbinom(n, 1, plogis(0.5 * x1))
y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
dat <- data.frame(x1, z, y)

ps_mod <- glm(z ~ x1, data = dat, family = binomial())
wts <- wt_ate(ps_mod)
#> ℹ Using exposure variable "z" from the propensity score model
#> ℹ Treating `.exposure` as binary
outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
result <- ipw(ps_mod, outcome_mod)

tidy(result)
#> # A tibble: 5 × 6
#>   term    contrast estimate std.error statistic  p.value
#>   <chr>   <chr>       <dbl>     <dbl>     <dbl>    <dbl>
#> 1 mean    0           0.440    0.0506      8.69 3.69e-18
#> 2 mean    1           0.582    0.0487     11.9  7.14e-33
#> 3 rd      1 vs 0      0.142    0.0702      2.03 4.27e- 2
#> 4 log(rr) 1 vs 0      0.280    0.142       1.97 4.87e- 2
#> 5 log(or) 1 vs 0      0.573    0.287       2.00 4.55e- 2

# A 90% interval, with the ratio measures on their natural scale
tidy(result, conf.int = TRUE, conf.level = 0.9, exponentiate = TRUE)
#> # A tibble: 5 × 8
#>   term  contrast estimate std.error statistic  p.value conf.low conf.high
#>   <chr> <chr>       <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 mean  0           0.440    0.0506      8.69 3.69e-18   0.357      0.523
#> 2 mean  1           0.582    0.0487     11.9  7.14e-33   0.502      0.662
#> 3 rd    1 vs 0      0.142    0.0702      2.03 4.27e- 2   0.0268     0.258
#> 4 rr    1 vs 0      1.32     0.142       1.97 4.87e- 2   1.05       1.67 
#> 5 or    1 vs 0      1.77     0.287       2.00 4.55e- 2   1.11       2.84 

# The outcome model's coefficients, with the standard errors the joint
# estimation of the weights and the outcome implies
tidy(result, conf.int = TRUE, effects = "conditional")
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 (Intercept)   -0.242     0.205     -1.18  0.239   -0.645      0.161
#> 2 z              0.573     0.287      2.00  0.0455   0.0115     1.14 
```
