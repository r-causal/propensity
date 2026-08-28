# Tidy a pooled inverse probability weighted result

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
estimates of a
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
result as a tibble using the column names broom conventions use. It
holds the rows the results themselves report: for a binary or
categorical exposure one row per exposure level, under the term
`"mean"`, and then one row per effect measure per contrast of levels;
for a continuous exposure one row per reported coefficient. The rows
arrive in the order the pooled result stores them and nothing is
dropped.

Those columns are the ones the pooled result's own coercion surface
reports, so this method is that surface read as a tibble rather than a
second assembly of the same table. The covariance of the effects the
frame attaches is the one thing that does not travel: it is an attribute
rather than a column, and a tidied table is its columns.

The pooling has already happened by the time this method sees the
result. Nothing here re-pools and no argument changes an estimate: they
say how the table already in the object is reported, which is which of
the two readings it is, whether an interval comes with it, the level
that interval is reported at, and the scale the ratio rows are put on.

A pooled result reports its effects in one of two readings, and
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
one the result records unless `effects` names the other for the call.
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
pools both readings of the results it is given, so naming one here
reports a table the pooling already built rather than pooling again. The
marginal reading is the table of pooled causal contrasts described
above. The conditional reading is the outcome models' coefficient
surface, pooled the same way: one row per coefficient, from the block of
the joint estimation that carries the uncertainty of having estimated
the weights from the same data. The two readings return the same columns
in the same order, with one exception: the `contrast` column belongs to
the marginal reading. A binary and a categorical pool both carry it
there, naming the exposure level each mean row belongs to and the pair
of levels each effect measure compares, and the conditional reading of
either returns the same table one column narrower.

## Usage

``` r
# S3 method for class 'ipw_pooled'
tidy(
  x,
  conf.int = FALSE,
  conf.level = NULL,
  exponentiate = FALSE,
  ...,
  effects = NULL
)
```

## Arguments

- x:

  An `ipw_pooled` object, as returned by
  [`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html).

- conf.int:

  Logical. Should the confidence interval bounds be returned in the
  `conf.low` and `conf.high` columns? Defaults to `FALSE`.

- conf.level:

  The confidence level of the interval. `NULL`, the default, reports the
  bounds at the level the pooled result stored, which is the level its
  own intervals were built at. A number between 0 and 1 rebuilds them at
  that level instead, from the estimate and its standard error referred
  to t on the pooled degrees of freedom. Not used when `conf.int` is
  `FALSE`, though it must still be a valid level.

- exponentiate:

  Logical. Should the estimate and its bounds be exponentiated on the
  rows reported on the log scale? Defaults to `FALSE`. This behaves
  exactly as it does on the pooled result's own coercion surface: in the
  marginal reading the `log(rr)` and `log(or)` rows have their estimate
  and bounds exponentiated and their `term` relabeled to `rr` and `or`,
  and in the conditional reading the outcome model's link settles the
  scale for every row at once and no term is relabeled. In both readings
  the standard error, the test statistic, and the p-value describe the
  log scale and stay there, which is where the inference was done.

- ...:

  These dots are for future extensions and must be empty.

- effects:

  The reading to report, either `"marginal"` or `"conditional"`. `NULL`,
  the default, reports the reading the pooled result records; naming a
  reading reports that one for the one call and leaves the result as it
  is. The marginal reading reports the pooled causal contrasts; the
  conditional reading reports the pooled coefficients of the outcome
  models, which is the reading `exponentiate` above puts on the natural
  scale every row at once.
  [`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  and
  [`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
  move a pooled result between the two readings for good rather than for
  a call.

  A reading the pooling could not compute is refused rather than
  reported under the other one's name. The conditional reading needs the
  covariance the joint estimation of the weights and the outcome
  implies, and a linearization fit stacks no such system and records no
  such block, so a set of those fits pools the marginal reading alone. A
  result pooled before both readings were kept holds only the one it was
  pooled in, and asking it for the other is refused the same way.

  `"fixed"` is not among the values this method takes, though
  [`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw.md)
  on a single result accepts it.
  [`mice::pool()`](https://amices.org/mice/reference/pool.html) asks for
  that reading from the tidier of each imputation's own fit, which is
  where the name arrives and where the tolerance for it belongs; nothing
  tidies a pooled result on its behalf.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per pooled estimate and the columns:

- `term`:

  The effect measure, such as `"mean"` for a counterfactual mean, or
  `"rd"`, `"log(rr)"`, `"log(or)"`, `"diff"`, `"slope"`, or `"coef"` for
  the measures built from them, or the coefficient name in the
  conditional reading.

- `contrast`:

  What the row is about: the exposure level a `"mean"` row belongs to,
  such as `"0"` or `"b"`; the pair of levels an effect measure compares,
  such as `"1 vs 0"` or `"b vs a"`; or the coefficient name on a surface
  whose rows are named after their coefficients. Present only in the
  marginal reading, and there only where a row reports one: the
  conditional reading names each coefficient in `term` and returns no
  `contrast` column.

- `estimate`:

  The pooled point estimate.

- `std.error`:

  The pooled standard error.

- `statistic`:

  The test statistic, the estimate over its standard error.

- `df`:

  The pooled degrees of freedom the statistic is referred to. These
  carry the Barnard-Rubin small-sample adjustment, so the statistic is
  referred to t on them rather than to the normal.

- `p.value`:

  The two-sided p-value of that statistic.

- `conf.low`, `conf.high`:

  The interval bounds. Present only when `conf.int` is `TRUE`.

## See also

[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
for the pooling,
[`glance()`](https://r-causal.github.io/propensity/reference/glance.ipw_pooled.md)
for a one-row summary of the pooled fit,
[`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
and
[`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.html)
for the reading a pooled result records, and the Multiple imputation
section of
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for the workflow they belong to.

## Examples

``` r
set.seed(2024)
n <- 150
x1 <- rnorm(n)
z <- rbinom(n, 1, plogis(0.3 * x1))
y <- rbinom(n, 1, plogis(-0.4 + 0.9 * z + 0.5 * x1))
dat <- data.frame(x1, z, y)
dat$x1[rbinom(n, 1, 0.2) == 1] <- NA

imp <- mice::mice(dat, m = 2, print = FALSE, seed = 1)
fits <- with(imp, {
  ps <- glm(z ~ x1, family = binomial())
  w <- wt_ate(ps)
  om <- glm(y ~ z, family = quasibinomial(), weights = w)
  ipw(ps, om)
})
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary

tidy(pool_ipw(fits))
#> # A tibble: 5 × 7
#>   term    contrast estimate std.error statistic    df  p.value
#>   <chr>   <chr>       <dbl>     <dbl>     <dbl> <dbl>    <dbl>
#> 1 mean    0           0.412    0.0595      6.91  136. 1.68e-10
#> 2 mean    1           0.564    0.0551     10.2   138. 1.29e-18
#> 3 rd      1 vs 0      0.152    0.0810      1.88  134. 6.25e- 2
#> 4 log(rr) 1 vs 0      0.315    0.174       1.80  134. 7.33e- 2
#> 5 log(or) 1 vs 0      0.614    0.332       1.85  134. 6.68e- 2

tidy(pool_ipw(fits), conf.int = TRUE, exponentiate = TRUE)
#> # A tibble: 5 × 9
#>   term  contrast estimate std.error statistic    df  p.value conf.low conf.high
#>   <chr> <chr>       <dbl>     <dbl>     <dbl> <dbl>    <dbl>    <dbl>     <dbl>
#> 1 mean  0           0.412    0.0595      6.91  136. 1.68e-10  0.294       0.529
#> 2 mean  1           0.564    0.0551     10.2   138. 1.29e-18  0.455       0.673
#> 3 rd    1 vs 0      0.152    0.0810      1.88  134. 6.25e- 2 -0.00803     0.312
#> 4 rr    1 vs 0      1.37     0.174       1.80  134. 7.33e- 2  0.970       1.93 
#> 5 or    1 vs 0      1.85     0.332       1.85  134. 6.68e- 2  0.958       3.56 

# The outcome models' coefficients, pooled the same way
tidy(pool_ipw(fits), effects = "conditional")
#> # A tibble: 2 × 6
#>   term        estimate std.error statistic    df p.value
#>   <chr>          <dbl>     <dbl>     <dbl> <dbl>   <dbl>
#> 1 (Intercept)   -0.358     0.246     -1.46  136.  0.148 
#> 2 z              0.614     0.332      1.85  134.  0.0668
```
