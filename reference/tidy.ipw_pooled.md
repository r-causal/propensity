# Tidy a pooled inverse probability weighted result

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
estimates of a
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
result as a tibble using the column names broom conventions use. There
is one row per effect measure, and one row per effect measure per
contrast for a categorical exposure, in the order the pooled result
stores them. Nothing is dropped.

Those columns are the ones the pooled result's own coercion surface
reports, so this method is that surface read as a tibble rather than a
second assembly of the same table. The covariance of the effects the
frame attaches is the one thing that does not travel: it is an attribute
rather than a column, and a tidied table is its columns.

The pooling has already happened by the time this method sees the
result. Nothing here re-pools and no argument changes an estimate: they
say how the table already in the object is reported, which is whether an
interval comes with it, the level that interval is reported at, and the
scale the ratio rows are put on.

## Usage

``` r
# S3 method for class 'ipw_pooled'
tidy(x, conf.int = FALSE, conf.level = NULL, exponentiate = FALSE, ...)
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

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per pooled estimate and the columns:

- `term`:

  The effect measure, such as `"rd"`, `"log(rr)"`, `"log(or)"`,
  `"diff"`, or `"slope"`, or the coefficient name in the conditional
  reading.

- `contrast`:

  The contrast the row reports, such as `"b vs a"`. Categorical
  exposures only.

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
for a one-row summary of the pooled fit, and the Multiple imputation
section of
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for the workflow the two belong to.

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
#> # A tibble: 3 × 6
#>   term    estimate std.error statistic    df p.value
#>   <chr>      <dbl>     <dbl>     <dbl> <dbl>   <dbl>
#> 1 rd         0.152    0.0810      1.88  134.  0.0625
#> 2 log(rr)    0.315    0.174       1.80  134.  0.0733
#> 3 log(or)    0.614    0.332       1.85  134.  0.0668

tidy(pool_ipw(fits), conf.int = TRUE, exponentiate = TRUE)
#> # A tibble: 3 × 8
#>   term  estimate std.error statistic    df p.value conf.low conf.high
#>   <chr>    <dbl>     <dbl>     <dbl> <dbl>   <dbl>    <dbl>     <dbl>
#> 1 rd       0.152    0.0810      1.88  134.  0.0625 -0.00803     0.312
#> 2 rr       1.37     0.174       1.80  134.  0.0733  0.970       1.93 
#> 3 or       1.85     0.332       1.85  134.  0.0668  0.958       3.56 
```
