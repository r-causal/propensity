# Tidy an inverse probability weighted result

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns the
estimates of an
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
result as a tibble using the column names broom conventions use. There
is one row per effect measure, and one row per effect measure per
comparison for a categorical exposure, in the order the result stores
them. Nothing is dropped.

The values are the ones the result already holds:
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) renames and
selects rather than re-estimating anything. The one exception is the
confidence interval, which is rebuilt from the estimate and its standard
error when the requested `conf.level` differs from the level the result
was fit at.

## Usage

``` r
# S3 method for class 'ipw'
tidy(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...)
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
  itself reports. Ignored when `conf.int` is `FALSE`.

- exponentiate:

  Logical. Should the estimate and its bounds be exponentiated on the
  rows reported on the log scale? Defaults to `FALSE`. This behaves
  exactly as it does in
  [`as.data.frame()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.html):
  the `log(rr)` and `log(or)` rows have their estimate and bounds
  exponentiated and their `term` relabeled to `rr` and `or`, while every
  other row is left alone. The standard error, the test statistic, and
  the p-value describe the log scale estimate and stay there.

- ...:

  These dots are for future extensions and must be empty.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per estimate and the columns:

- `term`:

  The effect measure, such as `"rd"`, `"log(rr)"`, `"log(or)"`,
  `"diff"`, or `"slope"`.

- `comparison`:

  The contrast the row reports, such as `"b vs a"`. Categorical
  exposures only.

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

## See also

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for the estimator,
[`glance()`](https://r-causal.github.io/propensity/reference/glance.ipw.md)
for a one-row summary of the fit,
[`augment()`](https://r-causal.github.io/propensity/reference/augment.ipw.md)
for its per-observation columns, and
[`as.data.frame()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.html)
for the result's own columns.

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
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
result <- ipw(ps_mod, outcome_mod)

tidy(result)
#> # A tibble: 3 × 5
#>   term    estimate std.error statistic p.value
#>   <chr>      <dbl>     <dbl>     <dbl>   <dbl>
#> 1 rd         0.142    0.0702      2.03  0.0427
#> 2 log(rr)    0.280    0.142       1.97  0.0487
#> 3 log(or)    0.573    0.287       2.00  0.0455

# A 90% interval, with the ratio measures on their natural scale
tidy(result, conf.int = TRUE, conf.level = 0.9, exponentiate = TRUE)
#> # A tibble: 3 × 7
#>   term  estimate std.error statistic p.value conf.low conf.high
#>   <chr>    <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 rd       0.142    0.0702      2.03  0.0427   0.0268     0.258
#> 2 rr       1.32     0.142       1.97  0.0487   1.05       1.67 
#> 3 or       1.77     0.287       2.00  0.0455   1.11       2.84 
```
