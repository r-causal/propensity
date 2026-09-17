# Tabulate Truncated Weights Over a Grid of Bounds

`wt_trunc_sensitivity()` truncates a set of weights with
[`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md)
at each bound in a grid and tabulates what each truncation does to them:
the bound as given and as applied, the number of weights moved, and the
range and mean of the weights that result. The first row describes the
weights before any truncation, so every other row can be read against
it. This is the weight side of the sensitivity table Cole and Hernán
(2008) recommend for choosing a bound.

## Usage

``` r
wt_trunc_sensitivity(
  .weights,
  method = c("pctl", "wt", "count"),
  lower = NULL,
  upper = NULL
)
```

## Arguments

- .weights:

  A [psw](https://r-causal.github.io/propensity/reference/psw.md)
  vector, such as one returned by
  [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
  or a numeric vector of weights. Propensity scores, including those
  returned by
  [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
  [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
  and
  [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md),
  are refused with an error of class `propensity_type_error`, and
  weights that have already been truncated with an error of class
  `propensity_already_modified_error`.

- method:

  The rule that reads each bound, one of `"pctl"` (the default), `"wt"`,
  or `"count"`, as in
  [`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md).
  The `"adaptive"` method sets a single bound from the number of
  weights, so it has no grid to vary and is not accepted.

- lower:

  An optional grid of lower bounds, read according to `method` as
  [`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md)
  reads `lower`. A single value is used with every upper bound;
  otherwise `lower` must have the same length as `upper`, and the two
  are paired element by element.

- upper:

  The grid of upper bounds, read according to `method` as
  [`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md)
  reads `upper`. For `"pctl"` it defaults to
  `c(0.99, 0.975, 0.95, 0.90)`; for `"wt"` and `"count"` it is required.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row for the untruncated weights followed by one row per bound, and the
columns:

- `lower`, `upper`: the bounds as given, `NA` in the first row and
  `lower` `NA` when no lower bound was given.

- `lower_value`, `upper_value`: the bounds as applied, on the scale of
  the weights, recorded by
  [`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md).
  Both are `NA` in the first row, and `lower_value` is `NA` when no
  lower bound was given.

- `n_truncated`: the number of weights moved, 0 in the first row.

- `min`, `max`, `mean`: the smallest, largest, and mean weight present.

- `range_ratio`: `max / min`.

## Details

Each row after the first is computed by calling
[`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md)
with that row's bounds, so it describes exactly the weights that call
returns. Every bound in the grid must be one
[`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md)
accepts: a single invalid bound refuses the whole grid, with the error
[`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md)
would raise for it, rather than leaving a gap in the table. The rows
follow the grid in the order given.

The summaries are taken over the weights that are present, so missing
weights, such as those of units set aside by
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
take no part in them and are never counted as truncated. A grid over
weights of which none are present is refused with an error of class
`propensity_range_error`. Weights can be zero, and a zero weight makes
`range_ratio` infinite.

Weights that have already been truncated are refused, because the first
row would then describe weights that were already bounded rather than
the weights before truncation. Pass the weights as they were before
[`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md)
instead.

### What is left to the caller

The table does not include effect estimates or effective sample sizes.
propensity fits no outcome models, so to add estimates to the table, map
over the grid, truncate the weights at each bound with
[`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md),
and fit the outcome model with each set of weights. For the effective
sample size, add the weights before and after truncation to a data frame
as separate columns and pass them to `halfmoon::check_ess()`.

## References

Cole, S. R., & Hernán, M. A. (2008). Constructing inverse probability
weights for marginal structural models. *American Journal of
Epidemiology*, 168(6), 656–664.

## See also

[`wt_trunc()`](https://r-causal.github.io/propensity/reference/wt_trunc.md)
to truncate the weights at a chosen bound.

## Examples

``` r
set.seed(1)
n <- 200
x <- rnorm(n)
z <- rbinom(n, 1, plogis(2 * x))
fit <- glm(z ~ x, family = binomial)
w <- wt_ate(fit, .exposure = z)
#> ℹ Treating `.exposure` as binary

# The default percentile grid. The first row is the untruncated weights.
wt_trunc_sensitivity(w)
#> # A tibble: 5 × 9
#>   lower  upper lower_value upper_value n_truncated   min   max range_ratio  mean
#>   <dbl>  <dbl>       <dbl>       <dbl>       <int> <dbl> <dbl>       <dbl> <dbl>
#> 1    NA NA              NA       NA              0  1.02 31.6        31.0   2.05
#> 2    NA  0.99           NA        7.94           2  1.02  7.94        7.79  1.91
#> 3    NA  0.975          NA        6.44           5  1.02  6.44        6.32  1.87
#> 4    NA  0.95           NA        4.31          10  1.02  4.31        4.23  1.79
#> 5    NA  0.9            NA        3.33          20  1.02  3.33        3.27  1.72

# A grid of absolute bounds, and a lower tail bounded at every one of them
wt_trunc_sensitivity(w, method = "wt", upper = c(20, 10, 5))
#> # A tibble: 4 × 9
#>   lower upper lower_value upper_value n_truncated   min   max range_ratio  mean
#>   <dbl> <dbl>       <dbl>       <dbl>       <int> <dbl> <dbl>       <dbl> <dbl>
#> 1    NA    NA          NA          NA           0  1.02  31.6       31.0   2.05
#> 2    NA    20          NA          20           1  1.02  20         19.6   1.99
#> 3    NA    10          NA          10           2  1.02  10          9.81  1.93
#> 4    NA     5          NA           5           8  1.02   5          4.91  1.82
wt_trunc_sensitivity(w, method = "wt", lower = 1.05, upper = c(20, 10, 5))
#> # A tibble: 4 × 9
#>   lower upper lower_value upper_value n_truncated   min   max range_ratio  mean
#>   <dbl> <dbl>       <dbl>       <dbl>       <int> <dbl> <dbl>       <dbl> <dbl>
#> 1 NA       NA       NA             NA           0  1.02  31.6       31.0   2.05
#> 2  1.05    20        1.05          20          12  1.05  20         19.0   1.99
#> 3  1.05    10        1.05          10          13  1.05  10          9.52  1.93
#> 4  1.05     5        1.05           5          19  1.05   5          4.76  1.83

# The largest one, three, and ten weights
wt_trunc_sensitivity(w, method = "count", upper = c(1, 3, 10))
#> # A tibble: 4 × 9
#>   lower upper lower_value upper_value n_truncated   min   max range_ratio  mean
#>   <dbl> <dbl>       <dbl>       <dbl>       <int> <dbl> <dbl>       <dbl> <dbl>
#> 1    NA    NA          NA       NA              0  1.02 31.6        31.0   2.05
#> 2    NA     1          NA       12.5            1  1.02 12.5        12.3   1.95
#> 3    NA     3          NA        7.58           3  1.02  7.58        7.43  1.90
#> 4    NA    10          NA        4.29          10  1.02  4.29        4.21  1.79
```
