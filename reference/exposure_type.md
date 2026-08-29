# What a set of weights records about the exposure

Weights carry two records describing the exposure they were built for.
`exposure_type()` returns the type of exposure the weight function was
given, and `density_meta()` returns the record left by weights whose
value is a ratio of densities.

Both describe the exposure rather than the units, so neither holds a
length of its own: both survive subsetting, arithmetic, and anything
else that changes the number of weights. Combining two sets of weights
that record either of them differently drops the record they disagree
on, with a warning of class `propensity_metadata_conflict_warning`. A
result that no longer records an exposure type is not described as
continuous, so it drops any density record along with it.

## Usage

``` r
exposure_type(wt)

density_meta(wt)

# S3 method for class 'propensity_density_meta'
format(x, ...)

# S3 method for class 'propensity_density_meta'
print(x, ...)
```

## Arguments

- wt:

  A `psw` or `causal_wts` object.

- x:

  A density record, as returned by `density_meta()`.

- ...:

  Not used.

## Value

- `exposure_type()`: a single string, or `NULL`.

- `density_meta()`: a list of class `propensity_density_meta` with the
  elements `density`, `numerator`, `sigma`, and `sigma_value`, or
  `NULL`.

- [`print()`](https://rdrr.io/r/base/print.html) and
  [`format()`](https://rdrr.io/r/base/format.html): the record,
  invisibly, and a character vector of one line per element.

## Details

`exposure_type()` returns `"binary"`, `"categorical"`, or
`"continuous"`. It returns `NULL` for weights that were built without an
exposure to describe, such as those written with
[`psw()`](https://r-causal.github.io/propensity/reference/psw.md)
directly.

`density_meta()` returns a record of what the ratio was built from, and
`NULL` for weights that are not a ratio of densities, which is every set
of weights for a binary or categorical exposure:

- `density`, the specification of the conditional density family, as
  built by
  [`dens_normal()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
  and its relatives.

- `numerator`, what stabilized the weights: `"marginal"` for the
  marginal density of the exposure, `"integrated"` for the conditional
  density marginalized over the units, `"score"` for a
  `stabilization_score` the caller supplied, and `"none"` for weights
  that were not stabilized.

- `sigma`, where the residual spread of the conditional density came
  from: `"pooled"` for the pooled residual standard deviation, and
  `"supplied"` for a `.sigma` the caller gave.

- `sigma_value`, the single spread the caller supplied, and `NULL` for a
  pooled spread and for one supplied per observation. A spread that is
  one number is a constant the weights can be rebuilt from, which is
  what
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  needs of it; a spread that changes with the observation is not, so the
  record holds where it came from and nothing more.

## See also

[`psw()`](https://r-causal.github.io/propensity/reference/psw.md) for
the weight vector class and the rest of what it records, and
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for the functions that write these records.

## Examples

``` r
set.seed(1)
ps <- runif(20, 0.2, 0.8)
trt <- rbinom(20, 1, ps)

# A binary exposure is weighted by propensity scores rather than by a
# density, so there is no density to record.
exposure_type(wt_ate(ps, trt))
#> ℹ Treating `.exposure` as binary
#> [1] "binary"
density_meta(wt_ate(ps, trt))
#> ℹ Treating `.exposure` as binary
#> NULL

dose <- rnorm(20)
mu <- 0.3 * ps
w <- wt_ate(mu, dose, exposure_type = "continuous")

exposure_type(w)
#> [1] "continuous"
density_meta(w)
#> density:   normal
#> numerator: marginal
#> sigma:     pooled
```
