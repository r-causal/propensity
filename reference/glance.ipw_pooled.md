# Glance at a pooled inverse probability weighted result

[`glance()`](https://generics.r-lib.org/reference/glance.html) describes
a
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
result rather than its estimates: one row naming the estimand, counting
the observations and the results that were pooled, and reporting the
complete-data degrees of freedom the small-sample adjustment used. A
pooled result reporting several effect measures, or several contrasts of
a categorical exposure, still returns exactly one row.

## Usage

``` r
# S3 method for class 'ipw_pooled'
glance(x, ...)
```

## Arguments

- x:

  An `ipw_pooled` object, as returned by
  [`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html).

- ...:

  These dots are for future extensions and must be empty.

## Value

A one-row [tibble](https://tibble.tidyverse.org/reference/tibble.html)
with the columns:

- `estimand`:

  The causal estimand every pooled result targeted.

- `nobs`:

  The smallest number of observations any pooled result was estimated
  from. The imputed datasets are the same size, so this is that size
  unless a result was fitted on a subset of one. Reported by
  [`stats::nobs()`](https://rdrr.io/r/stats/nobs.html).

- `m`:

  The number of results pooled.

- `dfcom`:

  The complete-data degrees of freedom the Barnard-Rubin adjustment
  used.

## See also

[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.html)
for the pooling,
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw_pooled.md)
for the pooled estimates, and the Multiple imputation section of
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

glance(pool_ipw(fits))
#> # A tibble: 1 × 4
#>   estimand  nobs     m dfcom
#>   <chr>    <int> <int> <dbl>
#> 1 ate        150     2   141
```
