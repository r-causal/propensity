# Check if propensity scores have been refit

`is_refit()` tests whether `x` is a `ps_trim` object whose propensity
model has been refit on the retained (non-trimmed) observations via
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md).

## Usage

``` r
is_refit(x)
```

## Arguments

- x:

  An object to test (typically a
  [ps_trim](https://r-causal.github.io/propensity/reference/ps_trim.md)
  vector).

## Value

A single `TRUE` or `FALSE`.

## See also

[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md)
to refit a propensity model after trimming,
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
to trim propensity scores.

## Examples

``` r
set.seed(2)
n <- 30
x <- rnorm(n)
z <- rbinom(n, 1, plogis(0.4 * x))
fit <- glm(z ~ x, family = binomial)
ps <- predict(fit, type = "response")

trimmed <- ps_trim(ps, lower = 0.2, upper = 0.8)
is_refit(trimmed)
#> [1] FALSE

refit <- ps_refit(trimmed, fit)
is_refit(refit)
#> [1] TRUE
```
