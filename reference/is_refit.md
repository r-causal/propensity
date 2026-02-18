# Check if an object has been refit

**`is_refit()`** is an S3 generic that returns `TRUE` if its argument
represents a
[ps_trim](https://r-causal.github.io/propensity/reference/ps_trim.md)
object (or a weighting object) that has had the propensity model refit
on the retained subset.

## Usage

``` r
is_refit(x)
```

## Arguments

- x:

  An R object (e.g. a
  [ps_trim](https://r-causal.github.io/propensity/reference/ps_trim.md)
  or [psw](https://r-causal.github.io/propensity/reference/psw.md)).

## Value

A logical scalar (`TRUE` or `FALSE`).
