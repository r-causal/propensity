# Check if object is truncated

`is_ps_truncated()` is an S3 generic that returns `TRUE` if its argument
represents a `ps_trunc` object or `psw` object created from truncated
propensity scores. `is_ps_truncated()` is a question about whether the
propensity scores *have* been truncated, as opposed to
[`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md),
which is a question about which *units* have been truncated.

## Usage

``` r
is_ps_truncated(x)
```

## Arguments

- x:

  An object.

## Value

A logical scalar (`TRUE` or `FALSE`).
