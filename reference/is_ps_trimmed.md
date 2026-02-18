# Check if object is trimmed

`is_ps_trimmed()` is an S3 generic that returns `TRUE` if its argument
represents a `ps_trim` object or `psw` object created from trimmed
propensity scores. `is_ps_trimmed()` is a question about whether or not
the propensity scores *have* been trimmed, as opposed to
[`is_unit_trimmed()`](https://r-causal.github.io/propensity/reference/is_unit_trimmed.md),
which is a question about which *units* have been trimmed.

## Usage

``` r
is_ps_trimmed(x)
```

## Arguments

- x:

  An object.

## Value

A logical scalar (`TRUE` or `FALSE`).
