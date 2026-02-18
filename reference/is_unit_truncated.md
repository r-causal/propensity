# Check if units have been truncated

[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md)
is an S3 generic that returns a vector of `TRUE` or `FALSE`,
representing if the element has been truncated. `is_unit_truncated()` is
a question about which *units* have been truncated, as opposed to
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md),
which is a question about whether the propensity scores *have* been
truncated.

## Usage

``` r
is_unit_truncated(x)
```

## Arguments

- x:

  An object.

## Value

A logical vector.
