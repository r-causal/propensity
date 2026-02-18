# Check if units have been trimmed

`is_unit_trimmed()` is an that vector of `TRUE` or `FALSE` values,
representing if the unit was trimmed. `is_unit_trimmed()` is a question
about which *units* have been trimmed, as opposed to
[`is_ps_trimmed()`](https://r-causal.github.io/propensity/reference/is_ps_trimmed.md),
which is a question about whether or not the propensity scores *have*
been trimmed.

## Usage

``` r
is_unit_trimmed(x)
```

## Arguments

- x:

  An object.

## Value

A logical scalar (`TRUE` or `FALSE`).
