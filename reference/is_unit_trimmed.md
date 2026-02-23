# Identify which units were trimmed

`is_unit_trimmed()` returns a logical vector indicating which
observations were removed by trimming. This is a per-unit query, as
opposed to
[`is_ps_trimmed()`](https://r-causal.github.io/propensity/reference/is_ps_trimmed.md),
which tests whether the object has been trimmed at all.

## Usage

``` r
is_unit_trimmed(x)
```

## Arguments

- x:

  A `ps_trim` object created by
  [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md).

## Value

A logical vector the same length as `x`, where `TRUE` marks a trimmed
unit.

## See also

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
for trimming propensity scores,
[`is_ps_trimmed()`](https://r-causal.github.io/propensity/reference/is_ps_trimmed.md)
to test whether an object has been trimmed,
[`ps_trim_meta()`](https://r-causal.github.io/propensity/reference/ps_trim_meta.md)
to retrieve full trimming metadata.

## Examples

``` r
ps <- c(0.05, 0.3, 0.6, 0.95)
trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)

is_unit_trimmed(trimmed)
#> [1]  TRUE FALSE FALSE  TRUE

# Use to subset data to retained observations
kept <- !is_unit_trimmed(trimmed)
ps[kept]
#> [1] 0.3 0.6
```
