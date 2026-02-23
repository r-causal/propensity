# Identify which units were truncated

`is_unit_truncated()` returns a logical vector indicating which
observations had their propensity scores modified by truncation. Use
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md)
to test whether an object has been truncated at all.

## Usage

``` r
is_unit_truncated(x)
```

## Arguments

- x:

  A `ps_trunc` object created by
  [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md).

## Value

A logical vector the same length as `x` (or number of rows for matrix
input). `TRUE` marks observations whose values were winsorized.

## See also

[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md),
[`ps_trunc_meta()`](https://r-causal.github.io/propensity/reference/ps_trunc_meta.md)

## Examples

``` r
ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
is_unit_truncated(ps_t)
#> [1]  TRUE FALSE FALSE FALSE  TRUE
```
