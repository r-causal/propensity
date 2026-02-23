# Test whether propensity scores have been truncated

`is_ps_truncated()` returns `TRUE` if `x` is a `ps_trunc` object or a
`psw` object derived from truncated propensity scores. Use
[`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md)
to find out *which* observations were modified.

## Usage

``` r
is_ps_truncated(x)
```

## Arguments

- x:

  An object.

## Value

A single `TRUE` or `FALSE`.

## See also

[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
[`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md),
[`ps_trunc_meta()`](https://r-causal.github.io/propensity/reference/ps_trunc_meta.md)

## Examples

``` r
ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
is_ps_truncated(ps)
#> [1] FALSE

ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
is_ps_truncated(ps_t)
#> [1] TRUE
```
