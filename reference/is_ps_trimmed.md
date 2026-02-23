# Test whether propensity scores have been trimmed

`is_ps_trimmed()` returns `TRUE` if `x` is a `ps_trim` object or a `psw`
object created from trimmed propensity scores, and `FALSE` otherwise.
This tests whether the *object* carries trimming information, not which
individual units were trimmed; see
[`is_unit_trimmed()`](https://r-causal.github.io/propensity/reference/is_unit_trimmed.md)
for that.

## Usage

``` r
is_ps_trimmed(x)
```

## Arguments

- x:

  An object to test.

## Value

A logical scalar (`TRUE` or `FALSE`).

## See also

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
for trimming propensity scores,
[`is_unit_trimmed()`](https://r-causal.github.io/propensity/reference/is_unit_trimmed.md)
to identify which units were trimmed,
[`ps_trim_meta()`](https://r-causal.github.io/propensity/reference/ps_trim_meta.md)
to retrieve full trimming metadata.

## Examples

``` r
ps <- c(0.05, 0.3, 0.6, 0.95)
trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)

is_ps_trimmed(trimmed)
#> [1] TRUE
is_ps_trimmed(ps)
#> [1] FALSE
```
