# Check if propensity scores are calibrated

`is_ps_calibrated()` tests whether `x` is a calibrated propensity score
object (class `ps_calib`) or a `psw` object derived from calibrated
scores.

## Usage

``` r
is_ps_calibrated(x)
```

## Arguments

- x:

  An object to test.

## Value

A single `TRUE` or `FALSE`.

## See also

[`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
to calibrate propensity scores.

## Examples

``` r
ps <- runif(100)
exposure <- rbinom(100, 1, ps)

is_ps_calibrated(ps)
#> [1] FALSE

calibrated <- ps_calibrate(ps, exposure, smooth = FALSE)
#> ℹ Setting focal level to 1
is_ps_calibrated(calibrated)
#> [1] TRUE
```
