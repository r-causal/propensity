# Check if propensity scores are calibrated

`is_ps_calibrated()` tests whether `x` was calibrated, rather than
whether it is still a `ps_calib` object: it is `TRUE` for a calibrated
propensity score, for weights derived from one, and for scores that were
calibrated before
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
or
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
bounded them, each of which records the calibration in its own metadata.

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
is_ps_calibrated(calibrated)
#> [1] TRUE
```
