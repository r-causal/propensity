# Calibrate propensity scores

`ps_calibrate()` adjusts estimated propensity scores so they better
reflect true treatment probabilities. This can improve the accuracy of
inverse probability weights derived from a misspecified propensity score
model.

## Usage

``` r
ps_calibrate(
  ps,
  .exposure,
  method = c("logistic", "isoreg"),
  smooth = TRUE,
  .focal_level = NULL,
  .reference_level = NULL,
  .treated = NULL,
  .untreated = NULL
)
```

## Arguments

- ps:

  A numeric vector of propensity scores between 0 and 1. Must not
  already be calibrated.

- .exposure:

  A binary vector of observed treatment assignments, the same length as
  `ps`.

- method:

  Calibration method. One of:

  `"logistic"`

  :   (Default) Logistic calibration, also called Platt scaling. Fits a
      logistic regression of `.exposure` on `ps`, yielding a smooth,
      parametric correction. Works well with small samples and when the
      bias in `ps` is approximately monotone.

  `"isoreg"`

  :   Isotonic regression. Fits a non-parametric, monotone step
      function. More flexible than logistic calibration because it makes
      no distributional assumption, but needs larger samples for stable
      estimates.

- smooth:

  Logical. When `method = "logistic"`, controls the form of the
  calibration model. If `TRUE` (default), fits a GAM with a spline on
  `ps` via [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html); if
  `FALSE`, fits a simple logistic regression via
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html). Ignored when
  `method = "isoreg"`.

- .focal_level:

  The value of `.exposure` representing the focal group (typically the
  treated group). If `NULL` (default), coding is determined
  automatically.

- .reference_level:

  The value of `.exposure` representing the reference group (typically
  the control group). If `NULL` (default), coding is determined
  automatically.

- .treated:

  **\[deprecated\]** Use `.focal_level` instead.

- .untreated:

  **\[deprecated\]** Use `.reference_level` instead.

## Value

A `ps_calib` vector the same length as `ps`. The attribute
`ps_calib_meta` stores calibration metadata (method and whether
smoothing was applied). Use
[`is_ps_calibrated()`](https://r-causal.github.io/propensity/reference/is_ps_calibrated.md)
to test whether an object has been calibrated.

## Details

Calibration is useful when the propensity score model is correctly
specified in terms of variable selection but produces probabilities that
are systematically too high or too low. Unlike
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
and
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
which handle extreme scores by removing or bounding them, calibration
reshapes the entire distribution of scores.

**Choosing a method:**

- Use `"logistic"` (the default) as a first choice. It is stable and
  fast, and the `smooth = TRUE` option adds flexibility via a spline.

- Use `"isoreg"` when you suspect a non-smooth or irregular relationship
  between estimated and true probabilities and have a sufficiently large
  sample.

The calibrated scores are returned as a `ps_calib` object, which can be
passed directly to weight functions such as
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md).

## References

Platt, J. (1999). Probabilistic outputs for support vector machines and
comparisons to regularized likelihood methods. *Advances in Large Margin
Classifiers*, 61–74.

Zadrozny, B., & Elkan, C. (2002). Transforming classifier scores into
accurate multiclass probability estimates. *Proceedings of the Eighth
ACM SIGKDD International Conference on Knowledge Discovery and Data
Mining*, 694–699.
[doi:10.1145/775047.775151](https://doi.org/10.1145/775047.775151)

## See also

[`is_ps_calibrated()`](https://r-causal.github.io/propensity/reference/is_ps_calibrated.md)
to test for calibrated scores;
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
and
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
for alternative approaches to extreme propensity scores;
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and other weight functions that accept `ps_calib` objects.

## Examples

``` r
# Simulate data
set.seed(42)
ps <- runif(200)
exposure <- rbinom(200, 1, ps)

# Logistic calibration without smoothing (simple Platt scaling)
cal <- ps_calibrate(ps, exposure, smooth = FALSE)
#> ℹ Setting focal level to 1
cal
#> <ps_calib[200]; method=logistic>
#>  [1] 0.8768754 0.8875980 0.2784575 0.8280651 0.6675032 0.5319951 0.7570779
#>  [8] 0.1604955 0.6830053 0.7291950
#> # ... with 190 more values

# Use calibrated scores to calculate weights
wt_ate(cal, exposure)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = ate; calibrated}[200]>
#>   [1]  1.140413  1.126636  3.591212  1.207634  1.498120  1.879717  4.116546
#>   [8]  1.191179  1.464118  1.371375  1.855251  1.347954  1.128055  1.334694
#>  [15]  1.873492  1.124922  1.104636  1.176539  2.079333  2.376276  1.147607
#>  [22]  1.194797  1.099587  9.255849  1.150056  2.111273  1.625281  1.146443
#>  [29]  1.813578  1.202353  1.319373  1.227172  1.619236  1.407268 10.590004
#>  [36]  1.205271  1.105925  4.728705  1.145858  1.572381  1.595167  1.772408
#>  [43]  1.121790  1.106935  1.758141  1.115151  1.159179  1.502218  1.108219
#>  [50]  2.805228  1.480542  1.511163  1.649763  1.256712  1.122644  4.298055
#>  [57]  1.422451  1.226540  1.343593  2.112310  3.348883  1.102432  1.288466
#>  [64]  2.416129  1.189910  1.246501  3.775976  1.209850  3.548601  4.201312
#>  [71]  1.124970  1.196402  1.279265  1.945603  4.910192  1.347561  1.106196
#>  [78]  1.584041  1.899052  1.103132  1.658348  1.212932  1.541116  3.044056
#>  [85]  1.267491  1.715517  1.302617  1.155398  1.152281  1.421620  3.261444
#>  [92]  1.102497  1.269325  1.129032  1.133530  1.324601  3.084413  2.115670
#>  [99]  1.310064  1.553123  1.535242  1.280267  1.279501  1.621643  1.123516
#> [106]  1.112496  1.316044  1.325880  1.814291  1.103468  1.579973  5.960174
#> [113]  1.299399  1.835609  1.814183  2.237296  1.103041  2.877072  2.749961
#> [120]  5.782649  1.535360  1.687422  1.683636  1.634153  1.347076  2.564296
#> [127]  1.137578  1.112516  1.302365  3.946616  7.762300  1.594854  1.522340
#> [134]  1.126454  1.189213  1.663814  1.216527  1.173481  4.547330  1.541815
#> [141]  1.203794  1.148552  2.135436  1.263130  1.325454  1.220759  5.436839
#> [148]  1.122225  3.502822  5.892699  1.347525  1.460171  1.263813  1.637690
#> [155]  3.381630  1.267489  1.244673  1.117167  1.192109  1.416832  1.127966
#> [162]  2.314895  1.599584  4.917740  2.225077  1.235420  2.201430  3.245145
#> [169]  1.175468  1.242677  1.331237  1.691377  2.431857  1.949609  2.345311
#> [176]  6.186653  1.213263  1.626478  1.245416  4.622547  1.138312  1.178849
#> [183]  1.445244  1.340694  1.312530  4.276192  1.138410  5.051990  1.189997
#> [190]  1.388812  4.960086  1.257408  1.186110  1.186298  1.143134  1.130987
#> [197]  2.206125  6.802145  1.311197  1.328810

# Isotonic regression calibration
cal_iso <- ps_calibrate(ps, exposure, method = "isoreg")
#> ℹ Setting focal level to 1

if (rlang::is_installed("mgcv")) {
  # Logistic calibration with spline smoothing (default)
  cal_smooth <- ps_calibrate(ps, exposure)
}
#> ℹ Setting focal level to 1
```
