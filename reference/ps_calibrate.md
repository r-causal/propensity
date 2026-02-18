# Calibrate propensity scores

This function calibrates propensity scores to improve their accuracy
using either Platt scaling (logistic regression) or isotonic regression.
It preserves the attributes of causal weight objects when applicable.

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

  Numeric vector of propensity scores between 0 and 1

- .exposure:

  A binary vector of treatment assignments

- method:

  Calibration method:

  `"logistic"`

  :   (Default) Logistic calibration (also known as Platt scaling).
      Assumes a sigmoid relationship between observed and true
      probabilities. Best when: propensity scores follow a logistic
      pattern but are systematically biased. Provides smooth, parametric
      calibration. Faster and more stable with small samples.

  `"isoreg"`

  :   Isotonic regression calibration. Uses a non-parametric monotonic
      transformation. Best when: the relationship between observed and
      true probabilities is non-linear or when you want to preserve the
      rank order without assuming a specific functional form. More
      flexible but requires larger samples for stable estimates.

- smooth:

  Logical. For `method = "logistic"`, whether to use a smoothed logistic
  spline model (`smooth = TRUE`, default) or simple logistic regression
  (`smooth = FALSE`). When `TRUE`, uses
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) with spline
  smoothing. When `FALSE`, uses
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html). Ignored for
  `method = "isoreg"`.

- .focal_level:

  The value representing the focal group (typically treatment). If not
  provided, `ps_calibrate()` will attempt to automatically determine the
  coding.

- .reference_level:

  The value representing the reference group (typically control). If not
  provided, `ps_calibrate()` will attempt to automatically determine the
  coding.

- .treated:

  **\[deprecated\]** Use `.focal_level` instead.

- .untreated:

  **\[deprecated\]** Use `.reference_level` instead.

## Value

A calibrated propensity score object (`ps_calib`)

## Examples

``` r
# Generate example data
ps <- runif(100)
exposure <- rbinom(100, 1, ps)

# Logistic calibration with smoothing (default)
calibrated_smooth <- ps_calibrate(ps, exposure)
#> ℹ Setting focal level to 1

# Logistic calibration without smoothing (simple logistic regression)
calibrated_simple <- ps_calibrate(ps, exposure, smooth = FALSE)
#> ℹ Setting focal level to 1

# Isotonic regression
calibrated_iso <- ps_calibrate(ps, exposure, method = "isoreg")
#> ℹ Setting focal level to 1
```
