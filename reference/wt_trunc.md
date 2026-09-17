# Truncate (Winsorize) Propensity Score Weights

`wt_trunc()` bounds extreme weights at a limit read on the scale of the
weights themselves, replacing each weight beyond the limit with the
limit. No unit is removed. Where
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
bounds the propensity scores the weights are built from, `wt_trunc()`
bounds the weights after they are built, so it applies to weights for
any exposure type, including the density-ratio weights of a continuous
exposure, which have no propensity score to bound.

## Usage

``` r
wt_trunc(
  .weights,
  method = c("adaptive", "wt", "pctl", "count"),
  lower = NULL,
  upper = NULL,
  ...
)
```

## Arguments

- .weights:

  A [psw](https://r-causal.github.io/propensity/reference/psw.md)
  vector, such as one returned by
  [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
  or a numeric vector of weights. A numeric vector is returned as a
  [psw](https://r-causal.github.io/propensity/reference/psw.md) vector
  with no estimand. Propensity scores, including those returned by
  [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
  [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
  and
  [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md),
  are refused with an error of class `propensity_type_error`: build
  weights from them first.

- method:

  The rule that sets the bound:

  - `"adaptive"` (default): an upper bound of \\\sqrt{n} \log(n) / 5\\
    (Gruber et al., 2022), where \\n\\ is the number of weights present.

  - `"wt"`: bounds given as weight values.

  - `"pctl"`: bounds read at sample quantiles of the weights.

  - `"count"`: bounds that move a given number of the most extreme
    weights.

- lower, upper:

  The bounds, read according to `method`:

  - `"adaptive"`: not used. Supplying either is ignored with a warning.

  - `"wt"`: `upper` is required and must be a positive finite number.
    `lower` is optional and must be a finite number of at least 0 and
    below `upper`.

  - `"pctl"`: `upper` is required and must lie in (0.5, 1). `lower` is
    optional and must lie in (0, 0.5).

  - `"count"`: `upper` is required and must be a whole number of at
    least 1; the `upper` largest weights are bounded at the next
    largest. `lower` is optional and bounds the `lower` smallest weights
    at the next smallest in the same way.

- ...:

  Not used. Any argument supplied here is an error, so that a misspelled
  bound is not silently ignored.

## Value

A [psw](https://r-causal.github.io/propensity/reference/psw.md) vector
of the same length as `.weights`, with the bounded values, the estimand
labeled as truncated, and a record of the truncation.

## Details

### Methods

`"adaptive"` is the default. Its bound depends on the sample size alone
and loosens as the sample grows: at \\n = 1000\\ it is 43.7. Gruber et
al. (2022) derived it for the propensity score of a binary treatment,
choosing a bound that loosens with the sample so that its bias shrinks.
Whether the rule carries over to the density-ratio weights of a
continuous exposure has not been established in the literature. Here
\\n\\ counts the weights that are present, so for weights built from a
trimmed propensity score, whose trimmed units are `NA`, it is the number
of units kept by the trim. The bound is below 1 for \\n \le 6\\ (0.72 at
\\n = 5\\), so for a sample that small `"adaptive"` pulls every weight
above that value down to a constant below 1.

`"wt"` winsorizes at the values given, and `"count"` at the order
statistic next to the weights it bounds. Under `"count"`, a weight tied
with that order statistic is not moved, so fewer weights than `upper`
may change. The two counts must leave at least two weights between them,
and `upper` alone must leave at least one weight below it.

`"pctl"` reads its bounds with
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) at its
default type, as `ps_trunc(method = "pctl")` does. Because such a rule
moves a fixed share of the units however the weights are distributed, it
is best used to check how sensitive an analysis is to its bound rather
than as the bound itself. An `upper` below 0.5 is refused rather than
read as the other tail: write `upper = 0.9` to bound the upper tail at
the 0.9 quantile, or `lower = 0.1` to bound the lower tail at the 0.1
quantile.

`"adaptive"` and `"pctl"` need at least two weights present, and every
method leaves a missing weight missing.

### One-sided by default

Each method bounds only the upper tail unless `lower` is supplied, and
`"adaptive"`, which takes no `lower`, always bounds only the upper tail.
The largest weights are the ones that dominate an estimate; small
weights contribute little to it, so bounding them changes little.
`lower` bounds the lower tail when that is wanted.

### What is recorded

The result keeps every attribute of `.weights`: the stabilization status
and any numerator model, the density record of a continuous exposure,
the attributes of a categorical exposure, and the records left by
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
or
[`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md).
It gains a `wt_truncated` flag, read by
[`is_wt_truncated()`](https://r-causal.github.io/propensity/reference/is_wt_truncated.md),
and a record of the method, the bounds as given and as applied, and the
positions of the weights moved, read by
[`is_unit_wt_truncated()`](https://r-causal.github.io/propensity/reference/is_wt_truncated.md).
The record is positional and is dropped by any operation that changes
the length of the weights; see
[psw](https://r-causal.github.io/propensity/reference/psw.md).

The estimand gains the label `"; weights truncated"`, as in
`"ate; weights truncated"`. A truncated weight targets a different
quantity from the weight before truncation: any fixed bound introduces a
bias (Ma and Wang, 2020; Gruber et al., 2022), and an interval computed
with the truncated weights held fixed is an interval for that quantity.
A weight truncation does not describe the same operation as truncating
the propensity score, which is labeled `"; truncated"` and reported by
[`is_ps_truncated()`](https://r-causal.github.io/propensity/reference/is_ps_truncated.md).

Weights that have already been truncated are returned unchanged with a
warning of class `propensity_already_modified_warning`. One record
describes one truncation, so to truncate at a different bound, call
`wt_trunc()` on the original weights.

### Weights from a trimmed, truncated, or calibrated score

Weights built from a score modified by
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
or
[`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
can be truncated as well. The result carries both records, and the
labels stack, as in `"ate; trimmed; weights truncated"`. The package
does not prevent combining these tools.

### Effect estimation

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
does not accept weights whose values were truncated, because its
standard errors do not account for the bound. An analysis with truncated
weights needs an interval built another way: one computed with the
weights held fixed, which conditions on the bound applied, or a
bootstrap that repeats the truncation in every resample, which a
percentile bound in particular needs, since that bound is itself
estimated from the weights.

To see how much a truncation changed the effective sample size, add the
weights before and after it to a data frame as separate columns and pass
both to `halfmoon::check_ess()`.

### Truncating weights or trimming the density

For a continuous exposure there are two ways to hold down extreme
weights. `wt_trunc()` bounds the weights once they are built and keeps
every unit.
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
with `method = "density"` or `method = "resid"` sets aside the units
whose observed dose is implausible under the dose model, and
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md)
refits the model on the units that remain; see **Trimming a dose model**
in
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md).
Both change the estimand. A fixed bound introduces a bias that does not
vanish however large the sample, while a bound that loosens with the
sample, as `"adaptive"` does, leaves the estimand alone in the limit. A
trim describes only the units it keeps. `wt_trunc()` applies to weights
from any density family, and so does `"density"`; `"resid"` accepts the
normal, t, and Laplace families only.
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
refuses the weights from either tool.

A simulation study run for this package (samples of 1000, 1000
replications, a correctly specified mean for the dose) found that which
tool helps depends on whether the density family fits the residuals. The
coverage figures below are for a sandwich interval that holds the
weights fixed, unless stated. With normal residuals, a normal density,
and an effective sample size of 64% of the sample, the `"adaptive"`
bound moved almost no weight and covered 0.933, as the untruncated
weights did under the same interval. The M-estimation interval on the
untruncated weights covered 0.863, so the gain in coverage there comes
from the interval, not the bound. What the bound bought was precision: a
root mean squared error of 0.084 against 0.096, for a bias of 0.015 (3%
of the effect) against 0.012. At an effective sample size of 53% it
covered 0.847 against 0.839 untruncated, and at 81% it changed almost
nothing (it applied in 0.2% of replications, leaving coverage and RMSE
unchanged). A density trim at `lower = 0.01` covered 0.843 at 64%, below
the untruncated weights.

With t(4) or Laplace residuals, or a residual spread that varied with a
covariate, the normal density fit some units badly. The bound then
helped little (coverage 0.822, 0.906, and 0.644, against 0.788, 0.881,
and 0.618 untruncated), while the 1% density trim covered **0.989**,
**0.980**, and **0.979**, because it removed the units that density fit
worst. The trim's intervals were computed on the trimmed sample and were
conservative, and its bias was measured against the effect in the full
population, so part of that bias is the change of estimand. The
`"resid"` method was not simulated; these results carry over to it only
because a bound on the standardized residual is a floor on the
conditional density.

Neither tool is the first step. Under the same heavy-tailed residuals,
`dens_t(4)` at its default scale covered 0.926 to 0.939 with the
M-estimation interval, whose standard errors were close to right, with
about a third of the root mean squared error of the bounded normal
weights, and kept every unit. So first check the density family against
the residuals, for instance by comparing the weights to those of a
[`dens_kernel()`](https://r-causal.github.io/propensity/reference/dens_normal.md)
density or a heavier tail. If the weights are still heavy under a family
that fits, bound them with `wt_trunc()`. Trim the dose model when the
units at the low-density end are ones the analysis should not describe,
and state the population the trimmed analysis describes.

## References

Gruber, S., Phillips, R. V., Lee, H., & van der Laan, M. J. (2022).
Data-adaptive selection of the propensity score truncation level for
inverse-probability-weighted and targeted maximum likelihood estimators
of marginal point treatment effects. *American Journal of Epidemiology*,
191(9), 1640–1651.

Ma, X., & Wang, J. (2020). Robust inference using inverse probability
weighting. *Journal of the American Statistical Association*, 115(532),
1851–1860.

## See also

[`is_wt_truncated()`](https://r-causal.github.io/propensity/reference/is_wt_truncated.md)
and
[`is_unit_wt_truncated()`](https://r-causal.github.io/propensity/reference/is_wt_truncated.md)
to read the record,
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
to bound propensity scores instead, and
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
to remove units.

## Examples

``` r
set.seed(1)
n <- 50
x <- rnorm(n)
z <- rbinom(n, 1, plogis(2 * x))
fit <- glm(z ~ x, family = binomial)
w <- wt_ate(fit, .exposure = z)
#> ℹ Treating `.exposure` as binary

# The adaptive bound, on the upper tail
w_adaptive <- wt_trunc(w)
w_adaptive
#> <psw{estimand = ate; weights truncated}[50]>
#>  [1] 1.212205 2.278488 1.151519 5.532436 2.010846 1.155264 1.783866 1.523330
#>  [9] 2.470781 1.355873 5.532436 2.090212 1.213994 1.016444 1.280802 1.541310
#> [17] 2.763810 1.375881 1.457933 1.660358 3.556024 1.487681 2.523988 1.023638
#> [25] 1.633358 1.531639 1.452808 1.054486 1.269447 1.876674 1.192720 1.493157
#> [33] 1.920466 2.873950 1.063360 1.298292 1.308405 1.528920 4.421230 1.502801
#> [41] 3.239700 1.386973 1.559374 1.701168 1.191948 1.186242 1.955335 1.498482
#> [49] 1.485625 3.404809
#> truncation: adaptive (upper 5.53), 2 of 50 weights truncated
is_wt_truncated(w_adaptive)
#> [1] TRUE
which(is_unit_wt_truncated(w_adaptive))
#> [1]  4 11

# A bound given as a weight, on one tail or both
wt_trunc(w, method = "wt", upper = 5)
#> <psw{estimand = ate; weights truncated}[50]>
#>  [1] 1.212205 2.278488 1.151519 5.000000 2.010846 1.155264 1.783866 1.523330
#>  [9] 2.470781 1.355873 5.000000 2.090212 1.213994 1.016444 1.280802 1.541310
#> [17] 2.763810 1.375881 1.457933 1.660358 3.556024 1.487681 2.523988 1.023638
#> [25] 1.633358 1.531639 1.452808 1.054486 1.269447 1.876674 1.192720 1.493157
#> [33] 1.920466 2.873950 1.063360 1.298292 1.308405 1.528920 4.421230 1.502801
#> [41] 3.239700 1.386973 1.559374 1.701168 1.191948 1.186242 1.955335 1.498482
#> [49] 1.485625 3.404809
#> truncation: wt (upper 5), 2 of 50 weights truncated
wt_trunc(w, method = "wt", lower = 1.05, upper = 5)
#> <psw{estimand = ate; weights truncated}[50]>
#>  [1] 1.212205 2.278488 1.151519 5.000000 2.010846 1.155264 1.783866 1.523330
#>  [9] 2.470781 1.355873 5.000000 2.090212 1.213994 1.050000 1.280802 1.541310
#> [17] 2.763810 1.375881 1.457933 1.660358 3.556024 1.487681 2.523988 1.050000
#> [25] 1.633358 1.531639 1.452808 1.054486 1.269447 1.876674 1.192720 1.493157
#> [33] 1.920466 2.873950 1.063360 1.298292 1.308405 1.528920 4.421230 1.502801
#> [41] 3.239700 1.386973 1.559374 1.701168 1.191948 1.186242 1.955335 1.498482
#> [49] 1.485625 3.404809
#> truncation: wt (lower 1.05, upper 5), 4 of 50 weights truncated

# The 99th percentile, as a sensitivity check
wt_trunc(w, method = "pctl", upper = 0.99)
#> <psw{estimand = ate; weights truncated}[50]>
#>  [1] 1.212205 2.278488 1.151519 8.127027 2.010846 1.155264 1.783866 1.523330
#>  [9] 2.470781 1.355873 7.639676 2.090212 1.213994 1.016444 1.280802 1.541310
#> [17] 2.763810 1.375881 1.457933 1.660358 3.556024 1.487681 2.523988 1.023638
#> [25] 1.633358 1.531639 1.452808 1.054486 1.269447 1.876674 1.192720 1.493157
#> [33] 1.920466 2.873950 1.063360 1.298292 1.308405 1.528920 4.421230 1.502801
#> [41] 3.239700 1.386973 1.559374 1.701168 1.191948 1.186242 1.955335 1.498482
#> [49] 1.485625 3.404809
#> truncation: pctl 0.99 (upper 8.13), 1 of 50 weights truncated

# The three largest weights, bounded at the fourth largest
wt_trunc(w, method = "count", upper = 3)
#> <psw{estimand = ate; weights truncated}[50]>
#>  [1] 1.212205 2.278488 1.151519 3.556024 2.010846 1.155264 1.783866 1.523330
#>  [9] 2.470781 1.355873 3.556024 2.090212 1.213994 1.016444 1.280802 1.541310
#> [17] 2.763810 1.375881 1.457933 1.660358 3.556024 1.487681 2.523988 1.023638
#> [25] 1.633358 1.531639 1.452808 1.054486 1.269447 1.876674 1.192720 1.493157
#> [33] 1.920466 2.873950 1.063360 1.298292 1.308405 1.528920 3.556024 1.502801
#> [41] 3.239700 1.386973 1.559374 1.701168 1.191948 1.186242 1.955335 1.498482
#> [49] 1.485625 3.404809
#> truncation: count 3 (upper 3.56), 3 of 50 weights truncated

# Weights for a continuous exposure
a <- 1 + 0.5 * x + rnorm(n)
dose <- lm(a ~ x)
w_dose <- wt_ate(dose, .exposure = a)
#> ℹ Treating `.exposure` as continuous
wt_trunc(w_dose)
#> <psw{estimand = ate; weights truncated; stabilized}[50]>
#>  [1] 0.8785574 0.9114861 0.6906516 0.6455340 1.0463966 0.5023458 0.8706861
#>  [8] 0.6045779 3.1339880 1.0193753 0.4955401 0.8296752 0.6312966 0.8002803
#> [15] 0.6251934 0.8386485 1.2895826 0.4850001 0.5878021 0.7993661 0.5804047
#> [22] 2.8493299 0.9049329 0.0812741 1.1002735 0.8494276 0.7988081 0.2190744
#> [29] 0.8113444 1.1431540 0.2551815 1.0066451 0.7982113 0.8629044 5.5324360
#> [36] 0.6677786 0.6846820 1.6158321 1.5777330 0.9343494 0.7724890 0.7433665
#> [43] 0.9745417 0.7112151 0.6488288 0.5576057 0.9507051 0.9397858 0.7858439
#> [50] 0.8673159
#> density:   normal
#> numerator: marginal
#> sigma:     pooled
#> truncation: adaptive (upper 5.53), 1 of 50 weights truncated

# Weights from a trimmed propensity score carry both records
trimmed <- ps_refit(ps_trim(fit, lower = 0.05, upper = 0.95), fit)
wt_trunc(wt_ate(trimmed, .exposure = z), method = "wt", upper = 5)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate; trimmed; weights truncated}[50]>
#>  [1] 1.218631 2.270252 1.157053 5.000000 2.008570 1.160864 1.785672 1.528346
#>  [9] 2.463665 1.363270 5.000000 2.090773 1.220441       NA 1.286673 1.548425
#> [17] 2.742373 1.381740 1.463429 1.663916 3.518597 1.492978 2.509412       NA
#> [25] 1.637246 1.538800 1.460229 1.057521 1.276421 1.876934 1.198080 1.500474
#> [33] 1.919935 2.849167 1.066708 1.305456 1.315623 1.536094 4.353572 1.507984
#> [41] 3.203032 1.394423 1.564062 1.704189 1.198117 1.192332 1.954148 1.503699
#> [49] 1.492965 3.372186
#> truncation: wt (upper 5), 2 of 50 weights truncated
```
