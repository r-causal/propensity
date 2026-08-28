# Getting Started with propensity

This vignette walks through the core propensity score weighting
workflow: fitting a propensity score model, calculating weights, and
estimating causal effects with
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).
We’ll also cover what to do when propensity scores are extreme.

## Setup

``` r

library(propensity)
```

We’ll work with a simulated dataset throughout. There are two
confounders (`x1` and `x2`), a binary exposure (`z`), and a binary
outcome (`y`):

``` r

set.seed(42)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n)
z <- rbinom(n, 1, plogis(0.5 * x1 + 0.3 * x2))
y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1 + 0.2 * x2))
dat <- data.frame(x1, z, y, x2)
```

Both `x1` and `x2` affect treatment and outcome, so we need to adjust
for them.

## Basic workflow

### Step 1: Fit a propensity score model

Start with a model for treatment assignment. Here we use logistic
regression:

``` r

ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
```

### Step 2: Calculate weights and fit a weighted outcome model

Pass the fitted model directly to
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
to get ATE weights. It pulls out the fitted values and exposure for you:

``` r

wts <- wt_ate(ps_mod)
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
outcome_mod <- glm(y ~ z, data = dat, family = binomial(), weights = wts)
#> Warning in eval(family$initialize): non-integer #successes in a binomial glm!
```

[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
returns a `psw` object, which is just a numeric vector with some extra
metadata attached:

``` r

estimand(wts)
#> [1] "ate"
is_stabilized(wts)
#> [1] FALSE
```

You can also pass propensity scores as a plain numeric vector. In that
case you need to supply the exposure too:

``` r

ps <- fitted(ps_mod)
wt_ate(ps, dat$z)
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[100]>
#>   [1]  1.237569  1.962759  2.211732  1.312977  1.974772  1.918957  3.413991  1.844849  1.223426
#>  [10]  2.048453  1.409967  1.189795  1.283684  2.580633  1.439961  1.771951  1.627989  3.438494
#>  [19]  1.092310  1.379591  1.414973  1.142879  2.132832  2.539924  1.264028  1.584122  1.614753
#>  [28]  1.115628  2.235160  1.641530  1.598952  1.767794  1.494051  2.039262  3.465881  1.174226
#>  [37]  1.511863  1.832668  1.135144  2.045876  2.067593  2.960898  1.724205  2.807457  1.296458
#>  [46]  1.487979  1.433057  3.287998  2.085343  2.000254  1.845028  1.286187  1.207434  2.360698
#>  [55]  1.840088  1.704295  1.642486  2.362152 12.582758  2.974447  1.677742  1.704949  2.553764
#>  [64]  1.438721  1.711034  1.227343  1.812465  1.409825  1.518867  3.314572  1.404951  1.799540
#>  [73]  2.354036  1.941761  1.909359  1.731474  2.080547  2.731912  1.606549  3.350612  1.327948
#>  [82]  2.103802  2.178471  2.018730  3.813295  1.864473  2.078958  1.959235  1.747083  1.907159
#>  [91]  3.853789  1.584359  2.693732  1.644175  1.286716  1.788770  3.037240  1.416308  1.474800
#> [100]  1.619529
```

### Step 3: Estimate causal effects

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
takes the propensity score model and the weighted outcome model and
returns causal effect estimates. By default, the standard errors come
from M-estimation: the propensity score and outcome models are stacked
into a single system of estimating equations, so the uncertainty of
estimating the propensity scores is carried into the standard errors
rather than ignored. Set `se_method = "linearization"` to use the
influence-function method instead:

``` r

result <- ipw(ps_mod, outcome_mod)
result
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z, family = binomial(), data = dat, weights = wts) 
#> 
#> Marginal estimates:
#>                estimate  std.err      z ci.lower ci.upper conf.level   p.value    
#> mean 0         0.321135 0.069594 4.6144  0.18473  0.45754       0.95 3.942e-06 ***
#> mean 1         0.641133 0.075977 8.4385  0.49222  0.79004       0.95 < 2.2e-16 ***
#> rd 1 vs 0      0.319997 0.103584 3.0892  0.11698  0.52302       0.95  0.002007 ** 
#> log(rr) 1 vs 0 0.691374 0.248114 2.7865  0.20508  1.17767       0.95  0.005328 ** 
#> log(or) 1 vs 0 1.328843 0.461758 2.8778  0.42381  2.23387       0.95  0.004005 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Choosing an estimand

Each estimand targets a different population:

| Estimand | Target population | Function |
|----|----|----|
| ATE | Entire study population | [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md) |
| ATT | Treated (focal) group | [`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md) |
| ATU | Untreated (reference) group | [`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md) |
| ATO | Overlap population | [`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md) |
| ATM | Matched population | [`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md) |
| Entropy | Entropy-tilted population | [`wt_entropy()`](https://r-causal.github.io/propensity/reference/wt_ate.md) |

[`wt_atc()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
is an alias for
[`wt_atu()`](https://r-causal.github.io/propensity/reference/wt_ate.md).

ATE is the most common choice. ATT and ATU narrow the question to the
treated or untreated, respectively. ATO, ATM, and entropy weights target
overlap populations – they produce bounded weights by construction,
which makes them a good option when propensity scores are extreme (more
on that below).

To switch estimands, just swap the weight function:

``` r

wts_ate <- wt_ate(ps_mod)
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
wts_att <- wt_att(ps_mod)
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
wts_ato <- wt_ato(ps_mod)
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
```

## Handling extreme weights

Propensity scores near 0 or 1 produce large weights that can blow up
your variance. The [`summary()`](https://rdrr.io/r/base/summary.html)
method gives a quick look at the weight distribution:

``` r

summary(wts_ate)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.092   1.440   1.780   2.028   2.111  12.583
```

If you see a very large maximum or high variance, you have a few
options.

### Overlap estimands

The easiest fix is to use an estimand with bounded weights.
[`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and
[`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
down-weight observations where overlap is poor:

``` r

summary(wt_ato(ps_mod))
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.08451 0.30539 0.43830 0.43370 0.52629 0.92053
summary(wt_atm(ps_mod))
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.09231 0.43965 0.78036 0.70946 1.00000 1.00000
```

The trade-off is that you’re now targeting a different population.

### Trimming

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
drops observations with extreme propensity scores by setting them to
`NA`. The `"ps"` method uses fixed thresholds (by default, 0.1 and 0.9):

``` r

ps_trimmed <- ps_trim(ps, method = "ps")
```

The `"adaptive"` method (Crump et al., 2009) finds a data-driven
threshold:

``` r

ps_trimmed_adapt <- ps_trim(ps, method = "adaptive")
```

You can inspect the result with a few helpers:

``` r

# Confirm the object has been trimmed
is_ps_trimmed(ps_trimmed)
#> [1] TRUE

# Which observations were removed?
sum(is_unit_trimmed(ps_trimmed))
#> [1] 2

# View trimming metadata (method, cutoffs, etc.)
ps_trim_meta(ps_trimmed)
#> $method
#> [1] "ps"
#> 
#> $lower
#> [1] 0.1
#> 
#> $upper
#> [1] 0.9
#> 
#> $keep_idx
#>   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  20  21  22  23  24  25  26 
#>   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  20  21  22  23  24  25  26 
#>  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51 
#>  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51 
#>  52  53  54  55  56  57  58  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77 
#>  52  53  54  55  56  57  58  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77 
#>  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 
#>  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 
#> 
#> $trimmed_idx
#> [1] 19 59
#> 
#> $n_obs
#> [1] 100
```

Use `!is_unit_trimmed()` to subset your data down to the retained
observations:

``` r

retained <- !is_unit_trimmed(ps_trimmed)
dat_trimmed <- dat[retained, ]
```

After trimming, you should refit the propensity score model on the
retained sample so the scores reflect the trimmed population:

``` r

ps_refitted <- ps_refit(ps_trimmed, ps_mod)
is_refit(ps_refitted)
#> [1] TRUE
```

Then pass the refitted scores to the weight function as usual:

``` r

wts_trimmed <- wt_ate(ps_refitted, dat$z)
#> ℹ Treating `.exposure` as binary
summary(wts_trimmed)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.     NAs 
#>   1.073   1.386   1.726   1.970   2.157   4.724       2
```

See
[`?ps_trim`](https://r-causal.github.io/propensity/reference/ps_trim.md)
for other trimming methods, including percentile-based (`"pctl"`),
preference score (`"pref"`), and common range (`"cr"`).

### Truncation

Truncation is similar to trimming but keeps all observations – it just
clips extreme scores to specified bounds:

``` r

ps_truncated <- ps_trunc(ps, lower = 0.05, upper = 0.95)
```

[`is_unit_truncated()`](https://r-causal.github.io/propensity/reference/is_unit_truncated.md)
tells you which observations were clipped:

``` r

is_ps_truncated(ps_truncated)
#> [1] TRUE
sum(is_unit_truncated(ps_truncated))
#> [1] 0
ps_trunc_meta(ps_truncated)
#> $method
#> [1] "ps"
#> 
#> $lower_bound
#> [1] 0.05
#> 
#> $upper_bound
#> [1] 0.95
#> 
#> $truncated_idx
#> integer(0)
#> 
#> $n_obs
#> [1] 100
```

``` r

wts_truncated <- wt_ate(ps_truncated, dat$z)
#> ℹ Treating `.exposure` as binary
summary(wts_truncated)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.092   1.440   1.780   2.028   2.111  12.583
```

### Which approach?

These aren’t mutually exclusive. In general: overlap estimands like
[`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
are the easiest path if your research question allows it. Trimming
(followed by
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md))
is the standard choice when you need ATE but have near-violations of
positivity. Truncation is a lighter touch when you want to keep the full
sample.

## Interpreting results

### Binary outcomes

For binary outcomes,
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
returns the marginal risk under each exposure level, then three effect
measures built from those risks: the risk difference, log risk ratio,
and log odds ratio. The `contrast` column names the level a mean row
belongs to and the pair of levels an effect measure compares:

``` r

result
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z, family = binomial(), data = dat, weights = wts) 
#> 
#> Marginal estimates:
#>                estimate  std.err      z ci.lower ci.upper conf.level   p.value    
#> mean 0         0.321135 0.069594 4.6144  0.18473  0.45754       0.95 3.942e-06 ***
#> mean 1         0.641133 0.075977 8.4385  0.49222  0.79004       0.95 < 2.2e-16 ***
#> rd 1 vs 0      0.319997 0.103584 3.0892  0.11698  0.52302       0.95  0.002007 ** 
#> log(rr) 1 vs 0 0.691374 0.248114 2.7865  0.20508  1.17767       0.95  0.005328 ** 
#> log(or) 1 vs 0 1.328843 0.461758 2.8778  0.42381  2.23387       0.95  0.004005 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) pulls the
estimates into a data frame:

``` r

as.data.frame(result)
#>      term contrast  estimate  std.error statistic      p.value
#> 1    mean        0 0.3211354 0.06959360  4.614439 3.941578e-06
#> 2    mean        1 0.6411327 0.07597704  8.438507 3.214169e-17
#> 3      rd   1 vs 0 0.3199973 0.10358436  3.089243 2.006669e-03
#> 4 log(rr)   1 vs 0 0.6913736 0.24811365  2.786520 5.327738e-03
#> 5 log(or)   1 vs 0 1.3288426 0.46175800  2.877790 4.004714e-03
```

Use `exponentiate = TRUE` to get risk ratios and odds ratios on their
natural scale. The standard errors, z-statistics, and p-values stay on
the log scale:

``` r

as.data.frame(result, exponentiate = TRUE)
#>   term contrast  estimate  std.error statistic      p.value
#> 1 mean        0 0.3211354 0.06959360  4.614439 3.941578e-06
#> 2 mean        1 0.6411327 0.07597704  8.438507 3.214169e-17
#> 3   rd   1 vs 0 0.3199973 0.10358436  3.089243 2.006669e-03
#> 4   rr   1 vs 0 1.9964559 0.24811365  2.786520 5.327738e-03
#> 5   or   1 vs 0 3.7766699 0.46175800  2.877790 4.004714e-03
```

### Continuous outcomes

For continuous outcomes,
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
returns the marginal mean under each exposure level and the difference
between them. Use [`lm()`](https://rdrr.io/r/stats/lm.html) for the
outcome model:

``` r

y_cont <- 2 + 0.8 * z + 0.3 * x1 + 0.2 * x2 + rnorm(n)
dat$y_cont <- y_cont
outcome_cont <- lm(y_cont ~ z, data = dat, weights = wts)
ipw(ps_mod, outcome_cont)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x1 + x2, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: lm(formula = y_cont ~ z, data = dat, weights = wts) 
#> 
#> Marginal estimates:
#>             estimate std.err      z ci.lower ci.upper conf.level   p.value    
#> mean 0       1.90947 0.13958 13.681  1.63591   2.1830       0.95 < 2.2e-16 ***
#> mean 1       2.83684 0.16397 17.301  2.51546   3.1582       0.95 < 2.2e-16 ***
#> diff 1 vs 0  0.92737 0.20395  4.547  0.52763   1.3271       0.95 5.443e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Next steps

The examples above all use binary exposures. propensity also handles
continuous and categorical treatments.

### Continuous exposures

For continuous exposures, weights use density ratios. Stabilization is
usually a good idea here:

``` r

# Fit a model for the continuous exposure
ps_cont <- glm(continuous_exposure ~ x1 + x2, data = dat, family = gaussian())

# Stabilized weights (strongly recommended for continuous exposures)
wts_cont <- wt_ate(ps_cont, stabilize = TRUE)
```

### Categorical exposures

For multi-level treatments, pass a matrix or data frame of predicted
probabilities with one column per level:

``` r

# Multinomial propensity scores (one column per treatment level)
ps_matrix <- predict(multinom_model, type = "probs")
wt_ate(ps_matrix, exposure, exposure_type = "categorical")

# ATT and ATU require specifying the focal level
wt_att(ps_matrix, exposure, .focal_level = "treated")
```

### Calibration

[`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
adjusts propensity scores so they better reflect treatment
probabilities. Where trimming and truncation deal with the tails,
calibration reshapes the whole distribution. It supports logistic
calibration (the default) and isotonic regression:

``` r

ps_calibrated <- ps_calibrate(ps, dat$z, method = "logistic", smooth = FALSE)
is_ps_calibrated(ps_calibrated)

wts_calibrated <- wt_ate(ps_calibrated, dat$z)
```

### Censoring weights

[`wt_cens()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
calculates inverse probability of censoring weights for survival or
longitudinal analyses:

``` r

# Model the probability of being uncensored
cens_mod <- glm(uncensored ~ x1 + x2, data = dat, family = binomial())
wts_cens <- wt_cens(cens_mod)

# Censoring weights use the same formula as ATE weights
estimand(wts_cens) # "uncensored"
```

### Learning more

See the function reference for details:

- [`?wt_ate`](https://r-causal.github.io/propensity/reference/wt_ate.md)
  – Weight calculation for all estimands
- [`?ps_trim`](https://r-causal.github.io/propensity/reference/ps_trim.md),
  [`?ps_trunc`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
  [`?ps_calibrate`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
  – Handling extreme propensity scores
- [`?ipw`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  – Inverse probability weighted estimation
