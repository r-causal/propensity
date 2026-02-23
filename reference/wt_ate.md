# Calculate propensity score weights

Compute inverse probability weights for causal inference under different
estimands. Each function targets a different population:

- `wt_ate()`: **Average Treatment Effect** – the full population.

- `wt_att()`: **Average Treatment Effect on the Treated** – the treated
  (focal) group.

- `wt_atu()`: **Average Treatment Effect on the Untreated** – the
  untreated (reference) group. `wt_atc()` is an alias.

- `wt_atm()`: **Average Treatment Effect for the Evenly Matchable** –
  units with the most overlap.

- `wt_ato()`: **Average Treatment Effect for the Overlap Population** –
  weights proportional to overlap.

- `wt_entropy()`: **Entropy-weighted Average Treatment Effect** – an
  entropy-balanced population.

- `wt_cens()`: **Inverse probability of censoring weights** – uses the
  same formula as `wt_ate()` but labels the estimand `"uncensored"`. Use
  these to adjust for censoring in survival analysis, not for treatment
  weighting.

`.propensity` accepts a numeric vector of predicted probabilities, a
`data.frame` of per-level probabilities, a fitted `glm` object, or a
modified propensity score created by
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md),
or
[`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md).

All functions return a
[`psw`](https://r-causal.github.io/propensity/reference/psw.md) object –
a numeric vector that tracks the estimand, stabilization status, and any
trimming or truncation applied.

## Usage

``` r
wt_ate(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_ate(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_att(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_att(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_atu(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_atu(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_atm(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_atm(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_ato(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_ato(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_entropy(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_entropy(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)

wt_atc(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical"),
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

wt_cens(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
)

# S3 method for class 'data.frame'
wt_cens(
  .propensity,
  .exposure,
  .sigma = NULL,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .focal_level = NULL,
  .reference_level = NULL,
  stabilize = FALSE,
  stabilization_score = NULL,
  ...,
  .propensity_col = NULL,
  .treated = NULL,
  .untreated = NULL
)
```

## Arguments

- .propensity:

  Propensity scores in one of several forms:

  - A **numeric vector** of predicted probabilities (binary/continuous).

  - A **data frame** or matrix with one column per exposure level
    (categorical), or two columns for binary (see `.propensity_col`).

  - A fitted **`glm`** object – fitted values are extracted
    automatically.

  - A modified propensity score created by
    [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
    [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
    [`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md),
    or
    [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md).

- .exposure:

  The exposure (treatment) variable. For binary exposures, a numeric 0/1
  vector, logical, or two-level factor. For categorical exposures, a
  factor or character vector. For continuous exposures, a numeric
  vector. Optional when `.propensity` is a `glm` object (extracted from
  the model).

- .sigma:

  Numeric vector of observation-level standard deviations for continuous
  exposures (e.g., `influence(model)$sigma`). Extracted automatically
  when `.propensity` is a `glm` object.

- exposure_type:

  Type of exposure: `"auto"` (default), `"binary"`, `"categorical"`, or
  `"continuous"`. `"auto"` detects the type from `.exposure`.

- .focal_level:

  The value of `.exposure` representing the focal (treated) group. For
  binary exposures, defaults to the higher value. Required for
  `wt_att()` and `wt_atu()` with categorical exposures.

- .reference_level:

  The value of `.exposure` representing the reference (control) group.
  Automatically detected if not supplied.

- stabilize:

  If `TRUE`, multiply weights by an estimate of the marginal treatment
  probability (binary) or density (continuous). Only supported by
  `wt_ate()` and `wt_cens()`. See **Stabilization** in Details.

- stabilization_score:

  Optional numeric value to use as the stabilization multiplier instead
  of the default (the marginal mean of `.exposure`). Ignored when
  `stabilize = FALSE`.

- ...:

  These dots are for future extensions and must be empty.

- .treated:

  **\[deprecated\]** Use `.focal_level` instead.

- .untreated:

  **\[deprecated\]** Use `.reference_level` instead.

- .propensity_col:

  Column to use when `.propensity` is a data frame with a binary
  exposure. Accepts a column name (quoted or unquoted) or numeric index.
  Defaults to the second column. Ignored for categorical exposures,
  where all columns are used.

## Value

A [`psw`](https://r-causal.github.io/propensity/reference/psw.md) vector
(a double vector with class `psw`) carrying these attributes:

- `estimand`: character, e.g. `"ate"`, `"att"`, `"uncensored"`.

- `stabilized`: logical, whether stabilization was applied.

- `trimmed`: logical, whether the propensity scores were trimmed.

- `truncated`: logical, whether the propensity scores were truncated.

- `calibrated`: logical, whether the propensity scores were calibrated.

## Details

### Exposure types

All weight functions support binary exposures. `wt_ate()` and
`wt_cens()` also support continuous exposures. All except `wt_cens()`
support categorical exposures.

- **Binary**: `.exposure` is a two-level vector (e.g., 0/1, logical, or
  a two-level factor). `.propensity` is a numeric vector of P(treatment
  \| X).

- **Categorical**: `.exposure` is a factor or character vector with 3+
  levels. `.propensity` must be a matrix or data frame with one column
  per level, where rows sum to 1.

- **Continuous**: `.exposure` is a numeric vector. `.propensity` is a
  vector of conditional means (fitted values). Weights use a normal
  density ratio; stabilization is strongly recommended.

- **Auto** (default): Detects the exposure type from `.exposure`.

### Stabilization

Setting `stabilize = TRUE` multiplies the base weight by an estimate of
P(A) (binary) or f_A(A) (continuous), reducing variance. When no
`stabilization_score` is supplied, the marginal mean of `.exposure` is
used. Stabilization is supported for ATE and censoring weights
(`wt_ate()` and `wt_cens()`) and is strongly recommended for continuous
exposures.

### Handling extreme weights

Extreme weights signal positivity violations, poor model fit, or limited
overlap. You can address them by:

- Choosing an overlap-focused estimand (`wt_ato()`, `wt_atm()`,
  `wt_entropy()`), which down-weight units in regions of poor overlap.

- Trimming
  ([`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md))
  or truncating
  ([`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md))
  propensity scores before computing weights.

- Calibrating weights with
  [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md).

- Stabilizing ATE weights (`stabilize = TRUE`).

See the [halfmoon](https://CRAN.R-project.org/package=halfmoon) package
for weight diagnostics and visualization.

### Weight formulas

#### Binary exposures

For binary treatments (\\A \in \\0, 1\\\\), with propensity score \\e(X)
= P(A=1 \mid X)\\:

- **ATE**: \\w = \frac{A}{e(X)} + \frac{1-A}{1-e(X)}\\

- **ATT**: \\w = A + \frac{(1-A) \cdot e(X)}{1-e(X)}\\

- **ATU**: \\w = \frac{A \cdot (1-e(X))}{e(X)} + (1-A)\\

- **ATM**: \\w = \frac{\min(e(X), 1-e(X))}{A \cdot e(X) + (1-A) \cdot
  (1-e(X))}\\

- **ATO**: \\w = A \cdot (1-e(X)) + (1-A) \cdot e(X)\\

- **Entropy**: \\w = \frac{h(e(X))}{A \cdot e(X) + (1-A) \cdot
  (1-e(X))}\\, where \\h(e) = -\[e \log(e) + (1-e) \log(1-e)\]\\

#### Continuous exposures

Weights use the density ratio \\w = f_A(A) / f\_{A\|X}(A \mid X)\\,
where \\f_A\\ is the marginal density and \\f\_{A\|X}\\ is the
conditional density (both assumed normal). Only `wt_ate()` and
`wt_cens()` support continuous exposures.

#### Categorical exposures

For \\K\\-level treatments, weights take the tilting-function form \\w_i
= h(\mathbf{e}\_i) / e\_{i,Z_i}\\, where \\e\_{i,Z_i}\\ is the
propensity for unit \\i\\'s observed level and \\h(\cdot)\\ depends on
the estimand:

- **ATE**: \\h(\mathbf{e}) = 1\\

- **ATT**: \\h(\mathbf{e}) = e\_{\text{focal}}\\

- **ATU**: \\h(\mathbf{e}) = 1 - e\_{\text{focal}}\\

- **ATM**: \\h(\mathbf{e}) = \min(e_1, \ldots, e_K)\\

- **ATO**: \\h(\mathbf{e}) = \bigl(\sum_k 1/e_k\bigr)^{-1}\\

- **Entropy**: \\h(\mathbf{e}) = -\sum_k e_k \log(e_k)\\

## References

Barrett, M., D'Agostino McGowan, L., & Gerke, T. *Causal Inference in
R*. <https://www.r-causal.org/>

Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the
propensity score in observational studies for causal effects.
*Biometrika*, 70(1), 41–55.

Li, L., & Greene, T. (2013). A weighting analogue to pair matching in
propensity score analysis. *The International Journal of Biostatistics*,
9(2), 215–234. (ATM weights)

Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates
via propensity score weighting. *Journal of the American Statistical
Association*, 113(521), 390–400. (ATO weights)

Zhou, Y., Matsouaka, R. A., & Thomas, L. (2020). Propensity score
weighting under limited overlap and model misspecification. *Statistical
Methods in Medical Research*, 29(12), 3721–3756. (Entropy weights)

Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous
treatments. In *Applied Bayesian Modeling and Causal Inference from
Incomplete-Data Perspectives* (pp. 73–84).

Austin, P. C., & Stuart, E. A. (2015). Moving towards best practice when
using inverse probability of treatment weighting (IPTW). *Statistics in
Medicine*, 34(28), 3661–3679.

## See also

- [`psw()`](https://r-causal.github.io/propensity/reference/psw.md) for
  the returned weight vector class.

- [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
  [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
  [`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md),
  and
  [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
  for modifying propensity scores before weighting.

- [`ipw()`](https://r-causal.github.io/propensity/reference/ipw.md) for
  inverse-probability-weighted estimation of causal effects.

## Examples

``` r
# -- Binary exposure, numeric propensity scores ----------------------
set.seed(123)
ps <- runif(100, 0.1, 0.9)
trt <- rbinom(100, 1, ps)

wt_ate(ps, trt)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = ate}[100]>
#>   [1] 1.492675 1.368655 1.745754 5.165661 1.173194 7.328950 2.094172 1.228599
#>   [9] 1.847923 1.870179 7.433103 1.861044 1.557495 2.262990 1.223002 1.219720
#>  [17] 1.422212 7.482363 1.568225 1.157940 1.232086 1.528485 1.632905 1.116800
#>  [25] 1.601115 3.001420 1.868276 1.738182 1.495501 1.278267 1.148871 5.612908
#>  [33] 2.878230 3.793252 1.135965 2.073670 3.410265 3.661309 2.820518 1.399190
#>  [41] 1.272653 1.759439 1.757406 1.653101 4.505402 1.267499 1.401399 1.896705
#>  [49] 1.455134 1.271840 1.158299 1.830697 1.352924 1.246136 1.822296 1.360961
#>  [57] 1.253173 1.423191 1.225436 1.665474 1.582048 1.213405 2.455942 1.469523
#>  [65] 1.330297 1.847790 1.336806 1.333490 1.359668 1.824369 1.421302 1.657339
#>  [73] 3.013373 1.111729 2.082235 1.381397 1.677439 1.694293 2.621656 1.232906
#>  [81] 3.391031 1.576182 2.303524 1.368819 1.222930 1.811313 1.126170 1.227836
#>  [89] 5.240410 4.165936 1.257160 1.606473 2.667996 1.598960 2.806635 1.333605
#>  [97] 1.377723 1.211939 1.899058 2.037508
wt_att(ps, trt)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = att}[100]>
#>   [1] 0.4926755 1.0000000 0.7457538 4.1656608 1.0000000 1.0000000 1.0941724
#>   [8] 1.0000000 1.0000000 0.8701789 6.4331026 0.8610445 1.0000000 1.2629898
#>  [15] 0.2230018 1.0000000 0.4222125 1.0000000 0.5682254 1.0000000 1.0000000
#>  [22] 1.0000000 1.0000000 1.0000000 1.0000000 2.0014200 1.0000000 1.0000000
#>  [29] 0.4955011 0.2782671 1.0000000 4.6129081 1.8782298 2.7932516 0.1359647
#>  [36] 1.0000000 2.4102647 1.0000000 1.0000000 0.3991897 0.2726533 0.7594392
#>  [43] 0.7574058 0.6531012 1.0000000 0.2674992 0.4013989 0.8967053 0.4551341
#>  [50] 1.0000000 0.1582988 0.8306973 1.0000000 0.2461361 1.0000000 0.3609610
#>  [57] 0.2531726 1.0000000 1.0000000 0.6654737 1.0000000 0.2134045 1.0000000
#>  [64] 0.4695226 1.0000000 0.8477904 1.0000000 1.0000000 1.0000000 0.8243692
#>  [71] 1.0000000 1.0000000 2.0133726 0.1117285 1.0000000 0.3813969 0.6774393
#>  [78] 1.0000000 1.0000000 0.2329063 1.0000000 1.0000000 1.0000000 1.0000000
#>  [85] 0.2229300 0.8113126 1.0000000 1.0000000 4.2404103 1.0000000 0.2571604
#>  [92] 1.0000000 1.0000000 1.0000000 1.0000000 0.3336052 1.0000000 0.2119390
#>  [99] 0.8990583 1.0375079
wt_atu(ps, trt)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = atu}[100]>
#>   [1] 1.0000000 0.3686554 1.0000000 1.0000000 0.1731942 6.3289497 1.0000000
#>   [8] 0.2285990 0.8479233 1.0000000 1.0000000 1.0000000 0.5574953 1.0000000
#>  [15] 1.0000000 0.2197205 1.0000000 6.4823626 1.0000000 0.1579396 0.2320863
#>  [22] 0.5284847 0.6329051 0.1167996 0.6011153 1.0000000 0.8682760 0.7381824
#>  [29] 1.0000000 1.0000000 0.1488715 1.0000000 1.0000000 1.0000000 1.0000000
#>  [36] 1.0736701 1.0000000 2.6613092 1.8205180 1.0000000 1.0000000 1.0000000
#>  [43] 1.0000000 1.0000000 3.5054016 1.0000000 1.0000000 1.0000000 1.0000000
#>  [50] 0.2718404 1.0000000 1.0000000 0.3529239 1.0000000 0.8222956 1.0000000
#>  [57] 1.0000000 0.4231912 0.2254357 1.0000000 0.5820478 1.0000000 1.4559422
#>  [64] 1.0000000 0.3302967 1.0000000 0.3368064 0.3334905 0.3596676 1.0000000
#>  [71] 0.4213022 0.6573389 1.0000000 1.0000000 1.0822347 1.0000000 1.0000000
#>  [78] 0.6942927 1.6216558 1.0000000 2.3910308 0.5761821 1.3035242 0.3688192
#>  [85] 1.0000000 1.0000000 0.1261698 0.2278362 1.0000000 3.1659355 1.0000000
#>  [92] 0.6064733 1.6679958 0.5989600 1.8066347 1.0000000 0.3777228 1.0000000
#>  [99] 1.0000000 1.0000000
wt_atm(ps, trt)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = atm}[100]>
#>   [1] 0.4926755 0.3686554 0.7457538 1.0000000 0.1731942 1.0000000 1.0000000
#>   [8] 0.2285990 0.8479233 0.8701789 1.0000000 0.8610445 0.5574953 1.0000000
#>  [15] 0.2230018 0.2197205 0.4222125 1.0000000 0.5682254 0.1579396 0.2320863
#>  [22] 0.5284847 0.6329051 0.1167996 0.6011153 1.0000000 0.8682760 0.7381824
#>  [29] 0.4955011 0.2782671 0.1488715 1.0000000 1.0000000 1.0000000 0.1359647
#>  [36] 1.0000000 1.0000000 1.0000000 1.0000000 0.3991897 0.2726533 0.7594392
#>  [43] 0.7574058 0.6531012 1.0000000 0.2674992 0.4013989 0.8967053 0.4551341
#>  [50] 0.2718404 0.1582988 0.8306973 0.3529239 0.2461361 0.8222956 0.3609610
#>  [57] 0.2531726 0.4231912 0.2254357 0.6654737 0.5820478 0.2134045 1.0000000
#>  [64] 0.4695226 0.3302967 0.8477904 0.3368064 0.3334905 0.3596676 0.8243692
#>  [71] 0.4213022 0.6573389 1.0000000 0.1117285 1.0000000 0.3813969 0.6774393
#>  [78] 0.6942927 1.0000000 0.2329063 1.0000000 0.5761821 1.0000000 0.3688192
#>  [85] 0.2229300 0.8113126 0.1261698 0.2278362 1.0000000 1.0000000 0.2571604
#>  [92] 0.6064733 1.0000000 0.5989600 1.0000000 0.3336052 0.3777228 0.2119390
#>  [99] 0.8990583 1.0000000
wt_ato(ps, trt)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = ato}[100]>
#>   [1] 0.3300620 0.2693559 0.4271815 0.8064139 0.1476262 0.8635548 0.5224844
#>   [8] 0.1860648 0.4588520 0.4652918 0.8654667 0.4626673 0.3579435 0.5581067
#>  [15] 0.1823397 0.1801400 0.2968702 0.8663524 0.3623366 0.1363971 0.1883685
#>  [22] 0.3457573 0.3875945 0.1045842 0.3754354 0.6668244 0.4647472 0.4246864
#>  [29] 0.3313278 0.2176909 0.1295806 0.8218392 0.6525642 0.7363739 0.1196909
#>  [36] 0.5177632 0.7067676 0.7268737 0.6454552 0.2853006 0.2142400 0.4316371
#>  [43] 0.4309795 0.3950764 0.7780442 0.2110449 0.2864273 0.4727700 0.3127781
#>  [50] 0.2137378 0.1366649 0.4537601 0.2608601 0.1975194 0.4512416 0.2652251
#>  [57] 0.2020253 0.2973537 0.1839637 0.3995702 0.3679078 0.1758725 0.5928243
#>  [64] 0.3195069 0.2482880 0.4588131 0.2519485 0.2500884 0.2645261 0.4518654
#>  [71] 0.2964199 0.3966231 0.6681459 0.1004998 0.5197467 0.2760951 0.4038532
#>  [78] 0.4097832 0.6185617 0.1889083 0.7051044 0.3655555 0.5658826 0.2694433
#>  [85] 0.1822917 0.4479142 0.1120344 0.1855591 0.8091752 0.7599579 0.2045566
#>  [92] 0.3775185 0.6251868 0.3745935 0.6437014 0.2501529 0.2741646 0.1748760
#>  [99] 0.4734232 0.5092044
wt_entropy(ps, trt)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = entropy}[100]>
#>   [1] 0.9466884 0.7974021 1.1914845 2.5383087 0.4910630 2.9202760 1.4494516
#>   [8] 0.5903003 1.2746181 1.2917997 2.9354392 1.2847853 1.0158386 1.5532690
#>  [15] 0.5808315 0.5752272 0.8649741 2.9425314 1.0267968 0.4612870 0.5961433
#>  [22] 0.9855373 1.0902253 0.3741728 1.0595945 1.9101181 1.2903427 1.1850225
#>  [29] 0.9498151 0.6697735 0.4429918 2.6301699 1.8588898 2.1880067 0.4161138
#>  [36] 1.4360497 2.0632803 2.1467853 1.8339424 0.8365617 0.6611694 1.2030532
#>  [43] 1.2013433 1.1091726 2.3850345 0.6531902 0.8393281 1.3118818 0.9040828
#>  [50] 0.6599161 0.4620024 1.2611029 0.7765170 0.6192602 1.2544407 0.7872501
#>  [57] 0.6305931 0.8661618 0.5849629 1.1205924 1.0407229 0.5643262 1.6597604
#>  [64] 0.9206519 0.7455600 1.2745145 0.7545814 0.7499980 0.7855318 1.2560894
#>  [71] 0.8638679 1.1130997 1.9149498 0.3626234 1.4416708 0.8139569 1.1315051
#>  [78] 1.1466626 1.7427820 0.5975110 2.0565825 1.0348389 1.5766261 0.7976169
#>  [85] 0.5807093 1.2456605 0.3950015 0.5890166 2.5542633 2.2959657 0.6369461
#>  [92] 1.0648287 1.7647932 1.0574806 1.8278454 0.7501570 0.8092153 0.5617751
#>  [99] 1.3136430 1.4119476

# Stabilized ATE weights (reduces variance)
wt_ate(ps, trt, stabilize = TRUE)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = ate; stabilized}[100]>
#>   [1] 0.7612645 0.6706411 0.8903344 2.6344870 0.5748651 3.5911853 1.0680279
#>   [8] 0.6020135 0.9054824 0.9537912 3.7908823 0.9491327 0.7631727 1.1541248
#>  [15] 0.6237309 0.5976630 0.7253284 3.6663577 0.7997950 0.5673904 0.6037223
#>  [22] 0.7489575 0.8001235 0.5472318 0.7845465 1.5307242 0.9154552 0.8517094
#>  [29] 0.7627055 0.6519162 0.5629470 2.8625831 1.4678972 1.9345583 0.5793420
#>  [36] 1.0160984 1.7392350 1.7940415 1.3820538 0.7135867 0.6490532 0.8973140
#>  [43] 0.8962770 0.8430816 2.2076468 0.6464246 0.7147134 0.9673197 0.7421184
#>  [50] 0.6232018 0.5907324 0.9336556 0.6629327 0.6355294 0.8929248 0.6940901
#>  [57] 0.6391180 0.6973637 0.6004635 0.8493916 0.7752034 0.6188363 1.2034117
#>  [64] 0.7494566 0.6518454 0.9423731 0.6550351 0.6534103 0.6662371 0.9304283
#>  [71] 0.6964381 0.8120960 1.5368200 0.5669815 1.0202950 0.7045124 0.8554940
#>  [78] 0.8302034 1.2846113 0.6287822 1.6616051 0.7723292 1.1287269 0.6707214
#>  [85] 0.6236943 0.9237694 0.5518232 0.6016397 2.6726093 2.0413084 0.6411518
#>  [92] 0.7871719 1.3073180 0.7834904 1.3752510 0.6801387 0.6750841 0.6180889
#>  [99] 0.9685198 1.0391291

# Inspect the result
w <- wt_ate(ps, trt)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
estimand(w)
#> [1] "ate"
summary(w)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.112   1.317   1.591   2.044   2.047   7.482 

# -- Overlap-focused estimands handle extreme PS better --------------
ps_extreme <- c(0.01, 0.02, 0.98, 0.99, rep(0.5, 4))
trt_extreme <- c(0, 0, 1, 1, 0, 1, 0, 1)

max(wt_ate(ps_extreme, trt_extreme))
#> ℹ Treating `.exposure` as binary
#> [1] 2
max(wt_ato(ps_extreme, trt_extreme))
#> ℹ Treating `.exposure` as binary
#> [1] 0.5

# -- From a fitted GLM -----------------------------------------------
x1 <- rnorm(100)
x2 <- rnorm(100)
trt2 <- rbinom(100, 1, plogis(0.5 * x1 + 0.3 * x2))
ps_model <- glm(trt2 ~ x1 + x2, family = binomial)

# Exposure is extracted from the model automatically
wt_ate(ps_model)
#> ℹ Using exposure variable "trt2" from GLM model
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = ate}[100]>
#>   [1] 1.772299 2.816364 1.919886 2.059945 1.590514 2.019637 1.592626 3.675866
#>   [9] 1.731064 1.526937 2.194737 1.585679 1.604308 2.035997 1.836701 3.038426
#>  [17] 1.926720 2.527851 2.956648 1.456686 1.951071 2.308830 2.122857 1.859348
#>  [25] 1.333779 1.789153 2.089527 2.008161 1.809376 1.845397 1.236823 2.436497
#>  [33] 1.817852 2.352419 1.300776 1.637110 3.144659 2.924006 4.687625 1.406029
#>  [41] 2.396310 2.163061 2.510037 1.356974 1.430312 2.334509 3.102993 2.202745
#>  [49] 1.310198 1.607962 2.624988 2.535210 1.793952 1.647484 2.569296 2.061068
#>  [57] 3.046296 1.936586 1.517183 1.624578 2.826262 2.825176 1.554039 1.098856
#>  [65] 1.637046 1.507218 2.581387 1.554926 2.385541 1.698365 1.811373 1.670399
#>  [73] 1.760472 1.365580 1.911171 1.539045 2.208360 1.724933 1.617692 1.874702
#>  [81] 1.406166 3.813823 1.982021 2.553269 1.982884 1.957207 1.440986 1.560568
#>  [89] 1.603633 1.885138 1.584638 2.218637 1.630055 1.585865 1.872610 1.274698
#>  [97] 1.393716 1.394134 1.800220 1.778728

# -- Data frame input ------------------------------------------------
ps_df <- data.frame(
  control = c(0.9, 0.7, 0.3, 0.1),
  treated = c(0.1, 0.3, 0.7, 0.9)
)
exposure <- c(0, 0, 1, 1)
wt_ate(ps_df, exposure)
#> ℹ Treating `.exposure` as binary
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[4]>
#> [1] 1.111111 1.428571 1.428571 1.111111
wt_ate(ps_df, exposure, .propensity_col = "treated")
#> ℹ Treating `.exposure` as binary
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[4]>
#> [1] 1.111111 1.428571 1.428571 1.111111

# -- Censoring weights -----------------------------------------------
cens_ps <- runif(50, 0.6, 0.95)
cens_ind <- rbinom(50, 1, cens_ps)
wt_cens(cens_ps, cens_ind)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> <psw{estimand = uncensored}[50]>
#>  [1] 1.116412 1.447081 1.644465 1.108844 1.394494 1.217840 1.263824 3.727943
#>  [9] 1.503634 1.333104 1.218367 1.278290 1.203561 1.169213 1.298047 1.361360
#> [17] 1.660049 1.663092 1.054599 2.760604 3.945306 1.174730 1.162941 1.104729
#> [25] 1.230385 1.156509 1.115406 1.227788 1.139434 2.551015 1.340626 1.103498
#> [33] 1.209817 1.082758 2.943590 1.407391 1.134443 1.056929 1.402849 2.502609
#> [41] 1.055677 1.533923 1.620017 1.235981 3.477837 1.274250 1.138436 1.602240
#> [49] 1.185878 1.072937
estimand(wt_cens(cens_ps, cens_ind))  # "uncensored"
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
#> [1] "uncensored"
```
