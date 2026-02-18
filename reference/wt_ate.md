# Calculate Propensity Score Weights for Causal Inference

This family of functions computes propensity score weights for various
causal estimands:

- **ATE** (Average Treatment Effect)

- **ATT** (Average Treatment Effect on the Treated)

- **ATU** (Average Treatment Effect on the Untreated, sometimes called
  the **ATC**, where the "C" stands for "control"). `wt_atc()` is
  provided as an alias for `wt_atu()`

- **ATM** (Average Treatment Effect for the Evenly Matchable)

- **ATO** (Average Treatment Effect for the Overlap population)

- **Entropy** (Average Treatment Effect for the Entropy-weighted
  population)

- **Censoring weights** can be calculated using `wt_cens()`, which uses
  the same formula as ATE weights but with estimand "uncensored". These
  are useful for handling censoring in survival analysis

  The propensity score can be provided as a numeric vector of predicted
  probabilities, as a `data.frame` where each column represents the
  predicted probability for a level of the exposure, or as a fitted GLM
  object. They can also be propensity score objects created by
  [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
  [`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md),
  or
  [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)

  The returned weights are encapsulated in a `psw` object, which is a
  numeric vector with additional attributes that record the estimand,
  and whether the weights have been stabilized, trimmed, or truncated.

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

  Either a numeric vector of predicted probabilities, a `data.frame`
  where each column corresponds to a level of the exposure, or a fitted
  GLM object. For data frames, the second column is used by default for
  binary exposures unless specified otherwise with `.propensity_col`.
  For GLM objects, fitted values are extracted automatically.

- .exposure:

  The exposure variable. For binary exposures, a vector of 0s and 1s;
  for continuous exposures, a numeric vector. When `.propensity` is a
  GLM object, this argument is optional and will be extracted from the
  model if not provided.

- .sigma:

  For continuous exposures, a numeric vector of standard errors used
  with [`dnorm()`](https://rdrr.io/r/stats/Normal.html). For example,
  this can be derived from the influence measures of a model (e.g.,
  `influence(model)$sigma`).

- exposure_type:

  Character string specifying the type of exposure. Options are
  `"auto"`, `"binary"`, `"categorical"`, and `"continuous"`. Defaults to
  `"auto"`, which detects the type automatically.

- .focal_level:

  For binary exposures, the value representing the focal group
  (typically the treatment group). For categorical exposures with ATT or
  ATU estimands, specifies the focal category. Must be one of the levels
  of the exposure variable. Required for `wt_att()` and `wt_atu()` with
  categorical exposures.

- .reference_level:

  For binary exposures, the value representing the reference group
  (typically the control group). If not provided, it is automatically
  detected.

- stabilize:

  Logical indicating whether to stabilize the weights. For ATE weights,
  stabilization multiplies the weight by either the mean of `.exposure`
  or the supplied `stabilization_score`. Note: stabilization is only
  supported for ATE and continuous exposures.

- stabilization_score:

  Optional numeric value for stabilizing the weights (e.g., a predicted
  value from a regression model without predictors). Only used when
  `stabilize` is `TRUE`.

- ...:

  Reserved for future expansion. Not currently used.

- .treated:

  **\[deprecated\]** Use `.focal_level` instead.

- .untreated:

  **\[deprecated\]** Use `.reference_level` instead.

- .propensity_col:

  With a binary exposure, when `.propensity` is a data frame, specifies
  which column to use for propensity scores. Can be a column name
  (quoted or unquoted) or a numeric index. Defaults to the second column
  if available, otherwise the first. For categorical exposures, the
  entire data frame is used as a matrix of propensity scores.

## Value

A `psw` object (a numeric vector) with additional attributes:

- **estimand**: A description of the estimand (e.g., "ate", "att").

- **stabilized**: A logical flag indicating if stabilization was
  applied.

- **trimmed**: A logical flag indicating if the weights are based on
  trimmed propensity scores.

- **truncated**: A logical flag indicating if the weights are based on
  truncated propensity scores.

## Details

### Theoretical Background

Propensity score weighting is a method for estimating causal effects by
creating a pseudo-population where the exposure is independent of
measured confounders. The propensity score, \\e(X)\\, is the probability
of receiving treatment given observed covariates \\X\\. By weighting
observations inversely proportional to their propensity scores, we can
balance the distribution of covariates between treatment groups. Other
weights allow for different target populations.

### Mathematical Formulas

#### Binary Exposures

For binary treatments (\\A = 0\\ or \\1\\), the weights are:

- **ATE**: \\w = \frac{A}{e(X)} + \frac{1-A}{1-e(X)}\\

- **ATT**: \\w = A + \frac{(1-A) \cdot e(X)}{1-e(X)}\\

- **ATU**: \\w = \frac{A \cdot (1-e(X))}{e(X)} + (1-A)\\

- **ATM**: \\w = \frac{\min(e(X), 1-e(X))}{A \cdot e(X) + (1-A) \cdot
  (1-e(X))}\\

- **ATO**: \\w = A \cdot (1-e(X)) + (1-A) \cdot e(X)\\

- **Entropy**: \\w = \frac{h(e(X))}{A \cdot e(X) + (1-A) \cdot
  (1-e(X))}\\, where \\h(e) = -\[e \cdot \log(e) + (1-e) \cdot
  \log(1-e)\]\\

#### Continuous Exposures

For continuous treatments, weights use the density ratio: \\w =
\frac{f_A(A)}{f\_{A\|X}(A\|X)}\\, where \\f_A\\ is the marginal density
of \\A\\ and \\f\_{A\|X}\\ is the conditional density given \\X\\.

#### Categorical Exposures

For categorical treatments with \\K\\ levels, weights use a tilting
function approach: \\w_i = \frac{h(e_i)}{e\_{i,Z_i}}\\, where
\\e\_{i,Z_i}\\ is the propensity score for unit \\i\\'s observed
treatment level, and \\h(e_i)\\ is a tilting function that depends on
the estimand:

- **ATE**: \\h(e) = 1\\

- **ATT**: \\h(e) = e\_{focal}\\ (propensity score for the focal
  category)

- **ATU**: \\h(e) = 1 - e\_{focal}\\ (complement of focal category
  propensity)

- **ATM**: \\h(e) = \min(e_1, ..., e_K)\\

- **ATO**: \\h(e) = 1 / \sum_k(1/e_k)\\ (reciprocal of harmonic mean
  denominator)

- **Entropy**: \\h(e) = -\sum_k\[e_k \cdot \log(e_k)\]\\ (entropy of
  propensity scores)

### Exposure Types

The functions support different types of exposures:

- **`binary`**: For dichotomous treatments (e.g. 0/1).

- **`continuous`**: For numeric exposures. Here, weights are calculated
  via the normal density using
  [`dnorm()`](https://rdrr.io/r/stats/Normal.html).

- **`categorical`**: For exposures with more than 2 categories. Requires
  `.propensity` to be a matrix or data frame with columns representing
  propensity scores for each category.

- **`auto`**: Automatically detects the exposure type based on
  `.exposure`.

### Stabilization

For ATE weights, stabilization can improve the performance of the
estimator by reducing variance. When `stabilize` is `TRUE` and no
`stabilization_score` is provided, the weights are multiplied by the
mean of `.exposure`. Alternatively, if a `stabilization_score` is
provided, it is used as the multiplier. Stabilized weights have the
form: \\w_s = f_A(A) \times w\\, where \\f_A(A)\\ is the marginal
probability or density.

### Weight Properties and Diagnostics

Extreme weights can indicate:

- Positivity violations (near 0 or 1 propensity scores)

- Poor model specification

- Lack of overlap between treatment groups

See the halfmoon package for tools to diagnose and visualize weights.

You can address extreme weights in several ways. The first is to modify
the target population: use trimming, truncation, or alternative
estimands (ATM, ATO, entropy). Another technique that can help is
stabilization, which reduces variance of the weights.

### Trimmed and Truncated Weights

In addition to the standard weight functions, versions exist for trimmed
and truncated propensity score weights created by
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
and
[`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md).
These variants calculate the weights using modified propensity scores
(trimmed or truncated) and update the estimand attribute accordingly.

## References

For detailed guidance on causal inference in R, see [*Causal Inference
in R*](https://www.r-causal.org/) by Malcolm Barrett, Lucy D'Agostino
McGowan, and Travis Gerke.

### Foundational Papers

Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the
propensity score in observational studies for causal effects.
*Biometrika*, 70(1), 41-55.

### Estimand-Specific Methods

Li, L., & Greene, T. (2013). A weighting analogue to pair matching in
propensity score analysis. *The International Journal of Biostatistics*,
9(2), 215-234. (ATM weights)

Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates
via propensity score weighting. *Journal of the American Statistical
Association*, 113(521), 390-400. (ATO weights)

Zhou, Y., Matsouaka, R. A., & Thomas, L. (2020). Propensity score
weighting under limited overlap and model misspecification. *Statistical
Methods in Medical Research*, 29(12), 3721-3756. (Entropy weights)

### Continuous Exposures

Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous
treatments. *Applied Bayesian Modeling and Causal Inference from
Incomplete-Data Perspectives*, 226164, 73-84.

### Practical Guidance

Austin, P. C., & Stuart, E. A. (2015). Moving towards best practice when
using inverse probability of treatment weighting (IPTW) using the
propensity score to estimate causal treatment effects in observational
studies. *Statistics in Medicine*, 34(28), 3661-3679.

## See also

- [`psw()`](https://r-causal.github.io/propensity/reference/psw.md) for
  details on the structure of the returned weight objects.

- [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
  [`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md),
  and
  [`ps_refit()`](https://r-causal.github.io/propensity/reference/ps_refit.md)
  for handling extreme weights.

- [`ps_calibrate()`](https://r-causal.github.io/propensity/reference/ps_calibrate.md)
  for calibrating weights.

## Examples

``` r
## Basic Usage with Binary Exposures

# Simulate a simple dataset
set.seed(123)
n <- 100
propensity_scores <- runif(n, 0.1, 0.9)
treatment <- rbinom(n, 1, propensity_scores)

# Calculate different weight types
weights_ate <- wt_ate(propensity_scores, treatment)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
weights_att <- wt_att(propensity_scores, treatment)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
weights_atu <- wt_atu(propensity_scores, treatment)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1

# With explicit focal and reference levels
weights_att_explicit <- wt_att(propensity_scores, treatment,
                               .focal_level = 1, .reference_level = 0)
#> ℹ Treating `.exposure` as binary
weights_atm <- wt_atm(propensity_scores, treatment)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
weights_ato <- wt_ato(propensity_scores, treatment)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
weights_entropy <- wt_entropy(propensity_scores, treatment)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1

# Compare weight distributions
summary(weights_ate)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.112   1.317   1.591   2.044   2.047   7.482 
summary(weights_ato)  # Often more stable than ATE
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1005  0.2406  0.3713  0.3976  0.5113  0.8664 

## Stabilized Weights

# Stabilization reduces variance
weights_ate_stab <- wt_ate(propensity_scores, treatment, stabilize = TRUE)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1

## Handling Extreme Propensity Scores

# Create data with positivity violations
ps_extreme <- c(0.01, 0.02, 0.98, 0.99, rep(0.5, 4))
trt_extreme <- c(0, 0, 1, 1, 0, 1, 0, 1)

# Standard ATE weights can be extreme
wt_extreme <- wt_ate(ps_extreme, trt_extreme)
#> ℹ Treating `.exposure` as binary
# Very large!
max(wt_extreme)
#> [1] 2

# ATO weights are bounded
wt_extreme_ato <- wt_ato(ps_extreme, trt_extreme)
#> ℹ Treating `.exposure` as binary
# Much more reasonable
max(wt_extreme_ato)
#> [1] 0.5
# but they target a different population
estimand(wt_extreme_ato) # "ato"
#> [1] "ato"

## Working with Data Frames

# Example with custom data frame
ps_df <- data.frame(
  control = c(0.9, 0.7, 0.3, 0.1),
  treated = c(0.1, 0.3, 0.7, 0.9)
)
exposure <- c(0, 0, 1, 1)

# Uses second column by default (treated probabilities)
wt_ate(ps_df, exposure)
#> ℹ Treating `.exposure` as binary
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[4]>
#> [1] 1.111111 1.428571 1.428571 1.111111

# Explicitly specify column by name
wt_ate(ps_df, exposure, .propensity_col = "treated")
#> ℹ Treating `.exposure` as binary
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[4]>
#> [1] 1.111111 1.428571 1.428571 1.111111

# Or by position
wt_ate(ps_df, exposure, .propensity_col = 2)
#> ℹ Treating `.exposure` as binary
#> ℹ Treating `.exposure` as binary
#> <psw{estimand = ate}[4]>
#> [1] 1.111111 1.428571 1.428571 1.111111

## Working with GLM Objects

# Fit a propensity score model
set.seed(123)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n)
treatment <- rbinom(n, 1, plogis(0.5 * x1 + 0.3 * x2))

ps_model <- glm(treatment ~ x1 + x2, family = binomial)

# Use GLM directly for weight calculation
weights_from_glm <- wt_ate(ps_model, treatment)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1

# Or omit the exposure argument (it will be extracted from the GLM)
weights_from_glm_auto <- wt_ate(ps_model)
#> ℹ Using exposure variable "treatment" from GLM model
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
```
