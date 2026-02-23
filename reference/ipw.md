# Inverse Probability Weighted Estimation

`ipw()` is a bring-your-own-model (BYOM) inverse probability weighted
estimator for causal inference. You supply a fitted propensity score
model and a fitted weighted outcome model; `ipw()` computes causal
effect estimates with standard errors that correctly account for the
two-step estimation process.

`ipw()` currently supports binary exposures with binary or continuous
outcomes. For binary outcomes, it returns the risk difference, log risk
ratio, and log odds ratio. For continuous outcomes, it returns the
difference in means.

## Usage

``` r
ipw(
  ps_mod,
  outcome_mod,
  .data = NULL,
  estimand = NULL,
  ps_link = NULL,
  conf_level = 0.95
)

# S3 method for class 'ipw'
as.data.frame(x, row.names = NULL, optional = NULL, exponentiate = FALSE, ...)
```

## Arguments

- ps_mod:

  A fitted propensity score model of class
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html), typically a
  logistic regression with the exposure as the left-hand side of the
  formula.

- outcome_mod:

  A fitted weighted outcome model of class
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html), with the outcome as
  the dependent variable and propensity score weights supplied via the
  `weights` argument. The weights should be created with a propensity
  weight function such as
  [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md).

- .data:

  A data frame containing the exposure, outcome, and covariates. If
  `NULL` (the default), `ipw()` extracts data from the model objects.
  Supply `.data` explicitly if the outcome model formula contains
  transformations that prevent extraction of the exposure variable from
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html).

- estimand:

  A character string specifying the causal estimand: `"ate"`, `"att"`,
  `"ato"`, or `"atm"`. If `NULL`, the estimand is inferred from the
  weights in `outcome_mod`. Auto-detection requires weights created with
  [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
  [`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
  [`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
  or
  [`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md).

- ps_link:

  A character string specifying the link function used in the propensity
  score model: `"logit"`, `"probit"`, or `"cloglog"`. Defaults to the
  link used by `ps_mod`.

- conf_level:

  Confidence level for intervals. Default is `0.95`.

- x:

  An `ipw` object.

- row.names, optional, ...:

  Passed to
  [`base::as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- exponentiate:

  If `TRUE`, exponentiate the log risk ratio and log odds ratio to
  produce risk ratios and odds ratios on their natural scale. The
  confidence interval bounds are also exponentiated. Standard errors, z
  statistics, and p-values remain on the log scale. Default is `FALSE`.

## Value

An S3 object of class `ipw` with the following components:

- `estimand`:

  The causal estimand: one of `"ate"`, `"att"`, `"ato"`, or `"atm"`.

- `ps_mod`:

  The fitted propensity score model.

- `outcome_mod`:

  The fitted outcome model.

- `estimates`:

  A data frame with one row per effect measure and the following
  columns: `effect` (the measure name), `estimate` (point estimate),
  `std.err` (standard error), `z` (z-statistic), `ci.lower` and
  `ci.upper` (confidence interval bounds), `conf.level`, and `p.value`.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the `estimates` component as a data frame. When `exponentiate = TRUE`,
the `log(rr)` and `log(or)` rows are transformed: point estimates and
confidence limits are exponentiated and the effect labels become `"rr"`
and `"or"`.

## Workflow

`ipw()` is designed around a three-step workflow:

1.  Fit a propensity score model (e.g., logistic regression of exposure
    on confounders).

2.  Calculate propensity score weights for your estimand (e.g.,
    [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md))
    and fit a weighted outcome model.

3.  Pass both models to `ipw()` to obtain causal effect estimates with
    correct standard errors.

You are responsible for specifying and fitting both models. `ipw()`
handles the variance estimation.

## Effect measures

For binary outcomes ([`stats::glm()`](https://rdrr.io/r/stats/glm.html)
with `family = binomial()`), `ipw()` returns three effect measures:

- `rd`: Risk difference (marginal risk in exposed minus unexposed)

- `log(rr)`: Log risk ratio

- `log(or)`: Log odds ratio

For continuous outcomes
([`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) with
`family = gaussian()`), only the difference in means (`diff`) is
returned.

Use [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) with
`exponentiate = TRUE` to obtain risk ratios and odds ratios on their
natural scale.

## Variance estimation

Standard errors are computed via linearization, which correctly accounts
for the uncertainty introduced by estimating propensity scores. This
avoids the known problem of underestimated standard errors that arises
from treating estimated weights as fixed. See Kostouraki et al. (2024)
for details.

## References

Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of the
inverse probability-of-treatment weighting estimator: A tutorial for
different types of propensity score weights. *Statistics in Medicine*.
2024;43(13):2672–2694.
[doi:10.1002/sim.10078](https://doi.org/10.1002/sim.10078)

## See also

[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
[`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for calculating propensity score weights.

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`ps_trunc()`](https://r-causal.github.io/propensity/reference/ps_trunc.md)
for handling extreme propensity scores before weighting.

## Examples

``` r
# Simulate data with a confounder, binary exposure, and binary outcome
set.seed(123)
n <- 200
x1 <- rnorm(n)
z <- rbinom(n, 1, plogis(0.5 * x1))
y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
dat <- data.frame(x1, z, y)

# Step 1: Fit a propensity score model
ps_mod <- glm(z ~ x1, data = dat, family = binomial())

# Step 2: Calculate ATE weights and fit a weighted outcome model
wts <- wt_ate(ps_mod)
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1
outcome_mod <- glm(y ~ z, data = dat, family = binomial(), weights = wts)
#> Warning: non-integer #successes in a binomial glm!

# Step 3: Estimate causal effects with correct standard errors
result <- ipw(ps_mod, outcome_mod)
result
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Propensity Score Model:
#>   Call: glm(formula = z ~ x1, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z, family = binomial(), data = dat, weights = wts) 
#> 
#> Estimates:
#>         estimate std.err       z ci.lower ci.upper conf.level   p.value    
#> rd       0.14230 0.07038 2.02194   0.0044  0.28025       0.95 0.0431831 *  
#> log(rr)  0.28031 0.10770 2.60262   0.0692  0.49141       0.95 0.0092513 ** 
#> log(or)  0.57339 0.16200 3.53950   0.2559  0.89090       0.95 0.0004009 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Exponentiate log-RR and log-OR to get RR and OR
as.data.frame(result, exponentiate = TRUE)
#>   effect  estimate   std.err        z    ci.lower  ci.upper conf.level
#> 1     rd 0.1423042 0.0703802 2.021935 0.004361509 0.2802468       0.95
#> 2     rr 1.3235458 0.1077045 2.602624 1.071669109 1.6346215       0.95
#> 3     or 1.7742759 0.1619980 3.539503 1.291600480 2.4373286       0.95
#>        p.value
#> 1 0.0431831009
#> 2 0.0092513440
#> 3 0.0004008815

# Continuous outcome example
y_cont <- 2 + 0.8 * z + 0.3 * x1 + rnorm(n)
dat$y_cont <- y_cont
outcome_cont <- lm(y_cont ~ z, data = dat, weights = wts)
ipw(ps_mod, outcome_cont)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Propensity Score Model:
#>   Call: glm(formula = z ~ x1, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: lm(formula = y_cont ~ z, data = dat, weights = wts) 
#> 
#> Estimates:
#>      estimate std.err       z ci.lower ci.upper conf.level   p.value    
#> diff  0.90057 0.13695 6.57579   0.6322    1.169       0.95 4.839e-11 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
