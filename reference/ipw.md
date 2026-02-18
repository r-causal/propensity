# Inverse Probability Weights for Causal Inference

`ipw()` is a bring-your-own-model (BYOM) inverse probability weighted
estimator. `ipw()` accepts a propensity score model and a weighted
outcome model that you have already fit. The purpose of `ipw()` is to
capture the uncertainty inherent to this two-step process and calculate
the correct standard errors for each estimate. Currently, `ipw()`
supports binary exposures and either binary or continuous outcomes. For
binary outcomes, `ipw()` calculates the marginal risk difference, log
risk ratio, and log odds ratio. For continuous outcomes, `ipw()` only
calculates the marginal difference in means.

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
  logistic regression with the exposure as the dependent variable.

- outcome_mod:

  A fitted, weighted outcome model of class
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html), with the outcome as
  the dependent variable.

- .data:

  A data frame containing the exposure, outcome, and covariates. If
  `NULL`, `ipw()` will try to extract the data from `ps_mod` and
  `outcome_mod`.

- estimand:

  A character string specifying the causal estimand: `ate`, `att`,
  `ato`, or `atm`. If `NULL`, the function attempts to infer it from
  existing weights in `outcome_mod`, assuming they were calculated with
  [`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
  [`wt_att()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
  [`wt_atm()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
  or
  [`wt_ato()`](https://r-causal.github.io/propensity/reference/wt_ate.md).

- ps_link:

  A character string specifying the link function for the propensity
  score model: `logit`, `probit`, or `cloglog`. Defaults to whatever
  link was used by `ps_mod`.

- conf_level:

  Numeric. Confidence level for intervals (default `0.95`).

- x:

  an `ipw` object

- row.names, optional, ...:

  Passed to
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- exponentiate:

  Logical. Should the log-RR and log-OR be exponentiated?

## Value

An S3 object of class `ipw` containing:

- `estimand`: One of `"ate"`, `"att"`, `"ato"`, or `"atm"`.

- `ps_mod`: The fitted propensity score model.

- `outcome_mod`: The fitted outcome model.

- `estimates`: A data frame of point estimates, standard errors,
  z-statistics, confidence intervals, and p-values.

## Details

The function constructs inverse probability weights based on the chosen
`estimand`, then uses these weights (or extracts them from
`outcome_mod`) to compute effect measures:

- `rd`: Risk difference

- `log(rr)`: Log risk ratio

- `log(or)`: Log odds ratio

  For a linear outcome model (using
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) with
  `family = gaussian()`), only the difference in means (`diff`) is
  returned.

**Variance Estimation**

The variance is estimated via linearization, which provide variance
estimates for IPW that correctly account for the uncertainty in
estimation of the propensity scores. For more details on various types
of propensity score weights and their corresponding variance estimators,
see:

- *Kostouraki A, Hajage D, Rachet B, et al. On variance estimation of
  the inverse probability-of-treatment weighting estimator: A tutorial
  for different types of propensity score weights. Statistics in
  Medicine. 2024; 43(13): 2672-2694. doi:
  [10.1002/sim.10078](https://doi.org/10.1002/sim.10078)*

## See also

[`stats::glm()`](https://rdrr.io/r/stats/glm.html),
[`stats::lm()`](https://rdrr.io/r/stats/lm.html)

## Examples

``` r
set.seed(123)
n <- 100
# confounder
x1 <- rnorm(n)
# exposure
z  <- rbinom(n, 1, plogis(0.2 * x1))
# binary outcome
y  <- rbinom(n, 1, plogis(1 + 2*z + 0.5*x1))

dat <- data.frame(x1, z, y)

# fit a propensity score model (exposure ~ x1)
ps_mod <- glm(z ~ x1, data = dat, family = binomial())

# calculate weights for ATE
ps <- predict(ps_mod, type = "response")
wts <- wt_ate(ps, z)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1

# fit an outcome model (binary y ~ z) using IPW
outcome_mod <- glm(y ~ z, data = dat, family = binomial(), weights = wts)
#> Warning: non-integer #successes in a binomial glm!

# run IPW
ipw_res <- ipw(ps_mod, outcome_mod)

ipw_res
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
#>         estimate  std.err        z ci.lower ci.upper conf.level   p.value    
#> rd       0.16311 0.076670 2.127466   0.0128  0.31338       0.95   0.03338 *  
#> log(rr)  0.19720 0.080683 2.444134   0.0391  0.35533       0.95   0.01452 *  
#> log(or)  1.24122 0.170217 7.291963   0.9076  1.57484       0.95 3.055e-13 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Convert to a data frame with exponentiated RR/OR
ipw_res_df <- as.data.frame(ipw_res, exponentiate = TRUE)
ipw_res_df
#>   effect  estimate    std.err        z   ci.lower  ci.upper conf.level
#> 1     rd 0.1631125 0.07666988 2.127466 0.01284233 0.3133827       0.95
#> 2     rr 1.2179865 0.08068260 2.444134 1.03983716 1.4266572       0.95
#> 3     or 3.4598274 0.17021737 7.291963 2.47836430 4.8299620       0.95
#>        p.value
#> 1 3.338141e-02
#> 2 1.452002e-02
#> 3 3.055334e-13
```
