# Print a result whose standard errors are a diagnostic

`se_method = "robust"` reports the sandwich the weighted outcome model
computes for itself, which treats the estimated weights as known.
Printing such a result writes what
[print()](https://rdrr.io/r/base/print.html) writes for any
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
result and then one line naming the method, so that the number beside
each estimate is not read as one accounting for the propensity score
model. See **Standard errors as a diagnostic** in
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for what the diagnostic does and does not cover.

## Usage

``` r
# S3 method for class 'ipw_diagnostic_se'
print(x, ...)

# S3 method for class 'ipw_diagnostic_se'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'ipw_diagnostic_se'
tidy(x, ...)
```

## Arguments

- x:

  An `ipw` result fit with `se_method = "robust"`.

- ...:

  Passed to the next method.

- row.names, optional:

  Passed to the next method.

## Value

`x`, invisibly.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
what it returns for any
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
result, carrying an `ipw_se_diagnostic` attribute naming the method.

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns what
it returns for any
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
result, carrying an `ipw_se_diagnostic` attribute naming the method.

## See also

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for the estimator and the standard error methods.

## Examples

``` r
set.seed(2)
n <- 300
x <- rnorm(n)
z <- rbinom(n, 1, plogis(0.4 * x))
y <- rbinom(n, 1, plogis(-0.3 + 0.7 * z + 0.5 * x))
dat <- data.frame(x, z, y)

ps_mod <- glm(z ~ x, data = dat, family = binomial())
wts <- wt_ate(ps_mod)
#> ℹ Using exposure variable "z" from the propensity score model
#> ℹ Treating `.exposure` as binary
outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

print(ipw(ps_mod, outcome_mod, se_method = "robust"))
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z, family = quasibinomial(), data = dat, weights = wts) 
#> 
#> Marginal estimates:
#>                estimate  std.err       z ci.lower ci.upper conf.level   p.value
#> mean 0         0.398460 0.041088  9.6978 0.317929  0.47899       0.95 < 2.2e-16
#> mean 1         0.585255 0.042736 13.6947 0.501494  0.66902       0.95 < 2.2e-16
#> rd 1 vs 0      0.186795 0.059284  3.1509 0.070601  0.30299       0.95  0.001628
#> log(rr) 1 vs 0 0.384441 0.126353  3.0426 0.136794  0.63209       0.95  0.002345
#> log(or) 1 vs 0 0.756270 0.245729  3.0777 0.274649  1.23789       0.95  0.002086
#>                   
#> mean 0         ***
#> mean 1         ***
#> rd 1 vs 0      ** 
#> log(rr) 1 vs 0 ** 
#> log(or) 1 vs 0 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Standard errors: robust, a diagnostic that treats the weights as known
```
