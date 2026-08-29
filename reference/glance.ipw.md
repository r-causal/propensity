# Glance at an inverse probability weighted result

[`glance()`](https://generics.r-lib.org/reference/glance.html) describes
an
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
result rather than its estimates: one row naming the estimand and
counting the observations and the residual degrees of freedom of the
system the standard errors came from. A fit reporting several effect
measures, or several contrasts of a categorical exposure, still returns
exactly one row.

Under M-estimation that system is the stacked estimating equations,
which hold the propensity score model, the outcome model, and the effect
measures at once. Its residual degrees of freedom are the observations
it was solved on less the parameters it solves for, so a multinomial
propensity score model leaves fewer of them than a binary one does on
the same data. Linearization stacks nothing and records no parameter
count, so the observations are the outcome model's and there is no count
to subtract from them.

The columns and their types are the same on every route
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
takes, so the rows of several results stack into one table.

## Usage

``` r
# S3 method for class 'ipw'
glance(x, ...)
```

## Arguments

- x:

  An `ipw` object, as returned by
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).

- ...:

  These dots are for future extensions and must be empty.

## Value

A one-row [tibble](https://tibble.tidyverse.org/reference/tibble.html)
with the columns:

- `estimand`:

  The estimand the weights target, such as `"ate"`.

- `nobs`:

  The number of observations the outcome model was fit on, which is also
  what the stacked estimating equations are solved on under
  M-estimation. Reported by
  [`stats::nobs()`](https://rdrr.io/r/stats/nobs.html).

- `df.residual`:

  The residual degrees of freedom of the stacked estimating equations,
  `nobs` less the number of parameters the system solves for. `NA` under
  linearization, which records no parameter count. Reported by
  [`stats::df.residual()`](https://rdrr.io/r/stats/df.residual.html).

## See also

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for the estimator,
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw.md)
for its estimates, and
[`augment()`](https://r-causal.github.io/propensity/reference/augment.ipw.md)
for its per-observation columns.

## Examples

``` r
set.seed(123)
n <- 200
x1 <- rnorm(n)
z <- rbinom(n, 1, plogis(0.5 * x1))
y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
dat <- data.frame(x1, z, y)

ps_mod <- glm(z ~ x1, data = dat, family = binomial())
wts <- wt_ate(ps_mod)
#> ℹ Using exposure variable "z" from the propensity score model
#> ℹ Treating `.exposure` as binary
outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
result <- ipw(ps_mod, outcome_mod)

glance(result)
#> # A tibble: 1 × 3
#>   estimand  nobs df.residual
#>   <chr>    <int>       <int>
#> 1 ate        200         191
```
