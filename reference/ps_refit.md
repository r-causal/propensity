# Refit the Propensity Score Model on Retained Observations

Takes a `ps_trim` object and the original model used to calculate the
propensity score, then:

1.  Retrieves data from the model (or from `.df` argument if provided)

2.  Subsets rows to the non‐trimmed indices

3.  Refits the model

4.  Predicts new propensity scores for all rows (trimmed rows -\> `NA`)

5.  Returns a new `ps_trim` object with `refit = TRUE`.

## Usage

``` r
ps_refit(trimmed_ps, model, .data = NULL, ...)
```

## Arguments

- trimmed_ps:

  A `ps_trim` object (same length as data, NAs for trimmed).

- model:

  The fitted model used to get the original PS (e.g. a glm).

- .data:

  Optional. A data frame. If `NULL`, we try to retrieve from `model`.

- ...:

  Additional arguments passed to
  [`update()`](https://rdrr.io/r/stats/update.html).

## Value

A new `ps_trim` object with updated propensity scores and
`ps_trim_meta(x)$refit` set to `TRUE`.

## See also

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md),
[`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md),
[`is_ps_trimmed()`](https://r-causal.github.io/propensity/reference/is_ps_trimmed.md)

## Examples

``` r
set.seed(2)
n <- 30
x <- rnorm(n)
z <- rbinom(n, 1, plogis(0.4 * x))
fit <- glm(z ~ x, family = binomial)
ps <- predict(fit, type = "response")

# trim and refit
refit <- ps_trim(ps, lower = .2, upper = .8) |>
  ps_refit(fit)

is_refit(refit)
#> [1] TRUE
```
