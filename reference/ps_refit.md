# Refit a Propensity Score Model on Retained Observations

Re-estimates a propensity score model using only the observations
retained after trimming. This is the recommended intermediate step
between
[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
and weight calculation (e.g.
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)):

**[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
-\> `ps_refit()` -\> `wt_*()`**

Trimming changes the target population by removing observations with
extreme propensity scores. Refitting the model on the retained subset
produces propensity scores that better reflect this population,
improving both model fit and downstream weight estimation. Weight
functions warn if a trimmed propensity score has not been refit.

## Usage

``` r
ps_refit(trimmed_ps, model, .data = NULL, ...)
```

## Arguments

- trimmed_ps:

  A `ps_trim` object returned by
  [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md).

- model:

  The original fitted model used to estimate the propensity scores (e.g.
  a [glm](https://rdrr.io/r/stats/glm.html) or
  [multinom](https://rdrr.io/pkg/nnet/man/multinom.html) object). The
  model is refit via [update()](https://rdrr.io/r/stats/update.html) on
  the retained subset.

- .data:

  A data frame. If `NULL` (the default), the data are extracted from
  `model` via [model.frame()](https://rdrr.io/r/stats/model.frame.html).

- ...:

  Additional arguments passed to
  [update()](https://rdrr.io/r/stats/update.html).

## Value

A `ps_trim` object with re-estimated propensity scores for retained
observations and `NA` for trimmed observations. Use
[`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md)
to confirm refitting was applied.

## See also

[`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md)
for the trimming step,
[`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md)
to check refit status,
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
and other weight functions for the next step in the pipeline.

## Examples

``` r
set.seed(2)
n <- 200
x <- rnorm(n)
z <- rbinom(n, 1, plogis(0.4 * x))

# fit a propensity score model
ps_model <- glm(z ~ x, family = binomial)
ps <- predict(ps_model, type = "response")

# trim -> refit -> weight pipeline
trimmed <- ps_trim(ps, lower = 0.1, upper = 0.9)
refit <- ps_refit(trimmed, ps_model)
wts <- wt_ate(refit, .exposure = z)
#> ℹ Treating `.exposure` as binary
#> ℹ Setting focal level to 1

is_refit(refit)
#> [1] TRUE
```
