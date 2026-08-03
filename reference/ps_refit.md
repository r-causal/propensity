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
  Refitting reads the retained positions out of the trimming record, so
  an object whose record was dropped or no longer covers it raises an
  error of class `propensity_missing_meta_error`; see
  [`ps_trim()`](https://r-causal.github.io/propensity/reference/ps_trim.md).

- model:

  The original fitted model used to estimate the propensity scores (e.g.
  a [glm](https://rdrr.io/r/stats/glm.html) or
  [multinom](https://rdrr.io/pkg/nnet/man/multinom.html) object). The
  model is refit via [update()](https://rdrr.io/r/stats/update.html) on
  the retained subset.

- .data:

  A data frame with one row per observation in `trimmed_ps`, in the same
  order. If `NULL` (the default), the data are recovered from `model`:
  its [model.frame()](https://rdrr.io/r/stats/model.frame.html) when
  that already holds every variable the refit reads, and otherwise the
  data the model names, restricted by row name to the rows the model
  analyzed. A model fit without a data argument names none, and its
  variables are read out of the formula's environment instead. A formula
  that transforms a term, such as `z ~ log(x)` or a spline basis, stores
  that term already computed, so only the underlying variables let the
  transformation be recomputed from the retained rows. Pass `.data` when
  the data the model was fit on can no longer be reached.

- ...:

  Additional arguments passed to
  [update()](https://rdrr.io/r/stats/update.html).

## Value

A `ps_trim` object with re-estimated propensity scores for retained
observations and `NA` for trimmed observations. Use
[`is_refit()`](https://r-causal.github.io/propensity/reference/is_refit.md)
to confirm refitting was applied.

## Details

### Composing with a `subset`

A `subset` in the original call has already chosen the sample the
propensity scores are about, and the trimming record indexes that sample
rather than every row the data carry. Refitting narrows that sample
further, to the retained rows, so the original `subset` is dropped from
the call rather than put to work a second time on rows it was never
about. A `subset` passed through `...` is an instruction of its own and
is honored.

### Arguments read from outside the formula

`weights`, `offset`, and `na.action` in the original call are
re-evaluated against the retained rows. A `weights` or `offset` naming a
column of the data the model was fit on is read from that column and
follows the retained rows, whether the data are recovered from `model`
or passed to `.data`. A vector held outside the data cannot follow them:
it keeps the length it had and raises an error about differing variable
lengths.

Scores predicted from a fit with `na.action = na.exclude` are padded
back to the full length of the data, so they describe more observations
than the fit read and the trimming record indexes a sample the model
never analyzed. `ps_refit()` refuses such scores. Trim scores from a fit
whose `na.action` drops those rows instead.

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

is_refit(refit)
#> [1] TRUE
```
