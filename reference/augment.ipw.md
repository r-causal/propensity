# Augment an inverse probability weighted result with per-observation columns

[`augment()`](https://generics.r-lib.org/reference/augment.html) works
per observation rather than per estimate: the data an
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
result was produced from, carried through in full, with the propensity
score, the weights, the fitted values, and the residuals attached as
dot-prefixed columns. Nothing is dropped and no column of the source
frame is changed.

The added columns are read from the two models the result holds, so they
describe the fit rather than the frame they are attached to. `data` says
which frame to carry them on, which is how broom uses the argument: it
names the data the fit was produced from, not new data to predict on.
Supplying the modeling data rather than leaving the default is how
covariates the outcome formula left out arrive beside the fit's own
columns.

The weights are reported once, as `.weights`, which is a deliberate
departure from broom. A broom
[`augment()`](https://generics.r-lib.org/reference/augment.html) method
carries the model frame's `(weights)` column through and adds no weights
column of its own; here the weights are the central per-observation
quantity of the analysis, and the
[`psw()`](https://r-causal.github.io/propensity/reference/psw.md) vector
reports them with the estimand they target attached, which the model
frame's plain numeric copy of the same values does not. The default
frame therefore leaves `(weights)` out. A frame supplied through `data`
is the caller's own and is carried as it arrives, so a weight column
they hold stays where they put it, whatever it is named.

A column is an answer per observation only while the fit it is read from
and the frame it is carried on describe the same observations, so three
things are refused rather than reported. A result whose propensity score
model and outcome model were fit to different rows is refused, because
the propensity score of one observation would arrive beside the fitted
value of another. Two models of the same data most often part this way
over missing values, each dropping the rows a variable of its own is
missing on. What is compared is the number of observations each model
produced an answer for and, when the outcome model's frame names the
exposure, the exposure values of the two model frames position by
position, reading a factor and the numbers its labels spell as one
encoding of one exposure. Exposure values that disagree prove the two
models hold different observations; exposure values that agree prove
nothing, since two different sets of rows are free to carry the same
sequence of values, and two encodings that cannot be read onto each
other, such as a factor of `"a"` and `"b"` against a recoding of it as 0
and 1, prove nothing either way. The check is one-sided in that
direction on purpose: it refuses only what it can prove, so a fit whose
models do line up is never refused over the labels the rows of either
frame carry or the encoding its exposure is written in. A frame that
already holds a column this call would add is refused as well, because
the fit's column has nowhere to go: added beside the frame's, the two
would be told apart by nothing, and written over it, the frame's data
would be gone. Which names are in the way depends on the fit and the
frame, `.resid` being among them only for a frame it would be added to
and the propensity columns being the ones the propensity score model
produced. Last, a result whose outcome model was fit without weights is
refused under the class
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
refuses the same fact with: `.weights` has nothing to report for such a
fit.

## Usage

``` r
# S3 method for class 'ipw'
augment(x, data = NULL, ...)
```

## Arguments

- x:

  An `ipw` object, as returned by
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html).

- data:

  A data frame with one row for each observation the fit used, or
  `NULL`, the default, to use the outcome model's own model frame. That
  frame holds the response and the terms of the outcome formula, and
  both are kept. The `(weights)` column
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  records is the one thing left out of it, because those weights are
  reported as `.weights`.

- ...:

  These dots are for future extensions and must be empty.

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per observation, holding every column of the source frame in its own
order followed by:

- `.propensity`:

  The propensity score, as the propensity score model predicts it: the
  probability of exposure for a binary exposure, and the conditional
  mean of the exposure for a continuous one. This is the propensity
  score the weight functions take, so feeding it back with the exposure
  returns the weights the outcome model was fit with.

- `.propensity_<level>`:

  For a categorical exposure, which has a probability for every level
  rather than a single number, one column per level in place of
  `.propensity`, named for the level and ordered as the propensity score
  model orders its levels.

- `.weights`:

  The weights the outcome model was fit with, as the
  [`psw()`](https://r-causal.github.io/propensity/reference/psw.md)
  vector they were supplied as, so the estimand they record travels with
  them.

- `.fitted`:

  The outcome model's fitted values, on the response scale.

- `.resid`:

  The observed outcome less `.fitted`. Present only when the source
  frame holds the outcome the model was fit on, which the default frame
  always does. A factor or logical outcome is differenced on the 0/1
  scale its fitted values are on.

## See also

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
for the estimator,
[`tidy()`](https://r-causal.github.io/propensity/reference/tidy.ipw.md)
for its estimates, and
[`glance()`](https://r-causal.github.io/propensity/reference/glance.ipw.md)
for a one-row summary of the fit.

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
#> ℹ Using exposure variable "z" from GLM model
#> ℹ Treating `.exposure` as binary
outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
result <- ipw(ps_mod, outcome_mod)

augment(result)
#> # A tibble: 200 × 6
#>        y     z .propensity   .weights .fitted .resid
#>    <int> <int>       <dbl> <psw{ate}>   <dbl>  <dbl>
#>  1     1     1       0.471   2.123269   0.582  0.418
#>  2     1     0       0.495   1.978251   0.440  0.560
#>  3     0     0       0.620   2.629856   0.440 -0.440
#>  4     0     0       0.516   2.065896   0.440 -0.440
#>  5     1     1       0.520   1.922573   0.582  0.418
#>  6     0     1       0.630   1.586777   0.582 -0.582
#>  7     0     0       0.544   2.191512   0.440 -0.440
#>  8     1     0       0.421   1.728108   0.440  0.560
#>  9     0     0       0.462   1.858725   0.440 -0.440
#> 10     1     1       0.479   2.087063   0.582  0.418
#> # ℹ 190 more rows

# The modeling data carries through the covariates the outcome formula left
# out
augment(result, data = dat)
#> # A tibble: 200 × 7
#>         x1     z     y .propensity   .weights .fitted .resid
#>      <dbl> <int> <int>       <dbl> <psw{ate}>   <dbl>  <dbl>
#>  1 -0.560      1     1       0.471   2.123269   0.582  0.418
#>  2 -0.230      0     1       0.495   1.978251   0.440  0.560
#>  3  1.56       0     0       0.620   2.629856   0.440 -0.440
#>  4  0.0705     0     0       0.516   2.065896   0.440 -0.440
#>  5  0.129      1     1       0.520   1.922573   0.582  0.418
#>  6  1.72       1     0       0.630   1.586777   0.582 -0.582
#>  7  0.461      0     0       0.544   2.191512   0.440 -0.440
#>  8 -1.27       0     1       0.421   1.728108   0.440  0.560
#>  9 -0.687      0     0       0.462   1.858725   0.440 -0.440
#> 10 -0.446      1     1       0.479   2.087063   0.582  0.418
#> # ℹ 190 more rows
```
