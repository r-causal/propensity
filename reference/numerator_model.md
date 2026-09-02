# The model a set of weights was stabilized on

`numerator_model()` returns the fitted model supplied to `stabilize`
that estimated the numerator of a set of weights, and `NULL` for weights
whose numerator was not estimated by one.

## Usage

``` r
numerator_model(wt)
```

## Arguments

- wt:

  A `psw` or `causal_wts` object.

## Value

The fitted model supplied to `stabilize`, or `NULL`.

## Details

A model supplied to
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)'s
`stabilize` argument estimates the numerator of the weights, and what it
estimates follows the exposure. A
[`binomial()`](https://rdrr.io/r/stats/family.html) fit of a binary
exposure reports the probability of the level each unit took; an
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) fit of
a categorical exposure reports one probability per level, and the
numerator is the column named for the level each unit took; a model of a
dose reports the conditional mean the density family is read at, and the
numerator is that density. The model itself is recorded rather than the
numerator it evaluates to, because
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
rebuilds that numerator at every value of its parameter vector, which
takes the model's design and its coefficients.

Where the record is kept differs by exposure type, and this accessor is
the one place to read it either way. Weights for a continuous exposure
are a ratio of densities and keep the model inside the density record
[`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md)
returns; weights for a binary or categorical exposure are a ratio of
probabilities, leave no density record, and keep the model on its own.
Both come back here.

The record describes the numerator rather than the units, so it holds no
length of its own and survives subsetting, arithmetic, and anything else
that changes the number of weights. Two sets of weights record the same
numerator when their models write the same formula and were fit to the
same coefficients; combining two that record different ones drops the
record, with a warning of class `propensity_metadata_conflict_warning`.
Weights stabilized without a model were divided by a numerator no model
estimated, which is a different numerator again: combining those with
weights stabilized on a model drops the model the same way, and the
result is stabilized on nothing it can name.

A `stabilization_score` is a numerator the caller computed rather than
one a model estimated, so weights stabilized on a score record no model,
and
[`stabilization_score()`](https://r-causal.github.io/propensity/reference/psw.md)
returns `NULL` for weights stabilized on a model.

## See also

[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for the `stabilize` argument that writes this record,
[`stabilization_score()`](https://r-causal.github.io/propensity/reference/psw.md)
for the numerator a caller supplies as a vector, and
[`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md)
for the rest of what a continuous exposure's weights record.

## Examples

``` r
set.seed(1)
x <- rnorm(50)
v <- rbinom(50, 1, 0.5)
z <- rbinom(50, 1, plogis(0.4 * x - 0.7 * v))

ps_mod <- glm(z ~ x + v, family = binomial())
num_mod <- glm(z ~ v, family = binomial())

w <- wt_ate(ps_mod, stabilize = num_mod)
#> ℹ Using exposure variable "z" from the propensity score model
#> ℹ Treating `.exposure` as binary
numerator_model(w)
#> 
#> Call:  glm(formula = z ~ v, family = binomial())
#> 
#> Coefficients:
#> (Intercept)            v  
#>     -0.3365      -0.1335  
#> 
#> Degrees of Freedom: 49 Total (i.e. Null);  48 Residual
#> Null Deviance:       67.3 
#> Residual Deviance: 67.25     AIC: 71.25

# Weights stabilized on the marginal probability record no model.
numerator_model(wt_ate(ps_mod, stabilize = TRUE))
#> ℹ Using exposure variable "z" from the propensity score model
#> ℹ Treating `.exposure` as binary
#> NULL
```
