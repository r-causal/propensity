# Product weights for a joint intervention on two treatments

`wt_joint()` builds the weight for a joint intervention on two
treatments by multiplying the two component weight vectors. A joint
exposure needs a weight for the pair, and that weight factorizes
sequentially: the density of the first treatment given the covariates,
times the density of the second given the first treatment and the
covariates. `wt_joint()` is the product, plus the checks that make the
product a joint weight rather than an arbitrary one.

`is_joint_wt()` reports whether a set of weights was built this way, and
`joint_wt_meta()` reads back what each component was.

## Usage

``` r
wt_joint(w_a, w_e, exposure_type = NULL)

is_joint_wt(x)

joint_wt_meta(x)
```

## Arguments

- w_a, w_e:

  The two component weight vectors, as
  [`psw()`](https://r-causal.github.io/propensity/reference/psw.md)
  objects. Both must target the `ate`, both must be the same length, and
  a component weighting a continuous exposure must be stabilized. The
  order is the factorization's order: `w_a` weights the first treatment
  and `w_e` weights the second, whose model conditions on the first.

- exposure_type:

  A character vector of length two naming each component's exposure
  type, one of `"binary"`, `"categorical"`, or `"continuous"`, in the
  order the components were given. Defaults to the types the components
  record, which is what a weight function writes on the weights it
  builds. A value given here is used instead of what they record, and a
  component recording no type, such as one assembled by hand, is refused
  unless both types are named.

- x:

  An object to test, or the weights to read the record from.

## Value

`wt_joint()` returns a
[`psw()`](https://r-causal.github.io/propensity/reference/psw.md) vector
of the elementwise product, carrying a `joint_wt_meta` record.

`is_joint_wt()` returns a single logical.

`joint_wt_meta()` returns the record as a list of five elements, with a
sixth on a record that has dropped a stabilization score, each one per
component and in the order the components were given, or `NULL` for
weights that are not a product:

- `exposure_type`:

  The two components' exposure types, in order.

- `stabilized`:

  Whether each component was stabilized, in order.

- `density`:

  Each component's density record, as
  [`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md)
  returns it, in order. A component weighting a discrete exposure is
  `NULL`, since its weights are no ratio of densities.

- `numerator_model`:

  Each component's numerator model, as
  [`numerator_model()`](https://r-causal.github.io/propensity/reference/numerator_model.md)
  returns it, in order, for the components that have nowhere else to
  keep one. A component weighting a continuous exposure keeps its model
  inside its own `density` record, so this element is `NULL` for it, and
  so it is for a component stabilized on anything other than a fitted
  model. A record written before this element existed holds none at all,
  which says the same thing as a record whose components each record no
  model.

- `stabilization_score`:

  Each component's stabilization score, as
  [`stabilization_score()`](https://r-causal.github.io/propensity/reference/psw.md)
  returns it, in order. A component stabilized on anything other than a
  score the caller supplied is `NULL`. A record written before this
  element existed holds none at all, which says the same thing as a
  record whose components each record no score.

- `score_dropped`:

  Whether each component's stabilization score was dropped by an
  operation that changed the length of the weights, in order. A drop
  empties the component's `stabilization_score` slot, which is the slot
  a component stabilized on anything else has, so a component marked
  `TRUE` here was stabilized on a score that is gone rather than on a
  numerator an estimator can rebuild. A record holding no such element,
  which is what a product nothing has shortened carries, says no
  component's score was dropped.

The record names the components rather than the observations, so it
survives subsetting and every other operation that keeps the vector a
[`psw()`](https://r-causal.github.io/propensity/reference/psw.md). A
stabilization score is the exception: it holds one value per
observation, so a score that no longer describes the observations the
result holds is dropped the way the score on the weights themselves is,
and a score of one value is carried at any length. Such a drop is
recorded in `score_dropped`, so a component stabilized on a score that
is gone stays distinguishable from one stabilized on the marginal
numerator, and
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
names it rather than standing a numerator in for it silently. Arithmetic
that leaves the result unstabilized reduces the record's numerator side,
as described under **What the product records**.

## The factorization

The weight for a joint intervention on `A` and `E` is \\1 / \[f(A \| L)
f(E \| A, L)\]\\. The second factor conditions on the first treatment.
The product of two *marginal* weights, \\1 / \[f(A \| L) f(E \| L)\]\\,
is a different quantity, and it is not the joint weight for any data in
which `E` depends on `A`.

Nothing downstream can tell the two apart. Both are ordinary vectors of
positive numbers, both fit an outcome model without complaint, and both
return estimates with standard errors. That is why
[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
refuses a second model that does not condition on the first treatment,
and why the dependence should be modeled flexibly rather than as a
single additive term: an additive term satisfies the factorization but
may model it badly, and a badly modeled dependence biases the joint
weight in a way that is equally invisible.

## Stabilization

A continuous component must be stabilized. The unstabilized density
ratio \\1 / f(A \| L)\\ has a heavy right tail on its own, and
multiplying it by a second weight inherits that tail, leaving the
product with no usable variance. A continuous component is stabilized
unless it was built with `stabilize = FALSE`, so the requirement is met
by default. A binary or categorical component needs no stabilization and
is accepted either way.

Any numerator that stabilizes the ratio satisfies the requirement, and
any density the ratio is read in is accepted: a component built with
`numerator = "integrated"`, with a heavier-tailed `.density`, or with a
`stabilization_score` of its own is stabilized as much as the default
marginal one is. What each component was built with is recorded per
component rather than merged, so the product carries both answers.

## What the product records

The result is a
[`psw()`](https://r-causal.github.io/propensity/reference/psw.md) with
estimand `"ate"`, since the product of two `ate` weights targets the
joint `ate`. It is marked stabilized when *either* component is, because
a stabilizing numerator on either factor is a stabilizing numerator on
the product; `joint_wt_meta()` keeps the per-component truth, so the
coarse flag hides nothing.

The product records nothing about trimming, truncation, or calibration
of the propensity scores the components were built from. Those records
name the observations of one component, and the product is not that
component.

Arithmetic that leaves a result unstabilized, such as multiplying the
product by unstabilized weights, reduces the numerator side of the
record it carries: the components read as unstabilized, their numerator
models and stabilization scores are dropped, and each density record's
numerator fields read as the record of no numerator. Such a result was
divided by no numerator, so a record naming one would describe a
different vector. What each component's weights divide by stays, and so
does the component structure.

## See also

[`joint_wt_models()`](https://r-causal.github.io/propensity/reference/joint_wt_models.md)
for the container recording the two treatment models, and
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
for building the components. These are used by
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
methods for joint treatment models.

## Examples

``` r
set.seed(4)
n <- 200
x1 <- rnorm(n)
a <- rbinom(n, 1, plogis(0.3 * x1))
e <- rbinom(n, 1, plogis(-0.2 + 0.5 * x1 - 0.8 * a))
dat <- data.frame(x1, a, e)

# The second model conditions on the first treatment, and does so flexibly
mod_a <- glm(a ~ x1, data = dat, family = binomial())
mod_e <- glm(e ~ a * x1, data = dat, family = binomial())

w <- wt_joint(wt_ate(mod_a), wt_ate(mod_e))
#> ℹ Using exposure variable "a" from the propensity score model
#> ℹ Treating `.exposure` as binary
#> ℹ Using exposure variable "e" from the propensity score model
#> ℹ Treating `.exposure` as binary
head(w)
#> <psw{estimand = ate}[6]>
#> [1] 4.486306 4.196643 4.279964 5.496731 2.693151 2.638138
estimand(w)
#> [1] "ate"
is_joint_wt(w)
#> [1] TRUE
joint_wt_meta(w)
#> $exposure_type
#> [1] "binary" "binary"
#> 
#> $stabilized
#> [1] FALSE FALSE
#> 
#> $density
#> $density[[1]]
#> NULL
#> 
#> $density[[2]]
#> NULL
#> 
#> 
#> $numerator_model
#> $numerator_model[[1]]
#> NULL
#> 
#> $numerator_model[[2]]
#> NULL
#> 
#> 
#> $stabilization_score
#> $stabilization_score[[1]]
#> NULL
#> 
#> $stabilization_score[[2]]
#> NULL
#> 
#> 

# A continuous component must be stabilized, which it is by default
d <- 0.5 + 0.6 * x1 - 0.7 * a + rnorm(n)
mod_d <- lm(d ~ a * x1, data = dat)
w_d <- wt_ate(mod_d)
#> ℹ Using exposure variable "d" from the propensity score model
#> ℹ Treating `.exposure` as continuous
# Each component records the exposure type it weights, so the product knows
# which of them is the dose without being told
joint_wt_meta(wt_joint(wt_ate(mod_a), w_d))
#> ℹ Using exposure variable "a" from the propensity score model
#> ℹ Treating `.exposure` as binary
#> $exposure_type
#> [1] "binary"     "continuous"
#> 
#> $stabilized
#> [1] FALSE  TRUE
#> 
#> $density
#> $density[[1]]
#> NULL
#> 
#> $density[[2]]
#> density:   normal
#> numerator: marginal
#> sigma:     pooled
#> 
#> 
#> $numerator_model
#> $numerator_model[[1]]
#> NULL
#> 
#> $numerator_model[[2]]
#> NULL
#> 
#> 
#> $stabilization_score
#> $stabilization_score[[1]]
#> NULL
#> 
#> $stabilization_score[[2]]
#> NULL
#> 
#> 
```
