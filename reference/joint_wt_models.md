# Record the two treatment models of a joint exposure

`joint_wt_models()` records the two fitted treatment models a joint
weight is built from, in the order the joint density factorizes: the
first treatment's model, then the second treatment's model, which
conditions on the first. `is_joint_wt_models()` reports whether an
object is such a record.

The container exists to refuse. It checks that the pair describes a
sequential factorization at all, which is the one thing about a joint
weight that cannot be detected from the weight vector afterwards.

## Usage

``` r
joint_wt_models(...)

is_joint_wt_models(x)
```

## Arguments

- ...:

  Exactly two fitted treatment models, each named for the treatment it
  fits. The order is the factorization's order. Supported models are a
  binomial [`stats::glm()`](https://rdrr.io/r/stats/glm.html) for a
  binary treatment, a
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) for a
  categorical one, and an
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html) or gaussian
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) for a continuous
  one.

- x:

  An object to test.

## Value

`joint_wt_models()` returns an S3 object of class `joint_wt_models`, a
list with three fields:

- `names`:

  The two treatment names, in the factorization's order.

- `models`:

  The two fitted models, named for their treatments.

- `exposure_type`:

  Each model's exposure type, named for its treatment: `"binary"`,
  `"categorical"`, or `"continuous"`.

`is_joint_wt_models()` returns a single logical.

## What is checked

Each argument is named for the treatment its model fits, and the name
has to be that model's response: the factorization check looks for the
first treatment by name in the second model's formula, so a model
recorded under the wrong name would be checked for the wrong variable.

The second model must condition on the first treatment. The product of
two marginal models, \\f(A \| L) f(E \| L)\\, is not the joint weight
\\f(A \| L) f(E \| A, L)\\ for any data in which the second treatment
depends on the first, and nothing downstream can tell the two apart.
What is checked is whether the second model's right-hand side *mentions*
the first treatment, so `e ~ factor(a) * x1` and `e ~ x1 + a:x1` both
pass, as does the additive `e ~ a + x1`.

Passing is not the same as modeling the dependence well. An additive
term satisfies the factorization and may still model it badly, and a
badly modeled dependence biases the joint weight invisibly. Model the
dependence on the first treatment flexibly, with an interaction or a
spline, unless you have a reason to believe it is additive.

Neither may the two models condition on each other. A sequential
factorization has a first factor that is marginal in the second
treatment, so a pair that each read the other's treatment are not the
two factors of any factorization, in either order.

## See also

[`wt_joint()`](https://r-causal.github.io/propensity/reference/wt_joint.md)
for the product weight itself. These are used by
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

models <- joint_wt_models(
  a = glm(a ~ x1, data = dat, family = binomial()),
  e = glm(e ~ a * x1, data = dat, family = binomial())
)
models$names
#> [1] "a" "e"
models$exposure_type
#>        a        e 
#> "binary" "binary" 
is_joint_wt_models(models)
#> [1] TRUE
```
