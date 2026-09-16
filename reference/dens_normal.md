# Density specifications for continuous exposures

Weights for a continuous exposure are a ratio of densities. Both
densities are evaluated on the standardized residual \\z_i = (A_i -
\mu_i) / \sigma\\, where \\A_i\\ is the exposure, \\\mu_i\\ the fitted
conditional mean, and \\\sigma\\ the residual spread, so a specification
describes a density on a standardized scale rather than on the scale of
the exposure itself.

These constructors name the family and record the parameters that
identify it. Pass one to the `.density` argument of
[`wt_ate()`](https://r-causal.github.io/propensity/reference/wt_ate.md)
or
[`wt_cens()`](https://r-causal.github.io/propensity/reference/wt_ate.md),
which also accept the strings `"normal"`, `"laplace"`, and `"kernel"`
for the families that need no parameters, and a bare function of one
argument, which is wrapped with `dens_fn()`.

- `dens_normal()` is the standard normal density, the default.

- `dens_laplace()` is the standard Laplace density, \\\exp(-\|z\|) /
  2\\, which puts more mass in the tails than the normal. Its scale is
  estimated under the Laplace itself, where it is the mean absolute
  residual, or by the root mean square of the residuals, which is the
  spread every family without an estimator of its own is read at.

- `dens_t()` is Student's t density with `df` degrees of freedom,
  heavier tailed still, and heavier the smaller `df` is. Its scale is
  estimated under the t itself, by maximum likelihood, or by the root
  mean square of the residuals.

- `dens_kernel()` is a kernel density estimate of the standardized
  residuals, fit with
  [`stats::density()`](https://rdrr.io/r/stats/density.html) and
  interpolated to each observation. It assumes no family at all, at the
  cost of a density that is not a smooth function of the model's
  parameters: weights built from it have no closed-form standard error,
  and
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
  reports none for them at all. Bootstrap the whole fit by hand to put
  an interval around such an estimate.

- `dens_fn()` is a density you write yourself.

## Usage

``` r
dens_normal()

dens_laplace(sigma_method = c("mle", "rms"))

dens_t(df, sigma_method = c("mle", "rms"))

dens_kernel(bw = "nrd0", adjust = 1, kernel = "gaussian", n = 512)

dens_fn(f)
```

## Arguments

- sigma_method:

  How the spread of the conditional density is estimated from the
  residuals of the propensity score model, either `"rms"`, the root mean
  square, or `"mle"`, the scale estimated under the family itself, which
  both `dens_laplace()` and `dens_t()` take by default. See the section
  below.

- df:

  Degrees of freedom for Student's t, a single positive, finite number.

- bw:

  The bandwidth passed to
  [`stats::density()`](https://rdrr.io/r/stats/density.html): a single
  positive number, or the name of one of its selection rules (`"nrd0"`,
  `"nrd"`, `"ucv"`, `"bcv"`, `"SJ"`, `"SJ-ste"`, or `"SJ-dpi"`).

- adjust:

  A single positive number the bandwidth is multiplied by, as in
  [`stats::density()`](https://rdrr.io/r/stats/density.html). Values
  above 1 smooth the estimate further.

- kernel:

  The smoothing kernel, one of the kernels
  [`stats::density()`](https://rdrr.io/r/stats/density.html) accepts.

- n:

  The number of grid points
  [`stats::density()`](https://rdrr.io/r/stats/density.html) evaluates,
  at least 2. The estimate is interpolated between them, so a larger `n`
  follows the shape of the residuals more closely.

- f:

  A function of one argument, the standardized residual, returning one
  non-negative, finite density value for each element it is given.

## Value

An object of class `propensity_density`: a list with the elements
`family`, the name of the family; `params`, the parameters that identify
it, which is an empty list for a family that takes none; and `fn`, the
function that evaluates the density, which is `NULL` for
`dens_kernel()`. `dens_t()` and `dens_laplace()` carry a fourth element,
`sigma_method`, the name of the estimator their scale is read with.

## The spread of the conditional density

Both densities are evaluated on a residual standardized by a spread, and
for the conditional density that spread is estimated from the residuals
of the propensity score model. A family with no scale estimator of its
own takes the root mean square of those residuals, which is
`sigma_method = "rms"` and is the maximum likelihood estimator of the
spread of a normal density.

It is not the maximum likelihood estimator of the scale of a
heavier-tailed family. The scale of such a family is smaller than its
standard deviation, and the root mean square is pulled outward by the
large residuals a heavy tail produces, which is what the family was
chosen to accommodate. `sigma_method = "mle"` estimates the scale under
the family itself, and both `dens_laplace()` and `dens_t()` take it by
default, since the spread the density is divided by is the scale
parameter of the family reading it. For `dens_laplace()` it is the mean
absolute residual. For `dens_t()` it is the root of \\\sum_i
\left\[(\nu + 1) r_i^2 / (\nu \sigma^2 + r_i^2)\right\] = n\\, where
\\r_i\\ is the residual and \\\nu\\ is `df`. Each residual enters that
sum through a bounded term, so a residual far out in the tail moves the
estimate by less than it moves the root mean square. Prefer it when the
residuals are heavy tailed, which is the case these families are for;
the two estimators answer the same question, and answer it alike, as
`df` grows and the t approaches the normal.

`dens_normal()` offers no such choice, because the root mean square is
already its maximum likelihood estimator. Neither do `dens_kernel()` and
`dens_fn()`: a kernel is fit to the standardized residuals and divided
by the same spread, so the spread cancels up to the grid the estimate is
interpolated on, and a density you write yourself names no family an
estimator could be derived from.

The choice describes both densities of the ratio. The marginal density
that stabilizes the weights is the exposure's own, read at the
exposure's mean and at the spread the same estimator gives for it, so
the two halves of the ratio are densities of the same width. A scale
estimated by maximum likelihood is recorded by
[`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md)
as `sigma = "mle"`, and
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
estimates it alongside the propensity score model's coefficients,
solving the equation it is the root of as part of its stacked system so
that the standard errors account for it. Supplying a `.sigma` says the
spread is a number of your own rather than one estimated from the
residuals, so the two cannot be given together.

## Examples

``` r
dens_normal()
#> <density: normal>

dens_t(df = 4)
#> <density: t(df = 4)>

dens_t(df = 4, sigma_method = "mle")
#> <density: t(df = 4)>

dens_laplace(sigma_method = "rms")
#> <density: laplace>

dens_kernel(adjust = 1.5)
#> <density: kernel(bw = "nrd0", adjust = 1.5, kernel = "gaussian", n = 512)>

dens_fn(function(z) stats::dt(z, df = 4))
#> <density: function>
```
