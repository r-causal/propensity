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
  2\\, which puts more mass in the tails than the normal.

- `dens_t()` is Student's t density with `df` degrees of freedom,
  heavier tailed still, and heavier the smaller `df` is. Its scale is
  estimated by the root mean square of the residuals, as every other
  family's spread is, or by maximum likelihood under the t itself.

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

dens_laplace()

dens_t(df, sigma_method = c("rms", "mle"))

dens_kernel(bw = "nrd0", adjust = 1, kernel = "gaussian", n = 512)

dens_fn(f)
```

## Arguments

- df:

  Degrees of freedom for Student's t, a single positive, finite number.

- sigma_method:

  How the spread of the conditional density is estimated from the
  residuals of the propensity score model, either `"rms"` or `"mle"`.
  See the section below.

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
`dens_kernel()`. `dens_t()` carries a fourth element, `sigma_method`,
the name of the estimator its scale is read with.

## The spread of a t density

Both densities are evaluated on a residual standardized by a spread, and
for the conditional density that spread is estimated from the residuals
of the propensity score model. Every family takes the root mean square
of those residuals, which is `sigma_method = "rms"`, the default, and is
the maximum likelihood estimator of the spread of a normal density.

It is not the maximum likelihood estimator of the scale of a t. The
scale of a t is smaller than its standard deviation by a factor that
depends on `df`, and the root mean square is pulled outward by the large
residuals a heavy tail produces, which is what the family was chosen to
accommodate. `sigma_method = "mle"` estimates the scale under the t
itself, as the root of \\\sum_i \left\[(\nu + 1) r_i^2 / (\nu \sigma^2 +
r_i^2)\right\] = n\\, where \\r_i\\ is the residual and \\\nu\\ is `df`.
Each residual enters that sum through a bounded term, so a residual far
out in the tail moves the estimate by less than it moves the root mean
square. Prefer it when the residuals are heavy tailed, which is the case
the t family is for; the two estimators answer the same question, and
answer it alike, as `df` grows and the t approaches the normal.

The choice describes the conditional density alone. The marginal density
that stabilizes the weights is the exposure's own, read at the
exposure's mean and root mean square, whatever the conditional spread
was estimated by. A scale estimated by maximum likelihood is recorded by
[`density_meta()`](https://r-causal.github.io/propensity/reference/exposure_type.md)
as `sigma = "mle"`, and
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.html)
estimates it alongside the propensity score model's coefficients,
solving the equation above as part of its stacked system so that the
standard errors account for it. Supplying a `.sigma` says the spread is
a number of your own rather than one estimated from the residuals, so
the two cannot be given together.

## Examples

``` r
dens_normal()
#> <density: normal>

dens_t(df = 4)
#> <density: t(df = 4)>

dens_t(df = 4, sigma_method = "mle")
#> <density: t(df = 4)>

dens_kernel(adjust = 1.5)
#> <density: kernel(bw = "nrd0", adjust = 1.5, kernel = "gaussian", n = 512)>

dens_fn(function(z) stats::dt(z, df = 4))
#> <density: function>
```
