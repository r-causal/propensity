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
  heavier tailed still, and heavier the smaller `df` is.

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

dens_t(df)

dens_kernel(bw = "nrd0", adjust = 1, kernel = "gaussian", n = 512)

dens_fn(f)
```

## Arguments

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
`family`, `params`, and `fn`. `fn` is the function that evaluates the
density, and is `NULL` for `dens_kernel()`.

## Examples

``` r
dens_normal()
#> <density: normal>

dens_t(df = 4)
#> <density: t(df = 4)>

dens_kernel(adjust = 1.5)
#> <density: kernel(bw = "nrd0", adjust = 1.5, kernel = "gaussian", n = 512)>

dens_fn(function(z) stats::dt(z, df = 4))
#> <density: function>
```
