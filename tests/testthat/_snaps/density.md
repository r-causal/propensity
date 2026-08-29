# the constructors say what they refuse

    Code
      expr
    Condition <propensity_density_error>
      Error in `dens_t()`:
      ! `df` must be a single positive, finite number.
      x It is -1.
      i Supply one number, greater than zero.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `dens_t()`:
      ! `df` must be a single positive, finite number.
      x It has length 2.
      i Supply one number, greater than zero.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `dens_t()`:
      ! `df` must be a single positive, finite number.
      x It has class <character>.
      i Supply one number, greater than zero.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `dens_kernel()`:
      ! `bw` must be a single positive, finite number.
      x It is 0.
      i Supply one number, greater than zero.

---

    Code
      expr
    Condition <rlang_error>
      Error in `dens_kernel()`:
      ! `bw` must be one of "nrd0", "nrd", "ucv", "bcv", "SJ", "SJ-ste", or "SJ-dpi", not "silverman".

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `dens_kernel()`:
      ! `adjust` must be a single positive, finite number.
      x It has class <character>.
      i Supply one number, greater than zero.

---

    Code
      expr
    Condition <rlang_error>
      Error in `dens_kernel()`:
      ! `kernel` must be one of "gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", or "optcosine", not "gausian".
      i Did you mean "gaussian"?

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `dens_kernel()`:
      ! `n` must be a single whole number of at least 2.
      x It is 512.5.
      i `n` is the number of points `stats::density()` evaluates the estimate at before it is interpolated.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `dens_fn()`:
      ! `f` must be a function.
      x It has class <character>.
      i Supply a function of the standardized residual, such as `function(z) dnorm(z)`.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `dens_fn()`:
      ! `f` must be a function of at least one argument.
      x It takes no arguments.
      i The density is called on the standardized residuals, so its first argument is the vector it is evaluated at.

# as_density_spec() says what it refuses

    Code
      expr
    Condition <rlang_error>
      Error:
      ! `.density` must be one of "normal", "laplace", or "kernel", not "kernal".
      i Did you mean "kernel"?

---

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` must be a density specification, a function, or a string.
      x It has class <numeric>.
      i Supply "normal", "laplace", or "kernel", one of `dens_normal()`, `dens_laplace()`, `dens_t()`, `dens_kernel()`, or `dens_fn()`, or a function of the standardized residual.

# check_density_values() says what went wrong

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` must return a numeric vector.
      x It returned <character>.
      i A density returns one number for each standardized residual it is given.

---

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` must return one value for each observation.
      x It returned 2 values for 4 observations.
      i The density is evaluated on the whole vector of standardized residuals at once, so it must be vectorized over its argument.

---

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` must not return missing values.
      x It returned 1 missing value.
      i The standardized residuals run from -0.213 to 0.982.
      i The density must be defined at every standardized residual, and a kernel density is defined only over the range it was fit on.

---

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` must return finite, non-negative values.
      x It returned 2 values that are negative or not finite.
      i It failed at standardized residuals from -0.108 to 0.469.
      i A density is non-negative everywhere. Check that `.density` returns a density rather than a log density or a distribution function.

---

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` must not be zero everywhere.
      x Every one of the 4 values it returned is zero.
      i The standardized residuals run from -0.213 to 0.982.
      i Every weight built from a density that is zero at every observation is infinite. Check that the family matches the scale of the standardized residuals.

---

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.numerator` must return finite, non-negative values.
      x It returned 1 value that is negative or not finite.
      i It failed at the standardized residual 0.982.
      i A density is non-negative everywhere. Check that `.numerator` returns a density rather than a log density or a distribution function.

# a kernel that cannot be fit says why

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` cannot fit a kernel on standardized residuals that do not vary.
      x They are every one 0.5, so the estimate has no range to be fit over.
      i A model that fits the exposure exactly, or an exposure that is constant, leaves residuals a kernel cannot smooth. Use a parametric family, such as `dens_normal()`.

---

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` cannot fit a kernel on fewer than two standardized residuals.
      x It was given 1 to fit on.
      i The bandwidth is chosen from the spread of the residuals, which takes at least two of them. Use a parametric family, such as `dens_normal()`, on a sample this small.

---

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` cannot fit a kernel on missing standardized residuals.
      x 1 of them is missing.
      i A missing exposure or fitted value leaves a missing residual. Drop those observations before weighting.

---

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` cannot fit a kernel over a missing range.
      x The range a kernel is fit over is read from the standardized residuals, and at least one of them is missing.
      i A missing exposure or fitted value leaves a missing residual. Drop those observations before weighting.

# as_density_spec() says what it refuses more than one of

    Code
      expr
    Condition <propensity_density_error>
      Error:
      ! `.density` must be a single string.
      x It has length 2.
      i Supply "normal", "laplace", or "kernel".

# a refusal names the function the density was supplied to

    Code
      expr
    Condition <rlang_error>
      Error in `weight_from_density()`:
      ! `.density` must be one of "normal", "laplace", or "kernel", not "kernal".
      i Did you mean "kernel"?

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `weight_from_density()`:
      ! `.density` must return finite, non-negative values.
      x It returned 5 values that are negative or not finite.
      i It failed at standardized residuals from -0.213 to 1.16.
      i A density is non-negative everywhere. Check that `.density` returns a density rather than a log density or a distribution function.

---

    Code
      expr
    Condition <propensity_density_error>
      Error in `weight_from_density()`:
      ! `.density` must return finite, non-negative values.
      x It returned 2 values that are negative or not finite.
      i It failed at standardized residuals from -0.213 to -0.108.
      i A density is non-negative everywhere. Check that `.density` returns a density rather than a log density or a distribution function.

