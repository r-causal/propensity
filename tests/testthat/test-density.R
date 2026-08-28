# Tests for the density specifications that choose the conditional density
# family for continuous-exposure weights: the exported `dens_*()` constructors,
# the internal coercion `as_density_spec()`, the shared output check
# `check_density_values()`, and `density_eval()`, which evaluates a spec at the
# standardized residuals.
#
# The parametric families are measured against the stats functions they are
# named for, and the kernel family against a `stats::density()` fit and its
# interpolation written out by hand, so each test compares a spec against the
# density it claims to be rather than against itself.
#
# `as_density_spec()`, `check_density_values()`, and `density_eval()` are
# internal, and are called here by name the way the rest of the suite calls
# internals.

density_z <- function(n = 40, seed = 2024) {
  withr::with_seed(seed, stats::rnorm(n))
}

# The kernel evaluation as `density_eval()` is specified to perform it: one
# `stats::density()` fit on `fit_on`, bounded by `range`, interpolated to `z`.
kernel_oracle <- function(
  z,
  fit_on = z,
  range = range(z),
  bw = "nrd0",
  adjust = 1,
  kernel = "gaussian",
  n = 512
) {
  fit <- stats::density(
    fit_on,
    bw = bw,
    adjust = adjust,
    kernel = kernel,
    n = n,
    from = range[1],
    to = range[2]
  )

  stats::approxfun(fit$x, fit$y)(z)
}

# ---- constructors -----------------------------------------------------------

test_that("dens_normal() is the standard normal density", {
  z <- density_z()
  spec <- dens_normal()

  expect_s3_class(spec, "propensity_density")
  expect_type(spec, "list")
  expect_named(spec, c("family", "params", "fn"))
  expect_identical(spec$family, "normal")
  expect_identical(spec$params, list())
  expect_type(spec$fn, "closure")
  expect_equal(spec$fn(z), stats::dnorm(z))
})

test_that("dens_laplace() is the standard Laplace density", {
  z <- density_z()
  spec <- dens_laplace()

  expect_s3_class(spec, "propensity_density")
  expect_named(spec, c("family", "params", "fn"))
  expect_identical(spec$family, "laplace")
  expect_identical(spec$params, list())
  expect_equal(spec$fn(z), exp(-abs(z)) / 2)
})

test_that("dens_t() records its degrees of freedom", {
  z <- density_z()
  spec <- dens_t(df = 4)

  expect_s3_class(spec, "propensity_density")
  expect_named(spec, c("family", "params", "fn"))
  expect_identical(spec$family, "t")
  expect_named(spec$params, "df")
  expect_equal(spec$params$df, 4)
  expect_equal(spec$fn(z), stats::dt(z, df = 4))
  expect_equal(dens_t(df = 12)$fn(z), stats::dt(z, df = 12))
})

test_that("dens_kernel() records its stats::density() controls", {
  spec <- dens_kernel()

  expect_s3_class(spec, "propensity_density")
  expect_named(spec, c("family", "params", "fn"))
  expect_identical(spec$family, "kernel")
  expect_named(spec$params, c("bw", "adjust", "kernel", "n"))
  expect_identical(spec$params$bw, "nrd0")
  expect_equal(spec$params$adjust, 1)
  expect_identical(spec$params$kernel, "gaussian")
  expect_equal(spec$params$n, 512)
  expect_null(spec$fn)

  tuned <- dens_kernel(bw = 0.5, adjust = 1.5, kernel = "epanechnikov", n = 64)

  expect_equal(tuned$params$bw, 0.5)
  expect_equal(tuned$params$adjust, 1.5)
  expect_identical(tuned$params$kernel, "epanechnikov")
  expect_equal(tuned$params$n, 64)
})

test_that("dens_fn() carries the user's function", {
  z <- density_z()
  f <- function(r) stats::dt(r, df = 7)
  spec <- dens_fn(f)

  expect_s3_class(spec, "propensity_density")
  expect_named(spec, c("family", "params", "fn"))
  expect_identical(spec$family, "function")
  expect_identical(spec$params, list())
  expect_equal(spec$fn(z), stats::dt(z, df = 7))
})

# ---- constructor validation -------------------------------------------------

test_that("dens_t() requires a single positive finite df", {
  expect_error(dens_t(df = -1), class = "propensity_density_error")
  expect_error(dens_t(df = 0), class = "propensity_density_error")
  expect_error(dens_t(df = Inf), class = "propensity_density_error")
  expect_error(dens_t(df = NA_real_), class = "propensity_density_error")
  expect_error(dens_t(df = c(2, 4)), class = "propensity_density_error")
  expect_error(dens_t(df = "4"), class = "propensity_density_error")
  expect_error(dens_t(df = NULL), class = "propensity_density_error")
})

test_that("a propensity_density_error is a propensity_error", {
  expect_error(dens_t(df = -1), class = "propensity_error")
})

test_that("dens_kernel() requires a usable bandwidth", {
  expect_error(dens_kernel(bw = -1), class = "propensity_density_error")
  expect_error(dens_kernel(bw = 0), class = "propensity_density_error")
  expect_error(
    dens_kernel(bw = c(0.5, 1)),
    class = "propensity_density_error"
  )
  expect_error(dens_kernel(bw = NA_real_), class = "propensity_density_error")

  # A name that is not one of the rules is refused against the rules, however
  # the refusal is worded.
  expect_error(dens_kernel(bw = "silverman"), regexp = "nrd0")

  for (rule in c("nrd0", "nrd", "ucv", "bcv", "SJ", "SJ-ste", "SJ-dpi")) {
    expect_identical(dens_kernel(bw = rule)$params$bw, rule)
  }
})

test_that("dens_kernel() requires a single positive adjust", {
  expect_error(dens_kernel(adjust = -1), class = "propensity_density_error")
  expect_error(dens_kernel(adjust = 0), class = "propensity_density_error")
  expect_error(
    dens_kernel(adjust = c(1, 2)),
    class = "propensity_density_error"
  )
  expect_error(
    dens_kernel(adjust = NA_real_),
    class = "propensity_density_error"
  )
  expect_error(dens_kernel(adjust = "1"), class = "propensity_density_error")
})

test_that("dens_kernel() matches kernel against stats::density()", {
  kernels <- eval(formals(stats::density.default)$kernel)

  for (kernel in kernels) {
    expect_identical(dens_kernel(kernel = kernel)$params$kernel, kernel)
  }

  # A typo is corrected against the list of valid names rather than passed on
  # to `stats::density()`.
  expect_error(dens_kernel(kernel = "gausian"), regexp = "gaussian")
  expect_error(dens_kernel(kernel = "nonesuch"), regexp = "epanechnikov")
})

test_that("dens_kernel() requires an integer-valued n of at least 2", {
  expect_error(dens_kernel(n = 1), class = "propensity_density_error")
  expect_error(dens_kernel(n = 0), class = "propensity_density_error")
  expect_error(dens_kernel(n = -512), class = "propensity_density_error")
  expect_error(dens_kernel(n = 512.5), class = "propensity_density_error")
  expect_error(dens_kernel(n = c(256, 512)), class = "propensity_density_error")
  expect_error(dens_kernel(n = NA_integer_), class = "propensity_density_error")
  expect_error(dens_kernel(n = "512"), class = "propensity_density_error")

  expect_equal(dens_kernel(n = 2)$params$n, 2)
  expect_equal(dens_kernel(n = 1024L)$params$n, 1024)
})

test_that("dens_fn() requires a function of at least one argument", {
  expect_error(dens_fn(1), class = "propensity_density_error")
  expect_error(dens_fn("dnorm"), class = "propensity_density_error")
  expect_error(dens_fn(NULL), class = "propensity_density_error")
  expect_error(
    dens_fn(function() 1),
    class = "propensity_density_error"
  )
})

# ---- as_density_spec() ------------------------------------------------------

test_that("as_density_spec() maps strings to the matching spec", {
  z <- density_z()

  from_string <- as_density_spec("normal")
  expect_identical(from_string$family, "normal")
  expect_identical(from_string$params, list())
  expect_equal(from_string$fn(z), stats::dnorm(z))

  from_string <- as_density_spec("laplace")
  expect_identical(from_string$family, "laplace")
  expect_equal(from_string$fn(z), exp(-abs(z)) / 2)

  from_string <- as_density_spec("kernel")
  expect_identical(from_string$family, "kernel")
  expect_identical(from_string$params, dens_kernel()$params)
  expect_null(from_string$fn)
})

test_that("as_density_spec() corrects a mistyped string", {
  expect_error(as_density_spec("kernal"), regexp = "kernel")
  expect_error(as_density_spec("gaussian"), regexp = "normal")
  expect_error(as_density_spec("dt_4"), regexp = "laplace")
})

test_that("as_density_spec() passes a spec through unchanged", {
  for (spec in list(
    dens_normal(),
    dens_laplace(),
    dens_t(df = 4),
    dens_kernel(adjust = 1.5),
    dens_fn(function(r) stats::dnorm(r))
  )) {
    expect_identical(as_density_spec(spec), spec)
  }
})

test_that("as_density_spec() wraps a bare function", {
  z <- density_z()
  spec <- as_density_spec(function(r) stats::dt(r, df = 3))

  expect_s3_class(spec, "propensity_density")
  expect_identical(spec$family, "function")
  expect_equal(spec$fn(z), stats::dt(z, df = 3))

  # The wrap is `dens_fn()`, so it refuses what `dens_fn()` refuses.
  expect_error(
    as_density_spec(function() 1),
    class = "propensity_density_error"
  )
})

test_that("as_density_spec() refuses anything else", {
  expect_error(
    as_density_spec(1),
    class = "propensity_density_error"
  )
  expect_error(
    as_density_spec(NULL),
    class = "propensity_density_error"
  )
  expect_error(
    as_density_spec(NA),
    class = "propensity_density_error"
  )
  expect_error(
    as_density_spec(list(family = "normal")),
    class = "propensity_density_error"
  )
})

# ---- check_density_values() -------------------------------------------------

test_that("check_density_values() accepts a usable density", {
  expect_no_error(check_density_values(c(0.1, 0.2, 0.3), 3))
  # A zero is allowed as long as the density is not zero everywhere.
  expect_no_error(check_density_values(c(0, 0.2), 2))
})

test_that("check_density_values() refuses a malformed result", {
  expect_error(
    check_density_values("a", 1),
    class = "propensity_density_error"
  )
  expect_error(
    check_density_values(list(1, 2), 2),
    class = "propensity_density_error"
  )
  expect_error(
    check_density_values(c(0.1, NA, 0.3), 3),
    class = "propensity_density_error"
  )
  expect_error(
    check_density_values(c(0.1, 0.2), 3),
    class = "propensity_density_error"
  )
  expect_error(
    check_density_values(numeric(0), 3),
    class = "propensity_density_error"
  )
})

test_that("check_density_values() refuses non-finite or negative values", {
  expect_error(
    check_density_values(c(0.1, Inf), 2),
    class = "propensity_density_error"
  )
  expect_error(
    check_density_values(c(0.1, -Inf), 2),
    class = "propensity_density_error"
  )
  expect_error(
    check_density_values(c(0.1, -0.2), 2),
    class = "propensity_density_error"
  )
})

test_that("check_density_values() refuses an everywhere-zero density", {
  expect_error(
    check_density_values(c(0, 0, 0), 3),
    class = "propensity_density_error"
  )
})

test_that("check_density_values() names the argument it is checking", {
  expect_error(
    check_density_values(c(-1, 1), 2, what = ".numerator"),
    regexp = "numerator"
  )
})

# ---- density_eval() ---------------------------------------------------------

test_that("density_eval() reproduces the parametric densities", {
  z <- density_z()

  expect_equal(
    density_eval(dens_normal(), z),
    stats::dnorm(z)
  )
  expect_equal(
    density_eval(dens_laplace(), z),
    exp(-abs(z)) / 2
  )
  expect_equal(
    density_eval(dens_t(df = 4), z),
    stats::dt(z, df = 4)
  )
  expect_equal(
    density_eval(dens_fn(function(r) stats::dt(r, df = 9)), z),
    stats::dt(z, df = 9)
  )
})

test_that("density_eval() ignores fit_on and range for parametric families", {
  z <- density_z()

  expect_equal(
    density_eval(
      dens_normal(),
      z,
      fit_on = density_z(seed = 99),
      range = c(-10, 10)
    ),
    stats::dnorm(z)
  )
})

test_that("density_eval() fits the kernel with the spec's controls", {
  z <- density_z(n = 200)

  expect_equal(
    density_eval(dens_kernel(), z),
    kernel_oracle(z)
  )
  expect_equal(
    density_eval(dens_kernel(adjust = 1.5), z),
    kernel_oracle(z, adjust = 1.5)
  )
  expect_equal(
    density_eval(dens_kernel(bw = 0.4, n = 64), z),
    kernel_oracle(z, bw = 0.4, n = 64)
  )
  expect_equal(
    density_eval(dens_kernel(kernel = "epanechnikov"), z),
    kernel_oracle(z, kernel = "epanechnikov")
  )

  # The controls are not decoration: a wider bandwidth is a different density.
  expect_false(isTRUE(all.equal(
    density_eval(dens_kernel(adjust = 3), z),
    density_eval(dens_kernel(), z)
  )))
})

test_that("density_eval() fits the kernel on fit_on over the given range", {
  fit_on <- density_z(n = 200)
  grid <- seq(min(fit_on), max(fit_on), length.out = 50)
  range <- range(c(fit_on, grid))

  values <- density_eval(
    dens_kernel(),
    grid,
    fit_on = fit_on,
    range = range
  )

  expect_equal(values, kernel_oracle(grid, fit_on = fit_on, range = range))
  expect_length(values, length(grid))
  expect_false(anyNA(values))
  expect_true(all(is.finite(values)))
  expect_true(all(values >= 0))
})

test_that("density_eval() covers the whole supplied range without NA", {
  fit_on <- density_z(n = 200)
  range <- c(min(fit_on) - 1, max(fit_on) + 1)
  grid <- seq(range[1], range[2], length.out = 101)

  values <- density_eval(
    dens_kernel(),
    grid,
    fit_on = fit_on,
    range = range
  )

  expect_length(values, length(grid))
  expect_false(anyNA(values))
  expect_true(all(is.finite(values)))
})

test_that("density_eval() refuses a point outside the kernel's range", {
  fit_on <- density_z(n = 200)

  expect_error(
    density_eval(
      dens_kernel(),
      c(fit_on, max(fit_on) + 1),
      fit_on = fit_on,
      range = range(fit_on)
    ),
    class = "propensity_density_error"
  )
})

test_that("density_eval() checks what a user function returns", {
  z <- density_z(n = 10)

  expect_error(
    density_eval(dens_fn(function(r) -stats::dnorm(r)), z),
    class = "propensity_density_error"
  )
  expect_error(
    density_eval(dens_fn(function(r) rep(0, length(r))), z),
    class = "propensity_density_error"
  )
  expect_error(
    density_eval(dens_fn(function(r) 1), z),
    class = "propensity_density_error"
  )
  expect_error(
    density_eval(dens_fn(function(r) rep(NA_real_, length(r))), z),
    class = "propensity_density_error"
  )
  expect_error(
    density_eval(dens_fn(function(r) rep(Inf, length(r))), z),
    class = "propensity_density_error"
  )
  expect_error(
    density_eval(
      dens_fn(function(r) as.character(stats::dnorm(r))),
      z
    ),
    class = "propensity_density_error"
  )
})

# ---- print and format -------------------------------------------------------

# The printed forms are pinned as exact text here rather than by snapshot: the
# constructors do not exist yet in this round, so a snapshot would record the
# absence rather than the format.

test_that("format() writes the family and its parameters", {
  expect_identical(format(dens_normal()), "normal")
  expect_identical(format(dens_laplace()), "laplace")
  expect_identical(format(dens_t(df = 4)), "t(df = 4)")
  expect_identical(
    format(dens_kernel()),
    'kernel(bw = "nrd0", adjust = 1, kernel = "gaussian", n = 512)'
  )
  expect_identical(
    format(dens_kernel(bw = 0.5, adjust = 1.5, kernel = "cosine", n = 64)),
    'kernel(bw = 0.5, adjust = 1.5, kernel = "cosine", n = 64)'
  )
  expect_identical(format(dens_fn(function(r) stats::dnorm(r))), "function")
})

test_that("print() wraps the format in angle brackets", {
  expect_output(print(dens_normal()), "<density: normal>", fixed = TRUE)
  expect_output(print(dens_laplace()), "<density: laplace>", fixed = TRUE)
  expect_output(print(dens_t(df = 4)), "<density: t(df = 4)>", fixed = TRUE)
  expect_output(
    print(dens_kernel()),
    '<density: kernel(bw = "nrd0", adjust = 1, kernel = "gaussian", n = 512)>',
    fixed = TRUE
  )
  expect_output(
    print(dens_fn(function(r) stats::dnorm(r))),
    "<density: function>",
    fixed = TRUE
  )
})
