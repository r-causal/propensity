# The conditional density family of continuous-exposure weights: the `.density`
# argument of `wt_ate()` and `wt_cens()`, the arithmetic each family produces,
# and what the resulting weights record about it.
#
# fixtures/continuous-weights-normal.rds holds continuous ATE weights as they
# were before `.density` existed, over three seeds and the six combinations of
# stabilization and residual spread. It was written by
# scratch/make-continuous-weights-fixture.R against the package as it stood
# then, and it is what holds the normal family, which is the default, to numbers
# it did not choose for itself. Regenerating it against a later version of the
# package would defeat that, so the file is checked in and the script is not
# rerun.

# One continuous-exposure problem the hand calculations are written against: an
# exposure, the fitted conditional means of a model for it, and a spread for
# each observation.
continuous_density_data <- local({
  set.seed(20250828)

  n <- 60
  x <- rnorm(n)
  exposure <- 0.6 * x + rnorm(n)

  list(
    n = n,
    x = x,
    exposure = exposure,
    mu = as.numeric(fitted(lm(exposure ~ x))),
    sigma_i = runif(n, 0.7, 1.3)
  )
})

# Weights for that problem, stabilized on the marginal density unless the test
# asks for something else.
continuous_density_wt <- function(..., stabilize = TRUE) {
  wt_ate(
    continuous_density_data$mu,
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = stabilize,
    ...
  )
}

# The pooled residual spread: the uncentered root mean square of the residuals,
# which is what the conditional density is spread by when no `.sigma` is given.
continuous_density_pooled <- function() {
  sqrt(mean(
    (continuous_density_data$exposure - continuous_density_data$mu)^2
  ))
}

# The standardized residual the conditional density is evaluated at.
continuous_density_z <- function(sigma = continuous_density_pooled()) {
  (continuous_density_data$exposure - continuous_density_data$mu) / sigma
}

# The standardized marginal value the numerator density is evaluated at, and the
# population standard deviation it is standardized by.
continuous_density_sd_a <- function() {
  exposure <- continuous_density_data$exposure

  sqrt(mean((exposure - mean(exposure))^2))
}

continuous_density_z_a <- function() {
  exposure <- continuous_density_data$exposure

  (exposure - mean(exposure)) / continuous_density_sd_a()
}

# A kernel density estimate of `z`, fit over the range of the values it is fit
# on and read back at `z` by linear interpolation: what `stats::density()` and
# `stats::approxfun()` do when they are called by hand.
continuous_density_kde <- function(
  z,
  fit_on = z,
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
    from = min(fit_on),
    to = max(fit_on)
  )

  stats::approxfun(fit$x, fit$y)(z)
}

# ---- the fixture ------------------------------------------------------------

test_that("continuous weights still equal the recorded normal weights", {
  fixture <- readRDS(test_path("fixtures", "continuous-weights-normal.rds"))

  for (problem in fixture) {
    for (case in problem$cases) {
      # Some of the recorded cases were built with a spread small enough to put
      # the weights past the finite-variance boundary. The recorded values are
      # what this test is about, so the report that says so is muffled.
      weights <- muffle_variance_warning(wt_ate(
        problem$mu,
        problem$exposure,
        .sigma = case$sigma,
        exposure_type = "continuous",
        stabilize = case$stabilize
      ))

      expect_equal(
        as.numeric(weights),
        case$weights,
        tolerance = 1e-12,
        info = paste(problem$seed, case$stabilize, case$sigma_kind)
      )
    }
  }
})

test_that("the normal family reproduces the recorded weights exactly", {
  fixture <- readRDS(test_path("fixtures", "continuous-weights-normal.rds"))

  for (problem in fixture) {
    for (case in problem$cases) {
      by_string <- muffle_variance_warning(wt_ate(
        problem$mu,
        problem$exposure,
        .sigma = case$sigma,
        exposure_type = "continuous",
        stabilize = case$stabilize,
        .density = "normal"
      ))

      by_spec <- muffle_variance_warning(wt_ate(
        problem$mu,
        problem$exposure,
        .sigma = case$sigma,
        exposure_type = "continuous",
        stabilize = case$stabilize,
        .density = dens_normal()
      ))

      label <- paste(problem$seed, case$stabilize, case$sigma_kind)

      expect_equal(
        as.numeric(by_string),
        case$weights,
        tolerance = 1e-12,
        info = label
      )
      expect_equal(
        as.numeric(by_spec),
        case$weights,
        tolerance = 1e-12,
        info = label
      )
    }
  }
})

# ---- hand calculations ------------------------------------------------------

test_that("Student's t weights are the t density ratio computed by hand", {
  z <- continuous_density_z()
  sigma <- continuous_density_pooled()
  sd_a <- continuous_density_sd_a()

  f_den <- stats::dt(z, df = 4) / sigma
  f_num <- stats::dt(continuous_density_z_a(), df = 4) / sd_a

  stabilized <- continuous_density_wt(.density = dens_t(df = 4))
  expect_equal(as.numeric(stabilized), f_num / f_den, tolerance = 1e-12)

  unstabilized <- continuous_density_wt(
    .density = dens_t(df = 4),
    stabilize = FALSE
  )
  expect_equal(as.numeric(unstabilized), 1 / f_den, tolerance = 1e-12)

  # A t with many degrees of freedom is nearly the normal, and the weights it
  # gives are nearly the weights the default gives.
  heavy <- as.numeric(continuous_density_wt(.density = dens_t(df = 4)))
  light <- as.numeric(continuous_density_wt(.density = dens_t(df = 500)))
  normal <- as.numeric(continuous_density_wt())

  expect_equal(light, normal, tolerance = 1e-2)
  expect_false(isTRUE(all.equal(heavy, normal, tolerance = 1e-2)))
})

test_that("Laplace weights are the Laplace density ratio computed by hand", {
  z <- continuous_density_z()
  sigma <- continuous_density_pooled()
  sd_a <- continuous_density_sd_a()

  f_den <- exp(-abs(z)) / 2 / sigma
  f_num <- exp(-abs(continuous_density_z_a())) / 2 / sd_a

  by_string <- continuous_density_wt(.density = "laplace")
  by_spec <- continuous_density_wt(.density = dens_laplace())

  expect_equal(as.numeric(by_string), f_num / f_den, tolerance = 1e-12)
  expect_equal(as.numeric(by_spec), f_num / f_den, tolerance = 1e-12)

  unstabilized <- continuous_density_wt(
    .density = "laplace",
    stabilize = FALSE
  )
  expect_equal(as.numeric(unstabilized), 1 / f_den, tolerance = 1e-12)
})

test_that("a density the user writes is evaluated as written", {
  z <- continuous_density_z()
  sigma <- continuous_density_pooled()
  sd_a <- continuous_density_sd_a()

  # Deliberately not a family the package names, so nothing but the function
  # itself can produce these numbers.
  f <- function(z) stats::dlogis(z)

  f_den <- f(z) / sigma
  f_num <- f(continuous_density_z_a()) / sd_a

  bare <- continuous_density_wt(.density = f)
  wrapped <- continuous_density_wt(.density = dens_fn(f))

  expect_equal(as.numeric(bare), f_num / f_den, tolerance = 1e-12)
  expect_equal(as.numeric(wrapped), f_num / f_den, tolerance = 1e-12)

  # A function that writes out a named family gives that family's weights.
  as_function <- continuous_density_wt(
    .density = function(z) stats::dt(z, df = 4)
  )
  as_spec <- continuous_density_wt(.density = dens_t(df = 4))

  expect_equal(
    as.numeric(as_function),
    as.numeric(as_spec),
    tolerance = 1e-12
  )
})

test_that("the kernel denominator is the estimate stats::density gives", {
  z <- continuous_density_z()
  sigma <- continuous_density_pooled()

  f_den <- continuous_density_kde(z) / sigma

  unstabilized <- continuous_density_wt(
    .density = "kernel",
    stabilize = FALSE
  )
  expect_equal(as.numeric(unstabilized), 1 / f_den, tolerance = 1e-12)

  # The `stats::density()` controls reach the estimate: a bandwidth multiplied
  # by 1.5, a different kernel, and a coarser grid all show up in the weights.
  smoother <- continuous_density_wt(
    .density = dens_kernel(adjust = 1.5, kernel = "epanechnikov", n = 128),
    stabilize = FALSE
  )
  smoother_den <- continuous_density_kde(
    z,
    adjust = 1.5,
    kernel = "epanechnikov",
    n = 128
  ) /
    sigma

  expect_equal(as.numeric(smoother), 1 / smoother_den, tolerance = 1e-12)
  expect_false(
    isTRUE(all.equal(as.numeric(smoother), as.numeric(unstabilized)))
  )
})

test_that("the kernel numerator is fit on the standardized marginal values", {
  z <- continuous_density_z()
  z_a <- continuous_density_z_a()
  sigma <- continuous_density_pooled()
  sd_a <- continuous_density_sd_a()

  f_den <- continuous_density_kde(z) / sigma
  f_num <- continuous_density_kde(z_a) / sd_a

  stabilized <- continuous_density_wt(.density = "kernel")

  expect_equal(as.numeric(stabilized), f_num / f_den, tolerance = 1e-12)
})

test_that("a supplied spread scales the density it is given to", {
  scalar <- 0.85
  per_observation <- continuous_density_data$sigma_i
  sd_a <- continuous_density_sd_a()
  f_num <- stats::dt(continuous_density_z_a(), df = 4) / sd_a

  z_scalar <- continuous_density_z(scalar)
  f_den_scalar <- stats::dt(z_scalar, df = 4) / scalar

  expect_equal(
    as.numeric(continuous_density_wt(
      .sigma = scalar,
      .density = dens_t(df = 4)
    )),
    f_num / f_den_scalar,
    tolerance = 1e-12
  )

  z_each <- continuous_density_z(per_observation)
  f_den_each <- stats::dt(z_each, df = 4) / per_observation

  expect_equal(
    as.numeric(continuous_density_wt(
      .sigma = per_observation,
      .density = dens_t(df = 4)
    )),
    f_num / f_den_each,
    tolerance = 1e-12
  )

  # The numerator is the marginal density of the exposure, which a conditional
  # spread has nothing to say about, so only the denominator moves.
  expect_equal(
    as.numeric(continuous_density_wt(
      .sigma = per_observation,
      .density = dens_t(df = 4),
      stabilize = FALSE
    )),
    1 / f_den_each,
    tolerance = 1e-12
  )
})

test_that("a stabilization score replaces the numerator of any family", {
  score <- seq(0.4, 1.6, length.out = continuous_density_data$n)
  sigma <- continuous_density_pooled()
  f_den <- stats::dt(continuous_density_z(), df = 4) / sigma

  weights <- continuous_density_wt(
    .density = dens_t(df = 4),
    stabilization_score = score
  )

  expect_equal(as.numeric(weights), score / f_den, tolerance = 1e-12)
})

# ---- refusals ---------------------------------------------------------------

test_that("a density other than the default refuses a binary exposure", {
  set.seed(11)
  ps <- runif(20, 0.2, 0.8)
  trt <- rbinom(20, 1, ps)

  expect_error(
    wt_ate(ps, trt, exposure_type = "binary", .density = dens_t(df = 4)),
    class = "propensity_density_error"
  )
  expect_error(
    wt_ate(ps, trt, exposure_type = "binary", .density = "laplace"),
    class = "propensity_density_error"
  )
  expect_error(
    wt_cens(ps, trt, exposure_type = "binary", .density = "kernel"),
    class = "propensity_density_error"
  )

  # The message says what the exposure is being treated as, the way
  # `check_sigma()` does for `.sigma`.
  binary_error <- rlang::catch_cnd(
    wt_ate(ps, trt, exposure_type = "binary", .density = "laplace")
  )
  expect_match(conditionMessage(binary_error), "binary")
})

test_that("a density other than the default refuses a categorical exposure", {
  exposure <- factor(rep(c("a", "b", "c"), each = 4))
  ps <- matrix(
    rep(c(0.5, 0.3, 0.2), times = 12),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )

  expect_error(
    wt_ate(
      ps,
      exposure,
      exposure_type = "categorical",
      .density = dens_t(df = 4)
    ),
    class = "propensity_density_error"
  )

  categorical_error <- rlang::catch_cnd(
    wt_ate(ps, exposure, exposure_type = "categorical", .density = "kernel")
  )
  expect_match(conditionMessage(categorical_error), "categorical")
})

test_that("the default density is accepted for every exposure type", {
  set.seed(12)
  ps <- runif(20, 0.2, 0.8)
  trt <- rbinom(20, 1, ps)

  plain <- wt_ate(ps, trt, exposure_type = "binary")

  expect_equal(
    as.numeric(wt_ate(ps, trt, exposure_type = "binary", .density = "normal")),
    as.numeric(plain)
  )
  expect_equal(
    as.numeric(wt_ate(
      ps,
      trt,
      exposure_type = "binary",
      .density = dens_normal()
    )),
    as.numeric(plain)
  )

  # Weights that are not a ratio of densities record none, whether or not the
  # default was written out.
  expect_null(density_meta(
    wt_ate(ps, trt, exposure_type = "binary", .density = "normal")
  ))
})

test_that("the default density is accepted for a categorical exposure", {
  exposure <- factor(rep(c("a", "b", "c"), each = 4))
  ps <- matrix(
    rep(c(0.5, 0.3, 0.2), times = 12),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )

  plain <- wt_ate(ps, exposure, exposure_type = "categorical")

  expect_equal(
    as.numeric(wt_ate(
      ps,
      exposure,
      exposure_type = "categorical",
      .density = "normal"
    )),
    as.numeric(plain)
  )
  expect_equal(
    as.numeric(wt_ate(
      ps,
      exposure,
      exposure_type = "categorical",
      .density = dens_normal()
    )),
    as.numeric(plain)
  )

  expect_null(density_meta(
    wt_ate(ps, exposure, exposure_type = "categorical", .density = "normal")
  ))
})

test_that("a density other than the default refuses a detected binary exposure", {
  set.seed(15)
  ps <- runif(20, 0.2, 0.8)
  trt <- rbinom(20, 1, ps)

  # The type was not named, so the refusal is against the type detection
  # resolved the 0/1 exposure to rather than against anything the caller wrote.
  expect_error(
    wt_ate(ps, trt, .density = dens_t(df = 4)),
    class = "propensity_density_error"
  )
  expect_error(
    wt_cens(ps, trt, .density = "kernel"),
    class = "propensity_density_error"
  )

  detected <- rlang::catch_cnd(wt_ate(ps, trt, .density = "laplace"))
  expect_match(conditionMessage(detected), "binary")
})

test_that("the refusal of a density reads the same way for either type", {
  set.seed(16)
  ps <- runif(20, 0.2, 0.8)
  trt <- rbinom(20, 1, ps)

  categorical_exposure <- factor(rep(c("a", "b", "c"), each = 4))
  categorical_ps <- matrix(
    rep(c(0.5, 0.3, 0.2), times = 12),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )

  expect_propensity_error(
    wt_ate(ps, trt, exposure_type = "binary", .density = "laplace")
  )
  expect_propensity_error(
    wt_ate(
      categorical_ps,
      categorical_exposure,
      exposure_type = "categorical",
      .density = dens_t(df = 4)
    )
  )
  expect_propensity_error(wt_cens(ps, trt, .density = "kernel"))
})

test_that("a density that is not a density is refused", {
  expect_error(
    continuous_density_wt(.density = 1),
    class = "propensity_density_error"
  )
  expect_error(
    continuous_density_wt(.density = c("normal", "laplace")),
    class = "propensity_density_error"
  )
  expect_error(
    continuous_density_wt(.density = function() 1),
    class = "propensity_density_error"
  )

  # The output checks apply to a density however it was written.
  expect_error(
    continuous_density_wt(.density = function(z) 1),
    class = "propensity_density_error"
  )
  expect_error(
    continuous_density_wt(.density = function(z) -stats::dnorm(z)),
    class = "propensity_density_error"
  )

  # A mistyped family name is corrected by the argument matcher rather than
  # refused on the package's own terms.
  expect_error(
    continuous_density_wt(.density = "gaussian"),
    regexp = "must be one of"
  )
})

test_that("the density cannot be supplied by position", {
  # `.density` sits after the dots, so a value meant for it that is written
  # without a name is an unused argument rather than a silent collision with
  # `.sigma` or `exposure_type`.
  expect_error(
    wt_ate(
      continuous_density_data$mu,
      continuous_density_data$exposure,
      NULL,
      "continuous",
      NULL,
      NULL,
      TRUE,
      NULL,
      dens_t(df = 4)
    ),
    class = "rlib_error_dots_nonempty"
  )
})

# A multi-response linear model of the fixture's exposure. Its fitted values are
# a matrix of conditional means, one column for each response, rather than the
# vector of conditional means a continuous exposure is weighted from.
continuous_multi_response_fit <- function() {
  fit_data <- data.frame(
    x = continuous_density_data$x,
    dose = continuous_density_data$exposure,
    other = continuous_density_data$exposure - 0.5 * continuous_density_data$x
  )

  stats::lm(cbind(dose, other) ~ x, data = fit_data)
}

test_that("a continuous exposure refuses a matrix of conditional means", {
  # A continuous exposure has one conditional mean for each unit. A matrix with
  # as many rows as the exposure is long otherwise passes the length check, and
  # the conditional density is then evaluated over every cell of it, giving one
  # weight for each entry of the matrix rather than one for each unit.
  fit <- continuous_multi_response_fit()

  expect_error(
    wt_ate(
      stats::fitted(fit),
      continuous_density_data$exposure,
      exposure_type = "continuous"
    ),
    class = "propensity_ps_shape_error"
  )

  # The same conditional mean written as a vector is still weighted.
  expect_length(
    continuous_density_wt(),
    continuous_density_data$n
  )
})

test_that("a continuous exposure reads one column of conditional means", {
  # A one-column matrix and a one-dimensional array each hold one conditional
  # mean for each unit, which is the shape the weights are built from, so each
  # is read as the vector it is rather than refused for carrying a dimension.
  # Both weighted correctly before the shape guard was written, and what the
  # guard is for is a matrix holding more than one mean for each unit.
  oracle <- as.numeric(continuous_density_wt())

  one_column <- wt_ate(
    matrix(continuous_density_data$mu, ncol = 1),
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  expect_s3_class(one_column, "psw")
  expect_equal(as.numeric(one_column), oracle)

  one_dimension <- wt_ate(
    array(continuous_density_data$mu),
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  expect_s3_class(one_dimension, "psw")
  expect_equal(as.numeric(one_dimension), oracle)
})

test_that("a multi-response linear model is refused by the model route", {
  fit <- continuous_multi_response_fit()

  expect_error(
    wt_ate(fit, continuous_density_data$exposure, exposure_type = "continuous"),
    class = "propensity_ps_shape_error"
  )

  # Read without an exposure, the model route takes the fit's own response,
  # which is the response matrix. The refusal has to come before the lengths
  # are compared, so that the account is of the multi-response fit rather than
  # of a length the caller never wrote.
  expect_error(
    wt_ate(fit, exposure_type = "continuous"),
    class = "propensity_ps_shape_error"
  )
  expect_error(wt_ate(fit), class = "propensity_ps_shape_error")
})

test_that("the refusal of a matrix of conditional means reads plainly", {
  fit <- continuous_multi_response_fit()

  expect_propensity_error(
    wt_ate(
      stats::fitted(fit),
      continuous_density_data$exposure,
      exposure_type = "continuous"
    )
  )

  # Censoring weights reach the continuous formula through the same route, so
  # the refusal is the same one, named for the function the user called.
  expect_propensity_error(
    wt_cens(
      stats::fitted(fit),
      continuous_density_data$exposure,
      exposure_type = "continuous"
    )
  )

  expect_propensity_error(
    wt_ate(fit, continuous_density_data$exposure, exposure_type = "continuous")
  )

  expect_propensity_error(
    wt_cens(fit, continuous_density_data$exposure, exposure_type = "continuous")
  )
})

# ---- what the weights record ------------------------------------------------

test_that("the weights record the density family they were built from", {
  families <- list(
    normal = list(input = "normal", recorded = "normal"),
    laplace = list(input = dens_laplace(), recorded = "laplace"),
    t = list(input = dens_t(df = 4), recorded = "t(df = 4)"),
    kernel = list(
      input = dens_kernel(adjust = 1.5),
      recorded = 'kernel(bw = "nrd0", adjust = 1.5, kernel = "gaussian", n = 512)'
    ),
    fn = list(input = function(z) stats::dlogis(z), recorded = "function")
  )

  for (family in families) {
    weights <- continuous_density_wt(.density = family$input)
    record <- density_meta(weights)

    expect_s3_class(record, "propensity_density_meta")
    expect_s3_class(record$density, "propensity_density")
    expect_identical(format(record$density), family$recorded)
    expect_identical(exposure_type(weights), "continuous")

    # The rest of the record is what it was: the family says nothing about
    # what stabilized the weights or where the spread came from.
    expect_identical(record$numerator, "marginal")
    expect_identical(record$sigma, "pooled")
  }
})

test_that("the density record travels with the numerator and the spread", {
  score <- seq(0.4, 1.6, length.out = continuous_density_data$n)

  scored <- density_meta(continuous_density_wt(
    .density = dens_t(df = 4),
    stabilization_score = score
  ))
  expect_identical(scored$numerator, "score")
  expect_identical(format(scored$density), "t(df = 4)")

  unstabilized <- density_meta(continuous_density_wt(
    .density = "laplace",
    stabilize = FALSE
  ))
  expect_identical(unstabilized$numerator, "none")
  expect_identical(format(unstabilized$density), "laplace")

  supplied <- density_meta(continuous_density_wt(
    .density = "kernel",
    .sigma = 0.85
  ))
  expect_identical(supplied$sigma, "supplied")
  expect_identical(supplied$density$family, "kernel")
})

test_that("a density the user wrote is recorded as the object they wrote", {
  f <- function(z) stats::dlogis(z)
  record <- density_meta(continuous_density_wt(.density = f))

  expect_identical(record$density$family, "function")
  expect_identical(record$density$fn, f)
})

test_that("weights built from the same user density combine", {
  f <- function(z) stats::dlogis(z)

  first <- continuous_density_wt(.density = f)
  second <- continuous_density_wt(.density = f)

  combined <- expect_silent(c(first, second))

  expect_length(combined, 2 * continuous_density_data$n)
  expect_identical(exposure_type(combined), "continuous")
  expect_identical(format(density_meta(combined)$density), "function")
  expect_identical(density_meta(combined)$density$fn, f)
})

test_that("weights built from different densities drop the density record", {
  four <- continuous_density_wt(.density = function(z) stats::dt(z, df = 4))
  six <- continuous_density_wt(.density = function(z) stats::dt(z, df = 6))

  combined <- NULL
  expect_warning(
    combined <- c(four, six),
    class = "propensity_metadata_conflict_warning"
  )

  expect_s3_class(combined, "psw")
  expect_null(density_meta(combined))
  expect_identical(exposure_type(combined), "continuous")

  # Two named families disagree the same way.
  mixed <- NULL
  expect_warning(
    mixed <- c(
      continuous_density_wt(.density = dens_t(df = 4)),
      continuous_density_wt(.density = "laplace")
    ),
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(density_meta(mixed))
})

# ---- the other methods ------------------------------------------------------

test_that("the density reaches the weights through a fitted model", {
  exposure <- continuous_density_data$exposure
  x <- continuous_density_data$x
  model_data <- data.frame(exposure = exposure, x = x)
  fit <- glm(exposure ~ x, data = model_data, family = gaussian())

  from_model <- wt_ate(
    fit,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_t(df = 4)
  )
  from_numeric <- wt_ate(
    as.numeric(fitted(fit)),
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_t(df = 4)
  )

  expect_equal(
    as.numeric(from_model),
    as.numeric(from_numeric),
    tolerance = 1e-12
  )
  expect_identical(format(density_meta(from_model)$density), "t(df = 4)")
})

test_that("the density reaches the weights through a data frame", {
  scores <- data.frame(mu = continuous_density_data$mu)

  from_frame <- wt_ate(
    scores,
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "laplace"
  )

  expect_equal(
    as.numeric(from_frame),
    as.numeric(continuous_density_wt(.density = "laplace")),
    tolerance = 1e-12
  )
  expect_identical(format(density_meta(from_frame)$density), "laplace")
})

test_that("the density reaches the weights through a modified score", {
  set.seed(13)
  n <- 40
  x <- rnorm(n)
  exposure <- 0.5 * x + rnorm(n)
  model_data <- data.frame(exposure = exposure, x = x)
  fit <- glm(exposure ~ x, data = model_data, family = gaussian())

  # A score a modification can be applied to has to lie in (0, 1); it stands in
  # for the fitted conditional mean here.
  scores <- plogis(0.5 * x)

  refit <- ps_refit(
    ps_trim(scores, method = "ps", lower = 0.2, upper = 0.8),
    model = fit
  )
  truncated <- ps_trunc(scores, method = "ps", lower = 0.2, upper = 0.8)
  calibrated <- ps_calibrate(scores, rbinom(n, 1, scores))

  for (modified in list(refit, truncated, calibrated)) {
    weights <- wt_ate(
      modified,
      exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = dens_t(df = 4)
    )
    plain <- wt_ate(
      as.numeric(modified),
      exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = dens_t(df = 4)
    )

    expect_equal(as.numeric(weights), as.numeric(plain), tolerance = 1e-12)
    expect_identical(format(density_meta(weights)$density), "t(df = 4)")
  }
})

test_that("censoring weights take a density of their own", {
  weights <- wt_cens(
    continuous_density_data$mu,
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_t(df = 4)
  )

  expect_equal(
    as.numeric(weights),
    as.numeric(continuous_density_wt(.density = dens_t(df = 4))),
    tolerance = 1e-12
  )
  expect_identical(estimand(weights), "uncensored")
  expect_identical(format(density_meta(weights)$density), "t(df = 4)")

  # And through the model method, which routes to the numeric one.
  model_data <- data.frame(
    exposure = continuous_density_data$exposure,
    x = continuous_density_data$x
  )
  fit <- glm(exposure ~ x, data = model_data, family = gaussian())

  from_model <- wt_cens(
    fit,
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "laplace"
  )

  expect_identical(estimand(from_model), "uncensored")
  expect_identical(format(density_meta(from_model)$density), "laplace")
})

test_that("censoring weights carry a density through their other methods", {
  set.seed(17)
  n <- 40
  x <- rnorm(n)
  exposure <- 0.5 * x + rnorm(n)

  from_frame <- wt_cens(
    data.frame(mu = continuous_density_data$mu),
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "laplace"
  )

  expect_equal(
    as.numeric(from_frame),
    as.numeric(continuous_density_wt(.density = "laplace")),
    tolerance = 1e-12
  )
  expect_identical(estimand(from_frame), "uncensored")
  expect_identical(format(density_meta(from_frame)$density), "laplace")

  # A score a modification can be applied to has to lie in (0, 1); it stands in
  # for the fitted conditional mean here.
  scores <- plogis(0.5 * x)

  modified_scores <- list(
    ps_refit(
      ps_trim(scores, method = "ps", lower = 0.2, upper = 0.8),
      model = glm(
        exposure ~ x,
        data = data.frame(exposure = exposure, x = x),
        family = gaussian()
      )
    ),
    ps_trunc(scores, method = "ps", lower = 0.2, upper = 0.8),
    ps_calibrate(scores, rbinom(n, 1, scores))
  )

  for (modified in modified_scores) {
    weights <- wt_cens(
      modified,
      exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = dens_t(df = 4)
    )
    plain <- wt_cens(
      as.numeric(modified),
      exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = dens_t(df = 4)
    )

    expect_equal(as.numeric(weights), as.numeric(plain), tolerance = 1e-12)
    expect_identical(format(density_meta(weights)$density), "t(df = 4)")
  }
})

# ---- missing values ---------------------------------------------------------

# A unit with no standardized residual has no weight, and the density is asked
# only about the units that have one. That is what lets a kernel be fit at all
# when some of the residuals are missing, and it is why a family is never the
# reason a weight came back missing.

test_that("a missing fitted value leaves a weight the density was not asked about", {
  exposure <- continuous_density_data$exposure
  mu <- continuous_density_data$mu
  missing_at <- c(4L, 37L)
  mu[missing_at] <- NA_real_

  weights <- as.numeric(wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "kernel"
  ))

  expect_identical(which(is.na(weights)), missing_at)

  # The pooled spread, the kernel, and the range it is fit over all read the
  # residuals that exist and nothing else.
  present <- setdiff(seq_len(continuous_density_data$n), missing_at)
  sigma <- sqrt(mean((exposure[present] - mu[present])^2))
  z <- (exposure[present] - mu[present]) / sigma
  f_den <- continuous_density_kde(z) / sigma

  # The exposure is whole, so the marginal density is still fit on all of it.
  f_num <- continuous_density_kde(continuous_density_z_a()) /
    continuous_density_sd_a()

  expect_equal(weights[present], f_num[present] / f_den, tolerance = 1e-12)
})

test_that("a missing exposure is left out of the numerator as well", {
  exposure <- continuous_density_data$exposure
  mu <- continuous_density_data$mu
  missing_at <- 11L
  exposure[missing_at] <- NA_real_

  weights <- as.numeric(wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "kernel"
  ))

  expect_identical(which(is.na(weights)), missing_at)

  present <- setdiff(seq_len(continuous_density_data$n), missing_at)
  sigma <- sqrt(mean((exposure[present] - mu[present])^2))
  z <- (exposure[present] - mu[present]) / sigma
  f_den <- continuous_density_kde(z) / sigma

  # A missing exposure is missing from its own marginal density too, so the
  # moments it is read at and the kernel fit on it both skip that unit.
  mu_a <- mean(exposure[present])
  sd_a <- sqrt(mean((exposure[present] - mu_a)^2))
  f_num <- continuous_density_kde((exposure[present] - mu_a) / sd_a) / sd_a

  expect_equal(weights[present], f_num / f_den, tolerance = 1e-12)
})

test_that("a trimmed propensity score leaves the units it set aside missing", {
  set.seed(18)
  n <- 40
  x <- rnorm(n)
  exposure <- 0.5 * x + rnorm(n)

  # A score a modification can be applied to has to lie in (0, 1); it stands in
  # for the fitted conditional mean here. Trimmed and not refit, it carries a
  # missing value at every unit it set aside, which is the ordinary way a
  # standardized residual goes missing.
  trimmed <- ps_trim(plogis(0.9 * x), method = "ps", lower = 0.2, upper = 0.8)
  mu <- as.numeric(trimmed)
  missing_at <- which(is.na(mu))
  expect_gt(length(missing_at), 1)

  present <- setdiff(seq_len(n), missing_at)
  sigma <- sqrt(mean((exposure[present] - mu[present])^2))
  z <- (exposure[present] - mu[present]) / sigma
  mu_a <- mean(exposure)
  sd_a <- sqrt(mean((exposure - mu_a)^2))
  z_a <- (exposure - mu_a) / sd_a

  families <- list(
    list(
      input = dens_t(df = 4),
      g = function(values) stats::dt(values, df = 4)
    ),
    list(input = "kernel", g = continuous_density_kde)
  )

  for (family in families) {
    weights <- NULL
    expect_warning(
      weights <- wt_ate(
        trimmed,
        exposure,
        exposure_type = "continuous",
        stabilize = TRUE,
        .density = family$input
      ),
      class = "propensity_no_refit_warning"
    )
    weights <- as.numeric(weights)

    expect_identical(which(is.na(weights)), missing_at)

    f_den <- family$g(z) / sigma
    f_num <- family$g(z_a) / sd_a

    expect_equal(weights[present], f_num[present] / f_den, tolerance = 1e-12)
  }
})

test_that("no residual at all leaves every weight missing", {
  exposure <- continuous_density_data$exposure
  mu <- rep(NA_real_, continuous_density_data$n)

  # Every family answers alike, the kernel included. There is nothing for a
  # kernel to be fit on, but the density is not the reason these weights are
  # missing, so it is not asked and has nothing to object to.
  for (family in list("normal", dens_t(df = 4), "kernel")) {
    weights <- expect_silent(wt_ate(
      mu,
      exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = family
    ))

    expect_length(weights, continuous_density_data$n)
    expect_true(all(is.na(as.numeric(weights))))
  }
})

test_that("an infinite exposure or fitted value is refused where it arrives", {
  exposure <- continuous_density_data$exposure
  mu <- continuous_density_data$mu

  infinite_exposure <- exposure
  infinite_exposure[[3]] <- Inf
  infinite_mu <- mu
  infinite_mu[[5]] <- -Inf

  # A missing value is a unit with nothing to weight and leaves that unit's
  # weight missing. An infinite one is not that: the spread computed from it is
  # infinite, so every weight comes back missing however few units carry the
  # infinity, and under a kernel the residuals it leaves are refused as
  # residuals that do not vary, which they do. Neither report names the value
  # that caused it, so the value is refused where it arrives.
  for (family in list("normal", dens_t(df = 4), "kernel")) {
    expect_error(
      wt_ate(
        mu,
        infinite_exposure,
        exposure_type = "continuous",
        stabilize = TRUE,
        .density = family
      ),
      class = "propensity_density_error"
    )

    expect_error(
      wt_ate(
        infinite_mu,
        exposure,
        exposure_type = "continuous",
        stabilize = TRUE,
        .density = family
      ),
      class = "propensity_density_error"
    )
  }

  err <- expect_error(
    wt_ate(
      mu,
      infinite_exposure,
      exposure_type = "continuous",
      stabilize = TRUE
    ),
    class = "propensity_density_error"
  )
  message <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(message, "infinite", fixed = TRUE)
  expect_match(message, ".exposure", fixed = TRUE)

  expect_propensity_error(wt_ate(
    mu,
    infinite_exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  ))
  expect_propensity_error(wt_ate(
    infinite_mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  ))

  # A missing value still weights the way it always has: the units that have a
  # residual are weighted and the ones that do not are missing.
  missing_exposure <- exposure
  missing_exposure[[3]] <- NA_real_
  weights <- expect_silent(wt_ate(
    mu,
    missing_exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  ))
  expect_identical(which(is.na(as.numeric(weights))), 3L)
})

# ---- a density that reaches zero --------------------------------------------

# A density with no support beyond two standardized residuals: what a
# light-tailed family comes to at an outlier, and what a kernel estimate returns
# past the range it was fit on. The conditional density is what a continuous
# exposure's weights divide by, so a unit whose own exposure sits outside the
# support has a weight of infinity.
truncated_support_density <- function(limit = 2) {
  dens_fn(function(z) ifelse(abs(z) > limit, 0, dnorm(z)))
}

test_that("a conditional density of zero is refused rather than weighted as infinite", {
  exposure <- continuous_density_data$exposure
  mu <- continuous_density_data$mu

  zero_density_weights <- function(...) {
    wt_ate(
      mu,
      exposure,
      exposure_type = "continuous",
      .density = truncated_support_density(),
      ...
    )
  }

  # A propensity score of exactly 0 or 1 is refused on the binary path because
  # the weight it leaves is undefined. A conditional density of zero is the same
  # case for a continuous exposure, and is answered the same way whatever stands
  # over it: the marginal density, nothing at all, or the marginalized one.
  expect_error(
    zero_density_weights(stabilize = TRUE),
    class = "propensity_density_error"
  )
  expect_error(
    zero_density_weights(stabilize = FALSE),
    class = "propensity_density_error"
  )
  expect_error(
    zero_density_weights(stabilize = TRUE, numerator = "integrated"),
    class = "propensity_density_error"
  )

  # One unit of this fixture sits outside the support, and the refusal says
  # which one and what to do about it.
  err <- expect_error(
    zero_density_weights(stabilize = TRUE),
    class = "propensity_density_error"
  )
  message <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(message, "33", fixed = TRUE)
  expect_match(message, "heavier|dens_t|dens_kernel")
})

test_that("the zero-density refusal reads the way it is written", {
  expect_propensity_error(wt_ate(
    continuous_density_data$mu,
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = truncated_support_density()
  ))
})

test_that("censoring weights refuse a conditional density of zero as well", {
  expect_error(
    wt_cens(
      continuous_density_data$mu,
      continuous_density_data$exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = truncated_support_density()
    ),
    class = "propensity_density_error"
  )
})

test_that("a weight that is merely large is still a weight", {
  # A conditional density that is small everywhere rather than zero anywhere.
  # The weights it leaves are enormous, which is a fact about the fit and not an
  # arithmetic failure, so they are returned.
  weights <- expect_silent(wt_ate(
    continuous_density_data$mu,
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_fn(function(z) dnorm(z, sd = 0.2))
  ))

  values <- as.numeric(weights)
  expect_true(all(is.finite(values)))
  expect_gt(max(values), 1e6)
})

test_that("a zero in the numerator is a weight of zero, not a refusal", {
  # An exposure the model follows closely, with one unit far out in the
  # marginal distribution and no unit far out in the conditional one. The
  # marginal density is zero at that unit and the conditional density is
  # positive at every unit, so the ratio is a weight of zero rather than an
  # infinite one.
  x <- c(seq(-1, 1, length.out = 20), 8)
  mu <- 2 * x
  exposure <- mu + rep(c(-0.4, 0.4), length.out = length(x))

  weights <- expect_silent(wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = truncated_support_density()
  ))

  values <- as.numeric(weights)
  expect_true(all(is.finite(values)))
  expect_identical(which(values == 0), length(x))
  expect_true(all(values[-length(x)] > 0))
})

# ---- the finite-variance boundary -------------------------------------------

# Stabilized normal density-ratio weights have a finite second moment only while
# the marginal variance of the exposure stays below twice the variance of the
# conditional density. This builds an exposure sitting at a chosen point
# relative to that boundary: `share` is Var(mu) as a multiple of the conditional
# variance, so the marginal variance is `share + 1` times it and the boundary is
# at `share = 1`. The residuals are made orthogonal to the covariate so that the
# fitted model reproduces the decomposition exactly rather than to within
# sampling.
boundary_exposure <- function(share) {
  withr::local_seed(913)

  n <- 200
  x <- rnorm(n)
  x <- x - mean(x)
  e <- as.numeric(residuals(lm(rnorm(n) ~ x)))
  b <- sqrt(share * mean(e^2) / mean(x^2))
  exposure <- b * x + e

  list(n = n, x = x, exposure = exposure, mu = b * x)
}

# How many messages of one class a single call raises.
count_class_messages <- function(expr, class) {
  n <- 0L
  withCallingHandlers(
    expr,
    message = function(cnd) {
      if (inherits(cnd, class)) {
        n <<- n + 1L
        rlang::cnd_muffle(cnd)
      }
    }
  )

  n
}

test_that("stabilized normal weights report an exposure past the finite-variance boundary", {
  past <- boundary_exposure(1.4)

  wt_past <- function() {
    wt_ate(
      past$mu,
      past$exposure,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  }

  cnd <- expect_message(
    wt_past(),
    class = "propensity_density_variance_message"
  )
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "variance", fixed = TRUE)

  # One report for one call, whatever the weights are computed over.
  expect_identical(
    count_class_messages(wt_past(), "propensity_density_variance_message"),
    1L
  )

  # The weights are a diagnostic short of nothing: they are still returned, and
  # they are still finite.
  weights <- suppressMessages(wt_past())
  expect_length(weights, past$n)
  expect_true(all(is.finite(as.numeric(weights))))
})

test_that("the finite-variance report reads the way it is written", {
  past <- boundary_exposure(1.4)

  expect_propensity_warning(wt_ate(
    past$mu,
    past$exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  ))
})

test_that("the finite-variance boundary is read only where it holds", {
  below <- boundary_exposure(0.4)
  near <- boundary_exposure(0.9)
  past <- boundary_exposure(1.4)

  wt_continuous <- function(fixture, ...) {
    wt_ate(
      fixture$mu,
      fixture$exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      ...
    )
  }

  # Well inside the boundary, and just inside it. The second moment exists on
  # both sides of that distinction, and the report is a statement that it does
  # not, so it is made at the boundary rather than in front of it.
  expect_no_message(wt_continuous(below))
  expect_no_message(wt_continuous(near))

  # The boundary is a property of the normal family read against the exposure's
  # own marginal density. Unstabilized weights have no marginal density over
  # them, a heavier-tailed family has a different boundary, a marginalized
  # numerator is not the marginal density, and a score the caller supplied is
  # not a density at all.
  expect_no_message(wt_ate(
    past$mu,
    past$exposure,
    exposure_type = "continuous",
    stabilize = FALSE
  ))
  expect_no_message(wt_continuous(past, .density = dens_t(df = 4)))
  expect_no_message(wt_continuous(past, numerator = "integrated"))
  expect_no_message(wt_continuous(
    past,
    stabilization_score = rep(0.5, past$n)
  ))

  # One spread for each unit is not one conditional variance, so there is no
  # boundary to read the marginal variance against.
  expect_no_message(wt_continuous(
    past,
    .sigma = seq(0.3, 0.7, length.out = past$n)
  ))
})

# ---- against WeightIt -------------------------------------------------------

test_that("the denominator is the one WeightIt divides by", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")

  set.seed(14)
  n <- 400
  x <- rnorm(n)
  # A residual spread well away from one, so that a denominator missing its
  # factor of `1 / sigma` is a denominator on a visibly different scale.
  exposure <- 1.5 * x + rnorm(n, sd = 2)
  model_data <- data.frame(exposure = exposure, x = x)

  fit <- lm(exposure ~ x, data = model_data)
  residual <- as.numeric(residuals(fit))
  sigma <- sqrt(mean(residual^2))
  z <- residual / sigma

  families <- list(
    list(
      ours = "normal",
      theirs = "dnorm",
      g = stats::dnorm(z)
    ),
    list(
      ours = dens_t(df = 4),
      theirs = "dt_4",
      g = stats::dt(z, df = 4)
    ),
    list(
      ours = "laplace",
      theirs = "dlaplace",
      g = exp(-abs(z)) / 2
    )
  )

  for (family in families) {
    f_den <- family$g / sigma

    # WeightIt exposes no unstabilized weight, so the comparison of the
    # denominator itself is against the standardized residual computation
    # written out by hand.
    unstabilized <- wt_ate(
      as.numeric(fitted(fit)),
      exposure,
      exposure_type = "continuous",
      stabilize = FALSE,
      .density = family$ours
    )
    expect_equal(as.numeric(unstabilized), 1 / f_den, tolerance = 1e-12)

    # WeightIt's own weights carry its integrated numerator, so they are not
    # comparable to ours until `numerator = "integrated"` exists. What they do
    # pin now is that the denominator underneath them is this one: dividing
    # them out leaves a density in the exposure, which integrates to one over
    # the exposure's range less the mass in the tails beyond it. A denominator
    # on a different scale would leave something that does not.
    theirs <- WeightIt::weightit(
      exposure ~ x,
      data = model_data,
      method = "glm",
      density = family$theirs
    )
    implied <- as.numeric(theirs$weights) * f_den

    ordered <- order(exposure)
    height <- implied[ordered]
    width <- diff(exposure[ordered])
    mass <- sum(width * (utils::head(height, -1) + utils::tail(height, -1)) / 2)

    expect_equal(mass, 1, tolerance = 0.05)
  }
})

# ---- the integrated numerator -----------------------------------------------

# The number of points the integrated numerator marginalizes the conditional
# density over. It is an internal constant rather than an argument, so the hand
# calculations write it out once here.
continuous_numerator_grid <- 50L

# The integrated numerator, written out: average the conditional density over
# the units at each of 50 points spanning the exposure, then interpolate back to
# each observation with `stats::spline()`. `g` evaluates the density on a
# standardized residual.
continuous_integrated_numerator <- function(
  g,
  exposure = continuous_density_data$exposure,
  mu = continuous_density_data$mu,
  sigma
) {
  grid <- seq(
    min(exposure),
    max(exposure),
    length.out = continuous_numerator_grid
  )
  standardized <- outer(grid, mu, "-") / sigma
  on_grid <- rowMeans(matrix(
    g(as.vector(standardized)),
    nrow = continuous_numerator_grid
  )) /
    sigma

  stats::spline(grid, on_grid, xout = exposure, method = "fmm")$y
}

# The whole integrated ratio for a parametric family, denominator included.
continuous_integrated_wt <- function(
  g,
  exposure = continuous_density_data$exposure,
  mu = continuous_density_data$mu
) {
  sigma <- sqrt(mean((exposure - mu)^2))
  f_den <- g((exposure - mu) / sigma) / sigma
  f_num <- continuous_integrated_numerator(
    g,
    exposure = exposure,
    mu = mu,
    sigma = sigma
  )

  f_num / f_den
}

# The integrated ratio for a kernel, which is one estimate rather than two: the
# same fit answers for the standardized residuals and for the whole standardized
# grid, so it is fit over a range covering both.
continuous_integrated_kernel_wt <- function(
  exposure = continuous_density_data$exposure,
  mu = continuous_density_data$mu
) {
  sigma <- sqrt(mean((exposure - mu)^2))
  z <- (exposure - mu) / sigma

  grid <- seq(
    min(exposure),
    max(exposure),
    length.out = continuous_numerator_grid
  )
  standardized <- outer(grid, mu, "-") / sigma
  span <- range(c(z, as.vector(standardized)))

  fit <- stats::density(
    z,
    bw = "nrd0",
    adjust = 1,
    kernel = "gaussian",
    n = 512,
    from = span[1],
    to = span[2]
  )
  kde <- stats::approxfun(fit$x, fit$y)

  f_den <- kde(z) / sigma
  on_grid <- rowMeans(matrix(
    kde(as.vector(standardized)),
    nrow = continuous_numerator_grid
  )) /
    sigma
  f_num <- stats::spline(grid, on_grid, xout = exposure, method = "fmm")$y

  f_num / f_den
}

test_that("the integrated numerator is the grid marginalization by hand", {
  families <- list(
    list(input = "normal", g = function(z) stats::dnorm(z)),
    list(input = dens_t(df = 4), g = function(z) stats::dt(z, df = 4)),
    list(input = "laplace", g = function(z) exp(-abs(z)) / 2),
    list(
      input = function(z) stats::dlogis(z),
      g = function(z) stats::dlogis(z)
    )
  )

  for (family in families) {
    weights <- continuous_density_wt(
      .density = family$input,
      numerator = "integrated"
    )

    expect_equal(
      as.numeric(weights),
      continuous_integrated_wt(family$g),
      tolerance = 1e-12
    )
  }
})

test_that("the integrated numerator differs from the marginal one", {
  # The two coincide only when the fitted means are themselves normal, so a
  # family with heavier tails than the fitted means separates them. This is the
  # guard that the integrated route is doing its own arithmetic rather than
  # falling back on the marginal moments.
  marginal <- as.numeric(continuous_density_wt(.density = dens_t(df = 4)))
  integrated <- as.numeric(continuous_density_wt(
    .density = dens_t(df = 4),
    numerator = "integrated"
  ))

  expect_false(isTRUE(all.equal(marginal, integrated, tolerance = 1e-6)))
})

test_that("an integrated kernel is one estimate over the grid and the residuals", {
  weights <- continuous_density_wt(
    .density = "kernel",
    numerator = "integrated"
  )

  expect_equal(
    as.numeric(weights),
    continuous_integrated_kernel_wt(),
    tolerance = 1e-12
  )

  # The range the kernel is fit over has to cover the whole standardized grid,
  # which reaches further than the residuals do. A kernel fit only over the
  # residuals would interpolate to missing values on the grid, and the density
  # output check would refuse them, so the weights existing at all is part of
  # what this pins.
  expect_true(all(is.finite(as.numeric(weights))))
})

test_that("an intercept-only model gives integrated weights of exactly one", {
  exposure <- continuous_density_data$exposure
  mu <- rep(mean(exposure), continuous_density_data$n)

  # With nothing to condition on, the conditional density of the exposure is its
  # marginal density and the weights are one. Sending them through the grid and
  # the spline instead returns values near but not equal to one, off by about
  # 1e-05 on this problem, so the equality asked for here is the identity of the
  # numbers rather than a tolerance the grid could also meet.
  weights <- wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "normal",
    numerator = "integrated"
  )

  expect_identical(as.numeric(weights), rep(1, continuous_density_data$n))

  # Not only the normal family: any family read at the same value in the
  # numerator and the denominator cancels.
  heavy <- wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_t(df = 4),
    numerator = "integrated"
  )

  expect_identical(as.numeric(heavy), rep(1, continuous_density_data$n))

  # The fitted values an intercept-only model actually returns are not the same
  # double: the decomposition behind them leaves the last few bits of each one
  # to the arithmetic, so this problem has seven distinct values spanning
  # 4e-16. They are one number as far as the model is concerned, and the
  # weights they give are one as far as the package is concerned.
  fitted_mu <- as.numeric(fitted(lm(exposure ~ 1)))
  expect_gt(length(unique(fitted_mu)), 1)

  from_model <- wt_ate(
    fitted_mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "normal",
    numerator = "integrated"
  )

  expect_identical(as.numeric(from_model), rep(1, continuous_density_data$n))
})

test_that("integrated weights do not depend on the units of the exposure", {
  exposure <- continuous_density_data$exposure
  mu <- continuous_density_data$mu

  # A density ratio is a ratio of densities in the same units, so it is the
  # same number whether the dose is read in grams or in nanograms. Nothing
  # about the exposure changes here except the unit it is measured in.
  unscaled <- wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    numerator = "integrated"
  )

  scaled <- wt_ate(
    mu * 1e-9,
    exposure * 1e-9,
    exposure_type = "continuous",
    stabilize = TRUE,
    numerator = "integrated"
  )

  expect_equal(as.numeric(scaled), as.numeric(unscaled), tolerance = 1e-8)

  # The failure this guards against is silent: fitted means that vary by less
  # than the absolute floor of the constancy test are read as one number, and
  # every weight comes back as exactly one.
  expect_false(all(as.numeric(scaled) == 1))
})

test_that("integrated weights do not depend on where the exposure is centered", {
  exposure <- continuous_density_data$exposure
  mu <- continuous_density_data$mu
  offset <- 1e9

  # A density ratio is invariant to the origin the dose is measured from as
  # well as to the unit it is measured in: shifting the exposure and the fitted
  # means by the same constant leaves every standardized residual and every
  # standardized grid point where it was.
  unshifted <- wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    numerator = "integrated"
  )

  shifted <- wt_ate(
    mu + offset,
    exposure + offset,
    exposure_type = "continuous",
    stabilize = TRUE,
    numerator = "integrated"
  )

  # The comparison is to the precision doubles hold a shifted dose at, which is
  # about seven digits below the shift itself.
  expect_equal(
    as.numeric(shifted),
    as.numeric(unshifted),
    tolerance = 1e-5
  )

  # The failure this guards against is silent: a constancy test whose floor is
  # a fraction of the mean reads fitted means that differ by a unit or two at
  # an offset of a billion as one number, and returns a weight of one for every
  # unit of a model that conditions on a covariate.
  expect_false(all(as.numeric(shifted) == 1))
})

test_that("is_constant() reads a spread against the scale it is on", {
  mu <- continuous_density_data$mu
  sigma <- continuous_density_pooled()

  # The fitted means of a model with a covariate vary, and they vary just as
  # much relative to the spread of the exposure when both are rescaled.
  expect_false(is_constant(mu, scale = sigma))
  expect_false(is_constant(mu * 1e-9, scale = sigma * 1e-9))

  # Fitted means that really are one number are still read as one number.
  expect_true(is_constant(rep(mean(mu), length(mu)), scale = sigma))
  expect_true(is_constant(
    rep(mean(mu) * 1e-9, length(mu)),
    scale = sigma * 1e-9
  ))

  # An offset says nothing about how much the fitted means vary. Read against a
  # floor that is a fraction of the mean, means that differ by a unit or two at
  # an offset of a billion are one number, which they are not.
  expect_false(is_constant(mu + 1e9, scale = sigma))

  # What the offset term is for is still met: the fitted values of an
  # intercept-only model at that offset spread over the last bits of the
  # arithmetic that produced them, about 1e-5 here, and are one number as far
  # as the model is concerned.
  offset_intercept <- as.numeric(fitted(
    lm(I(continuous_density_data$exposure + 1e9) ~ 1)
  ))
  expect_gt(diff(range(offset_intercept)), 0)
  expect_true(is_constant(offset_intercept, scale = sigma))
})

test_that("an exposure that does not vary has no grid to be marginalized over", {
  # Every unit took the same dose, so the grid the conditional density is
  # averaged over collapses to a point and the interpolation has nowhere to run.
  # `stats::spline()` reaches that as a base warning about tied points; it is
  # refused here first, in terms of the exposure.
  n <- continuous_density_data$n
  exposure <- rep(2, n)
  mu <- seq(1.5, 2.5, length.out = n)

  expect_error(
    wt_ate(
      mu,
      exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      numerator = "integrated"
    ),
    class = "propensity_density_error"
  )

  constant_error <- rlang::catch_cnd(wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    numerator = "integrated"
  ))
  expect_match(conditionMessage(constant_error), "does not vary")
})

test_that("the integrated numerator reads only the units with a residual", {
  set.seed(24)
  n <- 50
  x <- rnorm(n)
  exposure <- 0.5 * x + rnorm(n)

  # A score a modification can be applied to has to lie in (0, 1); it stands in
  # for the fitted conditional mean here. Trimmed and not refit, it carries a
  # missing value at every unit it set aside.
  trimmed <- ps_trim(plogis(0.9 * x), method = "ps", lower = 0.2, upper = 0.8)
  mu <- as.numeric(trimmed)
  missing_at <- which(is.na(mu))
  expect_gt(length(missing_at), 1)

  present <- setdiff(seq_len(n), missing_at)

  # The ends of the grid are the smallest and largest exposure among the units
  # that have a residual, the average of the conditional density is an average
  # over those units, and the spread that standardizes it is theirs as well. So
  # the weights the units that are left carry are the weights of the problem
  # those units make on their own, and a unit with no residual has no weight.
  families <- list(
    list(
      input = dens_t(df = 4),
      oracle = function(exposure, mu) {
        continuous_integrated_wt(
          function(z) stats::dt(z, df = 4),
          exposure = exposure,
          mu = mu
        )
      }
    ),
    list(input = "kernel", oracle = continuous_integrated_kernel_wt)
  )

  for (family in families) {
    weights <- NULL
    expect_warning(
      weights <- wt_ate(
        trimmed,
        exposure,
        exposure_type = "continuous",
        stabilize = TRUE,
        .density = family$input,
        numerator = "integrated"
      ),
      class = "propensity_no_refit_warning"
    )
    weights <- as.numeric(weights)

    expect_identical(which(is.na(weights)), missing_at)
    expect_equal(
      weights[present],
      family$oracle(exposure[present], mu[present]),
      tolerance = 1e-12
    )
  }

  # When no unit has a residual there is nothing to read: no grid, no average,
  # and no weight anywhere. A kernel has nothing to be fit on either, and is not
  # asked, since the density is not the reason these weights are missing.
  none <- expect_silent(wt_ate(
    rep(NA_real_, n),
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "kernel",
    numerator = "integrated"
  ))

  expect_true(all(is.na(as.numeric(none))))
})

# ---- refusals of the numerator ----------------------------------------------

test_that("an integrated numerator refuses a binary exposure", {
  set.seed(21)
  ps <- runif(20, 0.2, 0.8)
  trt <- rbinom(20, 1, ps)

  expect_error(
    wt_ate(ps, trt, exposure_type = "binary", numerator = "integrated"),
    class = "propensity_numerator_error"
  )
  expect_error(
    wt_cens(ps, trt, exposure_type = "binary", numerator = "integrated"),
    class = "propensity_numerator_error"
  )

  # The type was not named here, so the refusal is against the type the
  # detection resolved the 0/1 exposure to.
  expect_error(
    wt_ate(ps, trt, numerator = "integrated"),
    class = "propensity_numerator_error"
  )

  binary_error <- rlang::catch_cnd(
    wt_ate(ps, trt, exposure_type = "binary", numerator = "integrated")
  )
  expect_match(conditionMessage(binary_error), "binary")
})

test_that("an integrated numerator refuses a categorical exposure", {
  exposure <- factor(rep(c("a", "b", "c"), each = 4))
  ps <- matrix(
    rep(c(0.5, 0.3, 0.2), times = 12),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )

  expect_error(
    wt_ate(
      ps,
      exposure,
      exposure_type = "categorical",
      numerator = "integrated"
    ),
    class = "propensity_numerator_error"
  )

  categorical_error <- rlang::catch_cnd(
    wt_ate(
      ps,
      exposure,
      exposure_type = "categorical",
      numerator = "integrated"
    )
  )
  expect_match(conditionMessage(categorical_error), "categorical")
})

test_that("the default numerator is accepted for every exposure type", {
  set.seed(22)
  ps <- runif(20, 0.2, 0.8)
  trt <- rbinom(20, 1, ps)

  exposure <- factor(rep(c("a", "b", "c"), each = 4))
  categorical_ps <- matrix(
    rep(c(0.5, 0.3, 0.2), times = 12),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )

  # Writing out the numerator the weights were getting anyway is not a reason to
  # refuse them, so the default holds for every type and is ignored outside the
  # continuous route.
  expect_equal(
    as.numeric(wt_ate(
      ps,
      trt,
      exposure_type = "binary",
      numerator = "marginal"
    )),
    as.numeric(wt_ate(ps, trt, exposure_type = "binary"))
  )
  expect_equal(
    as.numeric(wt_ate(
      categorical_ps,
      exposure,
      exposure_type = "categorical",
      numerator = "marginal"
    )),
    as.numeric(wt_ate(categorical_ps, exposure, exposure_type = "categorical"))
  )
  expect_null(density_meta(
    wt_ate(ps, trt, exposure_type = "binary", numerator = "marginal")
  ))
})

test_that("an integrated numerator refuses weights that are not stabilized", {
  # The integrated numerator is a numerator, so there is nothing for it to be
  # when the weights carry none.
  expect_error(
    continuous_density_wt(numerator = "integrated", stabilize = FALSE),
    class = "propensity_numerator_error"
  )
  expect_error(
    continuous_density_wt(
      .density = dens_t(df = 4),
      numerator = "integrated",
      stabilize = FALSE
    ),
    class = "propensity_numerator_error"
  )

  stabilize_error <- rlang::catch_cnd(
    continuous_density_wt(numerator = "integrated", stabilize = FALSE)
  )
  expect_match(conditionMessage(stabilize_error), "stabilize")
})

test_that("an integrated numerator refuses a stabilization score", {
  score <- seq(0.4, 1.6, length.out = continuous_density_data$n)

  # A score the caller supplies is the numerator, so it cannot also be
  # marginalized out of the conditional density.
  expect_error(
    continuous_density_wt(
      numerator = "integrated",
      stabilization_score = score
    ),
    class = "propensity_numerator_error"
  )

  score_error <- rlang::catch_cnd(
    continuous_density_wt(
      numerator = "integrated",
      stabilization_score = score
    )
  )
  expect_match(conditionMessage(score_error), "stabilization_score")
})

test_that("an integrated numerator refuses a supplied spread", {
  # The grid marginalization reads the conditional density at every unit's
  # fitted mean, which a spread the caller supplied rather than the model
  # estimated has no standing to describe.
  expect_error(
    continuous_density_wt(numerator = "integrated", .sigma = 0.85),
    class = "propensity_numerator_error"
  )
  expect_error(
    continuous_density_wt(
      numerator = "integrated",
      .sigma = continuous_density_data$sigma_i
    ),
    class = "propensity_numerator_error"
  )

  sigma_error <- rlang::catch_cnd(
    continuous_density_wt(numerator = "integrated", .sigma = 0.85)
  )
  expect_match(conditionMessage(sigma_error), "sigma")
})

test_that("a numerator that is not a numerator is refused", {
  # A mistyped value is corrected by the argument matcher rather than refused on
  # the package's own terms.
  expect_error(
    continuous_density_wt(numerator = "intergrated"),
    regexp = "must be one of"
  )
  expect_error(
    continuous_density_wt(numerator = "none"),
    regexp = "must be one of"
  )
})

test_that("the refusals of an integrated numerator read the way they should", {
  set.seed(25)
  ps <- runif(20, 0.2, 0.8)
  trt <- rbinom(20, 1, ps)

  categorical_exposure <- factor(rep(c("a", "b", "c"), each = 4))
  categorical_ps <- matrix(
    rep(c(0.5, 0.3, 0.2), times = 12),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )

  expect_propensity_error(
    wt_ate(ps, trt, exposure_type = "binary", numerator = "integrated")
  )
  expect_propensity_error(
    wt_ate(
      categorical_ps,
      categorical_exposure,
      exposure_type = "categorical",
      numerator = "integrated"
    )
  )
  expect_propensity_error(
    continuous_density_wt(numerator = "integrated", stabilize = FALSE)
  )
  expect_propensity_error(
    continuous_density_wt(
      numerator = "integrated",
      stabilization_score = seq(
        0.4,
        1.6,
        length.out = continuous_density_data$n
      )
    )
  )
  expect_propensity_error(
    continuous_density_wt(numerator = "integrated", .sigma = 0.85)
  )

  # And the two refusals that are about the problem rather than the arguments:
  # an exposure with no spread for the grid to run over, and an interpolated
  # numerator that dipped below zero.
  expect_propensity_error(
    wt_ate(
      seq(1.5, 2.5, length.out = continuous_density_data$n),
      rep(2, continuous_density_data$n),
      exposure_type = "continuous",
      stabilize = TRUE,
      numerator = "integrated"
    )
  )
  expect_propensity_error(
    continuous_density_wt(
      .density = function(r) stats::dnorm(r, sd = 0.1),
      numerator = "integrated"
    )
  )
})

# ---- positivity of the interpolated numerator -------------------------------

test_that("an interpolated numerator that dips below zero is refused", {
  # A density this narrow leaves the marginalized density on the grid close to
  # zero between sharp peaks, and a cubic spline through those points undershoots
  # below zero at three of the observed exposures. The refusal is the ordinary
  # output check on a density, reached through the non-finite-or-negative rule:
  # the denominator is positive everywhere here, so nothing else is wrong with
  # these weights.
  narrow <- function(r) stats::dnorm(r, sd = 0.1)

  expect_error(
    continuous_density_wt(.density = narrow, numerator = "integrated"),
    class = "propensity_density_error"
  )

  # The message names the numerator, so a reader can tell an interpolated
  # marginal density that went negative from a conditional one that did.
  numerator_error <- rlang::catch_cnd(
    continuous_density_wt(.density = narrow, numerator = "integrated")
  )
  expect_match(conditionMessage(numerator_error), "numerator")

  # The same density with the marginal numerator is fine, so the refusal is
  # about the interpolation rather than about the family.
  expect_silent(continuous_density_wt(.density = narrow))
})

# ---- what the integrated weights record -------------------------------------

test_that("integrated weights record the numerator they were built from", {
  for (family in list("normal", dens_t(df = 4), "laplace", "kernel")) {
    record <- density_meta(continuous_density_wt(
      .density = family,
      numerator = "integrated"
    ))

    expect_s3_class(record, "propensity_density_meta")
    expect_identical(record$numerator, "integrated")
    expect_identical(record$sigma, "pooled")
  }

  # And the marginal numerator still records itself, written out or not.
  expect_identical(
    density_meta(continuous_density_wt(numerator = "marginal"))$numerator,
    "marginal"
  )
})

test_that("weights that disagree on the numerator drop the density record", {
  # The two records agree on the family and the spread and differ only on the
  # numerator, and the pair is compatible on every field that decides the type:
  # both are stabilized and neither carries a score. So the combine is a
  # disagreement about the density record alone.
  marginal <- continuous_density_wt(.density = dens_t(df = 4))
  integrated <- continuous_density_wt(
    .density = dens_t(df = 4),
    numerator = "integrated"
  )

  combined <- NULL
  expect_warning(
    combined <- c(marginal, integrated),
    class = "propensity_metadata_conflict_warning"
  )

  expect_s3_class(combined, "psw")
  expect_length(combined, 2 * continuous_density_data$n)
  expect_null(density_meta(combined))

  # The exposure type is not what they disagreed on, so it survives.
  expect_identical(exposure_type(combined), "continuous")
})

test_that("weights that agree on the integrated numerator combine", {
  first <- continuous_density_wt(
    .density = "laplace",
    numerator = "integrated"
  )
  second <- continuous_density_wt(
    .density = "laplace",
    numerator = "integrated"
  )

  combined <- expect_silent(c(first, second))

  expect_length(combined, 2 * continuous_density_data$n)
  expect_identical(density_meta(combined)$numerator, "integrated")
  expect_identical(format(density_meta(combined)$density), "laplace")
})

# ---- the numerator through the other methods --------------------------------

test_that("the numerator reaches the weights through a fitted model", {
  exposure <- continuous_density_data$exposure
  x <- continuous_density_data$x
  model_data <- data.frame(exposure = exposure, x = x)
  fit <- glm(exposure ~ x, data = model_data, family = gaussian())

  from_model <- wt_ate(
    fit,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_t(df = 4),
    numerator = "integrated"
  )

  expect_equal(
    as.numeric(from_model),
    continuous_integrated_wt(function(z) stats::dt(z, df = 4)),
    tolerance = 1e-12
  )
  expect_identical(density_meta(from_model)$numerator, "integrated")
})

test_that("the numerator reaches the weights through a data frame", {
  from_frame <- wt_ate(
    data.frame(mu = continuous_density_data$mu),
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "laplace",
    numerator = "integrated"
  )

  expect_equal(
    as.numeric(from_frame),
    continuous_integrated_wt(function(z) exp(-abs(z)) / 2),
    tolerance = 1e-12
  )
  expect_identical(density_meta(from_frame)$numerator, "integrated")
})

test_that("censoring weights take an integrated numerator of their own", {
  weights <- wt_cens(
    continuous_density_data$mu,
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_t(df = 4),
    numerator = "integrated"
  )

  expect_equal(
    as.numeric(weights),
    continuous_integrated_wt(function(z) stats::dt(z, df = 4)),
    tolerance = 1e-12
  )
  expect_identical(estimand(weights), "uncensored")
  expect_identical(density_meta(weights)$numerator, "integrated")
})

# The three modified propensity scores the numerator has to reach through,
# built on one seeded problem. A score a modification can be applied to has to
# lie in (0, 1); it stands in for the fitted conditional mean here.
#
# The seed is chosen so that the scores all lie inside the trimming bounds and
# nothing is trimmed away: with no unit set aside, none of the three leaves a
# missing value, and the grid marginalization reads every unit. That is a
# property of this fixture rather than of trimming, and it keeps these tests
# about the numerator reaching each route. What the marginalization does when a
# unit has no residual is a separate question, and the trimmed weights already
# have tests of their own above.
continuous_modified_scores <- function() {
  set.seed(23)
  n <- 40
  x <- rnorm(n)
  exposure <- 0.5 * x + rnorm(n)
  fit <- glm(
    exposure ~ x,
    data = data.frame(exposure = exposure, x = x),
    family = gaussian()
  )

  scores <- plogis(0.5 * x)
  trimmed <- ps_trim(scores, method = "ps", lower = 0.2, upper = 0.8)
  stopifnot(!anyNA(as.numeric(trimmed)))

  list(
    n = n,
    x = x,
    exposure = exposure,
    scores = list(
      refit = ps_refit(trimmed, model = fit),
      truncated = ps_trunc(scores, method = "ps", lower = 0.2, upper = 0.8),
      calibrated = ps_calibrate(scores, rbinom(n, 1, scores))
    )
  )
}

test_that("the numerator reaches the weights through a modified score", {
  problem <- continuous_modified_scores()

  for (modified in problem$scores) {
    weights <- wt_ate(
      modified,
      problem$exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = dens_t(df = 4),
      numerator = "integrated"
    )

    expect_equal(
      as.numeric(weights),
      continuous_integrated_wt(
        function(z) stats::dt(z, df = 4),
        exposure = problem$exposure,
        mu = as.numeric(modified)
      ),
      tolerance = 1e-12
    )
    expect_identical(density_meta(weights)$numerator, "integrated")
  }
})

test_that("censoring weights carry a numerator through their other methods", {
  from_frame <- wt_cens(
    data.frame(mu = continuous_density_data$mu),
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = "laplace",
    numerator = "integrated"
  )

  expect_equal(
    as.numeric(from_frame),
    continuous_integrated_wt(function(z) exp(-abs(z)) / 2),
    tolerance = 1e-12
  )
  expect_identical(estimand(from_frame), "uncensored")
  expect_identical(density_meta(from_frame)$numerator, "integrated")

  model_data <- data.frame(
    exposure = continuous_density_data$exposure,
    x = continuous_density_data$x
  )
  fit <- glm(exposure ~ x, data = model_data, family = gaussian())

  from_model <- wt_cens(
    fit,
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_t(df = 4),
    numerator = "integrated"
  )

  expect_equal(
    as.numeric(from_model),
    continuous_integrated_wt(function(z) stats::dt(z, df = 4)),
    tolerance = 1e-12
  )
  expect_identical(estimand(from_model), "uncensored")
  expect_identical(density_meta(from_model)$numerator, "integrated")
})

test_that("censoring weights carry a numerator through a modified score", {
  problem <- continuous_modified_scores()

  # A modification is recorded on the estimand after the censoring weights have
  # named themselves, so each of the three reads as the estimand and then the
  # modification.
  estimands <- c(
    refit = "uncensored; trimmed",
    truncated = "uncensored; truncated",
    calibrated = "uncensored; calibrated"
  )

  for (modification in names(problem$scores)) {
    modified <- problem$scores[[modification]]

    weights <- wt_cens(
      modified,
      problem$exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = dens_t(df = 4),
      numerator = "integrated"
    )

    expect_equal(
      as.numeric(weights),
      continuous_integrated_wt(
        function(z) stats::dt(z, df = 4),
        exposure = problem$exposure,
        mu = as.numeric(modified)
      ),
      tolerance = 1e-12
    )
    expect_identical(estimand(weights), estimands[[modification]])
    expect_identical(density_meta(weights)$numerator, "integrated")
  }
})

# ---- against WeightIt -------------------------------------------------------

test_that("integrated weights are the weights WeightIt gives", {
  skip_on_cran()
  skip_if_not_installed("WeightIt", "2.0.0")

  set.seed(123)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  exposure <- 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  model_data <- data.frame(exposure = exposure, x1 = x1, x2 = x2)

  ps_mod <- glm(exposure ~ x1 + x2, data = model_data, family = gaussian())

  # WeightIt divides neither density by the residual spread, and we divide both,
  # so the factor cancels and the two ratios are the same number rather than
  # proportional to each other. Both are checked: the equality is what the
  # implementation owes, and the proportionality says the shape agrees even if
  # the scale ever drifts.
  families <- list(
    list(ours = "normal", theirs = NULL),
    list(ours = dens_t(df = 4), theirs = "dt_4"),
    list(ours = "laplace", theirs = "dlaplace"),
    list(ours = "kernel", theirs = "kernel")
  )

  for (family in families) {
    ours <- wt_ate(
      ps_mod,
      exposure,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = family$ours,
      numerator = "integrated"
    )

    theirs <- if (is.null(family$theirs)) {
      WeightIt::weightit(
        exposure ~ x1 + x2,
        data = model_data,
        method = "glm"
      )$weights
    } else {
      WeightIt::weightit(
        exposure ~ x1 + x2,
        data = model_data,
        method = "glm",
        density = family$theirs
      )$weights
    }
    theirs <- as.numeric(theirs)

    expect_equal(
      as.numeric(ours),
      theirs,
      tolerance = 1e-6,
      ignore_attr = TRUE
    )

    expect_equal(
      as.numeric(ours) / mean(as.numeric(ours)),
      theirs / mean(theirs),
      tolerance = 1e-6,
      ignore_attr = TRUE
    )
  }
})

# ---- the t scale by maximum likelihood --------------------------------------

# A continuous-exposure problem whose residuals come from a law with heavy
# tails, which is what the t family is for. The root mean square of these
# residuals is pulled out by the few large ones; the maximum likelihood scale of
# the same t is not.
continuous_t_data <- local({
  set.seed(20260830)

  n <- 400
  x <- rnorm(n)
  exposure <- 1 + 0.7 * x + rt(n, df = 3)

  list(
    x = x,
    exposure = exposure,
    mu = as.numeric(fitted(lm(exposure ~ x)))
  )
})

continuous_t_residuals <- continuous_t_data$exposure - continuous_t_data$mu

# Unstabilized weights for that problem, which are the conditional density and
# nothing else. The marginal density that stabilizes weights is read at the
# exposure's own moments, and what is under test here is the spread of the
# conditional one.
continuous_t_wt <- function(..., exposure = continuous_t_data$exposure) {
  wt_ate(
    continuous_t_data$mu,
    exposure,
    exposure_type = "continuous",
    stabilize = FALSE,
    ...
  )
}

test_that("the t scale by maximum likelihood is the root of the t score", {
  # The scale the weights were built at is not reported on its own, so it is
  # read back through them: the same family spread by a number supplied by hand
  # is the same ratio, and the two agree only if the package solved the equation
  # the oracle solved.
  scale <- t_scale_mle(continuous_t_residuals, df = 3)

  mle <- continuous_t_wt(.density = dens_t(3, sigma_method = "mle"))
  oracle <- continuous_t_wt(.density = dens_t(3), .sigma = scale)

  expect_equal(as.numeric(mle), as.numeric(oracle), tolerance = 1e-8)
})

test_that("a maximum likelihood scale parts from the root mean square in heavy tails", {
  rms <- continuous_t_wt(.density = dens_t(3))
  mle <- continuous_t_wt(.density = dens_t(3, sigma_method = "mle"))

  relative <- abs(as.numeric(mle) / as.numeric(rms) - 1)

  # Neither set of weights is a rescaling of the other: the residuals that pull
  # the root mean square out are the ones the two spreads disagree about most,
  # and they are the units whose weights the choice of spread decides.
  expect_gt(max(relative), 0.5)
  expect_gt(mean(relative), 0.1)
})

test_that("the two spreads meet as the t approaches the normal", {
  # A t with many degrees of freedom is the normal, and the maximum likelihood
  # scale of a normal is its root mean square, so on residuals that are normal
  # the two estimators answer the same question and answer it alike.
  residuals <- continuous_density_data$exposure - continuous_density_data$mu
  scale <- t_scale_mle(residuals, df = 1000)

  expect_equal(scale, continuous_density_pooled(), tolerance = 0.005)

  rms <- continuous_density_wt(stabilize = FALSE, .density = dens_t(1000))
  mle <- continuous_density_wt(
    stabilize = FALSE,
    .density = dens_t(1000, sigma_method = "mle")
  )

  expect_lt(max(abs(as.numeric(mle) / as.numeric(rms) - 1)), 0.01)
})

test_that("the t family is spread by the root mean square unless asked otherwise", {
  # The default is the estimator every family takes, so weights written the way
  # they have always been written are the weights they have always been.
  default <- continuous_t_wt(.density = dens_t(6))
  asked <- continuous_t_wt(.density = dens_t(6, sigma_method = "rms"))

  expect_identical(as.numeric(default), as.numeric(asked))
  expect_identical(density_meta(default)$sigma, "pooled")
  expect_identical(density_meta(asked)$sigma, "pooled")
})

test_that("weights record a maximum likelihood spread as a source of its own", {
  mle <- continuous_t_wt(.density = dens_t(6, sigma_method = "mle"))

  # A third source beside the pooled root mean square and a spread the caller
  # supplied. It is estimated from the residuals, as the pooled one is, but by a
  # different equation, and that is the equation `ipw()` has to solve to rebuild
  # these weights.
  expect_identical(density_meta(mle)$sigma, "mle")

  # It is a function of the data rather than a number the caller chose, so there
  # is no constant to keep, exactly as for the pooled spread.
  expect_null(density_meta(mle)$sigma_value)
})

test_that("a maximum likelihood scale spreads the conditional density alone", {
  # `.sigma` and `sigma_method` describe the same thing: the spread of the
  # conditional density the weights divide by. The marginal density that
  # stabilizes them is the exposure's own, read at the exposure's own mean and
  # root mean square, which neither of them has ever changed.
  scale <- t_scale_mle(continuous_t_residuals, df = 6)

  stabilized <- wt_ate(
    continuous_t_data$mu,
    continuous_t_data$exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_t(6, sigma_method = "mle")
  )
  oracle <- wt_ate(
    continuous_t_data$mu,
    continuous_t_data$exposure,
    .sigma = scale,
    exposure_type = "continuous",
    stabilize = TRUE,
    .density = dens_t(6)
  )

  expect_equal(as.numeric(stabilized), as.numeric(oracle), tolerance = 1e-8)
})

test_that("a maximum likelihood scale is fit on the residuals that are there", {
  exposure <- continuous_t_data$exposure
  exposure[c(4, 19)] <- NA

  scale <- t_scale_mle(exposure - continuous_t_data$mu, df = 4)

  mle <- continuous_t_wt(
    exposure = exposure,
    .density = dens_t(4, sigma_method = "mle")
  )
  oracle <- continuous_t_wt(
    exposure = exposure,
    .density = dens_t(4),
    .sigma = scale
  )

  # A unit with no exposure has no residual for the likelihood to read and no
  # weight to be given, and it is the only unit that comes back missing.
  expect_equal(as.numeric(mle), as.numeric(oracle), tolerance = 1e-8)
  expect_true(all(is.na(as.numeric(mle)[c(4, 19)])))
  expect_false(anyNA(as.numeric(mle)[-c(4, 19)]))
})

test_that("a spread supplied and a spread estimated under the density are refused together", {
  # `.sigma` says the spread is a number of the caller's own, and
  # `sigma_method = "mle"` says it is estimated from the residuals: two
  # instructions about the same quantity, one of which would go unread whichever
  # way the pairing was resolved.
  expect_propensity_error(
    continuous_t_wt(.density = dens_t(4, sigma_method = "mle"), .sigma = 0.9)
  )

  # The default estimator takes a supplied spread as every other family does, so
  # what is refused is the pairing rather than the family.
  expect_no_error(continuous_t_wt(.density = dens_t(4), .sigma = 0.9))
})

# A propensity score model that all but interpolates its exposure: every unit
# but one is reproduced to within a billionth, and the one that is not is far
# away. Residuals like these have a maximum likelihood scale of their own, and
# it sits orders of magnitude below their root mean square, which the single
# large residual is almost all of.
continuous_interpolating_data <- local({
  set.seed(20260831)

  mu <- rnorm(60)

  list(mu = mu, exposure = mu + c(rep(1e-9, 59), 4))
})

test_that("a maximum likelihood scale far below the root mean square is found", {
  mu <- continuous_interpolating_data$mu
  exposure <- continuous_interpolating_data$exposure
  scale <- t_scale_mle(exposure - mu, df = 4)

  # The scale the likelihood is maximized at is nowhere near the spread the
  # residuals average to, so an estimator that only looks near the root mean
  # square never reaches it.
  expect_lt(scale, sqrt(mean((exposure - mu)^2)) / 1e4)

  mle <- wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = FALSE,
    .density = dens_t(4, sigma_method = "mle")
  )
  oracle <- wt_ate(
    mu,
    exposure,
    exposure_type = "continuous",
    stabilize = FALSE,
    .density = dens_t(4),
    .sigma = scale
  )

  expect_equal(as.numeric(mle), as.numeric(oracle), tolerance = 1e-8)
})

test_that("a likelihood with no maximum at a positive scale is refused", {
  # Residuals a model reproduced exactly say nothing about the spread of the
  # density around it, and once enough of them are exactly zero the likelihood
  # rises the whole way down: there is no scale to return, and how many is
  # enough is a question the degrees of freedom answer.
  mu <- seq(0, 1, length.out = 10)
  exposure <- mu
  exposure[[3]] <- exposure[[3]] + 1

  expect_propensity_error(
    wt_ate(
      mu,
      exposure,
      exposure_type = "continuous",
      stabilize = FALSE,
      .density = dens_t(4, sigma_method = "mle")
    )
  )

  # The same residuals under a t with lighter tails still leave a maximum: how
  # many exact zeros are too many is a question of how readily the family
  # dismisses the residual that is there as an outlier, and a nearly normal t
  # does not dismiss it at all.
  expect_no_error(
    wt_ate(
      mu,
      exposure,
      exposure_type = "continuous",
      stabilize = FALSE,
      .density = dens_t(20, sigma_method = "mle")
    )
  )
})

# ---- a numerator model supplied to stabilize --------------------------------

# A conditional numerator for the same problem: the exposure fit on a baseline
# covariate the propensity score model does not read. The weights it stabilizes
# are f(A | V) / f(A | X), each density read at the residual spread of its own
# model rather than at one spread shared between them.
continuous_numerator_model <- local({
  v <- withr::with_seed(
    20250830,
    stats::rnorm(continuous_density_data$n)
  )
  exposure <- continuous_density_data$exposure

  stats::lm(exposure ~ v)
})

# The pooled residual spread of a fitted model, which is what its density is
# read at: the uncentered root mean square, the same estimator the conditional
# density uses.
continuous_model_rms <- function(model) {
  sqrt(mean(stats::residuals(model)^2))
}

test_that("a numerator model supplied to stabilize is the conditional density ratio", {
  exposure <- continuous_density_data$exposure
  fit <- continuous_numerator_model

  f_num <- stats::dnorm(
    exposure,
    as.numeric(stats::fitted(fit)),
    continuous_model_rms(fit)
  )
  f_den <- stats::dnorm(
    exposure,
    continuous_density_data$mu,
    continuous_density_pooled()
  )

  weights <- wt_ate(
    continuous_density_data$mu,
    exposure,
    exposure_type = "continuous",
    stabilize = fit
  )

  expect_s3_class(weights, "psw")
  expect_equal(as.numeric(weights), f_num / f_den, tolerance = 1e-12)

  # The conditional numerator is not the marginal one, so the weights it builds
  # are a different set of weights rather than the same ones under another name.
  marginal <- continuous_density_wt()
  expect_false(isTRUE(all.equal(
    as.numeric(weights),
    as.numeric(marginal),
    tolerance = 1e-8
  )))
})

test_that("the weights record the numerator model they were stabilized on", {
  weights <- wt_ate(
    continuous_density_data$mu,
    continuous_density_data$exposure,
    exposure_type = "continuous",
    stabilize = continuous_numerator_model
  )
  record <- density_meta(weights)

  expect_true(is_stabilized(weights))
  expect_identical(estimand(weights), "ate")
  expect_identical(exposure_type(weights), "continuous")

  # The model itself is recorded, not the numerator it evaluates to: `ipw()`
  # rebuilds the numerator at every value of theta, which takes the design and
  # the coefficients rather than the fitted values.
  expect_s3_class(record, "propensity_density_meta")
  expect_identical(record$numerator, "model")
  expect_identical(record$numerator_model, continuous_numerator_model)
  expect_identical(record$sigma, "pooled")
  expect_identical(format(record$density), "normal")

  # A model is a numerator rather than a score, so nothing is recorded as one.
  expect_null(stabilization_score(weights))

  # The record prints the numerator it holds the way it prints every other
  # field of the ratio.
  expect_true(any(grepl("model", format(record), fixed = TRUE)))
})

test_that("the numerator model leaves the other stabilizations alone", {
  exposure <- continuous_density_data$exposure
  sigma <- continuous_density_pooled()

  f_den <- stats::dnorm(exposure, continuous_density_data$mu, sigma)
  f_marginal <- stats::dnorm(
    exposure,
    mean(exposure),
    continuous_density_sd_a()
  )

  marginal <- continuous_density_wt()
  expect_equal(as.numeric(marginal), f_marginal / f_den, tolerance = 1e-12)
  expect_identical(density_meta(marginal)$numerator, "marginal")
  expect_null(density_meta(marginal)$numerator_model)

  unstabilized <- continuous_density_wt(stabilize = FALSE)
  expect_equal(as.numeric(unstabilized), 1 / f_den, tolerance = 1e-12)
  expect_identical(density_meta(unstabilized)$numerator, "none")
  expect_null(density_meta(unstabilized)$numerator_model)
  expect_false(is_stabilized(unstabilized))
})

test_that("a numerator model of a dose refuses a binary or a categorical exposure", {
  set.seed(4218)
  ps <- runif(20, 0.2, 0.8)
  trt <- rbinom(20, 1, ps)
  dose <- rnorm(20)
  fit <- lm(dose ~ ps)

  # Stabilizing on a fitted numerator is a statement about the exposure the
  # model reads. A binary exposure takes one too, but its numerator is a
  # probability rather than a density, so a least-squares fit of a dose names
  # nothing a binary exposure's weights could divide by and is refused as the
  # wrong family; see tests/testthat/test-weights-numerator-binary.R for the
  # binomial fit that is the right one.
  expect_error(
    wt_ate(ps, trt, exposure_type = "binary", stabilize = fit),
    class = "propensity_model_family_error"
  )
  expect_error(
    wt_ate(ps, trt, stabilize = fit),
    class = "propensity_model_family_error"
  )

  # A categorical exposure's weights are a ratio over more than two levels and
  # take no fitted numerator at all, which is the refusal an integrated
  # numerator gets for the same reason.
  exposure <- factor(rep(c("a", "b", "c"), each = 4))
  categorical_ps <- matrix(
    rep(c(0.5, 0.3, 0.2), times = 12),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )
  categorical_fit <- lm(rnorm(12) ~ seq_len(12))

  expect_error(
    wt_ate(
      categorical_ps,
      exposure,
      exposure_type = "categorical",
      stabilize = categorical_fit
    ),
    class = "propensity_numerator_error"
  )
})

test_that("a numerator model is read the way a propensity score model is", {
  exposure <- continuous_density_data$exposure

  # A fit whose spread changes with its fitted values describes no single
  # conditional density, and it is refused as such wherever it arrives: the
  # numerator model is held to the family the propensity score model is held to.
  binary_outcome <- as.numeric(exposure > 0)
  wrong_family <- glm(
    binary_outcome ~ continuous_density_data$x,
    family = binomial()
  )

  expect_error(
    wt_ate(
      continuous_density_data$mu,
      exposure,
      exposure_type = "continuous",
      stabilize = wrong_family
    ),
    class = "propensity_model_family_error"
  )

  # A model fit to a different number of observations has one numerator for
  # each of its own units and none for the units here, so it is refused rather
  # than recycled against them.
  short <- stats::lm(
    exposure[1:20] ~ continuous_density_data$x[1:20]
  )

  expect_error(
    wt_ate(
      continuous_density_data$mu,
      exposure,
      exposure_type = "continuous",
      stabilize = short
    ),
    class = "propensity_length_error"
  )
})

test_that("a numerator model and a stabilization score are two numerators", {
  # A score the caller computed is itself the numerator of the weights, and a
  # model is a second one; the two together are two instructions about the same
  # quantity, the way a supplied spread and an estimated one are.
  expect_error(
    wt_ate(
      continuous_density_data$mu,
      continuous_density_data$exposure,
      exposure_type = "continuous",
      stabilize = continuous_numerator_model,
      stabilization_score = rep(0.5, continuous_density_data$n)
    ),
    class = "propensity_numerator_error"
  )

  # An integrated numerator marginalizes the conditional density over the
  # units, which is a third numerator and not one a model can be read into.
  expect_error(
    wt_ate(
      continuous_density_data$mu,
      continuous_density_data$exposure,
      exposure_type = "continuous",
      stabilize = continuous_numerator_model,
      numerator = "integrated"
    ),
    class = "propensity_numerator_error"
  )
})

test_that("the numerator model refusals read as the refusals they are", {
  exposure <- continuous_density_data$exposure
  fit <- continuous_numerator_model

  # A model of a dose handed to a binary exposure is refused as the wrong
  # family rather than as a model where none is taken, so its wording belongs
  # with the binary route's own refusals.

  expect_propensity_error(
    wt_ate(
      continuous_density_data$mu,
      exposure,
      exposure_type = "continuous",
      stabilize = fit,
      stabilization_score = rep(0.5, continuous_density_data$n)
    )
  )

  expect_propensity_error(
    wt_ate(
      continuous_density_data$mu,
      exposure,
      exposure_type = "continuous",
      stabilize = fit,
      numerator = "integrated"
    )
  )

  expect_propensity_error(
    wt_ate(
      continuous_density_data$mu,
      exposure,
      exposure_type = "continuous",
      stabilize = stats::glm(
        as.numeric(exposure > 0) ~ continuous_density_data$x,
        family = binomial()
      )
    )
  )

  expect_propensity_error(
    wt_ate(
      continuous_density_data$mu,
      exposure,
      exposure_type = "continuous",
      stabilize = stats::lm(
        exposure[1:20] ~ continuous_density_data$x[1:20]
      )
    )
  )
})

test_that("stabilize takes one of the answers or a model and nothing else", {
  # Every consumer of the resolved argument asks `isTRUE()` of it, so a value
  # naming neither answer would otherwise be read as `FALSE` and the weights
  # would come back unstabilized without a word.
  expect_propensity_error(
    wt_ate(
      continuous_density_data$mu,
      continuous_density_data$exposure,
      exposure_type = "continuous",
      stabilize = "yes"
    )
  )

  expect_propensity_error(
    wt_ate(
      continuous_density_data$mu,
      continuous_density_data$exposure,
      exposure_type = "continuous",
      stabilize = NA
    )
  )
})

test_that("censoring weights reach the numerator model route too", {
  # `wt_cens()` builds its weights with the ATE formula, so a numerator model
  # arrives there through the same argument and is recorded the same way.
  exposure <- continuous_density_data$exposure
  fit <- continuous_numerator_model

  weights <- wt_cens(
    continuous_density_data$mu,
    exposure,
    exposure_type = "continuous",
    stabilize = fit
  )

  f_num <- stats::dnorm(
    exposure,
    as.numeric(stats::fitted(fit)),
    continuous_model_rms(fit)
  )
  f_den <- stats::dnorm(
    exposure,
    continuous_density_data$mu,
    continuous_density_pooled()
  )

  expect_identical(estimand(weights), "uncensored")
  expect_equal(as.numeric(weights), f_num / f_den, tolerance = 1e-12)
  expect_identical(density_meta(weights)$numerator, "model")
  expect_identical(density_meta(weights)$numerator_model, fit)
})
