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
      weights <- wt_ate(
        problem$mu,
        problem$exposure,
        .sigma = case$sigma,
        exposure_type = "continuous",
        stabilize = case$stabilize
      )

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
      by_string <- wt_ate(
        problem$mu,
        problem$exposure,
        .sigma = case$sigma,
        exposure_type = "continuous",
        stabilize = case$stabilize,
        .density = "normal"
      )

      by_spec <- wt_ate(
        problem$mu,
        problem$exposure,
        .sigma = case$sigma,
        exposure_type = "continuous",
        stabilize = case$stabilize,
        .density = dens_normal()
      )

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
