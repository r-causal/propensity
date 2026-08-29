# What `stabilize` resolves to when the caller does not name it. The default is
# read per exposure type rather than fixed once: a binary or categorical
# exposure is weighted unstabilized, as it always has been, and a continuous
# exposure is stabilized on its marginal density, which is the only form of the
# continuous weights the package recommends. An explicit `TRUE` or `FALSE` is
# honored wherever it was before, and an explicit `FALSE` on a continuous
# exposure still says that the weights it built are not the recommended ones.
#
# The resolution belongs to the numeric method, after the exposure type is
# known, so every route that funnels through it reads the same default. The
# model methods, the data frame method, and the modified-score methods are
# checked here for that reason rather than for the arithmetic they share with
# the numeric route, which is pinned in test-weights-continuous.R.

# ---- fixtures ---------------------------------------------------------------

# One problem per exposure type, all built from the same covariate so that the
# three defaults are read against comparable data: a dose and the fitted
# conditional means of a linear model for it, a binary exposure and its
# propensity score, and a three-level exposure and the matrix of level
# probabilities that go with it.
stabilize_default_data <- local({
  set.seed(20250903)

  n <- 120
  x <- rnorm(n)

  dose <- 0.7 * x + rnorm(n)
  mu <- as.numeric(fitted(lm(dose ~ x)))

  ps <- plogis(0.4 * x)
  trt <- rbinom(n, 1, ps)

  eta <- cbind(0, 0.5 * x, -0.3 * x)
  ps_matrix <- exp(eta) / rowSums(exp(eta))
  colnames(ps_matrix) <- c("a", "b", "c")
  levels3 <- factor(
    vapply(
      seq_len(n),
      function(i) sample(c("a", "b", "c"), 1, prob = ps_matrix[i, ]),
      character(1)
    ),
    levels = c("a", "b", "c")
  )

  list(
    n = n,
    x = x,
    dose = dose,
    mu = mu,
    ps = ps,
    trt = trt,
    ps_matrix = ps_matrix,
    levels3 = levels3
  )
})

# The conditional density the continuous weights divide by, spread by the
# pooled residual standard deviation, and the marginal density that stabilizes
# them, read at the exposure's own population moments.
stabilize_default_denominator <- function() {
  dat <- stabilize_default_data
  sigma <- sqrt(mean((dat$dose - dat$mu)^2))

  dnorm(dat$dose, mean = dat$mu, sd = sigma)
}

stabilize_default_numerator <- function() {
  dose <- stabilize_default_data$dose
  mu_a <- mean(dose)

  dnorm(dose, mean = mu_a, sd = sqrt(mean((dose - mu_a)^2)))
}

# ---- the continuous default -------------------------------------------------

test_that("a continuous exposure is stabilized when `stabilize` is not named", {
  dat <- stabilize_default_data

  # The recommended weights raise nothing to say so, so the alert is opted back
  # in to for this test and the call still has nothing to report.
  withr::local_options(propensity.quiet = FALSE)
  expect_no_message(
    weights <- wt_ate(dat$mu, dat$dose, exposure_type = "continuous")
  )

  expect_true(is_stabilized(weights))
  expect_identical(density_meta(weights)$numerator, "marginal")
  expect_equal(
    as.numeric(weights),
    stabilize_default_numerator() / stabilize_default_denominator(),
    tolerance = 1e-12
  )
})

test_that("the continuous default is the same weights as `stabilize = TRUE`", {
  dat <- stabilize_default_data

  expect_equal(
    wt_ate(dat$mu, dat$dose, exposure_type = "continuous"),
    wt_ate(dat$mu, dat$dose, exposure_type = "continuous", stabilize = TRUE)
  )
})

test_that("an explicit `stabilize = FALSE` still leaves a continuous exposure unstabilized", {
  dat <- stabilize_default_data

  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    weights <- wt_ate(
      dat$mu,
      dat$dose,
      exposure_type = "continuous",
      stabilize = FALSE
    ),
    "Using unstabilized weights for continuous exposures is not recommended."
  )

  expect_false(is_stabilized(weights))
  expect_identical(density_meta(weights)$numerator, "none")
  expect_equal(
    as.numeric(weights),
    1 / stabilize_default_denominator(),
    tolerance = 1e-12
  )
})

test_that("the continuous alert is raised once for one call", {
  dat <- stabilize_default_data

  withr::local_options(propensity.quiet = FALSE)
  messages <- character()
  withCallingHandlers(
    wt_ate(
      dat$mu,
      dat$dose,
      exposure_type = "continuous",
      stabilize = FALSE
    ),
    message = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }
  )

  expect_length(messages, 1L)
})

# ---- the binary and categorical defaults ------------------------------------

test_that("a binary exposure is unstabilized when `stabilize` is not named", {
  dat <- stabilize_default_data

  withr::local_options(propensity.quiet = FALSE)
  expect_no_message(
    weights <- wt_ate(dat$ps, dat$trt, exposure_type = "binary")
  )

  expect_false(is_stabilized(weights))
  expect_identical(
    weights,
    wt_ate(dat$ps, dat$trt, exposure_type = "binary", stabilize = FALSE)
  )
})

test_that("`stabilize = TRUE` on a binary exposure is unchanged", {
  dat <- stabilize_default_data

  weights <- wt_ate(
    dat$ps,
    dat$trt,
    exposure_type = "binary",
    stabilize = TRUE
  )

  p1 <- mean(dat$trt)
  expect_true(is_stabilized(weights))
  expect_equal(
    as.numeric(weights),
    dat$trt * p1 / dat$ps + (1 - dat$trt) * (1 - p1) / (1 - dat$ps),
    tolerance = 1e-12
  )
})

test_that("a categorical exposure is unstabilized when `stabilize` is not named", {
  dat <- stabilize_default_data

  withr::local_options(propensity.quiet = FALSE)
  expect_no_message(
    weights <- wt_ate(dat$ps_matrix, dat$levels3, exposure_type = "categorical")
  )

  expect_false(is_stabilized(weights))
  expect_identical(
    weights,
    wt_ate(
      dat$ps_matrix,
      dat$levels3,
      exposure_type = "categorical",
      stabilize = FALSE
    )
  )
})

test_that("`stabilize = TRUE` on a categorical exposure is unchanged", {
  dat <- stabilize_default_data

  weights <- wt_ate(
    dat$ps_matrix,
    dat$levels3,
    exposure_type = "categorical",
    stabilize = TRUE
  )

  expect_true(is_stabilized(weights))
  expect_false(identical(
    as.numeric(weights),
    as.numeric(wt_ate(
      dat$ps_matrix,
      dat$levels3,
      exposure_type = "categorical",
      stabilize = FALSE
    ))
  ))
})

# ---- censoring weights ------------------------------------------------------

test_that("censoring weights read the same per-type default", {
  dat <- stabilize_default_data

  withr::local_options(propensity.quiet = FALSE)
  expect_no_message(
    continuous <- wt_cens(dat$mu, dat$dose, exposure_type = "continuous")
  )

  expect_true(is_stabilized(continuous))
  expect_identical(estimand(continuous), "uncensored")
  expect_identical(density_meta(continuous)$numerator, "marginal")
  expect_equal(
    as.numeric(continuous),
    as.numeric(wt_cens(
      dat$mu,
      dat$dose,
      exposure_type = "continuous",
      stabilize = TRUE
    )),
    tolerance = 1e-12
  )

  binary <- wt_cens(dat$ps, dat$trt, exposure_type = "binary")
  expect_false(is_stabilized(binary))
  expect_identical(
    binary,
    wt_cens(dat$ps, dat$trt, exposure_type = "binary", stabilize = FALSE)
  )
})

test_that("an unstabilized continuous censoring weight still alerts", {
  dat <- stabilize_default_data

  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    weights <- wt_cens(
      dat$mu,
      dat$dose,
      exposure_type = "continuous",
      stabilize = FALSE
    ),
    "Using unstabilized weights for continuous exposures is not recommended."
  )

  expect_false(is_stabilized(weights))
  expect_identical(density_meta(weights)$numerator, "none")
})

# ---- the routes that funnel into the numeric method -------------------------

test_that("the model methods resolve the continuous default", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("MASS")

  dat <- stabilize_default_data
  model_data <- data.frame(dose = dat$dose, x = dat$x)

  fits <- list(
    glm = glm(dose ~ x, data = model_data, family = gaussian()),
    lm = lm(dose ~ x, data = model_data),
    gam = mgcv::gam(dose ~ s(x), data = model_data),
    rlm = MASS::rlm(dose ~ x, data = model_data)
  )

  for (route in names(fits)) {
    fit <- fits[[route]]

    weights <- wt_ate(fit, dat$dose, exposure_type = "continuous")

    expect_true(is_stabilized(weights), info = route)
    expect_identical(density_meta(weights)$numerator, "marginal", info = route)
    expect_equal(
      as.numeric(weights),
      as.numeric(wt_ate(
        fit,
        dat$dose,
        exposure_type = "continuous",
        stabilize = TRUE
      )),
      tolerance = 1e-12,
      info = route
    )
  }
})

test_that("the data frame method resolves the continuous default", {
  dat <- stabilize_default_data
  scores <- data.frame(mu = dat$mu)

  weights <- wt_ate(scores, dat$dose, exposure_type = "continuous")

  expect_true(is_stabilized(weights))
  expect_identical(density_meta(weights)$numerator, "marginal")
  expect_equal(
    as.numeric(weights),
    as.numeric(wt_ate(
      dat$mu,
      dat$dose,
      exposure_type = "continuous",
      stabilize = TRUE
    )),
    tolerance = 1e-12
  )
})

test_that("the modified-score methods resolve the continuous default", {
  dat <- stabilize_default_data

  # A score a modification can be applied to has to lie in (0, 1), so the
  # conditional mean these routes are handed is one that does. What is being
  # read here is the default, not the arithmetic, which is the same as the
  # numeric route's on the same numbers.
  scores <- plogis(0.5 * dat$x)
  fit <- glm(
    trt ~ x,
    data = data.frame(trt = dat$trt, x = dat$x),
    family = binomial()
  )

  modified <- list(
    trim = ps_refit(
      ps_trim(scores, method = "ps", lower = 0.2, upper = 0.8),
      model = fit
    ),
    trunc = ps_trunc(scores, method = "ps", lower = 0.2, upper = 0.8),
    calib = ps_calibrate(scores, dat$trt)
  )

  for (route in names(modified)) {
    ps <- modified[[route]]

    weights <- wt_ate(ps, dat$dose, exposure_type = "continuous")

    expect_true(is_stabilized(weights), info = route)
    expect_identical(density_meta(weights)$numerator, "marginal", info = route)
    expect_equal(
      as.numeric(weights),
      as.numeric(wt_ate(
        ps,
        dat$dose,
        exposure_type = "continuous",
        stabilize = TRUE
      )),
      tolerance = 1e-12,
      info = route
    )
  }
})

test_that("the modified-score methods leave a binary exposure unstabilized", {
  dat <- stabilize_default_data

  truncated <- ps_trunc(dat$ps, method = "ps", lower = 0.05, upper = 0.95)
  weights <- wt_ate(truncated, dat$trt, exposure_type = "binary")

  expect_false(is_stabilized(weights))
  expect_identical(
    weights,
    wt_ate(truncated, dat$trt, exposure_type = "binary", stabilize = FALSE)
  )
})

# ---- ipw() ------------------------------------------------------------------

test_that("ipw() reads a continuous fit built with the default as stabilized", {
  set.seed(20250904)
  n <- 400
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  y <- 1 + 0.6 * A + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  dat <- data.frame(x1 = x1, x2 = x2, A = A, y = y)

  ps_mod <- lm(A ~ x1 + x2, data = dat)
  fitted_ps <- as.double(fitted(ps_mod))

  wts_default <- wt_ate(fitted_ps, dat$A, exposure_type = "continuous")
  wts_stabilized <- wt_ate(
    fitted_ps,
    dat$A,
    exposure_type = "continuous",
    stabilize = TRUE
  )

  msm_default <- lm(y ~ A, data = dat, weights = wts_default)
  msm_stabilized <- lm(y ~ A, data = dat, weights = wts_stabilized)

  res_default <- ipw(ps_mod, msm_default)
  res_stabilized <- ipw(ps_mod, msm_stabilized)

  # The marginal numerator is two more parameters in the stacked estimating
  # equations, so a fit built on the default now carries them.
  theta_names <- names(coef(res_default$fit))
  expect_true("mu_a" %in% theta_names)
  expect_true("sigma2_a" %in% theta_names)

  expect_equal(res_default$estimates, res_stabilized$estimates)
})
