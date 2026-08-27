# Base method coverage for `ipw` results. `coef()`, `vcov()`, `confint()`,
# `nobs()`, `df.residual()`, and `weights()` are causalgenerics' methods on the
# shared result class, and they read only the fields that class contracts to
# hold. What propensity owes them is the covariance of the effects it reports,
# attached to the estimates table as the `ipw_vcov` attribute.
#
# That attribute is a square matrix of the reported effects, in the row order of
# `estimates` and labeled the way `print()` labels its rows: the effect measure
# and the contrast it names together, which for the counterfactual mean rows is
# the exposure level the mean belongs to. Its diagonal is the square
# of the standard errors the result already reports, which is what keeps the
# matrix and the estimates table describing the same fit. Its off-diagonal is the
# part that cannot be recovered from the estimates table, and it is a real
# quantity: the effect measures of one fit are transformations of the same
# marginal means and move together.
#
# M-estimation stacks the propensity score model, the outcome model, and the
# effect measures into one system, so every block of the resulting sandwich
# accounts for the others having been estimated from the same data. The
# covariance of a component model there is not the one its own fitting routine
# reports, and `outcome_mod` and `wt_mod` carry the joint blocks by way of the
# `ipw_model` subclass. Linearization solves no such system, so its component
# models stay exactly as they were fit and only the effects covariance is
# attached.
#
# Coverage spans every route that builds an `ipw` object: a binary glm exposure
# under both standard error methods, a categorical exposure through
# nnet::multinom, and a continuous exposure through lm.

# ---- data simulators ---------------------------------------------------------

sim_base_binary <- function(seed = 2024, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.6 * x2))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * x1 - 0.3 * x2))
  yc <- 1.5 + 0.8 * z + 0.6 * x1 - 0.4 * x2 + rnorm(n)
  data.frame(x1, x2, z, y, yc)
}

sim_base_categorical <- function(seed = 2024, n = 700) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  eta_b <- -0.2 + 0.5 * x1 + 0.3 * x2
  eta_c <- 0.1 - 0.4 * x1 + 0.6 * x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  pa <- 1 / denom
  pb <- exp(eta_b) / denom
  u <- runif(n)
  a <- factor(
    ifelse(u < pa, "a", ifelse(u < pa + pb, "b", "c")),
    levels = c("a", "b", "c")
  )
  y <- rbinom(
    n,
    1,
    plogis(-0.3 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.5 * x1)
  )
  data.frame(x1, x2, a, y)
}

sim_base_continuous <- function(seed = 2024, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  yc <- 1 + 0.6 * A + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  yb <- rbinom(n, 1, plogis(-0.5 + 0.4 * A + 0.3 * x1 - 0.2 * x2))
  data.frame(x1, x2, A, yc, yb)
}

# ---- reported labels ---------------------------------------------------------

# The rows a binary fit reports, in order: the counterfactual mean at each
# exposure level with the reference level first, then the effect measures, each
# naming the contrast it compares. The fixtures code the exposure 0/1, so those
# are the level names the labels are written with.
binary_effect_labels <- function(contrasts = c("rd", "log(rr)", "log(or)")) {
  c("mean 0", "mean 1", paste(contrasts, "1 vs 0"))
}

# ---- model fitting -----------------------------------------------------------

# Each fixture returns the weights beside the two models, because the weights are
# themselves one of the things the accessors report and comparing `weights()`
# against the model frame it is read from would assert nothing.

# Binary exposure: a logit propensity score model and an exposure-only weighted
# outcome model. The outcome model carries no covariate, so the same pair of
# models serves the M-estimation and the linearization paths, the latter being
# restricted to marginal outcome models. `estimand` picks the weighting scheme,
# so a fit reporting something other than the ATE is available here rather than
# through a second fixture. `stabilize` scales the ATE weights by the marginal
# exposure probability, which the stacked system solves for as a parameter of
# its own rather than holding fixed. Stabilization is defined for the ATE alone,
# so the ATT branch takes no such argument.
fit_base_binary_models <- function(
  dat,
  outcome_family = "binomial",
  estimand = "ate",
  stabilize = FALSE
) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    if (estimand == "att") {
      wt_att(ps_mod)
    } else {
      wt_ate(ps_mod, stabilize = stabilize)
    }
  )
  outcome_var <- if (outcome_family == "binomial") "y" else "yc"
  fmla <- stats::reformulate("z", response = outcome_var)
  outcome_mod <- if (outcome_family == "binomial") {
    glm(
      fmla,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  } else {
    lm(fmla, data = dat, weights = wts)
  }
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# Categorical exposure: a multinomial propensity score model and an
# exposure-only weighted outcome model. The tightened `reltol` matches the
# convention in the other categorical test files, keeping the fitted
# coefficients at the score root the stacked system re-solves. The ATT of a
# categorical exposure is defined against one level, so that fixture names the
# focal level the weights target.
fit_base_categorical_models <- function(dat, estimand = "ate") {
  ps_mod <- nnet::multinom(
    a ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_named <- unname(predict(ps_mod, type = "probs"))
  colnames(ps_named) <- ps_mod$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    if (estimand == "att") {
      wt_att(
        ps_named,
        dat$a,
        exposure_type = "categorical",
        .focal_level = "a"
      )
    } else {
      wt_ate(ps_named, dat$a, exposure_type = "categorical")
    }
  )
  outcome_mod <- glm(
    y ~ a,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# Continuous exposure: an lm propensity score model of the exposure, stabilized
# continuous ATE weights, and a weighted marginal structural model with the one
# exposure term the continuous path requires.
fit_base_continuous_models <- function(
  dat,
  outcome_family = "gaussian",
  msm_rhs = "A"
) {
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(fitted(ps_mod)),
      dat$A,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )
  outcome_var <- if (outcome_family == "binomial") "yb" else "yc"
  fmla <- stats::reformulate(msm_rhs, response = outcome_var)
  outcome_mod <- if (outcome_family == "binomial") {
    glm(
      fmla,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  } else {
    lm(fmla, data = dat, weights = wts)
  }
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# ---- expected shape ----------------------------------------------------------

# The labels every surface of a result keys its rows by. `print()` labels its
# rows this way and the accessors name their coefficients this way, so the tests
# read the labels off the estimates table rather than repeating a list of them
# that could agree with the wrong thing.
base_effect_labels <- function(result) {
  estimates <- result$estimates
  if (is.null(estimates[["contrast"]])) {
    estimates$effect
  } else {
    paste(estimates$effect, estimates$contrast)
  }
}

# The covariance contract, asserted against the result the matrix is attached to:
# a square numeric matrix of the reported effects, labeled on both margins by
# the effect labels, symmetric, and carrying the reported standard errors on its
# diagonal. Symmetry and the diagonal are compared rather than matched exactly,
# because a sandwich is assembled from a product of three matrices rather than
# read off a symmetric form, and the reported standard errors are square roots of
# the variances rather than the variances themselves.
expect_ipw_vcov_contract <- function(result) {
  labels <- base_effect_labels(result)
  covariance <- attr(result$estimates, "ipw_vcov", exact = TRUE)

  expect_true(is.matrix(covariance))
  expect_type(covariance, "double")
  expect_identical(dim(covariance), c(length(labels), length(labels)))
  expect_identical(dimnames(covariance), list(labels, labels))

  # `vcov()` reports the attribute itself rather than anything derived from it,
  # and `coef()` names its elements the way the matrix names its margins, so the
  # two are usable together without a reordering step in between.
  expect_identical(vcov(result), covariance)
  expect_identical(names(coef(result)), rownames(covariance))
  expect_identical(names(coef(result)), colnames(covariance))

  expect_true(all(is.finite(covariance)))
  expect_equal(covariance, t(covariance))
  expect_equal(unname(diag(covariance)), result$estimates$std.err^2)

  invisible(covariance)
}

# The accessors that read fields the result class has always held, asserted
# against the result and against the weights the outcome model was fit with.
# Nothing here consults the covariance, so a result that carries none still
# answers all of it.
expect_ipw_accessor_contract <- function(result, wts) {
  estimates <- result$estimates
  labels <- base_effect_labels(result)

  expect_identical(coef(result), stats::setNames(estimates$estimate, labels))

  # At the level the result was fit at, `confint()` reports the bounds the result
  # stores rather than a recomputation that merely agrees with them.
  interval <- confint(result, level = estimates$conf.level[[1]])
  expect_true(is.matrix(interval))
  expect_identical(rownames(interval), labels)
  expect_identical(unname(interval[, 1]), estimates$ci.lower)
  expect_identical(unname(interval[, 2]), estimates$ci.upper)

  # The two counts describe the fit rather than its estimates, which is what
  # `glance()` reports, so the two surfaces report the same numbers.
  glanced <- glance(result)
  expect_identical(nobs(result), glanced[["nobs"]])
  expect_identical(df.residual(result), glanced[["df.residual"]])

  # The weights arrive as they were stored, so the estimand they record and the
  # class that carries it travel with them rather than being flattened to a
  # numeric copy of their values.
  expect_identical(weights(result), wts)
  expect_s3_class(weights(result), "psw")
  expect_identical(estimand(weights(result)), result$estimand)

  invisible(result)
}

# The component model contract on the M-estimation path: the model the result
# holds is the model the user fit, with the block of the joint sandwich that
# belongs to it attached, so `vcov()` on it accounts for the rest of the system
# having been estimated from the same data. `naive` is the covariance the model's
# own fitting routine reports, which the corrected block must not be. The
# comparison is made on the parameters the two have in common, because the
# propensity score block of a continuous exposure carries one parameter the
# model itself has no coefficient for.
expect_ipw_model_contract <- function(model, labels, naive) {
  expect_s3_class(model, "ipw_model")

  covariance <- vcov(model)
  expect_true(is.matrix(covariance))
  expect_identical(dimnames(covariance), list(labels, labels))
  expect_true(all(is.finite(covariance)))
  expect_equal(covariance, t(covariance))
  expect_true(all(diag(covariance) > 0))

  shared <- rownames(naive)
  expect_false(isTRUE(all.equal(
    covariance[shared, shared],
    naive,
    check.attributes = FALSE
  )))

  invisible(covariance)
}

# ---- the covariance attribute, M-estimation ---------------------------------

test_that("a binary mestimation fit records the covariance of its effects", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  covariance <- expect_ipw_vcov_contract(res)
  expect_identical(rownames(covariance), binary_effect_labels())

  # The three effect measures are increasing transformations of the same pair of
  # marginal means, so they move together and that block of the matrix is not
  # the diagonal one the standard errors on their own would make.
  contrasts <- covariance[3:5, 3:5]
  expect_true(all(contrasts[upper.tri(contrasts)] > 0))

  # A mean row covaries with the contrasts built from it, and with the opposite
  # arm's mean, so no row of the matrix is the standard error alone.
  expect_true(all(covariance[1:2, 3:5] != 0))
})

test_that("a gaussian outcome under mestimation records one variance", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # A continuous outcome reports one effect measure, so the matrix is one by one
  # and has no covariance to hold. It is still a labeled matrix, because that is
  # what `vcov()` returns whatever the fit reports.
  covariance <- expect_ipw_vcov_contract(res)
  labels <- c("mean 0", "mean 1", "diff 1 vs 0")
  expect_identical(dimnames(covariance), list(labels, labels))
})

test_that("a binary att fit records the covariance of its effects", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat, estimand = "att")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # The estimand changes the weights and so the whole stacked system. The
  # covariance is of whatever effects the result reports, so it follows.
  expect_identical(res$estimand, "att")
  covariance <- expect_ipw_vcov_contract(res)
  expect_identical(rownames(covariance), binary_effect_labels())
})

test_that("a stabilized binary fit records the covariance of its effects", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat, stabilize = TRUE)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # Stabilizing the weights adds a parameter to the stacked system, so the
  # sandwich the covariance is read off is a wider matrix than the unstabilized
  # fit's. What the estimates table is owed is unchanged: the covariance of the
  # three effect measures it reports, carrying their standard errors.
  expect_true(is_stabilized(mods$wts))
  covariance <- expect_ipw_vcov_contract(res)
  expect_identical(rownames(covariance), binary_effect_labels())
})

test_that("a categorical fit labels its covariance by effect and contrast", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_base_categorical()
  mods <- fit_base_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  # One mean for each of three levels, then three effect measures for each of
  # two contrasts, and neither part of a contrast label identifies a row on its
  # own, so the label is both together.
  covariance <- expect_ipw_vcov_contract(res)
  labels <- c(
    "mean a",
    "mean b",
    "mean c",
    "rd b vs a",
    "log(rr) b vs a",
    "log(or) b vs a",
    "rd c vs a",
    "log(rr) c vs a",
    "log(or) c vs a"
  )
  expect_identical(rownames(covariance), labels)
  expect_identical(anyDuplicated(labels), 0L)

  # The contrasts share the reference level's marginal mean, so the covariance
  # spans the whole table rather than being block diagonal by contrast.
  across <- covariance[seq(4, 6), seq(7, 9)]
  expect_true(all(is.finite(across)))
  expect_true(all(across != 0))

  # The mean rows are in the covariance too, and each contrast covaries with the
  # means it is built from.
  expect_true(all(is.finite(covariance[seq(1, 3), seq(4, 9)])))
  expect_true(all(covariance[seq(1, 3), seq(4, 9)] != 0))
})

test_that("a categorical att fit records the covariance of its effects", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_base_categorical()
  mods <- fit_base_categorical_models(dat, estimand = "att")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_identical(res$estimand, "att")
  covariance <- expect_ipw_vcov_contract(res)
  expect_identical(dim(covariance), c(9L, 9L))
})

test_that("a continuous exposure fit records the variance of its slope", {
  skip_if_not_installed("deli")
  dat <- sim_base_continuous()
  mods <- fit_base_continuous_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  # A continuous exposure reports the one coefficient of its marginal structural
  # model, so the matrix is one by one and named for that effect measure.
  covariance <- expect_ipw_vcov_contract(res)
  expect_identical(dimnames(covariance), list("slope", "slope"))
})

test_that("a continuous logit fit records the variance of its log odds ratio", {
  skip_if_not_installed("deli")
  dat <- sim_base_continuous()
  mods <- fit_base_continuous_models(dat, outcome_family = "binomial")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  covariance <- expect_ipw_vcov_contract(res)
  expect_identical(dimnames(covariance), list("log(or)", "log(or)"))
})

# ---- the covariance attribute, linearization --------------------------------

test_that("a binary linearization fit records the covariance of its effects", {
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # This path stores no stacked fit, so the covariance has to come from the
  # influence functions the standard errors themselves come from. The contract it
  # answers is the same one.
  expect_null(res$fit)
  covariance <- expect_ipw_vcov_contract(res)
  expect_identical(rownames(covariance), binary_effect_labels())

  # Each of the three influence functions is a difference of the same two per-arm
  # influence values, with a positive coefficient on each, so every pair of
  # effect measures has a real covariance and a positive one. Zeros here would
  # report the effect measures as independent estimates of unrelated things.
  contrasts <- covariance[3:5, 3:5]
  off_diagonal <- contrasts[upper.tri(contrasts)]
  expect_length(off_diagonal, 3L)
  expect_true(all(is.finite(off_diagonal)))
  expect_true(all(off_diagonal > 0))

  # The per-arm means are reported rows of their own here too, built from the
  # same two influence values, so they covary with the contrasts.
  expect_true(all(is.finite(covariance[1:2, 3:5])))
  expect_true(all(covariance[1:2, 3:5] != 0))
})

test_that("a gaussian linearization fit records one variance", {
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  expect_null(res$fit)
  covariance <- expect_ipw_vcov_contract(res)
  labels <- c("mean 0", "mean 1", "diff 1 vs 0")
  expect_identical(dimnames(covariance), list(labels, labels))
})

# ---- accessor integration ----------------------------------------------------

test_that("the accessors read a binary mestimation fit", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  expect_ipw_accessor_contract(res, mods$wts)
  expect_identical(names(coef(res)), binary_effect_labels())

  # The counts are of the stacked system: every observation of the fixture went
  # into it, and the parameters it solves for are what the observations are
  # spent on.
  expect_identical(nobs(res), 600L)
  expect_identical(df.residual(res), 600L - as.integer(res$fit@n_params))
})

test_that("a stabilized fit spends an observation on its stabilizer", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  stabilized <- fit_base_binary_models(dat, stabilize = TRUE)
  unstabilized <- fit_base_binary_models(dat)
  expect_true(is_stabilized(stabilized$wts))
  expect_false(is_stabilized(unstabilized$wts))

  res <- ipw(
    stabilized$ps_mod,
    stabilized$outcome_mod,
    se_method = "mestimation"
  )
  base <- ipw(
    unstabilized$ps_mod,
    unstabilized$outcome_mod,
    se_method = "mestimation"
  )

  # The weights the result reports are the stabilized vector the outcome model
  # was fit with, so the stabilization travels with the fit rather than living
  # only in the fixture.
  expect_identical(weights(res), stabilized$wts)

  # The stabilizer is solved for rather than held fixed, so it is a parameter of
  # the stacked system and sits between the propensity score block and the
  # outcome block. Everything else the system solves for is what the
  # unstabilized fit solves for, in the same order.
  theta_names <- names(res$fit@theta)
  expect_identical(theta_names[[4]], "stab_pi")
  expect_identical(theta_names[-4], names(base$fit@theta))

  # An observation is spent on that parameter, so the counts describing the fit
  # report one fewer residual degree of freedom than the unstabilized fit of the
  # same data does.
  expect_identical(as.integer(res$fit@n_params), 11L)
  expect_identical(nobs(res), 600L)
  expect_identical(df.residual(res), 589L)
  expect_identical(df.residual(res), df.residual(base) - 1L)

  # The row that parameter answers is the mean of the centered exposure, whose
  # root is the marginal exposure probability, so the solved value is the
  # proportion exposed rather than anything the weights were scaled by
  # approximately. The solver reaches it to within two units in the last place, a
  # measured 1.1e-16 on this fixture, so the comparison carries a tolerance
  # rather than asking for the same double.
  expect_equal(res$fit@theta[["stab_pi"]], mean(dat$z), tolerance = 1e-12)
})

test_that("the accessors read a binary linearization fit", {
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  expect_ipw_accessor_contract(res, mods$wts)

  # There is no stacked system here, so there is no parameter count to spend the
  # observations against and the degrees of freedom are reported as missing
  # rather than invented. The observations are the outcome model's.
  expect_identical(nobs(res), 600L)
  expect_identical(df.residual(res), NA_integer_)
})

test_that("the accessors label a categorical fit by effect and contrast", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_base_categorical()
  mods <- fit_base_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_ipw_accessor_contract(res, mods$wts)

  labels <- paste(res$estimates$effect, res$estimates$contrast)
  expect_identical(names(coef(res)), labels)
  expect_identical(anyDuplicated(names(coef(res))), 0L)
  expect_identical(rownames(confint(res)), labels)

  # The count is of the fixture rather than of the six rows the effects table
  # holds, which the agreement with `glance()` above does not distinguish.
  expect_identical(nobs(res), 700L)
})

test_that("the accessors read a continuous exposure fit", {
  skip_if_not_installed("deli")
  dat <- sim_base_continuous()
  mods <- fit_base_continuous_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_ipw_accessor_contract(res, mods$wts)
  expect_identical(names(coef(res)), "slope")
  expect_identical(nobs(res), 600L)
})

test_that("confint() rebuilds the interval at a level the fit does not store", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # The interval at a level other than the stored one is the normal
  # approximation the result itself reports, which is the interval `tidy()`
  # builds for the same request, so the two surfaces agree to the bit.
  interval <- confint(res, level = 0.9)
  tidied <- tidy(res, conf.int = TRUE, conf.level = 0.9)
  expect_identical(colnames(interval), c("5 %", "95 %"))
  expect_identical(unname(interval[, 1]), tidied$conf.low)
  expect_identical(unname(interval[, 2]), tidied$conf.high)

  # A 90% interval is not the stored 95% one read back.
  expect_true(all(interval[, 1] > res$estimates$ci.lower))
  expect_true(all(interval[, 2] < res$estimates$ci.upper))
})

test_that("confint() reports the stored bounds of a fit built at any level", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    se_method = "mestimation",
    conf_level = 0.9
  )

  # The stored level is whatever the fit was built at, so a request for it
  # returns the result's own numbers.
  at_stored <- confint(res, level = 0.9)
  expect_identical(unname(at_stored[, 1]), res$estimates$ci.lower)
  expect_identical(unname(at_stored[, 2]), res$estimates$ci.upper)

  # The 0.95 default is the accessor's own rather than the fit's, so it is
  # recomputed rather than served from the stored 0.9 columns.
  at_default <- confint(res)
  expect_identical(colnames(at_default), c("2.5 %", "97.5 %"))
  expect_true(all(at_default[, 1] < res$estimates$ci.lower))
  expect_true(all(at_default[, 2] > res$estimates$ci.upper))
})

test_that("weights() returns the psw vector the outcome model was fit with", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat, estimand = "att")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # The weights are the analysis' central per-observation quantity, and the
  # estimand they target is part of what they are, so they come back as the
  # object they were supplied as rather than as its values.
  wts <- weights(res)
  expect_identical(wts, mods$wts)
  expect_s3_class(wts, "psw")
  expect_identical(estimand(wts), "att")
  expect_length(wts, 600L)
})

# ---- component models, M-estimation ------------------------------------------

test_that("a binary mestimation fit corrects its component covariances", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  expect_ipw_model_contract(
    res$outcome_mod,
    names(coef(mods$outcome_mod)),
    vcov(mods$outcome_mod)
  )
  expect_ipw_model_contract(
    res$wt_mod,
    names(coef(mods$ps_mod)),
    vcov(mods$ps_mod)
  )

  # The models the caller fit are their own and are left as they are: the
  # correction travels with the copies the result holds.
  expect_identical(class(mods$outcome_mod), c("glm", "lm"))
  expect_identical(class(mods$ps_mod), c("glm", "lm"))
  expect_null(attr(mods$outcome_mod, "ipw_vcov", exact = TRUE))

  # Wrapping adds a class and an attribute and changes nothing else, so the model
  # the result holds still does what a model does.
  expect_identical(coef(res$outcome_mod), coef(mods$outcome_mod))
  newdata <- data.frame(z = c(0, 1))
  expect_identical(
    predict(res$outcome_mod, newdata = newdata, type = "response"),
    predict(mods$outcome_mod, newdata = newdata, type = "response")
  )
})

test_that("a categorical fit names its propensity block by level and term", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_base_categorical()
  mods <- fit_base_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_ipw_model_contract(
    res$outcome_mod,
    names(coef(mods$outcome_mod)),
    vcov(mods$outcome_mod)
  )

  # A multinomial model has one coefficient per term per non-reference level, and
  # `nnet::multinom()` names the pair with a colon between them. The propensity
  # block of the stacked system holds the same parameters in the same order, so
  # it carries the same names.
  ps_labels <- rownames(vcov(mods$ps_mod))
  expect_identical(
    ps_labels,
    c("b:(Intercept)", "b:x1", "b:x2", "c:(Intercept)", "c:x1", "c:x2")
  )
  expect_ipw_model_contract(res$wt_mod, ps_labels, vcov(mods$ps_mod))

  # The propensity score model is read for its predictions elsewhere in the
  # package, and wrapping leaves that alone.
  expect_identical(
    predict(res$wt_mod, type = "probs"),
    predict(mods$ps_mod, type = "probs")
  )
})

test_that("a continuous fit carries the error variance in its propensity block", {
  skip_if_not_installed("deli")
  dat <- sim_base_continuous()
  mods <- fit_base_continuous_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_ipw_model_contract(
    res$outcome_mod,
    names(coef(mods$outcome_mod)),
    vcov(mods$outcome_mod)
  )

  # The weights of a continuous exposure are built from a conditional density,
  # which needs the spread of the propensity score model's errors as well as its
  # coefficients. That spread is a parameter of the stacked system, so the
  # propensity block is one row wider than the model's own covariance.
  ps_labels <- c(names(coef(mods$ps_mod)), "sigma2_d")
  expect_ipw_model_contract(res$wt_mod, ps_labels, vcov(mods$ps_mod))
})

# ---- component models, linearization -----------------------------------------

test_that("a linearization fit leaves its component models as they were fit", {
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # This path solves nothing jointly: the standard errors come from influence
  # functions computed after both models were fit, and there is no sandwich with
  # a block for either of them. There is therefore no corrected covariance to
  # attach, and a wrapper claiming one would be claiming a correction that was
  # never made.
  expect_false(inherits(res$outcome_mod, "ipw_model"))
  expect_false(inherits(res$wt_mod, "ipw_model"))
  expect_identical(class(res$outcome_mod), c("glm", "lm"))
  expect_identical(class(res$wt_mod), c("glm", "lm"))
  expect_null(attr(res$outcome_mod, "ipw_vcov", exact = TRUE))
  expect_null(attr(res$wt_mod, "ipw_vcov", exact = TRUE))

  expect_identical(vcov(res$outcome_mod), vcov(mods$outcome_mod))
  expect_identical(vcov(res$wt_mod), vcov(mods$ps_mod))
})

# ---- a result that carries no covariance -------------------------------------

test_that("a result with no covariance answers every accessor but vcov()", {
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # A result built before the covariance was part of the class contract, or by a
  # fitting package that has not adopted it. Its `fit` is a bare list, which is
  # what the contract allows a fitting package to store there and what the
  # generics that reach into it have nothing to dispatch on.
  estimates <- res$estimates
  attr(estimates, "ipw_vcov") <- NULL
  bare <- causalgenerics::new_ipw(
    estimand = res$estimand,
    wt_mod = mods$ps_mod,
    outcome_mod = mods$outcome_mod,
    estimates = estimates,
    se_method = "mestimation",
    fit = list(theta = estimates$estimate, vcov = diag(nrow(estimates)))
  )

  # The one thing that is missing is reported as missing, under a class that
  # names the kind of result it was missing from, rather than substituted for
  # with a diagonal built from the standard errors.
  expect_error(vcov(bare), class = "causalgenerics_no_vcov_ipw")

  # Nothing else is affected, because nothing else reads it.
  expect_identical(coef(bare), coef(res))
  expect_identical(confint(bare), confint(res))
  expect_identical(nobs(bare), nobs(res))
  expect_identical(df.residual(bare), NA_integer_)
  expect_identical(weights(bare), mods$wts)

  # `glance()` describes the fit, and a fit it can read nothing from leaves the
  # residual degrees of freedom missing rather than failing the whole summary.
  glanced <- glance(bare)
  expect_identical(nrow(glanced), 1L)
  expect_named(glanced, c("estimand", "nobs", "df.residual"))
  expect_identical(glanced[["estimand"]], "ate")
  expect_identical(glanced[["nobs"]], nobs(bare))
  expect_identical(glanced[["df.residual"]], NA_integer_)
})

# ---- the presentation mode ---------------------------------------------------

# A result reports its effects in one of two readings, recorded in the `effects`
# field the result class contracts to hold. The marginal reading is the causal
# contrast estimates every route reported before the field existed; the
# conditional reading presents the outcome model's coefficient surface. Both
# surfaces exist on every result, so the field says which one is presented rather
# than which one was computed, and `ipw()` takes an `effects` argument that says
# which one the result it builds should record.
#
# What propensity owes the two readings is the argument on every route and the
# corrected covariance the conditional one reports. The reading itself, the
# accessors that report it, and the generics that move a result between the two
# are causalgenerics'.

# One result from each of the four routes that builds an `ipw` object: a binary
# glm exposure under each standard error method, a categorical exposure through
# nnet::multinom, and a continuous exposure through lm. The tests that assert
# something of every route fit their models once here rather than once per route.
fit_base_all_routes <- function() {
  binary <- fit_base_binary_models(sim_base_binary())
  categorical <- fit_base_categorical_models(sim_base_categorical())
  continuous <- fit_base_continuous_models(sim_base_continuous())

  list(
    mestimation = ipw(
      binary$ps_mod,
      binary$outcome_mod,
      se_method = "mestimation"
    ),
    linearization = ipw(
      binary$ps_mod,
      binary$outcome_mod,
      se_method = "linearization"
    ),
    categorical = ipw(categorical$ps_mod, categorical$outcome_mod),
    continuous = ipw(continuous$ps_mod, continuous$outcome_mod)
  )
}

# The result an `effects` argument builds, against the same result built without
# the argument and moved to that reading afterwards. The field records which
# reading is presented rather than which one was computed, so naming a reading at
# construction settles the field and nothing else: the same estimates, the same
# covariance attached to them, and the same component models, wrapper included.
# Asserting only the field and the surface it selects would leave a construction
# that computed just the named reading, skipping the marginal estimates or the
# covariance because they are not the ones on show, passing every test here.
#
# `fit` is excluded because two solves are not comparable. The M-estimator deli
# returns carries the closures of the call that produced it, each with its own
# environment, so two solves of one system agree in every number they report and
# are identical in none of them. What the exclusion gives up is asserted against
# the default build elsewhere in this file.
expect_ipw_built_as <- function(result, expected) {
  expect_identical(result$effects, expected$effects)

  # Before the comparison below, which draws its fields from `expected` and would
  # otherwise pass over a field the construction dropped entirely.
  expect_identical(names(result), names(expected))

  comparable <- setdiff(names(expected), "fit")
  expect_identical(result[comparable], expected[comparable])

  invisible(result)
}

test_that("every route defaults to the marginal reading and round-trips", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  results <- fit_base_all_routes()

  # Marginal is what every route reported before the mode was a field, so a call
  # that names no mode still reports it.
  expect_identical(
    vapply(results, function(res) res$effects, character(1)),
    stats::setNames(rep("marginal", length(results)), names(results))
  )

  # A result says which readings it supports as well as which one it presents.
  # A binary exposure under either standard error method, a categorical one, and
  # a single-term continuous one all have both, which is what makes the flips
  # below answers rather than refusals.
  expect_identical(
    lapply(results, function(res) res$readings),
    stats::setNames(
      rep(list(c("marginal", "conditional")), length(results)),
      names(results)
    )
  )

  # Both readings exist on every result, so moving to the other one records the
  # move and reads nothing else, and moving back is the result that went in
  # rather than a rebuild of it.
  flipped <- lapply(results, causalgenerics::as_conditional)
  expect_identical(
    vapply(flipped, function(res) res$effects, character(1)),
    stats::setNames(rep("conditional", length(results)), names(results))
  )
  expect_identical(lapply(flipped, causalgenerics::as_marginal), results)

  # Asking for the reading a result already records says what asking once said.
  expect_identical(lapply(flipped, causalgenerics::as_conditional), flipped)
  expect_identical(lapply(results, causalgenerics::as_marginal), results)
})

test_that("a binary mestimation fit records the reading it was built in", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)

  base <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
  res <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    se_method = "mestimation",
    effects = "conditional"
  )

  # Building in the conditional reading is building the result and recording the
  # reading, so it is the default build moved to that reading and nothing else.
  # The stacked system is solved either way, so the corrected block still reaches
  # the component models and the marginal estimates are still there to present.
  expect_ipw_built_as(res, causalgenerics::as_conditional(base))
  expect_s3_class(res$outcome_mod, "ipw_model")
  expect_identical(vcov(res), vcov(res$outcome_mod))

  # The reading is a field of the result rather than something the accessors are
  # told at each call, so a result built in the conditional one reports it with
  # nothing named at the call site.
  expect_identical(coef(res), coef(mods$outcome_mod))

  # Naming the default is the other half of the argument, and it has to build the
  # result that leaving the argument out builds.
  named <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    se_method = "mestimation",
    effects = "marginal"
  )
  expect_ipw_built_as(named, causalgenerics::as_marginal(base))
  expect_identical(names(coef(named)), binary_effect_labels())
})

test_that("a binary linearization fit records the reading it was built in", {
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)

  base <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")
  res <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    se_method = "linearization",
    effects = "conditional"
  )

  # The argument records the reading on a path that solves no stacked system, so
  # it is threaded through the construction rather than through whatever the
  # standard errors were computed by.
  expect_ipw_built_as(res, causalgenerics::as_conditional(base))
  expect_null(res$fit)
  expect_identical(coef(res), coef(mods$outcome_mod))

  # Nothing was solved jointly here whichever reading was asked for, so this
  # stays the one route whose component models are exactly as they were fit. A
  # construction that wrapped them because the conditional reading wants a
  # covariance would be claiming a correction that was never computed.
  expect_false(inherits(res$outcome_mod, "ipw_model"))
  expect_false(inherits(res$wt_mod, "ipw_model"))
})

test_that("a categorical fit records the reading it was built in", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_base_categorical()
  mods <- fit_base_categorical_models(dat)

  res <- ipw(mods$ps_mod, mods$outcome_mod, effects = "conditional")

  # The conditional reading is of the outcome model, which has one coefficient
  # per non-reference exposure level rather than the six rows the marginal
  # reading of this fit reports.
  expect_identical(res$effects, "conditional")
  expect_identical(coef(res), coef(mods$outcome_mod))
  expect_identical(names(coef(res)), c("(Intercept)", "ab", "ac"))
})

test_that("a continuous fit records the reading it was built in", {
  skip_if_not_installed("deli")
  dat <- sim_base_continuous()
  mods <- fit_base_continuous_models(dat)

  base <- ipw(mods$ps_mod, mods$outcome_mod)
  res <- ipw(mods$ps_mod, mods$outcome_mod, effects = "conditional")

  # This route stacks the estimating equations too, so the same thing has to hold
  # of it: the argument records the reading and leaves the rest of the result
  # where the default build puts it, corrected block included.
  expect_ipw_built_as(res, causalgenerics::as_conditional(base))
  expect_s3_class(res$outcome_mod, "ipw_model")
  expect_identical(coef(res), coef(mods$outcome_mod))
  expect_identical(names(coef(res)), c("(Intercept)", "A"))
})

test_that("a dose basis fit records the conditional reading", {
  skip_if_not_installed("deli")
  dat <- sim_base_continuous()
  mods <- fit_base_continuous_models(dat, msm_rhs = c("A", "I(A^2)"))

  # An exposure entering the outcome model through several design columns has no
  # coefficient that is an effect, so there is no marginal reading of the fit to
  # present. The result records the conditional one and declares it as the only
  # reading it supports, which is what the refusals come from; the class is the
  # shared one every other result carries.
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  expect_identical(class(res), "ipw")
  expect_identical(res$effects, "conditional")
  expect_identical(res$readings, "conditional")

  # Naming the reading the default records builds the result the default builds,
  # which is what says the argument records a reading rather than selecting a
  # computation.
  named <- ipw(mods$ps_mod, mods$outcome_mod, effects = "conditional")
  expect_ipw_built_as(named, res)

  # The accessors answer from the outcome model with nothing named at the call
  # site, under the block the stacked system left for it.
  expect_s3_class(res$outcome_mod, "ipw_model")
  expect_identical(coef(res), coef(mods$outcome_mod))
  expect_identical(names(coef(res)), c("(Intercept)", "A", "I(A^2)"))
  expect_identical(vcov(res), vcov(res$outcome_mod))
  expect_identical(rownames(confint(res)), names(coef(mods$outcome_mod)))
})

test_that("a dose basis fit refuses the marginal reading everywhere", {
  skip_if_not_installed("deli")
  dat <- sim_base_continuous()
  mods <- fit_base_continuous_models(dat, msm_rhs = c("A", "I(A^2)"))
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  # Asking for the reading at construction is this package's refusal, since the
  # result that would record it is never built.
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, effects = "marginal"),
    class = "propensity_ipw_effects_error"
  )

  # Every surface a reading can be asked for at refuses the one this result does
  # not declare, under the classes the result class raises for it, so a caller
  # reaching for it from any of them is told the same thing rather than being
  # handed a table of coefficients under the name of an effect.
  cls <- "causalgenerics_unsupported_reading_marginal"
  expect_error(causalgenerics::as_marginal(res), class = cls)
  expect_error(coef(res, effects = "marginal"), class = cls)
  expect_error(vcov(res, effects = "marginal"), class = cls)
  expect_error(confint(res, effects = "marginal"), class = cls)
  expect_error(as.data.frame(res, effects = "marginal"), class = cls)
  expect_error(tidy(res, effects = "marginal"), class = cls)

  # The general class is what a caller catching any refusal of a reading writes
  # against, so it is on the condition beside the specific one.
  expect_s3_class(
    expect_error(coef(res, effects = "marginal")),
    "causalgenerics_unsupported_reading"
  )

  # The reading the result records is still asked for the ordinary way, and
  # asking for it says what recording it said.
  expect_identical(causalgenerics::as_conditional(res), res)
  expect_identical(coef(res, effects = "conditional"), coef(res))
  expect_identical(vcov(res, effects = "conditional"), vcov(res))

  # The accessors that describe the fit rather than a reading of it take no
  # reading at all, so they are untouched by any of this.
  expect_identical(nobs(res), nobs(mods$outcome_mod))
  expect_identical(weights(res), mods$wts)
  expect_identical(estimand(res), "ate")
})

test_that("an invalid effects value errors", {
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)

  err <- expect_error(
    ipw(
      mods$ps_mod,
      mods$outcome_mod,
      se_method = "linearization",
      effects = "both"
    ),
    class = "rlang_error"
  )

  # The message names the two readings, which is what tells a rejection of the
  # value apart from a rejection of the argument itself.
  expect_match(conditionMessage(err), "marginal", fixed = TRUE)
  expect_match(conditionMessage(err), "conditional", fixed = TRUE)
})

test_that("the effects argument reports the other reading for one call", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  expect_identical(res$effects, "marginal")
  expect_identical(coef(res, effects = "conditional"), coef(mods$outcome_mod))
  expect_identical(vcov(res, effects = "conditional"), vcov(res$outcome_mod))

  # Naming a reading at the call site answers in it and leaves the result where
  # it was, so the next call with nothing named answers in the stored one.
  expect_identical(res$effects, "marginal")
  expect_identical(
    coef(res),
    stats::setNames(res$estimates$estimate, base_effect_labels(res))
  )
})

test_that("a mestimation fit reports its corrected block conditionally", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
  conditional <- causalgenerics::as_conditional(res)

  # The covariance the conditional reading reports is the outcome block of the
  # stacked sandwich, which this path attaches to the model itself, rather than
  # anything derived from it here.
  covariance <- vcov(conditional)
  expect_identical(covariance, vcov(res$outcome_mod))
  expect_identical(
    dimnames(covariance),
    list(
      names(coef(mods$outcome_mod)),
      names(coef(mods$outcome_mod))
    )
  )

  # It is not the covariance the outcome model computed for itself, which treats
  # the estimated weights as fixed and reports an uncertainty the coefficients do
  # not have.
  expect_false(isTRUE(all.equal(
    covariance,
    vcov(mods$outcome_mod),
    check.attributes = FALSE
  )))

  # The limits are built from that block at every level, since the ones the
  # result stores belong to the effects the other reading reports.
  interval <- confint(conditional)
  expect_identical(rownames(interval), names(coef(mods$outcome_mod)))
  expect_identical(colnames(interval), c("2.5 %", "97.5 %"))
})

test_that("a continuous fit reports its corrected block conditionally", {
  skip_if_not_installed("deli")
  dat <- sim_base_continuous()
  mods <- fit_base_continuous_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  conditional <- causalgenerics::as_conditional(res)

  # This route stacks the estimating equations too, so its outcome model carries
  # the joint block and the conditional reading has a covariance to report.
  expect_s3_class(res$outcome_mod, "ipw_model")
  expect_identical(vcov(conditional), vcov(res$outcome_mod))
  expect_identical(
    rownames(vcov(conditional)),
    names(coef(mods$outcome_mod))
  )
})

test_that("a linearization fit has no covariance for the conditional reading", {
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")
  conditional <- causalgenerics::as_conditional(res)

  # Nothing was solved jointly here, so there is no corrected block to report and
  # the model's own covariance is not a substitute for one. The absence is
  # reported as an absence rather than filled in.
  expect_error(
    vcov(conditional),
    class = "causalgenerics_no_conditional_vcov"
  )
  expect_error(
    confint(conditional),
    class = "causalgenerics_no_conditional_vcov"
  )

  # The coefficients need no such block and come back in either reading.
  expect_identical(coef(conditional), coef(mods$outcome_mod))
})

test_that("a conditional mestimation result prints the model's coefficients", {
  skip_if_not_installed("deli")

  # testthat 3e pins the output width but not the number of significant digits,
  # and `printCoefmat()` wraps its table past 80 columns under a larger `digits`,
  # which splits the rows this test reads by position.
  withr::local_options(digits = 7)

  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
  out <- capture.output(print(causalgenerics::as_conditional(res)))

  # The reading is named twice, because the two readings are different tables of
  # different numbers: once beside the estimand and once over the table.
  expect_true(any(grepl(
    "Effects: conditional (outcome model)",
    out,
    fixed = TRUE
  )))
  expect_false(any(grepl("Marginal estimates:", out, fixed = TRUE)))

  header <- which(startsWith(out, "Conditional estimates (outcome model):"))
  expect_length(header, 1L)

  # The corrected block is there, so the coefficients are reported beside the
  # standard errors it implies rather than on their own.
  expect_true(grepl("Std. Error", out[[header + 1L]], fixed = TRUE))
  expect_false(any(grepl(
    "Standard errors are not reported",
    out,
    fixed = TRUE
  )))

  # The rows are the outcome model's coefficients rather than the effect
  # measures the marginal reading of the same result tabulates.
  labels <- names(coef(mods$outcome_mod))
  rows <- out[seq(header + 2L, header + 1L + length(labels))]
  expect_identical(sub(" .*$", "", rows), labels)
})

test_that("a conditional linearization result prints a note in place of errors", {
  withr::local_options(digits = 7)

  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")
  conditional <- causalgenerics::as_conditional(res)

  # `print()` is the view of whatever a caller is holding, so a result whose
  # conditional reading has no covariance to report is still printable. Refusing
  # here would leave a result that cannot be looked at at all.
  out <- expect_no_error(capture.output(print(conditional)))

  header <- which(startsWith(out, "Conditional estimates (outcome model):"))
  expect_length(header, 1L)

  # The coefficients on their own, under a note saying what is missing, rather
  # than beside the standard errors the model computed for itself: a column of
  # those under this heading would be read as the corrected ones.
  expect_identical(trimws(out[[header + 1L]]), "Estimate")
  expect_false(any(grepl("Std. Error", out, fixed = TRUE)))
  expect_true(any(grepl(
    "Standard errors are not reported",
    out,
    fixed = TRUE
  )))
  expect_length(grep("covariance", out, fixed = TRUE), 1L)
})
