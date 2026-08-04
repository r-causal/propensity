# Base method coverage for `ipw` results. `coef()`, `vcov()`, `confint()`,
# `nobs()`, `df.residual()`, and `weights()` are causalgenerics' methods on the
# shared result class, and they read only the fields that class contracts to
# hold. What propensity owes them is the covariance of the effects it reports,
# attached to the estimates table as the `ipw_vcov` attribute.
#
# That attribute is a square matrix of the reported effects, in the row order of
# `estimates` and labeled the way `print()` labels its rows: the effect measure
# on its own, or the effect measure and the comparison together where a
# categorical exposure reports one row per comparison. Its diagonal is the square
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

# ---- model fitting -----------------------------------------------------------

# Each fixture returns the weights beside the two models, because the weights are
# themselves one of the things the accessors report and comparing `weights()`
# against the model frame it is read from would assert nothing.

# Binary exposure: a logit propensity score model and an exposure-only weighted
# outcome model. The outcome model carries no covariate, so the same pair of
# models serves the M-estimation and the linearization paths, the latter being
# restricted to marginal outcome models. `estimand` picks the weighting scheme,
# so a fit reporting something other than the ATE is available here rather than
# through a second fixture.
fit_base_binary_models <- function(
  dat,
  outcome_family = "binomial",
  estimand = "ate"
) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wt_fun <- if (estimand == "att") wt_att else wt_ate
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_fun(ps_mod)
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
fit_base_continuous_models <- function(dat, outcome_family = "gaussian") {
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
  fmla <- stats::reformulate("A", response = outcome_var)
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
  if (is.null(estimates[["comparison"]])) {
    estimates$effect
  } else {
    paste(estimates$effect, estimates$comparison)
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
  expect_identical(rownames(covariance), c("rd", "log(rr)", "log(or)"))

  # The three effect measures are increasing transformations of the same pair of
  # marginal means, so they move together and the matrix is not the diagonal one
  # the standard errors on their own would make.
  expect_true(all(covariance[upper.tri(covariance)] > 0))
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
  expect_identical(dimnames(covariance), list("diff", "diff"))
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
  expect_identical(rownames(covariance), c("rd", "log(rr)", "log(or)"))
})

test_that("a categorical fit labels its covariance by effect and comparison", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_base_categorical()
  mods <- fit_base_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  # Three effect measures for each of two comparisons, and neither part of a
  # label identifies a row on its own, so the label is both together.
  covariance <- expect_ipw_vcov_contract(res)
  labels <- c(
    "rd b vs a",
    "log(rr) b vs a",
    "log(or) b vs a",
    "rd c vs a",
    "log(rr) c vs a",
    "log(or) c vs a"
  )
  expect_identical(rownames(covariance), labels)
  expect_identical(anyDuplicated(labels), 0L)

  # The comparisons share the reference level's marginal mean, so the covariance
  # spans the whole table rather than being block diagonal by comparison.
  across <- covariance[seq(1, 3), seq(4, 6)]
  expect_true(all(is.finite(across)))
  expect_true(all(across != 0))
})

test_that("a categorical att fit records the covariance of its effects", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_base_categorical()
  mods <- fit_base_categorical_models(dat, estimand = "att")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_identical(res$estimand, "att")
  covariance <- expect_ipw_vcov_contract(res)
  expect_identical(dim(covariance), c(6L, 6L))
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
  expect_identical(rownames(covariance), c("rd", "log(rr)", "log(or)"))

  # Each of the three influence functions is a difference of the same two per-arm
  # influence values, with a positive coefficient on each, so every pair of
  # effect measures has a real covariance and a positive one. Zeros here would
  # report the effect measures as independent estimates of unrelated things.
  off_diagonal <- covariance[upper.tri(covariance)]
  expect_length(off_diagonal, 3L)
  expect_true(all(is.finite(off_diagonal)))
  expect_true(all(off_diagonal > 0))
})

test_that("a gaussian linearization fit records one variance", {
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  expect_null(res$fit)
  covariance <- expect_ipw_vcov_contract(res)
  expect_identical(dimnames(covariance), list("diff", "diff"))
})

# ---- accessor integration ----------------------------------------------------

test_that("the accessors read a binary mestimation fit", {
  skip_if_not_installed("deli")
  dat <- sim_base_binary()
  mods <- fit_base_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  expect_ipw_accessor_contract(res, mods$wts)
  expect_identical(names(coef(res)), c("rd", "log(rr)", "log(or)"))

  # The counts are of the stacked system: every observation of the fixture went
  # into it, and the parameters it solves for are what the observations are
  # spent on.
  expect_identical(nobs(res), 600L)
  expect_identical(df.residual(res), 600L - as.integer(res$fit@n_params))
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

test_that("the accessors label a categorical fit by effect and comparison", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_base_categorical()
  mods <- fit_base_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_ipw_accessor_contract(res, mods$wts)

  labels <- paste(res$estimates$effect, res$estimates$comparison)
  expect_identical(names(coef(res)), labels)
  expect_identical(anyDuplicated(names(coef(res))), 0L)
  expect_identical(rownames(confint(res)), labels)
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
