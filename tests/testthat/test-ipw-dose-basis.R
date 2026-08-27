# Basis-matrix dose-response marginal structural models on the plain continuous
# path: an ordinary single-dose fit with stabilized continuous ate weights and
# an lm or glm marginal structural model whose exposure term is a basis such as
# `poly(e, 2)`, `splines::ns(e, 3)`, or `splines::bs(e, 3)`. No joint container,
# no `.by`; this is the route a dose-response curve has always taken, asked for
# a curve written as a basis rather than as `e + I(e^2)`.
#
# Three things follow from the boundary being variable membership. A basis term
# reads the exposure alone, so it is admitted however many columns it expands
# to; each of those columns carries a coefficient, so each is a row of the
# surface the stacked system estimates, named by the coefficient rather than by
# the term, since the two differ exactly here. A basis term reading a covariate
# reads a covariate, so it is admitted as one and contributes no row. A term
# reading both, such as `poly(e, 2):x1`, keeps the refusal every mixed term
# gets: its coefficient is a change in the effect per unit of the covariate, and
# no row could name the covariate value it holds at.
#
# What such a fit reports is the conditional reading and only that. An exposure
# entering through several columns has no coefficient that is an effect: a curve
# has a different slope at every dose, so no row of the table answers the
# question the marginal reading is asked. `ipw()` therefore records the
# conditional reading for these fits, says so once, and refuses the marginal one
# wherever it is asked for. Marginalizing over the dose is a separate estimand
# that this package does not compute, and the message names the marginaleffects
# package rather than leaving the user to assemble it from the coefficients.
#
# The one thing a basis changes about how the fit has to be called is `.data`. A
# model frame records the term, not the variables inside it, so `model.frame()`
# of a `poly(e, 2)` fit has no `e` column and the exposure cannot be read off
# it. That is already refused, with a message directing the user to `.data`, and
# supplying `.data` is the whole of the remedy: the outcome design is read from
# the fit rather than rebuilt, and the propensity design is rebuilt through the
# terms object, whose `predvars` records the basis each of `poly`, `ns`, and
# `bs` was fit with.

# ---- data simulator ---------------------------------------------------------

# One dose confounded by two covariates, with a curved response so a second
# basis coefficient is not zero, and a binary outcome carrying the same curve
# for the logit variant.
sim_dose_basis <- function(seed = 8301, n = 800) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  e <- 0.4 + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  y <- 1 + 0.5 * e + 0.2 * e^2 + 0.4 * x1 + rnorm(n)
  yb <- rbinom(n, 1, plogis(-0.3 + 0.5 * e - 0.15 * e^2))
  data.frame(x1, x2, e, y, yb)
}

quiet_dose_wt <- function(expr) {
  withr::with_options(list(propensity.quiet = TRUE), expr)
}

# ---- model fitting ----------------------------------------------------------

# The single-dose route with whatever marginal structural model a test wants to
# hand it. The weights are the stabilized continuous ate weights this path
# documents, kept as a `psw` object so the estimand reaches the outcome model
# frame. The quasibinomial fit tightens its tolerance so its coefficients sit at
# the weighted maximum likelihood estimate well below the comparison tolerance,
# and uses `quasibinomial()` because non-integer successes warn under
# `binomial()`.
fit_dose_basis <- function(
  dat,
  msm_rhs,
  outcome_family = c("gaussian", "binomial"),
  ps_rhs = c("x1", "x2")
) {
  outcome_family <- match.arg(outcome_family)

  ps_mod <- lm(stats::reformulate(ps_rhs, response = "e"), data = dat)
  wts <- quiet_dose_wt(wt_ate(
    as.double(fitted(ps_mod)),
    dat$e,
    exposure_type = "continuous",
    stabilize = TRUE
  ))

  response <- if (outcome_family == "binomial") "yb" else "y"
  fmla <- stats::reformulate(msm_rhs, response = response)
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

# ---- shared expectations ----------------------------------------------------

# The column contract of the multi-row coefficient surface the stacked system
# estimates. It is the frame the result class contracts to hold rather than a
# reading the result presents: a basis puts a contrast column on the table
# naming the coefficient each row belongs to, placed where every other route
# places it, and puts no group column there, since nothing about a
# dose-response curve is evaluated in a stratum.
expect_dose_basis_columns <- function(estimates) {
  testthat::expect_identical(
    names(estimates),
    c(
      "effect",
      "contrast",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  testthat::expect_null(estimates[["group"]])

  invisible(estimates)
}

# The reading a basis fit presents, written once because the claim is the same
# for every basis. The result records the conditional reading and declares it
# as the only one it supports, which is where the refusals of the other one come
# from. Its accessors key by the outcome model's coefficient names in one order,
# and the covariance they report is the outcome block of the stacked sandwich
# rather than the one the weighted fit computed for itself while treating its
# weights as fixed.
#
# The coefficients come out of one fit of one model, so they covary; a surface
# assembled row by row would report zeros off the diagonal.
#
# The stored frame is checked against that same block. It is not a reading the
# result presents, but `pool_ipw()` reads it rather than the accessors, so its
# standard errors reach a pooled result whether or not anything else reports
# them, and nothing else pins them to a covariance.
expect_dose_basis_conditional <- function(res, outcome_mod) {
  labels <- names(coef(outcome_mod))

  testthat::expect_identical(res$effects, "conditional")
  testthat::expect_identical(res$readings, "conditional")
  testthat::expect_identical(coef(res), coef(outcome_mod))
  testthat::expect_identical(names(coef(res)), labels)
  testthat::expect_identical(dimnames(vcov(res)), list(labels, labels))
  testthat::expect_identical(rownames(confint(res)), labels)
  testthat::expect_identical(tidy(res)$term, labels)

  covariance <- vcov(res)
  testthat::expect_identical(covariance, vcov(res$outcome_mod))
  testthat::expect_true(all(is.finite(sqrt(diag(covariance)))))
  testthat::expect_true(all(diag(covariance) > 0))
  testthat::expect_equal(covariance, t(covariance), tolerance = 1e-12)

  entry <- covariance[2L, 3L]
  testthat::expect_true(is.finite(entry))
  testthat::expect_gt(abs(entry), 1e-8)

  testthat::expect_false(isTRUE(all.equal(
    covariance,
    vcov(outcome_mod),
    check.attributes = FALSE
  )))

  # The stored frame reports the standard error of each exposure coefficient,
  # which is the square root of that coefficient's diagonal entry in the same
  # block, read at the coefficient the row names rather than by position. A
  # frame built from the wrong rows of the sandwich, or from the covariance the
  # weighted fit computed for itself, differs here by more than rounding.
  estimates <- res$estimates
  testthat::expect_equal(
    unname(sqrt(diag(covariance))[estimates$contrast]),
    estimates$std.err,
    tolerance = 1e-12
  )
  testthat::expect_true(all(is.finite(estimates$std.err)))
  testthat::expect_true(all(estimates$std.err > 0))
  testthat::expect_true(all(estimates$ci.lower < estimates$ci.upper))

  invisible(res)
}

# One refusal of a reading a result does not declare, asserted by its classes
# rather than by its words, since the words are causalgenerics'. The specific
# class names the reading that was asked for and the general one is what a
# caller catching any such refusal writes against.
expect_unsupported_marginal <- function(expr) {
  err <- testthat::expect_error(
    expr,
    class = "causalgenerics_unsupported_reading_marginal"
  )
  testthat::expect_s3_class(err, "causalgenerics_unsupported_reading")

  invisible(err)
}

# Every surface that could report the marginal reading of a basis fit refuses
# it, so a caller who reaches for it from any of them is told the same thing.
# The refusal belongs to the result class rather than to this package: the fit
# declares the readings it supports and the shared accessors answer for them,
# which is what keeps every other result answering as it always did without a
# subclass here to carry refusals of its own. `tidy()` is propensity's own
# method, so it is asserted here too rather than assumed to defer.
expect_dose_basis_refuses_marginal <- function(res) {
  expect_unsupported_marginal(causalgenerics::as_marginal(res))
  expect_unsupported_marginal(coef(res, effects = "marginal"))
  expect_unsupported_marginal(vcov(res, effects = "marginal"))
  expect_unsupported_marginal(confint(res, effects = "marginal"))
  expect_unsupported_marginal(tidy(res, effects = "marginal"))
  expect_unsupported_marginal(as.data.frame(res, effects = "marginal"))

  invisible(res)
}

# ---- one row per basis coefficient ------------------------------------------

test_that("a polynomial basis reports one row per basis coefficient", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  est <- res$estimates

  expect_dose_basis_columns(est)
  expect_identical(nrow(est), 2L)
  expect_identical(est$effect, c("coef", "coef"))

  # Named by the coefficient, not by the term: one term produced both columns,
  # so the term names neither of them on its own.
  expect_identical(est$contrast, c("poly(e, 2)1", "poly(e, 2)2"))
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[c("poly(e, 2)1", "poly(e, 2)2")]),
    tolerance = 1e-10
  )
  expect_dose_basis_conditional(res, fx$outcome_mod)

  # The stacked system for a stabilized single-dose fit: three propensity score
  # coefficients, the conditional variance the density ratio needs, the two
  # moments the stabilizing numerator needs, and the outcome coefficients. Only
  # the last block moves, and it moves by the basis dimension.
  expect_identical(as.integer(res$fit@n_params), 3L + 1L + 2L + 3L)

  # One exposure column is one slope, so a bare dose is the marginal reading it
  # always was, on the class it always had, declaring both readings as every
  # other result does.
  bare <- fit_dose_basis(dat, "e")
  bare_res <- ipw(bare$ps_mod, bare$outcome_mod)
  expect_identical(as.integer(bare_res$fit@n_params), 3L + 1L + 2L + 2L)
  expect_identical(bare_res$effects, "marginal")
  expect_identical(bare_res$readings, c("marginal", "conditional"))
  expect_identical(class(bare_res), "ipw")
})

test_that("a natural spline basis reports one row per basis coefficient", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()

  # No bare `e` term anywhere in this model. The boundary is variable
  # membership, not the presence of a linear term, so a model that is nothing
  # but a spline in the exposure is as reportable as one that is not.
  fx <- fit_dose_basis(dat, "splines::ns(e, 3)")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  est <- res$estimates

  expect_dose_basis_columns(est)
  expect_identical(nrow(est), 3L)
  expect_identical(est$effect, rep("coef", 3L))
  expect_identical(
    est$contrast,
    c("splines::ns(e, 3)1", "splines::ns(e, 3)2", "splines::ns(e, 3)3")
  )
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[est$contrast]),
    tolerance = 1e-10
  )
  expect_dose_basis_conditional(res, fx$outcome_mod)

  # Three basis columns where the bare dose has one, so the outcome block is
  # four coefficients wide and the system is two parameters wider than the
  # single-term fit's eight.
  expect_identical(as.integer(res$fit@n_params), 3L + 1L + 2L + 4L)
  expect_identical(nobs(res), 800L)
  expect_identical(df.residual(res), 800L - 10L)
})

test_that("a B-spline basis reports one row per basis coefficient", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "splines::bs(e, 3)")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  est <- res$estimates

  expect_dose_basis_columns(est)
  expect_identical(nrow(est), 3L)
  expect_identical(est$effect, rep("coef", 3L))
  expect_identical(
    est$contrast,
    c("splines::bs(e, 3)1", "splines::bs(e, 3)2", "splines::bs(e, 3)3")
  )
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[est$contrast]),
    tolerance = 1e-10
  )
  expect_dose_basis_conditional(res, fx$outcome_mod)

  # A cubic B-spline with no interior knots spans the same three columns a
  # natural spline of the same argument does, so the system is the same width.
  expect_identical(as.integer(res$fit@n_params), 3L + 1L + 2L + 4L)
})

# ---- the reading a basis fit reports ----------------------------------------

test_that("a basis fit announces the reading it reports", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)")

  # The announcement is the one thing about this route a user cannot infer from
  # the result they were handed: the surface looks like a coefficient table
  # either way, and what changed is which question it answers. It says which
  # reading was recorded, why there is no other, where the marginalization it
  # does not do belongs, and how to stop being told.
  withr::with_options(
    list(propensity.quiet = FALSE),
    expect_snapshot(res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat))
  )
})

test_that("naming the conditional reading silences the announcement", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)")

  # A caller who has named the reading has been told, so the message is a
  # default being explained rather than a fact being reported and it stops.
  named <- withr::with_options(
    list(propensity.quiet = FALSE),
    expect_no_message(
      ipw(fx$ps_mod, fx$outcome_mod, .data = dat, effects = "conditional")
    )
  )

  # Naming the reading the default records builds the result the default builds.
  # `fit` is excluded because two solves of one system agree in every number
  # they report and are identical in none of them: the M-estimator carries the
  # closures of the call that produced it.
  default <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  expect_identical(class(named), class(default))
  expect_identical(names(named), names(default))
  comparable <- setdiff(names(default), "fit")
  expect_identical(named[comparable], default[comparable])
})

test_that("a basis fit refuses the marginal reading", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)")

  # Asking for the marginal reading of a model that has none is a question
  # rather than a preference, so it is answered rather than quietly given the
  # other reading.
  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, .data = dat, effects = "marginal"),
    class = "propensity_ipw_effects_error"
  )
  expect_propensity_error(
    ipw(fx$ps_mod, fx$outcome_mod, .data = dat, effects = "marginal")
  )
})

test_that("the marginal reading is refused wherever it is asked for", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "splines::ns(e, 3)")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)

  # One refusal reached from every door, so a caller who goes around the
  # constructor gets the same answer the constructor gave.
  expect_dose_basis_refuses_marginal(res)

  # Asking for the reading the result already records is answered rather than
  # refused, which is what makes the refusal about the reading rather than about
  # the argument.
  expect_identical(causalgenerics::as_conditional(res), res)
  expect_identical(coef(res, effects = "conditional"), coef(res))
  expect_identical(tidy(res, effects = "conditional"), tidy(res))
})

test_that("a basis fit is an ipw result reporting its conditional reading", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)

  # The class is the shared one, with no subclass of this package's on it. What
  # marks the fit out is the set of readings it declares, which the result class
  # contracts to hold and its accessors answer from, so everything written
  # against an `ipw` result keeps working.
  expect_identical(class(res), "ipw")
  expect_identical(res$readings, "conditional")
  expect_identical(res$estimand, "ate")
  expect_identical(res$se_method, "mestimation")
  expect_s3_class(res$outcome_mod, "ipw_model")

  # The printed form is the conditional one every route prints, named over the
  # table as well as beside the estimand.
  withr::local_options(digits = 7)
  out <- capture.output(print(res))
  expect_true(any(grepl(
    "Effects: conditional (outcome model)",
    out,
    fixed = TRUE
  )))
  expect_false(any(grepl("Marginal estimates:", out, fixed = TRUE)))

  header <- which(startsWith(out, "Conditional estimates (outcome model):"))
  expect_length(header, 1L)
  labels <- names(coef(fx$outcome_mod))
  rows <- out[seq(header + 2L, header + 1L + length(labels))]
  expect_identical(substr(rows, 1L, nchar(labels)), labels)
})

# ---- the closed form ---------------------------------------------------------

test_that("the conditional reading is the outcome block of the stacked system", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)

  # An identity-link marginal structural model is the estimator: nothing is
  # standardized, so what the reading reports is exactly the weighted fit's
  # coefficient vector and the comparison is a closed form rather than an
  # agreement to within a g-computation. The intercept is a coefficient of that
  # model, so it is there, which is the difference between a coefficient surface
  # and the effects table a single-dose fit reports.
  beta <- coef(fx$outcome_mod)
  expect_equal(coef(res), beta, tolerance = 1e-10)
  expect_identical(names(coef(res))[[1L]], "(Intercept)")

  # The covariance is the block the stacked system leaves for the outcome
  # model, read out of the fitted variance object by position: three propensity
  # score coefficients, the conditional variance, and the two stabilizing
  # moments come first, and the outcome coefficients are the rest.
  covariance <- vcov(res$fit)
  outcome_block <- seq(7L, nrow(covariance))
  expect_equal(
    unname(vcov(res)),
    unname(covariance[outcome_block, outcome_block]),
    tolerance = 1e-12
  )

  # Row by row as well, so a reordering of the surface fails here rather than
  # passing on a vector that happens to line up.
  pick <- function(label) unname(coef(res)[[label]])
  expect_equal(pick("poly(e, 2)1"), unname(beta[["poly(e, 2)1"]]))
  expect_equal(pick("poly(e, 2)2"), unname(beta[["poly(e, 2)2"]]))
})

test_that("a logit marginal structural model reports the same rows on the odds scale", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)", outcome_family = "binomial")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  est <- res$estimates

  # The scale follows the outcome link, which is what the continuous path
  # already does with a single-term logit dose response. The rows the stacked
  # system estimates are otherwise the ones the gaussian fit reports, and the
  # contrast still names the coefficient.
  expect_dose_basis_columns(est)
  expect_identical(nrow(est), 2L)
  expect_identical(est$effect, c("log(or)", "log(or)"))
  expect_identical(est$contrast, c("poly(e, 2)1", "poly(e, 2)2"))
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[c("poly(e, 2)1", "poly(e, 2)2")]),
    tolerance = 1e-10
  )
  expect_dose_basis_conditional(res, fx$outcome_mod)
  expect_dose_basis_refuses_marginal(res)
  expect_identical(as.integer(res$fit@n_params), 3L + 1L + 2L + 3L)
})

test_that("the conditional reading of a basis fit is the outcome model's", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "splines::ns(e, 3)")
  conditional <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)

  # Read through the tidier, which is the surface the reading presents. The
  # conditional reading is the whole coefficient vector, so the intercept the
  # stored frame drops is back.
  tidied <- tidy(conditional)
  expect_identical(tidied$term, names(coef(fx$outcome_mod)))
  expect_equal(
    tidied$estimate,
    unname(coef(fx$outcome_mod)),
    tolerance = 1e-10
  )
  expect_false("contrast" %in% names(tidied))
  expect_false("group" %in% names(tidied))
})

# ---- pooling -----------------------------------------------------------------

test_that("pool_ipw() over basis fits pools the conditional reading", {
  skip_if_not_installed("deli")

  # Three fits of the same model to three datasets, which is the shape a set of
  # imputations arrives in. What is pooled is the reading the results record,
  # so a set of basis fits pools by coefficient.
  fits <- lapply(c(8301, 8302, 8303), function(seed) {
    dat <- sim_dose_basis(seed = seed)
    fx <- fit_dose_basis(dat, "poly(e, 2)")
    ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  })

  pooled <- pool_ipw(fits)
  coefficients <- names(coef(fits[[1L]]$outcome_mod))

  expect_identical(pooled$effects, "conditional")
  expect_identical(pooled$estimates$effect, coefficients)
  expect_identical(pooled$m, 3L)
  expect_true(all(is.finite(pooled$estimates$estimate)))
  expect_true(all(is.finite(pooled$estimates$std.err)))

  # The readings the results declare travel into the pooling, so a set of fits
  # with one reading between them pools that reading and records why the other
  # one is absent. Asking the pooled result for it is refused where the pooled
  # class is defined; what this package settles is the reading the pooling
  # starts from.
  expect_error(
    causalgenerics::as_marginal(pooled),
    class = "causalgenerics_pool_missing_surface_marginal"
  )
})

# ---- the .data requirement ---------------------------------------------------

test_that("a basis fit without .data directs the user to supply it", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()

  # A model frame records the term, so a basis fit has no exposure column in it
  # and the exposure cannot be read off the frame. The refusal is the one the
  # package already raises for a transformed outcome model, and the remedy it
  # names is the whole of what is needed.
  for (rhs in c("poly(e, 2)", "splines::ns(e, 3)", "splines::bs(e, 3)")) {
    fx <- fit_dose_basis(dat, rhs)
    expect_error(
      ipw(fx$ps_mod, fx$outcome_mod),
      class = "propensity_columns_exist_error",
      regexp = "Please specify"
    )
    expect_s3_class(ipw(fx$ps_mod, fx$outcome_mod, .data = dat), "ipw")
  }
})

test_that("the refusal without .data names the exposure and the remedy", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)")

  # The message names the column the frame holds instead of the exposure, which
  # is what tells a model that transforms the exposure from one that never reads
  # it. The first is what `.data` is for and the second is a model to refit, and
  # a message that said only that the exposure was missing left the reader to
  # work out which of the two they had.
  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod),
    class = "propensity_columns_exist_error",
    regexp = "poly\\(e, 2\\)"
  )
  expect_propensity_error(ipw(fx$ps_mod, fx$outcome_mod))

  # A basis of the exposure beside a basis of a covariate: only the term built
  # from the exposure is named, since the other one is not what the frame is
  # missing.
  beside <- fit_dose_basis(dat, "poly(e, 2) + splines::ns(x1, 3)")
  err <- expect_error(
    ipw(beside$ps_mod, beside$outcome_mod),
    class = "propensity_columns_exist_error"
  )
  expect_match(conditionMessage(err), "poly(e, 2)", fixed = TRUE)
  expect_false(grepl("ns(x1, 3)", conditionMessage(err), fixed = TRUE))
})

test_that("a spline propensity model rebuilds from .data as it was fit", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()

  # The one design that is rebuilt rather than read off a fit. `.data` is passed
  # through the terms object, whose `predvars` records the knots and the
  # orthogonalization the basis was fit with, so the rebuilt columns are the
  # fitted ones. A rebuild that refit the basis instead would give different
  # propensity scores, whose weights would fail the consistency preflight rather
  # than quietly changing the answer.
  fx <- fit_dose_basis(
    dat,
    "poly(e, 2)",
    ps_rhs = c("splines::ns(x1, 3)", "x2")
  )
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)

  expect_equal(
    res$estimates$estimate,
    unname(coef(fx$outcome_mod)[c("poly(e, 2)1", "poly(e, 2)2")]),
    tolerance = 1e-10
  )

  # Five propensity score coefficients now, and the rest of the system as
  # before.
  expect_identical(as.integer(res$fit@n_params), 5L + 1L + 2L + 3L)
  expect_dose_basis_conditional(res, fx$outcome_mod)
})

# ---- basis terms that are not exposure terms ---------------------------------

test_that("a basis term reading a covariate is admitted as a covariate", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()

  # `splines::ns(x1, 3)` reads a covariate and nothing else, so it is a
  # covariate term however many columns it expands to: it contributes no row and
  # the surface is the single-row one a bare dose reports, with no contrast
  # column to name a coefficient that is the only one.
  #
  # That single row keeps the word `slope`. The exposure enters through one
  # column, so its coefficient is the slope of the dose response everywhere and
  # the word is true of it; this is the vocabulary surface, and `slope` is
  # reserved for exactly the rows it is true of. The reading is marginal for the
  # same reason, and nothing is announced, since there is a row that answers the
  # question the marginal reading is asked.
  fx <- fit_dose_basis(dat, "e + splines::ns(x1, 3)")
  res <- withr::with_options(
    list(propensity.quiet = FALSE),
    expect_no_message(ipw(fx$ps_mod, fx$outcome_mod, .data = dat))
  )
  est <- res$estimates

  expect_identical(class(res), "ipw")
  expect_identical(res$effects, "marginal")
  expect_identical(res$readings, c("marginal", "conditional"))
  expect_identical(nrow(est), 1L)
  expect_identical(est$effect, "slope")
  expect_null(est[["contrast"]])
  expect_null(est[["group"]])
  expect_length(names(est), 8L)
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[["e"]]),
    tolerance = 1e-10
  )

  # The reading is the one every single-dose fit records, so it is reported
  # rather than refused and the result moves between the two as any other does.
  expect_identical(names(coef(res)), "slope")
  expect_identical(causalgenerics::as_marginal(res), res)
  expect_identical(
    coef(causalgenerics::as_conditional(res)),
    coef(fx$outcome_mod)
  )

  # The covariate basis is still in the outcome block of the stacked system,
  # which is five coefficients wide: an intercept, the dose, and three spline
  # columns.
  expect_identical(as.integer(res$fit@n_params), 3L + 1L + 2L + 5L)
})

test_that("an exposure basis beside a covariate basis reports only its own rows", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2) + splines::ns(x1, 3)")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  est <- res$estimates

  # Both terms expand to several columns and only one of them reads the
  # exposure, so the count of rows follows the exposure basis rather than the
  # width of the design. The exposure entering through two of them is what
  # settles the reading, so a covariate basis beside it changes nothing about
  # that either.
  expect_identical(nrow(est), 2L)
  expect_identical(est$contrast, c("poly(e, 2)1", "poly(e, 2)2"))
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[c("poly(e, 2)1", "poly(e, 2)2")]),
    tolerance = 1e-10
  )
  expect_dose_basis_conditional(res, fx$outcome_mod)
  expect_dose_basis_refuses_marginal(res)
  expect_identical(as.integer(res$fit@n_params), 3L + 1L + 2L + 6L)
})

test_that("a basis term mixing the exposure and a covariate is refused", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()

  # A basis expanded against a covariate contributes coefficients that are
  # changes in the dose response per unit of that covariate, so there is no one
  # effect for a row to report and no value of the covariate a label could name.
  # The refusal is the one every mixed term gets.
  interacted <- fit_dose_basis(dat, "poly(e, 2):x1")
  expect_error(
    ipw(interacted$ps_mod, interacted$outcome_mod, .data = dat),
    class = "propensity_ipw_msm_error"
  )

  crossed <- fit_dose_basis(dat, "poly(e, 2) * x1")
  expect_error(
    ipw(crossed$ps_mod, crossed$outcome_mod, .data = dat),
    class = "propensity_ipw_msm_error"
  )
})

test_that("the mixed-term refusal names the term and the way out", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2):x1")

  expect_propensity_error(ipw(fx$ps_mod, fx$outcome_mod, .data = dat))
})
