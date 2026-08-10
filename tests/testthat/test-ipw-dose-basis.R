# Basis-matrix dose-response marginal structural models on the plain continuous
# path: an ordinary single-dose fit with stabilized continuous ate weights and
# an lm or glm marginal structural model whose exposure term is a basis such as
# `poly(e, 2)`, `splines::ns(e, 3)`, or `splines::bs(e, 3)`. No joint container,
# no `.by`; this is the route a dose-response curve has always taken, asked for
# a curve written as a basis rather than as `e + I(e^2)`.
#
# Three things follow from the boundary being variable membership. A basis term
# reads the exposure alone, so it is admitted however many columns it expands
# to; each of those columns carries a coefficient, so each is a row, named by
# the coefficient rather than by the term, since the two differ exactly here. A
# basis term reading a covariate reads a covariate, so it is admitted as one and
# contributes no row. A term reading both, such as `poly(e, 2):x1`, keeps the
# refusal every mixed term gets: its coefficient is a change in the effect per
# unit of the covariate, and no row could name the covariate value it holds at.
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

# The column contract of a multi-row coefficient surface. A basis puts a
# contrast column on the table naming the coefficient each row reports, placed
# where every other route places it, and puts no group column there: nothing
# about a dose-response curve is evaluated in a stratum.
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

# The accessors a coefficient-shaped surface keeps in lockstep with the
# estimates frame. Written once because the claim is the same for every basis:
# the row labels are unique, every accessor keys by them in the same order, the
# reported standard errors are the diagonal of the reported covariance, and that
# covariance is a real one rather than a diagonal assembled row by row.
#
# `effect` is pinned here too, since it is one word repeated down the table and
# the same word for every basis. A surface whose rows are named by their
# coefficients reports `coef` at an identity link rather than `slope`: a basis
# coefficient is not a slope, and `slope` is reserved for the surfaces where a
# row really is one and says where it is evaluated.
expect_dose_basis_accessors <- function(res, effect = "coef") {
  estimates <- res$estimates
  labels <- paste(estimates$effect, estimates$contrast)

  testthat::expect_identical(estimates$effect, rep(effect, nrow(estimates)))
  testthat::expect_true(all(is.finite(estimates$std.err)))
  testthat::expect_true(all(estimates$std.err > 0))
  testthat::expect_true(all(estimates$ci.lower < estimates$ci.upper))
  testthat::expect_identical(anyDuplicated(labels), 0L)
  testthat::expect_identical(names(coef(res)), labels)
  testthat::expect_identical(dimnames(vcov(res)), list(labels, labels))
  testthat::expect_identical(rownames(confint(res)), labels)
  testthat::expect_equal(
    sqrt(diag(vcov(res))),
    estimates$std.err,
    ignore_attr = TRUE
  )

  # The basis coefficients come out of one fit of one model, so they covary; an
  # assembly that estimated each row on its own would report zeros here.
  covariance <- vcov(res)
  entry <- covariance[1L, 2L]
  testthat::expect_true(is.finite(entry))
  testthat::expect_gt(abs(entry), 1e-8)
  testthat::expect_equal(covariance, t(covariance), tolerance = 1e-12)

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
  expect_dose_basis_accessors(res)

  # The stacked system for a stabilized single-dose fit: three propensity score
  # coefficients, the conditional variance the density ratio needs, the two
  # moments the stabilizing numerator needs, and the outcome coefficients. Only
  # the last block moves, and it moves by the basis dimension.
  expect_identical(as.integer(res$fit@n_params), 3L + 1L + 2L + 3L)

  bare <- fit_dose_basis(dat, "e")
  bare_res <- ipw(bare$ps_mod, bare$outcome_mod)
  expect_identical(as.integer(bare_res$fit@n_params), 3L + 1L + 2L + 2L)
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
  expect_dose_basis_accessors(res)

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
  expect_dose_basis_accessors(res)

  # A cubic B-spline with no interior knots spans the same three columns a
  # natural spline of the same argument does, so the system is the same width.
  expect_identical(as.integer(res$fit@n_params), 3L + 1L + 2L + 4L)
})

# ---- the closed form ---------------------------------------------------------

test_that("the reported coefficients are the weighted marginal structural model's", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)

  # An identity-link marginal structural model is the estimator: nothing is
  # standardized, so each row is exactly a coefficient of the weighted fit and
  # the comparison is a closed form rather than an agreement to within a
  # g-computation.
  beta <- coef(fx$outcome_mod)
  expect_equal(
    res$estimates$estimate,
    unname(beta[c("poly(e, 2)1", "poly(e, 2)2")]),
    tolerance = 1e-10
  )

  # Row by row as well, so a reordering of the surface fails here rather than
  # passing on a vector that happens to line up.
  pick <- function(contrast) {
    res$estimates$estimate[res$estimates$contrast == contrast]
  }
  expect_equal(pick("poly(e, 2)1"), unname(beta[["poly(e, 2)1"]]))
  expect_equal(pick("poly(e, 2)2"), unname(beta[["poly(e, 2)2"]]))

  # The intercept is not a causal coefficient and is not reported, so the
  # surface is one row shorter than the coefficient vector.
  expect_identical(nrow(res$estimates), length(beta) - 1L)
})

test_that("a logit marginal structural model reports the same rows on the odds scale", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "poly(e, 2)", outcome_family = "binomial")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  est <- res$estimates

  # The scale follows the outcome link, which is what the continuous path
  # already does with a single-term logit dose response. The rows are otherwise
  # the ones the gaussian fit reports, and the contrast still names the
  # coefficient.
  expect_dose_basis_columns(est)
  expect_identical(nrow(est), 2L)
  expect_identical(est$effect, c("log(or)", "log(or)"))
  expect_identical(est$contrast, c("poly(e, 2)1", "poly(e, 2)2"))
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[c("poly(e, 2)1", "poly(e, 2)2")]),
    tolerance = 1e-10
  )
  expect_dose_basis_accessors(res, effect = "log(or)")
  expect_identical(as.integer(res$fit@n_params), 3L + 1L + 2L + 3L)
})

test_that("the conditional reading of a basis fit is the outcome model's", {
  skip_if_not_installed("deli")
  dat <- sim_dose_basis()
  fx <- fit_dose_basis(dat, "splines::ns(e, 3)")
  conditional <- causalgenerics::as_conditional(
    ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  )

  # Read through the tidier, which is the surface the reading presents; the
  # stored estimates frame is the marginal one whichever reading is recorded.
  # The conditional reading is the whole coefficient vector, so the intercept
  # the marginal surface drops is back.
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

  expect_propensity_error(ipw(fx$ps_mod, fx$outcome_mod))
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
  expect_dose_basis_accessors(res)
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
  # reserved for exactly the rows it is true of.
  fx <- fit_dose_basis(dat, "e + splines::ns(x1, 3)")
  res <- ipw(fx$ps_mod, fx$outcome_mod, .data = dat)
  est <- res$estimates

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
  # width of the design.
  expect_identical(nrow(est), 2L)
  expect_identical(est$contrast, c("poly(e, 2)1", "poly(e, 2)2"))
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[c("poly(e, 2)1", "poly(e, 2)2")]),
    tolerance = 1e-10
  )
  expect_dose_basis_accessors(res)
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
