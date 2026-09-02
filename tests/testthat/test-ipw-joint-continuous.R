# A joint intervention with a continuous component, and the marginal structural
# models the single-exposure-term validator now admits.
#
# Two things join here. The first is the two-model route with a continuous
# second treatment: `joint_wt_models(a = , e = )` where `e` is a dose. A dose has
# no cells, so there is no crossing to report counterfactual risks over and the
# surface is coefficient-shaped rather than cell-shaped. It reports the marginal
# structural model's own causal coefficients, which for an identity-link model
# are exactly the weighted fit's coefficients.
#
# The second is the relaxation that makes such a model fittable at all. The
# continuous path required a marginal structural model with exactly one exposure
# term, which refuses a dose-response curve as firmly as it refuses a model with
# two exposures in it. The requirement narrows rather than disappears: what is
# refused is a term reading something other than the exposure.
#
# ---- the vocabulary surface -------------------------------------------------
#
# For `y ~ a * e` with `a` binary and `e` a dose, three coefficients carry causal
# content and each is reported as a row:
#
#   effect   contrast          group            coefficient
#   diff     a: 1 vs 0         e = 0            a
#   slope    e: per unit       a = 0            e
#   diff     a: 1 vs 0         e + 1 vs e       a:e
#
# `effect` names the scale, keeping the vocabulary each treatment type already
# uses: a level contrast of the binary treatment is a `diff` under an identity
# link, a per-unit change in the dose is a `slope`, and both become `log(or)`
# under a logit link, which is what the continuous path already reports there.
# `contrast` names the treatment being varied and how, and `group` says where in
# the other treatment's range the row is evaluated. The interaction row keeps the
# discrete route's idiom, where the group is a comparison of the other
# treatment's levels; with a dose the comparison is a one-unit step, which for a
# linear model is the whole of it.
#
# ---- the coefficient surface ------------------------------------------------
#
# Those three readings hold of a model that is linear in each treatment and of
# no other, so a marginal structural model carrying a transformed or a basis
# treatment term reports a different surface rather than the same one with
# different numbers in it. The coefficient surface reports one row per
# treatment-reading coefficient, `contrast` carries the coefficient name exactly
# as the fit writes it, and there is no `group` column at all: the surface makes
# no claim about where a row is evaluated, because for a curve there is no one
# place. An interpretable dose response is built downstream from `coef()` and
# `vcov()`, which is what the covariance between the rows is for.
#
# `effect` names the scale and only the scale. A coefficient of an identity-link
# model is `coef`, a new word rather than `diff` or `slope`, since both of those
# carry an evaluated-at claim on the vocabulary surface and reusing one here
# would make the same word mean two things. Under a logit or a log link the
# existing words are already honest, since a coefficient there is a log odds
# ratio or a log risk ratio per unit of whatever column it multiplies, so
# `log(or)` and `log(rr)` are kept.
#
# Which surface a fit reports is decided by the model alone. A bare-term model
# reports the vocabulary surface, and every other treatment-reading model
# reports the coefficient surface. Two things are unchanged by this. A term
# reading a treatment and a covariate together is still refused, since its
# coefficient depends on the covariate and naming it after the column would
# report a quantity nothing pins down. And a bare-term model whose columns are
# coded some other way keeps its refusal rather than falling through to the
# coefficient surface: such a model reads as linear, and answering it with a
# differently shaped table would be a worse outcome than the guided refusal.

# ---- data simulator ---------------------------------------------------------

# A binary treatment and a dose that depends on it, sharing confounders. The
# outcome carries a genuine interaction, so the third row is not zero.
sim_joint_continuous <- function(seed = 7501, n = 800) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * x2))
  e <- 0.4 + 0.5 * x1 - 0.3 * x2 - 0.8 * a + rnorm(n)
  y <- 1 + 0.7 * a + 0.5 * e + 0.6 * a * e + rnorm(n)
  yb <- rbinom(n, 1, plogis(-0.4 + 0.6 * a + 0.4 * e + 0.5 * a * e))
  data.frame(x1, x2, a, e, y, yb)
}

# The same crossing with a dose that is positive everywhere, which a dose model
# with a log link needs: the link takes the log of the response to start its
# iteration, and a dose that can be zero or negative has no such start. The dose
# is the one above shifted up rather than exponentiated, so it stays
# conditionally normal with a single spread.
sim_joint_continuous_positive <- function(seed = 7501, n = 800) {
  dat <- sim_joint_continuous(seed = seed, n = n)
  dat$e <- dat$e + 5
  withr::local_seed(seed + 1L)
  dat$y <- 1 + 0.7 * dat$a + 0.5 * dat$e + 0.6 * dat$a * dat$e + rnorm(n)
  # yb came from the pre-shift dose and no caller reads it.
  dat$yb <- NULL
  dat
}

# The same crossing with a contaminated dose, which is what separates a robust
# dose model from the least-squares one fit to the same columns. Without the
# contamination the two agree and a system that stacked the wrong score would
# report the right coefficients anyway.
sim_joint_continuous_outliers <- function(
  seed = 7501,
  n = 800,
  n_outliers = 40
) {
  dat <- sim_joint_continuous(seed = seed, n = n)
  dat$e[seq_len(n_outliers)] <- dat$e[seq_len(n_outliers)] + 12
  withr::local_seed(seed + 2L)
  dat$y <- 1 + 0.7 * dat$a + 0.5 * dat$e + 0.6 * dat$a * dat$e + rnorm(n)
  dat$yb <- NULL
  dat
}

# A single dose with a curved response, for the validator relaxation.
sim_dose_curve <- function(seed = 7502, n = 800) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  e <- 0.4 + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  a <- rbinom(n, 1, plogis(0.2 * x1))
  y <- 1 + 0.5 * e + 0.2 * e^2 + 0.4 * x1 + rnorm(n)
  data.frame(x1, x2, a, e, y)
}

quiet_wt <- function(expr) {
  withr::with_options(list(propensity.quiet = TRUE), expr)
}

# ---- model fitting ----------------------------------------------------------

# The dose models the registry reads a conditional density from. A log-link fit
# is iterated to a far tighter tolerance than the default, so its coefficients
# are the root of the score stacked against them rather than the point the
# iteration happened to stop at; a robust fit is tightened for the same reason.
# The additive fit carries a smooth of the covariate, which is what makes it a
# `gam` rather than a `glm` written in another spelling.
fit_joint_dose <- function(dose_type, e_rhs, dat) {
  fmla <- stats::reformulate(e_rhs, response = "e")

  switch(
    dose_type,
    lm = lm(fmla, data = dat),
    glm_log = glm(
      fmla,
      data = dat,
      family = gaussian(link = "log"),
      control = glm.control(epsilon = 1e-14, maxit = 200)
    ),
    # `k` moves the psi's own threshold, which is where a caller writes one and
    # is not the `k2` that tunes the scale estimator. A fit at the default
    # threshold would be stacked correctly by a block that read either, so the
    # fixture is fit away from it.
    rlm = MASS::rlm(fmla, data = dat, acc = 1e-10, k = 2),
    gam = mgcv::gam(
      stats::reformulate(c("a", "s(x1)", "x2"), response = "e"),
      data = dat
    )
  )
}

# The two-model route with a continuous second treatment. The dose model
# conditions on the binary treatment, which `joint_wt_models()` requires of a
# continuous second component exactly as it requires it of a discrete one, and
# the dose weights are stabilized, which `wt_joint()` requires. Formulas are
# built here so the weights resolve when `lm()` and `glm()` look for them.
fit_joint_continuous <- function(
  dat,
  a_rhs = c("x1", "x2"),
  e_rhs = c("a", "x1", "x2"),
  outcome_rhs = "a * e",
  outcome_family = "gaussian",
  dose_score = NULL,
  dose_stabilize = TRUE,
  a_stabilize = NULL,
  a_score = NULL,
  dose_type = c("lm", "glm_log", "rlm", "gam"),
  .density = "normal",
  numerator = "marginal",
  .sigma = NULL
) {
  dose_type <- match.arg(dose_type)
  ps_a <- glm(
    stats::reformulate(a_rhs, response = "a"),
    data = dat,
    family = binomial()
  )
  ps_e <- fit_joint_dose(dose_type, e_rhs, dat)

  wts <- quiet_wt(wt_joint(
    # `NULL` leaves a binary component unstabilized, which is the default a
    # binary exposure resolves to, and a fitted model stabilizes it on the
    # probability of the level each unit took.
    wt_ate(ps_a, stabilize = a_stabilize, stabilization_score = a_score),
    wt_ate(
      as.double(fitted(ps_e)),
      dat$e,
      exposure_type = "continuous",
      # `TRUE` for the default numerator, and a fitted model for the one the
      # single-dose route estimates in its stacked system.
      stabilize = dose_stabilize,
      # A numerator the caller computed rather than one the estimator can
      # rebuild, which `dose_score` is here to supply.
      stabilization_score = dose_score,
      # The three things a dose's weights record about the ratio they are, each
      # of which the stacked system has to read back off the product.
      .density = .density,
      numerator = numerator,
      .sigma = .sigma
    ),
    exposure_type = c("binary", "continuous")
  ))

  outcome_var <- if (outcome_family == "binomial") "yb" else "y"
  fmla <- stats::reformulate(outcome_rhs, response = outcome_var)
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

  list(
    models = joint_wt_models(a = ps_a, e = ps_e),
    ps_a = ps_a,
    ps_e = ps_e,
    outcome_mod = outcome_mod,
    wts = wts
  )
}

# The same fixture with the binary treatment recoded as a two-level factor,
# which is the same intervention written a different way. Every model is refit
# from the recoded column, so the coding the outcome design carries is the one
# `ordered` and the session's contrast option imply, which is what the estimator
# has to read off the columns rather than off the formula.
fit_joint_factor <- function(dat, ordered = FALSE) {
  dat$a <- factor(
    ifelse(dat$a == 1, "yes", "no"),
    levels = c("no", "yes"),
    ordered = ordered
  )

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- lm(e ~ a + x1 + x2, data = dat)
  wts <- quiet_wt(wt_joint(
    wt_ate(ps_a),
    wt_ate(
      as.double(fitted(ps_e)),
      dat$e,
      exposure_type = "continuous",
      stabilize = TRUE
    ),
    exposure_type = c("binary", "continuous")
  ))

  list(
    models = joint_wt_models(a = ps_a, e = ps_e),
    outcome_mod = lm(y ~ a * e, data = dat, weights = wts),
    wts = wts
  )
}

# The single-dose route, with whatever marginal structural model a test wants to
# hand it. The weights are the stabilized continuous ate weights the path
# already documents.
fit_dose_msm <- function(dat, msm_rhs) {
  ps_mod <- lm(e ~ x1 + x2, data = dat)
  wts <- quiet_wt(wt_ate(
    as.double(fitted(ps_mod)),
    dat$e,
    exposure_type = "continuous",
    stabilize = TRUE
  ))
  outcome_mod <- lm(
    stats::reformulate(msm_rhs, response = "y"),
    data = dat,
    weights = wts
  )

  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# ---- the surface a continuous component reports -----------------------------

# The identity columns, written out against the conventions in the header.
joint_continuous_expected_rows <- function(scale = c("diff", "slope", "diff")) {
  data.frame(
    effect = scale,
    contrast = c("a: 1 vs 0", "e: per unit", "a: 1 vs 0"),
    group = c("e = 0", "a = 0", "e + 1 vs e"),
    stringsAsFactors = FALSE
  )
}

joint_continuous_labels <- function(estimates) {
  paste(estimates$effect, estimates$contrast, estimates$group)
}

# A coefficient row has no group to fold in, so its label is the two columns
# that name it.
joint_coefficient_labels <- function(estimates) {
  paste(estimates$effect, estimates$contrast)
}

# What every coefficient-surface fit holds, written once because it is the same
# claim each time: the column layout with no group in it, the identity columns,
# unique labels that every accessor keys by in the same order, standard errors
# that are the diagonal of the reported covariance, and a covariance that is a
# real one rather than a diagonal assembled a row at a time.
expect_joint_coefficient_surface <- function(res, contrast, effect = "coef") {
  estimates <- res$estimates
  labels <- joint_coefficient_labels(estimates)

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
  testthat::expect_identical(estimates$effect, rep(effect, length(contrast)))
  testthat::expect_identical(estimates$contrast, contrast)

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

  # Every row is a parameter of one system solved on one sample, so the rows
  # covary and a downstream dose curve has a covariance to build from.
  covariance <- vcov(res)
  testthat::expect_true(is.finite(covariance[1L, 2L]))
  testthat::expect_gt(abs(covariance[1L, 2L]), 1e-10)
  testthat::expect_equal(covariance, t(covariance), tolerance = 1e-12)

  invisible(res)
}

# The reported estimates are the weighted fit's own coefficients, keyed by the
# name each row carries rather than by position, so a reordering of the surface
# fails here rather than passing on a vector that happens to line up.
expect_joint_coefficient_estimates <- function(res, outcome_mod) {
  beta <- coef(outcome_mod)
  estimates <- res$estimates

  testthat::expect_equal(
    estimates$estimate,
    unname(beta[estimates$contrast]),
    tolerance = 1e-9
  )

  for (i in seq_len(nrow(estimates))) {
    testthat::expect_equal(
      estimates$estimate[[i]],
      unname(beta[[estimates$contrast[[i]]]]),
      tolerance = 1e-9
    )
  }

  invisible(res)
}

# ---- baselines --------------------------------------------------------------

test_that("the continuous fixture builds the container and the product weights", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)

  # The container records a binary treatment and a dose, and it accepted the
  # dose model because that model conditions on the binary treatment: the
  # factorization guardrail is about the factorization, not the response type.
  expect_true(is_joint_wt_models(fx$models))
  expect_identical(fx$models$names, c("a", "e"))
  expect_identical(
    fx$models$exposure_type,
    c(a = "binary", e = "continuous")
  )
  expect_true("a" %in% all.vars(formula(fx$ps_e)[[3]]))

  # The product is a joint weight whose continuous component is stabilized,
  # which `wt_joint()` requires and which the estimator has to rebuild.
  expect_true(is_joint_wt(fx$wts))
  expect_identical(estimand(fx$wts), "ate")
  expect_identical(
    joint_wt_meta(fx$wts)$exposure_type,
    c("binary", "continuous")
  )
  expect_identical(joint_wt_meta(fx$wts)$stabilized, c(FALSE, TRUE))
  expect_true(is_stabilized(fx$wts))
})

test_that("a single-term dose model still reports one slope with no contrast column", {
  dat <- sim_dose_curve()
  fx <- fit_dose_msm(dat, "e")
  res <- ipw(fx$ps_mod, fx$outcome_mod)

  # The shape the relaxation must leave alone. One exposure coefficient, one
  # row, and no contrast column: a column repeating one value down a
  # single-row table would read as a contrast that was named.
  expect_identical(nrow(res$estimates), 1L)
  expect_identical(res$estimates$effect, "slope")
  expect_null(res$estimates[["contrast"]])
  expect_length(names(res$estimates), 8L)
  expect_equal(
    res$estimates$estimate,
    unname(coef(fx$outcome_mod)[["e"]]),
    tolerance = 1e-8
  )

  # The stacked system for a stabilized single-dose fit: three propensity score
  # coefficients, the conditional variance the density ratio needs, the two
  # moments the stabilizing numerator needs, and the two outcome coefficients.
  expect_identical(as.integer(res$fit@n_params), 8L)
  expect_true("sigma2_d" %in% names(res$fit@theta))
})

test_that("a marginal structural model reading a covariate is still refused", {
  dat <- sim_dose_curve()

  # The relaxation narrows the requirement to terms that read the exposure and
  # nothing else, so a term reading a covariate keeps the refusal it has today:
  # the effect such a model reports depends on the covariate and there is no one
  # coefficient to name.
  covariate <- fit_dose_msm(dat, "e * x1")
  expect_error(
    ipw(covariate$ps_mod, covariate$outcome_mod),
    class = "propensity_ipw_msm_error"
  )
  expect_propensity_error(ipw(covariate$ps_mod, covariate$outcome_mod))

  # The case the guardrail exists for: two undeclared exposures handed to the
  # plain continuous route with one treatment's weights. Nothing declares a
  # crossing here, so the model has two exposures and no reported effect.
  undeclared <- fit_dose_msm(dat, "a * e")
  expect_error(
    ipw(undeclared$ps_mod, undeclared$outcome_mod),
    class = "propensity_ipw_msm_error"
  )
  expect_propensity_error(ipw(undeclared$ps_mod, undeclared$outcome_mod))
})

# ---- part 2: the validator relaxation ---------------------------------------

test_that("a polynomial dose model reports one row per dose coefficient", {
  dat <- sim_dose_curve()
  fx <- fit_dose_msm(dat, "e + I(e^2)")
  res <- ipw(fx$ps_mod, fx$outcome_mod)
  est <- res$estimates

  # More than one exposure coefficient, so the frame the stacked system
  # estimates gains a contrast column naming the term each one belongs to,
  # placed where every other route places it, and the scale word becomes `coef`.
  # Neither of these two coefficients is the slope of the dose response, since
  # the response has a different slope at every dose, so the surface says what
  # the numbers are rather than naming them after a quantity only one of them
  # could be. That is also why the result reports the conditional reading and
  # refuses the marginal one: there is no row of it that answers the question
  # the marginal reading is asked.
  expect_identical(
    names(est),
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
  expect_identical(nrow(est), 2L)
  expect_identical(est$effect, c("coef", "coef"))
  expect_identical(est$contrast, c("e", "I(e^2)"))
  expect_null(est[["group"]])

  # The reported coefficients are the weighted fit's coefficients: the marginal
  # structural model is the estimator, and the stacked system reports the
  # dose block of it rather than standardizing anything.
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[c("e", "I(e^2)")]),
    tolerance = 1e-8
  )

  expect_true(all(is.finite(est$std.err)))
  expect_true(all(est$std.err > 0))
  expect_true(all(est$ci.lower < est$ci.upper))

  # What the result reports is the conditional reading: the outcome model's own
  # coefficient vector, intercept included, under the block the stacked system
  # leaves for it.
  expect_identical(res$effects, "conditional")
  expect_identical(res$readings, "conditional")
  expect_identical(class(res), "ipw")
  expect_identical(coef(res), coef(fx$outcome_mod))
  expect_identical(names(coef(res)), c("(Intercept)", "e", "I(e^2)"))
  expect_identical(vcov(res), vcov(res$outcome_mod))

  # The frame the stacked system estimates carries the standard error of each
  # dose coefficient, which is that coefficient's diagonal entry in the same
  # block. Nothing else pins those numbers to a covariance, and `pool_ipw()`
  # reads them off the frame rather than through the accessors.
  expect_equal(
    unname(sqrt(diag(vcov(res)))[est$contrast]),
    est$std.err,
    tolerance = 1e-12
  )
  expect_error(
    causalgenerics::as_marginal(res),
    class = "causalgenerics_unsupported_reading_marginal"
  )

  # One outcome coefficient more than the single-term fit, and nothing else
  # about the system changes.
  expect_identical(as.integer(res$fit@n_params), 9L)
})

test_that("the relaxation admits any term that reads the exposure alone", {
  dat <- sim_dose_curve()

  # The boundary is variable membership, not the shape of the function: what
  # made a multi-term model unreportable was a term reading something other than
  # the exposure, and a curve in the exposure alone is as reportable as a
  # quadratic in it. Restricting to polynomials would also mean deciding by
  # inspection which expressions count as one.
  fx <- fit_dose_msm(dat, "e + sin(e)")
  res <- ipw(fx$ps_mod, fx$outcome_mod)

  expect_identical(nrow(res$estimates), 2L)
  expect_identical(res$estimates$effect, c("coef", "coef"))
  expect_identical(res$estimates$contrast, c("e", "sin(e)"))
  expect_identical(res$effects, "conditional")
  expect_equal(
    res$estimates$estimate,
    unname(coef(fx$outcome_mod)[c("e", "sin(e)")]),
    tolerance = 1e-8
  )
})

# ---- part 1: the joint route with a dose ------------------------------------

test_that("ipw() over a binary and a continuous treatment reports the vocabulary surface", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates

  expect_s3_class(res, "ipw")
  expect_identical(res$estimand, "ate")
  expect_identical(res$se_method, "mestimation")

  # Coefficient-shaped, not cell-shaped: a dose has no cells, so there are no
  # counterfactual risks to report and no `mean` rows.
  #
  # A bare-term model keeps this surface exactly as it is. The coefficient
  # surface exists for the models this vocabulary cannot describe, and a model
  # it can describe never reaches it, so these rows are a regression pin: the
  # words, the group claims, and the column layout all stay.
  expect_false(any(est$effect == "mean"))
  expected <- joint_continuous_expected_rows()
  expect_identical(nrow(est), 3L)
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)
  expect_identical(
    names(est),
    c(
      "effect",
      "contrast",
      "group",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
})

test_that("the reported coefficients are the weighted marginal structural model's", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)
  res <- ipw(fx$models, fx$outcome_mod)

  # The closed form. An identity-link marginal structural model is the
  # estimator: the binary treatment's effect at a dose of zero is its
  # coefficient, the dose slope among the untreated is its coefficient, and the
  # change in the binary treatment's effect per unit of dose is the interaction
  # coefficient. Nothing is standardized, so this is exact rather than to within
  # a g-computation.
  beta <- coef(fx$outcome_mod)
  expect_equal(
    res$estimates$estimate,
    unname(beta[c("a", "e", "a:e")]),
    tolerance = 1e-8
  )

  # Row by row, so a reordering of the surface fails here rather than passing on
  # a vector that happens to match.
  pick <- function(effect, contrast, group) {
    res$estimates$estimate[
      res$estimates$effect == effect &
        res$estimates$contrast == contrast &
        res$estimates$group == group
    ]
  }
  expect_equal(pick("diff", "a: 1 vs 0", "e = 0"), unname(beta[["a"]]))
  expect_equal(pick("slope", "e: per unit", "a = 0"), unname(beta[["e"]]))
  expect_equal(
    pick("diff", "a: 1 vs 0", "e + 1 vs e"),
    unname(beta[["a:e"]])
  )
})

test_that("a logit marginal structural model reports the same rows on the odds scale", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, outcome_family = "binomial")
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates

  # The scale follows the outcome link, which is what the continuous path
  # already does when it labels a logit dose response `log(or)`. The contrast
  # and the group say which coefficient a row is whatever the scale.
  expected <- joint_continuous_expected_rows(scale = rep("log(or)", 3))
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)

  beta <- coef(fx$outcome_mod)
  expect_equal(
    est$estimate,
    unname(beta[c("a", "e", "a:e")]),
    tolerance = 1e-8
  )
})

test_that("an additive marginal structural model reports two ungrouped rows", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, outcome_rhs = "a + e")
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates

  # With no interaction term the model forces one effect for the binary
  # treatment at every dose and one slope at either level of it, so neither row
  # is evaluated anywhere in particular and there is no interaction to report.
  # Bare terms, so the vocabulary surface, and these rows are pinned as they
  # are: `overall` is a group claim and not the absence of one, which is what
  # distinguishes this table from a coefficient-surface table of the same width.
  expect_identical(nrow(est), 2L)
  expect_identical(est$effect, c("diff", "slope"))
  expect_identical(est$contrast, c("a: 1 vs 0", "e: per unit"))
  expect_identical(est$group, c("overall", "overall"))
  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[c("a", "e")]),
    tolerance = 1e-8
  )
})

test_that("a covariate in a bare-term model contributes no vocabulary row", {
  dat <- sim_joint_continuous()

  # Adjusting the marginal structural model for a covariate is a different thing
  # from making the covariate a treatment. `x1` reads neither treatment, so it
  # contributes a coefficient the fit needs and the surface has nothing to say
  # about: the rows are the three the crossing reports and no more, and each of
  # them is still the coefficient it was without the covariate in the model.
  #
  # The design puts `x1` between the dose column and the interaction column, so
  # a surface that reported whatever was not the intercept would show it here,
  # in the middle of the table, rather than tidily at the end.
  fx <- fit_joint_continuous(dat, outcome_rhs = "a * e + x1")
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates

  expected <- joint_continuous_expected_rows()
  expect_identical(nrow(est), 3L)
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)
  expect_false("x1" %in% est$contrast)

  expect_equal(
    est$estimate,
    unname(coef(fx$outcome_mod)[c("a", "e", "a:e")]),
    tolerance = 1e-8
  )

  labels <- joint_continuous_labels(est)
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(names(coef(res)), labels)

  # The covariate's coefficient is inside the outcome block of the stacked
  # system even though no row reports it: five outcome coefficients against the
  # four a model without it has, and everything else unchanged.
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 5L)
})

# ---- standard errors and the stacked system ---------------------------------

test_that("a joint continuous fit reports a usable standard error for every row", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates

  expect_identical(nrow(est), 3L)
  expect_true(all(is.finite(est$std.err)))
  expect_true(all(est$std.err > 0))
  expect_true(all(est$ci.lower < est$ci.upper))
  expect_equal(sqrt(diag(vcov(res))), est$std.err, ignore_attr = TRUE)

  labels <- joint_continuous_labels(est)
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(names(coef(res)), labels)
  expect_identical(dimnames(vcov(res)), list(labels, labels))
  expect_identical(rownames(confint(res)), labels)

  # The three reported coefficients come out of one fit of one model, so they
  # covary; a stitched assembly would report zeros here.
  covariance <- vcov(res)
  entry <- covariance[labels[[1]], labels[[3]]]
  expect_true(is.finite(entry))
  expect_gt(abs(entry), 1e-8)
  expect_equal(covariance, t(covariance), tolerance = 1e-12)
})

test_that("the stacked system carries the dose model's score and its variance", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)
  res <- ipw(fx$models, fx$outcome_mod)

  # The dimension of the system, spelled out so a mismatch says which block is
  # wrong. Three coefficients for the binomial treatment model, whose score
  # carries no variance parameter of its own; four for the dose model plus the
  # conditional variance its density ratio needs; the two marginal moments the
  # stabilizing numerator needs, which the continuous component is required to
  # carry; and four outcome coefficients. The surface reports the outcome block
  # directly, so there are no marginal means and no contrast parameters.
  #
  # A gaussian outcome contributes no residual variance here: the weighted
  # least-squares score is the whole of what the outcome block solves, which is
  # what the single-dose route already does.
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 4L)
  expect_identical(nobs(res), 800L)
  expect_identical(df.residual(res), 800L - 14L)

  theta <- names(res$fit@theta)
  expect_true(any(grepl("sigma2", theta, fixed = TRUE)))
  expect_true(all(c("x1", "x2") %in% theta))

  # The dose model's dimension is what moves when the dose model changes, so the
  # count tracks the pair rather than carrying a fixed allowance for it.
  smaller <- fit_joint_continuous(dat, e_rhs = c("a", "x1"))
  res_smaller <- ipw(smaller$models, smaller$outcome_mod)
  expect_identical(as.integer(res_smaller$fit@n_params), 13L)
})

test_that("the estimator rebuilds the product weights the container implies", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)

  # The weights that fit the outcome model are the product of a binary inverse
  # probability weight and a stabilized density ratio, and the stacked system
  # rebuilds both from the propensity score blocks of theta on every evaluation.
  # A preflight that rebuilt them differently, or a stabilizing numerator the
  # system did not carry, would put the seed off the root and the fit would
  # either refuse or drift off the closed form above. That the fit returns and
  # the estimates equal the weighted coefficients exactly is what pins the
  # rebuild.
  res <- ipw(fx$models, fx$outcome_mod)
  expect_s3_class(res, "ipw")
  expect_equal(
    res$estimates$estimate,
    unname(coef(fx$outcome_mod)[c("a", "e", "a:e")]),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("the conditional reading of a joint continuous fit is the outcome model's", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)
  conditional <- causalgenerics::as_conditional(ipw(fx$models, fx$outcome_mod))

  # Read through the tidier, which is the surface the reading presents; the
  # stored estimates frame is the marginal one whichever reading is recorded.
  tidied <- tidy(conditional)
  expect_identical(tidied$term, names(coef(fx$outcome_mod)))
  expect_equal(
    tidied$estimate,
    unname(coef(fx$outcome_mod)),
    tolerance = 1e-8
  )
  expect_false("group" %in% names(tidied))
  expect_false("contrast" %in% names(tidied))
})

# ---- the dose models the registry reads -------------------------------------

# The stacked system's seed is the fits' own parameters, so every estimating
# equation in it has to be at its root there. A block written against the wrong
# score sits away from zero here, which is the same fault the estimates would
# show later as a propensity score model the user never fit.
expect_joint_root_seeded <- function(spec) {
  layout <- ipw_theta_layout(spec)
  mat <- build_ipw_psi(spec, layout)(layout$init)

  testthat::expect_true(is.matrix(mat))
  testthat::expect_identical(dim(mat), c(length(layout$init), spec$n))
  testthat::expect_false(anyNA(mat))
  testthat::expect_true(all(abs(rowSums(mat)) / spec$n < 1e-8))

  invisible(layout)
}

# The weights the system rebuilds from its seed against the ones the outcome
# model was fit with. The preflight compares the two at a relative tolerance of
# 1e-6 before solving anything; they agree far more closely than that when the
# dose's ratio is rebuilt from what the product records, so this holds the
# rebuild to the arithmetic rather than to the guard.
expect_joint_weights_at_init <- function(spec, wts) {
  layout <- ipw_theta_layout(spec)

  testthat::expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(wts),
    tolerance = 1e-12
  )

  invisible(layout)
}

# The message an expression's error carries, collapsed onto one line so that a
# refusal's wrapped text can be matched as it was written, and empty where the
# expression raised no error at all: what a message says is worth checking only
# once the class it came with has been.
joint_error_message <- function(expr) {
  cnd <- tryCatch(expr, error = function(e) e)

  if (!inherits(cnd, "condition")) {
    return("")
  }

  gsub("[[:space:]]+", " ", conditionMessage(cnd))
}

expect_joint_dose_estimates <- function(res, outcome_mod) {
  testthat::expect_equal(
    res$estimates$estimate,
    unname(coef(outcome_mod)[c("a", "e", "a:e")]),
    tolerance = 1e-8
  )
  testthat::expect_true(all(is.finite(res$estimates$std.err)))
  testthat::expect_true(all(res$estimates$std.err > 0))

  invisible(res)
}

test_that("a log-link dose model is stacked at its own score", {
  dat <- sim_joint_continuous_positive()
  fx <- fit_joint_continuous(dat, dose_type = "glm_log")

  # A gaussian dose model fit through a log link reads its conditional mean as
  # the exponential of the linear predictor, and the score its coefficients
  # solve is one the system can write, so the density the weights divide by is
  # a function of theta here as much as at the identity link. What made the
  # link unsupported was a block that was the least-squares one whatever the
  # link, which would have put the seed off the root.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  expect_joint_root_seeded(spec)
  expect_joint_weights_at_init(spec, fx$wts)

  res <- ipw(fx$models, fx$outcome_mod)
  expect_joint_dose_estimates(res, fx$outcome_mod)

  # The dimension is unchanged by the link: the dose block is its coefficients
  # and the conditional variance its density ratio divides by.
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 4L)
})

test_that("an integrated numerator on the dose stacks no marginal moments", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, numerator = "integrated")

  # The integrated numerator is the conditional density averaged over the
  # units, which is built from the dose block and the data alone, so there is
  # nothing left for a stabilization block to estimate. The product records
  # that the dose was stabilized either way, so the numerator it records is
  # what tells the two apart.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  layout <- expect_joint_weights_at_init(spec, fx$wts)
  expect_length(layout$idx$stab, 0L)
  expect_joint_root_seeded(spec)

  res <- ipw(fx$models, fx$outcome_mod)
  expect_joint_dose_estimates(res, fx$outcome_mod)
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 4L)
  expect_false(any(grepl("stab_", names(coef(res$fit)), fixed = TRUE)))
})

test_that("a heavier-tailed dose density is the one the system rebuilds", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, .density = dens_t(4))

  # Which family the ratio is read at is recorded per component on the product,
  # so the system rebuilds the dose weights as the t density they are rather
  # than as the normal the route assumed. A rebuild at the wrong family reaches
  # a different vector and the preflight refuses before anything is solved.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  expect_joint_weights_at_init(spec, fx$wts)
  expect_joint_root_seeded(spec)

  res <- ipw(fx$models, fx$outcome_mod)
  expect_joint_dose_estimates(res, fx$outcome_mod)
})

test_that("a dose scale fit by maximum likelihood is stacked at its own equation", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, .density = dens_t(4, sigma_method = "mle"))

  # The dose records which estimator its scale was read with, so the scale sits
  # in the dose block as the pooled spread does, and the row that estimates it
  # is the equation the scale the weights were built at is the root of. A block
  # stacked at the residual moment instead would be seeded away from its own
  # root and would rebuild a different vector of weights.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  layout <- expect_joint_weights_at_init(spec, fx$wts)
  expect_length(
    layout$idx$ps,
    length(coef(fx$ps_a)) + length(coef(fx$ps_e)) + 1L
  )
  expect_joint_root_seeded(spec)

  res <- ipw(fx$models, fx$outcome_mod)
  expect_joint_dose_estimates(res, fx$outcome_mod)
  expect_true("sigma2_e" %in% names(coef(res$fit)))
})

test_that("a fixed scalar sigma on the dose is a constant rather than a parameter", {
  dat <- sim_joint_continuous()

  # A spread the caller fixed is a known constant, so the dose block is its
  # coefficients alone and the system carries none of that number's
  # uncertainty. The value is away from the pooled residual standard deviation,
  # so weights rebuilt at the pooled spread would not match.
  fx <- fit_joint_continuous(dat, .sigma = 1.25)

  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  layout <- expect_joint_weights_at_init(spec, fx$wts)
  expect_length(
    layout$idx$ps,
    length(coef(fx$ps_a)) + length(coef(fx$ps_e))
  )
  expect_joint_root_seeded(spec)

  res <- ipw(fx$models, fx$outcome_mod)
  expect_joint_dose_estimates(res, fx$outcome_mod)
  expect_false("sigma2_e" %in% names(coef(res$fit)))
  expect_identical(as.integer(res$fit@n_params), 3L + 4L + 2L + 4L)
})

test_that("a robust dose model is stacked at the score it solves", {
  skip_if_not_installed("MASS")
  dat <- sim_joint_continuous_outliers()
  fx <- fit_joint_continuous(dat, dose_type = "rlm")

  # `rlm` minimizes a Huber loss rather than the sum of squares, so its
  # coefficients are the root of a different equation. The block has to be that
  # equation: seeded from the robust fit and solving the least-squares score,
  # the system drifts to the least-squares coefficients and reports a dose
  # model the user never fit.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  expect_joint_root_seeded(spec)
  expect_joint_weights_at_init(spec, fx$wts)

  res <- ipw(fx$models, fx$outcome_mod)
  theta <- coef(res$fit)
  dose_block <- theta[length(coef(fx$ps_a)) + seq_along(coef(fx$ps_e))]
  expect_equal(
    unname(dose_block),
    unname(coef(fx$ps_e)),
    tolerance = 1e-6
  )
  expect_joint_dose_estimates(res, fx$outcome_mod)

  # The contamination is what makes the distinction visible: the two fits
  # disagree in the first decimal place on this fixture, so a system that
  # solved the least-squares score would be caught above.
  ls_coefs <- coef(lm(e ~ a + x1 + x2, data = dat))
  expect_gt(max(abs(unname(coef(fx$ps_e)) - unname(ls_coefs))), 0.1)
})

test_that("an additive dose model is stacked at its penalized score", {
  skip_if_not_installed("mgcv")
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, dose_type = "gam")

  # A dose model is read through the registry the single-dose route reads, so an
  # additive fit reaches this route on the same terms: its design is the smooth
  # basis the fit reports, and the equation its coefficients solve is the
  # least-squares score less the penalty its smoothing parameters define. A
  # block written without that penalty would sit away from its root here.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  layout <- expect_joint_root_seeded(spec)

  # The rebuilt weights are held to a looser tolerance here than the other dose
  # models are held to. An additive fit's conditional mean comes back as its
  # basis times its coefficients rather than out of `fitted()`, and the two
  # agree to about 1e-11 rather than to the last bit; the preflight's own guard
  # is 1e-6.
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(fx$wts),
    tolerance = 1e-8
  )

  res <- ipw(fx$models, fx$outcome_mod)
  expect_joint_dose_estimates(res, fx$outcome_mod)

  theta <- coef(res$fit)
  dose_block <- theta[length(coef(fx$ps_a)) + seq_along(coef(fx$ps_e))]
  expect_equal(
    unname(dose_block),
    unname(coef(fx$ps_e)),
    tolerance = 1e-8
  )
})

test_that("ipw() reads a joint additive dose model's design once", {
  skip_if_not_installed("mgcv")
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, dose_type = "gam")

  # The dose component is read through the registry, whose entry evaluates the
  # smooth basis the fit's penalized score is checked at. That basis is the
  # design this route multiplies the dose block by as well, so one call needs
  # one of it. The binary component's design is a lookup and is not counted
  # here: `predict.gam()` is not what builds it.
  counted <- count_gam_designs(ipw(fx$models, fx$outcome_mod))

  expect_identical(counted$builds, 1L)

  # Counting the builds does not change what the call reports.
  expect_joint_dose_estimates(counted$value, fx$outcome_mod)
})

test_that("the joint additive route reports the estimates it reports today", {
  skip_if_not_installed("mgcv")
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, dose_type = "gam")

  # How often the design is built is arithmetic the answer must not see, so a
  # route that builds it once has to reproduce these to the digit.
  res <- ipw(fx$models, fx$outcome_mod)
  expect_equal(
    res$estimates$estimate,
    c(0.6282057290, 0.4153889358, 0.6694390103),
    tolerance = 1e-6
  )
  expect_equal(
    res$estimates$std.err,
    c(0.0747616013, 0.0445820487, 0.0546759779),
    tolerance = 1e-6
  )
})

test_that("ipw() refuses dose weights built from a kernel density", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, .density = "kernel")

  # The bandwidth of a kernel estimate is chosen from the residuals of the dose
  # model, so the weights are not a differentiable function of that model's
  # parameters. The weights are what they claim to be, which is why this is the
  # unavailable-method refusal rather than a mismatch.
  expect_error(
    ipw(fx$models, fx$outcome_mod),
    class = "propensity_ipw_se_method_unavailable_error"
  )

  msg <- joint_error_message(ipw(fx$models, fx$outcome_mod))
  expect_match(msg, "joint", fixed = TRUE)
  expect_match(msg, "wt_joint", fixed = TRUE)
  expect_no_match(msg, "wt_ate", fixed = TRUE)

  expect_propensity_error(ipw(fx$models, fx$outcome_mod))
})

test_that("ipw() refuses a dose psi it cannot write, under the component's name", {
  skip_if_not_installed("MASS")
  dat <- sim_joint_continuous_outliers()

  # A psi of the caller's own is the root of an equation this path does not
  # write, and the single-treatment route refuses it by name. Reached as the
  # second component of a joint intervention the refusal is the same one, and it
  # names the component rather than the container the two treatment models
  # arrived in, since `wt_mod` is not a model the reader can refit.
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  own_psi <- function(u, k = 1.5) pmin(1, k / abs(u))
  ps_e <- suppressWarnings(MASS::rlm(
    e ~ a + x1 + x2,
    data = dat,
    psi = own_psi,
    acc = 1e-10,
    maxit = 200
  ))
  wts <- quiet_wt(wt_joint(
    wt_ate(ps_a),
    wt_ate(
      as.double(fitted(ps_e)),
      dat$e,
      exposure_type = "continuous",
      stabilize = TRUE
    ),
    exposure_type = c("binary", "continuous")
  ))
  outcome_mod <- lm(y ~ a * e, data = dat, weights = wts)
  models <- joint_wt_models(a = ps_a, e = ps_e)

  expect_error(
    ipw(models, outcome_mod),
    class = "propensity_ipw_robust_psi_error"
  )

  msg <- joint_error_message(ipw(models, outcome_mod))
  expect_match(msg, "recognize", fixed = TRUE)
  expect_match(msg, "`e`", fixed = TRUE)
  expect_no_match(msg, "wt_mod", fixed = TRUE)

  expect_propensity_error(ipw(models, outcome_mod))
})

test_that("ipw() refuses a dose link it cannot write the score for", {
  dat <- sim_joint_continuous_positive()

  # The single-treatment route refuses the remaining gaussian links by name,
  # because the coefficients an IRLS iteration stops at under one of them are
  # not a tight enough root to seed the solve from. A dose reached as the second
  # component of a joint intervention is read by the same registry and is
  # refused the same way.
  for (link in c("inverse", "sqrt")) {
    ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
    ps_e <- glm(
      e ~ a + x1 + x2,
      data = dat,
      family = gaussian(link = link)
    )
    wts <- quiet_wt(wt_joint(
      wt_ate(ps_a),
      wt_ate(
        as.double(fitted(ps_e)),
        dat$e,
        exposure_type = "continuous",
        stabilize = TRUE
      ),
      exposure_type = c("binary", "continuous")
    ))
    outcome_mod <- lm(y ~ a * e, data = dat, weights = wts)
    models <- joint_wt_models(a = ps_a, e = ps_e)

    err <- expect_error(
      ipw(models, outcome_mod),
      class = "propensity_ipw_link_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, link, fixed = TRUE)

    # The refused model is the dose component, and the message names it. On
    # this route `wt_mod` is the container the two treatment models arrived in,
    # so a reader told to refit `wt_mod` is told to refit the wrong thing.
    expect_match(msg, "`e`", fixed = TRUE)
    expect_no_match(msg, "wt_mod", fixed = TRUE)
  }
})

test_that("the weights mismatch names the ratio the dose records", {
  dat <- sim_joint_continuous()

  # A dose component built from a model the container does not hold is a weight
  # the system cannot rebuild, and the mismatch is where a caller reads what it
  # tried. Naming the family and the numerator keeps those two off the list of
  # possible causes when they are the two the system already agreed with.
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- lm(e ~ a + x1 + x2, data = dat)
  stale <- lm(e ~ a + x1, data = dat)
  wts <- quiet_wt(wt_joint(
    wt_ate(ps_a),
    wt_ate(
      as.double(fitted(stale)),
      dat$e,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = dens_t(4)
    ),
    exposure_type = c("binary", "continuous")
  ))
  outcome_mod <- lm(y ~ a * e, data = dat, weights = wts)

  models <- joint_wt_models(a = ps_a, e = ps_e)
  expect_error(
    ipw(models, outcome_mod),
    class = "propensity_ipw_weights_mismatch_error"
  )

  msg <- joint_error_message(ipw(models, outcome_mod))
  expect_match(msg, "t(df = 4)", fixed = TRUE)
  expect_match(msg, "marginal", fixed = TRUE)

  # A joint specification targets the joint ate and resolves no focal level, so
  # a focal level is a cause that cannot apply here and naming it sends the
  # reader after a setting the route never read.
  expect_no_match(msg, "focal level", fixed = TRUE)

  expect_propensity_error(ipw(models, outcome_mod))
})

# ---- refusals ----------------------------------------------------------------

test_that("ipw() refuses an outcome model that does not read both treatments", {
  dat <- sim_joint_continuous()

  # The class is the one the package already uses for an outcome model that does
  # not read the exposure, and the message is matched as well: that class also
  # answers a request this route cannot serve at all, so the class alone would
  # be satisfied by a refusal about something else entirely.
  missing_dose <- fit_joint_continuous(dat, outcome_rhs = "a + x1")
  expect_error(
    ipw(missing_dose$models, missing_dose$outcome_mod),
    class = "propensity_ipw_exposure_error",
    regexp = "must contain the exposure"
  )

  missing_binary <- fit_joint_continuous(dat, outcome_rhs = "e + x1")
  expect_error(
    ipw(missing_binary$models, missing_binary$outcome_mod),
    class = "propensity_ipw_exposure_error",
    regexp = "must contain the exposure"
  )
})

test_that("a transformed dose term reports one row per treatment coefficient", {
  dat <- sim_joint_continuous()

  # A curve written without a bare dose term at all. Every treatment-reading
  # coefficient is a row, named by the coefficient rather than by the treatment,
  # and the rows come in the design's column order so the table reads as the
  # coefficient vector reads.
  fx <- fit_joint_continuous(dat, outcome_rhs = "a * sin(e)")
  res <- ipw(fx$models, fx$outcome_mod)

  expect_joint_coefficient_surface(res, c("a", "sin(e)", "a:sin(e)"))
  expect_joint_coefficient_estimates(res, fx$outcome_mod)
  expect_identical(nrow(res$estimates), 3L)

  # Three coefficients for the binomial treatment model, four for the dose model
  # plus the conditional variance its density ratio needs, the two marginal
  # moments the stabilizing numerator is built from, and the four outcome
  # coefficients. Only the outcome block moves when the marginal structural
  # model changes shape.
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 4L)
})

test_that("a covariate in a transformed model contributes no coefficient row", {
  dat <- sim_joint_continuous()

  # The rule that a covariate term contributes no row holds on this surface too,
  # and it has to be said here rather than inferred from the single-treatment
  # route, which reaches the same conclusion through different code.
  #
  # Naming rows by their coefficients makes the temptation concrete: the surface
  # reports coefficients, and every column of the design has one. What decides a
  # row is still whether the term reads a treatment. `x1` does not, so its
  # coefficient is estimated, is carried in the outcome block, and is reported
  # nowhere. Its column sits between the dose column and the interaction column,
  # so a surface that took everything but the intercept would show it in the
  # middle of the table and fail here rather than at the end of a list.
  fx <- fit_joint_continuous(dat, outcome_rhs = "a * sin(e) + x1")
  res <- ipw(fx$models, fx$outcome_mod)

  expect_joint_coefficient_surface(res, c("a", "sin(e)", "a:sin(e)"))
  expect_joint_coefficient_estimates(res, fx$outcome_mod)
  expect_identical(nrow(res$estimates), 3L)
  expect_false("x1" %in% res$estimates$contrast)

  # Five outcome coefficients, one of them the covariate's, against the three
  # rows the surface reports.
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 5L)
  expect_length(coef(fx$outcome_mod), 5L)
})

test_that("a curve beside the bare terms reports every treatment coefficient", {
  dat <- sim_joint_continuous()

  # The bare terms are still there and still contribute their coefficients; the
  # quadratic contributes one more. What the extra term costs is the vocabulary,
  # not the reporting: no row can be called the slope per unit once two columns
  # of the dose carry the response, so every row is named by its coefficient
  # instead and none of them claims a place it is evaluated at.
  fx <- fit_joint_continuous(dat, outcome_rhs = "a + e + I(e^2) + a:e")
  res <- ipw(fx$models, fx$outcome_mod)

  expect_joint_coefficient_surface(res, c("a", "e", "I(e^2)", "a:e"))
  expect_joint_coefficient_estimates(res, fx$outcome_mod)
  expect_identical(nrow(res$estimates), 4L)
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 5L)
})

test_that("a transformed term reading both treatments is named by its coefficient", {
  dat <- sim_joint_continuous()

  # `I(a * sin(e))` reads both treatments without being their bare interaction.
  # On a surface that named rows by what they mean there would be nowhere to put
  # it, since the model then holds two distinct interactions between the same
  # pair of treatments. On a surface that names rows by their coefficients there
  # is no difficulty: the row reports the coefficient of the column
  # `I(a * sin(e))` and says so, which is true and distinguishes it from the
  # bare interaction sitting beside it.
  fx <- fit_joint_continuous(dat, outcome_rhs = "a * e + I(a * sin(e))")
  res <- ipw(fx$models, fx$outcome_mod)

  expect_joint_coefficient_surface(
    res,
    c("a", "e", "I(a * sin(e))", "a:e")
  )
  expect_joint_coefficient_estimates(res, fx$outcome_mod)
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 5L)
})

test_that("a basis in the dose reports its basis and interaction coefficients", {
  dat <- sim_joint_continuous()

  # A basis is where naming rows by the coefficient stops being a refinement and
  # becomes the only option: one term contributes three columns of the dose and
  # three more of its interaction with the binary treatment, and the term names
  # none of them.
  fx <- fit_joint_continuous(dat, outcome_rhs = "a * splines::ns(e, 3)")
  res <- ipw(fx$models, fx$outcome_mod)

  expect_joint_coefficient_surface(
    res,
    c(
      "a",
      "splines::ns(e, 3)1",
      "splines::ns(e, 3)2",
      "splines::ns(e, 3)3",
      "a:splines::ns(e, 3)1",
      "a:splines::ns(e, 3)2",
      "a:splines::ns(e, 3)3"
    )
  )
  expect_joint_coefficient_estimates(res, fx$outcome_mod)
  expect_identical(nrow(res$estimates), 7L)
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 8L)
})

test_that("a basis joint fit needs no .data, unlike the single-treatment route", {
  dat <- sim_joint_continuous()

  # The single-treatment route reads the exposure off the outcome model frame
  # and so requires `.data` for a basis fit, whose frame holds the basis matrix
  # instead. This route reads each treatment off its own treatment model, so the
  # outcome frame is never asked for a column it does not have and the fit runs
  # from the container alone. Supplying `.data` anyway changes nothing.
  fx <- fit_joint_continuous(dat, outcome_rhs = "a * splines::ns(e, 3)")
  bare_call <- ipw(fx$models, fx$outcome_mod)
  with_data <- ipw(fx$models, fx$outcome_mod, .data = dat)

  expect_identical(bare_call$estimates$contrast, with_data$estimates$contrast)
  expect_equal(
    bare_call$estimates$estimate,
    with_data$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    bare_call$estimates$std.err,
    with_data$estimates$std.err,
    tolerance = 1e-10
  )
})

test_that("a logit joint dose model reports coefficient rows on the odds scale", {
  dat <- sim_joint_continuous()

  # The scale word follows the outcome link exactly as it does everywhere else.
  # Under a logit link a coefficient is a log odds ratio per unit of the column
  # it multiplies, whichever treatment that column reads, so the existing word
  # is already the honest one and the surface keeps it. Only the identity link
  # needed a word of its own.
  fx <- fit_joint_continuous(
    dat,
    outcome_rhs = "a * sin(e)",
    outcome_family = "binomial"
  )
  res <- ipw(fx$models, fx$outcome_mod)

  expect_joint_coefficient_surface(
    res,
    c("a", "sin(e)", "a:sin(e)"),
    effect = "log(or)"
  )
  expect_joint_coefficient_estimates(res, fx$outcome_mod)
})

test_that("the conditional reading of a transformed joint fit is the outcome model's", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, outcome_rhs = "a * sin(e)")
  conditional <- causalgenerics::as_conditional(
    ipw(fx$models, fx$outcome_mod)
  )

  # Read through the tidier, which is the surface the reading presents; the
  # stored estimates frame is the marginal one whichever reading is recorded.
  # The conditional reading is the whole coefficient vector, unchanged by the
  # marginal surface having become coefficient-shaped, so the intercept the
  # marginal surface drops is back and neither identity column is there.
  tidied <- tidy(conditional)
  expect_identical(tidied$term, names(coef(fx$outcome_mod)))
  expect_equal(
    tidied$estimate,
    unname(coef(fx$outcome_mod)),
    tolerance = 1e-9
  )
  expect_false("contrast" %in% names(tidied))
  expect_false("group" %in% names(tidied))
})

test_that("ipw() still refuses a term reading a treatment and a covariate", {
  dat <- sim_joint_continuous()

  # The boundary the coefficient surface moves is what a treatment term may look
  # like, not whether a term may read something other than a treatment. A
  # coefficient of a term reading a covariate is a change in an effect per unit
  # of that covariate, so a row named after the column would report a quantity
  # the table pins to no value of it.
  mixed <- fit_joint_continuous(dat, outcome_rhs = "a * e * x1")
  expect_error(
    ipw(mixed$models, mixed$outcome_mod),
    class = "propensity_ipw_msm_error"
  )
  expect_propensity_error(ipw(mixed$models, mixed$outcome_mod))

  dosewise <- fit_joint_continuous(dat, outcome_rhs = "a + e + e:x1")
  expect_error(
    ipw(dosewise$models, dosewise$outcome_mod),
    class = "propensity_ipw_msm_error"
  )
  expect_propensity_error(ipw(dosewise$models, dosewise$outcome_mod))
})

test_that("ipw() refuses a treatment column coded some other way", {
  dat <- sim_joint_continuous()

  # A bare term says which variables a column is built from and nothing about
  # how the column is coded, so a factor treatment under a coding other than
  # treatment contrasts reaches the surface with the term check satisfied. An
  # ordered two-level factor carries polynomial contrasts, whose column is the
  # contrast scaled and centered, so the fit runs to completion and reports
  # neither the effect at a dose of zero nor the slope at the reference level
  # under labels that claim both.
  ordered <- fit_joint_factor(dat, ordered = TRUE)
  expect_error(
    ipw(ordered$models, ordered$outcome_mod),
    class = "propensity_ipw_msm_error"
  )
  expect_propensity_error(ipw(ordered$models, ordered$outcome_mod))

  # The same hole reached through the session option rather than through the
  # vector, which is why the check reads the columns rather than the contrasts
  # a factor happens to carry.
  withr::with_options(list(contrasts = c("contr.sum", "contr.poly")), {
    summed <- fit_joint_factor(dat)
    expect_error(
      ipw(summed$models, summed$outcome_mod),
      class = "propensity_ipw_msm_error"
    )
    expect_propensity_error(ipw(summed$models, summed$outcome_mod))
  })
})

test_that("a factor treatment under treatment contrasts reports the numeric fit", {
  dat <- sim_joint_continuous()

  # Treatment contrasts code a two-level factor as the indicator the rows
  # describe, so the coding check admits it and the two fixtures are the same
  # intervention written two ways. Anchoring one fit against the other is what
  # keeps the check from being satisfied by a refusal of everything: the
  # estimates have to be the numeric fixture's, not merely present.
  numeric_fit <- ipw(
    fit_joint_continuous(dat)$models,
    fit_joint_continuous(dat)$outcome_mod
  )
  factored <- fit_joint_factor(dat)
  factor_fit <- ipw(factored$models, factored$outcome_mod)

  expect_identical(
    factor_fit$estimates$effect,
    numeric_fit$estimates$effect
  )
  expect_equal(
    factor_fit$estimates$estimate,
    numeric_fit$estimates$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    factor_fit$estimates$std.err,
    numeric_fit$estimates$std.err,
    tolerance = 1e-8
  )

  # The two identity columns that name a level are the one thing the fixtures do
  # not share, since a factor's levels are its own and the numeric fixture's are
  # 0 and 1.
  expect_identical(
    factor_fit$estimates$contrast,
    c("a: yes vs no", "e: per unit", "a: yes vs no")
  )
  expect_identical(
    factor_fit$estimates$group,
    c("e = 0", "a = no", "e + 1 vs e")
  )
})

test_that("the weights mismatch names a score the record does not keep", {
  dat <- sim_joint_continuous()

  # A product records the score each component was stabilized on, so a score is
  # a numerator the stacked system rebuilds exactly. What is left is a record
  # written before the slot existed, or one assembled by hand: it says the dose
  # was stabilized on a score and keeps no vector, so the system estimates the
  # dose's marginal moments and reaches a different product.
  #
  # The refusal has to name that cause. Reaching this message and reading
  # through a list of estimands, focal levels, and trimming, none of which
  # applies, is how a correct diagnosis goes unmade.
  fx <- fit_joint_continuous(
    dat,
    dose_score = dnorm(dat$e, mean(dat$e), stats::sd(dat$e))
  )

  wts <- strip_joint_scores(fx$wts)
  outcome_mod <- lm(y ~ a * e, data = dat, weights = wts)

  expect_error(
    ipw(fx$models, outcome_mod),
    class = "propensity_ipw_weights_mismatch_error",
    regexp = "stabilization_score"
  )
  expect_propensity_error(ipw(fx$models, outcome_mod))

  # The observation-level spread is the other cause a dose brings, and it was
  # reported for a single dose before it was reported for this one.
  expect_error(
    ipw(fx$models, outcome_mod),
    class = "propensity_ipw_weights_mismatch_error",
    regexp = "\\.sigma"
  )
})

test_that("a bare-term model with no intercept reports the coefficient surface", {
  dat <- sim_joint_continuous()
  dat$a <- factor(
    ifelse(dat$a == 1, "yes", "no"),
    levels = c("no", "yes")
  )

  # Dropping the intercept expands a factor treatment into one column per level,
  # so the first column is the indicator of the reference level rather than of
  # the focal one and no column of the fit is the 0/1 indicator the vocabulary's
  # rows are claims about. A model the vocabulary has no reading for is reported
  # on the coefficient surface rather than refused: the coefficients are the
  # weighted fit's own either way, and naming each row after the column it
  # multiplies says exactly what each one is without claiming a level contrast
  # none of them is.
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- lm(e ~ a + x1 + x2, data = dat)
  wts <- quiet_wt(wt_joint(
    wt_ate(ps_a),
    wt_ate(
      as.double(fitted(ps_e)),
      dat$e,
      exposure_type = "continuous",
      stabilize = TRUE
    ),
    exposure_type = c("binary", "continuous")
  ))
  models <- joint_wt_models(a = ps_a, e = ps_e)
  outcome_mod <- lm(y ~ a * e - 1, data = dat, weights = wts)

  expect_identical(
    colnames(model.matrix(outcome_mod)),
    c("ano", "ayes", "e", "ayes:e")
  )

  res <- ipw(models, outcome_mod)

  expect_joint_coefficient_surface(res, c("ano", "ayes", "e", "ayes:e"))
  expect_joint_coefficient_estimates(res, outcome_mod)
  expect_identical(nrow(res$estimates), 4L)

  # One row per treatment-reading column, and four of them where a model with an
  # intercept contributes three: the reference-level indicator is a column of
  # this fit and the surface reports it.
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 4L)
})

test_that("a numeric treatment with no intercept reports the coefficient surface", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)
  outcome_mod <- lm(y ~ a * e - 1, data = dat, weights = fx$wts)

  # A 0/1 numeric treatment contributes the same single column whether the
  # intercept is kept or dropped, so the columns here are the ones the
  # vocabulary describes and the coding check has nothing to catch. What the
  # forced zero moves is what the coefficients mean rather than what the columns
  # hold: the coefficient of `a` is the mean at `a = 1` and a dose of zero
  # rather than the difference between the treatment's levels there, so the row
  # the vocabulary would write for it names a contrast the fit does not carry.
  #
  # Which surface a fit reports is decided by whether the vocabulary reads it,
  # and no model without an intercept is one it reads, whatever its columns hold.
  expect_identical(
    colnames(model.matrix(outcome_mod)),
    c("a", "e", "a:e")
  )
  expect_identical(as.integer(attr(terms(outcome_mod), "intercept")), 0L)

  res <- ipw(fx$models, outcome_mod)

  expect_joint_coefficient_surface(res, c("a", "e", "a:e"))
  expect_joint_coefficient_estimates(res, outcome_mod)
  expect_identical(nrow(res$estimates), 3L)
  expect_identical(as.integer(res$fit@n_params), 3L + 5L + 2L + 3L)
})

test_that("a dose transformed in the formula is reported, not refused", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, outcome_rhs = "a * I(e / 10)")

  # Writing the transformation into the marginal structural model's own formula
  # is one of the two ways to work on a rescaled dose, and it is the one that
  # leaves the column the treatment models were fit to alone: the outcome model
  # reads that column and rescales it on its way into the design. The term is no
  # longer bare, so the fit reports the coefficient surface, where each row names
  # the column it multiplies and claims nothing about a step in the dose itself.
  res <- ipw(fx$models, fx$outcome_mod)

  expect_joint_coefficient_surface(res, c("a", "I(e/10)", "a:I(e/10)"))
  expect_joint_coefficient_estimates(res, fx$outcome_mod)
})

test_that("a dose column rescaled after the treatment models are fit says so", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)

  # The dose the estimator holds is the one the treatment models were fit to,
  # and a marginal structural model fit on a rescaled copy of that column is a
  # model of a different dose: its slope is per ten units of the dose the weights
  # were built for, and the weights themselves were built for the original. The
  # refusal is right, and what it has to say is about the dose. A reader who
  # rescaled a column has coded no factor and dropped no intercept, so advice
  # about contrasts and intercepts describes nothing they did.
  rescaled <- dat
  rescaled$e <- rescaled$e / 10
  outcome_mod <- lm(y ~ a * e, data = rescaled, weights = fx$wts)

  cnd <- tryCatch(ipw(fx$models, outcome_mod), error = identity)
  expect_s3_class(cnd, "propensity_ipw_msm_error")
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))

  # What the message has to say: the outcome model's dose column holds values
  # the models the weights came from were not fit to. Either name for those
  # models is admitted, since this route calls them treatment models elsewhere
  # and they are the propensity models a caller knows them as.
  expect_match(message, "same values", fixed = TRUE)
  expect_match(message, "treatment models|propensity")

  # And what it must not say. Scaling a dose is supported; only doing it to one
  # model and not the others is refused, so a remedy about factor coding is
  # wrong here, and the intercept is present in this fit besides.
  expect_no_match(message, "unordered factor", fixed = TRUE)
  expect_no_match(message, "no intercept", fixed = TRUE)

  expect_propensity_error(ipw(fx$models, outcome_mod))
})

test_that("a contrast-coded treatment keeps the contrast advice", {
  dat <- sim_joint_continuous()
  ordered <- fit_joint_factor(dat, ordered = TRUE)

  # The dose advice belongs to a dose column that does not hold what the
  # treatment models were fit to. Here that column is exactly what they were fit
  # to and only the factor's coding moved, so the remedy stays the one about
  # contrasts. The interaction is among the mismatched columns and reads the
  # dose, which is what makes this worth pinning: what decides the advice is
  # which column disagrees, not which variables its term names.
  cnd <- tryCatch(ipw(ordered$models, ordered$outcome_mod), error = identity)
  expect_s3_class(cnd, "propensity_ipw_msm_error")
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))

  expect_match(message, "unordered factor", fixed = TRUE)
  expect_no_match(message, "same values", fixed = TRUE)
})

test_that("ipw() refuses .by on a joint continuous fit", {
  dat <- sim_joint_continuous()
  dat$grp <- factor(ifelse(dat$x1 > 0, "hi", "lo"), levels = c("lo", "hi"))
  # The `.by` refusal fires before the marginal structural model's shape is
  # examined, which is what this fixture pins: the model here reads a covariate
  # in its treatment terms and would be refused on those grounds too, and the
  # error reported is nonetheless the one about `.by`.
  fx <- fit_joint_continuous(dat, outcome_rhs = "a * e * grp")

  # Effect modification of a joint intervention is a three-way question whether
  # the second treatment is a coin flip or a dose, so the refusal is the one the
  # other joint routes raise.
  expect_no_warning(expect_error(
    ipw(fx$models, fx$outcome_mod, .by = grp),
    class = "propensity_ipw_joint_by_error"
  ))
  expect_propensity_error(ipw(fx$models, fx$outcome_mod, .by = grp))
})

test_that("ipw() refuses linearization on a joint continuous fit", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)

  # Every reported coefficient is a parameter of the stacked system, and the
  # linearization path solves no such system.
  expect_error(
    ipw(fx$models, fx$outcome_mod, se_method = "linearization"),
    class = "propensity_ipw_joint_models_method_error"
  )
  expect_propensity_error(
    ipw(fx$models, fx$outcome_mod, se_method = "linearization")
  )
})

# ---- a binary component's own stabilizing numerator -------------------------
#
# A binary component may be stabilized as readily as a dose must be, and the
# numerator it carries is a factor of the product like any other. The system
# estimates it where the single-treatment route estimates it, in a block of its
# own, rather than rebuilding the component unstabilized and reporting the
# difference as a mismatch the caller did not cause.

test_that("the joint route stacks a marginal-stabilized binary component", {
  dat <- sim_joint_continuous()
  fits <- fit_joint_continuous(dat, a_stabilize = TRUE)

  # The product the container holds, written out: the binary factor is the
  # unstabilized one scaled by the marginal probability of the level each unit
  # took, and the dose factor is the density ratio it always was.
  ps_a <- as.numeric(fitted(fits$ps_a))
  p1 <- mean(dat$a)
  binary_wt <- ((dat$a * p1) / ps_a) + (((1 - dat$a) * (1 - p1)) / (1 - ps_a))
  dose_wt <- as.numeric(quiet_wt(wt_ate(
    as.double(fitted(fits$ps_e)),
    dat$e,
    exposure_type = "continuous",
    stabilize = TRUE
  )))

  expect_equal(as.numeric(fits$wts), binary_wt * dose_wt, tolerance = 1e-12)

  # The system rebuilds that product at its seed, which is what the preflight
  # compares, and every equation in it is at its root there.
  spec <- ipw_spec_joint_models(fits$models, fits$outcome_mod)
  expect_joint_weights_at_init(spec, fits$wts)
  expect_joint_root_seeded(spec)

  res <- ipw(fits$models, fits$outcome_mod)
  expect_s3_class(res, "ipw")
  expect_joint_dose_estimates(res, fits$outcome_mod)
})

test_that("a stabilized binary component's proportion is a parameter", {
  dat <- sim_joint_continuous()
  fits <- fit_joint_continuous(dat, a_stabilize = TRUE)
  unstabilized <- fit_joint_continuous(dat)

  res <- ipw(fits$models, fits$outcome_mod)
  res_unstabilized <- ipw(unstabilized$models, unstabilized$outcome_mod)

  # The binary numerator is one marginal proportion, so the system is exactly
  # one parameter wider than the same fit whose binary component carries no
  # numerator at all.
  expect_identical(
    as.integer(res$fit@n_params),
    as.integer(res_unstabilized$fit@n_params) + 1L
  )

  # Blocks are named for the component they belong to, since a joint system
  # carries two of them and a name saying only which role a parameter plays
  # would not say whose.
  theta <- coef(res$fit)
  stab <- theta[grepl("^stab_", names(theta))]
  expect_length(stab, 3L)
  expect_equal(
    unname(theta[["stab_a_pi"]]),
    mean(dat$a),
    tolerance = 1e-8
  )
})

test_that("the joint route stacks a dose stabilized on a numerator model", {
  # The single-dose route estimates a numerator model in its own stabilization
  # block. This route estimates each component's numerator in a block of its
  # own, so a dose whose weights record a fitted numerator is answered with the
  # weights the caller built rather than refused for want of somewhere to put
  # the model.
  dat <- sim_joint_continuous()
  num_mod <- lm(e ~ x2, data = dat)
  fits <- fit_joint_continuous(dat, dose_stabilize = num_mod)

  res <- muffle_coverage_warning(ipw(fits$models, fits$outcome_mod))
  expect_s3_class(res, "ipw")

  # The point estimates are the weighted marginal structural model's own
  # coefficients, which is the oracle every other fit on this route is held to.
  # A system that rebuilt the dose's numerator as anything but the model's would
  # not reach them at all: the preflight compares the weights it rebuilds
  # against the ones the outcome model was fit with.
  expect_equal(
    res$estimates$estimate,
    unname(coef(fits$outcome_mod)[c("a", "e", "a:e")]),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("the dose's numerator model is a block of the joint system", {
  dat <- sim_joint_continuous()
  num_mod <- lm(e ~ x2, data = dat)
  fits <- fit_joint_continuous(dat, dose_stabilize = num_mod)
  marginal <- fit_joint_continuous(dat)

  res <- muffle_coverage_warning(ipw(fits$models, fits$outcome_mod))
  res_marginal <- ipw(marginal$models, marginal$outcome_mod)

  # The default numerator is the dose's own two marginal moments. A fitted one
  # is its coefficients and the spread its density is read at, which is what the
  # single-dose route stacks for the same model, so the system is wider by the
  # difference between the two.
  expect_identical(
    as.integer(res$fit@n_params),
    as.integer(res_marginal$fit@n_params) + length(coef(num_mod)) + 1L - 2L
  )

  # The block is solved at the model it was given rather than carried at
  # whatever the seed was, with the spread after the coefficients as the
  # single-dose route orders them.
  theta <- coef(res$fit)
  stab <- theta[grepl("^stab_", names(theta))]
  expect_identical(length(stab), length(coef(num_mod)) + 1L)
  expect_equal(
    unname(stab[seq_along(coef(num_mod))]),
    unname(coef(num_mod)),
    tolerance = 1e-6
  )
})

test_that("the joint route stacks a binary component's numerator model", {
  # Both components may carry a fitted numerator, and each one is estimated in
  # its own block. A binary component's numerator is the probability of the
  # level each unit took, so the product the container holds is the dose's
  # density ratio times that probability over the binary denominator.
  dat <- sim_joint_continuous()
  num_a <- glm(a ~ x2, data = dat, family = binomial())
  fits <- fit_joint_continuous(dat, a_stabilize = num_a)

  ps_a <- as.numeric(fitted(fits$ps_a))
  p_a <- as.numeric(fitted(num_a))
  binary_wt <- ((dat$a / ps_a) + ((1 - dat$a) / (1 - ps_a))) *
    (dat$a * p_a + (1 - dat$a) * (1 - p_a))
  dose_wt <- as.numeric(quiet_wt(wt_ate(
    as.double(fitted(fits$ps_e)),
    dat$e,
    exposure_type = "continuous",
    stabilize = TRUE
  )))

  expect_equal(as.numeric(fits$wts), binary_wt * dose_wt, tolerance = 1e-12)

  res <- muffle_coverage_warning(ipw(fits$models, fits$outcome_mod))
  expect_s3_class(res, "ipw")
  expect_equal(
    res$estimates$estimate,
    unname(coef(fits$outcome_mod)[c("a", "e", "a:e")]),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res$estimates$std.err)))

  # The binary numerator's coefficients are parameters of the system, alongside
  # the dose's own two marginal moments.
  theta <- coef(res$fit)
  expect_identical(
    as.integer(sum(grepl("^stab_", names(theta)))),
    length(coef(num_a)) + 2L
  )
})

# ---- a component stabilized on a fixed score --------------------------------
#
# A numerator the caller computed is a vector rather than a model, and the
# product records it per component. The stacked system reads it back and holds
# it fixed, which is what the single-treatment routes do with a score: it
# multiplies the weights and estimates nothing.

test_that("the joint route rebuilds a dose stabilized on a fixed score", {
  dat <- sim_joint_continuous()

  # The numerator a caller would write by hand, whose `sd()` divides by n - 1
  # where the estimator's own moment divides by n. The two differ by more than
  # the preflight's tolerance, so a system that stood the marginal moments in
  # would be reported as a mismatch rather than reaching these estimates.
  score <- dnorm(dat$e, mean(dat$e), stats::sd(dat$e))
  fx <- fit_joint_continuous(dat, dose_score = score)

  # The product written out: the binary factor is the unstabilized one, and the
  # dose factor is the score over the conditional density.
  ps_a <- as.numeric(fitted(fx$ps_a))
  binary_wt <- (dat$a / ps_a) + ((1 - dat$a) / (1 - ps_a))
  mu <- as.numeric(fitted(fx$ps_e))
  sigma <- sqrt(mean(residuals(fx$ps_e)^2))
  dose_wt <- score / dnorm(dat$e, mu, sigma)

  expect_equal(as.numeric(fx$wts), binary_wt * dose_wt, tolerance = 1e-12)

  # The system rebuilds that product at its seed, which is what the preflight
  # compares, and every equation in it is at its root there.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  expect_joint_weights_at_init(spec, fx$wts)
  expect_joint_root_seeded(spec)

  res <- ipw(fx$models, fx$outcome_mod)
  expect_s3_class(res, "ipw")
  expect_equal(
    res$estimates$estimate,
    unname(coef(fx$outcome_mod)[c("a", "e", "a:e")]),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("a fixed score on the dose is no parameter of the system", {
  dat <- sim_joint_continuous()
  score <- dnorm(dat$e, mean(dat$e), stats::sd(dat$e))

  fixed <- fit_joint_continuous(dat, dose_score = score)
  marginal <- fit_joint_continuous(dat)

  res <- ipw(fixed$models, fixed$outcome_mod)
  res_marginal <- ipw(marginal$models, marginal$outcome_mod)

  # A score is a known multiplier rather than a quantity the system estimates,
  # so the dose's slice of the stabilization block is empty where the default
  # numerator's holds the exposure's two marginal moments.
  expect_identical(
    as.integer(res$fit@n_params),
    as.integer(res_marginal$fit@n_params) - 2L
  )
  expect_false(any(grepl("^stab_", names(coef(res$fit)))))
})

test_that("the joint route rebuilds a binary component stabilized on a score", {
  dat <- sim_joint_continuous()

  # A discrete component's score is the one the product used to record nothing
  # at all about: `stabilize = TRUE` and a score of the caller's left the same
  # record, so a score that differed from the marginal proportion was reported
  # as a mismatch naming a component the caller never stabilized by hand.
  score <- 0.3 + 0.4 * plogis(dat$x2)
  fx <- fit_joint_continuous(dat, a_stabilize = TRUE, a_score = score)

  ps_a <- as.numeric(fitted(fx$ps_a))
  binary_wt <- ((dat$a / ps_a) + ((1 - dat$a) / (1 - ps_a))) * score
  dose_wt <- as.numeric(quiet_wt(wt_ate(
    as.double(fitted(fx$ps_e)),
    dat$e,
    exposure_type = "continuous",
    stabilize = TRUE
  )))

  expect_equal(as.numeric(fx$wts), binary_wt * dose_wt, tolerance = 1e-12)

  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  expect_joint_weights_at_init(spec, fx$wts)
  expect_joint_root_seeded(spec)

  res <- ipw(fx$models, fx$outcome_mod)
  expect_equal(
    res$estimates$estimate,
    unname(coef(fx$outcome_mod)[c("a", "e", "a:e")]),
    tolerance = 1e-8
  )

  # The binary component estimates nothing for its numerator, so the only
  # stabilization parameters are the dose's own two marginal moments.
  theta <- coef(res$fit)
  expect_identical(sum(grepl("^stab_", names(theta))), 2L)
  expect_false(any(grepl("^stab_a_", names(theta))))
})

test_that("the weights mismatch names a component whose score was dropped", {
  dat <- sim_joint_continuous()

  # A score holds one value per observation, so an operation that changes the
  # length of the weights drops it. What is left is a component the record says
  # was stabilized and holds no vector for, which is what a component
  # stabilized on the marginal proportion looks like: the system stands the
  # marginal proportion in, reaches a different product, and the preflight is
  # where the caller finds out. The refusal has to name the component whose
  # numerator went missing, since the comparison sees the product alone.
  score <- 0.3 + 0.4 * plogis(dat$x2)
  fx <- fit_joint_continuous(dat, a_stabilize = TRUE, a_score = score)

  half <- seq_len(floor(nrow(dat) / 2))
  expect_warning(
    first <- fx$wts[half],
    class = "propensity_stabilization_score_warning"
  )
  expect_warning(
    second <- fx$wts[-half],
    class = "propensity_stabilization_score_warning"
  )

  # The halves hold every observation between them and each of them dropped the
  # same component's score, so the combined weights are the product the outcome
  # model was fit with, carrying a record that says the score is gone.
  wts <- vctrs::vec_c(first, second)

  expect_identical(joint_wt_meta(wts)$score_dropped, c(TRUE, FALSE))
  expect_equal(as.numeric(wts), as.numeric(fx$wts), tolerance = 1e-12)

  outcome_mod <- lm(y ~ a * e, data = dat, weights = wts)

  cnd <- expect_error(
    ipw(fx$models, outcome_mod),
    class = "propensity_ipw_weights_mismatch_error"
  )
  expect_match(conditionMessage(cnd), "stabilization_score")
  expect_match(conditionMessage(cnd), "`a`", fixed = TRUE)
})

test_that("a component's score left stale by a model frame is refused", {
  dat <- sim_joint_continuous()
  dat$w <- rev(dat$x1)
  dat$w[order(dat$e, decreasing = TRUE)[seq_len(10)]] <- NA
  kept <- !is.na(dat$w)

  # The other way the same drop arrives. Subsetting the weights empties the
  # slot and marks it, which the block above pins; a model frame instead drops
  # the incomplete rows in C and re-attaches the record whole, so the score
  # reaches `ipw()` at the length the product was built at rather than at the
  # length of the rows it is about to weight. Multiplying that score into a
  # density read over the rows that remain recycles, which base R reports in
  # terms of no argument the caller wrote, so the refusal comes first and names
  # the component whose score is the wrong length.
  dose_score <- dnorm(dat$e, mean(dat$e), stats::sd(dat$e))

  # The drop is reported where the outcome model restricts the rows, and is
  # asserted here so that the only condition left for the call below is the one
  # being pinned.
  expect_warning(
    fx <- fit_joint_continuous(
      dat,
      dose_score = dose_score,
      outcome_rhs = "a * e + w"
    ),
    class = "propensity_stabilization_score_warning"
  )

  seen <- character()
  cnd <- withCallingHandlers(
    tryCatch(
      ipw(fx$models, fx$outcome_mod, .data = dat[kept, ]),
      error = function(e) e
    ),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(seen, character())
  expect_s3_class(cnd, "propensity_ipw_stabilization_score_error")
  expect_match(conditionMessage(cnd), "stabilization_score")
  expect_match(conditionMessage(cnd), "`e`", fixed = TRUE)
})

test_that("a record that keeps no score is read as one that records none", {
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, a_stabilize = TRUE)

  # A product built before the slot existed holds a record one element short,
  # which says what a record whose components each record no score says. The
  # fit it supports is the fit it supported then.
  wts <- strip_joint_scores(fx$wts)
  outcome_mod <- lm(y ~ a * e, data = dat, weights = wts)

  res <- expect_silent(ipw(fx$models, outcome_mod))
  expect_equal(
    res$estimates$estimate,
    ipw(fx$models, fx$outcome_mod)$estimates$estimate,
    tolerance = 1e-10
  )
})

# ---- a dose model fit under case weights ------------------------------------
#
# The dose block this route stacks is the unweighted score of whichever class
# the registry read, so a dose model fit with prior case weights sits away from
# the root of the equation written for it and the solve, seeded at its
# coefficients, would walk off the fit the user has. The single-dose route
# refuses such a fit by name for every class it stacks, and this route reads the
# same field of the same fitted objects, so the refusal is the same one and is
# made of each class in turn.
#
# Where each class records prior weights differs: an `lm` and an `rlm` keep them
# in `weights`, and a `glm` and a `gam` in `prior.weights`. A guard reading one
# field would pass half of these fits through.

test_that("ipw() refuses a joint dose model fit with case weights", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("mgcv")
  dat <- sim_joint_continuous()
  dat$cw <- rep(c(1, 2), length.out = nrow(dat))

  weighted_dose <- function(dose_type) {
    switch(
      dose_type,
      lm = lm(e ~ a + x1 + x2, data = dat, weights = cw),
      glm = glm(
        e ~ a + x1 + x2,
        data = dat,
        family = gaussian(),
        weights = cw
      ),
      rlm = MASS::rlm(
        e ~ a + x1 + x2,
        data = dat,
        weights = dat$cw,
        acc = 1e-10,
        maxit = 200
      ),
      gam = mgcv::gam(e ~ a + s(x1) + x2, data = dat, weights = cw)
    )
  }

  for (dose_type in c("lm", "glm", "rlm", "gam")) {
    ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
    ps_e <- weighted_dose(dose_type)
    wts <- quiet_wt(wt_joint(
      wt_ate(ps_a),
      wt_ate(
        as.double(fitted(ps_e)),
        dat$e,
        exposure_type = "continuous",
        stabilize = TRUE
      ),
      exposure_type = c("binary", "continuous")
    ))
    outcome_mod <- lm(y ~ a * e, data = dat, weights = wts)
    models <- joint_wt_models(a = ps_a, e = ps_e)

    expect_error(
      ipw(models, outcome_mod),
      class = "propensity_ipw_ps_weights_error",
      label = dose_type
    )

    msg <- joint_error_message(ipw(models, outcome_mod))
    expect_match(msg, "`e`", fixed = TRUE)
    expect_no_match(msg, "wt_mod", fixed = TRUE)
  }
})

test_that("the joint dose weights refusal comes before the estimates do", {
  skip_if_not_installed("mgcv")
  dat <- sim_joint_continuous()
  dat$cw <- rep(c(1, 2), length.out = nrow(dat))

  # The failure this guard prevents is silent: without it the route returns a
  # full table of estimates whose dose block sits at the unweighted score. The
  # unweighted fit of the same columns is a different model, which is what makes
  # the drift real rather than notional.
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- lm(e ~ a + x1 + x2, data = dat, weights = cw)
  unweighted <- lm(e ~ a + x1 + x2, data = dat)
  expect_gt(max(abs(unname(coef(ps_e)) - unname(coef(unweighted)))), 1e-3)

  wts <- quiet_wt(wt_joint(
    wt_ate(ps_a),
    wt_ate(
      as.double(fitted(ps_e)),
      dat$e,
      exposure_type = "continuous",
      stabilize = TRUE
    ),
    exposure_type = c("binary", "continuous")
  ))
  outcome_mod <- lm(y ~ a * e, data = dat, weights = wts)

  expect_propensity_error(ipw(joint_wt_models(a = ps_a, e = ps_e), outcome_mod))
})

# ---- a component's numerator design rebuilt from `.data` --------------------
#
# A numerator model's design is one of the designs `ipw()` rebuilds when the
# caller supplies `.data`, and a component's numerator is no different from the
# single-treatment routes' one: it is rebuilt over the rows every model read,
# under the coding the fit recorded, and out of `.data` rather than out of a
# frame the fit may no longer keep. Read off the fit's own frame instead, it is
# a design over other rows than everything it is stacked with, and a column
# `.data` supplies as another type never reaches it at all.
#
# Both block builders are exercised here, since a product carries one numerator
# per component and the two are built by the routes their own exposures use.

# Whitespace in a cli-formatted message wraps where the console is narrow, so
# the message is flattened before anything is matched in it.
joint_numerator_ipw_message <- function(cnd) {
  gsub("[[:space:]]+", " ", conditionMessage(cnd))
}

# A crossing shaped like the one `sim_joint_continuous()` simulates, with two
# more baseline covariates: one only the binary component's numerator reads and
# one only the dose's does. `sim_joint_continuous()` has no column of either
# kind, and they are the columns every test below asks `.data` for.
sim_joint_continuous_numerator <- function(seed = 8811, n = 700) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  vb <- rnorm(n)
  vd <- rnorm(n)
  a <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * x2 + 0.3 * vb))
  e <- 0.4 + 0.5 * x1 - 0.3 * x2 - 0.8 * a + 0.4 * vd + rnorm(n)
  y <- 1 + 0.7 * a + 0.5 * e + 0.6 * a * e + rnorm(n)

  data.frame(x1, x2, vb, vd, a, e, y)
}

# The fixture above with a covariate only the outcome model reads, one of whose
# values is missing. The five fits then read one frame and keep different rows
# of it: the outcome model drops the incomplete row and the two treatment models
# and the two numerator models, which never read the column, keep it. A `.data`
# holding the frame all five were given therefore has a row to drop before any
# design is built over it.
sim_joint_continuous_numerator_gap <- function(seed = 8811, n = 700) {
  dat <- sim_joint_continuous_numerator(seed = seed, n = n)
  dat$w <- rev(dat$x1)
  dat$w[[11]] <- NA

  dat
}

# The product both of whose components carry a numerator model, with the
# numerator models handed back alongside the fits the route reads.
joint_numerator_data_fits <- function(dat, outcome_rhs = "a * e") {
  num_a <- glm(a ~ vb, data = dat, family = binomial())
  num_e <- lm(e ~ vd, data = dat)
  fits <- fit_joint_continuous(
    dat,
    outcome_rhs = outcome_rhs,
    a_stabilize = num_a,
    dose_stabilize = num_e
  )

  c(fits, list(num_a = num_a, num_e = num_e))
}

test_that("a component's numerator design is restricted to the rows .data keeps", {
  dat <- sim_joint_continuous_numerator_gap()
  fits <- joint_numerator_data_fits(dat, outcome_rhs = "a * e + w")
  kept <- !is.na(dat$w)

  # Supplying the frame the fits were given and supplying the rows `ipw()`
  # restricts it to are the same request, so they report the same thing. Each
  # component's numerator design is one of the designs that restriction is for.
  res_given <- muffle_coverage_warning(ipw(
    fits$models,
    fits$outcome_mod,
    .data = dat
  ))
  res_kept <- muffle_coverage_warning(ipw(
    fits$models,
    fits$outcome_mod,
    .data = dat[kept, ]
  ))

  expect_s3_class(res_given, "ipw")
  expect_equal(
    res_given$estimates$estimate,
    res_kept$estimates$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    res_given$estimates$std.err,
    res_kept$estimates$std.err,
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res_given$estimates$std.err)))
})

test_that("the stacked numerator blocks solve over the rows .data keeps", {
  dat <- sim_joint_continuous_numerator_gap()
  fits <- joint_numerator_data_fits(dat, outcome_rhs = "a * e + w")
  kept <- !is.na(dat$w)

  # Each component's numerator contributes a block that reads no parameter from
  # anywhere else in the stack: the binary one is its own binomial score, which
  # is exactly identified, and the dose's is its least-squares score together
  # with the second moment of its residuals, whose root is the normal equations
  # of the design the block carries followed by the mean square of the residuals
  # those coefficients leave. Both roots are unique, so the rows each design
  # carries are the rows its solved block is the fit over. The numerators
  # arrived fit to the whole frame and the system reads them over the rows the
  # outcome model kept, so what they solve to are the refits on those rows
  # rather than the coefficients they came with. Both differ by more than the
  # tolerance, which is what makes the pins say anything.
  res <- muffle_coverage_warning(ipw(
    fits$models,
    fits$outcome_mod,
    .data = dat
  ))
  refit_a <- glm(a ~ vb, data = dat[kept, ], family = binomial())
  refit_e <- lm(e ~ vd, data = dat[kept, ])

  theta <- coef(res$fit)
  names_a <- paste0("stab_a_", names(coef(refit_a)))
  names_e <- paste0("stab_e_", names(coef(refit_e)))

  expect_true(all(c(names_a, names_e, "stab_e_sigma2") %in% names(theta)))
  expect_equal(unname(theta[names_a]), unname(coef(refit_a)), tolerance = 1e-6)
  expect_equal(unname(theta[names_e]), unname(coef(refit_e)), tolerance = 1e-6)
  expect_equal(
    unname(theta[["stab_e_sigma2"]]),
    mean(residuals(refit_e)^2),
    tolerance = 1e-6
  )

  expect_false(isTRUE(all.equal(
    unname(coef(refit_a)),
    unname(coef(fits$num_a)),
    tolerance = 1e-6
  )))
  expect_false(isTRUE(all.equal(
    unname(coef(refit_e)),
    unname(coef(fits$num_e)),
    tolerance = 1e-6
  )))
})

# ---- a numerator covariate `.data` supplies as another type -----------------
#
# A numerator model is rebuilt from `.data` the way the other models are, so a
# column it alone reads is a column the rebuild can be given as the wrong type,
# or not given at all. Which fit a refusal names is this route's own vocabulary,
# where the models arrive under names the caller gave them, so what is pinned
# here is the column and the argument the frame arrived in.

# A three-level factor over the same rows, which nothing models. It is what a
# numeric numerator covariate is supplied as below.
joint_numerator_grouping <- function(dat) {
  factor(c("a", "b", "c")[1 + (rank(dat$x1) %% 3)], levels = c("a", "b", "c"))
}

test_that("a binary component's numerator covariate supplied as a factor is refused", {
  dat <- sim_joint_continuous_numerator()
  fits <- joint_numerator_data_fits(dat)

  # The type sweep compares the class the fit recorded for the column with the
  # class `.data` supplies and refuses the pair before any design is rebuilt.
  # What it heads off is a factor of three levels taking two design columns
  # where the number it stands in for took one.
  supplied <- dat
  supplied$vb <- joint_numerator_grouping(dat)

  err <- expect_error(
    muffle_coverage_warning(ipw(
      fits$models,
      fits$outcome_mod,
      .data = supplied
    )),
    class = "propensity_ipw_data_error"
  )

  message <- joint_numerator_ipw_message(err)
  expect_match(message, "vb", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)
})

test_that("a dose component's numerator covariate supplied as a factor is refused", {
  dat <- sim_joint_continuous_numerator()
  fits <- joint_numerator_data_fits(dat)

  supplied <- dat
  supplied$vd <- joint_numerator_grouping(dat)

  err <- expect_error(
    muffle_coverage_warning(ipw(
      fits$models,
      fits$outcome_mod,
      .data = supplied
    )),
    class = "propensity_ipw_data_error"
  )

  message <- joint_numerator_ipw_message(err)
  expect_match(message, "vd", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)
})

test_that("a binary component's numerator covariate absent from .data is refused", {
  dat <- sim_joint_continuous_numerator()
  fits <- joint_numerator_data_fits(dat)

  # The columns the rebuilds read are asked for before any of them runs, and a
  # numerator model's covariates are among them. Left out of the set, a column
  # only a numerator reads reaches `model.matrix()` as an object that is not
  # there.
  supplied <- dat
  supplied$vb <- NULL

  err <- expect_error(
    muffle_coverage_warning(ipw(
      fits$models,
      fits$outcome_mod,
      .data = supplied
    )),
    class = "propensity_columns_exist_error"
  )

  expect_match(joint_numerator_ipw_message(err), "vb", fixed = TRUE)
})

test_that("a dose component's numerator covariate absent from .data is refused", {
  dat <- sim_joint_continuous_numerator()
  fits <- joint_numerator_data_fits(dat)

  supplied <- dat
  supplied$vd <- NULL

  err <- expect_error(
    muffle_coverage_warning(ipw(
      fits$models,
      fits$outcome_mod,
      .data = supplied
    )),
    class = "propensity_columns_exist_error"
  )

  expect_match(joint_numerator_ipw_message(err), "vd", fixed = TRUE)
})

# ---- a dose component's integrated numerator read over its fit's rows -------
#
# A component's integrated numerator is the same reading the single-dose route
# takes: the conditional density averaged over the units, interpolated from a
# grid spanning their dose. The dose model was fit over one set of rows and the
# weights carry that reading, and a `.data` that leaves fewer rows changes both
# halves of it at once. The rebuild owes the fit's reading here as it does
# there, alongside the spreads and the moments the component's blocks are
# seeded at.

# The crossing `sim_joint_continuous_numerator()` simulates with a covariate
# only the outcome model reads, withheld at the ten largest doses. Withholding
# it at one arbitrary row moves the numerator's average without moving the
# grid; withholding it there moves the grid's upper end as well, so both halves
# of the reading are exercised.
sim_joint_continuous_grid_gap <- function(seed = 8811, n = 700) {
  dat <- sim_joint_continuous_numerator(seed = seed, n = n)
  dat$w <- rev(dat$x1)
  dat$w[order(dat$e, decreasing = TRUE)[seq_len(10)]] <- NA

  dat
}

test_that("a dose component's integrated numerator survives the rows .data drops", {
  dat <- sim_joint_continuous_grid_gap()
  kept <- !is.na(dat$w)
  fits <- fit_joint_continuous(
    dat,
    outcome_rhs = "a * e + w",
    numerator = "integrated"
  )

  # The grid's upper end moves under the restriction, which is what makes the
  # agreement below a statement about which rows the numerator was read over.
  expect_gt(max(dat$e), max(dat$e[kept]))

  res <- ipw(fits$models, fits$outcome_mod, .data = dat)
  res_kept <- ipw(fits$models, fits$outcome_mod, .data = dat[kept, ])

  expect_s3_class(res, "ipw")
  expect_equal(
    res$estimates$estimate,
    res_kept$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    res$estimates$std.err,
    res_kept$estimates$std.err,
    tolerance = 1e-10
  )
  expect_true(all(is.finite(res$estimates$std.err)))
  expect_true(all(res$estimates$std.err > 0))
})

test_that("the dose numerator the joint route rebuilds is the one its weights carry", {
  dat <- sim_joint_continuous_grid_gap()
  kept <- !is.na(dat$w)
  fits <- fit_joint_continuous(
    dat,
    outcome_rhs = "a * e + w",
    numerator = "integrated"
  )

  # An integrated numerator estimates nothing, so the component carries no
  # stabilization block and the only thing that can move its half of the
  # product is the reading itself.
  spec <- ipw_spec_joint_models(fits$models, fits$outcome_mod, .data = dat)
  layout <- expect_joint_weights_at_init(spec, as.double(fits$wts)[kept])
  expect_length(layout$idx$stab, 0L)
})
