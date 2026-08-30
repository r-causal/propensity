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
    wt_ate(ps_a),
    wt_ate(
      as.double(fitted(ps_e)),
      dat$e,
      exposure_type = "continuous",
      stabilize = TRUE,
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

test_that("ipw() refuses a dose model whose score it cannot write", {
  skip_if_not_installed("mgcv")
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, dose_type = "gam")

  # An additive model chooses how much to smooth by REML, and no equation
  # stacked here reproduces that choice, so what is unavailable is the standard
  # error method rather than the weights, which are exactly what they claim to
  # be. The single-dose route says the same thing about the same fit.
  expect_error(
    ipw(fx$models, fx$outcome_mod),
    class = "propensity_ipw_se_method_unavailable_error"
  )

  # Both routes point at a bootstrap the user writes, and this one says it in
  # its own terms: the weights a replicate rebuilds are the product weights of
  # two treatment models rather than the dose ratio of one.
  msg <- joint_error_message(ipw(fx$models, fx$outcome_mod))
  expect_match(msg, "joint", fixed = TRUE)
  expect_match(msg, "wt_joint", fixed = TRUE)
  expect_no_match(msg, "wt_ate", fixed = TRUE)

  expect_propensity_error(ipw(fx$models, fx$outcome_mod))
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

  # A redescending psi is the root of an equation this path does not write, and
  # the single-treatment route refuses it by name. Reached as the second
  # component of a joint intervention the refusal is the same one, and it names
  # the component rather than the container the two treatment models arrived
  # in, since `wt_mod` is not a model the reader can refit.
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- MASS::rlm(
    e ~ a + x1 + x2,
    data = dat,
    psi = MASS::psi.bisquare,
    acc = 1e-10
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

  expect_error(
    ipw(models, outcome_mod),
    class = "propensity_ipw_robust_psi_error"
  )

  msg <- joint_error_message(ipw(models, outcome_mod))
  expect_match(msg, "psi.bisquare", fixed = TRUE)
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

test_that("the weights mismatch names a fixed stabilization score", {
  dat <- sim_joint_continuous()

  # A product weight records that a component was stabilized and not which
  # numerator it was built with, so a dose component carrying one of its own is
  # a weight the stacked system cannot rebuild: it estimates the dose's marginal
  # moments and reaches a different vector. The numerator here is the one a
  # caller would write by hand, whose `sd()` divides by n - 1 where the
  # estimator's own moment divides by n, so the two differ by more than the
  # preflight's tolerance and by little else.
  #
  # The refusal has to name that cause. Reaching this message and reading
  # through a list of estimands, focal levels, and trimming, none of which
  # applies, is how a correct diagnosis goes unmade.
  fx <- fit_joint_continuous(
    dat,
    dose_score = dnorm(dat$e, mean(dat$e), stats::sd(dat$e))
  )

  expect_error(
    ipw(fx$models, fx$outcome_mod),
    class = "propensity_ipw_weights_mismatch_error",
    regexp = "stabilization_score"
  )
  expect_propensity_error(ipw(fx$models, fx$outcome_mod))

  # The observation-level spread is the other cause a dose brings, and it was
  # reported for a single dose before it was reported for this one.
  expect_error(
    ipw(fx$models, fx$outcome_mod),
    class = "propensity_ipw_weights_mismatch_error",
    regexp = "\\.sigma"
  )
})

test_that("a bare-term model with no intercept is refused, not errored", {
  dat <- sim_joint_continuous()
  dat$a <- factor(
    ifelse(dat$a == 1, "yes", "no"),
    levels = c("no", "yes")
  )

  # Dropping the intercept expands a factor treatment into one column per level,
  # so the first column is the indicator of the reference level rather than of
  # the focal one. The terms are bare, so the model reaches the vocabulary
  # surface, and it is the columns that say it cannot be reported there.
  #
  # What is pinned is that this arrives as a guided refusal. Reading the columns
  # is what makes it one: a surface that trusted the terms would have taken the
  # first column for the indicator and reported a coefficient that is not the
  # effect the row names.
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
  expect_error(
    ipw(models, outcome_mod),
    class = "propensity_ipw_msm_error"
  )

  # The remedy has to name the cause a reader of this fit is looking at. The
  # treatment here is already an unordered factor under treatment contrasts, so
  # the second half of the remedy describes what was done; only the dropped
  # intercept explains why that was not enough.
  expect_error(
    ipw(models, outcome_mod),
    regexp = "no intercept"
  )

  expect_propensity_error(ipw(models, outcome_mod))
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
