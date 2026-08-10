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
# ---- the coefficient-shaped surface -----------------------------------------
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
  outcome_family = "gaussian"
) {
  ps_a <- glm(
    stats::reformulate(a_rhs, response = "a"),
    data = dat,
    family = binomial()
  )
  ps_e <- lm(stats::reformulate(e_rhs, response = "e"), data = dat)

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
  skip_if_not_installed("deli")
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
  skip_if_not_installed("deli")
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
  skip_if_not_installed("deli")
  dat <- sim_dose_curve()
  fx <- fit_dose_msm(dat, "e + I(e^2)")
  res <- ipw(fx$ps_mod, fx$outcome_mod)
  est <- res$estimates

  # More than one exposure coefficient, so the table gains a contrast column
  # naming the term each one belongs to, placed where every other route places
  # it. The scale stays the one the outcome link implies.
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
  expect_identical(est$effect, c("slope", "slope"))
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
  expect_equal(sqrt(diag(vcov(res))), est$std.err, ignore_attr = TRUE)
  expect_identical(names(coef(res)), paste(est$effect, est$contrast))

  # One outcome coefficient more than the single-term fit, and nothing else
  # about the system changes.
  expect_identical(as.integer(res$fit@n_params), 9L)
})

test_that("the relaxation admits any term that reads the exposure alone", {
  skip_if_not_installed("deli")
  dat <- sim_dose_curve()

  # The boundary is variable membership, not the shape of the function: what
  # made a multi-term model unreportable was a term reading something other than
  # the exposure, and a curve in the exposure alone is as reportable as a
  # quadratic in it. Restricting to polynomials would also mean deciding by
  # inspection which expressions count as one.
  fx <- fit_dose_msm(dat, "e + sin(e)")
  res <- ipw(fx$ps_mod, fx$outcome_mod)

  expect_identical(nrow(res$estimates), 2L)
  expect_identical(res$estimates$contrast, c("e", "sin(e)"))
  expect_equal(
    res$estimates$estimate,
    unname(coef(fx$outcome_mod)[c("e", "sin(e)")]),
    tolerance = 1e-8
  )
})

# ---- part 1: the joint route with a dose ------------------------------------

test_that("ipw() over a binary and a continuous treatment reports the coefficient surface", {
  skip_if_not_installed("deli")
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat)
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates

  expect_s3_class(res, "ipw")
  expect_identical(res$estimand, "ate")
  expect_identical(res$se_method, "mestimation")

  # Coefficient-shaped, not cell-shaped: a dose has no cells, so there are no
  # counterfactual risks to report and no `mean` rows.
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
  skip_if_not_installed("deli")
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
  skip_if_not_installed("deli")
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
  skip_if_not_installed("deli")
  dat <- sim_joint_continuous()
  fx <- fit_joint_continuous(dat, outcome_rhs = "a + e")
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates

  # With no interaction term the model forces one effect for the binary
  # treatment at every dose and one slope at either level of it, so neither row
  # is evaluated anywhere in particular and there is no interaction to report.
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

# ---- standard errors and the stacked system ---------------------------------

test_that("a joint continuous fit reports a usable standard error for every row", {
  skip_if_not_installed("deli")
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
  skip_if_not_installed("deli")
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
  skip_if_not_installed("deli")
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
  skip_if_not_installed("deli")
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

# ---- refusals ----------------------------------------------------------------

test_that("ipw() refuses an outcome model that does not read both treatments", {
  skip_if_not_installed("deli")
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

test_that("ipw() refuses a treatment term that is not a bare one", {
  skip_if_not_installed("deli")
  dat <- sim_joint_continuous()

  # Every row of this surface says where in the other treatment's range it is
  # evaluated, and each of those readings holds only of a model that is linear
  # in each treatment. A transformed treatment term has no such reading, and it
  # lands in a row whose label then describes neither the coefficient reported
  # nor, where two transformed terms of one treatment land in the same row, one
  # row rather than two.
  curved <- fit_joint_continuous(dat, outcome_rhs = "a + e + I(e^2) + a:e")
  expect_error(
    ipw(curved$models, curved$outcome_mod),
    class = "propensity_ipw_msm_error"
  )
  expect_propensity_error(ipw(curved$models, curved$outcome_mod))

  # A curve written without a bare dose term at all is refused on the same
  # terms: the boundary is what a term is, not whether a linear one sits beside
  # it.
  transformed <- fit_joint_continuous(dat, outcome_rhs = "a * sin(e)")
  expect_error(
    ipw(transformed$models, transformed$outcome_mod),
    class = "propensity_ipw_msm_error"
  )
  expect_propensity_error(ipw(transformed$models, transformed$outcome_mod))

  # A transformed term reading both treatments is a second interaction between
  # them, which the interaction row cannot hold beside the one it reports.
  wrapped <- fit_joint_continuous(
    dat,
    outcome_rhs = "a * e + I(a * sin(e))"
  )
  expect_error(
    ipw(wrapped$models, wrapped$outcome_mod),
    class = "propensity_ipw_msm_error"
  )
  expect_propensity_error(ipw(wrapped$models, wrapped$outcome_mod))
})

test_that("ipw() refuses a treatment column coded some other way", {
  skip_if_not_installed("deli")
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
  skip_if_not_installed("deli")
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

test_that("ipw() refuses .by on a joint continuous fit", {
  skip_if_not_installed("deli")
  dat <- sim_joint_continuous()
  dat$grp <- factor(ifelse(dat$x1 > 0, "hi", "lo"), levels = c("lo", "hi"))
  # The modifier is interacted with both treatments, so the only thing wrong
  # with the request is the combination itself.
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
  skip_if_not_installed("deli")
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
