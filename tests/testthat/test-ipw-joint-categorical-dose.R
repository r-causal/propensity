# A joint intervention on a categorical treatment and a dose.
#
# The two-model route reports a dose from the marginal structural model's own
# coefficients, in a vocabulary written for a first treatment with two levels.
# This file is that vocabulary at more of them. A three-level treatment crossed
# with a dose reports one level contrast per non-reference level at a dose of
# zero, one dose slope at the reference arm, and one interaction row per
# non-reference level, which is the same three readings the binary crossing
# reports with the level contrasts run out over the levels there are:
#
#   effect   contrast          group            coefficient
#   diff     z: mid vs lo      e = 0            zmid
#   diff     z: hi vs lo       e = 0            zhi
#   slope    e: per unit       z = lo           e
#   diff     z: mid vs lo      e + 1 vs e       zmid:e
#   diff     z: hi vs lo       e + 1 vs e       zhi:e
#
# Two levels give back exactly the three rows the binary crossing reports, so
# the pinned binary suite is this surface's own degenerate case rather than a
# separate contract. Rows come out in the outcome design's column order, every
# row is one coefficient of the weighted fit, and a logit outcome relabels every
# row `log(or)` as it already does there: these are the marginal structural
# model's coefficients rather than contrasts of standardized means, and a
# coefficient of a logit model is a log odds ratio.
#
# Nothing else about the pairing is new. The stacked system carries the
# multinomial score for the categorical component and the linear score with a
# conditional variance for the dose, each written against its own component's
# type, and the stabilization slot resolves each numerator the same way. What
# this file pins there is that the two sides compose: a categorical numerator
# beside a dose numerator names its parameters for the component it belongs to
# and the system solves at the weights the outcome model was fit under.
#
# The one thing the vocabulary has no reading for is a dose entering the
# marginal structural model through a basis. That model is not linear in the
# dose, so no coefficient of it is a slope, and the surface reports it as this
# route already reports every other non-bare model: one row per
# treatment-reading coefficient, named by the coefficient, with no group column
# and no per-level rows invented for it.

# ---- data simulator ---------------------------------------------------------

# The level a unit takes, drawn from a row of category probabilities.
draw_categorical_levels <- function(probs, levels) {
  cumulative <- t(apply(probs, 1, cumsum))
  u <- stats::runif(nrow(probs))
  index <- 1L + rowSums(u > cumulative[, -ncol(cumulative), drop = FALSE])

  factor(levels[index], levels = levels)
}

# A three-level treatment assigned from the covariates and a dose that depends
# on it, which is the factorization the container asks for: the second model
# conditions on the first. The outcome carries a genuine interaction with each
# non-reference level, and the two interactions have opposite signs, so no
# reported row is zero and no two of them are the same number.
sim_joint_categorical_dose <- function(seed = 7601, n = 1200) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)

  eta_mid <- -0.1 + 0.5 * x1 + 0.3 * x2
  eta_hi <- -0.3 + 0.3 * x1 - 0.4 * x2
  denom <- 1 + exp(eta_mid) + exp(eta_hi)
  z <- draw_categorical_levels(
    cbind(1, exp(eta_mid), exp(eta_hi)) / denom,
    c("lo", "mid", "hi")
  )

  z_mid <- as.integer(z == "mid")
  z_hi <- as.integer(z == "hi")
  e <- 0.4 + 0.5 * x1 - 0.3 * x2 - 0.8 * z_mid + 0.5 * z_hi + rnorm(n)
  y <- 1 +
    0.7 * z_mid +
    0.3 * z_hi +
    0.5 * e +
    0.6 * z_mid * e -
    0.4 * z_hi * e +
    rnorm(n)
  yb <- rbinom(
    n,
    1,
    plogis(
      -0.4 +
        0.6 * z_mid +
        0.3 * z_hi +
        0.4 * e +
        0.5 * z_mid * e -
        0.3 * z_hi * e
    )
  )

  data.frame(x1, x2, z, e, y, yb)
}

# ---- model fitting ----------------------------------------------------------

# The two-model route over that frame. Every multinomial is fit inline against
# the frame this function holds, for the reason every multinomial in these files
# is: the fit keeps no model frame and the route rebuilds its design by
# re-evaluating the fitting call in the formula's own environment.
#
# The dose is always stabilized, which `wt_joint()` requires of a continuous
# component, and the categorical component's numerator is what moves between
# fixtures.
fit_joint_categorical_dose <- function(
  dat,
  z_rhs = c("x1", "x2"),
  e_rhs = c("z", "x1", "x2"),
  outcome_rhs = "z * e",
  outcome_family = c("gaussian", "binomial"),
  z_numerator = c("none", "marginal", "model"),
  dose_numerator = c("marginal", "model")
) {
  outcome_family <- match.arg(outcome_family)
  z_numerator <- match.arg(z_numerator)
  dose_numerator <- match.arg(dose_numerator)

  ps_z <- nnet::multinom(
    stats::reformulate(z_rhs, response = "z"),
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_e <- lm(stats::reformulate(e_rhs, response = "e"), data = dat)

  # A proper subset of each propensity score model's design, so the block a
  # numerator contributes is not the denominator's block under another name.
  num_z <- nnet::multinom(
    z ~ x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  num_e <- lm(e ~ x2, data = dat)

  # A multinomial over two levels is a logistic regression solved by another
  # optimizer, and every part of the package says so: `wt_ate()` refuses to
  # weight it as a categorical exposure and `joint_wt_models()` types it binary.
  # The fixture follows, so a two-level fit takes the binary weight the route
  # will read it under.
  stabilize <- switch(
    z_numerator,
    none = NULL,
    marginal = TRUE,
    model = num_z
  )
  two_level <- length(ps_z$lev) == 2L
  z_type <- if (two_level) "binary" else "categorical"

  w_z <- withr::with_options(
    list(propensity.quiet = TRUE),
    if (two_level) {
      wt_ate(
        ps_z,
        stabilize = if (identical(z_numerator, "none")) NULL else TRUE
      )
    } else {
      probabilities <- unname(stats::predict(ps_z, type = "probs"))
      colnames(probabilities) <- ps_z$lev

      wt_ate(
        probabilities,
        dat$z,
        exposure_type = "categorical",
        stabilize = stabilize
      )
    }
  )
  w_e <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(stats::fitted(ps_e)),
      dat$e,
      exposure_type = "continuous",
      stabilize = if (identical(dose_numerator, "model")) num_e else TRUE
    )
  )
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(w_z, w_e, exposure_type = c(z_type, "continuous"))
  )

  outcome_var <- if (identical(outcome_family, "binomial")) "yb" else "y"
  fmla <- stats::reformulate(outcome_rhs, response = outcome_var)
  outcome_mod <- if (identical(outcome_family, "binomial")) {
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
    models = joint_wt_models(z = ps_z, e = ps_e),
    ps_z = ps_z,
    ps_e = ps_e,
    num_z = num_z,
    num_e = num_e,
    w_z = w_z,
    w_e = w_e,
    outcome_mod = outcome_mod,
    wts = wts
  )
}

# ---- the reported surface ---------------------------------------------------

# The identity columns, written out against the conventions in the header rather
# than imported, so a change to the labeling shows up as a disagreement with a
# written-out contract.
joint_categorical_dose_expected_rows <- function(
  scale = c("diff", "diff", "slope", "diff", "diff")
) {
  data.frame(
    effect = scale,
    contrast = c(
      "z: mid vs lo",
      "z: hi vs lo",
      "e: per unit",
      "z: mid vs lo",
      "z: hi vs lo"
    ),
    group = c("e = 0", "e = 0", "z = lo", "e + 1 vs e", "e + 1 vs e"),
    stringsAsFactors = FALSE
  )
}

# The coefficients those rows name, in the outcome design's column order, which
# is the order the rows come out in.
joint_categorical_dose_coefficients <- c("zmid", "zhi", "e", "zmid:e", "zhi:e")

# Every reported row is one coefficient of the weighted marginal structural
# model, keyed by the coefficient it names rather than by position.
expect_joint_categorical_dose_estimates <- function(
  res,
  outcome_mod,
  coefs,
  tolerance = 1e-8
) {
  beta <- coef(outcome_mod)

  testthat::expect_equal(
    res$estimates$estimate,
    unname(beta[coefs]),
    tolerance = tolerance
  )
  testthat::expect_true(all(is.finite(res$estimates$std.err)))
  testthat::expect_true(all(res$estimates$std.err > 0))

  invisible(res)
}

# ---- the type gate ----------------------------------------------------------

test_that("a dose is admitted beside a categorical first treatment", {
  # The vocabulary a dose reports is written per non-reference level of the
  # first treatment, so a treatment with more than two of them has rows there
  # and the pair is no longer refused by type.
  expect_silent(
    check_ipw_joint_models_types(c(z = "categorical", e = "continuous"))
  )
  expect_silent(
    check_ipw_joint_models_types(c(a = "binary", e = "continuous"))
  )
})

test_that("a dose is still refused in the first slot and beside another dose", {
  # What remains refused is the position rather than the company: nothing here
  # carries the density of a first factor of the factorization, so a dose may
  # only be the second treatment whatever the first one is.
  cnd <- tryCatch(
    check_ipw_joint_models_types(c(e = "continuous", z = "categorical")),
    error = identity
  )
  expect_s3_class(cnd, "propensity_ipw_exposure_error")
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`e`", fixed = TRUE)
  expect_match(message, "first", fixed = TRUE)

  # The clause the retired refusal was written on. A first treatment with more
  # than two levels is reported now, so a message still saying so would describe
  # the package the tests below contradict.
  expect_false(grepl("only beside a binary first one", message, fixed = TRUE))

  cnd_pair <- tryCatch(
    check_ipw_joint_models_types(c(e1 = "continuous", e2 = "continuous")),
    error = identity
  )
  expect_s3_class(cnd_pair, "propensity_ipw_exposure_error")
})

test_that("ipw() estimates a categorical treatment crossed with a dose", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_dose()
  fx <- fit_joint_categorical_dose(dat)

  expect_identical(
    fx$models$exposure_type,
    c(z = "categorical", e = "continuous")
  )

  res <- ipw(fx$models, fx$outcome_mod)
  expect_s3_class(res, "ipw")
})

# ---- the vocabulary surface at three levels ---------------------------------

test_that("a three-level treatment crossed with a dose reports five rows", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_dose()
  fx <- fit_joint_categorical_dose(dat)

  # The seed rebuilds the weights the outcome model was fit under and every
  # equation sits at its root there, which is what the preflight compares before
  # anything is solved.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  layout <- ipw_theta_layout(spec)
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(fx$wts),
    tolerance = 1e-10
  )
  mat <- build_ipw_psi(spec, layout)(layout$init)
  expect_false(anyNA(mat))
  expect_true(all(abs(rowSums(mat)) / spec$n < 1e-8))

  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates
  expected <- joint_categorical_dose_expected_rows()

  expect_identical(nrow(est), 5L)
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)

  # A dose has no cells to standardize over, so the rows are the weighted fit's
  # own coefficients, read off the outcome block of theta.
  expect_joint_categorical_dose_estimates(
    res,
    fx$outcome_mod,
    joint_categorical_dose_coefficients
  )

  # Each row names a distinct quantity, which is what the accessors key by.
  labels <- paste(est$effect, est$contrast, est$group)
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(names(coef(res)), labels)
  expect_identical(dimnames(vcov(res)), list(labels, labels))
})

test_that("an additive three-level crossing reports three ungrouped rows", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_dose()
  fx <- fit_joint_categorical_dose(dat, outcome_rhs = "z + e")

  # With no interaction the model forces one contrast per non-reference level at
  # every dose and one slope in every arm, so no row is evaluated anywhere in
  # particular. `overall` is a group claim rather than the absence of one, which
  # is what separates this table from a coefficient-surface table of the same
  # width.
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates

  expect_identical(nrow(est), 3L)
  expect_identical(est$effect, c("diff", "diff", "slope"))
  expect_identical(
    est$contrast,
    c("z: mid vs lo", "z: hi vs lo", "e: per unit")
  )
  expect_identical(est$group, rep("overall", 3L))
  expect_joint_categorical_dose_estimates(
    res,
    fx$outcome_mod,
    c("zmid", "zhi", "e")
  )
})

test_that("a logit marginal structural model reports the five rows on the odds scale", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_dose()
  fx <- fit_joint_categorical_dose(dat, outcome_family = "binomial")

  # The scale follows the outcome link and the row's identity does not: a
  # coefficient of a logit marginal structural model is a log odds ratio,
  # admitted here where the standardized cell surface refuses it, because these
  # rows are the model's own coefficients.
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates
  expected <- joint_categorical_dose_expected_rows(scale = rep("log(or)", 5))

  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)
  expect_joint_categorical_dose_estimates(
    res,
    fx$outcome_mod,
    joint_categorical_dose_coefficients
  )
})

test_that("a treatment column coded some other way is refused at three levels", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_dose()
  ordered_dat <- dat
  ordered_dat$z <- factor(
    as.character(dat$z),
    levels = c("lo", "mid", "hi"),
    ordered = TRUE
  )

  # The rows name the coefficients of a model whose treatment columns are the
  # indicators of the non-reference levels. Polynomial contrasts leave the term
  # label untouched and rescale those columns, so the fit runs to completion and
  # reports numbers that are not the effects the rows name.
  fx <- fit_joint_categorical_dose(dat)
  bad <- lm(y ~ z * e, data = ordered_dat, weights = fx$wts)

  cnd <- tryCatch(ipw(fx$models, bad), error = identity)
  expect_s3_class(cnd, "propensity_ipw_msm_error")
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(cnd)),
    "z",
    fixed = TRUE
  )
})

# ---- two levels, which is the binary surface --------------------------------

test_that("a two-level treatment beside a dose keeps the three-row surface", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_dose()
  two_level <- dat[dat$z != "hi", ]
  two_level$z <- factor(as.character(two_level$z), levels = c("lo", "mid"))

  fx <- fit_joint_categorical_dose(two_level)

  # At two levels the per-level loop emits precisely the rows the binary
  # crossing reports: one contrast at a dose of zero, one slope at the reference
  # arm, and one interaction row. The vocabulary is the same one; what changed
  # is only how many non-reference levels it runs out over.
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates

  expect_identical(nrow(est), 3L)
  expect_identical(est$effect, c("diff", "slope", "diff"))
  expect_identical(
    est$contrast,
    c("z: mid vs lo", "e: per unit", "z: mid vs lo")
  )
  expect_identical(est$group, c("e = 0", "z = lo", "e + 1 vs e"))
  expect_joint_categorical_dose_estimates(
    res,
    fx$outcome_mod,
    c("zmid", "e", "zmid:e")
  )
})

# ---- a dose entering through a basis ----------------------------------------

test_that("a basis dose beside a categorical treatment reports coefficients", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_dose()
  fx <- fit_joint_categorical_dose(
    dat,
    outcome_rhs = "z * splines::ns(e, 3)"
  )

  # A model carrying a basis in the dose is not linear in it, so no coefficient
  # of it is a slope and the vocabulary has no reading for it. This route
  # reports such a model on its coefficient surface, which is what it already
  # does for a basis dose beside a binary treatment: one row per
  # treatment-reading coefficient, named by the coefficient, with no group
  # column and no per-level rows fabricated for a shape that has none.
  res <- ipw(fx$models, fx$outcome_mod)
  est <- res$estimates
  columns <- colnames(model.matrix(fx$outcome_mod))[-1]

  expect_identical(est$effect, rep("coef", length(columns)))
  expect_identical(est$contrast, columns)
  expect_null(est[["group"]])
  expect_true(all(is.finite(est$std.err)))

  # The conditional-reading convention the single-dose route applies to a basis
  # exposure is that route's own. This surface reports both readings and
  # announces nothing, so a caller naming the marginal one is answered rather
  # than refused.
  messages <- testthat::capture_messages(ipw(fx$models, fx$outcome_mod))
  expect_false(any(grepl("marginaleffects", messages, fixed = TRUE)))
  expect_identical(res$effects, "marginal")
  expect_s3_class(
    ipw(fx$models, fx$outcome_mod, effects = "marginal"),
    "ipw"
  )
})

# ---- the stabilization slot -------------------------------------------------

test_that("a categorical component's marginal numerator composes with the dose's", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_dose()
  fx <- fit_joint_categorical_dose(dat, z_numerator = "marginal")

  # Each component's numerator is resolved by its own type: the categorical one
  # is the k - 1 free level proportions and the dose's is its two marginal
  # moments, and the slice each takes is as wide as its own numerator.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  expect_identical(spec$stab$widths, c(2L, 2L))

  seed <- ipw_init_joint_models_stab(spec)
  expect_identical(
    names(seed),
    c("stab_z_mid", "stab_z_hi", "stab_e_mu", "stab_e_sigma2")
  )

  layout <- ipw_theta_layout(spec)
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(fx$wts),
    tolerance = 1e-10
  )

  res <- ipw(fx$models, fx$outcome_mod)
  theta <- coef(res$fit)
  expect_identical(
    names(theta)[grepl("^stab_", names(theta))],
    c("stab_z_mid", "stab_z_hi", "stab_e_mu", "stab_e_sigma2")
  )
  expect_equal(
    unname(theta[c("stab_z_mid", "stab_z_hi")]),
    c(mean(dat$z == "mid"), mean(dat$z == "hi")),
    tolerance = 1e-8
  )
  expect_joint_categorical_dose_estimates(
    res,
    fx$outcome_mod,
    joint_categorical_dose_coefficients
  )
})

test_that("both components' numerator models are blocks of the same system", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_dose()
  fx <- fit_joint_categorical_dose(
    dat,
    outcome_rhs = "z * e + x2",
    z_numerator = "model",
    dose_numerator = "model"
  )

  # A multinomial numerator's parameters are named level-major and term-minor
  # under the joint route's own prefix, and a dose numerator's are its
  # coefficients followed by the spread its density is read at. Two stabilized
  # components mean a name saying only which role a parameter plays would not
  # say whose, which is what the prefix is for.
  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  expect_identical(
    spec$stab$widths,
    c(length(as.vector(coef(fx$num_z))), length(coef(fx$num_e)) + 1L)
  )

  res <- ipw(fx$models, fx$outcome_mod)
  theta <- coef(res$fit)
  stab <- names(theta)[grepl("^stab_", names(theta))]

  expect_identical(
    stab,
    c(
      "stab_z_mid:(Intercept)",
      "stab_z_mid:x2",
      "stab_z_hi:(Intercept)",
      "stab_z_hi:x2",
      "stab_e_(Intercept)",
      "stab_e_x2",
      "stab_e_sigma2"
    )
  )
  expect_equal(
    unname(theta[c("stab_e_(Intercept)", "stab_e_x2")]),
    unname(coef(fx$num_e)),
    tolerance = 1e-6
  )
  # The categorical numerator's multinom stops short of its score's root, so
  # the solved weights sit about 2e-7 from the caller's, the same convergence
  # floor the stab pin above absorbs at 1e-6.
  expect_joint_categorical_dose_estimates(
    res,
    fx$outcome_mod,
    joint_categorical_dose_coefficients,
    tolerance = 1e-6
  )
})
