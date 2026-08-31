# Tests for ipw() with a continuous exposure: an lm (or gaussian-family glm)
# propensity score model of the exposure on covariates, stabilized continuous
# ATE weights, and a weighted marginal structural outcome model with exactly one
# exposure term. These exercise ipw() end to end the way a user would call it,
# pinning estimand detection, the estimates-table contract, point-estimate parity
# with the MSM coefficient, a bootstrap standard-error oracle, the stabilization
# variants, the estimand and MSM guards, the gaussian-glm routing, the fit
# object, and the print and as.data.frame surfaces.

# ---- data simulator ---------------------------------------------------------

sim_continuous <- function(seed = 2024, n = 800) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  yc <- 1 + 0.6 * A + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  yb <- rbinom(n, 1, plogis(-0.5 + 0.4 * A + 0.3 * x1 - 0.2 * x2))
  data.frame(x1, x2, A, yc, yb)
}

# The same simulation with an exposure that is positive everywhere, which a
# propensity model with a log link needs: the link takes the log of the
# response to start its iteration, and a dose that can be zero or negative has
# no such start. The dose is the one above shifted up rather than exponentiated,
# so it stays conditionally normal with a single spread and the density the
# weights divide by is the density the dose has. The columns are named as they
# are in `sim_continuous()` so that every formula the fixtures below write reads
# either simulation.
sim_continuous_positive <- function(seed = 2024, n = 800) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  yc <- 1 + 0.6 * A + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  yb <- rbinom(n, 1, plogis(-3 + 0.3 * A + 0.3 * x1 - 0.2 * x2))
  data.frame(x1, x2, A, yc, yb)
}

# ---- model fitting ----------------------------------------------------------

# Continuous ATE weights from a numeric fitted propensity, always computed
# silently. Stabilized weights are the recommended default for a continuous
# exposure; stabilize = FALSE emits an alert unless quieted. `.density`,
# `numerator`, and `.sigma` are the three choices a set of continuous weights
# records, and every one of them changes what `ipw()` has to rebuild.
continuous_weights <- function(
  fitted_ps,
  A,
  stabilize = TRUE,
  stab_score = NULL,
  .density = "normal",
  numerator = "marginal",
  .sigma = NULL
) {
  withr::with_options(
    list(propensity.quiet = TRUE),
    # A bootstrap resample of a fixture near the finite-variance boundary lands
    # on either side of it, so the report is muffled here rather than left to
    # arrive on whichever resamples happen to cross.
    muffle_variance_warning(wt_ate(
      fitted_ps,
      A,
      .sigma = .sigma,
      exposure_type = "continuous",
      stabilize = stabilize,
      stabilization_score = stab_score,
      .density = .density,
      numerator = numerator
    ))
  )
}

# Fit the propensity score model of A on the covariates, build continuous ATE
# weights from its fitted values, and fit a weighted MSM. `ps_type` selects the
# lm, the gaussian-family glm, the log-link gaussian glm, or the robust Huber
# form of the propensity model; the first two share fitted values and so produce
# identical weights, and the log-link one needs the positive exposure
# `sim_continuous_positive()` simulates. `msm_rhs` allows a multiple-term
# right-hand side for the MSM guard. `.density`, `numerator`, and `.sigma` are
# passed through to the weights. The weights are kept as a psw object so the
# estimand survives into the outcome model frame for detection. The
# quasibinomial MSM tightens its IRLS tolerance so its coefficients sit at the
# weighted MLE to well below the point-estimate comparison tolerance, and the
# log-link and robust propensity models tighten their own for the same reason:
# the stacked score is solved to a tighter root than the iteration's default
# stops at. `rlm` stops on the relative change in its coefficients, so `acc`
# rather than `epsilon` is what tightens it.
fit_continuous_models <- function(
  dat,
  ps_type = c("lm", "glm", "glm_log", "rlm"),
  outcome_family = c("gaussian", "binomial"),
  stabilize = TRUE,
  stab_score = NULL,
  msm_rhs = "A",
  .density = "normal",
  numerator = "marginal",
  .sigma = NULL
) {
  ps_type <- match.arg(ps_type)
  outcome_family <- match.arg(outcome_family)

  ps_mod <- switch(
    ps_type,
    lm = lm(A ~ x1 + x2, data = dat),
    glm = glm(A ~ x1 + x2, data = dat, family = gaussian()),
    glm_log = glm(
      A ~ x1 + x2,
      data = dat,
      family = gaussian(link = "log"),
      control = glm.control(epsilon = 1e-14, maxit = 200)
    ),
    rlm = MASS::rlm(A ~ x1 + x2, data = dat, acc = 1e-10)
  )
  fitted_ps <- as.double(fitted(ps_mod))
  wts <- continuous_weights(
    fitted_ps,
    dat$A,
    stabilize = stabilize,
    stab_score = stab_score,
    .density = .density,
    numerator = numerator,
    .sigma = .sigma
  )

  outcome_var <- if (outcome_family == "binomial") "yb" else "yc"
  msm_fmla <- stats::reformulate(msm_rhs, response = outcome_var)
  if (outcome_family == "binomial") {
    outcome_mod <- glm(
      msm_fmla,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  } else {
    outcome_mod <- lm(msm_fmla, data = dat, weights = wts)
  }

  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

estimates_columns <- c(
  "effect",
  "estimate",
  "std.err",
  "z",
  "ci.lower",
  "ci.upper",
  "conf.level",
  "p.value"
)

# ---- end-to-end -------------------------------------------------------------

test_that("ipw() runs continuous ate end to end and auto-detects the estimand", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_s3_class(res, "ipw")
  expect_equal(res$se_method, "mestimation")
  expect_equal(res$estimand, "ate")
  expect_false(is.null(res$fit))
})

# ---- point-estimate parity with the MSM exposure coefficient ----------------

test_that("ipw() continuous point estimate equals the MSM exposure coefficient", {
  dat <- sim_continuous()

  # identity link (lm) is exact
  mods_g <- fit_continuous_models(dat, outcome_family = "gaussian")
  res_g <- ipw(mods_g$ps_mod, mods_g$outcome_mod)
  expect_equal(nrow(res_g$estimates), 1L)
  expect_equal(
    res_g$estimates$estimate,
    unname(coef(mods_g$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )

  # logit link (quasibinomial) tightened to its weighted MLE
  mods_b <- fit_continuous_models(dat, outcome_family = "binomial")
  res_b <- ipw(mods_b$ps_mod, mods_b$outcome_mod)
  expect_equal(nrow(res_b$estimates), 1L)
  expect_equal(
    res_b$estimates$estimate,
    unname(coef(mods_b$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )
})

# ---- estimates table shape and effect labels --------------------------------

test_that("ipw() continuous labels the effect by the outcome link", {
  dat <- sim_continuous()

  # identity link: a single slope, the eight-column contract, no contrast
  mods_g <- fit_continuous_models(dat, outcome_family = "gaussian")
  res_g <- ipw(mods_g$ps_mod, mods_g$outcome_mod)
  expect_named(res_g$estimates, estimates_columns)
  expect_equal(nrow(res_g$estimates), 1L)
  expect_equal(res_g$estimates$effect, "slope")
  expect_true(all(res_g$estimates$conf.level == 0.95))

  # logit link: a single log odds ratio
  mods_b <- fit_continuous_models(dat, outcome_family = "binomial")
  res_b <- ipw(mods_b$ps_mod, mods_b$outcome_mod)
  expect_named(res_b$estimates, estimates_columns)
  expect_equal(nrow(res_b$estimates), 1L)
  expect_equal(res_b$estimates$effect, "log(or)")
})

# ---- standard-error oracle: nonparametric bootstrap -------------------------

# The M-estimation sandwich standard error, which accounts for estimating the
# propensity model, is validated against a seed-pinned nonparametric bootstrap
# that refits the propensity model, the weights, and the MSM on each resample.
# The two are asymptotically equivalent; a 15 percent relative band accommodates
# the finite-sample and Monte Carlo gap while keeping this a meaningful check.
test_that("ipw() continuous mestimation SE tracks a nonparametric bootstrap", {
  skip_on_cran()

  dat <- sim_continuous(seed = 2024, n = 600)
  mods <- fit_continuous_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  mest_se <- res$estimates$std.err

  boot_slope <- function(d) {
    ps <- lm(A ~ x1 + x2, data = d)
    fitted_ps <- as.double(fitted(ps))
    w <- continuous_weights(fitted_ps, d$A, stabilize = TRUE)
    msm <- lm(yc ~ A, data = d, weights = as.double(w))
    unname(coef(msm)[["A"]])
  }

  withr::local_seed(918)
  reps <- 400L
  n <- nrow(dat)
  boot <- vapply(
    seq_len(reps),
    function(i) {
      boot_slope(dat[sample.int(n, n, replace = TRUE), , drop = FALSE])
    },
    numeric(1)
  )
  boot_se <- stats::sd(boot)

  rel <- abs(mest_se - boot_se) / boot_se
  expect_lt(rel, 0.15)
})

# ---- stabilization variants -------------------------------------------------

test_that("ipw() continuous runs stabilized, unstabilized, and fixed-score weights", {
  dat <- sim_continuous()

  # stabilized: the recommended default for a continuous exposure
  mods_s <- fit_continuous_models(dat, stabilize = TRUE)
  res_s <- ipw(mods_s$ps_mod, mods_s$outcome_mod)
  expect_equal(res_s$estimand, "ate")
  expect_true(is.finite(res_s$estimates$estimate))
  expect_true(
    is.finite(res_s$estimates$std.err) && res_s$estimates$std.err > 0
  )

  # unstabilized: the ipw() call is silent and the SE is finite (typically
  # larger than the stabilized SE)
  mods_u <- fit_continuous_models(dat, stabilize = FALSE)
  res_u <- expect_silent(ipw(mods_u$ps_mod, mods_u$outcome_mod))
  expect_true(is.finite(res_u$estimates$estimate))
  expect_true(
    is.finite(res_u$estimates$std.err) && res_u$estimates$std.err > 0
  )

  # user-supplied stabilization score
  mods_f <- fit_continuous_models(dat, stabilize = TRUE, stab_score = 0.5)
  res_f <- ipw(mods_f$ps_mod, mods_f$outcome_mod)
  expect_true(is.finite(res_f$estimates$estimate))
  expect_true(
    is.finite(res_f$estimates$std.err) && res_f$estimates$std.err > 0
  )
})

# ---- a numerator model in the stack -----------------------------------------

# The three fits this section compares. All three weight the same propensity
# score model; what separates them is the numerator. `model` stabilizes on a
# fitted model of the dose on a baseline covariate, `score` stabilizes on the
# very same numerator handed over as a vector of numbers, and `marginal` takes
# the default. The first two build one and the same set of weights, so anything
# that differs between them is the sandwich reading the numerator differently
# rather than the weights being different.
continuous_numerator_fits <- function(dat) {
  num_mod <- lm(A ~ x1, data = dat)
  score <- stats::dnorm(
    dat$A,
    as.numeric(fitted(num_mod)),
    sqrt(mean(residuals(num_mod)^2))
  )

  list(
    num_mod = num_mod,
    score = score,
    model = fit_continuous_models(dat, stabilize = num_mod),
    score_fit = fit_continuous_models(
      dat,
      stabilize = TRUE,
      stab_score = score
    ),
    marginal = fit_continuous_models(dat)
  )
}

test_that("ipw() reports the same estimate from a numerator model and its score", {
  dat <- sim_continuous()
  fits <- continuous_numerator_fits(dat)

  # The two routes weight the units identically, which is what makes the
  # standard errors below comparable at all.
  expect_equal(
    as.numeric(fits$model$wts),
    as.numeric(fits$score_fit$wts),
    tolerance = 1e-12
  )

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)

  expect_equal(
    res_model$estimates$estimate,
    unname(coef(fits$model$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )
  expect_equal(
    res_model$estimates$estimate,
    res_score$estimates$estimate,
    tolerance = 1e-8
  )
})

test_that("ipw() stacks the numerator model it was given", {
  dat <- sim_continuous()
  fits <- continuous_numerator_fits(dat)

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)

  # A supplied score is a constant of the stacked system and a supplied model is
  # not: the model's coefficients and the spread its density is read at are
  # parameters solved alongside everything else.
  theta_model <- coef(res_model$fit)
  theta_score <- coef(res_score$fit)

  expect_equal(
    length(theta_model),
    length(theta_score) + length(coef(fits$num_mod)) + 1L
  )
  expect_equal(
    names(theta_model),
    c(
      names(coef(fits$model$ps_mod)),
      "sigma2_d",
      paste0("stab_", names(coef(fits$num_mod))),
      "sigma2_n",
      names(coef(fits$model$outcome_mod))
    )
  )

  # The block is solved at the model it was given rather than carried at
  # whatever the seed was.
  stab_block <- theta_model[grepl("^stab_", names(theta_model))]
  expect_equal(
    unname(stab_block),
    unname(coef(fits$num_mod)),
    tolerance = 1e-6
  )
  expect_equal(
    unname(theta_model[["sigma2_n"]]),
    mean(residuals(fits$num_mod)^2),
    tolerance = 1e-6
  )
})

test_that("a stacked numerator model reports a different standard error", {
  dat <- sim_continuous()
  fits <- continuous_numerator_fits(dat)

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)
  res_marginal <- ipw(fits$marginal$ps_mod, fits$marginal$outcome_mod)

  se_model <- res_model$estimates$std.err
  se_score <- res_score$estimates$std.err
  se_marginal <- res_marginal$estimates$std.err

  expect_true(is.finite(se_model) && se_model > 0)

  # The comparison against the score fit is only about the numerator while the
  # weights are the same weights, which is what this restates.
  expect_equal(
    as.numeric(fits$model$wts),
    as.numeric(fits$score_fit$wts),
    tolerance = 1e-12
  )

  # The score fit is the same weights with the numerator held fixed, so the
  # difference between the two standard errors is exactly the uncertainty of
  # having fit the numerator model.
  expect_false(isTRUE(all.equal(se_model, se_score, tolerance = 1e-6)))

  # The marginal stabilizer is a different numerator over the same denominator,
  # so it is a different set of weights and a different standard error again.
  expect_false(isTRUE(all.equal(se_model, se_marginal, tolerance = 1e-6)))
})

# ---- estimand support -------------------------------------------------------

test_that("ipw() continuous rejects estimands other than ate", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  # a valid non-ate estimand is unsupported for a continuous exposure
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, estimand = "att"),
    class = "propensity_ipw_estimand_error"
  )

  # an entirely invalid estimand string errors with the same class
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, estimand = "nonsense"),
    class = "propensity_ipw_estimand_error"
  )
})

# ---- MSM guard --------------------------------------------------------------

test_that("ipw() continuous rejects an MSM term that reads a covariate", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, msm_rhs = c("A", "A:x1"))

  # the requirement is on what a term reads rather than on how many terms there
  # are: a term reading a covariate contributes a coefficient that depends on
  # that covariate, so no row of the table names an effect. The error must
  # direct the user to the returned fit object for the full coefficient vector
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_msm_error",
    regexp = "fit"
  )

  # a curve in the exposure alone is admitted, and the stacked system estimates
  # one row per coefficient; the full contract for those rows is pinned in
  # test-ipw-joint-continuous.R
  curve <- fit_continuous_models(dat, msm_rhs = c("A", "I(A^2)"))
  res <- ipw(curve$ps_mod, curve$outcome_mod)
  expect_identical(nrow(res$estimates), 2L)
  expect_identical(res$estimates$contrast, c("A", "I(A^2)"))
})

# ---- the reading a curve reports --------------------------------------------

# A curve in the exposure is admitted, and what it reports is the outcome
# model's coefficient surface. No coefficient of a curve is the effect of the
# dose, since the response has a different slope at every dose, so there is no
# marginal reading of such a fit to present and `ipw()` records the conditional
# one instead of picking a row to call an effect. The boundary is how many
# design columns the exposure enters through, which is the same boundary that
# decides whether the stacked system's rows are named after their coefficients.

test_that("ipw() continuous reports the conditional reading for a curve", {
  dat <- sim_continuous()
  curve <- fit_continuous_models(dat, msm_rhs = c("A", "I(A^2)"))
  res <- ipw(curve$ps_mod, curve$outcome_mod)

  expect_identical(res$effects, "conditional")
  expect_identical(res$readings, "conditional")
  expect_identical(class(res), "ipw")

  # The reading is the outcome model's own coefficient vector, intercept
  # included, under the covariance the stacked system leaves for it.
  expect_identical(coef(res), coef(curve$outcome_mod))
  expect_identical(names(coef(res)), c("(Intercept)", "A", "I(A^2)"))
  expect_identical(vcov(res), vcov(res$outcome_mod))
  expect_identical(tidy(res)$term, names(coef(curve$outcome_mod)))
})

test_that("ipw() continuous refuses the marginal reading of a curve", {
  dat <- sim_continuous()
  curve <- fit_continuous_models(dat, msm_rhs = c("A", "I(A^2)"))

  expect_error(
    ipw(curve$ps_mod, curve$outcome_mod, effects = "marginal"),
    class = "propensity_ipw_effects_error"
  )

  # The result refuses it too, from the set of readings it declares, so a caller
  # who takes the conditional result and asks it for the other reading is told
  # the same thing the constructor told them.
  res <- ipw(curve$ps_mod, curve$outcome_mod)
  expect_error(
    causalgenerics::as_marginal(res),
    class = "causalgenerics_unsupported_reading_marginal"
  )
  expect_error(
    coef(res, effects = "marginal"),
    class = "causalgenerics_unsupported_reading_marginal"
  )
})

test_that("ipw() continuous leaves a single exposure term marginal", {
  dat <- sim_continuous()

  # One design column is one slope, and that slope is the dose response
  # everywhere, so the reading, the class, and the row are the ones this path
  # has always reported. A covariate beside the exposure changes none of it.
  for (rhs in list("A", c("A", "x1"))) {
    mods <- fit_continuous_models(dat, msm_rhs = rhs)
    res <- withr::with_options(
      list(propensity.quiet = FALSE),
      expect_no_message(ipw(mods$ps_mod, mods$outcome_mod))
    )

    expect_identical(res$effects, "marginal")
    expect_identical(res$readings, c("marginal", "conditional"))
    expect_identical(class(res), "ipw")
    expect_identical(res$estimates$effect, "slope")
    expect_identical(causalgenerics::as_marginal(res), res)
  }
})

# ---- gaussian-glm routing ---------------------------------------------------

test_that("ipw() routes a gaussian-family glm ps model identically to lm", {
  dat <- sim_continuous()

  ps_lm <- lm(A ~ x1 + x2, data = dat)
  ps_glm <- glm(A ~ x1 + x2, data = dat, family = gaussian())

  # The two propensity fits give identical coefficients but not bitwise
  # identical fitted values. Building the weights once, from a single fitted
  # vector, and weighting one outcome model with them keeps that difference out
  # of everything downstream, so what the two calls are compared on is the
  # routing of the propensity model itself rather than the last bits of two
  # separately fitted weight vectors amplified through the sandwich.
  fitted_ps <- as.double(fitted(ps_lm))
  wts <- continuous_weights(fitted_ps, dat$A)
  outcome_mod <- lm(yc ~ A, data = dat, weights = wts)

  res_lm <- ipw(ps_lm, outcome_mod)
  res_glm <- ipw(ps_glm, outcome_mod)

  expect_equal(res_glm$estimand, res_lm$estimand)
  expect_equal(
    res_glm$estimates$estimate,
    res_lm$estimates$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    res_glm$estimates$std.err,
    res_lm$estimates$std.err,
    tolerance = 1e-8
  )
})

# ---- fit object -------------------------------------------------------------

test_that("ipw() continuous fit theta matches the reported estimate", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  # the MSM exposure coefficient carries the outcome-model term name, unique
  # across theta since the propensity block uses the covariate names
  tidied <- generics::tidy(res$fit)
  a_rows <- tidied[tidied$term == "A", ]
  expect_equal(nrow(a_rows), 1L)
  expect_equal(a_rows$estimate, res$estimates$estimate, tolerance = 1e-8)
})

# ---- print and as.data.frame ------------------------------------------------

test_that("ipw() continuous print output is stable", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_snapshot(print(res))
})

test_that("as.data.frame(exponentiate = TRUE) relabels the continuous log odds ratio", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, outcome_family = "binomial")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  df <- as.data.frame(res, exponentiate = TRUE)
  expect_equal(nrow(df), 1L)
  expect_named(
    df,
    c("term", "estimate", "std.error", "statistic", "p.value")
  )
  expect_equal(df$term, "or")
})

# ---- offset guard -----------------------------------------------------------

test_that("ipw() continuous rejects an outcome model with an offset term", {
  dat <- sim_continuous()
  dat$off <- 0.3 * dat$x1

  ps_mod <- lm(A ~ x1 + x2, data = dat)
  fitted_ps <- as.double(fitted(ps_mod))
  wts <- continuous_weights(fitted_ps, dat$A, stabilize = TRUE)
  outcome_mod <- lm(yc ~ A + offset(off), data = dat, weights = wts)

  # The mestimation psi does not thread the offset into the MSM score block, so
  # the default path must reject an offset before any estimation and with no
  # message or other condition signaled first.
  expect_no_message(
    expect_error(
      ipw(ps_mod, outcome_mod),
      class = "propensity_ipw_offset_error",
      regexp = "offset"
    )
  )
})

# ---- linearization is unsupported for continuous ----------------------------

test_that("ipw() rejects linearization for a continuous exposure", {
  dat <- sim_continuous()

  # the lm propensity route
  mods_lm <- fit_continuous_models(dat, ps_type = "lm")
  expect_error(
    ipw(mods_lm$ps_mod, mods_lm$outcome_mod, se_method = "linearization"),
    class = "propensity_method_error"
  )

  # the gaussian-family glm propensity route rejects it identically
  mods_glm <- fit_continuous_models(dat, ps_type = "glm")
  expect_error(
    ipw(mods_glm$ps_mod, mods_glm$outcome_mod, se_method = "linearization"),
    class = "propensity_method_error"
  )
})

test_that("ipw() rejects robust standard errors for a continuous exposure", {
  dat <- sim_continuous()

  # A continuous exposure reports the marginal structural model's own
  # coefficients rather than counterfactual means of two cells, so there are no
  # Hajek means for a weights-fixed sandwich to describe. The refusal names the
  # method that was asked for.
  mods_lm <- fit_continuous_models(dat, ps_type = "lm")
  expect_error(
    ipw(mods_lm$ps_mod, mods_lm$outcome_mod, se_method = "robust"),
    class = "propensity_method_error",
    regexp = "robust"
  )

  mods_glm <- fit_continuous_models(dat, ps_type = "glm")
  expect_error(
    ipw(mods_glm$ps_mod, mods_glm$outcome_mod, se_method = "robust"),
    class = "propensity_method_error",
    regexp = "robust"
  )

  # One wording for both routes, pinned once.
  expect_propensity_error(
    ipw(mods_lm$ps_mod, mods_lm$outcome_mod, se_method = "robust")
  )
})

# ---- outcome-family validation ----------------------------------------------

test_that("ipw_spec_continuous rejects an unsupported outcome family", {
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)

  withr::local_seed(1)
  dat$ycount <- rpois(nrow(dat), exp(0.1 + 0.2 * dat$A))
  outcome_mod <- suppressWarnings(
    glm(ycount ~ A, data = dat, family = poisson(), weights = wts)
  )

  # Without the family guard a poisson marginal structural model is stacked as a
  # binomial score.
  expect_error(
    ipw_spec_continuous(ps_mod, outcome_mod),
    class = "propensity_ipw_family_error"
  )
})

# ---- continuous-path link validation ----------------------------------------
#
# The marginal structural model is rejected at entry with a classed error naming
# the supported links. Left to its downstream failure it misleads: a probit
# marginal structural model passes the family check but has no continuous effect
# label, so it errors late with a terse internal message. (A gaussian-identity
# glm ps model, already covered by the gaussian-glm routing test above, and a
# logit msm, covered by the effect-label test above, must keep working; the
# log-link msm below adds the one supported msm link not otherwise exercised.
# The propensity model's own link is the registry's business, and the log link
# it supports is pinned with the rest of the registry below.)

test_that("ipw() rejects a non-identity link on the continuous marginal structural model", {
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  msm <- suppressWarnings(
    glm(yb ~ A, data = dat, family = binomial(link = "probit"), weights = wts)
  )

  # The continuous effect label supports only identity, logit, and log. A probit
  # msm fails at entry with a message listing the supported links, rather than
  # reaching the terse "Unsupported outcome link" message downstream.
  expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_link_error",
    regexp = "identity"
  )
})

test_that("ipw() continuous works with a log-link marginal structural model", {
  dat <- sim_continuous()
  # Low baseline risk and a weak dose effect keep the log-binomial fitted
  # probabilities below 1 so the marginal structural model converges.
  withr::local_seed(11)
  dat$ylow <- rbinom(
    nrow(dat),
    1,
    exp(-1.6 + 0.15 * dat$A + 0.1 * dat$x1)
  )
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  msm <- suppressWarnings(glm(
    ylow ~ A,
    data = dat,
    family = binomial(link = "log"),
    weights = wts,
    start = c(-1.6, 0.15),
    control = glm.control(epsilon = 1e-14, maxit = 200)
  ))

  res <- ipw(ps_mod, msm)
  expect_equal(res$estimates$effect, "log(rr)")
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
})

# ---- ps_link on the continuous path -----------------------------------------
#
# ps_link is meaningful only for a binomial glm on the binary path, where it
# overrides the model's link. A continuous propensity model has no such link, so
# both continuous entries, the lm method and the gaussian-family branch of
# ipw.glm, must reject it rather than silently ignore it. The default
# ps_link = NULL is covered by every other test in this file, through both
# entries (see the end-to-end and gaussian-glm routing tests above).

test_that("ipw() rejects ps_link for an lm propensity model", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit"),
    class = "propensity_ipw_link_error",
    regexp = "ps_link"
  )
})

test_that("ipw() rejects ps_link for a gaussian glm propensity model", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, ps_type = "glm")

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit"),
    class = "propensity_ipw_link_error",
    regexp = "ps_link"
  )
})

test_that("the gaussian glm ps_link rejection deprecates nothing", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, ps_type = "glm")

  # `ps_link` is deprecated on the binary path, and the gaussian branch of
  # ipw.glm returns before that warning is reached. A model with no link to
  # override must be rejected outright rather than told to drop the argument
  # and to correct it in the same breath.
  #
  # The signal is forced unconditionally, so the empty record below is evidence
  # that this path raises nothing rather than that something else reached the
  # `ipw(ps_link)` id first.
  deprecations <- character()
  err <- with_always_deprecated(
    withCallingHandlers(
      expect_error(
        ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit"),
        class = "propensity_ipw_link_error"
      ),
      lifecycle_warning_deprecated = function(cnd) {
        deprecations <<- c(deprecations, conditionMessage(cnd))
        rlang::cnd_muffle(cnd)
      }
    )
  )

  expect_length(deprecations, 0)

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "continuous propensity score model", fixed = TRUE)
})

test_that("the continuous ps_link error explains why the argument does not apply", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit")
  )
})

# ---- factor outcome response ------------------------------------------------

test_that("continuous logistic MSM matches a factor outcome response to the numeric fit", {
  dat <- sim_continuous()
  dat$ybf <- factor(ifelse(dat$yb == 1, "yes", "no"), levels = c("no", "yes"))
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)

  num <- glm(
    yb ~ A,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = ctrl
  )
  fac <- suppressWarnings(
    glm(
      ybf ~ A,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = ctrl
    )
  )
  res_num <- ipw(ps_mod, num)$estimates
  res_fac <- ipw(ps_mod, fac)$estimates
  expect_equal(res_fac$estimate, res_num$estimate, tolerance = 1e-8)
  expect_equal(res_fac$std.err, res_num$std.err, tolerance = 1e-8)
})

# ---- ps design reconstruction from .data ------------------------------------

test_that("continuous reconstructs the ps design from .data when the fit frame is gone", {
  # The continuous path routes the ps design through ipw_extract_ps_design:
  # with the fitting data gone, no .data raises the guarded error and .data
  # rebuilds the design. The binary path matches this behavior; see
  # test-ipw-mestimation.R.
  dat <- sim_continuous()
  ps_gone <- lm(A ~ x1 + x2, data = dat, model = FALSE)
  wts <- continuous_weights(fitted(ps_gone), dat$A)
  om <- lm(yc ~ A, data = dat, weights = wts)
  ps_ref <- lm(A ~ x1 + x2, data = dat)
  dat_copy <- dat
  rm(dat)

  expect_error(ipw(ps_gone, om), class = "propensity_ipw_data_error")

  ref <- ipw(ps_ref, om)$estimates
  res <- ipw(ps_gone, om, .data = dat_copy)$estimates
  expect_equal(res, ref, tolerance = 1e-8)
})

test_that("continuous errors informatively when the outcome model frame is gone", {
  # The companion to the binary and categorical pins in test-ipw-mestimation.R.
  # The guard lives in the weight extraction every ipw() method runs first, so
  # the continuous route gets it by sharing that helper rather than by having
  # its own check. This pin exists so a refactor that moves this path off the
  # shared helper is caught rather than silently losing the guided error.
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  outcome_gone <- local({
    d_local <- dat
    w_local <- wts
    m <- lm(yc ~ A, data = d_local, weights = w_local, model = FALSE)
    rm(d_local, w_local)
    m
  })

  # the fixture is genuinely unreconstructable
  expect_error(stats::model.frame(outcome_gone))

  for (supplied in list(NULL, dat)) {
    err <- expect_error(
      ipw(ps_mod, outcome_gone, .data = supplied),
      class = "propensity_ipw_data_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "outcome_mod", fixed = TRUE)
  }
})

test_that("continuous rejects a .data whose covariate types skew the design", {
  # The continuous companion to the binary and categorical width pins. Same
  # rationale as the frame-gone pin above: the check lives in the shared
  # extraction, and this fixes that it covers this route.
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  dat_skew <- dat
  dat_skew$x1 <- cut(dat$x1, breaks = 3, labels = c("lo", "mid", "hi"))
  expect_equal(length(coef(mods$ps_mod)), 3)
  expect_equal(ncol(model.matrix(~ x1 + x2, data = dat_skew)), 4)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat_skew),
    class = "propensity_ipw_data_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, ".data", fixed = TRUE)
})

test_that("a log-transformed continuous outcome response through .data matches the model-frame route", {
  # The continuous companion to the binary pins in test-ipw-se-method.R. Same
  # rationale as the pins above: the response is read out of `.data` by the
  # shared extraction, so a transformation of a column has to be computed rather
  # than looked up under the name of the function wrapping it. The exponentiated
  # outcome keeps `log()` defined on every observation.
  dat <- sim_continuous()
  dat$yp <- exp(dat$yc)
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  outcome_mod <- lm(log(yp) ~ A, data = dat, weights = wts)

  expect_equal(
    ipw(ps_mod, outcome_mod, .data = dat)$estimates,
    ipw(ps_mod, outcome_mod)$estimates,
    tolerance = 1e-10
  )
})

# ---- a wrong-sized .data is a data problem -----------------------------------
#
# As on the binary and categorical paths, a `.data` whose row count disagrees
# with the fitted models is a data problem and reads as one only if it is
# reported as one. Without the guard it would reach the weight-consistency
# preflight, which talks about weights instead. The binary companion is "mestimation rejects a
# .data with too few rows" in test-ipw-mestimation.R.

test_that("ipw() continuous rejects a .data with too few rows", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat[-1, ]),
    class = "propensity_ipw_data_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "one row per observation", fixed = TRUE)
})

# ---- weight-consistency hint wording ----------------------------------------

test_that("the continuous weight-consistency error names no focal level", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)
  spec <- ipw_spec_continuous(mods$ps_mod, mods$outcome_mod)
  doubled <- 2 * as.double(model.weights(model.frame(mods$outcome_mod)))

  err <- expect_error(
    ipw_check_weight_consistency(spec, doubled),
    class = "propensity_ipw_weights_mismatch_error"
  )

  # A continuous exposure has no focal level, so the focal-level hint the binary
  # and categorical paths offer must be absent rather than reworded. The binary
  # wording is pinned by "the weight-consistency error names a focal level as a
  # possible cause" in test-ipw-mestimation.R.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_false(grepl("focal", msg, fixed = TRUE))

  # The spread of the conditional density is a continuous-only cause, and the
  # one users meet on the way from weights built with observation-level
  # standard deviations.
  expect_match(msg, ".sigma", fixed = TRUE)
})

test_that("ipw() refuses continuous weights built with an observation-level `.sigma`", {
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)

  # `ipw()` stacks a single pooled residual variance, so weights spread by the
  # model's leave-one-out standard deviations are not a function of the stacked
  # parameters at any value. The weights record the spread they were built with,
  # so the refusal names it directly rather than reaching the consistency check
  # and reporting a disagreement between two vectors.
  sigma_wts <- wt_ate(
    as.double(fitted(ps_mod)),
    dat$A,
    .sigma = influence(ps_mod)$sigma,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  outcome_mod <- lm(yc ~ A, data = dat, weights = sigma_wts)

  expect_error(
    ipw(ps_mod, outcome_mod),
    class = "propensity_ipw_sigma_error"
  )

  # The same weights without `.sigma` take the pooled default and are accepted,
  # so the refusal is about the spread rather than anything else in the fit.
  pooled_wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  pooled_outcome_mod <- lm(yc ~ A, data = dat, weights = pooled_wts)
  expect_s3_class(ipw(ps_mod, pooled_outcome_mod), "ipw")
})

test_that("the observation-level spread refusal names both spreads it takes", {
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  sigma_wts <- wt_ate(
    as.double(fitted(ps_mod)),
    dat$A,
    .sigma = influence(ps_mod)$sigma,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  outcome_mod <- lm(yc ~ A, data = dat, weights = sigma_wts)

  expect_propensity_error(ipw(ps_mod, outcome_mod))
})

# ---- propensity model class validation --------------------------------------
#
# The continuous path stacks the score of the propensity model the registry
# recognizes, so a class it does not recognize would be solved to the root of an
# equation the supplied model was not fit by. An lm subclass whose coefficients
# are not that root would therefore yield estimates for a propensity model the
# user never fit, and no downstream guard catches it: an lm subclass reaches the
# lm method by inheritance, and the weights built from its fitted values agree
# with the weights recomputed at the seeded init, so the weight-consistency
# preflight passes. Every class the registry was not written for is refused at
# entry, naming the class it was given; which refusal each one gets is pinned
# with the rest of the registry below. MASS::rlm is the one lm subclass with a
# score of its own, so it is stacked rather than refused, on the terms pinned in
# the robust section below; a gaussian mgcv::gam is stacked as the penalized fit
# it is, on the terms pinned in tests/testthat/test-ipw-gam.R.

# A fixture whose robust and least-squares propensity fits genuinely disagree:
# adding a block of large outliers to the exposure pulls the least-squares fit
# away from the robust fit, so the two imply different weights and different
# estimates. The outcome is regenerated from the contaminated exposure, under
# its own seed, so the marginal structural model remains coherent.
sim_continuous_outliers <- function(seed = 2024, n = 800, n_outliers = 40) {
  dat <- sim_continuous(seed = seed, n = n)
  dat$A[seq_len(n_outliers)] <- dat$A[seq_len(n_outliers)] + 12
  withr::local_seed(seed + 1L)
  dat$yc <- 1 + 0.6 * dat$A + 0.5 * dat$x1 - 0.3 * dat$x2 + rnorm(n)
  # yb came from the pre-contamination exposure and no caller reads it.
  dat$yb <- NULL
  dat
}

test_that("a least-squares propensity model still runs on the outlier fixture", {
  dat <- sim_continuous_outliers()

  # The control for the robust section below: on the contaminated fixture the
  # least-squares route is unchanged, so anything the robust route reports
  # differently comes from the fit rather than from the data.
  lm_mod <- lm(A ~ x1 + x2, data = dat)
  lm_wts <- continuous_weights(as.double(fitted(lm_mod)), dat$A)
  lm_msm <- lm(yc ~ A, data = dat, weights = lm_wts)
  expect_s3_class(ipw(lm_mod, lm_msm), "ipw")
})

# ---- arguments that fall into the dots --------------------------------------

test_that("ipw() lm rejects arguments that fall into the dots", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  # `.data` was the third positional argument in earlier releases; it now lands
  # in the dots, where it would be discarded without a signal.
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, dat),
    class = "rlib_error_dots_nonempty"
  )

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, nonsense = 42),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("ipw() lm accepts every argument supplied by name", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  baseline <- ipw(mods$ps_mod, mods$outcome_mod)
  named <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = dat,
    estimand = "ate",
    conf_level = 0.95,
    se_method = "mestimation"
  )

  expect_equal(named$estimates, baseline$estimates)
})

# ---- the propensity model's response has to be a column ---------------------
#
# The continuous path derives the exposure's name by deparsing the propensity
# model's left-hand side, so a left-hand side that is not a bare symbol gives
# back a vector of several names and everything downstream is indexed by it. A
# matrix response and a transformed response both do this, and neither was
# guarded here: without `.data` the call died on an internal length assertion
# about `.exposure_name`, and with `.data` it asked for a column named after the
# transforming function, which cannot exist. The binary path already rejects
# both shapes, and the contract is the same one for a continuous exposure: the
# propensity model's response is the exposure column itself.

continuous_lhs_data <- function() {
  dat <- sim_continuous()
  # a strictly positive copy of the exposure, so log() of it is defined, and a
  # precomputed column holding that transformation
  dat$Apos <- exp(dat$A)
  dat$A1 <- round(pmax(dat$A, 0.1) * 10)
  dat$A2 <- round(pmax(-dat$A, 0.1) * 10)
  dat
}

# Continuous ATE weights and an MSM to carry them, built from whatever exposure
# vector the propensity model models. The models are only ever rejected here, so
# the weights need to be well formed rather than matched to the fit.
continuous_lhs_outcome <- function(dat, ps_mod, exposure) {
  wts <- continuous_weights(
    as.double(fitted(ps_mod))[seq_along(exposure)],
    exposure
  )
  lm(yc ~ A, data = dat, weights = wts)
}

test_that("a matrix-response continuous propensity model is rejected on its shape", {
  dat <- continuous_lhs_data()
  ps_mod <- lm(cbind(A1, A2) ~ x1 + x2, data = dat)
  ps_ok <- lm(A ~ x1 + x2, data = dat)
  outcome_mod <- continuous_lhs_outcome(dat, ps_ok, dat$A)

  # the fixture is the shape it claims to be
  expect_true(is.matrix(stats::model.response(stats::model.frame(ps_mod))))

  err <- expect_error(
    ipw(ps_mod, outcome_mod, .data = dat),
    class = "propensity_ipw_response_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "matrix response", fixed = TRUE)

  # and the spec constructor, which callers can reach directly, rejects it the
  # same way rather than on the internal length assertion it used to
  err_spec <- expect_error(
    ipw_spec_continuous(ps_mod, outcome_mod, .data = dat),
    class = "propensity_ipw_response_error"
  )
  expect_false(grepl(
    ".exposure_name",
    conditionMessage(err_spec),
    fixed = TRUE
  ))
})

test_that("a transformed continuous propensity response is rejected on both routes", {
  dat <- continuous_lhs_data()
  ps_mod <- lm(log(Apos) ~ x1 + x2, data = dat)
  outcome_mod <- continuous_lhs_outcome(dat, ps_mod, log(dat$Apos))

  for (data_arg in list(NULL, dat)) {
    err <- expect_error(
      ipw(ps_mod, outcome_mod, .data = data_arg),
      class = "propensity_ipw_response_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "log(Apos)", fixed = TRUE)
    # the two errors it used to raise, one per route
    expect_false(grepl(".exposure_name", msg, fixed = TRUE))
    expect_false(grepl("\"log\" column", msg, fixed = TRUE))
  }
})

test_that("a transformed continuous propensity response is rejected in the spec", {
  dat <- continuous_lhs_data()
  ps_mod <- lm(log(Apos) ~ x1 + x2, data = dat)
  outcome_mod <- continuous_lhs_outcome(dat, ps_mod, log(dat$Apos))

  expect_error(
    ipw_spec_continuous(ps_mod, outcome_mod, .data = dat),
    class = "propensity_ipw_response_error"
  )
})

test_that("the lm and gaussian glm propensity routes reject a transformed response alike", {
  dat <- continuous_lhs_data()
  ps_lm <- lm(log(Apos) ~ x1 + x2, data = dat)
  ps_glm <- glm(log(Apos) ~ x1 + x2, data = dat, family = gaussian())
  outcome_mod <- continuous_lhs_outcome(dat, ps_lm, log(dat$Apos))

  # The two propensity models share fitted values, so they have to reject the
  # same shape with the same message. The gaussian glm reached the binary
  # guard's cbind wording, which named a matrix this model does not have.
  messages <- vapply(
    list(ps_lm, ps_glm),
    function(ps_mod) {
      err <- expect_error(
        ipw(ps_mod, outcome_mod, .data = dat),
        class = "propensity_ipw_response_error"
      )
      gsub("[[:space:]]+", " ", conditionMessage(err))
    },
    character(1)
  )

  expect_equal(messages[[1]], messages[[2]])
  expect_false(grepl("cbind", messages[[1]], fixed = TRUE))
})

test_that("the transformed continuous propensity response error reads in the user's terms", {
  dat <- continuous_lhs_data()
  ps_mod <- lm(log(Apos) ~ x1 + x2, data = dat)
  outcome_mod <- continuous_lhs_outcome(dat, ps_mod, log(dat$Apos))

  expect_snapshot(
    error = TRUE,
    ipw(ps_mod, outcome_mod, .data = dat)
  )
})

test_that("a plain continuous propensity response is still accepted", {
  dat <- continuous_lhs_data()
  mods <- fit_continuous_models(dat)

  expect_s3_class(ipw(mods$ps_mod, mods$outcome_mod, .data = dat), "ipw")
})

# ---- the contrast column a continuous exposure has no use for ---------------

test_that("ipw() continuous stores no contrast column under either name", {
  dat <- sim_continuous()

  # A continuous exposure has no pair of levels to name, so its rows are keyed
  # by the effect measure alone and the stored table carries no column for a
  # pair, under the canonical name or under the one it replaces.
  for (family in c("gaussian", "binomial")) {
    mods <- fit_continuous_models(dat, outcome_family = family)
    est <- ipw(mods$ps_mod, mods$outcome_mod)$estimates
    expect_false("contrast" %in% names(est))
    expect_false("comparison" %in% names(est))
  }
})

# ---- the model registry -----------------------------------------------------
#
# What `ipw()` can stack for a continuous exposure is the set of propensity
# model classes whose score it can write down: a plain lm, and a gaussian glm
# read through its link. Everything else is either refused by class, because the
# stacked score would describe a model the user never fit, or refused as a
# standard error the M-estimation path cannot produce, which is a different
# refusal with a different remedy.

test_that("ipw() stacks a log-link gaussian propensity model", {
  dat <- sim_continuous_positive()
  mods <- fit_continuous_models(dat, ps_type = "glm_log")

  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_equal(
    res$estimates$estimate,
    unname(coef(mods$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )

  # The propensity block leads theta, so the coefficients the sandwich solves
  # for are the first ones it reports, and they are the coefficients the user
  # fit: the stacked score is the log-link score rather than the least-squares
  # one, which would sit at a different root.
  theta <- coef(res$fit)
  ps_block <- theta[seq_along(coef(mods$ps_mod))]
  expect_equal(
    unname(ps_block),
    unname(coef(mods$ps_mod)),
    tolerance = 1e-6
  )
})

# ---- a robust Huber propensity model ----------------------------------------
#
# `MASS::rlm()` descends its own loss rather than the sum of squares, so its
# coefficients are the root of the Huber score read at the scale the fit settled
# on. The stacked system writes that score with deli, holding the fit's MAD
# scale fixed as a known constant: rlm clips the standardized residual and deli
# clips the raw one, so the two agree when deli's threshold is the psi's own
# constant times that scale. The psi carries its constant in its formals, which
# is where a caller's `k` is written and where the tuning constant of the scale
# estimator, `k2`, is not. The spread the density ratio divides by is the pooled
# residual root mean square, as it is for every other class; a caller who wants
# the robust scale there passes it as `.sigma`, which is the fixed-spread path.
#
# Only the Huber psi is written here, at whatever threshold the fit clipped at.
# A redescending psi and the MM method are the roots of other equations, and an
# unconverged fit is the root of none, so all three are refused rather than
# stacked at a point the user's fit does not occupy.

test_that("ipw() stacks a robust Huber propensity model", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()
  mods <- fit_continuous_models(dat, ps_type = "rlm")

  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_equal(
    res$estimates$estimate,
    unname(coef(mods$outcome_mod)[["A"]]),
    tolerance = 1e-6
  )

  # The propensity block leads theta, so the coefficients the sandwich solves
  # for are the first ones it reports, and they are the robust fit's rather than
  # the least-squares ones the same design would give.
  theta <- coef(res$fit)
  ps_block <- theta[seq_along(coef(mods$ps_mod))]
  expect_equal(unname(ps_block), unname(coef(mods$ps_mod)), tolerance = 1e-6)

  # The contamination is what makes that distinction visible: the two fits
  # disagree in the first decimal place on this fixture, so a stacked system
  # that solved the least-squares score would report a different propensity
  # model and be caught here.
  lm_coefs <- coef(lm(A ~ x1 + x2, data = dat))
  expect_gt(max(abs(unname(coef(mods$ps_mod)) - unname(lm_coefs))), 0.1)
})

test_that("the weights ipw() rebuilds at its seed match a robust fit's", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()
  mods <- fit_continuous_models(dat, ps_type = "rlm")

  # The weights a user would build straight from the fitted model, rather than
  # from its fitted values, are the ones the preflight has to reproduce: the
  # `rlm` method reads the exposure off the model's own response.
  direct <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(mods$ps_mod, exposure_type = "continuous", stabilize = TRUE)
  )
  expect_equal(as.double(direct), as.double(mods$wts), tolerance = 1e-12)

  spec <- ipw_spec_continuous(mods$ps_mod, mods$outcome_mod)
  layout <- ipw_theta_layout(spec)
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(direct),
    tolerance = 1e-12
  )
})

test_that("ipw() stacks a robust fit whose psi threshold was retuned", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # Passing `k` retunes the Huber psi rather than replacing it: `rlm` rewrites
  # the formals of the psi it was given, so the fit is still a Huber fit and is
  # still one the stacked system can write, at the threshold the caller chose.
  # A guard that recognized the psi by identity alone would refuse this fit, and
  # a score written at the fit's `k2` would sit at the root of an equation this
  # fit does not solve.
  ps_mod <- MASS::rlm(A ~ x1 + x2, data = dat, acc = 1e-10, k = 2)
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  res <- ipw(ps_mod, msm)

  expect_equal(
    res$estimates$estimate,
    unname(coef(msm)[["A"]]),
    tolerance = 1e-6
  )

  theta <- coef(res$fit)
  ps_block <- theta[seq_along(coef(ps_mod))]
  expect_equal(unname(ps_block), unname(coef(ps_mod)), tolerance = 1e-6)
})

# ---- the redescending psi functions -----------------------------------------
#
# `MASS::psi.bisquare` and `MASS::psi.hampel` are the other two psi functions
# `rlm` ships, and deli writes both as losses of its own: the bisquare is deli's
# "tukey" and Hampel is deli's "hampel". The translation is the one the Huber
# route already makes, applied per psi: `rlm` clips the residual divided by its
# scale estimate and deli clips the raw residual, so each of the psi's own
# constants is multiplied by the scale the fit settled on. Which constants those
# are differs by psi, and each is read off the formals `rlm` rewrote: `c` for
# the bisquare, and `a`, `b`, and `c` for Hampel.
#
# Both psi functions redescend, so the equation the stacked system writes has
# more than one root and the solve reports whichever one it is seeded at. The
# seed is the fit's own coefficients, so what the sandwich describes is the fit
# the user has, read locally: it is the covariance of that root and carries no
# claim about the others. The tests below pin the seeding, which is what makes
# the local reading the right one.
#
# `method = "MM"` stays refused. Its final step converges the bisquare at a
# fixed scale, so its coefficients do sit at a root of the score this route
# would write; what the route cannot write is how it got there. The
# high-breakdown start decides which root the fit finishes at and supplies the
# scale it clips at, and a sandwich read locally at those coefficients
# describes neither, so it is not the sampling behavior of an MM fit.

test_that("ipw() stacks a bisquare propensity model at its own coefficients", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  ps_mod <- MASS::rlm(
    A ~ x1 + x2,
    data = dat,
    psi = MASS::psi.bisquare,
    acc = 1e-10,
    maxit = 200
  )
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  res <- ipw(ps_mod, msm)

  expect_s3_class(res, "ipw")
  expect_equal(
    res$estimates$estimate,
    unname(coef(msm)[["A"]]),
    tolerance = 1e-6
  )

  # The parity that matters: the solved propensity block is the coefficient
  # vector `MASS::rlm()` reports, so the sandwich describes the fit the user
  # has rather than another root of the same redescending score.
  theta <- coef(res$fit)
  ps_block <- theta[seq_along(coef(ps_mod))]
  expect_equal(unname(ps_block), unname(coef(ps_mod)), tolerance = 1e-6)

  # The contamination is what makes the parity worth asserting: a bisquare fit
  # and a least-squares fit of the same columns disagree here, so a system that
  # stacked the wrong score would report a different propensity model.
  lm_coefs <- coef(lm(A ~ x1 + x2, data = dat))
  expect_gt(max(abs(unname(coef(ps_mod)) - unname(lm_coefs))), 0.1)

  expect_true(all(is.finite(res$estimates$std.err)))
  expect_gt(res$estimates$std.err, 0)
})

test_that("ipw() stacks a Hampel propensity model at its own coefficients", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  ps_mod <- MASS::rlm(
    A ~ x1 + x2,
    data = dat,
    psi = MASS::psi.hampel,
    acc = 1e-10,
    maxit = 200
  )
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  res <- ipw(ps_mod, msm)

  expect_equal(
    res$estimates$estimate,
    unname(coef(msm)[["A"]]),
    tolerance = 1e-6
  )

  theta <- coef(res$fit)
  ps_block <- theta[seq_along(coef(ps_mod))]
  expect_equal(unname(ps_block), unname(coef(ps_mod)), tolerance = 1e-6)

  lm_coefs <- coef(lm(A ~ x1 + x2, data = dat))
  expect_gt(max(abs(unname(coef(ps_mod)) - unname(lm_coefs))), 0.1)

  expect_true(all(is.finite(res$estimates$std.err)))
  expect_gt(res$estimates$std.err, 0)
})

test_that("ipw() reads a redescending psi's constants off the fit", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # Each psi names its constants differently and `rlm` records a caller's
  # choice by rewriting the psi's formals, so a translation that read Huber's
  # `k` for every psi would find nothing here, and one that read the defaults
  # rather than the formals would write the score of a fit nobody made. The
  # Hampel case is the one that pins the arity as well: three constants rather
  # than one.
  bisquare <- MASS::rlm(
    A ~ x1 + x2,
    data = dat,
    psi = MASS::psi.bisquare,
    c = 3,
    acc = 1e-10,
    maxit = 200
  )
  expect_equal(as.numeric(formals(bisquare$psi)$c), 3)

  hampel <- MASS::rlm(
    A ~ x1 + x2,
    data = dat,
    psi = MASS::psi.hampel,
    a = 1.5,
    b = 3,
    c = 6,
    acc = 1e-10,
    maxit = 200
  )
  expect_equal(
    as.numeric(c(
      formals(hampel$psi)$a,
      formals(hampel$psi)$b,
      formals(hampel$psi)$c
    )),
    c(1.5, 3, 6)
  )

  for (ps_mod in list(bisquare, hampel)) {
    wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
    msm <- lm(yc ~ A, data = dat, weights = wts)

    res <- ipw(ps_mod, msm)
    theta <- coef(res$fit)
    expect_equal(
      unname(theta[seq_along(coef(ps_mod))]),
      unname(coef(ps_mod)),
      tolerance = 1e-6
    )
  }
})

test_that("the MM method stays refused now that the bisquare is written", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # `rlm(method = "MM")` finishes on the bisquare, which this route can now
  # write, and it is still refused: the high-breakdown start decides which root
  # the fit finishes at and supplies the scale it clips at, neither of which the
  # stacked system writes, so a sandwich read locally at those coefficients
  # would not describe how the fit behaves across samples.
  ps_mod <- MASS::rlm(A ~ x1 + x2, data = dat, method = "MM")
  expect_equal(ipw_rlm_psi_name(ps_mod), "MASS::psi.bisquare")

  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_robust_psi_error"
  )
  expect_match(conditionMessage(err), "MM", fixed = TRUE)
})

test_that("ipw() refuses a robust fit whose psi it cannot write", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # A psi of the caller's own is the refusal that survives: `rlm` accepts any
  # function here, and one this route has no name for is one whose score it
  # cannot write. The three psi functions MASS ships are all written now, so the
  # refusal is reached through a custom psi rather than through the bisquare.
  own_psi <- function(u, k = 1.5) pmin(1, k / abs(u))
  ps_mod <- suppressWarnings(MASS::rlm(
    A ~ x1 + x2,
    data = dat,
    psi = own_psi,
    acc = 1e-10,
    maxit = 200
  ))
  expect_null(ipw_rlm_psi_name(ps_mod))

  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_robust_psi_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "recognize", fixed = TRUE)
  expect_match(msg, "bootstrap the whole fit yourself", fixed = TRUE)
})

test_that("ipw() refuses a robust bisquare fit that did not converge", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # The convergence guard is not psi-specific, and it matters more for a
  # redescending psi than for Huber: a fit stopped short sits at no root at all,
  # and a solve seeded there would walk to whichever root it reached.
  ps_mod <- suppressWarnings(MASS::rlm(
    A ~ x1 + x2,
    data = dat,
    psi = MASS::psi.bisquare,
    maxit = 1
  ))
  expect_false(ps_mod$converged)

  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_convergence_error"
  )
})

test_that("the registry translates each rlm psi to its deli loss and scale", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # The translation itself, read off the registry entry rather than through a
  # solve: which loss deli writes, and the constants on the raw-residual scale
  # deli clips at, which are the psi's own constants times the scale `rlm`
  # settled on.
  fits <- list(
    huber = MASS::rlm(A ~ x1 + x2, data = dat, acc = 1e-10, maxit = 200),
    tukey = MASS::rlm(
      A ~ x1 + x2,
      data = dat,
      psi = MASS::psi.bisquare,
      acc = 1e-10,
      maxit = 200
    ),
    hampel = MASS::rlm(
      A ~ x1 + x2,
      data = dat,
      psi = MASS::psi.hampel,
      acc = 1e-10,
      maxit = 200
    )
  )

  expected_k <- list(
    huber = as.numeric(formals(fits$huber$psi)$k) * fits$huber$s,
    tukey = as.numeric(formals(fits$tukey$psi)$c) * fits$tukey$s,
    hampel = as.numeric(c(
      formals(fits$hampel$psi)$a,
      formals(fits$hampel$psi)$b,
      formals(fits$hampel$psi)$c
    )) *
      fits$hampel$s
  )

  for (loss in names(fits)) {
    entry <- ipw_continuous_model(fits[[loss]])
    expect_true(entry$stackable, label = loss)
    expect_identical(entry$psi_loss, loss, label = loss)
    expect_equal(entry$psi_k, expected_k[[loss]], label = loss)
  }
})

test_that("the robust psi refusal names the psi it found and the remedy", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  own_psi <- function(u, k = 1.5) pmin(1, k / abs(u))
  ps_mod <- suppressWarnings(MASS::rlm(
    A ~ x1 + x2,
    data = dat,
    psi = own_psi,
    acc = 1e-10,
    maxit = 200
  ))
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  expect_propensity_error(ipw(ps_mod, msm))
})

test_that("ipw() refuses the MM method as an equation it cannot write", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # MM reaches its coefficients through a high-breakdown start whose root
  # selection and scale this system does not write, so the refusal is the same
  # one a psi the system cannot write gets.
  ps_mod <- MASS::rlm(A ~ x1 + x2, data = dat, method = "MM")
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_robust_psi_error"
  )
})

test_that("ipw() refuses the MM method named through a variable", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # `rlm` records the method as it was written rather than as it was resolved,
  # so a caller who named it through a variable leaves a symbol in the call. The
  # refusal evaluates it rather than reading the symbol as a method name.
  m <- "MM"
  ps_mod <- MASS::rlm(A ~ x1 + x2, data = dat, method = m)
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_robust_psi_error"
  )
  expect_match(conditionMessage(err), "MM", fixed = TRUE)
})

test_that("the MM refusal names the method and the psi it finishes on", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  ps_mod <- MASS::rlm(A ~ x1 + x2, data = dat, method = "MM")
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  expect_propensity_error(ipw(ps_mod, msm))
})

test_that("ipw() refuses the MM method whose psi it cannot name", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # The refusal reads the method off the call and the psi off the fitted object,
  # and it has a third thing to say when the first is MM and the second is a psi
  # it has no name for: the method alone is enough to refuse on, so the psi it
  # cannot name goes unmentioned. `rlm()` writes no such fit itself, because
  # `method = "MM"` replaces whatever psi it was handed with the bisquare, so
  # the pair is assembled here rather than fit.
  own_psi <- function(u, k = 1.5) pmin(1, k / abs(u))
  ps_mod <- suppressWarnings(MASS::rlm(
    A ~ x1 + x2,
    data = dat,
    psi = own_psi,
    acc = 1e-10,
    maxit = 200
  ))
  ps_mod$call$method <- "MM"
  expect_identical(ipw_rlm_method(ps_mod), "MM")
  expect_null(ipw_rlm_psi_name(ps_mod))

  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_robust_psi_error"
  )
  expect_propensity_error(ipw(ps_mod, msm))
})

test_that("ipw() refuses a robust fit that did not converge", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  # A fit stopped after one step is not at the root of the Huber score, so the
  # stacked system seeded from its coefficients would move away from them and
  # report a fit the user never saw. rlm says so itself, and the refusal reads
  # that flag rather than rediscovering it.
  ps_mod <- suppressWarnings(MASS::rlm(A ~ x1 + x2, data = dat, maxit = 1))
  expect_false(ps_mod$converged)

  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_convergence_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "maxit", fixed = TRUE)
})

test_that("the robust convergence refusal names the arguments that fix it", {
  skip_if_not_installed("MASS")
  dat <- sim_continuous_outliers()

  ps_mod <- suppressWarnings(MASS::rlm(A ~ x1 + x2, data = dat, maxit = 1))
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  expect_propensity_error(ipw(ps_mod, msm))
})

test_that("ipw() takes a robust fit's own scale as a fixed spread", {
  skip_if_not_installed("MASS")

  # The uncontaminated fixture, deliberately: `rlm` reports a scale that resists
  # outliers, so on the contaminated data the conditional density read at that
  # scale is far narrower than the residuals are and the ratio degenerates onto
  # a handful of units. The combination being pinned here is the fixed spread,
  # not the contamination.
  dat <- sim_continuous()
  ps_mod <- MASS::rlm(A ~ x1 + x2, data = dat, acc = 1e-10)
  wts <- continuous_weights(
    as.double(fitted(ps_mod)),
    dat$A,
    .sigma = ps_mod$s
  )
  msm <- lm(yc ~ A, data = dat, weights = wts)

  spec <- ipw_spec_continuous(ps_mod, msm)
  layout <- ipw_theta_layout(spec)
  expect_length(layout$idx$ps, ncol(spec$ps$X))
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(wts),
    tolerance = 1e-12
  )

  res <- ipw(ps_mod, msm)
  expect_false("sigma2_d" %in% names(coef(res$fit)))
  expect_equal(
    res$estimates$estimate,
    unname(coef(msm)[["A"]]),
    tolerance = 1e-6
  )
})

test_that("ipw() continuous refuses an lm subclass it does not know", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  # A subclass reaches the lm method by inheritance, and nothing downstream can
  # tell that its coefficients are not the least-squares root, so the registry
  # refuses any class it was not written for and says which ones it was.
  unknown <- structure(mods$ps_mod, class = c("mymodel", "lm"))

  err <- expect_error(
    ipw(unknown, mods$outcome_mod),
    class = "propensity_class_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "mymodel", fixed = TRUE)
  expect_match(msg, "lm", fixed = TRUE)
  expect_match(msg, "glm", fixed = TRUE)
})

test_that("the unknown-subclass error names the classes that are supported", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)
  unknown <- structure(mods$ps_mod, class = c("mymodel", "lm"))

  expect_propensity_error(ipw(unknown, mods$outcome_mod))
})

test_that("ipw() refuses a propensity model fit without a formula", {
  dat <- sim_continuous()

  # A robust fit through the matrix interface records no formula and no terms.
  # The weights read nothing but its fitted values, so it builds them, and then
  # every `ipw()` route deparses the left-hand side of the formula to name the
  # exposure and finds none there.
  fit <- MASS::rlm(cbind(1, dat$x1, dat$x2), dat$A)
  wts <- continuous_weights(as.numeric(fitted(fit)), dat$A)
  outcome_mod <- lm(yc ~ A, data = dat, weights = wts)

  err <- expect_error(
    ipw(fit, outcome_mod),
    class = "propensity_ipw_response_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))

  # The refusal names the interface the model has to be refit through, rather
  # than reporting a subscript of a formula that was never written.
  expect_match(msg, "formula", fixed = TRUE)
  expect_match(msg, "wt_mod", fixed = TRUE)

  expect_propensity_error(ipw(fit, outcome_mod))
})

test_that("ipw() refuses a propensity link it cannot write the score for", {
  dat <- sim_continuous_positive()

  # An identity link is least squares and a log link is the score deli writes
  # for a normal model of exp(X alpha). The remaining gaussian links are refused
  # by name: the coefficients an IRLS iteration stops at under one of them are
  # not a tight enough root to seed the solve from, so a stacked system built on
  # them would not sit where the user's fit does.
  for (link in c("inverse", "sqrt")) {
    ps_mod <- glm(A ~ x1 + x2, data = dat, family = gaussian(link = link))
    wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
    msm <- lm(yc ~ A, data = dat, weights = wts)

    err <- expect_error(ipw(ps_mod, msm), class = "propensity_ipw_link_error")
    expect_match(conditionMessage(err), link, fixed = TRUE)
  }
})

test_that("the continuous propensity-link error names the link and the remedy", {
  dat <- sim_continuous_positive()
  ps_mod <- glm(A ~ x1 + x2, data = dat, family = gaussian(link = "inverse"))
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  expect_propensity_error(ipw(ps_mod, msm))
})

# ---- a dose model of a family that is not gaussian --------------------------
#
# The exposure type a fitted propensity score model implies is read off its
# family: a gaussian one is a dose, and every other family is a binary
# propensity score. A dose model of a third family is neither, and reading it as
# the second of the two describes it wrongly: the refusal a reader gets names
# the links a binary propensity score is written for, which is a set the fit
# they have could not be refit into and a remedy for a model they did not fit.
#
# What such a fit needs said is the family it has and the family a continuous
# exposure needs, so the refusal is written on the family rather than on the
# link a family that reached the wrong route happens to carry. The links a
# binomial fit is refused for are unchanged: a binomial model is a binary
# propensity score whatever its link, and it keeps the refusal that names them.

test_that("ipw() refuses a non-gaussian glm dose model by its family", {
  dat <- sim_continuous_positive()

  ps_mod <- glm(A ~ x1 + x2, data = dat, family = Gamma())
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, msm),
    class = "propensity_model_family_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))

  # The family the fit has, and the family each exposure type needs. Naming the
  # binary links here is what the old refusal did and is what misdescribed the
  # fit, so the message must not fall back on them.
  expect_match(msg, "Gamma", fixed = TRUE)
  expect_match(msg, "gaussian", fixed = TRUE)
  expect_match(msg, "continuous", fixed = TRUE)
  expect_no_match(msg, "cloglog", fixed = TRUE)
})

test_that("ipw() refuses a non-gaussian gam dose model by its family", {
  skip_if_not_installed("mgcv")
  dat <- sim_continuous_positive()

  # A gam reaches the same routing through the glm method it inherits from, so
  # the fault and the refusal are shared.
  ps_mod <- mgcv::gam(A ~ s(x1) + x2, data = dat, family = Gamma())
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, msm),
    class = "propensity_model_family_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "Gamma", fixed = TRUE)
  expect_match(msg, "gaussian", fixed = TRUE)
  expect_no_match(msg, "cloglog", fixed = TRUE)
})

test_that("ipw() refuses a poisson dose model by its family", {
  dat <- sim_continuous_positive()

  # A count dose is the case that reaches this most naturally, and the family
  # named in the message is the one the fit carries rather than a family the
  # refusal assumed.
  withr::local_seed(31)
  dat$count <- rpois(nrow(dat), exp(0.2 * dat$x1 + 0.4))
  ps_mod <- glm(count ~ x1 + x2, data = dat, family = poisson())
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$count)
  msm <- lm(yc ~ count, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, msm),
    class = "propensity_model_family_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "poisson", fixed = TRUE)
  expect_match(msg, "gaussian", fixed = TRUE)
})

test_that("a binomial propensity model keeps its link refusal", {
  dat <- sim_continuous()

  # The family check must not swallow the refusal a binomial fit already gets.
  # A binomial model is a binary propensity score whatever its link, so an
  # unsupported link is still reported as a link.
  dat$z <- as.numeric(dat$A > stats::median(dat$A))
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial(link = "cauchit"))
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(as.double(fitted(ps_mod)), dat$z, exposure_type = "binary")
  )
  msm <- lm(yc ~ z, data = dat, weights = wts)

  expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_link_error"
  )
})

test_that("ipw() refuses kernel-density weights as a method it cannot produce", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, .density = "kernel")

  # A kernel bandwidth is chosen from the residuals rather than written as a
  # function of the parameters, so the weights are not differentiable in theta
  # and the sandwich has nothing to differentiate. What is unavailable is the
  # method, not the model: the weights are exactly what they claim to be, and
  # the refusal points at a bootstrap of the whole pipeline written by hand.
  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_se_method_unavailable_error"
  )
  expect_s3_class(err, "propensity_method_error")
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(err)),
    "bootstrap the whole fit yourself",
    fixed = TRUE
  )
})

test_that("the kernel-density refusal names the bandwidth as the reason", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, .density = "kernel")

  expect_propensity_error(ipw(mods$ps_mod, mods$outcome_mod))
})

test_that("ipw() stacks a density the user wrote", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(
    dat,
    .density = function(z) stats::dt(z, df = 5)
  )

  # A density the user supplies is a closure the weights carry, and the sandwich
  # calls it as given: nothing about it depends on the parameters except the
  # standardized residual it is read at.
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  expect_equal(
    res$estimates$estimate,
    unname(coef(mods$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )
})

# ---- densities and numerators end to end ------------------------------------
#
# The combinations the sandwich has to rebuild weights for. Each is a family and
# a numerator; the spread is pooled throughout, which is the case a fixed
# `.sigma` and an observation-level one are held against below.
continuous_density_cases <- list(
  list(.density = dens_t(4), numerator = "marginal"),
  list(.density = "laplace", numerator = "marginal"),
  list(.density = "normal", numerator = "integrated"),
  list(.density = dens_t(4), numerator = "integrated")
)

test_that("ipw() continuous point estimates match the MSM coefficient for every density", {
  dat <- sim_continuous()

  for (case in continuous_density_cases) {
    mods <- fit_continuous_models(
      dat,
      .density = case$.density,
      numerator = case$numerator
    )
    res <- ipw(mods$ps_mod, mods$outcome_mod)
    expect_equal(
      res$estimates$estimate,
      unname(coef(mods$outcome_mod)[["A"]]),
      tolerance = 1e-8
    )
    expect_true(is.finite(res$estimates$std.err) && res$estimates$std.err > 0)
  }
})

test_that("the weights ipw() rebuilds at its seed are the weights it was given", {
  dat <- sim_continuous()

  # The preflight compares the two at a relative tolerance of 1e-6 before
  # solving anything. They agree far more closely than that when the sandwich
  # reads the same density, numerator, and spread the weights recorded, so this
  # holds them to the arithmetic rather than to the guard.
  cases <- c(
    continuous_density_cases,
    list(list(.density = "normal", numerator = "marginal"))
  )
  for (case in cases) {
    mods <- fit_continuous_models(
      dat,
      .density = case$.density,
      numerator = case$numerator
    )
    spec <- ipw_spec_continuous(mods$ps_mod, mods$outcome_mod)
    layout <- ipw_theta_layout(spec)
    expect_equal(
      as.double(ipw_weights_at_init(spec, layout)),
      as.double(mods$wts),
      tolerance = 1e-12
    )
    expect_s3_class(ipw(mods$ps_mod, mods$outcome_mod), "ipw")
  }
})

test_that("a log-link propensity model rebuilds its weights at the seed too", {
  dat <- sim_continuous_positive()
  mods <- fit_continuous_models(dat, ps_type = "glm_log")

  spec <- ipw_spec_continuous(mods$ps_mod, mods$outcome_mod)
  layout <- ipw_theta_layout(spec)
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(mods$wts),
    tolerance = 1e-12
  )
})

test_that("an integrated numerator stacks no marginal moments", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, numerator = "integrated")

  # The numerator is the conditional density averaged over the units, which is
  # built from the propensity block and the data alone, so there is nothing left
  # for a stabilization block to estimate.
  spec <- ipw_spec_continuous(mods$ps_mod, mods$outcome_mod)
  layout <- ipw_theta_layout(spec)
  expect_length(layout$idx$stab, 0L)

  res <- ipw(mods$ps_mod, mods$outcome_mod)
  theta_names <- names(coef(res$fit))
  expect_false("mu_a" %in% theta_names)
  expect_false("sigma2_a" %in% theta_names)
})

# ---- the spread the weights were built with ---------------------------------

test_that("ipw() accepts a fixed scalar `.sigma` as a constant", {
  dat <- sim_continuous()

  # A spread the user fixed is a number rather than a parameter, so the ps block
  # is the coefficients alone and the conditional variance is the square of what
  # was supplied. The value is deliberately away from the pooled residual
  # standard deviation, so weights rebuilt at the pooled spread would not match.
  mods <- fit_continuous_models(dat, .sigma = 1.25)

  spec <- ipw_spec_continuous(mods$ps_mod, mods$outcome_mod)
  layout <- ipw_theta_layout(spec)
  expect_length(layout$idx$ps, ncol(spec$ps$X))
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(mods$wts),
    tolerance = 1e-12
  )

  res <- ipw(mods$ps_mod, mods$outcome_mod)
  expect_false("sigma2_d" %in% names(coef(res$fit)))
  expect_equal(
    res$estimates$estimate,
    unname(coef(mods$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )
})

# ---- what the weights record and what the spec reads ------------------------

test_that("ipw_spec_continuous reads the density, numerator, and spread off the weights", {
  dat <- sim_continuous()

  mods <- fit_continuous_models(
    dat,
    .density = dens_t(3),
    numerator = "integrated"
  )
  spec <- ipw_spec_continuous(mods$ps_mod, mods$outcome_mod)
  expect_true(density_specs_agree(spec$density, dens_t(3)))
  expect_equal(spec$numerator, "integrated")
  expect_equal(spec$sigma$kind, "pooled")
  expect_null(spec$sigma$value)

  fixed <- fit_continuous_models(dat, .sigma = 1.25)
  spec_fixed <- ipw_spec_continuous(fixed$ps_mod, fixed$outcome_mod)
  expect_equal(spec_fixed$sigma$kind, "fixed")
  expect_equal(spec_fixed$sigma$value, 1.25)
})

test_that("ipw() reads weights carrying no density record as a normal ratio", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  # Weights written by hand, and weights built before the record existed, say
  # only that they were stabilized. That is the normal marginal ratio, which is
  # what the weights were, so the fit is the fit the recorded weights give.
  plain <- psw(as.double(mods$wts), estimand = "ate", stabilized = TRUE)
  plain_mod <- lm(yc ~ A, data = dat, weights = plain)

  res_plain <- ipw(mods$ps_mod, plain_mod)
  res_recorded <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_equal(
    res_plain$estimates$estimate,
    res_recorded$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    res_plain$estimates$std.err,
    res_recorded$estimates$std.err,
    tolerance = 1e-10
  )
})

# ---- standard-error oracles for the new routes ------------------------------

test_that("the sandwich SE for a heavy-tailed integrated ratio tracks a bootstrap", {
  skip_on_cran()

  dat <- sim_continuous(seed = 2024, n = 600)
  mods <- fit_continuous_models(
    dat,
    .density = dens_t(4),
    numerator = "integrated"
  )
  mest_se <- ipw(mods$ps_mod, mods$outcome_mod)$estimates$std.err

  boot_slope <- function(d) {
    ps <- lm(A ~ x1 + x2, data = d)
    w <- continuous_weights(
      as.double(fitted(ps)),
      d$A,
      .density = dens_t(4),
      numerator = "integrated"
    )
    msm <- lm(yc ~ A, data = d, weights = as.double(w))
    unname(coef(msm)[["A"]])
  }

  withr::local_seed(918)
  reps <- 400L
  n <- nrow(dat)
  boot <- vapply(
    seq_len(reps),
    function(i) {
      boot_slope(dat[sample.int(n, n, replace = TRUE), , drop = FALSE])
    },
    numeric(1)
  )

  boot_se <- stats::sd(boot)
  expect_lt(abs(mest_se - boot_se) / boot_se, 0.15)
})

test_that("the sandwich SE for a log-link propensity model tracks a bootstrap", {
  skip_on_cran()

  dat <- sim_continuous_positive(seed = 2024, n = 600)
  mods <- fit_continuous_models(dat, ps_type = "glm_log")
  mest_se <- ipw(mods$ps_mod, mods$outcome_mod)$estimates$std.err

  boot_slope <- function(d) {
    ps <- glm(
      A ~ x1 + x2,
      data = d,
      family = gaussian(link = "log"),
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
    w <- continuous_weights(as.double(fitted(ps)), d$A)
    msm <- lm(yc ~ A, data = d, weights = as.double(w))
    unname(coef(msm)[["A"]])
  }

  withr::local_seed(4102)
  reps <- 400L
  n <- nrow(dat)
  boot <- vapply(
    seq_len(reps),
    function(i) {
      boot_slope(dat[sample.int(n, n, replace = TRUE), , drop = FALSE])
    },
    numeric(1)
  )

  boot_se <- stats::sd(boot)
  expect_lt(abs(mest_se - boot_se) / boot_se, 0.15)
})

test_that("the sandwich SE for a robust Huber propensity model tracks a bootstrap", {
  skip_if_not_installed("MASS")
  skip_on_cran()

  # The sandwich holds the fit's MAD scale fixed while the bootstrap re-estimates
  # it on every resample, so the two are not the same calculation. On the
  # contaminated fixture they still agree well inside the 15 percent band the
  # other continuous oracles use, which is what says the fixed scale costs the
  # standard error little: a prototype of this stack put the gap at 7 percent
  # here, and at 4 percent on the uncontaminated simulation.
  dat <- sim_continuous_outliers(seed = 2024, n = 600)
  mods <- fit_continuous_models(dat, ps_type = "rlm")
  mest_se <- ipw(mods$ps_mod, mods$outcome_mod)$estimates$std.err

  boot_slope <- function(d) {
    # A resample of a contaminated sample can need more than the default twenty
    # steps to reach `acc`, and a replicate that stopped short would be fit at a
    # different tolerance than the sample it is a resample of.
    ps <- MASS::rlm(A ~ x1 + x2, data = d, acc = 1e-10, maxit = 100)
    w <- continuous_weights(as.double(fitted(ps)), d$A)
    msm <- lm(yc ~ A, data = d, weights = as.double(w))
    unname(coef(msm)[["A"]])
  }

  withr::local_seed(918)
  reps <- 400L
  n <- nrow(dat)
  boot <- vapply(
    seq_len(reps),
    function(i) {
      boot_slope(dat[sample.int(n, n, replace = TRUE), , drop = FALSE])
    },
    numeric(1)
  )

  boot_se <- stats::sd(boot)
  expect_lt(abs(mest_se - boot_se) / boot_se, 0.15)
})

# ---- against WeightIt -------------------------------------------------------

# WeightIt computes the same two-step sandwich for a weighted outcome model, so
# its `glm_weightit(vcov = "asympt")` standard error is an independent
# implementation of what `ipw()` reports here. The comparison is only
# meaningful because the weights themselves are the same number: WeightIt's
# default numerator for a continuous exposure is the integrated one and its
# default density is normal, and both packages build the numerator from the
# same 50-point grid interpolated with `stats::spline(method = "fmm")`, so the
# weights agree to floating-point noise (pinned in detail in
# test-weights-continuous.R).
#
# The two sandwiches are assembled differently and still agree to about seven
# significant digits. WeightIt parameterizes the conditional variance as
# `log(s2)`, which is a reparameterization the sandwich for the outcome
# coefficients is invariant to; it folds the treatment score in through a
# correction term rather than stacking one joint set of estimating equations,
# which is the same sandwich for the outcome block algebraically; and it
# differentiates its psi with its own numerical routine while `ipw()` uses
# deli's central differences. What is left is the difference between two
# numerical derivatives of the same function at the same root. The observed
# relative gap at this fixture is under 1e-7 for every density, so the band
# below is set far tighter than a percentage-level check would be while still
# leaving three orders of magnitude of headroom for a different BLAS.
#
# The kernel density is left out: it makes the weights a non-smooth function of
# the parameters, so neither package offers M-estimation for it. The marginal
# numerator is left out too, because WeightIt has no equivalent of it and there
# is nothing to compare against; its standard errors are checked against the
# nonparametric bootstrap oracle earlier in this file.

test_that("ipw() continuous m-estimation standard errors match WeightIt", {
  skip_on_cran()
  skip_if_not_installed("WeightIt", "2.0.0")

  set.seed(123)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  A <- 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  yc <- 1 + 0.5 * A + x1 - 0.7 * x2 + rnorm(n)
  dat <- data.frame(x1 = x1, x2 = x2, A = A, yc = yc)

  densities <- list(
    list(ours = "normal", theirs = NULL),
    list(ours = dens_t(df = 4), theirs = "dt_4"),
    list(ours = "laplace", theirs = "dlaplace")
  )

  for (dens in densities) {
    mods <- fit_continuous_models(
      dat,
      ps_type = "glm",
      .density = dens$ours,
      numerator = "integrated"
    )
    ours <- ipw(mods$ps_mod, mods$outcome_mod)

    theirs_wt <- if (is.null(dens$theirs)) {
      WeightIt::weightit(A ~ x1 + x2, data = dat, method = "glm")
    } else {
      WeightIt::weightit(
        A ~ x1 + x2,
        data = dat,
        method = "glm",
        density = dens$theirs
      )
    }
    theirs <- WeightIt::glm_weightit(
      yc ~ A,
      data = dat,
      weightit = theirs_wt,
      vcov = "asympt"
    )

    # The premise of the comparison: both packages weight by the same numbers.
    expect_equal(
      as.numeric(mods$wts),
      as.numeric(theirs_wt$weights),
      tolerance = 1e-10,
      ignore_attr = TRUE
    )

    expect_equal(nrow(ours$estimates), 1L)
    expect_equal(ours$estimates$effect, "slope")
    expect_equal(
      ours$estimates$estimate,
      unname(coef(theirs)[["A"]]),
      tolerance = 1e-8
    )
    expect_equal(
      ours$estimates$std.err,
      sqrt(vcov(theirs)["A", "A"]),
      tolerance = 1e-4,
      ignore_attr = TRUE
    )
  }
})

# ---- a spread estimated by maximum likelihood -------------------------------

test_that("ipw() solves the scale equation of a t spread by maximum likelihood", {
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  scale <- t_scale_mle(dat$A - as.double(fitted(ps_mod)), df = 6)

  # The same weights twice: once with the scale the package estimates, once with
  # the number it estimates supplied as a constant. The weights are equal, so
  # the point estimates are, and what differs between the two fits is what the
  # sandwich accounts for.
  mle <- fit_continuous_models(dat, .density = dens_t(6, sigma_method = "mle"))
  fixed <- fit_continuous_models(dat, .density = dens_t(6), .sigma = scale)
  expect_equal(as.double(mle$wts), as.double(fixed$wts), tolerance = 1e-8)

  res_mle <- ipw(mle$ps_mod, mle$outcome_mod)
  res_fixed <- ipw(fixed$ps_mod, fixed$outcome_mod)

  # The scale is a parameter of the stacked system under the name the moment
  # equation gives it, and it is solved to the scale the weights were built at.
  # A spread held fixed is no parameter at all.
  expect_true("sigma2_d" %in% names(coef(res_mle$fit)))
  expect_false("sigma2_d" %in% names(coef(res_fixed$fit)))
  expect_equal(
    unname(coef(res_mle$fit)[["sigma2_d"]]),
    scale^2,
    tolerance = 1e-6
  )

  expect_equal(
    res_mle$estimates$estimate,
    res_fixed$estimates$estimate,
    tolerance = 1e-8
  )
  expect_gt(
    abs(res_mle$estimates$std.err / res_fixed$estimates$std.err - 1),
    1e-4
  )
})

test_that("ipw_spec_continuous reads a maximum likelihood spread off the weights", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(
    dat,
    .density = dens_t(6, sigma_method = "mle")
  )

  spec <- ipw_spec_continuous(mods$ps_mod, mods$outcome_mod)
  layout <- ipw_theta_layout(spec)

  # A spread estimated by maximum likelihood is estimated alongside the
  # coefficients, as the pooled one is, so the propensity block carries a row
  # for it beyond the columns of the design.
  expect_length(layout$idx$ps, ncol(spec$ps$X) + 1L)
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(mods$wts),
    tolerance = 1e-12
  )
})

test_that("ipw() refuses a numerator model whose score it cannot write", {
  skip_if_not_installed("MASS")

  # The numerator model joins the stack the way the propensity score model
  # does, so it is read through the same registry and refused on the same
  # terms, in the argument it arrived in. A robust fit descending a psi of the
  # caller's own is the case the registry refuses on both sides: its
  # coefficients are the root of a loss no equation stacked here writes.
  dat <- sim_continuous()
  own_psi <- function(u, k = 1.5) pmin(1, k / abs(u))
  num_mod <- suppressWarnings(MASS::rlm(
    A ~ x1,
    data = dat,
    psi = own_psi,
    acc = 1e-10,
    maxit = 200
  ))
  fits <- fit_continuous_models(dat, stabilize = num_mod)

  expect_propensity_error(ipw(fits$ps_mod, fits$outcome_mod))
})

test_that("ipw() refuses a numerator model fit with case weights", {
  # The stabilization block writes the numerator model's unweighted score, the
  # same score every other model in the stack is written at, so a fit made under
  # case weights is not at the root of the row seeded for it and the solve would
  # walk it off the fit the user has. The refusal names the argument the model
  # arrived in, since a reader told to refit the propensity score model would be
  # told to refit the wrong thing.
  dat <- sim_continuous()
  dat$case_wt <- rep(c(1, 2), length.out = nrow(dat))
  num_mod <- lm(A ~ x1, data = dat, weights = case_wt)
  fits <- fit_continuous_models(dat, stabilize = num_mod)

  expect_error(
    ipw(fits$ps_mod, fits$outcome_mod),
    class = "propensity_ipw_ps_weights_error"
  )
  expect_propensity_error(ipw(fits$ps_mod, fits$outcome_mod))
})

test_that("ipw() refuses a numerator model of another response", {
  # The block the numerator model contributes reads the exposure the propensity
  # score model reads, so a model of anything else would sit away from the root
  # of the row seeded for it and the solve would move it, reporting a numerator
  # nobody fit.
  dat <- sim_continuous()
  num_mod <- lm(yc ~ x1, data = dat)
  fits <- fit_continuous_models(dat, stabilize = num_mod)

  expect_propensity_error(ipw(fits$ps_mod, fits$outcome_mod))
})

test_that("ipw() refuses a numerator model fit to other observations", {
  # The models here were fit under `na.exclude` to a frame with a hole in `x1`,
  # so they analyzed the rows complete over it while the numerator model, which
  # reads a complete column, analyzed every row. The weights are built at the
  # padded length the fits predict at and are the weights the caller meant; what
  # the stacked system cannot do is multiply a design of one length by a block
  # of another.
  dat <- sim_continuous(n = 200)
  dat$x1[c(3, 17, 42)] <- NA

  ps_mod <- lm(A ~ x1 + x2, data = dat, na.action = stats::na.exclude)
  num_mod <- lm(A ~ x2, data = dat)
  wts <- continuous_weights(
    as.double(fitted(ps_mod)),
    dat$A,
    stabilize = num_mod
  )
  outcome_mod <- lm(
    yc ~ A,
    data = dat,
    weights = wts,
    na.action = stats::na.exclude
  )

  expect_propensity_error(ipw(ps_mod, outcome_mod))
})

test_that("ipw() refuses a numerator model with a dropped coefficient", {
  # The stabilization block multiplies the numerator model's coefficients
  # against its design, exactly as the propensity score block multiplies its
  # own, so a column the fit could not separate from the others leaves every
  # numerator undefined. Without the guard the failure arrives much later, as
  # weights that no longer match the ones the caller built, and the report
  # blames the record rather than the model.
  dat <- sim_continuous()
  dat$x1_again <- dat$x1
  num_mod <- lm(A ~ x1 + x1_again, data = dat)
  fits <- fit_continuous_models(dat, stabilize = num_mod)

  expect_error(
    ipw(fits$ps_mod, fits$outcome_mod),
    class = "propensity_ipw_rank_error"
  )
  expect_propensity_error(ipw(fits$ps_mod, fits$outcome_mod))
})
