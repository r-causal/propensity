# Tests for the user-facing se_method switch on ipw.glm and the modified-weight
# guards. The linearization reference estimates below were computed from the
# seeded fixture while ipw() was linearization-only and pin that path exactly.

# Seeded fixture shared across the file. The reference estimates hardcoded below
# were computed from this data while ipw() was linearization-only.
se_method_data <- function() {
  set.seed(2024)
  n <- 400
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.5 * x2))
  y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.4 * x1 - 0.3 * x2))
  data.frame(x1, x2, z, y)
}

se_method_ps_mod <- function(dat) {
  glm(z ~ x1 + x2, data = dat, family = binomial())
}

# Binary outcome model weighted with plain ATE weights.
se_method_outcome_ate <- function(dat, ps_mod) {
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)
  glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
}

# Outcome model weighted with ps_trim-derived weights (refit to avoid the
# no-refit warning).
se_method_outcome_trimmed <- function(dat, ps_mod) {
  withr::local_options(propensity.quiet = TRUE)
  ps <- predict(ps_mod, type = "response")
  trimmed <- ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7)
  refit <- ps_refit(trimmed, model = ps_mod)
  wts <- wt_ate(
    refit,
    .exposure = dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
}

# Outcome model weighted with ps_trunc-derived weights.
se_method_outcome_truncated <- function(dat, ps_mod) {
  withr::local_options(propensity.quiet = TRUE)
  ps <- predict(ps_mod, type = "response")
  truncated <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)
  wts <- wt_ate(
    truncated,
    .exposure = dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
}

# Outcome model weighted with ps_calibrate-derived weights.
se_method_outcome_calibrated <- function(dat, ps_mod) {
  withr::local_options(propensity.quiet = TRUE)
  ps <- predict(ps_mod, type = "response")
  calibrated <- ps_calibrate(ps, .exposure = dat$z, .focal_level = 1)
  wts <- wt_ate(
    calibrated,
    .exposure = dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
}

# Outcome model with an offset term supplied through the formula. The offset
# genuinely shifts the fit (off = 0.3 * x1), so dropping the guard would let
# ipw() converge to a root that disagrees with the fitted model.
se_method_outcome_offset <- function(dat, ps_mod) {
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)
  glm(
    y ~ z + offset(off),
    data = dat,
    family = quasibinomial(),
    weights = wts
  )
}

# The same offset supplied through the `offset` argument rather than the
# formula. terms() stores it separately, so its term labels are the exposure
# alone; only the model-frame offset column reveals it.
se_method_outcome_offset_arg <- function(dat, ps_mod) {
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)
  glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    offset = off
  )
}

# Linearization estimates for the seeded fixture, hardcoded so changes to the
# linearization path are detected.
se_method_lin_reference <- function() {
  data.frame(
    effect = c("rd", "log(rr)", "log(or)"),
    estimate = c(
      0.24308584586316151,
      0.59732185440835850,
      1.02197399113135590
    ),
    std.err = c(
      0.048386022223687541,
      0.126936312969173404,
      0.211670119599652012
    )
  )
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

# Seeded fixture with an added continuous outcome, so the linear-model cases
# mirror the covariate-adjusted repro. The continuous outcome is generated
# deterministically from the shared fixture columns.
se_method_data_cont <- function() {
  dat <- se_method_data()
  set.seed(2025)
  dat$yc <- 1 + 0.5 * dat$z + 0.7 * dat$x1 - 0.4 * dat$x2 + rnorm(nrow(dat))
  dat
}

# ATE weights shared by the outcome-model fixtures below.
se_method_ate_wts <- function(dat, ps_mod) {
  ps <- predict(ps_mod, type = "response")
  wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)
}

# Marginal (exposure-only) weighted linear outcome model.
se_method_outcome_marginal_lm <- function(dat, ps_mod) {
  wts <- se_method_ate_wts(dat, ps_mod)
  lm(yc ~ z, data = dat, weights = wts)
}

# Covariate-adjusted weighted linear outcome model. Its right-hand side carries
# a covariate beyond the exposure, so the linearization influence functions no
# longer match the g-computation estimator that is reported.
se_method_outcome_adjusted_lm <- function(dat, ps_mod) {
  wts <- se_method_ate_wts(dat, ps_mod)
  lm(yc ~ z + x1, data = dat, weights = wts)
}

# Covariate-adjusted weighted quasibinomial outcome model.
se_method_outcome_adjusted_glm <- function(dat, ps_mod) {
  wts <- se_method_ate_wts(dat, ps_mod)
  glm(y ~ z + x1, data = dat, family = quasibinomial(), weights = wts)
}

# Exposure-by-covariate interaction outcome model. Its expansion carries a bare
# covariate term (`x1`) alongside the exposure and interaction, so the model is
# adjusted.
se_method_outcome_interaction_lm <- function(dat, ps_mod) {
  wts <- se_method_ate_wts(dat, ps_mod)
  lm(yc ~ z * x1, data = dat, weights = wts)
}

# Adversarial adjusted model whose every right-hand-side term involves the
# exposure (`z` and `z:x1`) yet is still covariate-adjusted. A guard that only
# rejects terms not involving the exposure would wrongly admit this model.
se_method_outcome_exposure_interaction_lm <- function(dat, ps_mod) {
  wts <- se_method_ate_wts(dat, ps_mod)
  lm(yc ~ z + z:x1, data = dat, weights = wts)
}

# Exposure-only outcome model fit without an intercept. Its term labels are the
# exposure alone, so it passes a guard that only compares term labels, but the
# fit carries no reference level: the linear predictor at z = 0 is 0 for every
# unit.
se_method_outcome_no_intercept <- function(dat, ps_mod) {
  wts <- se_method_ate_wts(dat, ps_mod)
  glm(y ~ z - 1, data = dat, family = quasibinomial(), weights = wts)
}

# The same omission in a linear outcome model, whose mean under no exposure is
# pinned at 0 rather than 0.5.
se_method_outcome_no_intercept_lm <- function(dat, ps_mod) {
  wts <- se_method_ate_wts(dat, ps_mod)
  lm(yc ~ z - 1, data = dat, weights = wts)
}

test_that("ipw() defaults to the mestimation SE method and returns a fit", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  res <- ipw(ps_mod, outcome_mod, .data = dat)

  expect_equal(res$se_method, "mestimation")

  expect_false(is.null(res$fit))
  expect_false(is.null(stats::coef(res$fit)))
  expect_false(is.null(stats::vcov(res$fit)))

  # The estimates table keeps the shared eight-column contract.
  expect_named(res$estimates, estimates_columns)
})

test_that("ipw(se_method = 'linearization') matches the pre-change output", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")

  expect_equal(res$se_method, "linearization")
  expect_null(res$fit)

  reference <- se_method_lin_reference()
  expect_equal(res$estimates$effect, reference$effect)
  expect_equal(res$estimates$estimate, reference$estimate, tolerance = 1e-8)
  expect_equal(res$estimates$std.err, reference$std.err, tolerance = 1e-8)
})

test_that("the two SE methods share point estimates but differ in SEs", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  res_m <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "mestimation")
  res_l <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")

  # Point estimates are identical across methods.
  expect_equal(
    res_m$estimates$estimate,
    res_l$estimates$estimate,
    tolerance = 1e-8
  )

  # Standard errors are close but not identical: the two methods genuinely
  # differ, and the relative gap stays within the engine tests' 3 percent band.
  rel_diff <- abs(res_m$estimates$std.err - res_l$estimates$std.err) /
    res_l$estimates$std.err
  expect_true(all(rel_diff < 0.03))
  expect_false(isTRUE(all.equal(
    res_m$estimates$std.err,
    res_l$estimates$std.err
  )))
})

test_that("an invalid se_method value errors", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "not-a-method"),
    class = "rlang_error"
  )
})

test_that("trimmed weights error on both SE paths", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_trimmed(dat, ps_mod)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat),
    class = "propensity_ipw_trimmed_error",
    regexp = "trim"
  )
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_ipw_trimmed_error",
    regexp = "trim"
  )
})

test_that("truncated weights error on both SE paths", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_truncated(dat, ps_mod)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat),
    class = "propensity_ipw_truncated_error",
    regexp = "truncat"
  )
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_ipw_truncated_error",
    regexp = "truncat"
  )
})

test_that("calibrated weights error on both SE paths", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_calibrated(dat, ps_mod)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat),
    class = "propensity_ipw_calibrated_error",
    regexp = "calibrat"
  )
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_ipw_calibrated_error",
    regexp = "calibrat"
  )
})

test_that("an outcome model with an offset term errors before estimation", {
  dat <- se_method_data()
  dat$off <- 0.3 * dat$x1
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_offset(dat, ps_mod)

  # The default (mestimation) path threads no offset through the outcome-score
  # block, so it must reject an offset outright, before any estimation and with
  # no message or other condition signaled first.
  expect_no_message(
    expect_error(
      ipw(ps_mod, outcome_mod, .data = dat),
      class = "propensity_ipw_offset_error",
      regexp = "offset"
    )
  )
})

test_that("the modified-weight guard fires before estimand parsing", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_trimmed(dat, ps_mod)

  # Even with a deliberately wrong estimand argument, the trim guard must win:
  # the error class is the specific guard, not a generic estimand mismatch.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, estimand = "att"),
    class = "propensity_ipw_trimmed_error"
  )
})

test_that("mismatched weights error on both standard error paths", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)

  # Doubled weights: inconsistent with the propensity score model, but the psw
  # estimand attribute is still "ate" so estimand detection succeeds.
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts + wts
  )

  # Both preflights recompute the weights from the propensity score model and
  # detect the mismatch. The linearization path once accepted this call: it
  # predicts the propensity scores and takes the weights from the outcome model
  # without requiring the two to agree, so doubled weights changed the standard
  # errors and nothing said so.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "mestimation"),
    class = "propensity_ipw_weights_mismatch_error"
  )

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_ipw_weights_mismatch_error"
  )
})

test_that("a case-weighted binary propensity model errors on the mestimation path", {
  dat <- se_method_data()
  case_weights <- rep(c(1, 2), length.out = nrow(dat))
  ps_mod <- glm(
    z ~ x1 + x2,
    data = dat,
    family = binomial(),
    weights = case_weights
  )
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  # The stacked propensity score block is unweighted, so a case-weighted binary
  # ps model would not sit at the score root and the estimates would drift. The
  # mestimation path rejects it, matching the categorical and continuous paths.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "mestimation"),
    class = "propensity_ipw_ps_weights_error"
  )

  # The linearization path treats the weights as fixed and does not restack the
  # propensity score model, so a case-weighted ps model is not rejected there.
  expect_no_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
  )
})

test_that("stabilized ATE weights work end to end on the default path", {
  withr::local_options(propensity.quiet = TRUE)
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(
    ps,
    dat$z,
    exposure_type = "binary",
    .focal_level = 1,
    stabilize = TRUE
  )
  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  res <- ipw(ps_mod, outcome_mod, .data = dat)

  expect_equal(res$se_method, "mestimation")
  expect_false(is.null(res$fit))
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("as.data.frame(exponentiate = TRUE) matches across SE methods", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  res_m <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "mestimation")
  res_l <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")

  df_m <- as.data.frame(res_m, exponentiate = TRUE)
  df_l <- as.data.frame(res_l, exponentiate = TRUE)

  # Both paths relabel the ratio rows and share point estimates on the natural
  # scale (standard errors stay on the log scale and may differ by method).
  expect_equal(df_m$effect, c("rd", "rr", "or"))
  expect_equal(df_l$effect, c("rd", "rr", "or"))
  expect_equal(df_m$estimate, df_l$estimate, tolerance = 1e-8)
})

test_that("ipw() print output is stable per SE method", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  # Pins the default-path print output.
  expect_snapshot(print(ipw(ps_mod, outcome_mod, .data = dat)))

  # Pins the explicit linearization print output.
  expect_snapshot(
    print(ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"))
  )
})

test_that("a marginal linear outcome model works on the linearization path", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_marginal_lm(dat, ps_mod)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")

  expect_s3_class(res, "ipw")
  expect_equal(res$se_method, "linearization")
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("a marginal quasibinomial outcome model works on the linearization path", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")

  expect_s3_class(res, "ipw")
  expect_equal(res$se_method, "linearization")
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("a covariate-adjusted outcome model still works on the mestimation path", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_adjusted_lm(dat, ps_mod)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "mestimation")

  expect_s3_class(res, "ipw")
  expect_equal(res$se_method, "mestimation")
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("a covariate-adjusted linear outcome model errors on the linearization path", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_adjusted_lm(dat, ps_mod)

  expect_error(
    ipw(ps_mod, outcome_mod, se_method = "linearization"),
    class = "propensity_method_error",
    regexp = "mestimation"
  )
})

test_that("a covariate-adjusted linear outcome model errors on the linearization path with .data", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_adjusted_lm(dat, ps_mod)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error",
    regexp = "mestimation"
  )
})

test_that("a covariate-adjusted quasibinomial outcome model errors on the linearization path", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_adjusted_glm(dat, ps_mod)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error",
    regexp = "mestimation"
  )
})

test_that("an exposure-covariate interaction outcome model errors on the linearization path", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_interaction_lm(dat, ps_mod)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error",
    regexp = "mestimation"
  )
})

test_that("an outcome model whose every term involves the exposure still errors when adjusted", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_exposure_interaction_lm(dat, ps_mod)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error",
    regexp = "mestimation"
  )
})

test_that("a formula-offset outcome model errors on the linearization path", {
  dat <- se_method_data()
  dat$off <- 0.3 * dat$x1
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_offset(dat, ps_mod)

  # The offset shifts the g-computation point estimates while the linearization
  # influence functions know nothing about it, so reject it before estimating.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_ipw_offset_error",
    regexp = "offset"
  )
})

test_that("an offset-argument outcome model errors on the linearization path", {
  dat <- se_method_data()
  dat$off <- 0.3 * dat$x1
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_offset_arg(dat, ps_mod)

  # An offset supplied through the `offset` argument is invisible to the term
  # labels, so the guard must inspect the model frame to catch this route.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_ipw_offset_error",
    regexp = "offset"
  )
})

# ---- no-intercept outcome model on the linearization path -------------------
#
# `y ~ z - 1` has term labels equal to the exposure name, so the exposure-only
# guard admits it, but the model has no reference level: the linear predictor at
# z = 0 is 0 for every unit, so the g-computation mu0 is plogis(0) = 0.5 whatever
# the data say. The weighted score no longer ties mu0 to the Hajek mean the
# influence functions assume, and the reported effects drift far from the Hajek
# contrast while the standard errors stay narrow, with nothing signaled. The
# intercept-bearing counterpart is covered above by "a marginal quasibinomial
# outcome model works on the linearization path", which must keep passing.

test_that("a no-intercept outcome model errors on the linearization path", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_no_intercept(dat, ps_mod)

  # The message must name the missing intercept, which distinguishes this guard
  # from the covariate-adjusted rejection that shares the error class.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error",
    regexp = "[Ii]ntercept"
  )

  # and it must direct the user to the SE method that handles this model.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    regexp = "mestimation"
  )
})

test_that("a no-intercept outcome model errors on the linearization path without .data", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_no_intercept(dat, ps_mod)

  # The guard runs before either extraction route, so omitting .data must reject
  # the model just the same.
  expect_error(
    ipw(ps_mod, outcome_mod, se_method = "linearization"),
    class = "propensity_method_error"
  )
})

# Factor-coded exposure with the outcome model either intercept-bearing or
# saturated. The saturated form is the case the roxygen's "including the
# saturated factor codings" clause rests on. Local to this file: test files do
# not share each other's top-level definitions, and test-ipw-mestimation.R's
# equivalent fixture is not reachable here.
se_method_factor_models <- function(dat, intercept = TRUE) {
  dat$zf <- factor(dat$z, levels = c(0, 1), labels = c("no", "yes"))
  ps_mod <- glm(zf ~ x1 + x2, data = dat, family = binomial())
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_mod))
  outcome_mod <- glm(
    if (intercept) y ~ zf else y ~ 0 + zf,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, dat = dat)
}

test_that("a saturated factor no-intercept outcome model errors on the linearization path", {
  # The documented contract is that EVERY no-intercept outcome model is rejected
  # here, including the saturated factor codings the M-estimation path accepts.
  # The tests above use `y ~ z - 1`, whose cell mean under no exposure really is
  # pinned at the link's zero point. `y ~ 0 + zf` is the case the claim rests on
  # and the one that was unpinned: it estimates both cell means, is a
  # reparameterization of the with-intercept fit, and is still refused, because
  # the influence functions are derived for the intercept parameterization.
  dat <- se_method_data()
  mods <- se_method_factor_models(dat, intercept = FALSE)

  # the fixture is saturated rather than zero-pinned: one coefficient per level
  expect_equal(
    names(coef(mods$outcome_mod)),
    c("zfno", "zfyes")
  )

  expect_error(
    ipw(
      mods$ps_mod,
      mods$outcome_mod,
      .data = mods$dat,
      se_method = "linearization"
    ),
    class = "propensity_method_error",
    regexp = "[Ii]ntercept"
  )

  # and the with-intercept counterpart, the same fit reparameterized, runs
  with_intercept <- se_method_factor_models(dat, intercept = TRUE)
  expect_s3_class(
    ipw(
      with_intercept$ps_mod,
      with_intercept$outcome_mod,
      .data = with_intercept$dat,
      se_method = "linearization"
    ),
    "ipw"
  )
})

test_that("a no-intercept linear outcome model errors on the linearization path", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_no_intercept_lm(dat, ps_mod)

  # The guard reads the terms object, not the family, so a linear outcome model
  # is rejected on the same grounds.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error"
  )
})

test_that("the no-intercept rejection names the intercept and the SE method", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_no_intercept(dat, ps_mod)

  expect_snapshot(
    error = TRUE,
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
  )
})

# ---- probit and cloglog linearization variance ------------------------------
#
# The linearization ps correction adds H' A^{-1} s_i to the influence values,
# where s_i = x_i (Z_i - p_i) / (p_i (1 - p_i) g'(p_i)) is the score of a
# binomial GLM with link g. The link enters the weight derivatives, the Fisher
# information A, and the score itself: the canonical logit form
# x_i (Z_i - p_i) omits the factor 1 / (p (1 - p) g'(p)), which equals 1 only
# for logit. These reference tests pin the generalized score correction, so the
# probit and cloglog cases are the discriminating ones.

# g'(p) for the propensity score link, matching derive_link() in R/ipw.R.
lin_gprime <- function(p, link) {
  switch(
    link,
    logit = 1 / (p * (1 - p)),
    probit = 1 / dnorm(qnorm(p)),
    cloglog = -1 / ((1 - p) * log(1 - p))
  )
}

# Hand-coded generalized linearization RD standard error. `if_scale` is the score
# factor 1 / (p (1 - p) g'(p)) that turns x_i (Z_i - p_i) into the GLM score, and
# `info_scale` is the same factor entering the Fisher information A. The correct
# estimator uses the score factor in both places; passing if_scale = 1 reproduces
# the canonical-logit form, which differs from the correct value for probit and
# cloglog.
lin_se_rd <- function(X, p, z, y, w, dwdeta, if_scale, info_scale) {
  n <- length(z)
  n1 <- sum(w[z == 1])
  n0 <- sum(w[z == 0])
  mu1 <- weighted.mean(y[z == 1], w[z == 1])
  mu0 <- weighted.mean(y[z == 0], w[z == 0])
  l1u <- n / n1 * (w * z * (y - mu1))
  l0u <- n / n0 * (w * (1 - z) * (y - mu0))
  h1 <- colSums(X * (dwdeta * z * (y - mu1))) / n
  h0 <- colSums(X * (dwdeta * (1 - z) * (y - mu0))) / n
  info <- crossprod(X * sqrt(info_scale^2 * p * (1 - p))) / n
  if_beta <- t(solve(info, t(X * ((z - p) * if_scale))))
  l1 <- l1u + n / n1 * drop(if_beta %*% h1)
  l0 <- l0u + n / n0 * drop(if_beta %*% h0)
  sqrt(var(l1 - l0) / n)
}

# Design pieces for the RD linearization: the ps design and fitted scores, the
# score factor, and the ATE weight derivative dw/deta (treated -1/(p^2 g'(p)),
# untreated 1/((1 - p)^2 g'(p))).
lin_ate_pieces <- function(ps_mod, z) {
  p <- as.double(predict(ps_mod, type = "response"))
  link <- ps_mod$family$link
  gp <- lin_gprime(p, link)
  list(
    X = model.matrix(ps_mod),
    p = p,
    score_factor = 1 / (p * (1 - p) * gp),
    dwdeta = ifelse(z == 1, -1 / (p^2 * gp), 1 / ((1 - p)^2 * gp))
  )
}

# Correct RD linearization SE (score factor in both the influence and the
# information) from a fitted ps model and the outcome data. `score` is the
# per-observation stabilization score carried on the weights: the weights are
# base * score, so the score scales the weight derivative observation-wise while
# leaving the propensity score's own score factor and information untouched. The
# default of 1 leaves unstabilized analyses unchanged.
lin_se_correct <- function(ps_mod, z, y, w, score = 1) {
  pieces <- lin_ate_pieces(ps_mod, z)
  lin_se_rd(
    pieces$X,
    pieces$p,
    z,
    y,
    as.double(w),
    pieces$dwdeta * score,
    if_scale = pieces$score_factor,
    info_scale = pieces$score_factor
  )
}

# Seeded data whose exposure is generated on the requested link, so the fitted
# propensity score model is well specified. The confounding is moderate, so the
# propensity score fits converge without separation warnings and the ATE weights
# stay finite and stable for both probit and cloglog.
se_link_data <- function(link, seed = 4242, n = 4000) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  lp <- 0.7 * x1 - 0.5 * x2
  pz <- switch(
    link,
    logit = plogis(lp),
    probit = pnorm(lp),
    cloglog = 1 - exp(-exp(lp))
  )
  z <- rbinom(n, 1, pz)
  y <- rbinom(n, 1, plogis(-0.3 + 0.6 * z + 0.5 * x1 - 0.4 * x2))
  data.frame(x1, x2, z, y)
}

# Propensity score model on the requested link and a marginal quasibinomial
# outcome model weighted with ATE weights (tightened IRLS so the fit sits at the
# weighted MLE).
se_link_models <- function(dat, link) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial(link = link))
  ps <- as.double(predict(ps_mod, type = "response"))
  wts <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

test_that("the hand-coded linearization RD SE reproduces the logit path and reference", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)
  wts <- wt_ate(
    predict(ps_mod, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )

  # For a logit ps model the score factor is 1, so the generalized helper and the
  # package agree, validating the helper against the package implementation.
  hand <- lin_se_correct(ps_mod, dat$z, dat$y, wts)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
  rd_se <- res$estimates$std.err[res$estimates$effect == "rd"]
  expect_equal(hand, rd_se, tolerance = 1e-8)

  # and against the pinned linearization reference for this fixture
  reference <- se_method_lin_reference()
  expect_equal(
    hand,
    reference$std.err[reference$effect == "rd"],
    tolerance = 1e-8
  )
})

test_that("probit linearization RD SE matches the generalized score correction", {
  dat <- se_link_data("probit")
  mods <- se_link_models(dat, "probit")

  correct <- lin_se_correct(mods$ps_mod, dat$z, dat$y, mods$wts)
  res <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = dat,
    se_method = "linearization"
  )
  rd_se <- res$estimates$std.err[res$estimates$effect == "rd"]
  expect_equal(rd_se, correct, tolerance = 1e-8)
})

test_that("cloglog linearization RD SE matches the generalized score correction", {
  dat <- se_link_data("cloglog")
  mods <- se_link_models(dat, "cloglog")

  correct <- lin_se_correct(mods$ps_mod, dat$z, dat$y, mods$wts)
  res <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = dat,
    se_method = "linearization"
  )
  rd_se <- res$estimates$std.err[res$estimates$effect == "rd"]
  expect_equal(rd_se, correct, tolerance = 1e-8)
})

test_that("probit linearization RD SE agrees with mestimation", {
  dat <- se_link_data("probit")
  mods <- se_link_models(dat, "probit")

  lin <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = dat,
    se_method = "linearization"
  )$estimates
  mest <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = dat,
    se_method = "mestimation"
  )$estimates

  lin_se <- lin$std.err[lin$effect == "rd"]
  mest_se <- mest$std.err[mest$effect == "rd"]

  # The score-corrected linearization agrees with the sandwich to well under a
  # percent (about 0.02 percent here); dropping the score factor puts it about
  # 1.6 percent off, so this band separates the two.
  expect_lt(abs(lin_se - mest_se) / mest_se, 0.005)
})

# ---- outcome-family validation on the linearization path --------------------

test_that("the linearization path rejects an unsupported outcome family", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)

  # A low-rate count outcome; the marginal outcome model passes the exposure-only
  # guard, so the family rejection must fire even on a marginal model. Without it
  # the path returns a finite but meaningless log(or) from count-scale marginal
  # means.
  withr::local_seed(1)
  dat$ycount <- rpois(nrow(dat), exp(-0.7 + 0.3 * dat$z + 0.2 * dat$x1))
  outcome_mod <- suppressWarnings(
    glm(ycount ~ z, data = dat, family = quasipoisson(), weights = wts)
  )

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_ipw_family_error"
  )
})

# ---- factor outcome response on the linearization path ----------------------

test_that("linearization matches a factor outcome response with finite SEs", {
  dat <- se_method_data()
  dat$yf <- factor(ifelse(dat$y == 1, "yes", "no"), levels = c("no", "yes"))
  ps_mod <- se_method_ps_mod(dat)
  wts <- wt_ate(
    predict(ps_mod, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  num <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
  fac <- suppressWarnings(
    glm(yf ~ z, data = dat, family = quasibinomial(), weights = wts)
  )

  res_num <- ipw(
    ps_mod,
    num,
    .data = dat,
    se_method = "linearization"
  )$estimates

  # Y - mu must be computed on the 0/1 coding: on a raw factor it gives NA
  # std.err with "'-' not meaningful for factors" warnings. The linearization
  # branch extracts the outcome two ways that would fail identically:
  # `.data[[outcome_name]]` when .data is supplied, and fmla_extract_left_vctr
  # when it is not. Both return finite std.err equal to the numeric fit with no
  # warnings.
  res_data <- expect_no_warning(
    ipw(ps_mod, fac, .data = dat, se_method = "linearization")$estimates
  )
  res_nodata <- expect_no_warning(
    ipw(ps_mod, fac, se_method = "linearization")$estimates
  )

  expect_true(all(is.finite(res_data$std.err)))
  expect_true(all(is.finite(res_nodata$std.err)))
  expect_equal(res_data$estimate, res_num$estimate, tolerance = 1e-8)
  expect_equal(res_data$std.err, res_num$std.err, tolerance = 1e-8)
  expect_equal(res_nodata$estimate, res_num$estimate, tolerance = 1e-8)
  expect_equal(res_nodata$std.err, res_num$std.err, tolerance = 1e-8)
})

# ---- factor exposures on the linearization path -----------------------------
#
# The mestimation path recodes the exposure to 0/1 via
# as.double(exposure == exposure_values[[2]]) (the second sorted level), but the
# linearization path passes the raw extracted exposure into
# linearize_variables_for_wts (wts * Z, 1 - Z), ate_derivative (exposure == 1),
# and correct_for_ps ((exposure - ps)). A factor exposure then triggers
# "not meaningful for factors" warnings and, in these fixtures, a hard error
# ("'x' must be an array of at least two dimensions"); exposure == 1 is also
# FALSE for every unit, so the treated/control derivative branches would be
# misassigned. These tests pin parity with the numeric-0/1 analysis, isolating
# the exposure recoding: the factor and numeric arms share the same fitted
# propensity scores and the same weight object, differing only in whether the
# exposure is a factor or its 0/1 recode.

factor_exposure_data <- function(seed = 2024, n = 800) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x1))
  y <- rbinom(n, 1, plogis(-0.4 + 1.0 * z + 0.5 * x1))
  trt <- factor(
    ifelse(z == 1, "treated", "control"),
    levels = c("control", "treated")
  )
  trt_rev <- factor(
    ifelse(z == 1, "treated", "control"),
    levels = c("treated", "control")
  )
  data.frame(x1, z, y, trt, trt_rev)
}

# Build a factor-exposure arm and its numeric-0/1 counterpart. The 0/1 recode is
# the factor's second sorted level (the mestimation focal, coded-1, level), so
# the ps models share fitted values and the weights (built once from the factor
# with a matching focal level) are reused by both outcome models. `wt_fun`
# selects the estimand.
factor_numeric_arms <- function(dat, factor_col, wt_fun) {
  focal <- levels(dat[[factor_col]])[[2]]
  d <- dat
  d$znum <- as.double(d[[factor_col]] == focal)
  ps_fac <- glm(
    stats::reformulate("x1", factor_col),
    data = d,
    family = binomial()
  )
  ps_num <- glm(znum ~ x1, data = d, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_fun(
      predict(ps_fac, type = "response"),
      d[[factor_col]],
      exposure_type = "binary",
      .focal_level = focal
    )
  )
  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)
  om_fac <- suppressWarnings(glm(
    stats::reformulate(factor_col, "y"),
    data = d,
    family = quasibinomial(),
    weights = wts,
    control = ctrl
  ))
  om_num <- glm(
    y ~ znum,
    data = d,
    family = quasibinomial(),
    weights = wts,
    control = ctrl
  )
  list(
    ps_fac = ps_fac,
    om_fac = om_fac,
    ps_num = ps_num,
    om_num = om_num,
    dat = d
  )
}

test_that("linearization recodes a factor exposure to match the numeric fit", {
  arms <- factor_numeric_arms(factor_exposure_data(), "trt", wt_ate)
  ref <- ipw(arms$ps_num, arms$om_num, se_method = "linearization")$estimates

  # Both linearization extraction routes must recode the factor: no .data
  # (fmla_extract_left_vctr) and .data (the .data column).
  res_nodata <- expect_no_warning(
    ipw(arms$ps_fac, arms$om_fac, se_method = "linearization")$estimates
  )
  res_data <- expect_no_warning(
    ipw(
      arms$ps_fac,
      arms$om_fac,
      .data = arms$dat,
      se_method = "linearization"
    )$estimates
  )
  expect_true(all(is.finite(res_nodata$std.err)))
  expect_true(all(is.finite(res_data$std.err)))
  expect_equal(res_nodata, ref, tolerance = 1e-8)
  expect_equal(res_data, ref, tolerance = 1e-8)
})

test_that("linearization recodes a reversed-level factor exposure to match the numeric fit", {
  # Levels c("treated", "control"): the second sorted level is "control", so the
  # recode codes control as 1. A fix hardcoding a label or level order fails.
  arms <- factor_numeric_arms(factor_exposure_data(), "trt_rev", wt_ate)
  ref <- ipw(arms$ps_num, arms$om_num, se_method = "linearization")$estimates
  res <- expect_no_warning(
    ipw(arms$ps_fac, arms$om_fac, se_method = "linearization")$estimates
  )
  expect_true(all(is.finite(res$std.err)))
  expect_equal(res, ref, tolerance = 1e-8)
})

test_that("linearization recodes a factor exposure for the att estimand", {
  # att_derivative branches on exposure == 1, which is FALSE for every unit of a
  # factor exposure, so the recode must precede the derivative as well.
  arms <- factor_numeric_arms(factor_exposure_data(), "trt", wt_att)
  ref <- ipw(arms$ps_num, arms$om_num, se_method = "linearization")$estimates
  res <- expect_no_warning(
    ipw(arms$ps_fac, arms$om_fac, se_method = "linearization")$estimates
  )
  expect_true(all(is.finite(res$std.err)))
  expect_equal(res, ref, tolerance = 1e-8)
})

test_that("mestimation matches a factor exposure to the numeric fit", {
  skip_if_not_installed("deli")
  for (factor_col in c("trt", "trt_rev")) {
    arms <- factor_numeric_arms(factor_exposure_data(), factor_col, wt_ate)
    res_fac <- ipw(
      arms$ps_fac,
      arms$om_fac,
      se_method = "mestimation"
    )$estimates
    res_num <- ipw(
      arms$ps_num,
      arms$om_num,
      se_method = "mestimation"
    )$estimates
    expect_equal(res_fac, res_num, tolerance = 1e-8)
  }
})

# ---- atu and entropy rejection on the linearization path --------------------
#
# The linearization path supports only ate, att, ato, and atm. An atu or entropy
# request has to be rejected before derive_weights, whose rlang::arg_match
# raises a bare, unclassed rlang error attributed to the internal function and
# claiming the estimand "must be one of ate, att, ato, atm" (misleading: atu and
# entropy are valid estimands, just unsupported by this SE method). These tests
# pin a classed propensity_method_error directing to se_method = "mestimation",
# matching the categorical and continuous paths. The four supported estimands on
# linearization are covered by test-ipw.R ("variance works"); atu and entropy on
# mestimation by test-ipw-mestimation.R ("ipw_mestimation runs atu and entropy
# with finite SEs").

test_that("linearization rejects the atu estimand with a classed error", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_atu(ps_mod))
  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error",
    regexp = "mestimation"
  )
})

test_that("linearization rejects the entropy estimand with a classed error", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_entropy(ps_mod))
  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error",
    regexp = "mestimation"
  )
})

test_that("the linearization atu rejection names the SE method", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_atu(ps_mod))
  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  expect_snapshot(
    error = TRUE,
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
  )
})

# ---- more-than-two-level exposure on the linearization path ------------------

test_that("a more-than-two-level exposure on the linearization path aborts informatively", {
  # The linearization path routes the extracted exposure through
  # estimate_marginal_means, whose two-level abort places call = call inside the
  # c() message vector instead of passing it to abort, so a >2-level exposure
  # crashes with rlang's "message must be a character vector, not a list" rather
  # than the intended informative binary-only error. A .data column with three
  # exposure values reaches the abort.
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)
  dat3 <- dat
  dat3$z <- rep_len(0:2, nrow(dat))

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat3, se_method = "linearization"),
    class = "propensity_error",
    regexp = "binary"
  )
})

# ---- matrix-response (cbind) propensity model -------------------------------
#
# A ps model with a matrix LHS (cbind(successes, failures) binomial response) is
# not a binary exposure, and both SE paths reject it with the same error naming
# the matrix response. Left to their downstream guards the two diverge and
# mislead: the mestimation path reports case weights (the cbind totals become
# prior weights) and the linearization path reports an adjusted outcome model
# (the mangled length-3 exposure name c("cbind", "y1", "y0") never matches the
# outcome term). "matrix" is absent from both of those messages ("cbind" is not,
# since it appears as the mangled exposure name on the linearization path), so
# it is the discriminating regexp.

matrix_lhs_models <- function() {
  set.seed(2024)
  n <- 300
  x1 <- rnorm(n)
  y1 <- rbinom(n, 5, plogis(0.3 * x1))
  y0 <- 5 - y1
  z <- rbinom(n, 1, plogis(0.3 * x1))
  y <- rbinom(n, 1, plogis(-0.4 + z))
  dat <- data.frame(x1, y1, y0, z, y)
  ps_mod <- glm(cbind(y1, y0) ~ x1, data = dat, family = binomial())
  ps_bin <- glm(z ~ x1, data = dat, family = binomial())
  wts <- wt_ate(
    predict(ps_bin, type = "response"),
    z,
    exposure_type = "binary",
    .focal_level = 1
  )
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, dat = dat)
}

test_that("a matrix-response propensity model errors consistently on both SE paths", {
  m <- matrix_lhs_models()
  for (se in c("mestimation", "linearization")) {
    expect_error(
      ipw(m$ps_mod, m$outcome_mod, .data = m$dat, se_method = se),
      class = "propensity_error",
      regexp = "matrix"
    )
  }
})

test_that("the matrix-response propensity model error names the matrix response", {
  m <- matrix_lhs_models()
  expect_snapshot(
    error = TRUE,
    ipw(m$ps_mod, m$outcome_mod, .data = m$dat)
  )
})

test_that("the matrix-response propensity guard carries the response class", {
  # A matrix response is a statement about the shape of a model's left-hand
  # side, not about the semantics of the exposure. The guards that are about
  # exposure semantics keep propensity_ipw_exposure_error, so the two families
  # can be caught apart; see "mestimation rejects a binary outcome model that
  # omits the exposure without .data" in test-ipw-mestimation.R for one of them.
  m <- matrix_lhs_models()
  for (se in c("mestimation", "linearization")) {
    expect_error(
      ipw(m$ps_mod, m$outcome_mod, .data = m$dat, se_method = se),
      class = "propensity_ipw_response_error"
    )
  }
})

# ---- matrix-response (cbind) outcome model ----------------------------------
#
# The outcome model's counterpart, and unguarded. A cbind(successes, failures)
# response is twice the length of the exposure once the model frame flattens it,
# and what the user is told depends on the route: without .data the length
# reconciliation reports an outcome of length 600 against an exposure of length
# 300, and with .data the extraction asks for a column named "cbind", which does
# not exist and never will. Neither says the response is a matrix, and the
# second actively sends the user to add a column named after a function.
#
# The single-column case is what every other test in this file fits, so it needs
# no pin of its own here.

matrix_outcome_models <- function() {
  set.seed(2024)
  n <- 300
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x1))
  y1 <- rbinom(n, 5, plogis(-0.4 + z))
  y0 <- 5 - y1
  dat <- data.frame(x1, z, y1, y0)
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      predict(ps_mod, type = "response"),
      dat$z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  outcome_mod <- glm(
    cbind(y1, y0) ~ z,
    data = dat,
    family = binomial(),
    weights = wts
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, dat = dat)
}

test_that("a matrix-response outcome model errors on its shape without .data", {
  skip_if_not_installed("deli")
  m <- matrix_outcome_models()

  # the fixture is the shape it claims to be
  expect_true(is.matrix(stats::model.response(stats::model.frame(
    m$outcome_mod
  ))))

  for (se in c("mestimation", "linearization")) {
    err <- expect_error(
      ipw(m$ps_mod, m$outcome_mod, se_method = se),
      class = "propensity_ipw_response_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "matrix response", fixed = TRUE)
  }
})

test_that("a transformed single-column outcome response is not a matrix response", {
  skip_if_not_installed("deli")
  # `log(y)` and friends make the formula's left-hand side a call rather than a
  # symbol, which is what the propensity guard rejects and what an early version
  # of this guard rejected too. They are ordinary single-column responses and
  # the analysis is correct, so the only honest check is the built frame. The
  # reference fits the same model on a precomputed column, which is the same
  # model written differently and must give the same answer.
  set.seed(2024)
  n <- 300
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x1))
  y <- rlnorm(n, -0.4 + z + 0.3 * x1)
  dat <- data.frame(x1, z, y, log_y = log(y), sqrt_y = sqrt(y))
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      predict(ps_mod, type = "response"),
      dat$z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )

  for (se in c("mestimation", "linearization")) {
    for (pair in list(c("log(y)", "log_y"), c("sqrt(y)", "sqrt_y"))) {
      transformed <- lm(
        stats::reformulate("z", response = pair[[1]]),
        data = dat,
        weights = wts
      )
      precomputed <- lm(
        stats::reformulate("z", response = pair[[2]]),
        data = dat,
        weights = wts
      )

      expect_false(
        is.matrix(stats::model.response(stats::model.frame(transformed)))
      )
      expect_equal(
        ipw(ps_mod, transformed, se_method = se)$estimates$estimate,
        ipw(ps_mod, precomputed, se_method = se)$estimates$estimate,
        tolerance = 1e-10
      )
    }
  }
})

test_that("the matrix-response outcome model error names the outcome model", {
  skip_if_not_installed("deli")
  m <- matrix_outcome_models()

  expect_snapshot(
    error = TRUE,
    ipw(m$ps_mod, m$outcome_mod)
  )
})

test_that("a matrix-response outcome model errors on its shape with .data", {
  skip_if_not_installed("deli")
  # The .data route fails differently today, asking for a "cbind" column, so it
  # is pinned separately rather than folded into the test above.
  m <- matrix_outcome_models()

  for (se in c("mestimation", "linearization")) {
    err <- expect_error(
      ipw(m$ps_mod, m$outcome_mod, .data = m$dat, se_method = se),
      class = "propensity_ipw_response_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "matrix response", fixed = TRUE)
  }
})

# ---- per-observation stabilization scores on the linearization path ----------

# Marginal quasibinomial outcome model weighted with the supplied weights, fit
# with tightened IRLS so the fit sits at the weighted MLE and the hand-coded
# oracle applies exactly.
se_score_outcome <- function(dat, wts) {
  glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
}

# Propensity score model, score-carrying stabilized ATE weights, and the marginal
# outcome model they weight.
se_score_models <- function(dat, score) {
  ps_mod <- se_method_ps_mod(dat)
  ps <- as.double(predict(ps_mod, type = "response"))
  wts <- wt_ate(
    ps,
    dat$z,
    exposure_type = "binary",
    .focal_level = 1,
    stabilize = TRUE,
    stabilization_score = score
  )
  list(
    ps_mod = ps_mod,
    outcome_mod = se_score_outcome(dat, wts),
    wts = wts
  )
}

# RD standard error from the linearization path for the supplied weights.
se_score_lin_rd <- function(dat, ps_mod, wts) {
  res <- ipw(
    ps_mod,
    se_score_outcome(dat, wts),
    .data = dat,
    se_method = "linearization"
  )
  res$estimates$std.err[res$estimates$effect == "rd"]
}

test_that("linearization threads a covariate-correlated stabilization score through the weight derivatives", {
  dat <- se_method_data()
  # A score correlated with a propensity score covariate, so omitting it from the
  # weight derivatives misstates the propensity score correction rather than
  # cancelling out of it.
  score <- exp(0.5 * dat$x1)
  mods <- se_score_models(dat, score)

  correct <- lin_se_correct(
    mods$ps_mod,
    dat$z,
    dat$y,
    mods$wts,
    score = score
  )
  res <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = dat,
    se_method = "linearization"
  )
  rd_se <- res$estimates$std.err[res$estimates$effect == "rd"]
  expect_equal(rd_se, correct, tolerance = 1e-8)
})

test_that("linearization matches the hand-coded SE for a unit stabilization score", {
  dat <- se_method_data()
  # A recorded score of 1 leaves the weights and their derivatives at their
  # unstabilized values, so threading the score must not move this analysis.
  score <- rep(1, nrow(dat))
  mods <- se_score_models(dat, score)

  correct <- lin_se_correct(
    mods$ps_mod,
    dat$z,
    dat$y,
    mods$wts,
    score = score
  )
  res <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = dat,
    se_method = "linearization"
  )
  rd_se <- res$estimates$std.err[res$estimates$effect == "rd"]
  expect_equal(rd_se, correct, tolerance = 1e-8)
})

test_that("linearization default stabilization reproduces the unstabilized SE", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- as.double(predict(ps_mod, type = "response"))

  # The default stabilizer scales the treated weights by mean(z) and the
  # untreated weights by 1 - mean(z). Each group constant scales that group's
  # weights and that group's weight total identically, so it cancels out of the
  # Hajek ratio and out of the propensity score correction alike: the stabilized
  # analysis must report exactly the unstabilized standard error.
  stabilized <- wt_ate(
    ps,
    dat$z,
    exposure_type = "binary",
    .focal_level = 1,
    stabilize = TRUE
  )
  unstabilized <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)

  expect_equal(
    se_score_lin_rd(dat, ps_mod, stabilized),
    se_score_lin_rd(dat, ps_mod, unstabilized),
    tolerance = 1e-8
  )
})

test_that("linearization is unchanged by a scalar stabilization score", {
  # The companion to the default-stabilizer pin above, for a score supplied as a
  # single number rather than computed per group. One constant scales every
  # weight and the weight total identically, so it cancels from the Hajek ratio
  # and from the propensity score correction, and the whole estimates table has
  # to come back unchanged rather than only the risk difference standard error.
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- as.double(predict(ps_mod, type = "response"))

  unstabilized <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)
  scored <- function(score) {
    wt_ate(
      ps,
      dat$z,
      exposure_type = "binary",
      .focal_level = 1,
      stabilize = TRUE,
      stabilization_score = score
    )
  }
  lin_estimates <- function(wts) {
    ipw(
      ps_mod,
      se_score_outcome(dat, wts),
      .data = dat,
      se_method = "linearization"
    )$estimates
  }

  # A score of one leaves the weights themselves untouched, so that case is
  # exact; a score of 2.5 genuinely rescales them and has to cancel downstream.
  expect_identical(as.double(scored(1)), as.double(unstabilized))
  expect_false(identical(as.double(scored(2.5)), as.double(unstabilized)))
  expect_true(is_stabilized(scored(2.5)))

  reference <- lin_estimates(unstabilized)
  expect_identical(lin_estimates(scored(1)), reference)
  expect_equal(lin_estimates(scored(2.5)), reference, tolerance = 1e-12)
})

test_that("linearization is silent for weights carrying a vector stabilization score", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- predict(ps_mod, type = "response")

  # A per-observation score aligned with the weights: nothing about this call
  # changes the length of the weight vector, so nothing about the score is in
  # question and the call must say nothing about it.
  score <- ifelse(dat$z == 1, mean(dat$z), 1 - mean(dat$z))
  wts <- wt_ate(
    ps,
    dat$z,
    exposure_type = "binary",
    .focal_level = 1,
    stabilize = TRUE,
    stabilization_score = score
  )
  expect_length(stabilization_score(wts), nrow(dat))

  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts
  )

  res <- expect_no_warning(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
  )

  expect_s3_class(res, "ipw")
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
  expect_true(all(res$estimates$std.err > 0))
})

# ---- printCoefmat column formatting -----------------------------------------

test_that("print.ipw formats the z column as a test statistic, not a coefficient", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)
  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")

  # Reconstruct the numeric matrix print.ipw formats, then derive the two
  # candidate renderings. `fixed` marks estimate + std.err as the coefficient
  # pair and z as the test statistic; `current` is the off-by-one assignment
  # that marks std.err + z as the pair and ci.lower as the test statistic, under
  # which z prints at full precision rather than as a rounded test statistic.
  estimates <- res$estimates[-1]
  rownames(estimates) <- res$estimates$effect
  fixed <- capture.output(
    printCoefmat(estimates, has.Pvalue = TRUE, cs.ind = 1:2, tst.ind = 3)
  )
  current <- capture.output(
    printCoefmat(estimates, has.Pvalue = TRUE, cs.ind = 2:3, tst.ind = 4)
  )

  data_rows <- function(lines) grep("^(rd|log)", lines, value = TRUE)
  printed <- data_rows(capture.output(print(res)))

  # The two renderings genuinely differ (guards the discriminator), and the
  # printed table matches the correct assignment, not the off-by-one one.
  expect_false(identical(data_rows(fixed), data_rows(current)))
  expect_identical(printed, data_rows(fixed))
  expect_false(identical(printed, data_rows(current)))
})

# ---- linearization ps design extraction with .data ---------------------------
#
# The linearization path takes its propensity design from the shared
# ipw_extract_ps_design() helper, so it behaves like the M-estimation path when
# the data behind a propensity fit is gone: without .data the call raises the
# classed propensity_ipw_data_error that names .data as the remedy, and with
# .data the design is rebuilt and the estimates match a reconstructable
# reference fit. Redundant .data on a fit that needs none changes nothing. The
# parity pinned here is with the M-estimation tests in test-ipw-mestimation.R
# under "binary ps design extraction with .data".
#
# A .data whose row count disagrees with the fitted models is rejected. That
# guard exists because the design was once read from model.matrix(wt_mod) with
# .data ignored, so a shorter .data was recycled against the model-sized weights
# and the reported standard errors came out far too small with nothing signaled.

# A binary propensity score model fit with model = FALSE inside a scope whose
# fitting data is then gone. predict() still works from the stored fit, so the
# weights can be built, but model.matrix() cannot rebuild the design and .data is
# the only route to it. The inner name is distinct from any name in the calling
# scope, so the lookup fails on the whole parent chain.
se_method_ps_mod_gone <- function(dat) {
  local({
    d_local <- dat
    m <- glm(z ~ x1 + x2, data = d_local, family = binomial(), model = FALSE)
    rm(d_local)
    m
  })
}

test_that("linearization errors with the supply-.data hint when the ps fit frame is gone", {
  dat <- se_method_data()
  ps_gone <- se_method_ps_mod_gone(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_gone)

  # The lookup failure inside model.matrix() must surface as the guarded error
  # that names the missing data and directs the user to supply .data.
  expect_error(
    ipw(ps_gone, outcome_mod, se_method = "linearization"),
    class = "propensity_ipw_data_error"
  )
})

test_that("linearization reconstructs the ps design from .data when the fit frame is gone", {
  dat <- se_method_data()
  ps_gone <- se_method_ps_mod_gone(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_gone)

  # A reconstructable reference: the same formula on the same data, so the
  # coefficients are identical and only the retained model frame differs.
  ps_ref <- se_method_ps_mod(dat)
  dat_copy <- dat
  rm(dat)

  ref <- ipw(ps_ref, outcome_mod, se_method = "linearization")$estimates
  res <- ipw(
    ps_gone,
    outcome_mod,
    .data = dat_copy,
    se_method = "linearization"
  )$estimates
  expect_equal(res, ref, tolerance = 1e-8)
})

test_that("linearization estimates are unchanged when redundant .data is supplied", {
  # A normal, reconstructable ps fit: rebuilding the design from .data must
  # produce the same design and so the same estimates. Guards the rebuild against
  # drifting when .data is redundant.
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  no_data <- ipw(ps_mod, outcome_mod, se_method = "linearization")$estimates
  with_data <- ipw(
    ps_mod,
    outcome_mod,
    .data = dat,
    se_method = "linearization"
  )$estimates
  expect_equal(with_data, no_data, tolerance = 1e-8)
})

test_that("linearization errors when .data has fewer rows than the fitted models", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  # 100 rows of .data against 400-row fits. The weights still come from the
  # outcome model, so an unreconciled .data would be recycled against them and
  # report standard errors far too small; the mismatch must be rejected instead.
  expect_error(
    ipw(
      ps_mod,
      outcome_mod,
      .data = head(dat, 100),
      se_method = "linearization"
    ),
    class = "propensity_ipw_data_error",
    regexp = "length|rows"
  )
})

# A `.data` column whose type contradicts the propensity fit's recorded factor
# coding is rejected here for the same reason it is on the M-estimation path,
# and by the same shared extraction: the design is rebuilt under the coding the
# fit recorded, so a numeric column reaches `model.matrix()` with a contrast
# specification that cannot be applied to it. The M-estimation companions are
# "mestimation rejects a numeric .data column where the ps model fit a factor"
# and its character control in test-ipw-mestimation.R.
#
# The outcome model's factor covariates cannot be reached from here: an outcome
# model adjusted for anything beyond the exposure is refused on this path
# before any design is built, so that direction is pinned on the M-estimation
# path alone.

# A binary ps fit adjusted for a three-level factor covariate, with the
# exposure-only outcome model this path requires. Local to this file, as
# test-ipw-mestimation.R's equivalent fixture is not reachable here.
se_method_ps_factor_models <- function(dat) {
  dat$grp <- factor(
    ifelse(dat$x1 < -0.5, "lo", ifelse(dat$x1 < 0.5, "mid", "hi")),
    levels = c("lo", "mid", "hi")
  )
  ps_mod <- glm(z ~ x1 + grp, data = dat, family = binomial())
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_mod))
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, dat = dat)
}

test_that("linearization rejects a numeric .data column where the ps model fit a factor", {
  mods <- se_method_ps_factor_models(se_method_data())
  expect_identical(mods$ps_mod$xlevels$grp, c("lo", "mid", "hi"))

  dat_num <- mods$dat
  dat_num$grp <- as.numeric(mods$dat$grp)

  err <- expect_no_warning(
    expect_error(
      ipw(
        mods$ps_mod,
        mods$outcome_mod,
        .data = dat_num,
        se_method = "linearization"
      ),
      class = "propensity_ipw_data_error"
    )
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "\\bgrp\\b")
  expect_match(msg, "factor", fixed = TRUE)
  expect_false(grepl("contrasts apply only to factors", msg, fixed = TRUE))
})

test_that("linearization accepts a character .data column where the ps model fit a factor", {
  # The guard must not widen to every column whose type differs from the fit's.
  # `model.frame()` re-levels a character column against the recorded levels, so
  # the rebuilt design is the fitted one and the estimates match the factor
  # `.data` call, which is the working control for this seam.
  mods <- se_method_ps_factor_models(se_method_data())

  dat_chr <- mods$dat
  dat_chr$grp <- as.character(mods$dat$grp)

  chr <- expect_no_warning(
    ipw(
      mods$ps_mod,
      mods$outcome_mod,
      .data = dat_chr,
      se_method = "linearization"
    )
  )$estimates
  fac <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = mods$dat,
    se_method = "linearization"
  )$estimates
  expect_equal(chr, fac, tolerance = 1e-12)
})

# ---- the weights must be consistent with the propensity model ---------------
#
# The linearization path predicts the propensity scores from `wt_mod` and takes
# the weights from the outcome model fit, and nothing checks that the two belong
# together. The M-estimation path checks, by recomputing the weights the
# propensity scores imply and comparing them to the ones actually used, and the
# same check belongs here.
#
# What goes wrong when they disagree is confined to the standard errors. The
# point estimates come from the outcome model and its weights, so they are
# unmoved; the propensity scores enter only through the estimation correction.
# A wrong standard error that agrees to the last digit on the estimate beside it
# is exactly the kind of error nobody looks twice at, which is why these are
# rejections rather than warnings.
#
# The correct-`.data` and intact-weight cases are pinned already and are not
# repeated here: "linearization estimates are unchanged when redundant .data is
# supplied", "linearization default stabilization reproduces the unstabilized
# SE", "linearization is silent for weights carrying a vector stabilization
# score", and "linearization threads a covariate-correlated stabilization score
# through the weight derivatives".

test_that("linearization rejects a .data whose covariates disagree with the fitted weights", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  # Row-aligned and the right size, so nothing about the shape is wrong; only
  # the covariate values differ, which moves the predicted propensity scores
  # away from the ones the weights were built on. The M-estimation path rejects
  # this same call.
  perturbed <- dat
  perturbed$x1 <- perturbed$x1 + 1

  err <- expect_error(
    ipw(ps_mod, outcome_mod, .data = perturbed, se_method = "linearization"),
    class = "propensity_ipw_weights_mismatch_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(
    msg,
    "not consistent with the propensity score model",
    fixed = TRUE
  )
})

test_that("linearization rejects stabilized weights whose score was dropped", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- as.double(predict(ps_mod, type = "response"))
  score <- exp(0.5 * dat$x1)
  wts <- wt_ate(
    ps,
    dat$z,
    exposure_type = "binary",
    .focal_level = 1,
    stabilize = TRUE,
    stabilization_score = score
  )

  # Slicing a psw drops a per-observation score it can no longer align, and says
  # so, but leaves the weights marked stabilized. That is the state below: the
  # values are still the score-stabilized ones, and nothing records which score
  # produced them.
  sliced <- NULL
  expect_warning(
    sliced <- wts[1:200],
    class = "propensity_stabilization_score_warning"
  )
  expect_true(is_stabilized(sliced))
  expect_null(stabilization_score(sliced))

  # Reproduce that state at full length so the weights remain the ones this
  # propensity model implies and the dropped score is the only thing wrong.
  # Slicing itself cannot be used here: it would also shorten the weights, and
  # the row-count guard would reject the call before the weights were examined.
  degraded <- wts
  attr(degraded, "stabilization_score") <- NULL
  expect_true(is_stabilized(degraded))
  expect_null(stabilization_score(degraded))
  expect_identical(as.double(degraded), as.double(wts))

  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = degraded,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  # Without the score there is nothing to rebuild these weights from, and the
  # default group-constant stabilizer is not what produced them: it differs by
  # more than 8 on some observations. Reconstructing it silently is the defect.
  default_stabilized <- wt_ate(
    ps,
    dat$z,
    exposure_type = "binary",
    .focal_level = 1,
    stabilize = TRUE
  )
  expect_gt(max(abs(as.double(default_stabilized) - as.double(wts))), 1)

  err <- expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_ipw_weights_mismatch_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(
    msg,
    "not consistent with the propensity score model",
    fixed = TRUE
  )
})

# ---- ps_link is deprecated, and must agree with the model's own link --------
#
# On the linearization path `ps_link` reaches `derive_weights`, the score factor,
# and the correction matrix, none of which consult the fitted model. Naming a
# supported link that is not the one `wt_mod` was fit with therefore scales the
# estimation correction by the wrong derivative, and the only thing that moves is
# the standard error: the point estimates come from the outcome model and its
# weights and are unchanged to the last digit. The weight-consistency preflight
# cannot see it either, since the weights do not depend on the link.
#
# The M-estimation path refuses the same misuse for a different reason: there the
# link inverts the coefficient block into propensity scores, so the recomputed
# weights no longer reproduce the ones the outcome model was fit with and the
# misuse surfaces as a weights mismatch. That boundary is pinned below as well,
# so neither path drifts into the other's message.
#
# Since both paths accept only the model's own link, the one value `ps_link` can
# carry without changing anything is the value the default already resolves to,
# which leaves the argument as pure redundancy. Supplying it therefore raises a
# lifecycle deprecation warning, the matching value included; omitting it, or
# passing NULL, resolves the model's link silently as before. The warning is
# advisory, so the pins below assert every number against the call that omits the
# argument.
#
# The pins force `lifecycle_verbosity = "warning"` rather than relying on the
# default once-per-session behavior, so no pin depends on being the first to
# reach the `ipw(ps_link)` id.
#
# The legitimate probit path is covered by "probit linearization RD SE matches
# the generalized score correction" and "probit linearization RD SE agrees with
# mestimation" above; nothing here should disturb it.

se_link_fit <- function(dat, link) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial(link = link))
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      predict(ps_mod, type = "response"),
      dat$z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod)
}

# Collect the deprecation warnings a call signals, without letting them escape
# into the test's own condition record, and hand back the call's value alongside
# them so a pin can assert on both.
se_link_deprecations <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    lifecycle_warning_deprecated = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      rlang::cnd_muffle(cnd)
    }
  )
  list(value = value, messages = messages)
}

test_that("linearization rejects a ps_link that differs from the model's link", {
  dat <- se_method_data()
  withr::local_options(lifecycle_verbosity = "warning")

  # Both directions, so the guard cannot be written as a rule about logit.
  cases <- list(
    list(fit = "logit", ps_link = "probit"),
    list(fit = "logit", ps_link = "cloglog"),
    list(fit = "probit", ps_link = "logit")
  )

  for (case in cases) {
    mods <- se_link_fit(dat, case$fit)
    caught <- se_link_deprecations(
      expect_error(
        ipw(
          mods$ps_mod,
          mods$outcome_mod,
          ps_link = case$ps_link,
          se_method = "linearization"
        ),
        class = "propensity_ipw_link_error"
      )
    )

    # deprecated first, then rejected: the warning does not stand in for the
    # guard, and the guard does not swallow the warning
    expect_length(caught$messages, 1)

    # the message has to name both, or it cannot say what disagrees with what
    msg <- gsub("[[:space:]]+", " ", conditionMessage(caught$value))
    expect_match(msg, case$fit, fixed = TRUE)
    expect_match(msg, case$ps_link, fixed = TRUE)
  }
})

test_that("the ps_link mismatch error names both links", {
  dat <- se_method_data()
  mods <- se_link_fit(dat, "logit")

  # The deprecation warning that now precedes the rejection is pinned by
  # "linearization rejects a ps_link that differs from the model's link".
  # Silencing it here keeps this snapshot a record of the rejection message
  # alone, rather than one that also moves whenever the package version does.
  withr::local_options(lifecycle_verbosity = "quiet")

  expect_snapshot(
    error = TRUE,
    ipw(
      mods$ps_mod,
      mods$outcome_mod,
      ps_link = "probit",
      se_method = "linearization"
    )
  )
})

test_that("mestimation reports a mismatched ps_link as a weights mismatch", {
  skip_if_not_installed("deli")
  # The deliberate boundary between the two paths, pinned so neither drifts into
  # the other. The membership guard checks only that the link is one this package
  # supports, and a supported-but-wrong link passes it. What catches the misuse
  # on this path is the weight recomputation: the link inverts the coefficient
  # block into propensity scores, so a mismatched link produces weights that do
  # not reproduce the ones the outcome model was fit with. The linearization path
  # cannot detect it that way, because its weights do not depend on the link, and
  # rejects the same misuse with a link error instead; see "linearization rejects
  # a ps_link that differs from the model's link".
  dat <- se_method_data()
  mods <- se_link_fit(dat, "logit")

  expect_identical(mods$ps_mod$family$link, "logit")

  withr::local_options(lifecycle_verbosity = "warning")
  caught <- se_link_deprecations(
    expect_error(
      ipw(
        mods$ps_mod,
        mods$outcome_mod,
        ps_link = "probit",
        se_method = "mestimation"
      ),
      class = "propensity_ipw_weights_mismatch_error"
    )
  )

  # the deprecation reaches this path too, and the mismatch still follows it
  expect_length(caught$messages, 1)
})

test_that("linearization deprecates a ps_link equal to the model's link", {
  dat <- se_method_data()
  withr::local_options(lifecycle_verbosity = "warning")

  # Naming the link the model was already fit with asks for exactly what the
  # default gives, so it is deprecated rather than rejected. The warning is the
  # whole of the change: the estimates must still match the call that omits the
  # argument, value for value.
  for (link in c("logit", "probit")) {
    mods <- se_link_fit(dat, link)

    caught <- se_link_deprecations(
      ipw(
        mods$ps_mod,
        mods$outcome_mod,
        ps_link = link,
        se_method = "linearization"
      )
    )

    expect_length(caught$messages, 1)

    # the warning has to name the argument and the function it belongs to
    msg <- gsub("[[:space:]]+", " ", paste(caught$messages, collapse = " "))
    expect_match(msg, "ps_link", fixed = TRUE)
    expect_match(msg, "ipw()", fixed = TRUE)
    expect_match(msg, "deprecated", fixed = TRUE)

    # the baseline runs at the same verbosity, so a deprecation fired for the
    # default resolution would escape here and be recorded
    expect_identical(
      caught$value$estimates,
      ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")$estimates
    )
  }
})

test_that("mestimation deprecates a ps_link equal to the model's link", {
  skip_if_not_installed("deli")
  dat <- se_method_data()
  withr::local_options(lifecycle_verbosity = "warning")

  mods <- se_link_fit(dat, "logit")
  caught <- se_link_deprecations(
    ipw(
      mods$ps_mod,
      mods$outcome_mod,
      ps_link = "logit",
      se_method = "mestimation"
    )
  )

  expect_length(caught$messages, 1)
  expect_identical(
    caught$value$estimates,
    ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")$estimates
  )
})

test_that("omitting ps_link deprecates nothing on either path", {
  skip_if_not_installed("deli")
  dat <- se_method_data()
  withr::local_options(lifecycle_verbosity = "warning")

  mods <- se_link_fit(dat, "logit")

  # Resolving the model's own link is the default rather than a use of the
  # argument, so neither the omitted form nor an explicit NULL may warn, even at
  # the loudest verbosity.
  for (method in c("linearization", "mestimation")) {
    expect_no_warning(
      ipw(mods$ps_mod, mods$outcome_mod, se_method = method),
      class = "lifecycle_warning_deprecated"
    )
    expect_no_warning(
      ipw(mods$ps_mod, mods$outcome_mod, ps_link = NULL, se_method = method),
      class = "lifecycle_warning_deprecated"
    )
  }
})

# ---- focal level flipped against the coded-1 exposure level -----------------
#
# The binary path has no focal level of its own: it codes the second sorted
# exposure level as 1, and the estimand derivatives assume that coding. att
# weights targeting the first sorted level are numerically identical to atu
# weights targeting the second, so weights labelled att but built the other way
# round reach the linearization path and are corrected with the treated and
# control roles mirrored, producing standard errors for an estimand nobody asked
# for. The mestimation path rejects the same weights, but only through the
# weight-consistency preflight, whose hints never name the cause. Both paths must
# reject a recorded focal level that is not the level they code as 1.

se_method_data_factor <- function() {
  dat <- se_method_data()
  dat$zf <- factor(
    ifelse(dat$z == 1, "treat", "control"),
    levels = c("control", "treat")
  )
  dat
}

se_method_ps_mod_factor <- function(dat) {
  glm(zf ~ x1 + x2, data = dat, family = binomial())
}

# Weights targeting "control", the first sorted exposure level and the level the
# binary path treats as unexposed. The propensity score is flipped to 1 - e so
# the values are the correct weights for that focal level: the att weights below
# equal wt_atu(e, zf) exactly, and the atu weights equal wt_att(e, zf).
se_method_outcome_focal_control <- function(dat, ps_mod, wt_fun) {
  e_treat <- predict(ps_mod, type = "response")
  wts <- wt_fun(
    1 - e_treat,
    dat$zf,
    exposure_type = "binary",
    .focal_level = "control"
  )
  glm(y ~ zf, data = dat, family = quasibinomial(), weights = wts)
}

test_that("linearization rejects att weights focal on the first exposure level", {
  dat <- se_method_data_factor()
  ps_mod <- se_method_ps_mod_factor(dat)
  outcome_mod <- se_method_outcome_focal_control(dat, ps_mod, wt_att)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_focal_level_error",
    regexp = "control"
  )
})

test_that("mestimation rejects att weights focal on the first exposure level", {
  dat <- se_method_data_factor()
  ps_mod <- se_method_ps_mod_factor(dat)
  outcome_mod <- se_method_outcome_focal_control(dat, ps_mod, wt_att)

  # The weight-consistency preflight also rejects these weights, but its hints
  # describe a refit rather than the flipped focal level, so the focal guard must
  # run first.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "mestimation"),
    class = "propensity_focal_level_error",
    regexp = "control"
  )
})

test_that("mestimation rejects atu weights focal on the first exposure level", {
  # The atu mirror of the att case. The linearization path has no atu support to
  # reach, so this pair is pinned on the mestimation path alone.
  dat <- se_method_data_factor()
  ps_mod <- se_method_ps_mod_factor(dat)
  outcome_mod <- se_method_outcome_focal_control(dat, ps_mod, wt_atu)

  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "mestimation"),
    class = "propensity_focal_level_error",
    regexp = "control"
  )
})

test_that("the focal level rejection names both the recorded and coded levels", {
  dat <- se_method_data_factor()
  ps_mod <- se_method_ps_mod_factor(dat)
  outcome_mod <- se_method_outcome_focal_control(dat, ps_mod, wt_att)

  expect_snapshot(
    error = TRUE,
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
  )
})

test_that("att weights focal on the second exposure level estimate unchanged", {
  dat <- se_method_data_factor()
  ps_mod <- se_method_ps_mod_factor(dat)
  e_treat <- predict(ps_mod, type = "response")

  # "treat" is the level the binary path codes as 1, so recording it must leave
  # the estimates exactly where the same weights without a recorded focal level
  # put them.
  wts_focal <- wt_att(
    e_treat,
    dat$zf,
    exposure_type = "binary",
    .focal_level = "treat"
  )
  wts_plain <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_att(e_treat, dat$zf, exposure_type = "binary")
  )
  expect_equal(as.double(wts_focal), as.double(wts_plain))

  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)
  outcome_focal <- glm(
    y ~ zf,
    data = dat,
    family = quasibinomial(),
    weights = wts_focal,
    control = ctrl
  )
  outcome_plain <- glm(
    y ~ zf,
    data = dat,
    family = quasibinomial(),
    weights = wts_plain,
    control = ctrl
  )

  for (method in c("linearization", "mestimation")) {
    res_focal <- ipw(ps_mod, outcome_focal, .data = dat, se_method = method)
    res_plain <- ipw(ps_mod, outcome_plain, .data = dat, se_method = method)
    expect_equal(res_focal$estimates, res_plain$estimates, tolerance = 1e-8)
  }
})

# ---- separation in the propensity model on the linearization path -----------
#
# A propensity model whose covariates predict the exposure without error has no
# finite maximum likelihood estimate, and the weights it implies are undefined
# wherever a fitted probability reaches the boundary. The M-estimation path
# refuses such a fit outright.
#
# The linearization path has nothing that breaks. Its scores come from
# `predict(type = "response")`, which goes through the fitted family's inverse
# link, and that link clamps: it cannot return an exact 0 or 1. Every weight is
# then finite, no influence value divides by zero, and a design with no overlap
# at all yields an estimate with a small standard error beside it.
#
# Both paths refuse the same fit at the same threshold: the fitted linear
# predictors, put through the link's unclamped inverse, reach exactly 0 or
# exactly 1 for at least one observation. Anything short of saturation still
# runs. Graded overlap diagnostics are not this package's business, and the
# near-separated fixture below is here to keep them out of it.

se_separation_data <- function(seed = 2024, n = 400) {
  set.seed(seed)
  x1 <- rnorm(n)
  # x1 predicts z with no error, so the fit has no finite optimum
  z <- as.integer(x1 > 0)
  y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.4 * x1))
  data.frame(x1, z, y)
}

# A fit that is near separation but whose maximum likelihood estimate is finite.
# The linear predictors reach about 25 in absolute value, so the fitted
# probabilities come within 1e-11 of the boundary without touching it.
se_near_separation_data <- function(seed = 11, n = 400) {
  set.seed(seed)
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(8 * x1))
  y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.4 * x1))
  data.frame(x1, z, y)
}

# The propensity fit warns that fitted probabilities of zero or one occurred,
# and on the separated fixture it also exhausts its iterations; both are the
# fixture working as intended rather than anything under test. The iteration
# limit is raised so the linear predictor reaches the range where the unclamped
# inverse link saturates. The weights are built from the fitted model the way a
# user would, so they are exactly the ones this propensity model implies and
# the weight-consistency preflight has nothing to object to.
se_separation_models <- function(dat, link = "logit") {
  ps_mod <- suppressWarnings(glm(
    z ~ x1,
    data = dat,
    family = binomial(link),
    control = glm.control(maxit = 200, epsilon = 1e-14)
  ))
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_mod)
  )
  outcome_mod <- suppressWarnings(glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  ))
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# Count the fitted scores that saturate under the link's unclamped inverse,
# which is what the threshold is stated in terms of. Neither
# `binomial()$linkinv` nor `stats::make.link()$linkinv` can measure this: they
# are the same clamped function, and it never reaches the boundary.
se_n_saturated <- function(ps_mod) {
  inv <- switch(
    ps_mod$family$link,
    logit = stats::plogis,
    probit = stats::pnorm
  )
  e <- inv(predict(ps_mod))
  sum(e == 0 | e == 1)
}

test_that("linearization rejects a separated propensity model", {
  dat <- se_separation_data()
  mods <- se_separation_models(dat)

  # The premise. The fit is separated, `predict()` still hands back interior
  # scores because the family's inverse link clamps, and the weights the user
  # holds are finite and agree with the model. Nothing here breaks on its own.
  expect_gt(max(abs(predict(mods$ps_mod))), 100)
  ps <- as.double(predict(mods$ps_mod, type = "response"))
  expect_false(any(ps == 0 | ps == 1))
  expect_true(all(is.finite(as.double(mods$wts))))

  n_saturated <- se_n_saturated(mods$ps_mod)
  expect_gt(n_saturated, 0)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization"),
    class = "propensity_ipw_separation_error"
  )

  # The guard's own message, naming the count of saturated scores, rather than
  # a downstream failure that happens to arrive first. The weights are the ones
  # this model implies, so the weight-consistency preflight is not the thing
  # speaking here and must not be.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "separation", ignore.case = TRUE)
  expect_match(msg, as.character(n_saturated), fixed = TRUE)
  expect_false(inherits(err, "propensity_ipw_weights_mismatch_error"))
})

test_that("linearization rejects a separated propensity model given .data", {
  dat <- se_separation_data()
  mods <- se_separation_models(dat)
  n_saturated <- se_n_saturated(mods$ps_mod)

  # The scores are predicted from `.data` when it is supplied and from the
  # model frame when it is not, so the guard has to cover both routes.
  err <- expect_error(
    ipw(
      mods$ps_mod,
      mods$outcome_mod,
      .data = dat,
      se_method = "linearization"
    ),
    class = "propensity_ipw_separation_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, as.character(n_saturated), fixed = TRUE)

  # And that the supplied `.data` is what the guard reads, rather than the model
  # frame standing in for it: this fit converges and saturates nowhere on its own
  # frame, but the scores it gives on a `.data` whose covariate is doubled do
  # saturate, and the guard reaches that before the weight comparison does.
  near <- se_near_separation_data()
  near_mods <- se_separation_models(near)
  expect_identical(se_n_saturated(near_mods$ps_mod), 0L)

  doubled <- near
  doubled$x1 <- 2 * doubled$x1
  e_doubled <- stats::plogis(predict(near_mods$ps_mod, newdata = doubled))
  expect_gt(sum(e_doubled == 0 | e_doubled == 1), 0)

  expect_error(
    ipw(
      near_mods$ps_mod,
      near_mods$outcome_mod,
      .data = doubled,
      se_method = "linearization"
    ),
    class = "propensity_ipw_separation_error"
  )
})

test_that("both SE methods reject a separated propensity model alike", {
  skip_if_not_installed("deli")
  dat <- se_separation_data()
  mods <- se_separation_models(dat)

  errs <- lapply(c("linearization", "mestimation"), function(method) {
    expect_error(
      ipw(mods$ps_mod, mods$outcome_mod, se_method = method),
      class = "propensity_ipw_separation_error"
    )
  })

  # Same case and same threshold, so both paths count the same observations.
  n_saturated <- se_n_saturated(mods$ps_mod)
  for (err in errs) {
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "separation", ignore.case = TRUE)
    expect_match(msg, as.character(n_saturated), fixed = TRUE)
  }
})

test_that("linearization still runs on a near-separated propensity model", {
  dat <- se_near_separation_data()
  mods <- se_separation_models(dat)

  # The fixture's premise, asserted here so it cannot drift into saturation and
  # turn this test vacuous: the estimate is finite, the scores sit within 1e-9
  # of both boundaries, and not one of them reaches either.
  expect_true(mods$ps_mod$converged)
  e <- stats::plogis(predict(mods$ps_mod))
  expect_identical(se_n_saturated(mods$ps_mod), 0L)
  expect_lt(min(e), 1e-9)
  expect_lt(1 - max(e), 1e-9)

  # No overlap diagnostic fires short of saturation, so this returns.
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  expect_s3_class(res, "ipw")
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("linearization rejects a separated probit propensity model", {
  dat <- se_separation_data()
  mods <- se_separation_models(dat, link = "probit")

  # `pnorm` saturates around 8.3, far below the 36.7 `plogis` needs, so the same
  # fixture separates a probit fit at least as hard. The clamped inverse link
  # still returns interior scores, so the path has as little to break on here.
  ps <- as.double(predict(mods$ps_mod, type = "response"))
  expect_false(any(ps == 0 | ps == 1))
  expect_gt(se_n_saturated(mods$ps_mod), 0)
  expect_true(all(is.finite(as.double(mods$wts))))

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization"),
    class = "propensity_ipw_separation_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "separation", ignore.case = TRUE)
})

test_that("linearization leaves a separated cauchit fit to the link rejection", {
  dat <- se_separation_data()
  mods <- se_separation_models(dat, link = "cauchit")

  # Separated exactly as the logit and probit fixtures are, so the guard would
  # have something to say if it applied here.
  expect_gt(max(abs(predict(mods$ps_mod))), 100)

  # It does not apply. The guard measures saturation through the link's
  # unclamped inverse, and there is none for cauchit, so it steps aside rather
  # than aborting. What speaks instead is the rejection this path already made
  # of a link its influence functions are not derived for, unchanged and
  # carrying none of the package's own error classes.
  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization"),
    class = "rlang_error"
  )
  expect_match(conditionMessage(err), "cauchit", fixed = TRUE)
  expect_false(inherits(err, "propensity_ipw_separation_error"))
  expect_false(inherits(err, "propensity_error"))
})

test_that("the linearization separation error names the count and the model", {
  dat <- se_separation_data()
  mods <- se_separation_models(dat)

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")
  )
})
