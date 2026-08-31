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
# linearization path are detected. These are recorded values, so the tolerance
# they are read at is set by what varies between the platform that recorded them
# and the platform reading them: the trailing digits move with the BLAS the
# linear algebra runs through. A change to the linearization path itself moves
# these numbers by orders of magnitude more, so 1e-6 still catches what the pins
# exist to catch.
se_method_lin_reference <- function() {
  data.frame(
    effect = c("mean", "mean", "rd", "log(rr)", "log(or)"),
    contrast = c("0", "1", rep("1 vs 0", 3)),
    estimate = c(
      0.29744534143394004,
      0.54053118729710203,
      0.24308584586316151,
      0.59732185440835850,
      1.02197399113135590
    ),
    std.err = c(
      0.032268480487667904,
      0.036574156651138997,
      0.048386022223687541,
      0.126936312969173404,
      0.211670119599652012
    )
  )
}

estimates_columns <- c(
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
  expect_equal(res$estimates$contrast, reference$contrast)
  expect_equal(res$estimates$estimate, reference$estimate, tolerance = 1e-6)
  expect_equal(res$estimates$std.err, reference$std.err, tolerance = 1e-6)
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
  expect_equal(df_m$term, c("mean", "mean", "rd", "rr", "or"))
  expect_equal(df_l$term, c("mean", "mean", "rd", "rr", "or"))
  expect_equal(df_m$estimate, df_l$estimate, tolerance = 1e-8)
})

test_that("ipw() print output is stable per SE method", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  # The row labels name the effect measure and the contrast together, so the
  # table needs more than the 80 columns testthat pins or it wraps.
  withr::local_options(width = 120)

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
    tolerance = 1e-6
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
  # branch extracts the outcome two ways that would fail identically: the
  # response evaluated against `.data` when it is supplied, and
  # fmla_extract_left_vctr when it is not. Both return finite std.err equal to
  # the numeric fit with no warnings.
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

# The numeric content of an estimates frame, with the labels naming its rows set
# aside. A factor exposure and a 0/1 recode of it are the same fit reporting the
# same numbers, and each writes its labels with its own level names, so the
# comparison between the two arms is of the numbers alone. The labels are pinned
# separately, since which level names a fit writes is itself part of the
# contract.
factor_estimates_values <- function(estimates) {
  out <- estimates[setdiff(names(estimates), "contrast")]
  attr(out, "ipw_vcov") <- unname(attr(estimates, "ipw_vcov", exact = TRUE))
  out
}

# The rows a binary fit over `levels` reports, in order.
factor_row_labels <- function(levels, contrasts) {
  c(levels, rep(paste(levels[[2]], "vs", levels[[1]]), length(contrasts)))
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
  expect_equal(
    factor_estimates_values(res_nodata),
    factor_estimates_values(ref),
    tolerance = 1e-8
  )
  expect_equal(
    factor_estimates_values(res_data),
    factor_estimates_values(ref),
    tolerance = 1e-8
  )

  # The labels are written with the levels the fit codes the exposure on, so the
  # factor arm names its arms and the 0/1 arm names its codes.
  expect_identical(
    res_nodata$contrast,
    factor_row_labels(c("control", "treated"), c("rd", "log(rr)", "log(or)"))
  )
  expect_identical(res_nodata$contrast, res_data$contrast)
  expect_identical(
    ref$contrast,
    factor_row_labels(c("0", "1"), c("rd", "log(rr)", "log(or)"))
  )
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
  expect_equal(
    factor_estimates_values(res),
    factor_estimates_values(ref),
    tolerance = 1e-8
  )

  # The reference level is the fit's own first level, which here is "treated".
  expect_identical(
    res$contrast,
    factor_row_labels(c("treated", "control"), c("rd", "log(rr)", "log(or)"))
  )
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
  expect_equal(
    factor_estimates_values(res),
    factor_estimates_values(ref),
    tolerance = 1e-8
  )
})

test_that("mestimation matches a factor exposure to the numeric fit", {
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
    expect_equal(
      factor_estimates_values(res_fac),
      factor_estimates_values(res_num),
      tolerance = 1e-8
    )
    expect_identical(
      res_fac$contrast,
      factor_row_labels(
        levels(arms$dat[[factor_col]]),
        c("rd", "log(rr)", "log(or)")
      )
    )
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

test_that("a matrix response is named as one when the model frame is gone", {
  # The shape is read from the model frame where there is one, and a fit made
  # with model = FALSE whose data are gone has none. The formula still says
  # cbind, so the report stays the matrix one rather than falling to the wording
  # for a response that is one column written as an expression.
  m <- matrix_lhs_models()
  ps_gone <- local({
    d_local <- m$dat
    fit <- glm(
      cbind(y1, y0) ~ x1,
      data = d_local,
      family = binomial(),
      model = FALSE
    )
    rm(d_local)
    fit
  })
  # the fixture is the shape it claims to be: no frame to read
  expect_error(stats::model.frame(ps_gone))

  err <- expect_error(
    ipw(ps_gone, m$outcome_mod, .data = m$dat),
    class = "propensity_ipw_response_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "matrix response", fixed = TRUE)
  expect_false(grepl("an expression rather than", msg, fixed = TRUE))
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

# A lognormal outcome, the propensity model whose weights it carries, and a
# precomputed column for each transformation, so every transformed fit has a
# reference written on a bare symbol. Shared with the `.data` route pins below.
transformed_response_models <- function() {
  set.seed(2024)
  n <- 300
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x1))
  y <- rlnorm(n, -0.4 + z + 0.3 * x1)
  dat <- data.frame(
    x1,
    z,
    y,
    log_y = log(y),
    sqrt_y = sqrt(y),
    scale_y = as.vector(scale(y))
  )
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
  list(ps_mod = ps_mod, dat = dat, wts = wts)
}

transformed_response_outcome <- function(response, dat, wts) {
  lm(stats::reformulate("z", response = response), data = dat, weights = wts)
}

test_that("a transformed single-column outcome response is not a matrix response", {
  # `log(y)` and friends make the formula's left-hand side a call rather than a
  # symbol, which is what the propensity guard rejects and what an early version
  # of this guard rejected too. They are ordinary single-column responses and
  # the analysis is correct, so the only honest check is the built frame. The
  # reference fits the same model on a precomputed column, which is the same
  # model written differently and must give the same answer.
  m <- transformed_response_models()

  for (se in c("mestimation", "linearization")) {
    for (pair in list(c("log(y)", "log_y"), c("sqrt(y)", "sqrt_y"))) {
      transformed <- transformed_response_outcome(pair[[1]], m$dat, m$wts)
      precomputed <- transformed_response_outcome(pair[[2]], m$dat, m$wts)

      expect_false(
        is.matrix(stats::model.response(stats::model.frame(transformed)))
      )
      expect_equal(
        ipw(m$ps_mod, transformed, se_method = se)$estimates$estimate,
        ipw(m$ps_mod, precomputed, se_method = se)$estimates$estimate,
        tolerance = 1e-10
      )
    }
  }
})

test_that("a scale-transformed outcome response matches its precomputed column", {
  # The log and sqrt pair above leaves a plain vector in the model frame.
  # `scale()` leaves a one-column matrix there instead, so it is the shape the
  # precomputed cross-check has nothing to say about yet, and it is excluded
  # from that loop because it fails the matrix premise asserted there. Same
  # model written two ways, so the same estimates.
  m <- transformed_response_models()
  transformed <- transformed_response_outcome("scale(y)", m$dat, m$wts)
  precomputed <- transformed_response_outcome("scale_y", m$dat, m$wts)

  expect_true(is.matrix(stats::model.frame(transformed)[[1]]))
  expect_false(is.matrix(stats::model.frame(precomputed)[[1]]))

  for (se in c("mestimation", "linearization")) {
    expect_equal(
      ipw(m$ps_mod, transformed, se_method = se)$estimates$estimate,
      ipw(m$ps_mod, precomputed, se_method = se)$estimates$estimate,
      tolerance = 1e-10
    )
  }
})

test_that("the matrix-response outcome model error names the outcome model", {
  m <- matrix_outcome_models()

  expect_snapshot(
    error = TRUE,
    ipw(m$ps_mod, m$outcome_mod)
  )
})

test_that("a matrix-response outcome model errors on its shape with .data", {
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

# ---- transformed outcome responses on the .data route -----------------------
#
# A transformed left-hand side such as `log(y)` is an ordinary single-column
# response, and the tests above pin that it estimates correctly when the models
# carry their own frames. The `.data` route reads the response by name instead,
# so the name it asks for has to be the column the transformation reads rather
# than the function wrapping it. These pin the two routes against each other,
# and pin that a missing-column report only ever names a column the user could
# supply. The second contract has no test of its own: with a complete `.data`
# the transformed calls simply run, which is what the parity pin asserts.
#
# A binomial outcome with a transformed response is deliberately absent. The
# no-`.data` pins above cover only gaussian fits, so there is no pinned
# reference to compare a binomial one against.

test_that("a log-transformed outcome response through .data matches the model-frame route", {
  m <- transformed_response_models()
  outcome_mod <- transformed_response_outcome("log(y)", m$dat, m$wts)

  for (se in c("mestimation", "linearization")) {
    expect_equal(
      ipw(m$ps_mod, outcome_mod, .data = m$dat, se_method = se)$estimates,
      ipw(m$ps_mod, outcome_mod, se_method = se)$estimates,
      tolerance = 1e-10
    )
  }
})

test_that("a sqrt-transformed outcome response through .data matches the model-frame route", {
  m <- transformed_response_models()
  outcome_mod <- transformed_response_outcome("sqrt(y)", m$dat, m$wts)

  expect_equal(
    ipw(
      m$ps_mod,
      outcome_mod,
      .data = m$dat,
      se_method = "linearization"
    )$estimates,
    ipw(m$ps_mod, outcome_mod, se_method = "linearization")$estimates,
    tolerance = 1e-10
  )
})

test_that("a scale-transformed outcome response through .data matches the model-frame route", {
  m <- transformed_response_models()
  outcome_mod <- transformed_response_outcome("scale(y)", m$dat, m$wts)

  # `scale()` is the mestimation case on purpose: it leaves a one-column matrix
  # in the model frame, so the psi rebuild sees a shape the other transformations
  # never produce. Pin that premise, since the case is pointless without it.
  expect_true(is.matrix(stats::model.frame(outcome_mod)[[1]]))

  expect_equal(
    ipw(
      m$ps_mod,
      outcome_mod,
      .data = m$dat,
      se_method = "mestimation"
    )$estimates,
    ipw(m$ps_mod, outcome_mod, se_method = "mestimation")$estimates,
    tolerance = 1e-10
  )
})

test_that("a covariate missing from .data is reported as the missing column", {
  # Guard for the pins above: the report has to keep working for a column that
  # really is absent, not go quiet along with the spurious ones.
  m <- transformed_response_models()
  outcome_mod <- transformed_response_outcome("log_y", m$dat, m$wts)

  err <- expect_error(
    ipw(
      m$ps_mod,
      outcome_mod,
      .data = m$dat[setdiff(names(m$dat), "x1")],
      se_method = "mestimation"
    ),
    class = "propensity_columns_exist_error"
  )
  expect_match(conditionMessage(err), "x1", fixed = TRUE)
})

test_that("a .data missing the outcome column names the response, not its transformation", {
  m <- transformed_response_models()
  outcome_mod <- transformed_response_outcome("log(y)", m$dat, m$wts)
  without_y <- m$dat[setdiff(names(m$dat), "y")]

  for (se in c("mestimation", "linearization")) {
    err <- expect_error(
      ipw(m$ps_mod, outcome_mod, .data = without_y, se_method = se),
      class = "propensity_columns_exist_error"
    )
    msg <- conditionMessage(err)
    expect_match(msg, "\"y\"", fixed = TRUE)
    expect_no_match(msg, "log", fixed = TRUE)
  }
})

test_that("a response the formula's environment holds is still a missing column", {
  # The response is evaluated rather than looked up, with the formula's own
  # environment as the enclosure, so a `.data` missing the column reaches
  # whatever that environment happens to bind and estimates on it with nothing
  # signaled. What prevents that is the assert that runs first in the same
  # branch, over a required set the response's own variables are part of. The
  # pins above cannot see that: their formula environments hold nothing to find,
  # so dropping the response from the required set leaves them reporting a
  # missing object instead of a missing column and still failing. Only an
  # enclosure that answers shows the silent case, so one is stocked here.
  m <- transformed_response_models()

  # The response is bound in the formula's environment and nowhere the fit
  # reads: `model.frame()` finds `y` in `data` first, so the model is fit on the
  # real column and only the enclosure holds the decoy.
  env <- new.env(parent = environment())
  env$y <- rev(m$dat$y)
  fmla <- y ~ z
  environment(fmla) <- env
  outcome_mod <- lm(fmla, data = m$dat, weights = m$wts)

  without_y <- m$dat[setdiff(names(m$dat), "y")]

  # The premise: the enclosure really does answer for the response, and with
  # values the fit never saw, so an extraction that ran before the assert would
  # give an answer rather than fail.
  expect_identical(
    fmla_eval_left(outcome_mod, without_y),
    rev(m$dat$y)
  )

  for (se in c("mestimation", "linearization")) {
    err <- expect_error(
      ipw(m$ps_mod, outcome_mod, .data = without_y, se_method = se),
      class = "propensity_columns_exist_error"
    )
    expect_match(conditionMessage(err), "\"y\"", fixed = TRUE)
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
  estimates <- res$estimates[
    !names(res$estimates) %in%
      c("effect", "contrast")
  ]
  rownames(estimates) <- paste(
    res$estimates$effect,
    res$estimates$contrast
  )
  fixed <- capture.output(
    printCoefmat(estimates, has.Pvalue = TRUE, cs.ind = 1:2, tst.ind = 3)
  )
  current <- capture.output(
    printCoefmat(estimates, has.Pvalue = TRUE, cs.ind = 2:3, tst.ind = 4)
  )

  data_rows <- function(lines) grep("^(mean|rd|log)", lines, value = TRUE)
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

  # The comparison this route reports is the shared one, phrased the same way as
  # on the M-estimation route; its wording is pinned in test-ipw-mestimation.R.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "recomputed from `wt_mod`", fixed = TRUE)
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

  # The comparison this route reports is the shared one, phrased the same way as
  # on the M-estimation route; its wording is pinned in test-ipw-mestimation.R.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "recomputed from `wt_mod`", fixed = TRUE)
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
# The pins run under `with_always_deprecated()` rather than relying on the
# default once-per-session behavior, so no pin depends on being the first to
# reach the `ipw(ps_link)` id, nor on the environment the runner evaluates the
# file in.
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

  # Both directions, so the guard cannot be written as a rule about logit.
  cases <- list(
    list(fit = "logit", ps_link = "probit"),
    list(fit = "logit", ps_link = "cloglog"),
    list(fit = "probit", ps_link = "logit")
  )

  with_always_deprecated({
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

  with_always_deprecated({
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
})

test_that("linearization deprecates a ps_link equal to the model's link", {
  dat <- se_method_data()

  # Naming the link the model was already fit with asks for exactly what the
  # default gives, so it is deprecated rather than rejected. The warning is the
  # whole of the change: the estimates must still match the call that omits the
  # argument, value for value.
  with_always_deprecated({
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
        ipw(
          mods$ps_mod,
          mods$outcome_mod,
          se_method = "linearization"
        )$estimates
      )
    }
  })
})

test_that("mestimation deprecates a ps_link equal to the model's link", {
  dat <- se_method_data()
  mods <- se_link_fit(dat, "logit")

  with_always_deprecated({
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
})

test_that("omitting ps_link deprecates nothing on either path", {
  dat <- se_method_data()
  mods <- se_link_fit(dat, "logit")

  # Resolving the model's own link is the default rather than a use of the
  # argument, so neither the omitted form nor an explicit NULL may warn, even at
  # the loudest verbosity. Forcing the signal unconditionally is what keeps the
  # silence asserted here evidence of the default resolution rather than of a
  # deprecation another test already consumed.
  with_always_deprecated({
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

# The same weights against the codings a user is likelier to hold than a named
# factor. The guard compares levels as characters, so a recorded focal of `0`
# against a 0/1 column and one of `FALSE` against a logical column are the
# lower-level cases the factor fixture states in words: neither is the second
# sorted level, and both must be refused.
se_method_outcome_focal_lower <- function(dat, ps_mod, exposure_name, wt_fun) {
  exposure <- dat[[exposure_name]]
  lower <- sort(unique(exposure))[[1]]
  e_upper <- predict(ps_mod, type = "response")
  wts <- wt_fun(
    1 - e_upper,
    exposure,
    exposure_type = "binary",
    .focal_level = lower
  )
  glm(
    stats::reformulate(exposure_name, response = "y"),
    data = dat,
    family = quasibinomial(),
    weights = wts
  )
}

test_that("both paths reject att weights focal on 0 for a 0/1 exposure", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_focal_lower(dat, ps_mod, "z", wt_att)

  expect_identical(
    attr(
      stats::model.weights(stats::model.frame(outcome_mod)),
      "focal_category"
    ),
    0L
  )

  for (method in c("linearization", "mestimation")) {
    err <- expect_error(
      ipw(ps_mod, outcome_mod, .data = dat, se_method = method),
      class = "propensity_focal_level_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "\"0\"", fixed = TRUE)
    expect_match(msg, "\"1\"", fixed = TRUE)
  }
})

test_that("both paths reject atu weights focal on FALSE for a logical exposure", {
  dat <- se_method_data()
  dat$zl <- as.logical(dat$z)
  ps_mod <- glm(zl ~ x1 + x2, data = dat, family = binomial())
  outcome_mod <- se_method_outcome_focal_lower(dat, ps_mod, "zl", wt_atu)

  expect_identical(
    attr(
      stats::model.weights(stats::model.frame(outcome_mod)),
      "focal_category"
    ),
    FALSE
  )

  for (method in c("linearization", "mestimation")) {
    err <- expect_error(
      ipw(ps_mod, outcome_mod, .data = dat, se_method = method),
      class = "propensity_focal_level_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "\"FALSE\"", fixed = TRUE)
    expect_match(msg, "\"TRUE\"", fixed = TRUE)
  }
})

test_that("both paths reject a lower focal level without .data", {
  # The exposure comes from the propensity model's own frame when `.data` is
  # absent, so the guard has to reach the same values by that route.
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_focal_lower(dat, ps_mod, "z", wt_att)

  for (method in c("linearization", "mestimation")) {
    expect_error(
      ipw(ps_mod, outcome_mod, se_method = method),
      class = "propensity_focal_level_error"
    )
  }
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

test_that("linearization reports separation ahead of its estimand restriction", {
  dat <- se_separation_data()
  ps_mod <- suppressWarnings(glm(
    z ~ x1,
    data = dat,
    family = binomial(),
    control = glm.control(maxit = 200, epsilon = 1e-14)
  ))

  # atu and entropy have no linearization influence functions and are refused
  # for that. On a separated fit that refusal is the lesser truth: the model has
  # no finite maximum likelihood estimate, so no standard error method would
  # help, and a user told to switch to mestimation would only meet the same
  # separation there. The guard therefore runs first.
  for (est in c("atu", "entropy")) {
    wt_fun <- if (est == "atu") wt_atu else wt_entropy
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_fun(ps_mod)
    )
    outcome_mod <- suppressWarnings(glm(
      y ~ z,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    ))

    expect_error(
      ipw(ps_mod, outcome_mod, se_method = "linearization"),
      class = "propensity_ipw_separation_error"
    )
  }
})

test_that("the saturation guard leaves a missing score out of its count", {
  # A missing propensity score is missing rather than saturated, so it is left
  # out of the count. Counting it in would put an `NA` into the comparison the
  # guard branches on and stop the call with base R's report about a missing
  # value, on a fit whose scores never reached the boundary at all.
  #
  # No route through `ipw()` reaches this today. On the `.data` route a missing
  # covariate is refused by the completeness guard well before the scores are
  # predicted, and on the model-frame route the prediction is made over the
  # design the fit kept, which has the incomplete rows already dropped. The
  # exclusion is held against a change to either, so it is pinned where it is
  # rather than through a call that cannot produce it.
  expect_true(check_ipw_ps_saturation(c(2, NA, -2), link = "logit"))

  err <- expect_error(
    check_ipw_ps_saturation(c(50, NA, -2), link = "logit"),
    class = "propensity_ipw_separation_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "for 1 observation,", fixed = TRUE)
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

test_that("the linearization separation error reads in the user's terms", {
  dat <- se_separation_data()
  mods <- se_separation_models(dat)

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization"),
    transform = mask_saturated_count
  )
})

# ---- stabilization scores on the M-estimation path --------------------------

# RD standard error from the M-estimation path for the supplied weights.
se_score_mest_rd <- function(dat, ps_mod, wts) {
  res <- ipw(
    ps_mod,
    se_score_outcome(dat, wts),
    .data = dat,
    se_method = "mestimation"
  )
  res$estimates$std.err[res$estimates$effect == "rd"]
}

test_that("mestimation accepts a scalar stabilization score and reproduces the unstabilized estimates", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- as.double(predict(ps_mod, type = "response"))

  unstabilized <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)
  scored <- wt_ate(
    ps,
    dat$z,
    exposure_type = "binary",
    .focal_level = 1,
    stabilize = TRUE,
    stabilization_score = 2.5
  )

  mest_estimates <- function(wts) {
    ipw(
      ps_mod,
      se_score_outcome(dat, wts),
      .data = dat,
      se_method = "mestimation"
    )$estimates
  }

  # The weight-consistency comparator rebuilds the weights the propensity scores
  # imply from the recorded score, so a scalar score passes it rather than being
  # read as a mismatch with the default stabilizer.
  scored_estimates <- expect_silent(mest_estimates(scored))

  # One constant scales every weight and every weight total identically, so it
  # cancels from the Hajek means the stacked system solves. The tolerance is the
  # root solver's, not a modelling difference: the standard errors agree to
  # eight figures, and the p-values, which amplify that difference, to five.
  expect_equal(scored_estimates, mest_estimates(unstabilized), tolerance = 1e-5)
})

test_that("mestimation uses a per-observation stabilization score", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- as.double(predict(ps_mod, type = "response"))

  # A score correlated with a propensity score covariate does not cancel from
  # the Hajek means, so it has to move the standard error rather than wash out.
  score <- exp(0.5 * dat$x1)
  scored <- wt_ate(
    ps,
    dat$z,
    exposure_type = "binary",
    .focal_level = 1,
    stabilize = TRUE,
    stabilization_score = score
  )
  unstabilized <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)

  scored_mest <- expect_silent(se_score_mest_rd(dat, ps_mod, scored))
  unstabilized_mest <- se_score_mest_rd(dat, ps_mod, unstabilized)
  expect_gt(
    abs(scored_mest - unstabilized_mest) / unstabilized_mest,
    0.05
  )

  # The two standard error paths differ by their finite-sample corrections
  # rather than by what they make of the score, so they have to sit as close
  # together for the scored analysis as they do without a score at all.
  scored_gap <- abs(scored_mest - se_score_lin_rd(dat, ps_mod, scored)) /
    scored_mest
  unstabilized_gap <- abs(
    unstabilized_mest - se_score_lin_rd(dat, ps_mod, unstabilized)
  ) /
    unstabilized_mest
  expect_equal(scored_gap, unstabilized_gap, tolerance = 1e-2)
})

# ---- the shapes the propensity response guard rejects ------------------------
#
# The guard has two prongs and they describe two different models. A matrix
# response is the `cbind(successes, failures)` shape covered above. A left-hand
# side that is a call but evaluates to one column, `factor(z)` say, is not a
# matrix at all, and the cbind wording sent that user looking for a matrix they
# never wrote.

se_method_call_lhs_models <- function() {
  dat <- se_method_data()
  ps_mod <- glm(factor(z) ~ x1 + x2, data = dat, family = binomial())
  wts <- se_method_ate_wts(dat, se_method_ps_mod(dat))
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, dat = dat)
}

test_that("a call-form propensity response is rejected as an expression", {
  m <- se_method_call_lhs_models()

  # the fixture is the shape it claims to be: a call, and not a matrix
  expect_gt(length(as.character(formula(m$ps_mod)[[2]])), 1L)
  expect_false(is.matrix(stats::model.response(stats::model.frame(m$ps_mod))))

  for (se in c("mestimation", "linearization")) {
    err <- expect_error(
      ipw(m$ps_mod, m$outcome_mod, .data = m$dat, se_method = se),
      class = "propensity_ipw_response_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "factor(z)", fixed = TRUE)
    expect_false(grepl("matrix", msg, fixed = TRUE))
    expect_false(grepl("cbind", msg, fixed = TRUE))
  }
})

test_that("the call-form propensity response error reads in the user's terms", {
  m <- se_method_call_lhs_models()
  expect_snapshot(
    error = TRUE,
    ipw(m$ps_mod, m$outcome_mod, .data = m$dat)
  )
})

# ---- the outcome model the linearization path has nothing to compare ---------
#
# The first linearization outcome guard rejects any model whose terms are not
# the exposure alone. Two shapes reach it and they are opposites: a model that
# carries the exposure and more, and a model that carries the exposure not at
# all. `y ~ 1` has no terms whatever, and reporting it as "adjusted for terms
# beyond z" described a model the user did not fit.

se_method_outcome_intercept_only <- function(dat, ps_mod) {
  wts <- se_method_ate_wts(dat, ps_mod)
  glm(y ~ 1, data = dat, family = quasibinomial(), weights = wts)
}

# Adjusted for a covariate and missing the exposure at the same time. The
# missing exposure is the defect that has to be reported.
se_method_outcome_covariate_only <- function(dat, ps_mod) {
  wts <- se_method_ate_wts(dat, ps_mod)
  glm(y ~ x1, data = dat, family = quasibinomial(), weights = wts)
}

test_that("an intercept-only outcome model is reported as missing the exposure", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_intercept_only(dat, ps_mod)

  # the fixture is the shape it claims to be: no terms at all
  expect_length(attr(stats::terms(outcome_mod), "term.labels"), 0L)

  err <- expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "does not include the exposure", fixed = TRUE)
  expect_false(grepl("adjusted for terms beyond", msg, fixed = TRUE))
  expect_match(msg, "mestimation", fixed = TRUE)
})

test_that("an outcome model of the covariates alone is reported as missing the exposure", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_covariate_only(dat, ps_mod)

  err <- expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "does not include the exposure", fixed = TRUE)
  expect_false(grepl("adjusted for terms beyond", msg, fixed = TRUE))
})

test_that("a genuinely adjusted outcome model keeps the adjusted wording", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_adjusted_lm(dat, ps_mod)

  err <- expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization"),
    class = "propensity_method_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "adjusted for terms beyond", fixed = TRUE)
  expect_false(grepl("does not include the exposure", msg, fixed = TRUE))
})

test_that("the intercept-only outcome rejection reads in the user's terms", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_intercept_only(dat, ps_mod)

  expect_snapshot(
    error = TRUE,
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
  )
})

# ---- the contrast column a binary exposure names its one pair in -------------

test_that("ipw() binary names its contrast column contrast under both routes", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  # A binary exposure reports the counterfactual mean at each level and one pair
  # of levels, and names both in a column called `contrast`, the way a
  # categorical exposure does. Both routes are pinned here because linearization
  # is available for a binary exposure alone, a categorical or continuous
  # exposure refusing it, so this is the only fixture that can hold the two
  # routes to the same table shape.
  results <- list(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "mestimation"),
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
  )
  for (res in results) {
    expect_true("contrast" %in% names(res$estimates))
    expect_false("comparison" %in% names(res$estimates))
    expect_identical(
      names(res$estimates)[seq(1, 2)],
      c("effect", "contrast")
    )
    expect_identical(res$estimates$contrast, c("0", "1", rep("1 vs 0", 3)))
  }
})

# ---- se_method = "robust": the HC0 diagnostic -------------------------------
#
# The diagnostic route reports the sandwich the weighted outcome model computes
# for itself, with the weights entering as known constants: it is the
# linearization route with the correction for having estimated the propensity
# score dropped. The point estimates are therefore the linearization ones and
# only the standard errors move. The oracle below is built from the bread and
# the meat by hand rather than read off a package, so the route is checked
# against the definition of the thing it claims to report.

# The HC0 covariance of a weighted model's coefficients: the inverse of the
# weighted information, the outer product of the weighted score contributions,
# and the inverse of the information again. The prior weights enter as known
# constants, which is the whole of what this diagnostic leaves out. Written
# through the family's variance function and link derivative so the
# quasibinomial and the linear outcome models are read the same way.
hc0_coef_vcov <- function(model) {
  design <- stats::model.matrix(model)
  y <- as.numeric(stats::model.response(stats::model.frame(model)))
  w <- as.numeric(stats::weights(model))
  mu <- as.numeric(stats::fitted(model))

  if (inherits(model, "glm")) {
    dmu <- model$family$mu.eta(as.numeric(model$linear.predictors))
    variance <- model$family$variance(mu)
  } else {
    dmu <- rep(1, length(mu))
    variance <- rep(1, length(mu))
  }

  information <- crossprod(design, design * (w * dmu^2 / variance))
  score <- w * (y - mu) * dmu / variance
  meat <- crossprod(design, design * score^2)
  bread <- solve(information)

  bread %*% meat %*% bread
}

# The delta-method covariance of the reported effects, in the row order the
# estimates table reports them: the counterfactual mean at each exposure level,
# their difference, and, where the outcome model is not linear, the log risk
# ratio and the log odds ratio. The outcome model is the exposure alone with an
# intercept, so the two means are the two fitted cells.
hc0_effect_vcov <- function(model, linear = FALSE) {
  b <- unname(stats::coef(model))
  linkinv <- if (inherits(model, "glm")) model$family$linkinv else identity
  mu_eta <- if (inherits(model, "glm")) {
    model$family$mu.eta
  } else {
    function(eta) 1
  }

  mu0 <- linkinv(b[[1]])
  mu1 <- linkinv(b[[1]] + b[[2]])
  d0 <- mu_eta(b[[1]])
  d1 <- mu_eta(b[[1]] + b[[2]])

  jacobian <- rbind(
    c(d0, 0),
    c(d1, d1),
    c(d1 - d0, d1)
  )

  if (!linear) {
    odds0 <- mu0 * (1 - mu0)
    odds1 <- mu1 * (1 - mu1)
    jacobian <- rbind(
      jacobian,
      c(d1 / mu1 - d0 / mu0, d1 / mu1),
      c(d1 / odds1 - d0 / odds0, d1 / odds1)
    )
  }

  jacobian %*% hc0_coef_vcov(model) %*% t(jacobian)
}

test_that("ipw(se_method = 'robust') reports the linearization estimates with its own SEs", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  res_r <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "robust")
  res_l <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")

  expect_s3_class(res_r, "ipw")
  expect_equal(res_r$se_method, "robust")

  # No system is solved here, as none is solved on the linearization route.
  expect_null(res_r$fit)

  # The estimates table keeps the shared column contract and the shared rows.
  expect_named(res_r$estimates, estimates_columns)
  expect_identical(res_r$estimates$effect, res_l$estimates$effect)
  expect_identical(res_r$estimates$contrast, res_l$estimates$contrast)
  expect_equal(
    res_r$estimates$estimate,
    res_l$estimates$estimate,
    tolerance = 1e-10
  )

  # Only the standard errors move, and on this fixture they move up: the
  # correction the diagnostic drops is what removes the variance of having
  # estimated the propensity score from the ATE.
  expect_false(isTRUE(all.equal(
    res_r$estimates$std.err,
    res_l$estimates$std.err
  )))
  expect_true(all(res_r$estimates$std.err > res_l$estimates$std.err))
})

test_that("robust standard errors are the hand-built HC0 sandwich of the outcome model", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "robust")
  oracle <- hc0_effect_vcov(outcome_mod)

  expect_equal(res$estimates$std.err, sqrt(diag(oracle)), tolerance = 1e-8)

  # The meat is divided by n rather than by n - 1, so what is reported is the
  # HC0 sandwich itself rather than a finite-sample rescaling of it.
  n <- nrow(dat)
  expect_false(isTRUE(all.equal(
    res$estimates$std.err,
    sqrt(diag(oracle) * n / (n - 1))
  )))

  # The covariance the result attaches is that same sandwich, so the reported
  # standard errors and the block a caller reads off the frame describe one
  # variance rather than two.
  expect_equal(
    unname(attr(res$estimates, "ipw_vcov")),
    oracle,
    tolerance = 1e-8
  )

  # The test statistic is built from the reported standard error, as on every
  # other route.
  expect_equal(
    res$estimates$z,
    res$estimates$estimate / res$estimates$std.err,
    tolerance = 1e-10
  )
})

test_that("robust standard errors of a linear outcome model are the hand-built HC0 sandwich", {
  dat <- se_method_data_cont()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_marginal_lm(dat, ps_mod)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "robust")
  oracle <- hc0_effect_vcov(outcome_mod, linear = TRUE)

  # A linear outcome model reports the two means and their difference, the
  # ratio measures belonging to a risk.
  expect_identical(res$estimates$effect, c("mean", "mean", "diff"))
  expect_equal(res$estimates$std.err, sqrt(diag(oracle)), tolerance = 1e-8)
  expect_equal(
    unname(attr(res$estimates, "ipw_vcov")),
    oracle,
    tolerance = 1e-8
  )
})

test_that("a robust result marks its standard errors as diagnostic when printed", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  # The estimate table needs more than the 80 columns testthat pins or it wraps.
  withr::local_options(width = 120)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "robust")

  # The mark is a class on the result, which is what gives the printed output a
  # method of its own. The result is an `ipw` to everything else.
  expect_identical(class(res)[[1]], "ipw_diagnostic_se")
  expect_s3_class(res, "ipw")

  printed <- utils::capture.output(print(res))
  marked <- grepl("^Standard errors: robust", printed)
  expect_identical(sum(marked), 1L)
  expect_match(printed[marked], "diagnostic")

  # The two methods that account for the estimated weights carry neither the
  # class nor the line, so the mark reads as the exception it is.
  for (method in c("mestimation", "linearization")) {
    other <- ipw(ps_mod, outcome_mod, .data = dat, se_method = method)
    expect_false(inherits(other, "ipw_diagnostic_se"))
    expect_false(any(grepl(
      "^Standard errors:",
      utils::capture.output(print(other))
    )))
  }
})

test_that("tidy() and as.data.frame() of a robust result carry the diagnostic mark", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "robust")

  # The mark travels as an attribute rather than as a column, the table being
  # the one every other route reports.
  expect_identical(attr(as.data.frame(res), "ipw_se_diagnostic"), "robust")
  expect_identical(attr(tidy(res), "ipw_se_diagnostic"), "robust")
  expect_identical(
    attr(tidy(res, conf.int = TRUE), "ipw_se_diagnostic"),
    "robust"
  )
  expect_identical(
    attr(as.data.frame(res, exponentiate = TRUE), "ipw_se_diagnostic"),
    "robust"
  )

  reference <- ipw(
    ps_mod,
    outcome_mod,
    .data = dat,
    se_method = "linearization"
  )
  expect_identical(names(tidy(res)), names(tidy(reference)))

  for (method in c("mestimation", "linearization")) {
    other <- ipw(ps_mod, outcome_mod, .data = dat, se_method = method)
    expect_null(attr(as.data.frame(other), "ipw_se_diagnostic"))
    expect_null(attr(tidy(other), "ipw_se_diagnostic"))
  }
})

test_that("a robust result has no conditional reading", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  res <- ipw(ps_mod, outcome_mod, .data = dat, se_method = "robust")

  # The conditional reading reports the covariance the joint estimation of the
  # weights and the outcome implies, and this route estimates nothing jointly.
  # The refusal names the method the result records rather than naming the
  # linearization method the result was not fit with.
  expect_error(
    tidy(res, effects = "conditional"),
    class = "propensity_no_conditional_vcov_error",
    regexp = "robust"
  )
})

test_that("ipw() refuses .by with robust standard errors", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_ate(dat, ps_mod)

  # The stratum effects and their contrasts are parameters of the stacked
  # system, which this route no more solves than the linearization route does.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, .by = x2, se_method = "robust"),
    class = "propensity_ipw_by_method_error",
    regexp = "robust"
  )
})

test_that("robust rejects the atu estimand as linearization does", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_atu(ps_mod))
  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  # The diagnostic route is the linearization route with the correction
  # dropped, so it accepts the estimands that route accepts and no others,
  # leaving one set of requirements to learn rather than two. The refusal names
  # the method that was asked for.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "robust"),
    class = "propensity_method_error",
    regexp = "robust"
  )
})

test_that("robust rejects a covariate-adjusted outcome model as linearization does", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  outcome_mod <- se_method_outcome_adjusted_glm(dat, ps_mod)

  # The reported estimates are the Hajek means of a two-cell outcome model, and
  # an adjusted model reports something else, whichever variance is put beside
  # them.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "robust"),
    class = "propensity_method_error",
    regexp = "robust"
  )
})
