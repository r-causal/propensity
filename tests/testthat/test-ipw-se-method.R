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

test_that("mismatched weights error on the mestimation path only", {
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

  # The mestimation preflight recomputes the weights from the PS model and
  # detects the mismatch.
  expect_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "mestimation"),
    class = "propensity_ipw_weights_mismatch_error"
  )

  # The linearization path has no registry check, so mismatched weights are not
  # detectable there and estimation proceeds without error.
  expect_no_error(
    ipw(ps_mod, outcome_mod, .data = dat, se_method = "linearization")
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
