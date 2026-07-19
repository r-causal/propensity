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

# ---- probit and cloglog linearization variance ------------------------------
#
# The linearization ps correction adds H' A^{-1} s_i to the influence values,
# where s_i = x_i (Z_i - p_i) / (p_i (1 - p_i) g'(p_i)) is the score of a
# binomial GLM with link g. The package generalizes the link in the weight
# derivatives and in the Fisher information A, but hardcodes the score itself to
# the canonical logit form x_i (Z_i - p_i), omitting the factor
# 1 / (p (1 - p) g'(p)), which equals 1 only for logit. These reference tests pin
# the generalized score correction; the probit and cloglog cases fail until the
# correction is implemented.

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
# the canonical-logit form the package currently hardcodes, which is the value
# the probit and cloglog cases return today.
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
# information) from a fitted ps model and the outcome data.
lin_se_correct <- function(ps_mod, z, y, w) {
  pieces <- lin_ate_pieces(ps_mod, z)
  lin_se_rd(
    pieces$X,
    pieces$p,
    z,
    y,
    as.double(w),
    pieces$dwdeta,
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
  # package agree, validating the helper against the current implementation.
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

  # The corrected linearization agrees with the sandwich to well under a percent
  # (about 0.02 percent here); the current mis-scaled formula is about 1.6
  # percent off, so this band separates the two.
  expect_lt(abs(lin_se - mest_se) / mest_se, 0.005)
})

# ---- outcome-family validation on the linearization path --------------------

test_that("the linearization path rejects an unsupported outcome family", {
  dat <- se_method_data()
  ps_mod <- se_method_ps_mod(dat)
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)

  # A low-rate count outcome; the marginal outcome model passes the exposure-only
  # guard, so the family rejection must fire even on a marginal model. Today the
  # path returns a finite but meaningless log(or) from count-scale marginal means.
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

  # Today Y - mu is computed on the raw factor, giving NA std.err with
  # "'-' not meaningful for factors" warnings. The linearization branch extracts
  # the outcome two ways that fail identically: `.data[[outcome_name]]` when
  # .data is supplied, and fmla_extract_left_vctr when it is not. Both must
  # return finite std.err equal to the numeric fit with no warnings.
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
# request currently reaches derive_weights, whose rlang::arg_match raises a bare,
# unclassed rlang error attributed to the internal function and claiming the
# estimand "must be one of ate, att, ato, atm" (misleading: atu and entropy are
# valid estimands, just unsupported by this SE method). These tests pin a classed
# propensity_method_error directing to se_method = "mestimation", matching the
# categorical and continuous paths. The four supported estimands on linearization
# are covered by test-ipw.R ("variance works"); atu and entropy on mestimation by
# test-ipw-mestimation.R ("ipw_mestimation runs atu and entropy with finite SEs").

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
