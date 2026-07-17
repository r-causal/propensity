test_that("deli integration: M-estimation solves and returns sandwich variance", {
  skip_if_not_installed("deli")

  set.seed(1)
  y <- rnorm(50, mean = 3, sd = 2)

  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- deli::MEstimator(stacked_equations = psi, init = c(mu = 0))
  m <- deli::estimate(m)

  # the M-estimator for a mean solves to the sample mean
  expect_equal(unname(coef(m)), mean(y), tolerance = 1e-8)
  expect_named(coef(m), "mu")

  # empirical sandwich variance for a mean equals the finite-sample
  # variance of y divided by n
  se_sandwich <- sqrt(diag(vcov(m)))
  se_closed_form <- sqrt(mean((y - mean(y))^2) / length(y))
  expect_equal(unname(se_sandwich), se_closed_form, tolerance = 1e-6)
})

# ---- binary M-estimation engine --------------------------------------------
#
# Tests for the internal binary engine: the spec constructor ipw_spec_binary()
# (which builds the design-contract spec FROM FITTED MODELS, reusing the same
# extraction that ipw() uses), the preflight ipw_check_weight_consistency(), and
# the driver ipw_mestimation() (which runs the deli fit from the seeded init and
# returns list(estimates, fit)). Unlike the psi-layer tests in test-ipw-psi.R,
# these construct specs from fitted glm/lm models rather than hand-built lists.

# ---- data simulator ---------------------------------------------------------

sim_binary <- function(seed = 2024, n = 800) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.6 * x2))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * x1 - 0.3 * x2))
  yc <- 1.5 + 0.8 * z + 0.6 * x1 - 0.4 * x2 + rnorm(n)
  data.frame(x1, x2, z, y, yc)
}

# ---- model fitting: build the two models the engine consumes ----------------

# Fit the propensity score model and a weighted outcome model, returning both
# plus the psw weights. The weights are supplied to the outcome model as a psw
# object (unless strip_weights = TRUE) so the estimand and stabilization status
# survive into model.frame() and can be auto-detected by the spec constructor.
fit_binary_models <- function(
  dat,
  estimand,
  ps_link = "logit",
  outcome_family = "binomial",
  stabilize = FALSE,
  strip_weights = FALSE
) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial(link = ps_link))
  wt_fun <- switch(
    estimand,
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )
  wt_args <- list(ps_mod)
  if (estimand == "ate") {
    wt_args$stabilize <- stabilize
  }
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    do.call(wt_fun, wt_args)
  )
  model_wts <- if (strip_weights) as.double(wts) else wts

  outcome_var <- if (outcome_family == "binomial") "y" else "yc"
  out_fmla <- stats::reformulate(c("z", "x1"), response = outcome_var)
  if (outcome_family == "binomial") {
    # Tighten IRLS convergence so the reference outcome model sits at its MLE
    # to well below the 1e-8 point-estimate comparison tolerance. The default
    # glm tolerance can stop short of the root for some weight schemes (atm in
    # particular), while the M-estimator solves the score to machine precision.
    outcome_mod <- glm(
      out_fmla,
      data = dat,
      family = quasibinomial(),
      weights = model_wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  } else {
    outcome_mod <- lm(out_fmla, data = dat, weights = model_wts)
  }

  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# Direct plug-in reference: g-computation on the weighted outcome model, matching
# how estimate_marginal_means() constructs counterfactual predictions in R/ipw.R.
# Returns a named vector of the contrast point estimates.
plugin_contrasts <- function(outcome_mod, dat, outcome_family = "binomial") {
  d1 <- dat
  d0 <- dat
  d1$z <- 1
  d0$z <- 0
  mu1 <- mean(predict(outcome_mod, newdata = d1, type = "response"))
  mu0 <- mean(predict(outcome_mod, newdata = d0, type = "response"))
  if (outcome_family == "binomial") {
    c(
      rd = mu1 - mu0,
      "log(rr)" = log(mu1) - log(mu0),
      "log(or)" = qlogis(mu1) - qlogis(mu0)
    )
  } else {
    c(diff = mu1 - mu0)
  }
}

# ---- assertion helper -------------------------------------------------------

# Compare an engine estimates table to a named reference vector, matching by the
# effect label so row order is irrelevant.
expect_estimate_match <- function(estimates, reference, tolerance = 1e-8) {
  got <- estimates$estimate[match(names(reference), estimates$effect)]
  expect_equal(got, unname(reference), tolerance = tolerance)
}

# ---- point-estimate parity with ipw() ---------------------------------------

test_that("ipw_mestimation point estimates match ipw() for binomial outcomes", {
  skip_if_not_installed("deli")
  for (est in c("ate", "att", "ato", "atm")) {
    dat <- sim_binary()
    mods <- fit_binary_models(dat, est)
    # pin the linearization comparator so this stays a cross-method check once
    # mestimation becomes the ipw() default in a later issue
    ref <- ipw(
      mods$ps_mod,
      mods$outcome_mod,
      se_method = "linearization"
    )$estimates
    spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
    got <- ipw_mestimation(spec)$estimates
    expect_estimate_match(
      got,
      stats::setNames(ref$estimate, ref$effect),
      tolerance = 1e-8
    )
  }
})

test_that("ipw_mestimation point estimates match ipw() for gaussian outcomes", {
  skip_if_not_installed("deli")
  for (est in c("ate", "att", "ato", "atm")) {
    dat <- sim_binary()
    mods <- fit_binary_models(dat, est, outcome_family = "gaussian")
    ref <- ipw(
      mods$ps_mod,
      mods$outcome_mod,
      se_method = "linearization"
    )$estimates
    spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
    got <- ipw_mestimation(spec)$estimates
    expect_estimate_match(
      got,
      stats::setNames(ref$estimate, ref$effect),
      tolerance = 1e-8
    )
  }
})

# ---- standard-error cross-validation against linearization ------------------

# The linearization is a first-order (delta-method) approximation of the two-step
# variance, while M-estimation returns the exact empirical sandwich. The two are
# asymptotically equivalent but differ at finite samples by a realization-
# dependent amount. For most estimands that gap is well under 1 percent here, but
# the att log risk ratio sits around 2 percent at this sample size (the
# M-estimation value matches a nonparametric bootstrap; the linearization is the
# looser approximation). A 3 percent band keeps this a meaningful cross-method
# sanity check while accommodating that inherent gap.
test_that("ipw_mestimation std.err matches linearization within 3 percent", {
  skip_if_not_installed("deli")
  configs <- list(
    list(est = "ate", link = "logit", family = "binomial"),
    list(est = "ate", link = "probit", family = "binomial"),
    list(est = "ate", link = "cloglog", family = "binomial"),
    list(est = "att", link = "logit", family = "binomial"),
    list(est = "ato", link = "logit", family = "binomial"),
    list(est = "atm", link = "logit", family = "binomial"),
    list(est = "ate", link = "logit", family = "gaussian")
  )
  for (cfg in configs) {
    dat <- sim_binary()
    mods <- fit_binary_models(
      dat,
      cfg$est,
      ps_link = cfg$link,
      outcome_family = cfg$family
    )
    ref <- ipw(
      mods$ps_mod,
      mods$outcome_mod,
      se_method = "linearization"
    )$estimates
    spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
    got <- ipw_mestimation(spec)$estimates
    got_se <- got$std.err[match(ref$effect, got$effect)]
    rel <- abs(got_se - ref$std.err) / ref$std.err
    expect_true(
      all(rel < 0.03),
      label = paste0(
        cfg$est,
        "/",
        cfg$link,
        "/",
        cfg$family,
        " rel diff: ",
        paste(format(rel, digits = 3), collapse = ", ")
      )
    )
  }
})

# ---- atu and entropy: no ipw() comparator, use the direct plug-in -----------

test_that("ipw_mestimation estimates atu and entropy via plug-in with finite SEs", {
  skip_if_not_installed("deli")
  for (est in c("atu", "entropy")) {
    dat <- sim_binary()
    mods <- fit_binary_models(dat, est)
    spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
    got <- ipw_mestimation(spec)$estimates

    ref <- plugin_contrasts(mods$outcome_mod, dat)
    expect_estimate_match(got, ref, tolerance = 1e-8)

    expect_true(all(is.finite(got$std.err) & got$std.err > 0))
    expect_true(all(got$ci.lower < got$estimate & got$estimate < got$ci.upper))
  }
})

# ---- stabilized ate ---------------------------------------------------------

test_that("ipw_mestimation runs stabilized ate and matches the stabilized plug-in", {
  skip_if_not_installed("deli")
  dat <- sim_binary(seed = 99, n = 2000)
  mods <- fit_binary_models(dat, "ate", stabilize = TRUE)
  spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  got <- ipw_mestimation(spec)$estimates

  ref <- plugin_contrasts(mods$outcome_mod, dat)
  expect_estimate_match(got, ref, tolerance = 1e-8)
  expect_true(all(is.finite(got$std.err) & got$std.err > 0))
})

# ---- estimates table shape --------------------------------------------------

test_that("ipw_mestimation estimates table has the documented shape", {
  skip_if_not_installed("deli")
  dat <- sim_binary()

  mods <- fit_binary_models(dat, "ate")
  spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  got <- ipw_mestimation(spec, conf_level = 0.95)$estimates

  expect_named(
    got,
    c(
      "effect",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_equal(got$effect, c("rd", "log(rr)", "log(or)"))
  expect_true(all(got$conf.level == 0.95))

  # gaussian outcome collapses to a single difference row
  mods_g <- fit_binary_models(dat, "ate", outcome_family = "gaussian")
  spec_g <- ipw_spec_binary(mods_g$ps_mod, mods_g$outcome_mod)
  got_g <- ipw_mestimation(spec_g)$estimates
  expect_equal(got_g$effect, "diff")

  # a non-default confidence level rescales the interval half-width by qnorm
  got90 <- ipw_mestimation(spec, conf_level = 0.9)$estimates
  expect_true(all(got90$conf.level == 0.9))

  half95 <- (got$ci.upper - got$estimate) / got$std.err
  half90 <- (got90$ci.upper - got90$estimate) / got90$std.err
  expect_equal(half95, rep(qnorm(0.975), nrow(got)), tolerance = 1e-8)
  expect_equal(half90, rep(qnorm(0.95), nrow(got90)), tolerance = 1e-8)
  expect_true(all(half90 < half95))
})

# ---- fit object -------------------------------------------------------------

test_that("ipw_mestimation returns a converged deli fit with named theta", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  res <- ipw_mestimation(spec)
  fit <- res$fit

  co <- coef(fit)
  expect_true(all(is.finite(co)))
  expect_equal(length(co), length(ipw_theta_layout(spec)$init))

  nm <- names(co)
  expect_false(any(is.na(nm) | nm == ""))

  vc <- vcov(fit)
  expect_equal(dim(vc), c(length(co), length(co)))

  # the contrast rows of the solved theta equal the reported estimates
  got <- res$estimates
  expect_equal(unname(co[got$effect]), got$estimate, tolerance = 1e-8)
})

# ---- weight-consistency preflight -------------------------------------------

test_that("ipw_check_weight_consistency passes for consistent models", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  observed <- as.double(model.weights(model.frame(mods$outcome_mod)))
  expect_true(ipw_check_weight_consistency(spec, observed))
})

test_that("ipw_check_weight_consistency errors on weights from a different ps model", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  ps_a <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts_a <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_a)
  )
  outcome_mod <- glm(
    y ~ z + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts_a
  )
  # the spec's ps block comes from a different model than the one that produced
  # the outcome weights, so the recomputed weights will not match
  ps_b <- glm(z ~ x1, data = dat, family = binomial())
  spec <- ipw_spec_binary(ps_b, outcome_mod)
  observed <- as.double(model.weights(model.frame(outcome_mod)))
  expect_error(
    ipw_check_weight_consistency(spec, observed),
    class = "propensity_ipw_weights_mismatch_error"
  )
})

test_that("ipw_check_weight_consistency errors on manually doubled weights", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  doubled <- 2 * as.double(model.weights(model.frame(mods$outcome_mod)))
  expect_error(
    ipw_check_weight_consistency(spec, doubled),
    class = "propensity_ipw_weights_mismatch_error"
  )
})

test_that("ipw_spec_binary errors when the estimand cannot be determined", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  # plain numeric weights carry no estimand attribute, so with estimand = NULL
  # there is nothing to detect
  mods <- fit_binary_models(dat, "ate", strip_weights = TRUE)
  expect_error(
    ipw_spec_binary(mods$ps_mod, mods$outcome_mod, estimand = NULL),
    class = "propensity_error"
  )
})
