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
  strip_weights = FALSE,
  adjust = TRUE
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
  rhs <- if (adjust) c("z", "x1") else "z"
  out_fmla <- stats::reformulate(rhs, response = outcome_var)
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

# Fit the outcome model on the covariate alone, dropping the exposure, while
# keeping the ate weights the propensity score model implies. Used by the
# missing-exposure guard tests at the end of the file.
fit_outcome_without_exposure <- function(dat, wts) {
  glm(
    stats::reformulate("x1", response = "y"),
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
}

# Estimand tilt h(e) for the standardized g-computation plug-in. ate is the flat
# tilt (h = 1); a tilted estimand weights each unit's counterfactual predictions
# by h of its propensity score. Mirrors the tilt the M-estimation marginal-mean
# rows apply when standardizing to the target population.
plugin_tilt <- function(ps, estimand, n) {
  if (estimand == "ate") {
    return(rep(1, n))
  }
  if (is.null(ps)) {
    stop("plugin_tilt(): a tilted estimand requires propensity scores in `ps`.")
  }
  switch(
    estimand,
    att = ps,
    atu = 1 - ps,
    atm = pmin(ps, 1 - ps),
    ato = ps * (1 - ps),
    entropy = -ps * log(ps) - (1 - ps) * log(1 - ps)
  )
}

# Direct plug-in reference: g-computation on the weighted outcome model, matching
# how the M-estimation marginal-mean rows construct counterfactual predictions.
# The counterfactual predictions are standardized to the target population by the
# estimand tilt h(ps): for ate (h = 1, ps = NULL) this is the ordinary mean; for
# a tilted estimand it is the h-weighted mean. The outcome family is detected from
# the model. Returns a named vector of the contrast point estimates.
plugin_contrasts <- function(
  outcome_mod,
  data,
  exposure_name = "z",
  estimand = "ate",
  ps = NULL
) {
  d1 <- data
  d0 <- data
  d1[[exposure_name]] <- 1
  d0[[exposure_name]] <- 0
  m1 <- predict(outcome_mod, newdata = d1, type = "response")
  m0 <- predict(outcome_mod, newdata = d0, type = "response")
  h <- plugin_tilt(ps, estimand, nrow(data))
  mu1 <- weighted.mean(m1, h)
  mu0 <- weighted.mean(m0, h)
  linear <- if (inherits(outcome_mod, "glm")) {
    stats::family(outcome_mod)$family == "gaussian"
  } else {
    TRUE
  }
  if (linear) {
    c(diff = mu1 - mu0)
  } else {
    c(
      rd = mu1 - mu0,
      "log(rr)" = log(mu1) - log(mu0),
      "log(or)" = qlogis(mu1) - qlogis(mu0)
    )
  }
}

# ---- assertion helper -------------------------------------------------------

# Compare an engine estimates table to a named reference vector, matching by the
# effect label so row order is irrelevant.
expect_estimate_match <- function(estimates, reference, tolerance = 1e-8) {
  got <- estimates$estimate[match(names(reference), estimates$effect)]
  expect_equal(got, unname(reference), tolerance = tolerance)
}

# ---- point-estimate parity with the g-computation plug-in -------------------

test_that("ipw_mestimation ate point estimates match the g-computation plug-in for binomial outcomes", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  # The g-computation plug-in is the point-estimate oracle. Both standard error
  # methods report these same g-computation contrasts, so the plug-in pins the
  # engine estimates without routing through the linearization path, which is
  # restricted to exposure-only outcome models. ate is the flat tilt (h = 1);
  # the tilted estimands are pinned separately once standardization is in place.
  ref <- plugin_contrasts(mods$outcome_mod, dat, "z", estimand = "ate")
  spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  got <- ipw_mestimation(spec)$estimates
  expect_estimate_match(got, ref, tolerance = 1e-8)
})

test_that("ipw_mestimation ate point estimates match the g-computation plug-in for gaussian outcomes", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate", outcome_family = "gaussian")
  ref <- plugin_contrasts(mods$outcome_mod, dat, "z", estimand = "ate")
  spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  got <- ipw_mestimation(spec)$estimates
  expect_estimate_match(got, ref, tolerance = 1e-8)
})

# ---- standard-error cross-validation against linearization ------------------

# The linearization is a first-order (delta-method) approximation of the two-step
# variance, while M-estimation returns the exact empirical sandwich. The two are
# asymptotically equivalent but differ at finite samples by a realization-
# dependent amount. For most estimands that gap is well under 1 percent here, but
# the att log risk ratio sits around 2 percent at this sample size (the
# M-estimation value matches a nonparametric bootstrap; the linearization is the
# looser approximation). A 3 percent band keeps this a meaningful cross-method
# sanity check while accommodating that inherent gap. The outcome model is fit on
# the exposure alone, the domain where the linearization comparator is valid.
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
      outcome_family = cfg$family,
      adjust = FALSE
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

# ---- atu and entropy: engine runs with finite SEs ---------------------------

# The tilted point-estimate parity for these estimands with an adjusted outcome
# model is pinned in the skipped tilt-standardization test below; here the engine
# is only exercised for a converged fit with finite, positive standard errors and
# ordered intervals, which holds regardless of the marginal-mean standardization.
test_that("ipw_mestimation runs atu and entropy with finite SEs", {
  skip_if_not_installed("deli")
  for (est in c("atu", "entropy")) {
    dat <- sim_binary()
    mods <- fit_binary_models(dat, est)
    spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
    got <- ipw_mestimation(spec)$estimates

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

  ref <- plugin_contrasts(mods$outcome_mod, dat, "z", estimand = "ate")
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
  # (addressed by position, since theta names are not unique across blocks)
  got <- res$estimates
  contrast_idx <- ipw_theta_layout(spec)$idx$contrast
  expect_equal(unname(co[contrast_idx]), got$estimate, tolerance = 1e-8)
})

# ---- name collisions between covariates and contrast labels ------------------

test_that("ipw_mestimation addresses contrast rows by position, not by name", {
  skip_if_not_installed("deli")
  # A covariate literally named "rd" collides with the rd contrast label. The ps
  # and outcome coefficient blocks precede the contrast block in theta and keep
  # their term names, so name-based subsetting of coef(fit) would return the ps
  # coefficient on "rd" instead of the causal contrast. Positional indexing must
  # recover the true contrast.
  withr::local_seed(2024)
  n <- 800
  rd <- rnorm(n)
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * rd - 0.2 * x1))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * rd - 0.3 * x1))
  dat <- data.frame(rd, x1, z, y)

  ps_mod <- glm(z ~ rd + x1, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_mod)
  )
  outcome_mod <- glm(
    y ~ z + rd + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  spec <- ipw_spec_binary(ps_mod, outcome_mod)
  got <- ipw_mestimation(spec)$estimates

  ref <- plugin_contrasts(outcome_mod, dat, "z", estimand = "ate")
  expect_estimate_match(got, ref, tolerance = 1e-8)

  # the reported rd row is the contrast, not the ps model's coefficient on the
  # covariate that shares its name
  rd_estimate <- got$estimate[got$effect == "rd"]
  expect_false(isTRUE(all.equal(rd_estimate, unname(coef(ps_mod)["rd"]))))
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

test_that("the weight-consistency error names a focal level as a possible cause", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  doubled <- 2 * as.double(model.weights(model.frame(mods$outcome_mod)))

  # Hand-built weights carry no focal level for the guard in ipw() to read, so
  # the mismatch error is the only place a flipped focal level can be raised as
  # an explanation.
  expect_error(
    ipw_check_weight_consistency(spec, doubled),
    class = "propensity_ipw_weights_mismatch_error",
    regexp = "focal"
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

# ---- tilted marginal-mean standardization (att/atu/ato/atm/entropy) ----------
#
# The marginal-mean rows of the stacked system standardize to the estimand's
# target population, mu_a = sum_i h(e_i) m_a(x_i) / sum_i h(e_i). An unweighted
# mean of counterfactual predictions over the whole sample is correct only for
# ate; with a covariate-adjusted outcome model and a heterogeneous effect it
# reports the ate-type contrast for every estimand, which is why these tests
# combine both features. A saturated outcome model (y ~ z) has predictions
# constant within arm, so the tilt is a no-op on the contrast; that case is a
# live regression anchor.

# Heterogeneous-effect simulator with stored potential outcomes. The propensity
# e = plogis(1.5 x) and the effect y1 - y0 = 2 + 1.5 x are both increasing in x,
# so the tilted target populations (att up-weights high e, atu low e) carry
# genuinely different average effects than the ate.
sim_tilt <- function(seed = 501, n = 20000) {
  withr::local_seed(seed)
  x <- rnorm(n)
  e <- plogis(1.5 * x)
  z <- rbinom(n, 1, e)
  y0 <- x + rnorm(n)
  y1 <- y0 + (2 + 1.5 * x)
  y <- z * y1 + (1 - z) * y0
  data.frame(x, e, z, y0, y1, y)
}

test_that("mestimation diff matches the tilted oracle under a heterogeneous effect", {
  skip_if_not_installed("deli")
  dat <- sim_tilt()
  ps_mod <- glm(z ~ x, data = dat, family = binomial())
  e_hat <- as.double(predict(ps_mod, type = "response"))
  tau <- dat$y1 - dat$y0

  for (est in c("att", "atu", "ato")) {
    wt_fun <- switch(est, att = wt_att, atu = wt_atu, ato = wt_ato)
    wts <- withr::with_options(list(propensity.quiet = TRUE), wt_fun(ps_mod))
    om <- lm(y ~ z * x, data = dat, weights = wts)
    spec <- ipw_spec_binary(ps_mod, om)
    got <- ipw_mestimation(spec)$estimates
    est_diff <- got$estimate[got$effect == "diff"]
    se <- got$std.err[got$effect == "diff"]

    # exact: the engine reports the tilt-standardized g-computation from the
    # fitted model, weighted.mean(m1 - m0, h(e_hat))
    ref <- plugin_contrasts(om, dat, "z", estimand = est, ps = e_hat)[["diff"]]
    expect_equal(est_diff, ref, tolerance = 1e-6, label = est)

    # within sampling noise of the simulation oracle (h-tilted mean of the
    # simulated individual effects)
    oracle <- weighted.mean(tau, plugin_tilt(dat$e, est, nrow(dat)))
    expect_lt(abs(est_diff - oracle), 3 * se, label = est)

    # att and atu target populations are far from the ate
    if (est %in% c("att", "atu")) {
      ate_plain <- plugin_contrasts(om, dat, "z", estimand = "ate")[["diff"]]
      expect_gt(abs(est_diff - ate_plain), 10 * se, label = est)
    }
  }
})

test_that("mestimation att marginal means equal the tilt-standardized predictions at the solved theta", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_att(ps_mod))
  outcome_mod <- glm(
    y ~ z + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  spec <- ipw_spec_binary(ps_mod, outcome_mod)
  res <- ipw_mestimation(spec)

  co <- coef(res$fit)
  layout <- ipw_theta_layout(spec)
  th_ps <- co[layout$idx$ps]
  th_out <- co[layout$idx$out]
  th_mu <- co[layout$idx$mu]

  # rebuild e_hat and the counterfactual predictions from the solved blocks
  e_hat <- plogis(as.vector(spec$ps$X %*% th_ps))
  pred1 <- plogis(as.vector(spec$outcome$X_counterfactual$X1 %*% th_out))
  pred0 <- plogis(as.vector(spec$outcome$X_counterfactual$X0 %*% th_out))

  # att tilt is h(e) = e
  expect_equal(
    unname(th_mu[[1]]),
    weighted.mean(pred1, e_hat),
    tolerance = 1e-8
  )
  expect_equal(
    unname(th_mu[[2]]),
    weighted.mean(pred0, e_hat),
    tolerance = 1e-8
  )
})

test_that("mestimation tilted point estimates match the tilt-standardized plug-in for binomial outcomes", {
  skip_if_not_installed("deli")
  for (est in c("att", "atu", "atm", "ato", "entropy")) {
    dat <- sim_binary()
    mods <- fit_binary_models(dat, est)
    e_hat <- as.double(predict(mods$ps_mod, type = "response"))
    ref <- plugin_contrasts(
      mods$outcome_mod,
      dat,
      "z",
      estimand = est,
      ps = e_hat
    )
    spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
    got <- ipw_mestimation(spec)$estimates
    expect_estimate_match(got, ref, tolerance = 1e-8)
  }
})

test_that("mestimation att with a saturated outcome model matches the tilt-standardized plug-in", {
  skip_if_not_installed("deli")
  # A saturated outcome model has predictions constant within arm, so the att
  # tilt is a no-op on the contrast and the engine matches the tilt-standardized
  # plug-in whether or not the standardization is applied. This anchors the
  # exposure-only case.
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "att", adjust = FALSE)
  e_hat <- as.double(predict(mods$ps_mod, type = "response"))
  ref <- plugin_contrasts(
    mods$outcome_mod,
    dat,
    "z",
    estimand = "att",
    ps = e_hat
  )
  spec <- ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  got <- ipw_mestimation(spec)$estimates
  expect_estimate_match(got, ref, tolerance = 1e-8)
})

test_that("mestimation categorical att matches the tilt-standardized plug-in", {
  skip_if_not_installed("deli")
  skip_if_not_installed("nnet")
  withr::local_seed(913)
  n <- 900
  x <- rnorm(n)
  eta_b <- -0.3 + 0.9 * x
  eta_c <- -0.2 - 0.8 * x
  den <- 1 + exp(eta_b) + exp(eta_c)
  u <- runif(n)
  pa <- 1 / den
  pb <- exp(eta_b) / den
  a <- ifelse(u < pa, "a", ifelse(u < pa + pb, "b", "c"))
  a <- factor(a, levels = c("a", "b", "c"))
  y <- rbinom(
    n,
    1,
    plogis(-0.4 + 0.6 * (a == "b") + 0.9 * (a == "c") + 0.5 * x)
  )
  dat <- data.frame(x, a, y)

  # Tighten the multinom convergence so the fitted coefficients sit at the same
  # multinomial score root the stacked ee_mlogit block re-solves; the default
  # tolerance leaves a ~1e-6 gap between the fitted probabilities that tilt the
  # oracle and the re-solved probabilities the engine tilts by.
  ps_mod <- nnet::multinom(a ~ x, data = dat, trace = FALSE, reltol = 1e-13)
  ps_mat <- predict(ps_mod, type = "probs")
  focal <- "b"
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_att(ps_mat, dat$a, exposure_type = "categorical", .focal_level = focal)
  )
  outcome_mod <- glm(
    y ~ a + x,
    data = dat,
    family = quasibinomial(),
    weights = wts
  )
  spec <- ipw_spec_categorical(ps_mod, outcome_mod, .focal_level = focal)
  got <- ipw_mestimation(spec)$estimates

  # tilt-standardized plug-in with h = e_focal (the focal-level probability)
  levs <- levels(dat$a)
  h <- ps_mat[, match(focal, levs)]
  mu <- vapply(
    levs,
    function(l) {
      d <- dat
      d$a <- factor(l, levels = levs)
      weighted.mean(predict(outcome_mod, newdata = d, type = "response"), h)
    },
    numeric(1)
  )
  nonref <- levs[-1]
  ref_rd <- vapply(nonref, function(l) mu[[l]] - mu[[1]], numeric(1))
  got_rd <- got$estimate[got$effect == "rd"]
  expect_equal(got_rd, unname(ref_rd), tolerance = 1e-8)
})

# ---- SE oracle: nonparametric bootstrap of the g-computation estimator -------

# Bootstrap the tilt-standardized g-computation estimator for a given estimand,
# refitting both the propensity score and outcome models on each resample. The
# tilt uses the resample's own fitted propensity scores. Runs on the
# heterogeneous-effect fixture so the att target genuinely differs from the ate.
boot_gcomp_diffs <- function(dat, estimand, reps = 200, seed = 7) {
  n <- nrow(dat)
  withr::local_seed(seed)
  replicate(reps, {
    idx <- sample(n, replace = TRUE)
    d <- dat[idx, ]
    ps_b <- glm(z ~ x, data = d, family = binomial())
    e_b <- as.double(predict(ps_b, type = "response"))
    wt_fun <- switch(estimand, ate = wt_ate, att = wt_att)
    w_b <- withr::with_options(list(propensity.quiet = TRUE), wt_fun(ps_b))
    om_b <- lm(y ~ z * x, data = d, weights = w_b)
    plugin_contrasts(om_b, d, "z", estimand = estimand, ps = e_b)[["diff"]]
  })
}

test_that("mestimation ate SE matches a nonparametric bootstrap for an adjusted outcome model", {
  skip_if_not_installed("deli")
  dat <- sim_tilt(seed = 321, n = 1500)
  ps_mod <- glm(z ~ x, data = dat, family = binomial())
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_mod))
  outcome_mod <- lm(y ~ z * x, data = dat, weights = wts)
  spec <- ipw_spec_binary(ps_mod, outcome_mod)
  got <- ipw_mestimation(spec)$estimates
  eng_diff <- got$estimate[got$effect == "diff"]
  eng_se <- got$std.err[got$effect == "diff"]

  boot <- boot_gcomp_diffs(dat, "ate")
  boot_se <- sd(boot)

  # the sandwich SE and the bootstrap SE of the g-computation estimator agree to
  # within a generous band; the point estimate sits inside the bootstrap spread
  expect_lt(abs(eng_se - boot_se) / boot_se, 0.15)
  expect_lt(abs(eng_diff - mean(boot)), 3 * boot_se)
})

test_that("mestimation att SE matches a nonparametric bootstrap for an adjusted outcome model", {
  skip_if_not_installed("deli")
  dat <- sim_tilt(seed = 321, n = 1500)
  ps_mod <- glm(z ~ x, data = dat, family = binomial())
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_att(ps_mod))
  outcome_mod <- lm(y ~ z * x, data = dat, weights = wts)
  spec <- ipw_spec_binary(ps_mod, outcome_mod)
  got <- ipw_mestimation(spec)$estimates
  eng_diff <- got$estimate[got$effect == "diff"]
  eng_se <- got$std.err[got$effect == "diff"]

  boot <- boot_gcomp_diffs(dat, "att")
  boot_se <- sd(boot)

  # the engine must report the tilted estimator: its point estimate sits inside
  # the bootstrap spread of the tilted g-computation (the untilted engine reports
  # the ate-type contrast, many SEs away), and its SE matches the bootstrap
  expect_lt(abs(eng_diff - mean(boot)), 3 * boot_se)
  expect_lt(abs(eng_se - boot_se) / boot_se, 0.15)
})

# ---- outcome-family validation ----------------------------------------------
#
# Only an lm, a gaussian-identity glm, and a binomial or quasibinomial glm are
# supported as outcome models. The spec constructors classify the model with a
# two-way branch, gaussian/identity for an lm or gaussian glm and binomial/<link>
# otherwise, so without a validator ahead of it a poisson, quasipoisson, Gamma,
# or inverse-gaussian model is stacked as a binomial score and a non-identity
# gaussian link is stacked as identity. The rejection tests pin a classed error
# naming the unsupported family or link.

# Data with count and positive-continuous outcomes for the family tests, so an
# outcome model can be fit with poisson, quasipoisson, Gamma, inverse gaussian,
# or a non-identity gaussian link.
family_data <- function(seed = 2024, n = 500) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.6 * x2))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * x1 - 0.3 * x2))
  ycount <- rpois(n, exp(-0.6 + 0.3 * z + 0.2 * x1))
  ypos <- exp(0.2 + 0.3 * z + 0.2 * x1 + 0.3 * rnorm(n))
  data.frame(x1, x2, z, y, ycount, ypos)
}

# Logit ps model and an ATE-weighted outcome model of the requested family. The
# outcome fit is silenced because a binomial fit on fractional weights warns
# about non-integer successes; the family, not the fit, is under test.
family_binary_models <- function(dat, family, response) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, dat$z, exposure_type = "binary", .focal_level = 1)
  fmla <- stats::reformulate("z", response)
  outcome_mod <- if (identical(family, "lm")) {
    lm(fmla, data = dat, weights = wts)
  } else {
    suppressWarnings(glm(fmla, data = dat, family = family, weights = wts))
  }
  list(ps_mod = ps_mod, outcome_mod = outcome_mod)
}

test_that("ipw_spec_binary rejects unsupported outcome families and links", {
  dat <- family_data()
  cases <- list(
    list(label = "poisson", family = poisson(), response = "ycount"),
    list(label = "quasipoisson", family = quasipoisson(), response = "ycount"),
    list(label = "Gamma", family = Gamma(), response = "ypos"),
    list(
      label = "inverse.gaussian",
      family = inverse.gaussian(),
      response = "ypos"
    ),
    list(
      label = "gaussian(log)",
      family = gaussian(link = "log"),
      response = "ypos"
    )
  )
  for (case in cases) {
    mods <- family_binary_models(dat, case$family, case$response)
    expect_error(
      ipw_spec_binary(mods$ps_mod, mods$outcome_mod),
      class = "propensity_ipw_family_error",
      label = case$label
    )
  }
})

test_that("ipw_spec_binary error names the unsupported outcome family", {
  dat <- family_data()
  mods <- family_binary_models(dat, poisson(), "ycount")
  expect_snapshot(
    error = TRUE,
    ipw_spec_binary(mods$ps_mod, mods$outcome_mod)
  )
})

test_that("ipw_spec_binary accepts lm, gaussian-identity, binomial, and quasibinomial outcomes", {
  dat <- family_data()

  m_lm <- family_binary_models(dat, "lm", "y")
  spec_lm <- ipw_spec_binary(m_lm$ps_mod, m_lm$outcome_mod)
  expect_equal(spec_lm$outcome$family, "gaussian")
  expect_equal(spec_lm$outcome$link, "identity")

  m_g <- family_binary_models(dat, gaussian(), "y")
  spec_g <- ipw_spec_binary(m_g$ps_mod, m_g$outcome_mod)
  expect_equal(spec_g$outcome$family, "gaussian")
  expect_equal(spec_g$outcome$link, "identity")

  m_b <- family_binary_models(dat, binomial(), "y")
  spec_b <- ipw_spec_binary(m_b$ps_mod, m_b$outcome_mod)
  expect_equal(spec_b$outcome$family, "binomial")
  expect_equal(spec_b$outcome$link, "logit")

  m_qb <- family_binary_models(dat, quasibinomial(), "y")
  spec_qb <- ipw_spec_binary(m_qb$ps_mod, m_qb$outcome_mod)
  expect_equal(spec_qb$outcome$family, "binomial")
  expect_equal(spec_qb$outcome$link, "logit")
})

# ---- binary propensity score link validation --------------------------------
#
# `ps_link` documents exactly three links for a binary propensity score model:
# logit, probit, and cloglog, which is also the set the linearization path
# arg-matches in derive_link() and derive_score_factor(). The M-estimation path
# restricts the resolved link to the same three at the entry of
# ipw_spec_binary(). The guard belongs there because the resolved link is
# otherwise only ever consumed by ipw_inv_link(), which serves outcome-link
# reconstruction as well and so accepts identity and log too: an unsupported
# propensity model would be estimated without comment or abort from an internal
# frame with a message that names neither the argument nor the supported set.
#
# Controls for the three supported links on the M-estimation path are covered by
# "ipw_mestimation std.err matches linearization within 3 percent", which builds
# a spec with ipw_spec_binary() and solves it for each of logit, probit, and
# cloglog.

test_that("mestimation rejects a propensity score link outside the supported set", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  # Both entry points are asserted. The validation belongs at spec construction,
  # so ipw_spec_binary() must reject the link on its own rather than relying on
  # ipw() to reach it first.
  for (link in c("cauchit", "log", "identity")) {
    mods <- fit_binary_models(dat, "ate", ps_link = link)
    expect_error(
      ipw(mods$ps_mod, mods$outcome_mod),
      class = "propensity_ipw_link_error",
      label = link
    )
    expect_error(
      ipw_spec_binary(mods$ps_mod, mods$outcome_mod),
      class = "propensity_ipw_link_error",
      label = link
    )
  }
})

test_that("the propensity score link error names the link and the supported set", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate", ps_link = "cauchit")

  cnd <- rlang::catch_cnd(ipw(mods$ps_mod, mods$outcome_mod))
  msg <- conditionMessage(cnd)

  # The message must name the offending link and every supported link, so a user
  # reading it knows both what was rejected and what to fit instead. Asserted
  # outside the snapshot because snapshots are skipped under --as-cran.
  expect_match(msg, "cauchit", fixed = TRUE)
  expect_match(msg, "logit", fixed = TRUE)
  expect_match(msg, "probit", fixed = TRUE)
  expect_match(msg, "cloglog", fixed = TRUE)

  # The error must be attributed to ipw(), the user-facing entry point, not to
  # the internal frame that raises it.
  expect_identical(as.character(cnd$call[[1]]), "ipw")

  expect_snapshot(error = TRUE, ipw(mods$ps_mod, mods$outcome_mod))
})

test_that("the ps_link argument is validated against the supported set", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate", ps_link = "logit")

  # A supported override is still accepted.
  expect_s3_class(ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit"), "ipw")

  # The argument overrides the family link, so what must be validated is the
  # resolved link rather than the link the propensity model was fit with. An
  # unsupported override otherwise reaches the weight-consistency preflight,
  # which reports a weights mismatch and directs the user to refit the outcome
  # model, naming the wrong problem.
  for (link in c("log", "identity", "cauchit")) {
    expect_error(
      ipw(mods$ps_mod, mods$outcome_mod, ps_link = link),
      class = "propensity_ipw_link_error",
      label = link
    )
  }
})

# ---- factor and logical outcome responses -----------------------------------
#
# The spec constructors store the outcome through ipw_outcome_numeric(), where
# the raw response comes from model.frame or a .data column. A plain
# as.double() on a binomial glm's factor response (first level failure, others
# success) yields level codes 1/2 instead of 0/1, leaving the stacked binomial
# score with no root and failing the solve. A logical response already doubles
# to 0/1. These tests pin that a factor response matches the numeric-0/1 fit.
# quasibinomial is used throughout: its $y matches binomial for a factor response
# and it is the weighted-outcome convention. quasibinomial does not warn about
# non-integer successes on fractional weights (unlike binomial), so the fixture
# fits are wrapped in suppressWarnings only defensively.

response_binary_data <- function() {
  dat <- sim_binary()
  dat$yf <- factor(ifelse(dat$y == 1, "yes", "no"), levels = c("no", "yes"))
  dat$yl <- dat$y == 1
  dat
}

response_binary_setup <- function(dat) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- wt_ate(
    predict(ps_mod, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  list(ps_mod = ps_mod, wts = wts)
}

response_outcome <- function(response, dat, wts) {
  suppressWarnings(glm(
    stats::reformulate("z", response),
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  ))
}

test_that("mestimation matches a factor outcome response to the numeric fit", {
  skip_if_not_installed("deli")
  dat <- response_binary_data()
  s <- response_binary_setup(dat)
  num <- ipw(s$ps_mod, response_outcome("y", dat, s$wts))$estimates
  fac <- ipw(s$ps_mod, response_outcome("yf", dat, s$wts))$estimates
  expect_equal(fac$estimate, num$estimate, tolerance = 1e-8)
  expect_equal(fac$std.err, num$std.err, tolerance = 1e-8)
})

test_that("mestimation matches a logical outcome response to the numeric fit", {
  skip_if_not_installed("deli")
  dat <- response_binary_data()
  s <- response_binary_setup(dat)
  num <- ipw(s$ps_mod, response_outcome("y", dat, s$wts))$estimates
  lgl <- ipw(s$ps_mod, response_outcome("yl", dat, s$wts))$estimates
  expect_equal(lgl$estimate, num$estimate, tolerance = 1e-8)
  expect_equal(lgl$std.err, num$std.err, tolerance = 1e-8)
})

test_that("mestimation matches a factor outcome supplied through .data", {
  skip_if_not_installed("deli")
  dat <- response_binary_data()
  s <- response_binary_setup(dat)
  num <- ipw(s$ps_mod, response_outcome("y", dat, s$wts), .data = dat)$estimates
  fac <- ipw(
    s$ps_mod,
    response_outcome("yf", dat, s$wts),
    .data = dat
  )$estimates
  expect_equal(fac$estimate, num$estimate, tolerance = 1e-8)
  expect_equal(fac$std.err, num$std.err, tolerance = 1e-8)
})

test_that("mestimation matches a three-level factor outcome to the numeric fit", {
  skip_if_not_installed("deli")
  dat <- response_binary_data()
  # A three-level factor whose first level is the failure: glm maps every
  # non-first level to success, so the response is 0/1, but as.double(f) is
  # 1/2/3 and as.double(f) - 1 is 0/1/2. Only reading the converted response
  # (level != first) gives the numeric fit; this pins the .data-column rule.
  dat$yf3 <- factor(
    ifelse(dat$y == 0, "none", ifelse(dat$x1 > 0, "severe", "mild")),
    levels = c("none", "mild", "severe")
  )
  dat$ynum <- as.double(dat$yf3 != "none")
  s <- response_binary_setup(dat)
  num <- ipw(
    s$ps_mod,
    response_outcome("ynum", dat, s$wts),
    .data = dat
  )$estimates
  fac <- ipw(
    s$ps_mod,
    response_outcome("yf3", dat, s$wts),
    .data = dat
  )$estimates
  expect_equal(fac$estimate, num$estimate, tolerance = 1e-8)
  expect_equal(fac$std.err, num$std.err, tolerance = 1e-8)
})

# ---- binary ps design extraction with .data ---------------------------------
#
# ipw_spec_binary takes its propensity design from the shared
# ipw_extract_ps_design() helper, in parity with the categorical and continuous
# paths. So a binary ps model fit with model = FALSE whose fitting data is gone
# rebuilds its design from .data, and without .data raises the informative
# propensity_ipw_data_error rather than a raw "object not found". These tests
# pin the guarded, .data-aware behavior.

ps_design_data <- function(seed = 2024, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  cov <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * (cov == "b")))
  y <- rbinom(n, 1, plogis(-0.4 + 1.0 * z + 0.5 * x1))
  data.frame(x1, cov, z, y)
}

test_that("mestimation reconstructs the binary ps design from .data when the fit frame is gone", {
  skip_if_not_installed("deli")
  dat <- ps_design_data()
  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)
  # model = FALSE drops the stored model frame; removing the fitting data leaves
  # model.matrix(ps_mod) unable to reconstruct it, so .data must rebuild it.
  ps_gone <- glm(z ~ x1, data = dat, family = binomial(), model = FALSE)
  wts <- wt_ate(
    predict(ps_gone, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  out <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = ctrl
  )
  # reconstructable reference: an identical fit whose frame is retained
  ps_ref <- glm(z ~ x1, data = dat, family = binomial())
  dat_copy <- dat
  rm(dat)

  ref <- ipw(ps_ref, out, se_method = "mestimation")$estimates
  res <- ipw(
    ps_gone,
    out,
    .data = dat_copy,
    se_method = "mestimation"
  )$estimates
  expect_equal(res, ref, tolerance = 1e-8)
})

test_that("mestimation errors with the supply-.data hint when the binary ps fit frame is gone", {
  dat <- ps_design_data()
  ps_gone <- glm(z ~ x1, data = dat, family = binomial(), model = FALSE)
  wts <- wt_ate(
    predict(ps_gone, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  out <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
  rm(dat)

  # Without .data the design cannot be rebuilt, so the guarded error stands in
  # for the raw "object not found".
  expect_error(
    ipw(ps_gone, out, se_method = "mestimation"),
    class = "propensity_ipw_data_error"
  )
})

test_that("mestimation binary estimates are unchanged when redundant .data is supplied", {
  skip_if_not_installed("deli")
  # A normal, reconstructable ps fit with a factor covariate: supplying .data
  # must not change the estimates. The rebuild reproduces the same design
  # (xlevels and all), so this guards against drift when .data is redundant and
  # exercises the factor-covariate design rebuild.
  dat <- ps_design_data(seed = 11)
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- wt_ate(
    predict(ps_mod, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  out <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  no_data <- ipw(ps_mod, out, se_method = "mestimation")$estimates
  with_data <- ipw(
    ps_mod,
    out,
    .data = dat,
    se_method = "mestimation"
  )$estimates
  expect_equal(with_data, no_data, tolerance = 1e-8)
})

test_that("mestimation binary estimates are unchanged when .data re-levels a ps factor covariate", {
  skip_if_not_installed("deli")
  # .data with the same factor-covariate values but a reordered levels attribute
  # (a user re-leveling for presentation after fitting) stays row-aligned. The
  # ps design rebuild honors ps_mod$xlevels so the design matches the fit.
  dat <- ps_design_data(seed = 11)
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- wt_ate(
    predict(ps_mod, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  out <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  dat_relevel <- dat
  dat_relevel$cov <- factor(as.character(dat$cov), levels = c("c", "b", "a"))

  no_data <- ipw(ps_mod, out, se_method = "mestimation")$estimates
  releveled <- ipw(
    ps_mod,
    out,
    .data = dat_relevel,
    se_method = "mestimation"
  )$estimates
  expect_equal(releveled, no_data, tolerance = 1e-8)
})

# ---- a wrong-sized .data is a data problem ----------------------------------
#
# The linearization path already rejects a `.data` whose row count disagrees
# with the fitted models, naming the two counts and pointing at `.data`. The
# M-estimation paths size everything to `.data` instead, so the disagreement
# surfaces further downstream in the weight-consistency preflight, whose message
# is about weights and focal levels. That message sends someone to look at how
# they built their weights when what they actually passed was the wrong data
# frame. The row count is the more specific diagnosis and belongs first.
#
# The count check owns count mismatches only. A `.data` of the right size but
# the wrong content is not detectable from a row count and stays with the
# preflight, which recomputes the propensity scores and notices; that division
# is pinned below.
#
# The linearization behavior these mirror is pinned in test-ipw-se-method.R by
# "linearization rejects a .data whose row count disagrees with the fits".

test_that("mestimation rejects a .data with too few rows", {
  skip_if_not_installed("deli")
  dat <- ps_design_data()
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      predict(ps_mod, type = "response"),
      dat$z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  out <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(ps_mod, out, .data = dat[-1, ], se_method = "mestimation"),
    class = "propensity_ipw_data_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "one row per observation", fixed = TRUE)
})

test_that("mestimation rejects a .data with too many rows", {
  skip_if_not_installed("deli")
  # A longer .data is the same mistake as a shorter one and today produces the
  # same weights-mismatch message. The check is on inequality, not on being
  # short, so both directions must report the row counts.
  dat <- ps_design_data()
  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      predict(ps_mod, type = "response"),
      dat$z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  out <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(
      ps_mod,
      out,
      .data = rbind(dat, dat[1, ]),
      se_method = "mestimation"
    ),
    class = "propensity_ipw_data_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "one row per observation", fixed = TRUE)
})

test_that("a right-sized .data with the wrong content stays a weights problem", {
  skip_if_not_installed("deli")
  # The boundary between the two checks. Different data of the same size passes
  # any row count comparison, and the propensity scores recomputed from it no
  # longer reproduce the supplied weights, which is what the preflight exists to
  # catch. The new count check must not take this case over.
  dat <- ps_design_data(seed = 11)
  other <- ps_design_data(seed = 99)
  expect_identical(nrow(other), nrow(dat))

  ps_mod <- glm(z ~ x1 + cov, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      predict(ps_mod, type = "response"),
      dat$z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  out <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  expect_error(
    ipw(ps_mod, out, .data = other, se_method = "mestimation"),
    class = "propensity_ipw_weights_mismatch_error"
  )
})

# ---- design rebuilds must honor the fitted models' contrasts ----------------
#
# Every design the engine rebuilds is multiplied by coefficients from the model
# that produced it, and the multiply is positional. A model fit under
# non-default contrasts therefore has to be rebuilt under those same contrasts
# or the two disagree column for column. A contrast change is only a
# reparameterization: the fit, its fitted values, and every marginal quantity
# are unchanged, so the estimates must match the default-coded fit exactly.
# These pin against a default-coded fit computed in the test rather than against
# stored numbers.
#
# The exposure column and the covariate columns fare differently. model.frame()
# keeps the contrasts attribute on the columns it stores, so a rebuild that only
# replaces the exposure leaves a covariate's coding intact; the binary
# counterfactual rebuild is faithful for that reason and is pinned below as a
# regression. The .data rebuild of the ps design is not: model.matrix() rebuilds
# the factor against xlev, which drops the attribute.

contrast_design_data <- function(seed = 2024, n = 800) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  g <- factor(
    sample(c("p", "q", "r"), n, replace = TRUE),
    levels = c("p", "q", "r")
  )
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.6 * (g == "r")))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * x1 - 0.3 * (g == "q")))
  data.frame(x1, g, z, y)
}

# Propensity score model adjusted for the factor covariate, ate weights, and an
# exposure-only outcome model, so the covariate reaches the ps design rebuild
# and nothing else.
fit_contrast_models <- function(dat) {
  ps_mod <- glm(z ~ x1 + g, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_mod)
  )
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

test_that("mestimation accepts a sum-coded ps covariate rebuilt from .data", {
  skip_if_not_installed("deli")
  # The ps design rebuild drops the sum coding and rebuilds under treatment
  # coding, so the recomputed propensity scores no longer match the fit. The
  # weights supplied to the outcome model are correct; it is the rebuild that is
  # wrong, and today the mismatch is reported as though the user's weights were
  # at fault.
  dat <- contrast_design_data()
  dat_sum <- dat
  contrasts(dat_sum$g) <- contr.sum(3)

  default <- fit_contrast_models(dat)
  sum_coded <- fit_contrast_models(dat_sum)

  # a reparameterization: same fit, same weights, different coefficient basis
  expect_equal(
    names(coef(sum_coded$ps_mod)),
    c("(Intercept)", "x1", "g1", "g2")
  )
  expect_equal(
    as.double(sum_coded$wts),
    as.double(default$wts),
    tolerance = 1e-10
  )

  reference <- ipw(
    default$ps_mod,
    default$outcome_mod,
    .data = dat,
    se_method = "mestimation"
  )$estimates
  got <- ipw(
    sum_coded$ps_mod,
    sum_coded$outcome_mod,
    .data = dat_sum,
    se_method = "mestimation"
  )$estimates

  # As elsewhere in this section, the point estimates agree to machine precision
  # and the sandwich standard errors carry the numerical derivative's noise at a
  # different coefficient basis.
  expect_equal(got$estimate, reference$estimate, tolerance = 1e-12)
  expect_equal(got, reference, tolerance = 1e-6)
})

test_that("the .data ps design rebuild does not warn on a sum-coded covariate", {
  skip_if_not_installed("deli")
  # Rebuilding the factor against xlev drops its contrasts attribute, and R
  # signals that with a bare "contrasts dropped from factor g" warning. The
  # linearization estimates survive it because the two codings span the same
  # column space and the correction is a projection, but the warning is a real
  # signal that the rebuilt design is not the fitted one.
  dat <- contrast_design_data()
  dat_sum <- dat
  contrasts(dat_sum$g) <- contr.sum(3)
  mods <- fit_contrast_models(dat_sum)

  expect_no_warning(
    ipw(
      mods$ps_mod,
      mods$outcome_mod,
      .data = dat_sum,
      se_method = "linearization"
    )
  )
})

test_that("mestimation binary counterfactual designs keep a sum-coded outcome covariate", {
  skip_if_not_installed("deli")
  # Regression pin, correct today. The counterfactual rebuild replaces only the
  # exposure column, and model.frame() hands back the covariate with its
  # contrasts attribute intact, so the rebuilt design already matches the fit.
  # Threading contrasts through the rebuilds must not disturb this.
  dat <- contrast_design_data()
  dat_sum <- dat
  contrasts(dat_sum$g) <- contr.sum(3)
  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)

  fit_adjusted <- function(d) {
    ps_mod <- glm(z ~ x1 + g, data = d, family = binomial())
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_ate(ps_mod)
    )
    outcome_mod <- glm(
      y ~ z + g,
      data = d,
      family = quasibinomial(),
      weights = wts,
      control = ctrl
    )
    ipw(ps_mod, outcome_mod, se_method = "mestimation")$estimates
  }

  got <- fit_adjusted(dat_sum)
  reference <- fit_adjusted(dat)

  # The point estimates are what the counterfactual rebuild determines, and they
  # agree to machine precision. The sandwich standard errors carry a relative
  # difference around 1e-7, which is the numerical derivative evaluated at a
  # different coefficient basis rather than a coding mismatch; a real one moves
  # the estimates themselves by whole percentage points.
  expect_equal(got$estimate, reference$estimate, tolerance = 1e-12)
  expect_equal(got, reference, tolerance = 1e-6)
})

test_that("mestimation binary counterfactual designs keep a sum-coded factor exposure", {
  skip_if_not_installed("deli")
  # A two-level factor exposure carries contrasts of its own. Setting the
  # exposure column to one level subsets the factor, which drops the attribute,
  # so this is the binary counterpart of the categorical rebuild rather than of
  # the covariate case above: without the fit's coding the design is treatment
  # coded against sum-coded coefficients and the contrast changes sign.
  dat <- contrast_design_data()
  dat$zf <- factor(
    ifelse(dat$z == 1, "trt", "ctl"),
    levels = c("ctl", "trt")
  )
  dat_sum <- dat
  contrasts(dat_sum$zf) <- contr.sum(2)
  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)

  fit_factor_exposure <- function(d) {
    ps_mod <- glm(zf ~ x1, data = d, family = binomial())
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_ate(
        predict(ps_mod, type = "response"),
        d$zf,
        exposure_type = "binary",
        .focal_level = "trt"
      )
    )
    outcome_mod <- glm(
      y ~ zf + x1,
      data = d,
      family = quasibinomial(),
      weights = wts,
      control = ctrl
    )
    ipw(ps_mod, outcome_mod, se_method = "mestimation")$estimates
  }

  got <- fit_factor_exposure(dat_sum)
  reference <- fit_factor_exposure(dat)

  expect_equal(got$estimate, reference$estimate, tolerance = 1e-12)
  expect_equal(got, reference, tolerance = 1e-6)
})

test_that("mestimation is unaffected by a contrasts option changed after fitting", {
  skip_if_not_installed("deli")
  # The rebuilds used to consult the ambient contrasts option rather than the
  # fit, so a fit made under a non-default option and analyzed after it was
  # restored silently disagreed with itself. Reading the coding off the fitted
  # models makes the call independent of the option in force when it runs.
  dat <- contrast_design_data()

  under_option <- withr::with_options(
    list(contrasts = c("contr.sum", "contr.poly")),
    fit_contrast_models(dat)
  )
  expect_equal(
    names(coef(under_option$ps_mod)),
    c("(Intercept)", "x1", "g1", "g2")
  )
  expect_equal(getOption("contrasts")[["unordered"]], "contr.treatment")

  default <- fit_contrast_models(dat)

  expect_equal(
    ipw(
      under_option$ps_mod,
      under_option$outcome_mod,
      se_method = "mestimation"
    )$estimates$estimate,
    ipw(
      default$ps_mod,
      default$outcome_mod,
      se_method = "mestimation"
    )$estimates$estimate,
    tolerance = 1e-12
  )
})

# ---- separation in the propensity score model -------------------------------
#
# A propensity model whose covariates perfectly predict the exposure has no
# finite maximum likelihood estimate. `glm` stops when its convergence criterion
# is met rather than at an optimum, leaving a linear predictor in the tens of
# thousands, and reports that fitted probabilities of numerically zero or one
# occurred.
#
# `glm` itself never hands back an exact zero or one: `binomial()$linkinv`
# clamps its argument, so the largest fitted value it can produce is the double
# just below one. The M-estimation path does not go through that clamp. It
# rebuilds the propensity scores as `plogis(X %*% theta)`, and `plogis`
# saturates to exactly one above about 36.74, so the weights recomputed for the
# preflight contain 0/0 where the supplied weights were finite. What the user
# sees is the generic weight-consistency error, which asks them to check how
# they built their weights and to consider the focal level. Neither is the
# cause, and neither leads anywhere: the weights are the ones this propensity
# model implies, and the model is the problem.
#
# Only the `ate` weights go non-finite. The other estimands carry a tilt that
# vanishes at the saturation point, which cancels the singularity, so their
# recomputed weights stay finite and the preflight is satisfied. They break
# further in instead; that is recorded separately and is not what these pin.

separation_data <- function(seed = 2024, n = 400) {
  set.seed(seed)
  x1 <- rnorm(n)
  # x1 predicts z with no error, so the fit has no finite optimum
  z <- as.integer(x1 > 0)
  y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.4 * x1))
  data.frame(x1, z, y)
}

# The propensity fit warns that fitted probabilities of zero or one occurred,
# which is the fixture working as intended rather than anything under test. The
# iteration limit is raised so the linear predictor reaches the range where
# plogis saturates, which is the state these tests are about.
separation_models <- function(dat) {
  ps_mod <- suppressWarnings(glm(
    z ~ x1,
    data = dat,
    family = binomial(),
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

test_that("mestimation reports separation in the propensity model as its own error", {
  skip_if_not_installed("deli")
  dat <- separation_data()
  mods <- separation_models(dat)

  # The fixture is separated, and the weights the user holds are finite: the
  # mismatch is manufactured by the rebuild, not by anything they did.
  expect_gt(max(abs(predict(mods$ps_mod))), 100)
  expect_true(all(is.finite(as.double(mods$wts))))

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation"),
    class = "propensity_ipw_separation_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "separation", ignore.case = TRUE)
})

test_that("mestimation reports separation when .data is supplied", {
  skip_if_not_installed("deli")
  # The propensity scores are rebuilt from the same coefficients either way, so
  # supplying .data does not change what saturates.
  dat <- separation_data()
  mods <- separation_models(dat)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat, se_method = "mestimation"),
    class = "propensity_ipw_separation_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "separation", ignore.case = TRUE)
})

# ---- the outcome model must contain the exposure -----------------------------
#
# The counterfactual designs are built by setting the exposure to each level in
# the outcome model's data, so an outcome model fit on covariates alone yields
# identical X1 and X0 designs: without a guard every contrast would collapse to
# zero with a degenerate standard error rather than erroring.
# check_ipw_outcome_exposure() reads the model terms, so it fires on both
# extraction routes. The .data route was the degenerate one; the no-.data route
# reached check_exposure(), whose "Please specify .data" hint funnelled the user
# straight into it. Both routes are pinned below, so moving the guard into the
# .data branch alone would fail. A transformed exposure such as I(z^2) still
# names the exposure and is not this error; that case belongs to
# check_exposure() and is pinned in test-ipw.R.

test_that("mestimation rejects a binary outcome model that omits the exposure when .data is supplied", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  no_exposure <- fit_outcome_without_exposure(dat, mods$wts)

  expect_error(
    ipw(mods$ps_mod, no_exposure, .data = dat),
    class = "propensity_ipw_exposure_error"
  )
})

test_that("mestimation rejects a binary outcome model that omits the exposure without .data", {
  skip_if_not_installed("deli")
  # The guard reads the model terms, so it fires before the frame-based
  # check_exposure() on this route too. Without it the call reached
  # check_exposure(), whose "Please specify .data" hint sent the user to the
  # .data route and its silently degenerate fit.
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  no_exposure <- fit_outcome_without_exposure(dat, mods$wts)

  expect_error(
    ipw(mods$ps_mod, no_exposure),
    class = "propensity_ipw_exposure_error"
  )
})

test_that("the missing-exposure error names the exposure and directs the user to refit", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  no_exposure <- fit_outcome_without_exposure(dat, mods$wts)

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, no_exposure, .data = dat)
  )
})

test_that("mestimation accepts a binary outcome model containing the exposure when .data is supplied", {
  skip_if_not_installed("deli")
  # The guard must not fire on the ordinary case: the same call with the
  # exposure on the right-hand side runs through to an ipw object.
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")

  res <- ipw(mods$ps_mod, mods$outcome_mod, .data = dat)

  expect_s3_class(res, "ipw")
  expect_equal(res$estimand, "ate")
})

# ---- degenerate counterfactual designs --------------------------------------
#
# The counterfactual designs come from delete.response(terms(outcome_mod)) with
# the exposure set to each level in turn. A numeric exposure fit without an
# intercept, y ~ z - 1, leaves the design at z = 0 identically zero, so mu0 is
# inv_link(0) for every unit whatever the data say: 0.5 for a logit or probit
# outcome link and 0 for a linear model. Nothing about the fit signals it. On
# this fixture
# the risk difference reads 0.1306 against the correct 0.2347, with a standard
# error of 0.0255 against 0.0340, and the linear model reports a difference of
# 2.179 against the correct 0.788.
#
# The condition is a design that is identically zero, not a missing intercept.
# A saturated factor coding, y ~ 0 + zf, has dummy columns that sum to one at
# every level and reproduces the with-intercept fit; y ~ z + x1 - 1 still
# estimates mu0 from the covariate and is honest g-computation on the model as
# specified. Both must keep working.

# Two-level factor recoding of the numeric exposure, with the propensity score
# model and the ate weights that coding implies. A no-intercept outcome model
# on this coding is a saturated reparameterization, not a degenerate one.
fit_binary_factor_models <- function(dat, intercept = TRUE) {
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

# The weighted outcome model of a numeric binary exposure, refit without the
# intercept. Any `outcome_family` other than "binomial" fits the linear
# counterpart on the continuous outcome.
fit_outcome_no_intercept <- function(dat, wts, outcome_family = "binomial") {
  if (identical(outcome_family, "binomial")) {
    glm(
      y ~ z - 1,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  } else {
    lm(yc ~ z - 1, data = dat, weights = wts)
  }
}

test_that("mestimation rejects a binary outcome model whose unexposed design is identically zero", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  no_intercept <- fit_outcome_no_intercept(dat, mods$wts)

  err <- expect_error(
    ipw(mods$ps_mod, no_intercept),
    class = "propensity_error",
    regexp = "identically zero"
  )

  # The message must name the exposure level whose design is degenerate, so a
  # user can tell which marginal mean was pinned. Whitespace is normalized
  # because cli wraps the bullet.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "`z` to \"0\"", fixed = TRUE)
})

test_that("mestimation rejects a degenerate counterfactual design when .data is supplied", {
  skip_if_not_installed("deli")
  # Supplying .data rebuilds the counterfactual designs from the supplied
  # columns rather than from the model frame. The guard must fire on that route
  # too.
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  no_intercept <- fit_outcome_no_intercept(dat, mods$wts)

  expect_error(
    ipw(mods$ps_mod, no_intercept, .data = dat),
    class = "propensity_error",
    regexp = "identically zero"
  )
})

test_that("mestimation rejects a degenerate counterfactual design in a linear outcome model", {
  skip_if_not_installed("deli")
  # The guard reads the design, not the family, so the linear model whose mu0 is
  # pinned at 0 rather than 0.5 is rejected on the same grounds.
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate", outcome_family = "gaussian")
  no_intercept <- fit_outcome_no_intercept(dat, mods$wts, "gaussian")

  expect_error(
    ipw(mods$ps_mod, no_intercept),
    class = "propensity_error",
    regexp = "identically zero"
  )
})

test_that("the degenerate-design error names the pinned level and the remedy", {
  skip_if_not_installed("deli")
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  no_intercept <- fit_outcome_no_intercept(dat, mods$wts)

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, no_intercept)
  )
})

test_that("mestimation accepts a saturated no-intercept factor outcome model", {
  skip_if_not_installed("deli")
  # y ~ 0 + zf is a reparameterization of y ~ zf, not a model without a
  # baseline: its counterfactual designs are the level indicators, which are
  # never zero. It must give the same answer as the with-intercept fit.
  dat <- sim_binary()
  with_intercept <- fit_binary_factor_models(dat, intercept = TRUE)
  saturated <- fit_binary_factor_models(dat, intercept = FALSE)

  res_sat <- ipw(saturated$ps_mod, saturated$outcome_mod)
  res_int <- ipw(with_intercept$ps_mod, with_intercept$outcome_mod)

  expect_s3_class(res_sat, "ipw")
  expect_equal(res_sat$estimates$effect, res_int$estimates$effect)
  expect_equal(
    res_sat$estimates$estimate,
    res_int$estimates$estimate,
    tolerance = 1e-6
  )
  expect_equal(
    res_sat$estimates$std.err,
    res_int$estimates$std.err,
    tolerance = 1e-6
  )
})

test_that("mestimation accepts a no-intercept outcome model adjusted for a covariate", {
  skip_if_not_installed("deli")
  # y ~ z + x1 - 1 has no intercept, but the covariate still varies at z = 0, so
  # mu0 is estimated from the data rather than pinned. It is g-computation on
  # the model as specified and must run.
  dat <- sim_binary()
  mods <- fit_binary_models(dat, "ate")
  adjusted <- glm(
    y ~ z + x1 - 1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  res <- ipw(mods$ps_mod, adjusted)

  expect_s3_class(res, "ipw")
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(res$estimates$std.err > 0))
})
