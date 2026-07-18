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

# Estimand tilt h(e) for the standardized g-computation plug-in. ate is the flat
# tilt (h = 1); a tilted estimand weights each unit's counterfactual predictions
# by h of its propensity score. Mirrors the tilt the M-estimation marginal-mean
# rows apply once they standardize to the target population.
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
    weights = wts
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
# The marginal-mean rows of the stacked system currently take an unweighted mean
# of counterfactual predictions over the whole sample for every estimand, which
# is correct only for ate. For a tilted estimand the marginal means must be
# standardized to the target population, mu_a = sum_i h(e_i) m_a(x_i) / sum_i
# h(e_i). With a covariate-adjusted outcome model and a heterogeneous effect the
# current engine reports the ate-type contrast for every estimand. The tests
# below pin the tilt-standardized behavior; those that cannot pass until the
# standardization is in place carry the pending skip. A saturated outcome model
# (y ~ z) escapes the defect because predictions are constant within arm, so the
# tilt is a no-op on the contrast; that case is a live regression anchor.

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
  # tilt is a no-op on the contrast and the engine already matches the tilt-
  # standardized plug-in. This anchors the exposure-only case across the fix.
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
# The spec constructors classify the outcome model with a two-way branch:
# gaussian/identity for an lm or gaussian glm, binomial/<link> otherwise. That
# silently misclassifies poisson, quasipoisson, Gamma, and inverse-gaussian
# outcome models as a binomial score, and ignores the link of a gaussian glm
# (gaussian(link = "log") is stacked as identity). Only lm, gaussian-identity
# glm, and binomial or quasibinomial glm are supported. The rejection tests pin
# a classed error naming the unsupported family or link; they cannot pass until
# the shared validator exists.

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

# ---- factor and logical outcome responses -----------------------------------
#
# The spec constructors store as.double(outcome), where outcome is the raw
# response from model.frame or a .data column. A binomial glm with a factor
# response (first level failure, others success) then yields level codes 1/2
# instead of 0/1, so the stacked binomial score has no root and the solve fails.
# A logical response already doubles to 0/1. These tests pin that a factor
# response matches the numeric-0/1 fit; the fix will read the converted response.
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
  skip("pending factor outcome conversion")
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
  skip("pending factor outcome conversion")
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
  skip("pending factor outcome conversion")
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
