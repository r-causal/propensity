# Tests for ipw() with a multinomial (nnet::multinom) propensity score model and
# a categorical exposure. These exercise ipw() end to end on a fitted multinom
# propensity score model and a weighted outcome model, the way a user would call
# it, pinning estimand detection, the estimates-table contract, point-estimate
# parity with g-computation, a PSweight oracle, focal-level handling, and the
# print and as.data.frame surfaces.

# ---- data simulator ---------------------------------------------------------

sim_categorical <- function(seed = 2024, n = 700) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  eta_b <- -0.2 + 0.5 * x1 + 0.3 * x2
  eta_c <- 0.1 - 0.4 * x1 + 0.6 * x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  pa <- 1 / denom
  pb <- exp(eta_b) / denom
  u <- runif(n)
  lab <- ifelse(u < pa, "a", ifelse(u < pa + pb, "b", "c"))
  a <- factor(lab, levels = c("a", "b", "c"))
  y <- rbinom(
    n,
    1,
    plogis(-0.3 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.5 * x1)
  )
  yc <- 0.5 + 0.4 * (a == "b") + 0.9 * (a == "c") + 0.6 * x1 + rnorm(n)
  data.frame(x1, x2, a, y, yc)
}

# ---- model fitting ----------------------------------------------------------

# Tighten the multinomial fit so the seeded init in the M-estimator sits at the
# score root well below the 1e-8 point-estimate comparison tolerance, matching
# the reltol/maxit pattern used in test-ipw-psi.R.
fit_ps_multinom <- function(dat) {
  nnet::multinom(
    a ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
}

# Propensity score matrix in factor-level order, column names set to the levels.
# The levels come from the fitted model, which is the authority on the order its
# own columns arrive in and is what the package reads everywhere else.
ps_matrix_named <- function(ps_mod) {
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- ps_mod$lev
  ps
}

# Categorical weights for a given estimand, always computed silently.
categorical_weights <- function(
  ps_named,
  a,
  estimand,
  focal_level = NULL,
  stabilize = FALSE
) {
  wt_fun <- switch(
    estimand,
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )
  args <- list(ps_named, a, exposure_type = "categorical")
  if (estimand %in% c("att", "atu")) {
    args$.focal_level <- focal_level
  }
  if (estimand == "ate") {
    args$stabilize <- stabilize
  }
  withr::with_options(
    list(propensity.quiet = TRUE),
    do.call(wt_fun, args)
  )
}

# Fit a weighted outcome model. `covariates = TRUE` adds x1 (the covariate-
# adjusted g-computation setup); `covariates = FALSE` fits an exposure-only
# saturated model whose g-computation marginal means equal the Hajek weighted
# group means, the estimator PSweight reports without an outcome model. The
# binomial outcome uses a tightened IRLS tolerance so the fitted coefficients
# sit at the weighted MLE to well below the point-estimate tolerance.
fit_outcome <- function(
  dat,
  wts,
  outcome_family = "binomial",
  covariates = TRUE
) {
  outcome_var <- if (outcome_family == "binomial") "y" else "yc"
  rhs <- if (covariates) c("a", "x1") else "a"
  fmla <- stats::reformulate(rhs, response = outcome_var)
  if (outcome_family == "binomial") {
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
}

# Build the propensity score model, the categorical weights, and the weighted
# outcome model in one step, returning all three plus the level vector.
fit_categorical_models <- function(
  dat,
  estimand,
  focal_level = NULL,
  stabilize = FALSE,
  outcome_family = "binomial",
  covariates = TRUE,
  strip_weights = FALSE
) {
  ps_mod <- fit_ps_multinom(dat)
  ps_named <- ps_matrix_named(ps_mod)
  wts <- categorical_weights(
    ps_named,
    dat$a,
    estimand,
    focal_level = focal_level,
    stabilize = stabilize
  )
  model_wts <- if (strip_weights) as.double(wts) else wts
  outcome_mod <- fit_outcome(
    dat,
    model_wts,
    outcome_family = outcome_family,
    covariates = covariates
  )
  list(
    ps_mod = ps_mod,
    outcome_mod = outcome_mod,
    wts = wts,
    lev = levels(dat$a)
  )
}

# ---- plug-in reference ------------------------------------------------------

# Categorical estimand tilt h(e) for the standardized g-computation plug-in. ate
# is the flat tilt; the tilted estimands weight each unit's counterfactual
# predictions by h of the propensity score matrix (columns in level order).
# Mirrors ps_tilt_categorical() in R/ps_tilt.R.
categorical_plugin_tilt <- function(ps, estimand, lev, focal, n) {
  if (estimand == "ate") {
    return(rep(1, n))
  }
  focal_col <- if (!is.null(focal)) match(focal, lev) else NULL
  switch(
    estimand,
    att = ps[, focal_col],
    atu = 1 - ps[, focal_col],
    ato = 1 / rowSums(1 / ps),
    atm = do.call(pmin, lapply(seq_len(ncol(ps)), function(j) ps[, j])),
    entropy = -rowSums(ps * log(ps))
  )
}

# Direct g-computation on the weighted outcome model: predict the counterfactual
# marginal mean for each level, standardized to the estimand's tilted population
# by h(ps) (h = 1 for ate), then form per-non-reference-level contrasts against
# the first (reference) level. Returns a data frame keyed by effect and
# comparison so it can be matched to the ipw() estimates table irrespective of
# row order.
plugin_categorical <- function(
  outcome_mod,
  dat,
  lev,
  outcome_family,
  estimand = "ate",
  ps = NULL,
  focal = NULL
) {
  h <- categorical_plugin_tilt(ps, estimand, lev, focal, nrow(dat))
  mu <- vapply(
    lev,
    function(l) {
      d <- dat
      d$a <- factor(l, levels = lev)
      weighted.mean(predict(outcome_mod, newdata = d, type = "response"), h)
    },
    numeric(1)
  )
  ref <- mu[[1]]
  forms <- if (outcome_family == "binomial") {
    c("rd", "log(rr)", "log(or)")
  } else {
    "diff"
  }
  rows <- list()
  for (j in seq_along(lev)[-1]) {
    for (f in forms) {
      val <- switch(
        f,
        rd = mu[[j]] - ref,
        diff = mu[[j]] - ref,
        "log(rr)" = log(mu[[j]]) - log(ref),
        "log(or)" = stats::qlogis(mu[[j]]) - stats::qlogis(ref)
      )
      rows[[length(rows) + 1]] <- data.frame(
        effect = f,
        comparison = paste0(lev[[j]], " vs ", lev[[1]]),
        estimate = val,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# Match an estimates table to a reference table on the (effect, comparison) key
# and compare the point estimates.
expect_categorical_estimate_match <- function(
  estimates,
  reference,
  tolerance = 1e-8
) {
  key_got <- paste(estimates$effect, estimates$comparison)
  key_ref <- paste(reference$effect, reference$comparison)
  got <- estimates$estimate[match(key_ref, key_got)]
  expect_equal(got, reference$estimate, tolerance = tolerance)
}

cat_estimates_columns <- c(
  "effect",
  "comparison",
  "estimate",
  "std.err",
  "z",
  "ci.lower",
  "ci.upper",
  "conf.level",
  "p.value"
)

# ---- end-to-end ate ---------------------------------------------------------

test_that("ipw() runs categorical ate end to end and auto-detects the estimand", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")

  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_s3_class(res, "ipw")
  expect_equal(res$se_method, "mestimation")
  expect_equal(res$estimand, "ate")
  expect_false(is.null(res$fit))
})

# ---- estimates table shape --------------------------------------------------

test_that("ipw() categorical estimates table adds a comparison column", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()

  mods <- fit_categorical_models(dat, "ate")
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  est <- res$estimates

  # the binary eight-column contract plus a comparison column
  expect_named(est, cat_estimates_columns)

  # three-level exposure: two non-reference comparisons, each vs the reference
  expect_equal(unique(est$comparison), c("b vs a", "c vs a"))

  # binomial outcome: rd/log(rr)/log(or) per comparison, six rows for three levels
  expect_equal(nrow(est), 6L)
  expect_equal(est$effect, rep(c("rd", "log(rr)", "log(or)"), times = 2))
  expect_equal(
    est$comparison,
    rep(c("b vs a", "c vs a"), each = 3)
  )
  expect_true(all(est$conf.level == 0.95))

  # gaussian outcome: a single difference per comparison, two rows
  mods_g <- fit_categorical_models(dat, "ate", outcome_family = "gaussian")
  res_g <- ipw(mods_g$ps_mod, mods_g$outcome_mod)
  est_g <- res_g$estimates
  expect_equal(nrow(est_g), 2L)
  expect_equal(est_g$effect, c("diff", "diff"))
  expect_equal(est_g$comparison, c("b vs a", "c vs a"))
})

# ---- point-estimate parity with the plug-in g-computation -------------------

test_that("ipw() categorical point estimates match the plug-in g-computation", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  for (family in c("binomial", "gaussian")) {
    dat <- sim_categorical()
    mods <- fit_categorical_models(dat, "ate", outcome_family = family)
    res <- ipw(mods$ps_mod, mods$outcome_mod)
    ref <- plugin_categorical(mods$outcome_mod, dat, mods$lev, family)
    expect_categorical_estimate_match(res$estimates, ref, tolerance = 1e-8)
  }
})

# ---- PSweight oracle --------------------------------------------------------

# Compare risk differences and their standard errors per comparison against
# PSweight for ate (weight = "IPW") and ato (weight = "overlap"). PSweight
# reports Hajek weighted group means (no outcome model), so the outcome model
# here is exposure-only (a saturated y ~ a); its g-computation marginal means
# equal those Hajek means, making the two comparable. PSweight fits its own
# multinomial model internally, so the fits differ slightly: point estimates are
# compared at 1e-4 and standard errors at 5 percent relative.
test_that("ipw() categorical risk differences match PSweight for ate and ato", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  skip_if_not_installed("PSweight")

  dat <- sim_categorical()
  contrast <- rbind("b-a" = c(-1, 1, 0), "c-a" = c(-1, 0, 1))

  configs <- list(
    list(estimand = "ate", weight = "IPW"),
    list(estimand = "ato", weight = "overlap")
  )

  for (cfg in configs) {
    mods <- fit_categorical_models(dat, cfg$estimand, covariates = FALSE)
    res <- ipw(mods$ps_mod, mods$outcome_mod)

    psw_fit <- PSweight::PSweight(
      ps.formula = a ~ x1 + x2,
      yname = "y",
      data = dat,
      weight = cfg$weight,
      family = "binomial"
    )
    psw_sum <- summary(psw_fit, type = "DIF", contrast = contrast)
    psw_rd <- psw_sum$estimates[, "Estimate"]
    psw_se <- psw_sum$estimates[, "Std.Error"]

    rd_rows <- res$estimates[res$estimates$effect == "rd", ]
    rd_rows <- rd_rows[match(c("b vs a", "c vs a"), rd_rows$comparison), ]

    expect_equal(rd_rows$estimate, unname(psw_rd), tolerance = 1e-4)
    rel_se <- abs(rd_rows$std.err - psw_se) / psw_se
    expect_true(
      all(rel_se < 0.05),
      label = paste0(
        cfg$estimand,
        " SE rel diff: ",
        paste(format(rel_se, digits = 3), collapse = ", ")
      )
    )
  }
})

# ---- att and atu with a focal level -----------------------------------------

test_that("ipw() categorical att detects the focal level from the psw attribute", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  # att weights carry the focal level as the psw focal_category attribute
  mods <- fit_categorical_models(dat, "att", focal_level = "b")

  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_equal(res$estimand, "att")
  expect_equal(unique(res$estimates$comparison), c("b vs a", "c vs a"))
  ref <- plugin_categorical(
    mods$outcome_mod,
    dat,
    mods$lev,
    "binomial",
    estimand = "att",
    ps = ps_matrix_named(mods$ps_mod),
    focal = "b"
  )
  expect_categorical_estimate_match(res$estimates, ref, tolerance = 1e-8)
})

test_that("ipw() categorical atu accepts an explicit focal level argument", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "atu", focal_level = "c")

  res <- ipw(mods$ps_mod, mods$outcome_mod, .focal_level = "c")

  expect_equal(res$estimand, "atu")
  ref <- plugin_categorical(
    mods$outcome_mod,
    dat,
    mods$lev,
    "binomial",
    estimand = "atu",
    ps = ps_matrix_named(mods$ps_mod),
    focal = "c"
  )
  expect_categorical_estimate_match(res$estimates, ref, tolerance = 1e-8)
})

test_that("ipw() categorical att errors when the focal level is unavailable", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  # weights created with a focal level, then stripped to plain numeric so the
  # focal_category attribute is gone; with estimand = "att" supplied explicitly
  # and no .focal_level argument there is no focal level to be found
  mods <- fit_categorical_models(
    dat,
    "att",
    focal_level = "b",
    strip_weights = TRUE
  )

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, estimand = "att"),
    class = "propensity_focal_required_error"
  )
})

# ---- atm and entropy --------------------------------------------------------

test_that("ipw() categorical atm and entropy return finite, ordered intervals", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  for (estimand in c("atm", "entropy")) {
    dat <- sim_categorical()
    mods <- fit_categorical_models(dat, estimand)
    res <- ipw(mods$ps_mod, mods$outcome_mod)
    est <- res$estimates

    expect_true(all(is.finite(est$std.err) & est$std.err > 0))
    expect_true(all(est$ci.lower < est$estimate & est$estimate < est$ci.upper))
  }
})

# Finiteness says nothing about where the estimates sit. The atm and entropy
# target populations are not subsets of the rows, so tilted standardization is
# the only route to them: a flat mean of the counterfactual predictions would
# report the ate contrast under either name. Both fixtures adjust the outcome
# model for x1, which is also what the propensity score model reads, so the tilt
# genuinely moves the answer. The flat standardization sits 1.2e-2 away from the
# tilted one for atm and 3.4e-3 for entropy, and the ato tilt in place of the
# atm one moves the contrasts by 5.5e-3, all far outside the 1e-8 tolerance the
# pins below are compared at.

test_that("ipw() categorical atm and entropy marginal means are the tilt-standardized predictions", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  for (estimand in c("atm", "entropy")) {
    dat <- sim_categorical()
    mods <- fit_categorical_models(dat, estimand)
    res <- ipw(mods$ps_mod, mods$outcome_mod)

    ps <- ps_matrix_named(mods$ps_mod)
    h <- ps_tilt(ps, estimand)

    # The tilt is the one piece of arithmetic the reference reads from the
    # package, so it is checked against its definition rather than taken on
    # trust: the smallest level probability, and the entropy of the probability
    # vector.
    h_hand <- if (estimand == "atm") {
      pmin(ps[, 1], ps[, 2], ps[, 3])
    } else {
      -rowSums(ps * log(ps))
    }
    expect_equal(h, unname(h_hand), tolerance = 1e-12)

    # The marginal means are free parameters of the stacked system, reported
    # under mu_<level> in the solved theta, so they can be pinned directly
    # rather than through the contrasts they feed.
    mu_ref <- vapply(
      mods$lev,
      function(l) {
        d <- dat
        d$a <- factor(l, levels = mods$lev)
        weighted.mean(
          predict(mods$outcome_mod, newdata = d, type = "response"),
          h
        )
      },
      numeric(1)
    )
    mu_got <- coef(res$fit)[paste0("mu_", mods$lev)]

    expect_equal(unname(mu_got), unname(mu_ref), tolerance = 1e-8)
  }
})

test_that("ipw() categorical atm and entropy point estimates match the tilt-standardized plug-in", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  for (estimand in c("atm", "entropy")) {
    dat <- sim_categorical()
    mods <- fit_categorical_models(dat, estimand)
    res <- ipw(mods$ps_mod, mods$outcome_mod)

    expect_equal(res$estimand, estimand)
    ref <- plugin_categorical(
      mods$outcome_mod,
      dat,
      mods$lev,
      "binomial",
      estimand = estimand,
      ps = ps_matrix_named(mods$ps_mod)
    )
    expect_categorical_estimate_match(res$estimates, ref, tolerance = 1e-8)
  }
})

test_that("ipw() categorical atm and entropy match the tilt-standardized plug-in for a gaussian outcome", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The gaussian companion, which reports a single difference per comparison
  # rather than the three binomial effects. The outcome model carries an
  # exposure-by-covariate interaction, and that is what makes this arm sensitive
  # to the standardization: a linear model without one predicts a difference
  # that is the exposure coefficient at every unit, so its marginal contrast is
  # that coefficient whichever population it is standardized to, and the tilt
  # cannot be observed in it at all. The fixture's `yc ~ a + x1` was measured
  # agreeing with the flat standardization to 2e-15.
  for (estimand in c("atm", "entropy")) {
    dat <- sim_categorical()
    ps_mod <- fit_ps_multinom(dat)
    ps <- ps_matrix_named(ps_mod)
    wts <- categorical_weights(ps, dat$a, estimand)
    outcome_mod <- lm(yc ~ a * x1, data = dat, weights = wts)

    res <- ipw(ps_mod, outcome_mod)

    expect_equal(res$estimand, estimand)
    ref <- plugin_categorical(
      outcome_mod,
      dat,
      levels(dat$a),
      "gaussian",
      estimand = estimand,
      ps = ps
    )
    expect_categorical_estimate_match(res$estimates, ref, tolerance = 1e-8)
  }
})

# ---- stabilized ate ---------------------------------------------------------

test_that("ipw() runs stabilized categorical ate and matches the plug-in", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate", stabilize = TRUE)

  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_equal(res$estimand, "ate")
  ref <- plugin_categorical(mods$outcome_mod, dat, mods$lev, "binomial")
  expect_categorical_estimate_match(res$estimates, ref, tolerance = 1e-8)
  expect_true(all(is.finite(res$estimates$std.err) & res$estimates$std.err > 0))
})

# ---- linearization is unsupported for categorical ---------------------------

test_that("ipw() rejects linearization for a categorical exposure", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization"),
    class = "propensity_method_error"
  )
})

# ---- print and as.data.frame ------------------------------------------------

test_that("ipw() categorical print output is stable", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_snapshot(print(res))
})

test_that("as.data.frame(exponentiate = TRUE) relabels ratios per comparison", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  df <- as.data.frame(res, exponentiate = TRUE)

  # the ratio rows are relabelled while the comparison column is preserved
  expect_equal(df$effect, rep(c("rd", "rr", "or"), times = 2))
  expect_equal(df$comparison, rep(c("b vs a", "c vs a"), each = 3))
})

# ---- a transformed exposure in the outcome formula ---------------------------
#
# The two extraction routes treat a transformed exposure differently, and the
# difference is the intended behavior rather than an inconsistency.
#
# With `.data`, each counterfactual design is built by setting the exposure
# column to a level and rebuilding the outcome design, so `model.matrix`
# re-evaluates the transformation at the level being set. A term reading
# `as.numeric(a == "c")` becomes 1 in the design for level "c" and 0 elsewhere,
# which is what g-computation on that model means.
#
# Without `.data` there is no exposure column to set: the model frame holds the
# transformed term under its own label and the exposure cannot be read back
# through it, so the call errors and says to supply `.data`.
#
# The dual-indicator no-intercept coding is pinned running by "ipw() categorical
# accepts a no-intercept outcome model adjusted for a covariate"; what is pinned
# here is that the re-evaluation produces the same answer as writing the same
# design without a transformation.

test_that("a transformed exposure with .data matches the untransformed equivalent", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)

  # One indicator per non-reference level, with an intercept, is the treatment
  # coding `y ~ a + x1` writes for itself. The two designs hold the same numbers
  # in the same order, so the only difference is that one route has to
  # re-evaluate the transformation at each counterfactual level.
  transformed <- glm(
    y ~ as.numeric(a == "b") + as.numeric(a == "c") + x1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = ctrl
  )
  factor_coded <- glm(
    y ~ a + x1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = ctrl
  )

  # the two designs really are the same numbers, attributes aside
  expect_equal(
    unname(model.matrix(transformed)[,]),
    unname(model.matrix(factor_coded)[,])
  )

  transformed_estimates <- ipw(
    mods$ps_mod,
    transformed,
    .data = dat
  )$estimates
  factor_estimates <- ipw(mods$ps_mod, factor_coded)$estimates

  expect_equal(
    transformed_estimates$estimate,
    factor_estimates$estimate,
    tolerance = 1e-12
  )
  expect_equal(transformed_estimates, factor_estimates, tolerance = 1e-10)
})

test_that("a transformed exposure without .data errors and asks for it", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The other half of the contract. The exposure is present in the formula, so
  # the guard requiring it is satisfied; what fails is that the model frame
  # holds only the transformed terms, which nothing on this route recomputes.
  # It used to be reported as the exposure column missing from that frame.
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  transformed <- glm(
    y ~ as.numeric(a == "b") + as.numeric(a == "c") + x1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  # the exposure name is genuinely absent from the model frame
  expect_false("a" %in% names(stats::model.frame(transformed)))

  err <- expect_error(
    ipw(mods$ps_mod, transformed),
    class = "propensity_ipw_exposure_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "\"a\"", fixed = TRUE)
  expect_match(msg, "as.numeric(a == \"b\")", fixed = TRUE)
  expect_match(msg, ".data", fixed = TRUE)
})

# ---- a matrix-response outcome model -----------------------------------------
#
# The response-shape guard has one call per exposure type, so a fix applied only
# to the binary path would look complete. The binary companions are "a
# matrix-response outcome model errors on its shape without .data" and its
# .data counterpart in test-ipw-se-method.R.

test_that("ipw() categorical rejects a matrix-response outcome model", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  dat$y1 <- rbinom(nrow(dat), 5, 0.4)
  dat$y0 <- 5 - dat$y1
  ps_mod <- fit_ps_multinom(dat)
  ps <- ps_matrix_named(ps_mod)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps, dat$a, exposure_type = "categorical")
  )
  matrix_outcome <- glm(
    cbind(y1, y0) ~ a,
    data = dat,
    family = binomial(),
    weights = wts
  )

  err <- expect_error(
    ipw(ps_mod, matrix_outcome),
    class = "propensity_ipw_response_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "matrix response", fixed = TRUE)
})

# ---- a .data whose covariates rebuild a differently shaped design ------------
#
# The propensity design is rebuilt from `.data` while the coefficients come from
# the fit, and the multinomial coefficient block is named by pairing each level
# with each design column. If `.data` gives a covariate a different type than the
# fit saw, the rebuild can produce a different number of columns, and the pairing
# produces more names than there are coefficients. What surfaces is a raw error
# about a names attribute, which says nothing about `.data` or about the column
# that changed.
#
# A change in the number of columns is not the only thing that matters, and the
# pin below records the measurement that settled it. A two-level factor supplied
# where the fit saw a numeric column gives one dummy column, the same width the
# numeric took, but the dummy is zero and one whatever the numeric held. It is
# the same numbers only when the fitted column was itself coded zero and one,
# which the rebuild cannot see and the width cannot decide, so the column's type
# is what the guard reads.

test_that("ipw() categorical rejects a .data whose covariate types skew the design", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")

  # x2 was numeric in the fit, so the design had one column for it. As a
  # three-level factor it expands to two, and the design gains a column the
  # coefficients have no counterpart for.
  dat_skew <- dat
  dat_skew$x2 <- factor(
    ifelse(dat$x2 == 1, "hi", ifelse(dat$x1 > 0, "mid", "lo")),
    levels = c("lo", "mid", "hi")
  )
  expect_equal(ncol(model.matrix(~ x1 + x2, data = dat)), 3)
  expect_equal(ncol(model.matrix(~ x1 + x2, data = dat_skew)), 4)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat_skew),
    class = "propensity_ipw_data_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, ".data", fixed = TRUE)
})

test_that("the skewed-design error names the column and both types", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  dat_skew <- dat
  dat_skew$x2 <- factor(
    ifelse(dat$x2 == 1, "hi", ifelse(dat$x1 > 0, "mid", "lo")),
    levels = c("lo", "mid", "hi")
  )

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat_skew)
  )
})

test_that("ipw() categorical reports a .data missing a model covariate", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The column check once covered only the exposure and the outcome, so a
  # covariate missing from .data reached model.matrix and raised a raw
  # object-not-found naming the variable but nothing else.
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  dat_missing <- dat
  dat_missing$x2 <- NULL

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat_missing),
    class = "propensity_columns_exist_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "x2", fixed = TRUE)
})

test_that("ipw() categorical rejects a .data covariate recoded to the same width", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The design keeps its width here, so nothing about its shape can decide the
  # case. The two-level factor gives one dummy column of zeros and ones, which
  # reproduces the fitted numbers only because this fixture's x2 is itself coded
  # zero and one; the same recoding of a covariate coded one and two rebuilds a
  # design the model was never fit on. Where that lands depends on which model
  # reads the covariate: an outcome-model covariate moves the estimates with
  # nothing signaled, while a propensity-model covariate moves the recomputed
  # weights and used to surface at the weight-consistency preflight, whose
  # message is about how the weights were built. The rebuild cannot tell the two
  # recodings apart, so the type decides.
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")

  dat_recoded <- dat
  dat_recoded$x2 <- factor(
    ifelse(dat$x2 == 1, "hi", "lo"),
    levels = c("lo", "hi")
  )
  expect_equal(ncol(model.matrix(~ x1 + x2, data = dat_recoded)), 3)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat_recoded),
    class = "propensity_ipw_data_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "\\bx2\\b")
  expect_match(msg, "factor", fixed = TRUE)
})

test_that("a log-transformed categorical outcome response through .data matches the model-frame route", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The categorical companion to the binary pins in test-ipw-se-method.R. The
  # response is read out of `.data` by the shared extraction, so a transformation
  # of a column has to be computed rather than looked up under the name of the
  # function wrapping it. The exponentiated outcome keeps `log()` defined on
  # every observation.
  dat <- sim_categorical()
  dat$yp <- exp(dat$yc)
  ps_mod <- fit_ps_multinom(dat)
  wts <- categorical_weights(ps_matrix_named(ps_mod), dat$a, "ate")
  outcome_mod <- lm(log(yp) ~ a + x1, data = dat, weights = wts)

  expect_equal(
    ipw(ps_mod, outcome_mod, .data = dat)$estimates,
    ipw(ps_mod, outcome_mod)$estimates,
    tolerance = 1e-10
  )
})

# ---- a wrong-sized .data is a data problem -----------------------------------
#
# The categorical path sizes the exposure, the outcome, and the designs to
# `.data`, so a `.data` that does not line up with the fitted models is only
# noticed later, by the weight-consistency preflight, whose message is about
# weights and focal levels rather than about the data frame that was passed. The
# binary companion is "mestimation rejects a .data with too few rows" in
# test-ipw-mestimation.R.

test_that("ipw() categorical rejects a .data with too few rows", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat[-1, ]),
    class = "propensity_ipw_data_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "one row per observation", fixed = TRUE)
})

# ---- .data when the multinom model frame is unavailable ----------------------

test_that("ipw() categorical uses .data when the ps model frame is gone", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()

  # nnet::multinom stores no model frame, so its design must be rebuilt from the
  # fitting call. Fit it inline against a binding that is then removed, leaving
  # the model frame unreconstructable, and keep a copy of the data for .data.
  # Fitting inline (not through a helper) keeps the fitting frame from retaining
  # the data through the helper's arguments.
  ps_mod <- nnet::multinom(
    a ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_named <- ps_matrix_named(ps_mod)
  wts <- categorical_weights(ps_named, dat$a, "ate")
  outcome_mod <- glm(
    y ~ a + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  dat_copy <- dat
  rm(dat)

  # Without .data the propensity score design cannot be reconstructed, so the
  # error names the cause and directs the user to supply .data.
  expect_error(
    ipw(ps_mod, outcome_mod),
    class = "propensity_ipw_data_error"
  )

  # Supplying .data rebuilds the design and matches the default extraction path.
  res_data <- ipw(ps_mod, outcome_mod, .data = dat_copy)
  expect_s3_class(res_data, "ipw")
  expect_equal(res_data$estimand, "ate")
  ref <- plugin_categorical(
    outcome_mod,
    dat_copy,
    levels(dat_copy$a),
    "binomial"
  )
  expect_categorical_estimate_match(res_data$estimates, ref, tolerance = 1e-8)
})

# ---- guard against a weighted multinom --------------------------------------

test_that("ipw() rejects a multinom fit with case weights", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  ps_mod_wtd <- nnet::multinom(
    a ~ x1 + x2,
    data = dat,
    weights = runif(nrow(dat), 0.5, 2),
    trace = FALSE
  )
  mods <- fit_categorical_models(dat, "ate")

  expect_error(
    ipw(ps_mod_wtd, mods$outcome_mod),
    class = "propensity_ipw_ps_weights_error"
  )
})

# ---- guard against an invalid focal level -----------------------------------

test_that("ipw() rejects a focal level that is not an exposure level", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "att", focal_level = "b")

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .focal_level = "zzz"),
    class = "propensity_focal_level_error"
  )
})

# ---- guard against an estimand the weights record that is not one ------------

test_that("ipw() rejects weights recording an estimand it does not know", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  # `psw()` records the estimand it is given, so weights built by hand reach the
  # categorical tilt with a value it has no branch for
  banana_wts <- psw(as.double(mods$wts), estimand = "banana")
  outcome_mod <- fit_outcome(dat, banana_wts)

  err <- expect_error(
    ipw(mods$ps_mod, outcome_mod),
    class = "propensity_ipw_estimand_error"
  )
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(err)),
    "banana",
    fixed = TRUE
  )
})

# ---- outcome-family validation ----------------------------------------------

test_that("ipw_spec_categorical rejects an unsupported outcome family", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical()
  ps_mod <- fit_ps_multinom(dat)
  ps_named <- ps_matrix_named(ps_mod)
  wts <- categorical_weights(ps_named, dat$a, "ate")

  withr::local_seed(1)
  dat$ycount <- rpois(nrow(dat), exp(0.1 + 0.3 * (dat$a == "c") + 0.2 * dat$x1))
  outcome_mod <- suppressWarnings(
    glm(ycount ~ a, data = dat, family = poisson(), weights = wts)
  )

  # Without the family guard a poisson outcome model is stacked as a binomial
  # score.
  expect_error(
    ipw_spec_categorical(ps_mod, outcome_mod),
    class = "propensity_ipw_family_error"
  )
})

# ---- factor outcome response ------------------------------------------------

test_that("categorical mestimation matches a factor outcome response to the numeric fit", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  dat$yf <- factor(ifelse(dat$y == 1, "yes", "no"), levels = c("no", "yes"))
  ps_mod <- fit_ps_multinom(dat)
  wts <- categorical_weights(ps_matrix_named(ps_mod), dat$a, "ate")
  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)

  num <- glm(
    y ~ a,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = ctrl
  )
  fac <- suppressWarnings(
    glm(
      yf ~ a,
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

# ---- length-2 .focal_level --------------------------------------------------

test_that("ipw() categorical errors informatively on a length-2 .focal_level", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "att", focal_level = "b")

  # The length assertion runs ahead of `!focal_level %in% levs`, which sits
  # inside an `&&` and would otherwise raise a raw base error ("'length = 2' in
  # coercion to 'logical(1)'") that does not name .focal_level. The length-1
  # case is covered by "ipw() categorical atu accepts an explicit focal level
  # argument".
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .focal_level = c("b", "c")),
    class = "propensity_error",
    regexp = "focal_level"
  )
})

# ---- non-NULL ps_link -------------------------------------------------------

test_that("ipw() categorical rejects a non-NULL ps_link", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")

  # ps_link is meaningful only for a binomial glm on the binary path; a multinom
  # propensity model must reject it rather than silently ignore it. The default
  # ps_link = NULL is covered by "ipw() runs categorical ate end to end and
  # auto-detects the estimand".
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit"),
    class = "propensity_error",
    regexp = "ps_link"
  )
})

# ---- offset guard -----------------------------------------------------------

test_that("ipw() categorical rejects an outcome model with an offset term", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  dat$off <- 0.3 * dat$x1
  mods <- fit_categorical_models(dat, "ate")

  # The offset guard (check_ipw_offset) is wired at ipw_spec_categorical entry;
  # the stacked outcome score does not thread an offset. Both offset routes must
  # be rejected on the categorical path.
  om_formula <- glm(
    y ~ a + offset(off),
    data = dat,
    family = quasibinomial(),
    weights = mods$wts
  )
  expect_error(
    ipw(mods$ps_mod, om_formula),
    class = "propensity_ipw_offset_error"
  )

  om_arg <- glm(
    y ~ a,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    offset = off
  )
  expect_error(
    ipw(mods$ps_mod, om_arg),
    class = "propensity_ipw_offset_error"
  )
})

# ---- the outcome model must contain the exposure -----------------------------
#
# As on the binary path, supplying .data bypasses the check_exposure() guard, so
# an outcome model fit without the exposure would produce one identical
# counterfactual design per level and collapse every contrast to zero. The
# categorical path rejects it with the same classed error.

test_that("ipw() categorical rejects an outcome model that omits the exposure when .data is supplied", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  no_exposure <- glm(
    stats::reformulate("x1", response = "y"),
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  expect_error(
    ipw(mods$ps_mod, no_exposure, .data = dat),
    class = "propensity_ipw_exposure_error"
  )
})

test_that("ipw() categorical accepts an outcome model containing the exposure when .data is supplied", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The guard must not fire on the ordinary case: the same call with the
  # exposure on the right-hand side runs through to an ipw object.
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")

  res <- ipw(mods$ps_mod, mods$outcome_mod, .data = dat)

  expect_s3_class(res, "ipw")
  expect_equal(res$estimand, "ate")
})

# ---- exposure levels come from the fitted model, not from .data --------------
#
# nnet::multinom records its training levels in ps_mod$lev and lays its
# coefficient rows out in that order. The exposure column in .data carries no
# such guarantee: a character column resolves alphabetically and a factor can be
# releveled. The indicator matrix, the counterfactual designs, and the reference
# level must follow the fitted model's ordering rather than whatever ordering
# .data happens to imply. The factor-.data control is already pinned by "ipw()
# categorical uses .data when the ps model frame is gone" and by "ipw()
# categorical accepts an outcome model containing the exposure when .data is
# supplied".

# sim_categorical() builds levels c("a", "b", "c"), which already agree with the
# alphabetical ordering a character column would produce. Relabel them so the
# fitted order (low, mid, high) and the alphabetical order (high, low, mid)
# disagree, leaving the values and the model fit otherwise unchanged.
sim_categorical_nonalpha <- function(...) {
  dat <- sim_categorical(...)
  nonalpha <- c("low", "mid", "high")
  dat$a <- factor(nonalpha[as.integer(dat$a)], levels = nonalpha)
  dat
}

test_that("ipw() categorical resolves a character .data exposure against the fitted levels", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical_nonalpha()
  mods <- fit_categorical_models(dat, "ate")

  dat_chr <- dat
  dat_chr$a <- as.character(dat_chr$a)

  # the two orderings genuinely disagree, so the fixture exercises the contract
  expect_equal(mods$ps_mod$lev, c("low", "mid", "high"))
  expect_equal(levels(as.factor(dat_chr$a)), c("high", "low", "mid"))

  res_fac <- ipw(mods$ps_mod, mods$outcome_mod, .data = dat)
  res_chr <- ipw(mods$ps_mod, mods$outcome_mod, .data = dat_chr)

  # same values, same model: the character column must give the same answer,
  # including the reference level the comparisons are formed against
  expect_equal(res_chr$estimates, res_fac$estimates, tolerance = 1e-8)
  expect_equal(
    unique(res_chr$estimates$comparison),
    c("mid vs low", "high vs low")
  )
})

test_that("ipw() categorical ignores a releveled .data exposure ordering", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical_nonalpha()
  mods <- fit_categorical_models(dat, "ate")

  # the same values in a different factor ordering: releveling a column cannot
  # change which level the fitted model treats as the reference
  dat_relevel <- dat
  dat_relevel$a <- factor(
    as.character(dat$a),
    levels = c("high", "low", "mid")
  )

  res_fac <- ipw(mods$ps_mod, mods$outcome_mod, .data = dat)
  res_relevel <- ipw(mods$ps_mod, mods$outcome_mod, .data = dat_relevel)

  expect_equal(res_relevel$estimates, res_fac$estimates, tolerance = 1e-8)
  expect_equal(
    unique(res_relevel$estimates$comparison),
    c("mid vs low", "high vs low")
  )
})

test_that("ipw() categorical rejects a .data exposure value the model never saw", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical_nonalpha()
  mods <- fit_categorical_models(dat, "ate")

  dat_unknown <- dat
  dat_unknown$a <- as.character(dat_unknown$a)
  dat_unknown$a[[1]] <- "unknown"

  # Without level resolution, the extra level reached the theta layout and
  # raised a raw base error about a names-length mismatch. It must be a classed
  # propensity error that names the value the model cannot score.
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat_unknown),
    class = "propensity_error",
    regexp = "unknown"
  )
})

# The tests above hand the exposure to ipw() through `.data`. Without `.data`
# the exposure is read back out of the propensity score model's frame, which
# returns a character column as character, carrying no levels of its own.
# ipw_categorical_exposure_factor() is what makes it a factor, and it imposes
# ps_mod$lev rather than deriving an order from the values. The alphabetical
# order is therefore fixed at fit time, by the as.factor() multinom applies to
# its response before recording $lev, and sim_categorical() labels its levels
# a, b, c, so that order matches the factor fit's: same propensity scores, same
# weights, same reference level. The two calls must therefore agree, which pins
# the character column as a supported way to fit the models rather than
# something that happens to work.

test_that("ipw() categorical accepts models fit on a character exposure column", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  dat_chr <- dat
  dat_chr$a <- as.character(dat$a)
  expect_type(dat_chr$a, "character")

  fit_both_models <- function(d) {
    ps_mod <- fit_ps_multinom(d)
    ps_named <- unname(predict(ps_mod, type = "probs"))
    # the fitted model's own level order, which is the only order available
    # when the source column carries no levels of its own
    colnames(ps_named) <- ps_mod$lev
    wts <- categorical_weights(ps_named, d$a, "ate")
    list(ps_mod = ps_mod, outcome_mod = fit_outcome(d, wts))
  }

  mods_fac <- fit_both_models(dat)
  mods_chr <- fit_both_models(dat_chr)
  expect_equal(mods_chr$ps_mod$lev, mods_fac$ps_mod$lev)

  res_fac <- ipw(mods_fac$ps_mod, mods_fac$outcome_mod)
  res_chr <- ipw(mods_chr$ps_mod, mods_chr$outcome_mod)

  expect_equal(res_chr$estimates, res_fac$estimates, tolerance = 1e-8)
  expect_equal(unique(res_chr$estimates$comparison), c("b vs a", "c vs a"))
})

# ---- weight-consistency hint wording ----------------------------------------

test_that("the categorical weight-consistency error omits the binary focal convention", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  spec <- ipw_spec_categorical(mods$ps_mod, mods$outcome_mod)
  doubled <- 2 * as.double(model.weights(model.frame(mods$outcome_mod)))

  err <- expect_error(
    ipw_check_weight_consistency(spec, doubled),
    class = "propensity_ipw_weights_mismatch_error"
  )

  # The focal-level hint is worded per exposure type. A categorical exposure
  # takes its focal level from the weights or from .focal_level and has no
  # sorted-level convention, so the binary sentence must not appear here. The
  # binary wording is pinned by "the weight-consistency error names a focal
  # level as a possible cause" in test-ipw-mestimation.R.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(
    msg,
    "Weights built with a different `.focal_level`",
    fixed = TRUE
  )
  expect_false(grepl("second sorted level", msg, fixed = TRUE))
})

# ---- degenerate counterfactual designs --------------------------------------
#
# One counterfactual design is built per level by setting the exposure factor to
# that level and rebuilding the outcome design. A numeric no-intercept coding of
# the exposure, y ~ as.numeric(a == "c") - 1, leaves the designs at "a" and at
# "b" identically zero, so both marginal means are pinned at plogis(0) = 0.5:
# the "b vs a" contrast collapses to a rounding-error zero and "c vs a" is
# measured against a constant rather than against the data, with nothing
# signaled.
#
# The condition is a design that is identically zero, not a missing intercept.
# The saturated factor coding y ~ a - 1 is a reparameterization whose designs
# are the level indicators, and it must keep working.

test_that("ipw() categorical rejects an outcome model whose counterfactual design is identically zero", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  # The transformed term is unreadable from the model frame, so this coding
  # reaches the counterfactual designs only through .data.
  collapsed <- glm(
    y ~ as.numeric(a == "c") - 1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(mods$ps_mod, collapsed, .data = dat),
    class = "propensity_error",
    regexp = "identically zero"
  )

  # The message must name every level whose design is degenerate. Whitespace is
  # normalized because cli wraps the bullet.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "`a` to \"a\" and \"b\"", fixed = TRUE)
})

test_that("the categorical degenerate-design error names the pinned levels", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  collapsed <- glm(
    y ~ as.numeric(a == "c") - 1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, collapsed, .data = dat)
  )
})

test_that("ipw() categorical accepts a saturated no-intercept outcome model", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # y ~ a - 1 is a reparameterization of y ~ a, not a model without a baseline:
  # its counterfactual designs are the level indicators, which are never zero.
  # It must give the same answer as the with-intercept fit.
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)
  saturated <- glm(
    y ~ a - 1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = ctrl
  )
  with_intercept <- glm(
    y ~ a,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = ctrl
  )

  res_sat <- ipw(mods$ps_mod, saturated)
  res_int <- ipw(mods$ps_mod, with_intercept)

  expect_s3_class(res_sat, "ipw")
  expect_equal(res_sat$estimates$effect, res_int$estimates$effect)
  expect_equal(res_sat$estimates$comparison, res_int$estimates$comparison)
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

test_that("ipw() categorical accepts a no-intercept outcome model adjusted for a covariate", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The covariate is what keeps this model out of both guards, which is why the
  # exposure is coded as numeric indicators rather than as the factor. Dropping
  # an intercept does not cost a factor its baseline: R answers y ~ a + x1 - 1
  # with full dummy coding, so the designs there are nonzero whether or not x1
  # is present and the case reduces to the saturated test above. With one
  # indicator per non-reference level the reference design is the indicators at
  # zero, so y ~ ind_b + ind_c - 1 alone would trip the zero-design guard at
  # "a"; x1 is the only column keeping that design alive.
  #
  # The same designs are also pairwise distinct, "a" at the origin, "b" and "c"
  # one step along their own axis, so this doubles as the pin that the
  # indistinguishable-design guard leaves a legitimate reparameterization alone.
  #
  # This is the categorical companion to "mestimation accepts a no-intercept
  # outcome model adjusted for a covariate" in test-ipw-mestimation.R, where
  # numeric z gives the binary path the same reference-at-the-origin design.
  # The transformed terms are unreadable from the model frame, so this coding
  # reaches the counterfactual designs only through .data.
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  adjusted <- glm(
    y ~ as.numeric(a == "b") + as.numeric(a == "c") + x1 - 1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  res <- ipw(mods$ps_mod, adjusted, .data = dat)

  expect_s3_class(res, "ipw")
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
  expect_true(all(res$estimates$std.err > 0))
})

# ---- indistinguishable counterfactual designs -------------------------------
#
# A counterfactual design can be free of zero columns and still fail to separate
# two levels. A single numeric indicator for one level, such as
# y ~ as.numeric(a == "c") + x1 - 1, leaves the designs at "a" and at "b"
# identical element for element: the indicator is zero at both and the covariate
# is untouched. Every unit's counterfactual prediction is then the same number at
# both levels, so the "b vs a" contrast is zero to rounding error and its
# standard error is NaN. The outcome model cannot represent a difference between
# those two levels, so per-level contrasts against them are not estimable and the
# call must be rejected instead of reported.
#
# The condition is pairwise equality of the designs, not the absence of a
# covariate or of an intercept. A coding that carries one indicator per
# non-reference level separates every pair and must keep working; that negative
# case is pinned by "ipw() categorical accepts a no-intercept outcome model
# adjusted for a covariate" above, which uses exactly that coding.

test_that("ipw() categorical rejects an outcome model whose counterfactual designs are pairwise identical", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  # The covariate keeps every design nonzero, so the zero-design guard has
  # nothing to catch; only the "a" and "b" designs coinciding is wrong here. The
  # transformed term is unreadable from the model frame, so this coding reaches
  # the counterfactual designs only through .data.
  indistinguishable <- glm(
    y ~ as.numeric(a == "c") + x1 - 1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(mods$ps_mod, indistinguishable, .data = dat),
    class = "propensity_ipw_exposure_error"
  )

  # The message must name the levels the model cannot tell apart. Whitespace is
  # normalized because cli wraps the bullet.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "identical counterfactual designs", fixed = TRUE)
  expect_match(msg, "`a` to \"a\" and \"b\"", fixed = TRUE)
})

test_that("the categorical indistinguishable-design error names the pinned levels", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  indistinguishable <- glm(
    y ~ as.numeric(a == "c") + x1 - 1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, indistinguishable, .data = dat)
  )
})

test_that("ipw() categorical rejects pairwise-identical counterfactual designs under an intercept", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # Adding the intercept back changes nothing about the defect: the "a" and "b"
  # designs still agree in every column, they are just no longer zero anywhere.
  # The same error must fire, which is what shows the guard keys on pairwise
  # equality rather than on a missing intercept.
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  indistinguishable <- glm(
    y ~ as.numeric(a == "c") + x1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(mods$ps_mod, indistinguishable, .data = dat),
    class = "propensity_ipw_exposure_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "identical counterfactual designs", fixed = TRUE)
  expect_match(msg, "`a` to \"a\" and \"b\"", fixed = TRUE)
})

# ---- a numeric-coded exposure in the outcome model --------------------------
#
# nnet::multinom is happy to model a numeric exposure, recording its levels as
# the character forms of the values. The counterfactual designs are then built
# by setting the exposure to a factor of those levels, so model.matrix expands
# it into dummy columns. An outcome model that holds the same exposure as a
# single numeric term has one slot for it in its coefficient vector, and the two
# no longer line up: the design gains a column per non-reference level while the
# coefficients do not.
#
# The two codings also mean different things. A numeric term fits one slope
# across the levels and assumes they are equally spaced; the counterfactual
# designs assume nothing of the sort. Rather than reconcile them, the outcome
# model must carry the exposure as a factor, which is what the error says.
#
# The guard's negative case is every factor-exposure test in this file, starting
# with "ipw() runs categorical ate end to end and auto-detects the estimand".

sim_numeric_categorical <- function(seed = 2024, n = 700) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  z <- sample(0:2, n, replace = TRUE)
  y <- rbinom(n, 1, plogis(-0.3 + 0.4 * (z == 1) + 0.8 * (z == 2) + 0.5 * x1))
  data.frame(x1, z, y)
}

# Propensity score model, categorical weights, and a weighted outcome model that
# holds the exposure as the bare numeric column. `double` recodes the exposure to
# double first, so the same fixture serves the integer and double cases.
fit_numeric_categorical_models <- function(dat, double = FALSE) {
  if (double) {
    dat$z <- as.numeric(dat$z)
  }
  ps_mod <- nnet::multinom(
    z ~ x1,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- ps_mod$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps, dat$z, exposure_type = "categorical")
  )
  outcome_mod <- glm(
    y ~ z + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts, dat = dat)
}

test_that("ipw() categorical rejects an outcome model holding the exposure as a numeric term", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_numeric_categorical()
  mods <- fit_numeric_categorical_models(dat)

  # Without a guard this reaches the estimating equations and dies on a raw
  # base "non-conformable arguments" from X %*% beta, with nothing naming the
  # exposure or the coding that caused it.
  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_exposure_error"
  )

  # Whitespace is normalized because cli wraps the bullet. Only the remedy is
  # pinned here; the full wording is pinned by the snapshot.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "as a factor", fixed = TRUE)
})

test_that("the categorical numeric-exposure error names the exposure and the remedy", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_numeric_categorical()
  mods <- fit_numeric_categorical_models(dat)

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, mods$outcome_mod)
  )
})

test_that("ipw() categorical rejects a numeric-term exposure when .data is supplied", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # Supplying .data changes how the exposure column is recovered but not how the
  # outcome model codes it, so the same error must fire on this route too.
  dat <- sim_numeric_categorical()
  mods <- fit_numeric_categorical_models(dat)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat),
    class = "propensity_ipw_exposure_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "as a factor", fixed = TRUE)
})

test_that("ipw() categorical rejects a double-coded exposure in the outcome model", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The defect is the numeric coding, not the storage type. A double column
  # reaches the model frame as the same single term an integer column does, so
  # it must be rejected the same way.
  dat <- sim_numeric_categorical()
  mods <- fit_numeric_categorical_models(dat, double = TRUE)
  expect_type(mods$dat$z, "double")

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_exposure_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "as a factor", fixed = TRUE)
})

# ---- both models must share the exposure's factor level order ---------------
#
# The counterfactual designs are built by setting the exposure to
# factor(level, levels = ps_mod$lev), so their dummy columns are laid out in the
# propensity score model's level order. The coefficients they multiply come from
# the outcome model's own coding, and the multiply is positional: nothing
# consults the column names. Fitting the outcome model on a releveled copy of
# the same exposure keeps the same level set and the same fit, and changes only
# the order, which is enough to pair each design column with the wrong
# coefficient. The estimates that come back are wrong, sometimes sign-flipped,
# with no warning and no NaN to notice.
#
# The same mismatch is reachable without anyone releveling anything: an outcome
# model fit on a character exposure column gets its levels from model.frame,
# which sorts them alphabetically. When the propensity score model was fit on a
# factor whose order is not alphabetical, the two disagree.
#
# The aligned case is pinned throughout this file. Two tests pin it specifically
# for a non-alphabetical level order, so a guard that compared against sorted
# levels rather than against the fitted order would fail them: "ipw() categorical
# resolves a character .data exposure against the fitted levels" and "ipw()
# categorical ignores a releveled .data exposure ordering".

test_that("ipw() categorical rejects an outcome model fit on a releveled exposure", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")

  # Same units, same weights, same fit: relevel() only reorders the factor, so
  # the outcome model's fitted values are unchanged and the weight-consistency
  # preflight has nothing to object to. Only the coefficient order moves.
  dat_relevel <- dat
  dat_relevel$a <- relevel(dat$a, ref = "c")
  releveled <- glm(
    y ~ a + x1,
    data = dat_relevel,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  # the fixture is the mismatch it claims to be
  expect_equal(mods$ps_mod$lev, c("a", "b", "c"))
  expect_equal(releveled$xlevels$a, c("c", "a", "b"))

  err <- expect_error(
    ipw(mods$ps_mod, releveled),
    class = "propensity_ipw_exposure_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "level order", fixed = TRUE)
})

test_that("the categorical level-order error names both orders", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  dat_relevel <- dat
  dat_relevel$a <- relevel(dat$a, ref = "c")
  releveled <- glm(
    y ~ a + x1,
    data = dat_relevel,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, releveled)
  )
})

test_that("ipw() categorical rejects a releveled outcome exposure when .data is supplied", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # Supplying .data changes how the exposure column is recovered, not how the
  # outcome model's coefficients are ordered, so this route is equally wrong and
  # must be rejected too.
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  dat_relevel <- dat
  dat_relevel$a <- relevel(dat$a, ref = "c")
  releveled <- glm(
    y ~ a + x1,
    data = dat_relevel,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(mods$ps_mod, releveled, .data = dat_relevel),
    class = "propensity_ipw_exposure_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "level order", fixed = TRUE)
})

test_that("ipw() categorical rejects a character outcome exposure ordered against the fitted levels", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # Nobody relevels anything here. The propensity score model is fit on a factor
  # ordered low, mid, high; the outcome model is fit on the same column stored as
  # character, which model.frame sorts to high, low, mid. The orders disagree and
  # the estimates are silently wrong. The .data route fails the same way.
  dat <- sim_categorical_nonalpha()
  mods <- fit_categorical_models(dat, "ate")
  dat_chr <- dat
  dat_chr$a <- as.character(dat$a)
  chr_outcome <- glm(
    y ~ a + x1,
    data = dat_chr,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  # the orders genuinely disagree, so the fixture exercises the contract
  expect_equal(mods$ps_mod$lev, c("low", "mid", "high"))
  expect_equal(chr_outcome$xlevels$a, c("high", "low", "mid"))

  err <- expect_error(
    ipw(mods$ps_mod, chr_outcome),
    class = "propensity_ipw_exposure_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "level order", fixed = TRUE)
})

test_that("the level-order guard passes a formula-transformed exposure through", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # y ~ factor(a) + x1 records its levels under the term label rather than under
  # the exposure name, so the guard has nothing to compare and must stay out of
  # the way. What this path does instead is owned elsewhere and deliberately not
  # pinned here: the assertion is only that the level-order guard is not what
  # stops it, whether it errors or runs.
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")
  wrapped <- glm(
    y ~ factor(a) + x1,
    data = dat,
    family = quasibinomial(),
    weights = mods$wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  expect_false("a" %in% names(wrapped$xlevels))

  msg <- tryCatch(
    {
      ipw(mods$ps_mod, wrapped, .data = dat)
      ""
    },
    error = function(e) gsub("[[:space:]]+", " ", conditionMessage(e))
  )
  expect_false(grepl("level order", msg, fixed = TRUE))
})

# ---- counterfactual rebuilds must honor the outcome model's contrasts -------
#
# The counterfactual designs are rebuilt by replacing the exposure column with
# factor(level, levels = levs), a bare factor carrying no contrasts attribute,
# and calling model.matrix without passing the fitted model's contrasts. An
# outcome model fit under non-default contrasts is therefore rebuilt under
# treatment coding while its coefficients stay in the coding it was fit under,
# and the positional multiply pairs them up wrongly.
#
# Changing a factor's contrasts is only a reparameterization. The fit, its
# fitted values, and every marginal mean are unchanged, so the estimates must
# equal those from the default-coded fit on the same data with the same weights.
# The reference is computed in the test rather than stored, so it cannot drift.
#
# Contrasts do not touch the level order, so the level-order guard passes these
# models: xlevels still records the levels in the order they were declared.

test_that("ipw() categorical is unchanged by a sum-coded exposure", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  dat_sum <- dat
  contrasts(dat_sum$a) <- contr.sum(3)

  default <- fit_categorical_models(dat, "ate")
  sum_coded <- fit_categorical_models(dat_sum, "ate")

  # a reparameterization, and one the level-order guard has no quarrel with
  expect_equal(names(coef(sum_coded$outcome_mod))[2:3], c("a1", "a2"))
  expect_equal(sum_coded$outcome_mod$xlevels$a, c("a", "b", "c"))
  expect_equal(
    unname(fitted(sum_coded$outcome_mod)),
    unname(fitted(default$outcome_mod)),
    tolerance = 1e-10
  )

  reference <- ipw(default$ps_mod, default$outcome_mod)$estimates
  frame_route <- ipw(sum_coded$ps_mod, sum_coded$outcome_mod)$estimates
  data_route <- ipw(
    sum_coded$ps_mod,
    sum_coded$outcome_mod,
    .data = dat_sum
  )$estimates

  # The point estimates are what the counterfactual rebuild determines, and they
  # agree to machine precision. The sandwich standard errors differ in their
  # eighth significant figure, which is the numerical derivative taken at a
  # different coefficient basis rather than a coding mismatch; a real one
  # reverses the sign of a contrast.
  expect_equal(frame_route$estimate, reference$estimate, tolerance = 1e-12)
  expect_equal(data_route$estimate, reference$estimate, tolerance = 1e-12)
  expect_equal(frame_route, reference, tolerance = 1e-6)
  expect_equal(data_route, reference, tolerance = 1e-6)
})

test_that("ipw() categorical is unchanged by an ordered-factor exposure", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # An ordered factor carries polynomial contrasts by default, so this reaches
  # the same defect without anyone setting a contrasts attribute. multinom takes
  # an ordered exposure and records the same levels in the same order, the
  # weight layer accepts it, and every exposure-coding guard passes, so nothing
  # stands between the poly-coded coefficients and the treatment-coded rebuild.
  dat <- sim_categorical()
  dat_ord <- dat
  dat_ord$a <- ordered(as.character(dat$a), levels = c("a", "b", "c"))

  default <- fit_categorical_models(dat, "ate")
  ord <- fit_categorical_models(dat_ord, "ate")

  expect_s3_class(dat_ord$a, "ordered")
  expect_equal(ord$ps_mod$lev, c("a", "b", "c"))
  expect_equal(ord$outcome_mod$xlevels$a, c("a", "b", "c"))
  expect_equal(names(coef(ord$outcome_mod))[2:3], c("a.L", "a.Q"))

  got <- ipw(ord$ps_mod, ord$outcome_mod)$estimates
  reference <- ipw(default$ps_mod, default$outcome_mod)$estimates

  # As above: the estimates match to machine precision once the rebuild carries
  # the polynomial coding, and the residual is sandwich noise at a different
  # coefficient basis.
  expect_equal(got$estimate, reference$estimate, tolerance = 1e-12)
  expect_equal(got, reference, tolerance = 1e-6)
})

# ---- arguments that fall into the dots --------------------------------------

test_that("ipw() multinom rejects arguments that fall into the dots", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "ate")

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

  # `.focal_level` is a multinom-only argument, so a misspelling of it has no
  # name to match on any method.
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, focal_level = "b"),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("ipw() multinom accepts every argument supplied by name", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "att", focal_level = "b")

  baseline <- ipw(mods$ps_mod, mods$outcome_mod)
  named <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = dat,
    estimand = "att",
    conf_level = 0.95,
    se_method = "mestimation",
    .focal_level = "b"
  )

  expect_equal(named$estimates, baseline$estimates)
})
