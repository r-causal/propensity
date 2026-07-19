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
ps_matrix_named <- function(ps_mod, dat) {
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- levels(dat$a)
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
  ps_named <- ps_matrix_named(ps_mod, dat)
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
# Mirrors ipw_categorical_tilt() in R/ipw-psi.R.
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
    ps = ps_matrix_named(mods$ps_mod, dat),
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
    ps = ps_matrix_named(mods$ps_mod, dat),
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
  ps_named <- ps_matrix_named(ps_mod, dat)
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

# ---- outcome-family validation ----------------------------------------------

test_that("ipw_spec_categorical rejects an unsupported outcome family", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical()
  ps_mod <- fit_ps_multinom(dat)
  ps_named <- ps_matrix_named(ps_mod, dat)
  wts <- categorical_weights(ps_named, dat$a, "ate")

  withr::local_seed(1)
  dat$ycount <- rpois(nrow(dat), exp(0.1 + 0.3 * (dat$a == "c") + 0.2 * dat$x1))
  outcome_mod <- suppressWarnings(
    glm(ycount ~ a, data = dat, family = poisson(), weights = wts)
  )

  # A poisson outcome model is silently stacked as a binomial score today.
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
  wts <- categorical_weights(ps_matrix_named(ps_mod, dat), dat$a, "ate")
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
  skip("pending focal_level length assertion")
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_categorical()
  mods <- fit_categorical_models(dat, "att", focal_level = "b")

  # Today a length-2 .focal_level reaches `!focal_level %in% levs` inside an `&&`,
  # raising a raw base error ("'length = 2' in coercion to 'logical(1)'") that
  # does not name .focal_level. It must become a classed error naming the
  # argument. The length-1 case is covered by "ipw() categorical atu accepts an
  # explicit focal level argument".
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .focal_level = c("b", "c")),
    class = "propensity_error",
    regexp = "focal_level"
  )
})
