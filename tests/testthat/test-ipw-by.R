# Effect modification through `.by` for a binary exposure on the M-estimation
# path. `ipw(wt_mod, outcome_mod, .by = v)` reports three kinds of row:
#
#   * the overall rows the result already reports, keyed by the group label
#     "overall" and still carrying all three effect measures for a binary
#     outcome;
#   * one row per effect measure per stratum of the modifier, keyed by a
#     "var = value" group label, whose counterfactual risks are the tilt-weighted
#     g-computation means taken within that stratum;
#   * the effect-modification rows, each non-reference stratum's effect minus the
#     reference stratum's, keyed by a "var = value vs var = value" group label.
#
# The stratum rows and the effect-modification rows are reported on the rd and
# log(rr) scales for a binary outcome and on the diff scale for a continuous
# one. The odds ratio is noncollapsible, so it is reported for the overall rows
# alone.
#
# The oracle for every point estimate here is the same hand computation: set the
# exposure column of the fitted outcome model's data to each of its values in
# turn, leaving the modifier and the covariates as each unit has them, average
# the predictions over the units of one stratum weighted by the estimand tilt,
# and read the effect measures off the pair of averages. `by_stratum_effects()`
# does that and nothing else, and the first test in the file anchors it against
# the ungrouped estimates a fit already reports, so a disagreement below is a
# disagreement about strata rather than about the oracle.

# ---- data simulators --------------------------------------------------------

# Binary exposure, binary modifier, and both a binomial and a gaussian outcome
# whose effect differs across the strata of the modifier. The modifier is a
# covariate of the propensity score model as well, which is the usual case: a
# variable worth reporting an effect within is usually a variable worth
# adjusting for.
sim_by <- function(seed = 913, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  v <- factor(rbinom(n, 1, 0.5))
  high <- as.numeric(v == "1")
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.5 * high))
  y <- rbinom(
    n,
    1,
    plogis(-0.4 + 0.9 * z + 0.5 * x1 - 0.3 * high + 0.8 * z * high)
  )
  yc <- 1.5 + 0.8 * z + 0.6 * x1 - 0.4 * high + 1.2 * z * high + rnorm(n)
  data.frame(x1, v, z, y, yc)
}

# Heterogeneous-effect simulator with stored potential outcomes, in the style of
# sim_tilt() in test-ipw-mestimation.R. The individual effect y1 - y0 is 1 in the
# reference stratum and 2.5 in the other, so the difference the effect-
# modification row estimates is exactly 1.5 in the population that generated the
# sample.
sim_by_het <- function(seed = 905, n = 8000) {
  withr::local_seed(seed)
  x <- rnorm(n)
  v <- factor(rbinom(n, 1, 0.5))
  high <- as.numeric(v == "1")
  e <- plogis(0.8 * x - 0.5 * high)
  z <- rbinom(n, 1, e)
  y0 <- x - 0.5 * high + rnorm(n)
  y1 <- y0 + 1 + 1.5 * high
  y <- z * y1 + (1 - z) * y0
  data.frame(x, v, e, z, y0, y1, y)
}

# Continuous exposure, for the guard that refuses `.by` there. The modifier is a
# covariate of the propensity score model so that it is available to be named at
# all and the refusal is about the exposure rather than about a variable that
# could not be found.
sim_by_continuous <- function(seed = 924, n = 300) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  v <- factor(rbinom(n, 1, 0.5))
  a <- 0.5 + 0.8 * x1 - 0.4 * (v == "1") + rnorm(n)
  yc <- 1 + 0.6 * a + 0.5 * x1 + rnorm(n)
  data.frame(x1, v, a, yc)
}

# ---- model fitting ----------------------------------------------------------

# Fit the propensity score model and a weighted outcome model. The weights stay
# a psw object so the estimand survives into the outcome model frame and is
# detected there. `outcome_rhs` and `ps_rhs` carry the right-hand sides the guard
# tests need, and the binomial outcome model tightens its IRLS tolerance so its
# coefficients sit at the weighted maximum likelihood estimate well below the
# tolerance the point-estimate comparisons use.
fit_by_models <- function(
  dat,
  estimand = "ate",
  outcome_family = "binomial",
  outcome_rhs = "z * v",
  ps_rhs = c("x1", "v")
) {
  ps_mod <- glm(
    stats::reformulate(ps_rhs, response = "z"),
    data = dat,
    family = binomial()
  )
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    switch(estimand, ate = wt_ate(ps_mod), att = wt_att(ps_mod))
  )

  outcome_var <- if (outcome_family == "binomial") "y" else "yc"
  fmla <- stats::reformulate(outcome_rhs, response = outcome_var)
  outcome_mod <- if (outcome_family == "binomial") {
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

  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# ---- point-estimate oracle --------------------------------------------------

# Tilt-weighted g-computation within one stratum. `keep` is the stratum
# indicator; the counterfactual predictions come from setting the exposure column
# to each of its values across the whole sample, so a unit's modifier and
# covariates are the ones it has, and the averaging is what restricts the result
# to a stratum. `tilt` is the estimand tilt h(e) evaluated at the fitted
# propensity scores, omitted for ate, whose tilt is the constant one.
by_stratum_effects <- function(
  outcome_mod,
  data,
  keep,
  tilt = NULL,
  exposure_name = "z"
) {
  d1 <- data
  d0 <- data
  d1[[exposure_name]] <- 1
  d0[[exposure_name]] <- 0
  m1 <- predict(outcome_mod, newdata = d1, type = "response")
  m0 <- predict(outcome_mod, newdata = d0, type = "response")

  h <- if (is.null(tilt)) rep(1, nrow(data)) else tilt
  w <- h * as.numeric(keep)
  mu1 <- weighted.mean(m1, w)
  mu0 <- weighted.mean(m0, w)

  linear <- !inherits(outcome_mod, "glm") ||
    stats::family(outcome_mod)$family == "gaussian"
  if (linear) {
    return(c(diff = mu1 - mu0))
  }

  c(
    rd = mu1 - mu0,
    "log(rr)" = log(mu1) - log(mu0),
    "log(or)" = qlogis(mu1) - qlogis(mu0)
  )
}

# The fitted propensity scores an estimand tilts by. Only att is needed here;
# ate tilts by nothing.
by_fitted_ps <- function(ps_mod) {
  as.double(predict(ps_mod, type = "response"))
}

# One value out of the estimates table, addressed by the pair of columns that
# names a row. A request that matches no row returns an empty vector, which the
# comparison that reads it reports as the mismatch it is.
by_estimate <- function(estimates, effect, group) {
  estimates$estimate[estimates$effect == effect & estimates$group == group]
}

by_std_err <- function(estimates, effect, group) {
  estimates$std.err[estimates$effect == effect & estimates$group == group]
}

# ---- the oracle, anchored against the ungrouped estimates -------------------

test_that("the stratum oracle reproduces the ungrouped estimates over the whole sample", {
  skip_if_not_installed("deli")
  dat <- sim_by()
  whole <- rep(TRUE, nrow(dat))

  # A stratum that holds every unit is the whole sample, so the oracle has to
  # return what the fit already reports. This is what makes the stratum
  # comparisons below tests of the strata rather than of the hand computation.
  mods <- fit_by_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  ref <- by_stratum_effects(mods$outcome_mod, dat, whole)
  expect_equal(
    res$estimates$estimate,
    unname(ref[res$estimates$effect]),
    tolerance = 1e-8
  )

  # The same with a tilt, so that the att stratum comparisons rest on an anchored
  # oracle too.
  att <- fit_by_models(dat, estimand = "att")
  res_att <- ipw(att$ps_mod, att$outcome_mod)
  ref_att <- by_stratum_effects(
    att$outcome_mod,
    dat,
    whole,
    tilt = by_fitted_ps(att$ps_mod)
  )
  expect_equal(
    res_att$estimates$estimate,
    unname(ref_att[res_att$estimates$effect]),
    tolerance = 1e-8
  )
})

test_that("an ungrouped binary fit names no subgroups", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  # The frame a fit reports without `.by` carries no group column at all, rather
  # than one repeating a single value down the table, which is the shape `.by`
  # has to leave alone.
  expect_identical(res$estimates$effect, c("rd", "log(rr)", "log(or)"))
  expect_null(res$estimates[["group"]])
  expect_identical(names(coef(res)), c("rd", "log(rr)", "log(or)"))
  expect_false("group" %in% names(as.data.frame(res)))
})

# ---- the shape of a grouped result ------------------------------------------

test_that("a .by fit orders its rows overall, then by stratum, then by stratum contrast", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  # The group column sits immediately after effect, where a categorical exposure
  # puts its contrast column, and the rest of the frame keeps its order.
  expect_identical(
    names(res$estimates),
    c(
      "effect",
      "group",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_identical(
    res$estimates$effect,
    c(
      "rd",
      "log(rr)",
      "log(or)",
      "rd",
      "log(rr)",
      "rd",
      "log(rr)",
      "rd",
      "log(rr)"
    )
  )
  expect_identical(
    res$estimates$group,
    c(
      rep("overall", 3),
      rep("v = 0", 2),
      rep("v = 1", 2),
      rep("v = 1 vs v = 0", 2)
    )
  )
  expect_identical(res$estimand, "ate")
})

test_that("a .by fit reports no odds ratio outside its overall rows", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  # An odds ratio is noncollapsible, so the package reports one for the sample as
  # a whole and reports effect modification on the collapsible scales alone.
  odds <- res$estimates$effect == "log(or)"
  expect_identical(sum(odds), 1L)
  expect_identical(res$estimates$group[odds], "overall")
})

test_that("a .by fit leaves the overall rows it already reported unchanged", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  baseline <- ipw(mods$ps_mod, mods$outcome_mod)
  grouped <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  overall <- grouped$estimates[grouped$estimates$group == "overall", ]
  expect_identical(overall$effect, baseline$estimates$effect)
  expect_equal(
    overall$estimate,
    baseline$estimates$estimate,
    tolerance = 1e-8
  )

  # The stratum parameters and the stratum contrasts are stacked after the
  # overall ones and nothing already in the system depends on them, so the block
  # of the sandwich the overall rows read is the block the ungrouped fit read.
  expect_equal(overall$std.err, baseline$estimates$std.err, tolerance = 1e-6)
})

test_that(".by = NULL reports the frame an ungrouped fit reports", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  baseline <- ipw(mods$ps_mod, mods$outcome_mod)
  explicit <- ipw(mods$ps_mod, mods$outcome_mod, .by = NULL)

  # The default has to be the absence of the feature rather than a grouping over
  # one stratum, so the whole frame is compared, the covariance it carries
  # included.
  expect_identical(explicit$estimates, baseline$estimates)
})

# ---- point estimates --------------------------------------------------------

test_that("a .by ate fit reports the stratum-specific risk differences and log risk ratios", {
  skip_if_not_installed("deli")
  dat <- sim_by()
  mods <- fit_by_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  for (level in levels(dat$v)) {
    ref <- by_stratum_effects(mods$outcome_mod, dat, dat$v == level)
    group <- paste0("v = ", level)
    expect_equal(
      by_estimate(res$estimates, "rd", group),
      ref[["rd"]],
      tolerance = 1e-8,
      label = paste("rd", group)
    )
    expect_equal(
      by_estimate(res$estimates, "log(rr)", group),
      ref[["log(rr)"]],
      tolerance = 1e-8,
      label = paste("log(rr)", group)
    )
  }
})

test_that("a .by ate fit contrasts each stratum against the reference stratum", {
  skip_if_not_installed("deli")
  dat <- sim_by()
  mods <- fit_by_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  ref0 <- by_stratum_effects(mods$outcome_mod, dat, dat$v == "0")
  ref1 <- by_stratum_effects(mods$outcome_mod, dat, dat$v == "1")
  group <- "v = 1 vs v = 0"

  expect_equal(
    by_estimate(res$estimates, "rd", group),
    ref1[["rd"]] - ref0[["rd"]],
    tolerance = 1e-8
  )
  expect_equal(
    by_estimate(res$estimates, "log(rr)", group),
    ref1[["log(rr)"]] - ref0[["log(rr)"]],
    tolerance = 1e-8
  )
})

test_that("a .by att fit weights each stratum mean by the tilt within the stratum", {
  skip_if_not_installed("deli")
  dat <- sim_by()
  mods <- fit_by_models(dat, estimand = "att")
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)
  tilt <- by_fitted_ps(mods$ps_mod)

  expect_identical(res$estimand, "att")

  # The att stratum mean is the counterfactual prediction averaged over the units
  # of the stratum with weight e(x), which is the estimand tilt multiplied by the
  # stratum indicator.
  refs <- lapply(levels(dat$v), function(level) {
    by_stratum_effects(mods$outcome_mod, dat, dat$v == level, tilt = tilt)
  })
  names(refs) <- levels(dat$v)

  for (level in levels(dat$v)) {
    group <- paste0("v = ", level)
    expect_equal(
      by_estimate(res$estimates, "rd", group),
      refs[[level]][["rd"]],
      tolerance = 1e-8,
      label = paste("rd", group)
    )
    expect_equal(
      by_estimate(res$estimates, "log(rr)", group),
      refs[[level]][["log(rr)"]],
      tolerance = 1e-8,
      label = paste("log(rr)", group)
    )
  }

  expect_equal(
    by_estimate(res$estimates, "rd", "v = 1 vs v = 0"),
    refs[["1"]][["rd"]] - refs[["0"]][["rd"]],
    tolerance = 1e-8
  )
})

test_that("a .by fit on a continuous outcome reports one diff row per stratum", {
  skip_if_not_installed("deli")
  dat <- sim_by()
  mods <- fit_by_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  expect_identical(res$estimates$effect, rep("diff", 4))
  expect_identical(
    res$estimates$group,
    c("overall", "v = 0", "v = 1", "v = 1 vs v = 0")
  )

  ref0 <- by_stratum_effects(mods$outcome_mod, dat, dat$v == "0")
  ref1 <- by_stratum_effects(mods$outcome_mod, dat, dat$v == "1")
  expect_equal(
    by_estimate(res$estimates, "diff", "v = 0"),
    ref0[["diff"]],
    tolerance = 1e-8
  )
  expect_equal(
    by_estimate(res$estimates, "diff", "v = 1"),
    ref1[["diff"]],
    tolerance = 1e-8
  )

  # A linear outcome model of the form y ~ z * v predicts a within-stratum
  # difference of beta_z in the reference stratum and beta_z + beta_zv in the
  # other, whatever the covariate values are, so the stratum contrast is the
  # fitted interaction coefficient exactly.
  expect_equal(
    by_estimate(res$estimates, "diff", "v = 1 vs v = 0"),
    unname(coef(mods$outcome_mod)[["z:v1"]]),
    tolerance = 1e-8
  )
})

test_that("a .by contrast recovers the simulated difference in effects", {
  skip_if_not_installed("deli")
  dat <- sim_by_het()
  ps_mod <- glm(z ~ x + v, data = dat, family = binomial())
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_mod))
  outcome_mod <- lm(y ~ z * v + x, data = dat, weights = wts)
  res <- ipw(ps_mod, outcome_mod, .by = v)

  group <- "v = 1 vs v = 0"
  got <- by_estimate(res$estimates, "diff", group)
  se <- by_std_err(res$estimates, "diff", group)

  # The individual effects are stored by the simulator, so the quantity the
  # stratum contrast estimates is available exactly rather than as a target the
  # fit is compared against loosely.
  tau <- dat$y1 - dat$y0
  oracle <- mean(tau[dat$v == "1"]) - mean(tau[dat$v == "0"])
  expect_equal(oracle, 1.5, tolerance = 1e-12)

  # The outcome model is linear, so the contrast the fit reports is the fitted
  # interaction coefficient, which is 1.4934 on this fixture. The comparison
  # against the simulation is therefore a comparison of a known number against a
  # known number, and the tolerances only have to cover sampling noise.
  expect_equal(got, unname(coef(outcome_mod)[["z:v1"]]), tolerance = 1e-8)
  expect_lt(abs(got - oracle), 0.1)
  expect_lt(abs(got - oracle), 3 * se)
})

# ---- standard errors and labels ---------------------------------------------

test_that("a .by fit reports a usable standard error for every row", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)
  est <- res$estimates

  expect_true(all(is.finite(est$std.err)))
  expect_true(all(est$std.err > 0))
  expect_true(all(is.finite(est$ci.lower)))
  expect_true(all(is.finite(est$ci.upper)))
  expect_true(all(est$ci.lower < est$ci.upper))

  # The stratum rows and the stratum contrasts are stacked parameters, so the
  # sandwich covers them and the standard errors the table reports are its
  # diagonal.
  expect_equal(sqrt(diag(vcov(res))), est$std.err, ignore_attr = TRUE)
})

test_that("a .by fit labels its coefficients, covariance, and printed rows alike", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  # A binary exposure names no contrasts, so a row's identity is the effect
  # measure and the subgroup, in that order, which is what the causalgenerics
  # labelling contract pastes.
  labels <- paste(res$estimates$effect, res$estimates$group)
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(names(coef(res)), labels)
  expect_identical(dimnames(vcov(res)), list(labels, labels))
  expect_identical(rownames(confint(res)), labels)
  expect_identical(
    dimnames(attr(res$estimates, "ipw_vcov", exact = TRUE)),
    list(labels, labels)
  )

  printed <- capture.output(print(res))
  expect_true(all(vapply(
    labels,
    function(label) any(startsWith(printed, label)),
    logical(1)
  )))
})

test_that("the conditional reading of a .by fit is the ungrouped fit's", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  plain <- causalgenerics::as_conditional(ipw(mods$ps_mod, mods$outcome_mod))
  grouped <- causalgenerics::as_conditional(
    ipw(mods$ps_mod, mods$outcome_mod, .by = v)
  )

  # `.by` adds stacked parameters after the outcome block and nothing in that
  # block depends on them, so the coefficient surface and the covariance the
  # conditional reading reports are the ones the ungrouped fit reports.
  expect_identical(coef(grouped), coef(plain))
  expect_equal(vcov(grouped), vcov(plain), tolerance = 1e-8)
})

# ---- the tidier surfaces ----------------------------------------------------

test_that("as.data.frame() heads its subgroup column after the term column", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  tbl <- as.data.frame(res)
  expect_identical(
    names(tbl),
    c("term", "group", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(tbl$term, res$estimates$effect)
  expect_identical(tbl$group, res$estimates$group)

  # An interval is appended to the table rather than rearranging it, so the
  # subgroup column keeps its place.
  bounded <- as.data.frame(res, conf.int = TRUE)
  expect_identical(
    names(bounded),
    c(
      "term",
      "group",
      "estimate",
      "std.error",
      "statistic",
      "p.value",
      "conf.low",
      "conf.high"
    )
  )
})

test_that("tidy() carries the subgroup column", {
  skip_if_not_installed("deli")
  mods <- fit_by_models(sim_by())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  tidied <- tidy(res)
  expect_s3_class(tidied, "tbl_df")
  expect_identical(
    names(tidied),
    c("term", "group", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(tidied$group, res$estimates$group)
})

# ---- refusals ---------------------------------------------------------------
#
# Every refusal below is raised before the outcome model is inspected for an
# exposure-by-modifier term, so a fixture that would also draw that warning
# raises the error alone. `expect_no_warning()` pins that order: a request the
# function will not answer is refused rather than diagnosed first.

test_that(".by refuses linearization standard errors", {
  dat <- sim_by()
  mods <- fit_by_models(dat, outcome_rhs = "z")

  expect_no_warning(expect_error(
    ipw(
      mods$ps_mod,
      mods$outcome_mod,
      se_method = "linearization",
      .by = v
    ),
    class = "propensity_ipw_by_method_error"
  ))
  expect_propensity_error(
    ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization", .by = v)
  )
})

test_that(".by refuses a continuous exposure", {
  skip_if_not_installed("deli")
  dat <- sim_by_continuous()
  ps_mod <- lm(a ~ x1 + v, data = dat)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(fitted(ps_mod)),
      dat$a,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )
  msm <- lm(yc ~ a, data = dat, weights = wts)

  expect_no_warning(expect_error(
    ipw(ps_mod, msm, .by = v),
    class = "propensity_ipw_by_exposure_error"
  ))
  expect_propensity_error(ipw(ps_mod, msm, .by = v))
})

test_that(".by refuses a modifier with missing values", {
  skip_if_not_installed("deli")
  dat <- sim_by()
  dat$w <- factor(rep(c("lo", "hi"), length.out = nrow(dat)))
  dat$w[1] <- NA
  mods <- fit_by_models(dat)

  # The modifier is read from `.data` when it is supplied, so a column neither
  # model was fit with can still name the strata, and a missing value in it names
  # a stratum that is not one.
  expect_no_warning(expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat, .by = w),
    class = "propensity_ipw_by_missing_error"
  ))
  expect_propensity_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat, .by = w)
  )
})

test_that(".by refuses a numeric modifier", {
  skip_if_not_installed("deli")
  dat <- sim_by()
  # The outcome model carries the exposure-by-modifier term, so the only thing
  # wrong with the request is the type of the modifier.
  mods <- fit_by_models(dat, outcome_rhs = "z * x1")

  expect_no_warning(expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .by = x1),
    class = "propensity_ipw_by_type_error"
  ))
  expect_propensity_error(ipw(mods$ps_mod, mods$outcome_mod, .by = x1))
})

test_that(".by refuses a stratum missing an exposure level", {
  skip_if_not_installed("deli")
  dat <- sim_by()
  # Every unit of stratum "b" is exposed, so that stratum holds no unexposed arm
  # to contrast against and its effect is not estimable from the data. The
  # propensity score model is fit without the modifier, since a covariate that
  # separates the exposure would leave that model with no finite fit and report
  # something else.
  dat$g <- factor(ifelse(dat$z == 1 & dat$x1 > 0, "b", "a"))
  mods <- fit_by_models(dat, ps_rhs = "x1", outcome_rhs = "z + g")

  expect_no_warning(expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat, .by = g),
    class = "propensity_ipw_by_levels_error"
  ))
  # The remedy is a coarser modifier, and the message has to say so: the fit
  # cannot be rescued by refitting either model.
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat, .by = g),
    regexp = "coars"
  )
  expect_propensity_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat, .by = g)
  )
})

test_that(".by warns when the outcome model has no exposure-by-modifier term", {
  skip_if_not_installed("deli")
  dat <- sim_by()
  mods <- fit_by_models(dat, outcome_rhs = "z + v")

  # An outcome model with no interaction forces the same effect on every
  # stratum, which is almost never what a request for stratum-specific effects
  # means. It is a warning rather than an error: the model as specified is what
  # the g-computation reports, and the fit still returns.
  res <- expect_warning(
    ipw(mods$ps_mod, mods$outcome_mod, .by = v),
    class = "propensity_ipw_by_interaction_warning"
  )
  expect_s3_class(res, "ipw")
  expect_true("group" %in% names(res$estimates))

  expect_propensity_warning(ipw(mods$ps_mod, mods$outcome_mod, .by = v))
})
