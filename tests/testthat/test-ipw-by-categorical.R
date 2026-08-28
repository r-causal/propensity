# Effect modification through `.by` for a categorical exposure fit through
# nnet::multinom. The binary contract in test-ipw-by.R carries over with one
# addition: a categorical exposure already reports one contrast per
# non-reference level, so every reported row is now named by three things, the
# effect measure, the contrast, and the subgroup, and the reported rows are the
# crossing of the contrasts with the groups.
#
# The rows come in the same three blocks a binary fit reports:
#
#   * the overall rows, unchanged from the ungrouped categorical fit, all three
#     effect measures for each contrast, under the group label "overall";
#   * one block per stratum of the modifier, each contrast on the rd and log(rr)
#     scales, under a "var = value" label;
#   * one block per non-reference stratum against the reference stratum, the
#     same measures, under a "var = value vs var = value" label.
#
# Within every block the ordering is the one an ungrouped categorical fit
# already uses: contrast-major, effect-minor, so a contrast's measures sit
# together. The blocks themselves run group-major, so the full key runs
# group, then contrast, then effect.
#
# The oracle is the binary one generalized: predict each unit's outcome at every
# exposure level, average the predictions over the units of one stratum weighted
# by the estimand tilt, and form each non-reference level's contrast against the
# reference level from that stratum's means. `by_cat_stratum_effects()` does
# that, and the first test anchors it against the estimates an ungrouped
# categorical fit already reports.

# ---- data simulator ---------------------------------------------------------

# Three-level exposure, binary modifier, and both a binomial and a gaussian
# outcome whose exposure effect differs across the strata of the modifier. The
# modifier is a covariate of the multinomial propensity score model as well.
# `v_levels` declares the modifier's levels, so a caller can declare one no unit
# carries.
sim_by_categorical <- function(seed = 617, n = 700, v_levels = c("0", "1")) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  v <- factor(rbinom(n, 1, 0.5), levels = v_levels)
  high <- as.numeric(v == "1")

  eta_b <- -0.2 + 0.5 * x1 + 0.4 * high
  eta_c <- 0.1 - 0.4 * x1 - 0.3 * high
  denom <- 1 + exp(eta_b) + exp(eta_c)
  pa <- 1 / denom
  pb <- exp(eta_b) / denom
  u <- runif(n)
  a <- factor(
    ifelse(u < pa, "a", ifelse(u < pa + pb, "b", "c")),
    levels = c("a", "b", "c")
  )

  y <- rbinom(
    n,
    1,
    plogis(
      -0.3 +
        0.4 * (a == "b") +
        0.8 * (a == "c") +
        0.5 * x1 -
        0.2 * high +
        0.7 * (a == "b") * high -
        0.5 * (a == "c") * high
    )
  )
  yc <- 0.5 +
    0.4 * (a == "b") +
    0.9 * (a == "c") +
    0.6 * x1 -
    0.3 * high +
    0.8 * (a == "b") * high +
    rnorm(n)

  data.frame(x1, v, a, y, yc)
}

# ---- model fitting ----------------------------------------------------------

# The outcome model every point-estimate comparison is made under, named once so
# the anchor and the comparisons it anchors cannot drift apart.
#
# The covariate is what makes those comparisons discriminating, for the reason
# the binary file records. Under `y ~ a * v` the counterfactual predictions are
# constant within a stratum, so the tilt-weighted stratum mean, the untilted
# stratum mean, and the mean over the stratum's observed-arm units agree to
# about 7e-15, and an implementation that dropped the tilt or averaged the wrong
# population would match the oracle anyway. With `x1` in the model those
# candidates separate by 2e-3 to 7e-2, orders above the tolerance the
# comparisons carry.
by_cat_outcome_rhs <- "a * v + x1"

# The tolerance every point-estimate comparison in this file carries, the anchor
# included. It is looser than the 1e-8 the binary file uses because the
# multinomial stack solves less tightly: measured against the oracle on the
# ungrouped fit, the reported estimates agree to between 1.7e-8 and 8.6e-8
# relative, which is where this fixture's solve stops rather than a disagreement
# about what is reported.
#
# It is still four orders tighter than what the comparisons have to tell apart.
# The closest wrong candidate, the untilted stratum mean under att, differs from
# the correct one by 2e-3 on estimates of order 0.02 to 0.3, which is at least
# 8e-2 relative.
by_cat_tolerance <- 1e-6

# Fit the multinomial propensity score model, categorical weights for the
# estimand, and a weighted outcome model. The multinomial fit tightens its
# convergence so the seeded init sits at the score root well below the
# point-estimate tolerance, matching the pattern in test-ipw-categorical.R. The
# weights stay a psw object so the estimand, and the focal level an att fit
# carries, survive into the outcome model frame and are detected there.
fit_by_cat_models <- function(
  dat,
  estimand = "ate",
  focal_level = "b",
  outcome_family = "binomial",
  outcome_rhs = by_cat_outcome_rhs,
  ps_rhs = c("x1", "v")
) {
  ps_mod <- nnet::multinom(
    stats::reformulate(ps_rhs, response = "a"),
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_named <- unname(predict(ps_mod, type = "probs"))
  colnames(ps_named) <- ps_mod$lev

  args <- list(ps_named, dat$a, exposure_type = "categorical")
  if (identical(estimand, "att")) {
    args$.focal_level <- focal_level
  }
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    do.call(switch(estimand, ate = wt_ate, att = wt_att), args)
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

  list(
    ps_mod = ps_mod,
    outcome_mod = outcome_mod,
    wts = wts,
    ps = ps_named,
    lev = ps_mod$lev
  )
}

# ---- point-estimate oracle --------------------------------------------------

# The tilt a categorical estimand standardizes to, or NULL for ate, whose tilt
# is the constant one. Only att is needed here, whose tilt is the focal level's
# propensity score column.
by_cat_tilt <- function(mods, estimand, focal_level = "b") {
  if (identical(estimand, "ate")) {
    return(NULL)
  }
  mods$ps[, match(focal_level, mods$lev)]
}

by_cat_transform <- function(form, mu_hi, mu_lo) {
  switch(
    form,
    rd = mu_hi - mu_lo,
    diff = mu_hi - mu_lo,
    "log(rr)" = log(mu_hi) - log(mu_lo),
    "log(or)" = qlogis(mu_hi) - qlogis(mu_lo)
  )
}

# The counterfactual mean at every exposure level, taken over one stratum. The
# predictions come from setting the exposure column to each level across the
# whole sample, so a unit's modifier and covariates are the ones it has, and the
# averaging is what restricts the result to a stratum. `keep` is the stratum
# indicator and `tilt` the estimand tilt, omitted for ate.
by_cat_stratum_means <- function(
  outcome_mod,
  data,
  keep,
  lev,
  tilt = NULL,
  exposure_name = "a"
) {
  h <- if (is.null(tilt)) rep(1, nrow(data)) else tilt
  w <- h * as.numeric(keep)

  vapply(
    lev,
    function(level) {
      d <- data
      d[[exposure_name]] <- factor(level, levels = lev)
      stats::weighted.mean(
        predict(outcome_mod, newdata = d, type = "response"),
        w
      )
    },
    numeric(1)
  )
}

# One stratum's contrasts, each non-reference level against the reference level,
# on the measures the rows being compared report. Returned as a frame keyed by
# effect and contrast so it can be matched to the estimates table whatever order
# either is in.
by_cat_stratum_effects <- function(
  outcome_mod,
  data,
  keep,
  lev,
  tilt = NULL,
  forms = c("rd", "log(rr)"),
  exposure_name = "a"
) {
  mu <- by_cat_stratum_means(
    outcome_mod,
    data,
    keep,
    lev,
    tilt = tilt,
    exposure_name = exposure_name
  )
  ref <- mu[[1]]

  # The counterfactual mean at each level is a reported row of its own, keyed by
  # the level rather than by a pair of them.
  rows <- list(data.frame(
    effect = "mean",
    contrast = lev,
    estimate = unname(mu),
    stringsAsFactors = FALSE
  ))
  for (j in seq_along(lev)[-1]) {
    for (form in forms) {
      rows[[length(rows) + 1L]] <- data.frame(
        effect = form,
        contrast = paste(lev[[j]], "vs", lev[[1]]),
        estimate = by_cat_transform(form, mu[[j]], ref),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

# ---- reading the estimates table --------------------------------------------

# The rows of a grouped estimates table belonging to one group, in table order.
by_cat_group_rows <- function(estimates, group) {
  estimates[estimates$group == group, , drop = FALSE]
}

# Compare the rows of one group against a reference frame keyed by effect and
# contrast, matching on the key so row order is not what is being tested here.
expect_by_cat_match <- function(
  estimates,
  group,
  reference,
  tolerance = by_cat_tolerance
) {
  rows <- by_cat_group_rows(estimates, group)
  key_got <- paste(rows$effect, rows$contrast)
  key_ref <- paste(reference$effect, reference$contrast)
  got <- rows$estimate[match(key_ref, key_got)]
  expect_equal(got, reference$estimate, tolerance = tolerance, label = group)
}

# The contrast rows of an oracle frame. A stratum reports the counterfactual
# means as well as the contrasts over them, but the effect-modification rows
# compare two strata's contrasts and have no mean of their own.
by_cat_contrast_rows <- function(reference) {
  reference[reference$effect != "mean", , drop = FALSE]
}

# The labels every surface of a grouped categorical result keys its rows by.
by_cat_labels <- function(estimates) {
  paste(estimates$effect, estimates$contrast, estimates$group)
}

# ---- the oracle, anchored against the ungrouped estimates -------------------

test_that("the categorical stratum oracle reproduces the ungrouped estimates over the whole sample", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_by_categorical()
  whole <- rep(TRUE, nrow(dat))
  all_forms <- c("rd", "log(rr)", "log(or)")

  # A stratum that holds every unit is the whole sample, so the oracle has to
  # return what an ungrouped fit already reports. That is what makes the stratum
  # comparisons below tests of the strata rather than of the hand computation,
  # and it is done under the covariate-bearing outcome model those comparisons
  # use so the anchoring argument covers the model they are made under.
  mods <- fit_by_cat_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  ref <- by_cat_stratum_effects(
    mods$outcome_mod,
    dat,
    whole,
    mods$lev,
    forms = all_forms
  )
  key_got <- paste(res$estimates$effect, res$estimates$contrast)
  key_ref <- paste(ref$effect, ref$contrast)
  expect_equal(
    res$estimates$estimate[match(key_ref, key_got)],
    ref$estimate,
    tolerance = by_cat_tolerance
  )

  # The same with a tilt, so the att stratum comparisons rest on an anchored
  # oracle too.
  att <- fit_by_cat_models(dat, estimand = "att")
  res_att <- ipw(att$ps_mod, att$outcome_mod)
  ref_att <- by_cat_stratum_effects(
    att$outcome_mod,
    dat,
    whole,
    att$lev,
    tilt = by_cat_tilt(att, "att"),
    forms = all_forms
  )
  key_att <- paste(res_att$estimates$effect, res_att$estimates$contrast)
  expect_equal(
    res_att$estimates$estimate[match(key_ref, key_att)],
    ref_att$estimate,
    tolerance = by_cat_tolerance
  )
})

test_that("an ungrouped categorical fit names no subgroups", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  # The frame a categorical fit reports without `.by` names contrasts and no
  # subgroups, which is the shape `.by` has to leave alone.
  expect_identical(
    res$estimates$effect,
    c(rep("mean", 3), rep(c("rd", "log(rr)", "log(or)"), times = 2))
  )
  expect_identical(
    res$estimates$contrast,
    c(c("a", "b", "c"), rep(c("b vs a", "c vs a"), each = 3))
  )
  expect_null(res$estimates[["group"]])
  expect_identical(
    names(coef(res)),
    paste(
      res$estimates$effect,
      res$estimates$contrast
    )
  )
  expect_false("group" %in% names(as.data.frame(res)))
})

test_that(".by = NULL reports the frame an ungrouped categorical fit reports", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  baseline <- ipw(mods$ps_mod, mods$outcome_mod)
  explicit <- ipw(mods$ps_mod, mods$outcome_mod, .by = NULL)

  # The default has to be the absence of the feature rather than a grouping over
  # one stratum, so the whole frame is compared, the covariance it carries
  # included.
  expect_identical(explicit$estimates, baseline$estimates)
})

# ---- the shape of a grouped categorical result -------------------------------

test_that("a .by categorical fit crosses its contrasts with its groups", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)
  est <- res$estimates

  # A row is named by the effect measure, the contrast, and the subgroup, in
  # that order, which is the order the labelling contract pastes them and so the
  # order the columns sit in: `group` goes after `contrast`, not before it.
  expect_identical(
    names(est),
    c(
      "effect",
      "contrast",
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

  # The whole-sample means lead, then the whole-sample contrasts, then the
  # stratum means, then the stratum contrasts. Blocks run group-major; within a
  # contrast block the ordering is the one an ungrouped categorical fit uses,
  # contrast-major and effect-minor, so a contrast's measures sit together.
  levs <- c("a", "b", "c")
  contrasts <- c("b vs a", "c vs a")
  overall <- c("rd", "log(rr)", "log(or)")
  stratum <- c("rd", "log(rr)")
  strata <- c("v = 0", "v = 1")
  groups <- c("v = 0", "v = 1", "v = 1 vs v = 0")

  expect_identical(
    est$effect,
    c(
      rep("mean", length(levs)),
      rep(overall, times = length(contrasts)),
      rep("mean", length(levs) * length(strata)),
      rep(rep(stratum, times = length(contrasts)), times = length(groups))
    )
  )
  expect_identical(
    est$contrast,
    c(
      levs,
      rep(contrasts, each = length(overall)),
      rep(levs, times = length(strata)),
      rep(rep(contrasts, each = length(stratum)), times = length(groups))
    )
  )
  expect_identical(
    est$group,
    c(
      rep(
        "overall",
        length(levs) + length(overall) * length(contrasts)
      ),
      rep(strata, each = length(levs)),
      rep(groups, each = length(stratum) * length(contrasts))
    )
  )
  expect_identical(nrow(est), 27L)
  expect_identical(res$estimand, "ate")
})

test_that("a .by categorical fit reports no odds ratio outside its overall rows", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  # An odds ratio is noncollapsible, so the package reports one per contrast for
  # the sample as a whole and reports effect modification on the collapsible
  # scales alone. There is one per contrast and no more.
  odds <- res$estimates$effect == "log(or)"
  expect_identical(sum(odds), 2L)
  expect_identical(res$estimates$group[odds], rep("overall", 2))
  expect_identical(res$estimates$contrast[odds], c("b vs a", "c vs a"))
})

test_that("a .by categorical fit leaves the overall rows it already reported unchanged", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  baseline <- ipw(mods$ps_mod, mods$outcome_mod)
  grouped <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  overall <- by_cat_group_rows(grouped$estimates, "overall")
  expect_identical(overall$effect, baseline$estimates$effect)
  expect_identical(overall$contrast, baseline$estimates$contrast)
  expect_equal(
    overall$estimate,
    baseline$estimates$estimate,
    tolerance = by_cat_tolerance
  )

  # The stratum parameters and the stratum contrasts are stacked after the
  # whole-sample ones and nothing already in the system depends on them, so the
  # block of the sandwich the overall rows read is the block the ungrouped fit
  # read.
  expect_equal(overall$std.err, baseline$estimates$std.err, tolerance = 1e-6)
})

# ---- point estimates ---------------------------------------------------------

test_that("a .by categorical ate fit reports each contrast within each stratum", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_by_categorical()
  mods <- fit_by_cat_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  # The covariate in the outcome model is what makes this a test of the
  # population each mean is taken over: the counterfactual predictions vary
  # within a stratum, so restricting the average to the stratum's units gives an
  # answer that averaging over the whole sample, or over the stratum's
  # observed-arm units, does not.
  for (level in levels(dat$v)) {
    ref <- by_cat_stratum_effects(
      mods$outcome_mod,
      dat,
      dat$v == level,
      mods$lev
    )
    expect_by_cat_match(res$estimates, paste0("v = ", level), ref)
  }
})

test_that("a .by categorical ate fit contrasts each stratum against the reference stratum", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_by_categorical()
  mods <- fit_by_cat_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  ref0 <- by_cat_stratum_effects(mods$outcome_mod, dat, dat$v == "0", mods$lev)
  ref1 <- by_cat_stratum_effects(mods$outcome_mod, dat, dat$v == "1", mods$lev)

  # The stratum contrast is taken per contrast and per scale: the same
  # comparison of exposure levels, differenced between the two strata.
  em <- by_cat_contrast_rows(ref1)
  em$estimate <- em$estimate - by_cat_contrast_rows(ref0)$estimate
  expect_by_cat_match(res$estimates, "v = 1 vs v = 0", em)
})

test_that("a .by categorical att fit weights each stratum mean by the tilt within the stratum", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_by_categorical()
  mods <- fit_by_cat_models(dat, estimand = "att")
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)
  tilt <- by_cat_tilt(mods, "att")

  expect_identical(res$estimand, "att")

  # The att stratum mean averages the counterfactual predictions over the units
  # of the stratum with weight e_b(x), the focal level's propensity score, which
  # is the estimand tilt multiplied by the stratum indicator. The covariate in
  # the outcome model is what separates that from the untilted stratum mean;
  # without it the tilt divides out of every contrast.
  refs <- lapply(levels(dat$v), function(level) {
    by_cat_stratum_effects(
      mods$outcome_mod,
      dat,
      dat$v == level,
      mods$lev,
      tilt = tilt
    )
  })
  names(refs) <- levels(dat$v)

  for (level in levels(dat$v)) {
    expect_by_cat_match(res$estimates, paste0("v = ", level), refs[[level]])
  }

  em <- by_cat_contrast_rows(refs[["1"]])
  em$estimate <- em$estimate - by_cat_contrast_rows(refs[["0"]])$estimate
  expect_by_cat_match(res$estimates, "v = 1 vs v = 0", em)
})

test_that("a .by categorical fit on a continuous outcome reports one diff row per contrast", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_by_categorical()
  mods <- fit_by_cat_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)
  est <- res$estimates

  expect_identical(
    est$effect,
    c(
      rep("mean", 3),
      rep("diff", 2),
      rep("mean", 6),
      rep("diff", 6)
    )
  )
  expect_identical(
    est$contrast,
    c(
      c("a", "b", "c"),
      c("b vs a", "c vs a"),
      rep(c("a", "b", "c"), times = 2),
      rep(c("b vs a", "c vs a"), times = 3)
    )
  )
  expect_identical(
    est$group,
    c(
      rep("overall", 5),
      rep(c("v = 0", "v = 1"), each = 3),
      rep(c("v = 0", "v = 1", "v = 1 vs v = 0"), each = 2)
    )
  )

  # A linear outcome model has each level's counterfactual mean differing from
  # the reference level's by a constant within a stratum, so no comparison made
  # under one can pin the weights the stratum mean averages with. What this
  # pins is the row shape and the scale set; the binary-outcome tests above
  # carry the pin on the averaging population and the tilt.
  for (level in levels(dat$v)) {
    ref <- by_cat_stratum_effects(
      mods$outcome_mod,
      dat,
      dat$v == level,
      mods$lev,
      forms = "diff"
    )
    expect_by_cat_match(est, paste0("v = ", level), ref)
  }
})

# ---- standard errors and labels ----------------------------------------------

test_that("a .by categorical fit reports a usable standard error for every row", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
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

test_that("a .by categorical fit's covariance couples the groups and contrasts it reports", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)
  covariance <- vcov(res)

  # Every reported row is a function of one multinomial propensity score model
  # and one outcome model, fit to one sample and solved as one stacked system,
  # so rows covary across contrasts as well as across groups. A per-stratum
  # refit assembled into a block-diagonal matrix would report exact zeros in
  # these entries while matching its own standard errors everywhere, so the
  # off-diagonal is what tells the two apart.
  #
  # The floor is absolute and far below any covariance shared parameters
  # generate here. What it excludes is the structural zero, not a small number.
  couples <- list(
    c("rd b vs a v = 0", "rd b vs a v = 1"),
    c("rd b vs a v = 0", "rd c vs a v = 0"),
    c("rd b vs a overall", "rd b vs a v = 1"),
    c("rd b vs a v = 1", "rd b vs a v = 1 vs v = 0")
  )
  for (pair in couples) {
    entry <- covariance[pair[[1]], pair[[2]]]
    expect_true(is.finite(entry), label = paste(pair, collapse = " with "))
    expect_gt(abs(entry), 1e-8, label = paste(pair, collapse = " with "))
  }

  expect_equal(covariance, t(covariance), tolerance = 1e-12)
})

test_that("a .by categorical fit labels its rows by effect, contrast, and subgroup", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  labels <- by_cat_labels(res$estimates)
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

test_that("the conditional reading of a .by categorical fit is the ungrouped fit's", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  plain <- causalgenerics::as_conditional(ipw(mods$ps_mod, mods$outcome_mod))
  grouped <- causalgenerics::as_conditional(
    ipw(mods$ps_mod, mods$outcome_mod, .by = v)
  )

  # `.by` adds stacked parameters after the outcome block and nothing in that
  # block depends on them, so the coefficient surface and the covariance the
  # conditional reading reports are the ones the ungrouped fit reports.
  expect_identical(coef(grouped), coef(plain))
  expect_equal(vcov(grouped), vcov(plain), tolerance = by_cat_tolerance)
})

# ---- the tidier surfaces -----------------------------------------------------

test_that("as.data.frame() heads its subgroup column after the contrast column", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  tbl <- as.data.frame(res)
  expect_identical(
    names(tbl),
    c(
      "term",
      "contrast",
      "group",
      "estimate",
      "std.error",
      "statistic",
      "p.value"
    )
  )
  expect_identical(tbl$term, res$estimates$effect)
  expect_identical(tbl$contrast, res$estimates$contrast)
  expect_identical(tbl$group, res$estimates$group)
})

test_that("tidy() carries the contrast and the subgroup columns", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- fit_by_cat_models(sim_by_categorical())
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  tidied <- tidy(res)
  expect_s3_class(tidied, "tbl_df")
  expect_identical(
    names(tidied),
    c(
      "term",
      "contrast",
      "group",
      "estimate",
      "std.error",
      "statistic",
      "p.value"
    )
  )
  expect_identical(tidied$contrast, res$estimates$contrast)
  expect_identical(tidied$group, res$estimates$group)
})

# ---- the modifier ------------------------------------------------------------

test_that("a .by categorical fit drops a modifier level no unit carries", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The modifier declares a third level the simulator never draws. An empty
  # stratum has no mean to standardize and no contrast to report, and the fits
  # themselves drop the level, so reading the modifier the way the models were
  # fit means dropping it here too. That is not worth a warning: nothing about
  # the request is wrong.
  dat <- sim_by_categorical(v_levels = c("0", "1", "2"))
  expect_identical(levels(dat$v), c("0", "1", "2"))
  expect_identical(sum(dat$v == "2"), 0L)

  mods <- fit_by_cat_models(dat)
  res <- expect_no_warning(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat, .by = v)
  )

  expect_identical(
    unique(res$estimates$group),
    c("overall", "v = 0", "v = 1", "v = 1 vs v = 0")
  )
  expect_identical(nrow(res$estimates), 27L)
})

test_that(".by refuses a stratum missing one of a categorical exposure's levels", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_by_categorical()
  # Stratum "b" holds the exposure levels "a" and "b" and no unit of level "c",
  # so the contrast of "c" against the reference level is not estimable there
  # even though the stratum holds two of the three levels. The outcome model is
  # additive in the modifier because interacting it would leave the design rank
  # deficient at the level the stratum does not hold, and the propensity score
  # model omits it because a covariate separating the exposure would leave that
  # model with no finite fit.
  dat$g <- factor(ifelse(dat$a == "c" | dat$x1 > 0, "a", "b"))
  expect_identical(sum(dat$a == "c" & dat$g == "b"), 0L)
  expect_gt(sum(dat$a == "b" & dat$g == "b"), 0L)

  mods <- fit_by_cat_models(dat, outcome_rhs = "a + g + x1", ps_rhs = "x1")

  expect_no_warning(expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .by = g),
    class = "propensity_ipw_by_levels_error"
  ))
  # The remedy is a coarser modifier, and the message has to say so: neither
  # model can be refit to supply a comparison the data do not hold.
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .by = g),
    regexp = "coars"
  )
  expect_propensity_error(ipw(mods$ps_mod, mods$outcome_mod, .by = g))
})

test_that(".by warns when a categorical outcome model has no exposure-by-modifier term", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_by_categorical()
  mods <- fit_by_cat_models(dat, outcome_rhs = "a + v + x1")

  # An outcome model with no interaction forces the same set of contrasts on
  # every stratum, up to the covariate distribution each one has. It is a
  # warning rather than an error: the fit is an honest reading of the model that
  # was supplied, and it still returns.
  expect_warning(
    res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v),
    class = "propensity_ipw_by_interaction_warning"
  )
  expect_s3_class(res, "ipw")
  expect_identical(nrow(res$estimates), 27L)

  expect_propensity_warning(ipw(mods$ps_mod, mods$outcome_mod, .by = v))
})
