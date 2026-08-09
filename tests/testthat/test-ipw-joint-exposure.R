# Reporting for a joint exposure: two treatments intervened on at once.
#
# `causalgenerics::joint_exposure()` crosses two discrete treatments into one
# categorical exposure and records the crossing on the vector. The result is a
# factor, so the whole categorical path already accepts it and reports what it
# reports for any four-level exposure: each cell against the reference cell.
# Those rows are correct and they are not the question a joint intervention
# asks. What the declaration buys is a surface written in terms of the two
# treatments rather than of the four cells:
#
#   * the four counterfactual risks Pr[Y^(a, e)] as rows of their own, under the
#     effect label "mean", each keyed by the cell it belongs to;
#   * the simple effects, each treatment's effect within a fixed level of the
#     other, keyed by the contrast naming the treatment and the group naming the
#     level it is held at, including the comparisons that are not against the
#     reference cell;
#   * the additive and multiplicative interaction, the difference between the
#     two simple effects of the first treatment.
#
# The oracle is the categorical one specialized to the cells: predict each
# unit's outcome with the exposure set to each cell in turn, average over the
# whole sample weighted by the estimand tilt, and build every reported quantity
# from the four means. `joint_cell_means()` does that, and the first tests
# anchor it against the vs-reference estimates the current fit already reports.

# ---- data simulator ---------------------------------------------------------

# Two binary treatments, the second depending on the first and both on a
# covariate, with an outcome carrying a genuine interaction between them. Both a
# binomial and a gaussian outcome are simulated so the reported scales can be
# pinned on each.
sim_joint <- function(seed = 4210, n = 700) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  a <- rbinom(n, 1, plogis(0.3 * x1))
  e <- rbinom(n, 1, plogis(-0.2 + 0.5 * x1 - 0.4 * a))
  y <- rbinom(n, 1, plogis(-0.5 + 0.7 * a + 0.5 * e + 0.6 * x1 + 0.9 * a * e))
  yc <- 1 + 0.6 * a + 0.4 * e + 0.5 * x1 + 0.8 * a * e + rnorm(n)
  # A modifier for the `.by` refusal: a joint exposure and a modifier together
  # are a three-way question this reporting does not answer.
  grp <- factor(ifelse(x1 > 0, "hi", "lo"), levels = c("lo", "hi"))
  data.frame(x1, a, e, y, yc, grp)
}

# The cells, in the order the crossing declares them: the first component varies
# fastest, and the reference cell crosses the two reference levels.
joint_cells <- c("a = 0, e = 0", "a = 1, e = 0", "a = 0, e = 1", "a = 1, e = 1")

# The tolerance every point-estimate comparison in this file carries, the anchor
# included, following the constant the categorical `.by` file sets for the same
# reason: the multinomial stack solves to a floor of its own, and a comparison
# written tighter than that floor fails on arithmetic rather than on contract.
#
# Measured on this fixture, the reported vs-reference estimates agree with the
# oracle to 7.8e-9 relative. Every joint row combines two or four of those
# mean-level gaps, so the interaction rows are expected around 6e-8 relative in
# the green state, which is already past 1e-8.
#
# It stays four orders below what the comparisons have to tell apart. The
# closest wrong candidate, the Hajek weighted mean of the observed outcomes in a
# cell, differs from the correct standardized mean by 5.8e-3 on means of order
# 0.4 to 0.8, which is about 1e-2 relative.
joint_tolerance <- 1e-6

# ---- model fitting ----------------------------------------------------------

# Build the multinomial propensity score model over the cells, categorical
# weights, and a weighted outcome model. `exposure_name` selects the column the
# models read, so the same builder serves the joint exposure and the plain
# factor carrying the same cells. The multinomial fit tightens its convergence
# so the seeded init sits at the score root well below the point-estimate
# tolerance.
#
# The default outcome model carries a covariate, and that is what makes the
# point-estimate comparisons discriminating, for the reason the `.by` files
# record. Under `y ~ joint` alone the counterfactual predictions are constant
# within a cell, so the whole-sample standardized mean, the mean over the units
# observed in that cell, and the Hajek weighted mean of the observed outcomes
# there agree to 4e-15, and an implementation averaging the wrong population
# would match the oracle anyway. With `x1` in the model those candidates
# separate by 5.8e-3 and 4.7e-2, orders above the tolerance the comparisons
# carry.
fit_joint_models <- function(
  dat,
  estimand = "ate",
  focal_level = NULL,
  outcome_family = "binomial",
  exposure_name = "joint",
  outcome_rhs = NULL
) {
  ps_mod <- nnet::multinom(
    stats::reformulate("x1", response = exposure_name),
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_named <- unname(predict(ps_mod, type = "probs"))
  colnames(ps_named) <- ps_mod$lev

  args <- list(ps_named, dat[[exposure_name]], exposure_type = "categorical")
  if (identical(estimand, "att")) {
    args$.focal_level <- focal_level
  }
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    do.call(switch(estimand, ate = wt_ate, att = wt_att), args)
  )

  if (is.null(outcome_rhs)) {
    outcome_rhs <- paste(exposure_name, "+ x1")
  }
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

# The data with the joint exposure declared, and with a plain factor carrying
# exactly the same cells under exactly the same labels. The two differ in the
# declaration and in nothing else, which is what makes the comparison between
# them a test of the declaration.
joint_data <- function(dat = sim_joint()) {
  dat$joint <- causalgenerics::joint_exposure(a = dat$a, e = dat$e)
  dat$plain <- factor(as.character(dat$joint), levels = levels(dat$joint))
  dat
}

# ---- point-estimate oracle --------------------------------------------------

# The four counterfactual means, standardized over the whole sample and weighted
# by the estimand tilt. `tilt` is NULL for ate, whose tilt is the constant one.
#
# The counterfactual column is written as a plain factor rather than assigned
# into the joint column. Setting every row to one cell leaves the other three
# unpopulated, and a joint exposure that loses cells gives up its declaration
# and warns that it did; the plain factor carries the same levels and gives the
# same design.
joint_cell_means <- function(
  outcome_mod,
  data,
  levels,
  tilt = NULL,
  exposure_name = "joint"
) {
  h <- if (is.null(tilt)) rep(1, nrow(data)) else tilt

  vapply(
    levels,
    function(cell) {
      d <- data
      d[[exposure_name]] <- factor(cell, levels = levels)
      stats::weighted.mean(
        predict(outcome_mod, newdata = d, type = "response"),
        h
      )
    },
    numeric(1)
  )
}

joint_transform <- function(form, mu_hi, mu_lo) {
  switch(
    form,
    rd = mu_hi - mu_lo,
    diff = mu_hi - mu_lo,
    "log(rr)" = log(mu_hi) - log(mu_lo),
    "log(or)" = qlogis(mu_hi) - qlogis(mu_lo)
  )
}

# Everything the joint surface reports, built from the four cell means. `mu` is
# indexed in cell order, so `mu[[1]]` is the reference cell (a = 0, e = 0),
# `mu[[2]]` is a = 1 e = 0, `mu[[3]]` is a = 0 e = 1, and `mu[[4]]` is
# a = 1 e = 1.
#
# A simple effect fixes one treatment and contrasts the other. An interaction is
# the difference of the first treatment's two simple effects, which on the log
# scale is the log of the ratio of the two risk ratios.
joint_expected_estimates <- function(mu, forms = c("rd", "log(rr)")) {
  simple <- list(
    list(contrast = "a: 1 vs 0", group = "e = 0", hi = 2L, lo = 1L),
    list(contrast = "a: 1 vs 0", group = "e = 1", hi = 4L, lo = 3L),
    list(contrast = "e: 1 vs 0", group = "a = 0", hi = 3L, lo = 1L),
    list(contrast = "e: 1 vs 0", group = "a = 1", hi = 4L, lo = 2L)
  )

  rows <- list(
    data.frame(
      effect = "mean",
      contrast = joint_cells,
      group = "overall",
      estimate = unname(mu),
      stringsAsFactors = FALSE
    )
  )

  for (cell in simple) {
    rows[[length(rows) + 1L]] <- data.frame(
      effect = forms,
      contrast = cell$contrast,
      group = cell$group,
      estimate = vapply(
        forms,
        function(f) joint_transform(f, mu[[cell$hi]], mu[[cell$lo]]),
        numeric(1)
      ),
      stringsAsFactors = FALSE
    )
  }

  rows[[length(rows) + 1L]] <- data.frame(
    effect = forms,
    contrast = "a: 1 vs 0",
    group = "e = 1 vs e = 0",
    estimate = vapply(
      forms,
      function(f) {
        joint_transform(f, mu[[4]], mu[[3]]) -
          joint_transform(f, mu[[2]], mu[[1]])
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

# The labels every surface of a joint result keys its rows by.
joint_labels <- function(estimates) {
  paste(estimates$effect, estimates$contrast, estimates$group)
}

# ---- the declaration and the oracle, anchored ------------------------------

test_that("joint_exposure() crosses the two treatments into the cells the fixture reports", {
  dat <- joint_data()

  # The fixture rests on the declaration, so what it declares is pinned before
  # anything is read through it.
  expect_true(causalgenerics::is_joint_exposure(dat$joint))
  expect_identical(levels(dat$joint), joint_cells)
  expect_identical(
    causalgenerics::joint_components(dat$joint),
    list(a = c("0", "1"), e = c("0", "1"))
  )
  expect_identical(causalgenerics::joint_reference(dat$joint), joint_cells[[1]])

  # Every cell is populated, which the constructor requires and which the
  # reporting below needs a mean for.
  expect_true(all(table(dat$joint) > 0))

  # The plain factor carries the same cells under the same labels and declares
  # nothing, which is what makes the comparison between the two a test of the
  # declaration rather than of the data.
  expect_false(causalgenerics::is_joint_exposure(dat$plain))
  expect_identical(levels(dat$plain), levels(dat$joint))
  expect_identical(as.character(dat$plain), as.character(dat$joint))
})

test_that("the cell-mean oracle reproduces the vs-reference estimates a categorical fit reports", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  # The plain factor, not the declared exposure. Each cell against the reference
  # cell is what the categorical path reports for a four-level exposure, and it
  # is built from the same four means the joint surface is built from, so it
  # anchors the oracle. Anchoring against the declared exposure would anchor
  # against whatever the joint surface reports, which is the thing under test:
  # it holds only while the declaration falls through to the categorical path,
  # and every match below would go missing the moment it stops doing so.
  mods <- fit_joint_models(dat, exposure_name = "plain")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  mu <- joint_cell_means(
    mods$outcome_mod,
    dat,
    mods$lev,
    exposure_name = "plain"
  )
  expect_identical(names(mu), joint_cells)

  forms <- c("rd", "log(rr)", "log(or)")
  reference <- do.call(
    rbind,
    lapply(seq(2, 4), function(j) {
      data.frame(
        effect = forms,
        contrast = paste(joint_cells[[j]], "vs", joint_cells[[1]]),
        estimate = vapply(
          forms,
          function(f) joint_transform(f, mu[[j]], mu[[1]]),
          numeric(1)
        ),
        stringsAsFactors = FALSE
      )
    })
  )

  key_got <- paste(res$estimates$effect, res$estimates$contrast)
  key_ref <- paste(reference$effect, reference$contrast)
  expect_equal(
    res$estimates$estimate[match(key_ref, key_got)],
    reference$estimate,
    tolerance = joint_tolerance
  )
})

test_that("a plain factor over the same cells reports the vs-reference rows", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat, exposure_name = "plain")
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  est <- res$estimates

  # Without the declaration there is no joint reporting: three effect measures
  # for each of three non-reference cells, no cell means, and no subgroups.
  expect_identical(
    names(est),
    c(
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
  )
  expect_identical(nrow(est), 9L)
  expect_identical(est$effect, rep(c("rd", "log(rr)", "log(or)"), times = 3))
  expect_identical(
    est$contrast,
    rep(paste(joint_cells[-1], "vs", joint_cells[[1]]), each = 3)
  )
  expect_null(est[["group"]])
  expect_false(any(est$effect == "mean"))
})

# ---- the joint surface -------------------------------------------------------

test_that("a joint exposure reports the four counterfactual cell risks", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  est <- res$estimates

  # The rows the declaration adds first: a counterfactual risk per cell, under
  # the effect label "mean" and keyed by the cell in the contrast column. Their
  # absence is what fails here while the reporting is the plain one.
  expect_true("group" %in% names(est))
  means <- est[est$effect == "mean", , drop = FALSE]
  expect_identical(nrow(means), 4L)
  expect_identical(means$contrast, joint_cells)
  expect_identical(means$group, rep("overall", 4))

  mu <- joint_cell_means(mods$outcome_mod, dat, mods$lev)
  expect_equal(means$estimate, unname(mu), tolerance = joint_tolerance)

  # A counterfactual risk is a risk, so it lies where one lies, and it comes
  # out of the same stacked system as everything else.
  expect_true(all(means$estimate > 0 & means$estimate < 1))
  expect_true(all(is.finite(means$std.err)))
  expect_true(all(means$std.err > 0))
})

test_that("a joint exposure orders its rows cells, then simple effects, then interaction", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  est <- res$estimates

  # A row is named by the effect measure, the contrast, and the subgroup, in
  # that order, which is the order the labelling contract pastes them and so the
  # order the columns sit in.
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

  expected <- joint_expected_estimates(
    joint_cell_means(mods$outcome_mod, dat, mods$lev)
  )
  expect_identical(nrow(est), nrow(expected))
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)

  # Fourteen rows: four cell means, two measures for each of four simple
  # effects, and two for the interaction.
  expect_identical(nrow(est), 14L)

  # The odds ratio is reported nowhere on this surface. It is noncollapsible, so
  # neither a difference of two of them nor a simple effect reported beside one
  # says what it appears to.
  expect_false(any(est$effect == "log(or)"))
})

test_that("a joint exposure reports each component's effect within the other's levels", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  est <- res$estimates

  mu <- joint_cell_means(mods$outcome_mod, dat, mods$lev)
  expected <- joint_expected_estimates(mu)
  simple <- expected[
    expected$group %in% c("e = 0", "e = 1", "a = 0", "a = 1"),
  ]

  # Two of the four are not contrasts against the reference cell: the first
  # treatment's effect when the second is set to 1 contrasts two non-reference
  # cells, and so does the second treatment's effect when the first is set to 1.
  # Those are the rows the vs-reference reporting cannot express at all.
  for (i in seq_len(nrow(simple))) {
    row <- simple[i, ]
    got <- est$estimate[
      est$effect == row$effect &
        est$contrast == row$contrast &
        est$group == row$group
    ]
    expect_equal(
      got,
      row$estimate,
      tolerance = joint_tolerance,
      label = paste(row$effect, row$contrast, row$group)
    )
  }

  # Written out for the one row the contract names explicitly, so the identity
  # the labels claim is pinned against arithmetic rather than against the same
  # helper that built the table above.
  expect_equal(
    est$estimate[
      est$effect == "rd" & est$contrast == "a: 1 vs 0" & est$group == "e = 1"
    ],
    unname(mu[[4]] - mu[[3]]),
    tolerance = joint_tolerance
  )
})

test_that("a joint exposure reports the additive and multiplicative interaction", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  est <- res$estimates

  mu <- joint_cell_means(mods$outcome_mod, dat, mods$lev)
  group <- "e = 1 vs e = 0"

  # The additive interaction is the difference of the two risk differences,
  # which written out over the cells is mu11 - mu10 - mu01 + mu00.
  additive <- est$estimate[
    est$effect == "rd" & est$contrast == "a: 1 vs 0" & est$group == group
  ]
  expect_equal(
    additive,
    unname(mu[[4]] - mu[[2]] - mu[[3]] + mu[[1]]),
    tolerance = joint_tolerance
  )

  # The multiplicative interaction is the same combination of the logs, which is
  # the log of the ratio of the two risk ratios.
  multiplicative <- est$estimate[
    est$effect == "log(rr)" & est$contrast == "a: 1 vs 0" & est$group == group
  ]
  expect_equal(
    multiplicative,
    unname(log(mu[[4]]) - log(mu[[2]]) - log(mu[[3]]) + log(mu[[1]])),
    tolerance = joint_tolerance
  )

  # Interaction is symmetric in the two treatments: the difference between the
  # first treatment's two simple effects is the difference between the second
  # treatment's two. That identity is arithmetic on the same four means, so it
  # is pinned between the two framings of the oracle and to the last bits;
  # asserting it against a reported row instead would pin the solver's residual
  # at a tolerance no solve has to meet.
  a_framing <- (mu[[4]] - mu[[3]]) - (mu[[2]] - mu[[1]])
  e_framing <- (mu[[4]] - mu[[2]]) - (mu[[3]] - mu[[1]])
  expect_equal(unname(a_framing), unname(e_framing), tolerance = 1e-12)
  expect_equal(
    unname(
      (log(mu[[4]]) - log(mu[[3]])) - (log(mu[[2]]) - log(mu[[1]]))
    ),
    unname(
      (log(mu[[4]]) - log(mu[[2]])) - (log(mu[[3]]) - log(mu[[1]]))
    ),
    tolerance = 1e-12
  )

  # Reported once, under the first treatment's framing, rather than twice.
  expect_identical(sum(est$group == group), 2L)

  interaction_rows <- est[est$group == group, , drop = FALSE]
  expect_true(all(is.finite(interaction_rows$std.err)))
  expect_true(all(interaction_rows$std.err > 0))
  expect_true(all(interaction_rows$ci.lower < interaction_rows$ci.upper))
})

test_that("a joint exposure reports a usable standard error for every row", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  est <- res$estimates

  # Read against the row count the joint surface reports, so this is a statement
  # about that surface rather than one any four-level exposure would satisfy.
  expect_identical(nrow(est), 14L)
  expect_true(all(is.finite(est$std.err)))
  expect_true(all(est$std.err > 0))
  expect_true(all(is.finite(est$ci.lower)))
  expect_true(all(is.finite(est$ci.upper)))
  expect_true(all(est$ci.lower < est$ci.upper))
  expect_equal(sqrt(diag(vcov(res))), est$std.err, ignore_attr = TRUE)
})

test_that("a joint exposure's covariance couples the interaction with the cells it is built from", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  covariance <- vcov(res)

  # Every reported row is a function of the same four cell means, solved as one
  # stacked system on one sample. The interaction is a combination of all four,
  # so it covaries with each simple effect it is built from, and two simple
  # effects sharing a cell covary with each other. An assembly that computed the
  # interaction separately and stitched it in would report exact zeros here
  # while matching its own standard errors everywhere.
  #
  # The floor is absolute and far below any covariance shared parameters
  # generate. What it excludes is the structural zero, not a small number.

  couples <- list(
    c("rd a: 1 vs 0 e = 1 vs e = 0", "rd a: 1 vs 0 e = 1"),
    c("rd a: 1 vs 0 e = 1 vs e = 0", "rd a: 1 vs 0 e = 0"),
    c("rd a: 1 vs 0 e = 0", "rd e: 1 vs 0 a = 0"),
    c("mean a = 1, e = 1 overall", "rd a: 1 vs 0 e = 1 vs e = 0")
  )
  for (pair in couples) {
    where <- paste(pair, collapse = " with ")
    # Asked for by name before the entry is read, so a surface that does not
    # report these rows says which label is missing rather than failing on a
    # subscript, and the read is skipped rather than erroring after it.
    present <- all(pair %in% rownames(covariance))
    expect_true(present, label = where)
    if (!present) {
      next
    }
    entry <- covariance[pair[[1]], pair[[2]]]
    expect_true(is.finite(entry), label = where)
    expect_gt(abs(entry), 1e-8, label = where)
  }

  expect_equal(covariance, t(covariance), tolerance = 1e-12)
})

test_that("a joint exposure on a continuous outcome reports means and differences", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  est <- res$estimates

  mu <- joint_cell_means(mods$outcome_mod, dat, mods$lev)
  expected <- joint_expected_estimates(mu, forms = "diff")

  # A continuous outcome reports the counterfactual mean per cell and a
  # difference in means everywhere a binary outcome reports two measures, so the
  # surface is four rows plus one per simple effect plus the interaction.
  expect_identical(nrow(est), 9L)
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)
  expect_equal(est$estimate, expected$estimate, tolerance = joint_tolerance)
})

test_that("the declaration is what turns on the joint reporting", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  declared <- ipw(
    fit_joint_models(dat)$ps_mod,
    fit_joint_models(dat)$outcome_mod
  )
  plain_mods <- fit_joint_models(dat, exposure_name = "plain")
  plain <- ipw(plain_mods$ps_mod, plain_mods$outcome_mod)

  # The two fits read the same cells from the same data under the same labels
  # and differ only in whether the crossing was declared. The plain one reports
  # the vs-reference rows and the declared one reports the joint surface, so the
  # frames differ in shape rather than only in wording.
  expect_identical(nrow(plain$estimates), 9L)
  expect_identical(nrow(declared$estimates), 14L)
  expect_null(plain$estimates[["group"]])
  expect_false(is.null(declared$estimates[["group"]]))
})

# ---- surfaces ----------------------------------------------------------------

test_that("a joint fit labels its rows by effect, contrast, and cell or subgroup", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  labels <- joint_labels(res$estimates)
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

test_that("the tidiers of a joint fit head the group column after the contrast", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  frame <- as.data.frame(res)
  expect_identical(
    names(frame),
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
  expect_identical(frame$term, res$estimates$effect)
  expect_identical(frame$contrast, res$estimates$contrast)
  expect_identical(frame$group, res$estimates$group)

  tidied <- tidy(res)
  expect_s3_class(tidied, "tbl_df")
  expect_identical(names(tidied), names(frame))
  expect_identical(tidied[["group"]], res$estimates$group)

  # There are no ratio rows on the log scale to move, since the surface reports
  # no odds ratio and its log risk ratios are the simple effects and the
  # interaction. Exponentiating relabels those and leaves the cell means alone.
  exponentiated <- as.data.frame(res, exponentiate = TRUE)
  expect_identical(
    exponentiated$term[exponentiated$term %in% c("mean", "rr")],
    c(rep("mean", 4), rep("rr", sum(res$estimates$effect == "log(rr)")))
  )
})

test_that("the conditional reading of a joint fit is the plain coefficient surface", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(dat)
  conditional <- causalgenerics::as_conditional(
    ipw(mods$ps_mod, mods$outcome_mod)
  )

  # The declaration changes what is reported, not what is fitted. The outcome
  # model is a model of the four cells, and the conditional reading is its
  # coefficients with the covariance the joint estimation implies, exactly as it
  # would be for a plain factor.
  expect_identical(coef(conditional), coef(mods$outcome_mod))
  expect_identical(
    rownames(vcov(conditional)),
    names(coef(mods$outcome_mod))
  )
  expect_null(conditional$estimates[["group"]])
})

# ---- refusals ----------------------------------------------------------------

test_that("a joint exposure refuses a focal estimand", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  mods <- fit_joint_models(
    dat,
    estimand = "att",
    focal_level = "a = 1, e = 1"
  )
  expect_identical(estimand(mods$wts), "att")

  # The joint surface standardizes every cell mean to one population, and a
  # focal estimand names one cell as that population. The cell means would then
  # be the effects among the units treated at one corner of the crossing, which
  # is a coherent quantity and not the one the interaction rows are built to
  # combine. Reporting it under these labels would say something else.
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_joint_estimand_error"
  )
  # The remedy is the estimand the surface is defined for, and the message has
  # to name it.
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    regexp = "ate"
  )
  expect_propensity_error(ipw(mods$ps_mod, mods$outcome_mod))

  # The estimand the surface is defined for still works, which is what makes the
  # refusal a restriction rather than a wall.
  ate <- fit_joint_models(dat)
  expect_s3_class(ipw(ate$ps_mod, ate$outcome_mod), "ipw")
})

test_that("a joint exposure refuses .by", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- joint_data()
  # The outcome model carries the exposure-by-modifier term, so the only thing
  # wrong with the request is the combination itself.
  mods <- fit_joint_models(dat, outcome_rhs = "joint * grp + x1")

  # Effect modification of a joint intervention is a three-way question: the
  # interaction between two treatments, within the levels of a third variable.
  # The surface reports neither that nor a sensible projection of it, so the
  # combination is refused rather than answered with one of the two readings.
  expect_no_warning(expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .by = grp),
    class = "propensity_ipw_joint_by_error"
  ))
  expect_propensity_error(ipw(mods$ps_mod, mods$outcome_mod, .by = grp))

  # Each on its own is supported, which is what the refusal is narrowing.
  plain_mods <- fit_joint_models(
    dat,
    exposure_name = "plain",
    outcome_rhs = "plain * grp + x1"
  )
  expect_s3_class(
    ipw(plain_mods$ps_mod, plain_mods$outcome_mod, .by = grp),
    "ipw"
  )
})
