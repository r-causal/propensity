# The counterfactual mean at each exposure level, reported as a row of its own.
#
# Every binary and categorical `ipw()` result already estimates those means: the
# contrasts it reports are transformations of them, and on the M-estimation path
# they are parameters of the stacked system, so the sandwich already covers them.
# What was missing was reporting them, which left a reader with the difference
# between two risks and no way to say what either risk was.
#
# The rows lead the table, one per level in level order with the reference level
# first, under the effect label "mean" and keyed by the level they belong to. The
# contrast rows that follow name the pair they compare in the same column, so a
# binary result reads the way a categorical one does rather than leaving the one
# comparison it makes unnamed.
#
# What is pinned here is that contract: the row set and its order, the labels,
# the agreement between the two standard error routes, and that every surface
# built from the result carries the new rows. The numbers themselves are checked
# against a g-computation plug-in and against the contrasts the same table
# reports, rather than against recorded values, so nothing here has to be
# re-recorded when a fixture changes.

# ---- fixtures ----------------------------------------------------------------

sim_level_means_binary <- function(seed = 4321, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.4 * x1 - 0.5 * x2))
  y <- rbinom(n, 1, plogis(-0.3 + 0.9 * z + 0.6 * x1 - 0.4 * x2))
  yc <- 1.2 + 0.7 * z + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  v <- factor(ifelse(x2 == 1, "hi", "lo"), levels = c("lo", "hi"))
  data.frame(x1, x2, z, y, yc, v)
}

sim_level_means_categorical <- function(seed = 4322, n = 900) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  eta_b <- -0.2 + 0.5 * x1
  eta_c <- 0.1 - 0.4 * x1
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
    plogis(-0.3 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.5 * x1)
  )
  data.frame(x1, a, y)
}

# A binary exposure with an exposure-only outcome model, which is the one shape
# both standard error routes accept, so the two can be compared on one fixture.
fit_level_means_binary <- function(dat, outcome_family = "binomial") {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_mod)
  )
  outcome_mod <- if (identical(outcome_family, "binomial")) {
    glm(
      y ~ z,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  } else {
    lm(yc ~ z, data = dat, weights = wts)
  }
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# The modifier enters the outcome model, so a `.by` request has the
# exposure-by-modifier term it warns about the absence of.
fit_level_means_by <- function(dat) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_mod)
  )
  outcome_mod <- glm(
    y ~ z * v + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

fit_level_means_categorical <- function(dat) {
  ps_mod <- nnet::multinom(
    a ~ x1,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      predict(ps_mod, type = "probs"),
      dat$a,
      exposure_type = "categorical"
    )
  )
  outcome_mod <- glm(
    y ~ a,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# ---- oracles -----------------------------------------------------------------

# The counterfactual mean at one exposure level: the outcome model's predictions
# with the exposure column set to that level for every unit, averaged over the
# whole sample. This is the ate standardization, which the fixtures here all use.
# A tilted estimand standardizes over h(e) instead, and the mean rows it reports
# are pinned against a tilt-weighted oracle in test-ipw-by.R and
# test-ipw-by-categorical.R rather than a second time here.
level_mean_plugin <- function(outcome_mod, data, exposure_name, value) {
  d <- data
  d[[exposure_name]] <- value
  mean(predict(outcome_mod, newdata = d, type = "response"))
}

# The rows of an estimates table under the effect label "mean".
mean_rows <- function(estimates) {
  estimates[estimates$effect == "mean", , drop = FALSE]
}

contrast_rows <- function(estimates) {
  estimates[estimates$effect != "mean", , drop = FALSE]
}

# ---- the row set a binary result reports -------------------------------------

test_that("a binary mestimation fit leads with one mean row per exposure level", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  est <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")$estimates

  # The means come first, in level order with the reference level leading, and
  # each is keyed by the level it belongs to rather than by a pair of levels.
  expect_identical(
    est$effect,
    c("mean", "mean", "rd", "log(rr)", "log(or)")
  )
  expect_identical(est$contrast, c("0", "1", rep("1 vs 0", 3)))

  # The contrast column leads the frame, ahead of the numeric columns, the way a
  # categorical result already writes it.
  expect_identical(names(est)[seq(1, 2)], c("effect", "contrast"))
})

test_that("a binary gaussian outcome reports its means beside the difference", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat, outcome_family = "gaussian")
  est <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")$estimates

  # A continuous outcome reports one contrast rather than three, and the mean
  # rows are the counterfactual means on the outcome's own scale, so nothing
  # confines them to (0, 1).
  expect_identical(est$effect, c("mean", "mean", "diff"))
  expect_identical(est$contrast, c("0", "1", "1 vs 0"))
  expect_true(all(est$estimate[seq(1, 2)] > 1))
})

test_that("a binary contrast row is exactly the transform of the means beside it", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  est <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")$estimates

  mu <- mean_rows(est)$estimate
  names(mu) <- mean_rows(est)$contrast
  effects <- contrast_rows(est)$estimate
  names(effects) <- contrast_rows(est)$effect

  # A counterfactual risk is a probability, so both mean rows lie strictly
  # inside the unit interval and the transforms below are all defined.
  expect_true(all(mu > 0 & mu < 1))

  expect_equal(effects[["rd"]], mu[["1"]] - mu[["0"]], tolerance = 1e-8)
  expect_equal(
    effects[["log(rr)"]],
    log(mu[["1"]]) - log(mu[["0"]]),
    tolerance = 1e-8
  )
  expect_equal(
    effects[["log(or)"]],
    qlogis(mu[["1"]]) - qlogis(mu[["0"]]),
    tolerance = 1e-8
  )
})

test_that("the binary mean rows are the mu block of the stacked system", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # The stacked system carries the marginal means as parameters named `mu1` and
  # `mu0`. The reported rows are those parameters, reordered so the reference
  # level leads, rather than a second computation of the same quantity.
  theta <- stats::coef(res$fit)
  reported <- mean_rows(res$estimates)
  expect_identical(reported$contrast, c("0", "1"))
  expect_equal(
    reported$estimate,
    unname(c(theta[["mu0"]], theta[["mu1"]])),
    tolerance = 1e-12
  )
})

test_that("the binary mean rows match a g-computation plug-in", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  reported <- mean_rows(res$estimates)$estimate
  plugin <- vapply(
    c(0, 1),
    function(value) {
      level_mean_plugin(mods$outcome_mod, dat, "z", value)
    },
    numeric(1)
  )
  expect_equal(reported, plugin, tolerance = 1e-8)
})

# ---- the two standard error routes agree -------------------------------------

test_that("the two SE routes report the same mean rows", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  mest <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
  lin <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  expect_identical(lin$estimates$effect, mest$estimates$effect)
  expect_identical(lin$estimates$contrast, mest$estimates$contrast)
  expect_identical(nrow(mean_rows(mest$estimates)), 2L)
  expect_identical(nrow(mean_rows(lin$estimates)), 2L)

  # The estimates are the same g-computation on both routes, so they agree to
  # the numerical accuracy of the solve.
  expect_equal(
    mean_rows(lin$estimates)$estimate,
    mean_rows(mest$estimates)$estimate,
    tolerance = 1e-8
  )

  # The standard errors are not the same computation: linearization is a
  # first-order approximation of the variance the sandwich returns exactly, so
  # they agree to a few digits rather than to the last one.
  mest_se <- mean_rows(mest$estimates)$std.err
  lin_se <- mean_rows(lin$estimates)$std.err
  expect_true(all(is.finite(mest_se)))
  expect_true(all(is.finite(lin_se)))
  expect_true(all(mest_se > 0))
  expect_true(all(lin_se > 0))
  expect_true(all(abs(lin_se - mest_se) / mest_se < 0.03))
})

test_that("a linearization fit reports its mean rows with usable inference", {
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  est <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    se_method = "linearization"
  )$estimates
  means <- mean_rows(est)

  # This route stores no stacked fit, so the mean rows and their inference come
  # from the influence functions the contrast rows already came from.
  expect_identical(nrow(means), 2L)
  expect_true(all(is.finite(means$z)))
  expect_true(all(means$ci.lower < means$estimate))
  expect_true(all(means$estimate < means$ci.upper))

  # The interval is the normal approximation the rest of the table is built on.
  half <- (means$ci.upper - means$estimate) / means$std.err
  expect_equal(half, rep(qnorm(0.975), 2), tolerance = 1e-8)
})

# ---- the covariance covers the new rows --------------------------------------

test_that("the stored covariance covers the binary mean rows on both routes", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  labels <- c(
    "mean 0",
    "mean 1",
    paste(c("rd", "log(rr)", "log(or)"), "1 vs 0")
  )

  for (se in c("mestimation", "linearization")) {
    res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = se)
    covariance <- attr(res$estimates, "ipw_vcov", exact = TRUE)

    expect_identical(dim(covariance), c(5L, 5L), label = se)
    expect_identical(dimnames(covariance), list(labels, labels), label = se)

    # The diagonal is the standard errors the table reports, which is what keeps
    # the block and the table describing one fit.
    expect_equal(
      sqrt(diag(covariance)),
      res$estimates$std.err,
      tolerance = 1e-8,
      ignore_attr = TRUE
    )

    # The two arms are estimated from one sample through one pair of models, so
    # their means covary rather than being independent estimates.
    expect_true(is.finite(covariance[["mean 0", "mean 1"]]))
    expect_true(covariance[["mean 0", "mean 1"]] != 0)

    # A risk difference is the difference of the two means, so it moves against
    # the reference arm's mean and with the exposed arm's.
    expect_lt(covariance[["mean 0", "rd 1 vs 0"]], 0)
    expect_gt(covariance[["mean 1", "rd 1 vs 0"]], 0)
  }
})

# ---- the accessors and tidiers carry the new rows ----------------------------

test_that("the accessors report the binary mean rows", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
  labels <- c(
    "mean 0",
    "mean 1",
    paste(c("rd", "log(rr)", "log(or)"), "1 vs 0")
  )

  expect_identical(names(coef(res)), labels)
  expect_equal(unname(coef(res)), res$estimates$estimate, tolerance = 1e-12)
  expect_identical(dimnames(vcov(res)), list(labels, labels))
  expect_identical(rownames(confint(res)), labels)

  bounds <- confint(res)
  expect_equal(bounds[, 1], res$estimates$ci.lower, ignore_attr = TRUE)
  expect_equal(bounds[, 2], res$estimates$ci.upper, ignore_attr = TRUE)
})

test_that("tidy() reports the binary mean rows under term and contrast", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  tidied <- tidy(res, conf.int = TRUE)
  expect_identical(sum(tidied$term == "mean"), 2L)
  expect_identical(tidied$term, res$estimates$effect)
  expect_identical(tidied$contrast, res$estimates$contrast)
  expect_identical(tidied$estimate, res$estimates$estimate)
  expect_identical(tidied$conf.low, res$estimates$ci.lower)
})

test_that("tidy(exponentiate = TRUE) leaves the mean rows alone", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_binary(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  plain <- tidy(res, conf.int = TRUE)
  expo <- tidy(res, conf.int = TRUE, exponentiate = TRUE)

  # A counterfactual risk is not a ratio, so exponentiating the ratio rows must
  # not touch it and must not rename it to something it is not.
  is_mean <- plain$term == "mean"
  expect_identical(expo$term[is_mean], rep("mean", 2))
  expect_identical(expo$contrast[is_mean], plain$contrast[is_mean])
  expect_identical(expo$estimate[is_mean], plain$estimate[is_mean])
  expect_identical(expo$conf.low[is_mean], plain$conf.low[is_mean])
  expect_identical(expo$conf.high[is_mean], plain$conf.high[is_mean])

  # The ratio rows still move, so the agreement above is with rows that were
  # left alone rather than with a table nothing happened to.
  ratio <- plain$term %in% c("log(rr)", "log(or)")
  expect_identical(expo$term[ratio], c("rr", "or"))
  expect_equal(
    expo$estimate[ratio],
    exp(plain$estimate[ratio]),
    tolerance = 1e-12
  )
})

# ---- a categorical exposure --------------------------------------------------

test_that("a categorical fit reports one mean row per level in level order", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_level_means_categorical()
  mods <- fit_level_means_categorical(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  est <- res$estimates

  expect_identical(
    est$effect,
    c(rep("mean", 3), rep(c("rd", "log(rr)", "log(or)"), times = 2))
  )
  expect_identical(
    est$contrast,
    c(c("a", "b", "c"), rep(c("b vs a", "c vs a"), each = 3))
  )

  # The mean rows are the categorical mu block, one parameter per level, in the
  # order the fit codes the levels.
  theta <- stats::coef(res$fit)
  expect_equal(
    mean_rows(est)$estimate,
    unname(theta[c("mu_a", "mu_b", "mu_c")]),
    tolerance = 1e-12
  )
})

test_that("the categorical mean rows match a g-computation plug-in", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_level_means_categorical()
  mods <- fit_level_means_categorical(dat)
  est <- ipw(mods$ps_mod, mods$outcome_mod)$estimates

  plugin <- vapply(
    levels(dat$a),
    function(level) {
      level_mean_plugin(
        mods$outcome_mod,
        dat,
        "a",
        factor(level, levels = levels(dat$a))
      )
    },
    numeric(1)
  )
  expect_equal(mean_rows(est)$estimate, unname(plugin), tolerance = 1e-8)
})

test_that("each categorical contrast is the transform of the two means it names", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_level_means_categorical()
  mods <- fit_level_means_categorical(dat)
  est <- ipw(mods$ps_mod, mods$outcome_mod)$estimates

  mu <- mean_rows(est)$estimate
  names(mu) <- mean_rows(est)$contrast
  contrasts <- contrast_rows(est)

  for (level in c("b", "c")) {
    label <- paste(level, "vs a")
    rows <- contrasts[contrasts$contrast == label, , drop = FALSE]
    values <- rows$estimate
    names(values) <- rows$effect
    expect_equal(
      values[["rd"]],
      mu[[level]] - mu[["a"]],
      tolerance = 1e-8,
      label = label
    )
    expect_equal(
      values[["log(rr)"]],
      log(mu[[level]]) - log(mu[["a"]]),
      tolerance = 1e-8,
      label = label
    )
  }
})

test_that("a categorical fit keys its covariance by the level each mean belongs to", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_level_means_categorical()
  mods <- fit_level_means_categorical(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  covariance <- attr(res$estimates, "ipw_vcov", exact = TRUE)

  expect_identical(dim(covariance), c(9L, 9L))
  expect_identical(
    rownames(covariance)[seq(1, 3)],
    c("mean a", "mean b", "mean c")
  )
  expect_identical(anyDuplicated(rownames(covariance)), 0L)
  expect_equal(
    sqrt(diag(covariance)),
    res$estimates$std.err,
    tolerance = 1e-8,
    ignore_attr = TRUE
  )

  # Both contrasts are taken against the reference level, so both covary with
  # the reference level's mean.
  expect_true(covariance[["mean a", "rd b vs a"]] != 0)
  expect_true(covariance[["mean a", "rd c vs a"]] != 0)
})

# ---- a `.by` request ---------------------------------------------------------

test_that("a .by fit reports the means overall and within every stratum", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_by(dat)
  est <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)$estimates
  means <- mean_rows(est)

  # Two overall means and two within each of the two strata. The stratum
  # contrasts compare two strata's effects and have no mean of their own.
  expect_identical(nrow(means), 6L)
  expect_identical(
    means$group,
    c(rep("overall", 2), rep("v = lo", 2), rep("v = hi", 2))
  )
  expect_identical(means$contrast, rep(c("0", "1"), times = 3))
  expect_false(any(means$group == "v = hi vs v = lo"))

  # The overall means lead the table and the stratum means follow the overall
  # contrasts, so a block of means is never split by the contrasts built from it.
  expect_identical(which(est$effect == "mean"), c(1L, 2L, 6L, 7L, 8L, 9L))
})

test_that("a .by stratum contrast is the difference of that stratum's means", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_by(dat)
  est <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)$estimates

  for (group in c("overall", "v = lo", "v = hi")) {
    rows <- est[est$group == group, , drop = FALSE]
    mu <- rows$estimate[rows$effect == "mean"]
    rd <- rows$estimate[rows$effect == "rd"]
    expect_equal(rd, mu[[2]] - mu[[1]], tolerance = 1e-8, label = group)
  }
})

test_that("a .by fit gives every mean row a usable standard error", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_by(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)
  means <- mean_rows(res$estimates)

  expect_identical(nrow(means), 6L)
  expect_true(all(is.finite(means$std.err)))
  expect_true(all(means$std.err > 0))

  # The stratum means are stacked parameters, so their standard errors are the
  # diagonal of the same sandwich the rest of the table is read off.
  expect_equal(
    sqrt(diag(vcov(res))),
    res$estimates$std.err,
    tolerance = 1e-8,
    ignore_attr = TRUE
  )

  # A stratum mean is estimated from that stratum's units alone, so it carries
  # more uncertainty than the whole-sample mean of the same arm.
  overall <- means$std.err[means$group == "overall"]
  strata <- means$std.err[means$group != "overall"]
  expect_true(all(strata > min(overall)))
})

test_that("a .by fit labels every mean row by level and stratum together", {
  skip_if_not_installed("deli")
  dat <- sim_level_means_binary()
  mods <- fit_level_means_by(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, .by = v)

  labels <- paste(
    res$estimates$effect,
    res$estimates$contrast,
    res$estimates$group
  )
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(names(coef(res)), labels)
  expect_true(all(c("mean 0 v = lo", "mean 1 v = hi") %in% labels))
})
