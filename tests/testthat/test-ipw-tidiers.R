# Tests for the broom tidiers on `ipw` results. `tidy()` reshapes the object's
# `estimates` table into the column names the tidymodels broom developer guide
# specifies, one row per estimate, adding interval bounds only when asked for
# them. The values themselves are the object's own: the tidier renames and
# selects, and recomputes only when a requested confidence level differs from
# the one the fit stored.
#
# `glance()` describes the fit rather than its estimates: a single row naming the
# estimand, the exposure type, the standard error method, the two models, and the
# number of observations the outcome model saw. That row carries the same columns
# with the same types on every route.
#
# Coverage spans every route that builds an `ipw` object: a binary glm exposure
# under both standard error methods, a categorical exposure through
# nnet::multinom, and a continuous exposure through lm.

# ---- data simulators ---------------------------------------------------------

sim_tidy_binary <- function(seed = 2024, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.6 * x2))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * x1 - 0.3 * x2))
  yc <- 1.5 + 0.8 * z + 0.6 * x1 - 0.4 * x2 + rnorm(n)
  data.frame(x1, x2, z, y, yc)
}

sim_tidy_categorical <- function(seed = 2024, n = 700) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  eta_b <- -0.2 + 0.5 * x1 + 0.3 * x2
  eta_c <- 0.1 - 0.4 * x1 + 0.6 * x2
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
  data.frame(x1, x2, a, y)
}

sim_tidy_continuous <- function(seed = 2024, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  yc <- 1 + 0.6 * A + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  yb <- rbinom(n, 1, plogis(-0.5 + 0.4 * A + 0.3 * x1 - 0.2 * x2))
  data.frame(x1, x2, A, yc, yb)
}

# ---- model fitting -----------------------------------------------------------

# Binary exposure: a logit propensity score model and an exposure-only weighted
# outcome model. The outcome model carries no covariate so the same pair of
# models serves the M-estimation and the linearization paths, the latter being
# restricted to marginal outcome models. `estimand` picks the weighting scheme,
# so a fit reporting something other than the ATE is available here rather than
# through a second fixture.
fit_tidy_binary_models <- function(
  dat,
  outcome_family = "binomial",
  estimand = "ate"
) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wt_fun <- if (estimand == "att") wt_att else wt_ate
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_fun(ps_mod)
  )
  outcome_var <- if (outcome_family == "binomial") "y" else "yc"
  fmla <- stats::reformulate("z", response = outcome_var)
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
  list(ps_mod = ps_mod, outcome_mod = outcome_mod)
}

# Categorical exposure: a multinomial propensity score model and an
# exposure-only weighted outcome model. The tightened `reltol` matches the
# convention in the other categorical test files, keeping the fitted
# coefficients at the score root the stacked system re-solves.
fit_tidy_categorical_models <- function(dat) {
  ps_mod <- nnet::multinom(
    a ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_named <- unname(predict(ps_mod, type = "probs"))
  colnames(ps_named) <- ps_mod$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_named, dat$a, exposure_type = "categorical")
  )
  outcome_mod <- glm(
    y ~ a,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(ps_mod = ps_mod, outcome_mod = outcome_mod)
}

# Continuous exposure: an lm propensity score model of the exposure, stabilized
# continuous ATE weights, and a weighted marginal structural model with the one
# exposure term the continuous path requires. `ps_model` chooses between the two
# propensity model classes that path accepts, an `lm` and a gaussian `glm`, which
# agree on fitted values and so produce the same weights.
fit_tidy_continuous_models <- function(
  dat,
  outcome_family = "gaussian",
  ps_model = "lm"
) {
  ps_mod <- if (ps_model == "glm") {
    glm(A ~ x1 + x2, data = dat, family = gaussian())
  } else {
    lm(A ~ x1 + x2, data = dat)
  }
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(fitted(ps_mod)),
      dat$A,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )
  outcome_var <- if (outcome_family == "binomial") "yb" else "yc"
  fmla <- stats::reformulate("A", response = outcome_var)
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
  list(ps_mod = ps_mod, outcome_mod = outcome_mod)
}

# ---- expected shape ----------------------------------------------------------

# The broom column contract. `comparison` is the categorical extra column and
# sits immediately after `term`, following the precedent broom sets with the
# `y.level` column of `tidy.multinom`. The interval bounds close the frame.
tidy_names <- function(conf_int = FALSE, comparison = FALSE) {
  nms <- c("term", "estimate", "std.error", "statistic", "p.value")
  if (comparison) {
    nms <- append(nms, "comparison", after = 1L)
  }
  if (conf_int) {
    nms <- c(nms, "conf.low", "conf.high")
  }
  nms
}

# Assert the whole contract against the object the tidier was called on: a
# tibble of the documented columns whose values are the `estimates` columns
# themselves, in the order the object stores them and with nothing dropped.
expect_tidy_contract <- function(
  tidied,
  result,
  conf_int = FALSE,
  comparison = FALSE
) {
  estimates <- result$estimates
  expect_s3_class(tidied, c("tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_named(tidied, tidy_names(conf_int, comparison))
  expect_identical(nrow(tidied), nrow(estimates))
  expect_identical(tidied$term, estimates$effect)
  expect_identical(tidied$estimate, estimates$estimate)
  expect_identical(tidied$std.error, estimates$std.err)
  expect_identical(tidied$statistic, estimates$z)
  expect_identical(tidied$p.value, estimates$p.value)
  if (comparison) {
    expect_identical(tidied$comparison, estimates$comparison)
  }
  if (conf_int) {
    expect_identical(tidied$conf.low, estimates$ci.lower)
    expect_identical(tidied$conf.high, estimates$ci.upper)
  }
  invisible(tidied)
}

# The glance column contract, as a named vector of the storage types the columns
# hold. Every route reports these columns, in this order, at these types.
glance_types <- function() {
  c(
    estimand = "character",
    exposure_type = "character",
    se_method = "character",
    wt_model = "character",
    outcome_model = "character",
    nobs = "integer"
  )
}

# Assert the whole contract against the result glanced at: one row of the
# documented columns, each describing the fit rather than restating a constant.
# `exposure_type` and `nobs` are the caller's expectations because the result
# records neither, the first being read off the propensity model and the second
# off the outcome model.
expect_glance_contract <- function(glanced, result, exposure_type, nobs) {
  expect_s3_class(glanced, c("tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_named(glanced, names(glance_types()))
  expect_identical(
    unname(vapply(glanced, typeof, character(1))),
    unname(glance_types())
  )
  expect_identical(nrow(glanced), 1L)

  expect_identical(glanced$estimand, result$estimand)
  expect_identical(glanced$exposure_type, exposure_type)
  expect_identical(glanced$se_method, result$se_method)
  expect_identical(glanced$wt_model, class(result$wt_mod)[[1]])
  expect_identical(glanced$outcome_model, class(result$outcome_mod)[[1]])
  expect_identical(glanced$nobs, nobs)
  expect_identical(glanced$nobs, as.integer(stats::nobs(result$outcome_mod)))

  invisible(glanced)
}

# ---- binary exposure, M-estimation ------------------------------------------

test_that("tidy() returns the broom columns for a binary mestimation fit", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  tidied <- tidy(res)

  # conf.int defaults to FALSE, so no bounds appear, and the effect labels
  # arrive on the log scale because exponentiate defaults to FALSE
  expect_tidy_contract(tidied, res)
  expect_identical(tidied$term, c("rd", "log(rr)", "log(or)"))
})

test_that("tidy(conf.int = TRUE) returns the stored bounds for a binary fit", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # conf.level defaults to 0.95, the level the fit was built at, so the bounds
  # are the stored ones rather than a recomputation that merely agrees with them
  tidied <- tidy(res, conf.int = TRUE)
  expect_tidy_contract(tidied, res, conf_int = TRUE)
})

test_that("tidy() reports the diff row for a gaussian binary outcome", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  tidied <- tidy(res, conf.int = TRUE)
  expect_tidy_contract(tidied, res, conf_int = TRUE)
  expect_identical(tidied$term, "diff")

  # a table with no ratio rows is untouched by exponentiate
  expect_identical(tidy(res, conf.int = TRUE, exponentiate = TRUE), tidied)
})

# ---- binary exposure, linearization ------------------------------------------

test_that("tidy() returns the broom columns for a binary linearization fit", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # The linearization path stores no deli fit, so the tidier has only the
  # estimates table to read. Its contract is the same as the M-estimation one.
  expect_null(res$fit)
  expect_tidy_contract(tidy(res), res)
  expect_tidy_contract(tidy(res, conf.int = TRUE), res, conf_int = TRUE)
})

# ---- categorical exposure ----------------------------------------------------

test_that("tidy() keeps the comparison column after term for a categorical fit", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_tidy_categorical()
  mods <- fit_tidy_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_tidy_contract(tidy(res), res, comparison = TRUE)

  tidied <- tidy(res, conf.int = TRUE)
  expect_tidy_contract(tidied, res, conf_int = TRUE, comparison = TRUE)

  # every effect of every comparison keeps its own row, in the stored order
  expect_identical(tidied$term, rep(c("rd", "log(rr)", "log(or)"), times = 2))
  expect_identical(tidied$comparison, rep(c("b vs a", "c vs a"), each = 3))
})

# ---- continuous exposure -----------------------------------------------------

test_that("tidy() returns the single slope row for a continuous exposure fit", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_continuous()
  mods <- fit_tidy_continuous_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  tidied <- tidy(res, conf.int = TRUE)

  # a continuous exposure carries no comparison column
  expect_tidy_contract(tidied, res, conf_int = TRUE)
  expect_identical(tidied$term, "slope")
})

test_that("tidy() returns the log odds ratio row for a continuous logit fit", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_continuous()
  mods <- fit_tidy_continuous_models(dat, outcome_family = "binomial")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  tidied <- tidy(res, conf.int = TRUE)
  expect_tidy_contract(tidied, res, conf_int = TRUE)
  expect_identical(tidied$term, "log(or)")
})

# ---- confidence level --------------------------------------------------------

test_that("tidy() recomputes the interval at a non-default conf.level", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
  estimates <- res$estimates

  tidied <- tidy(res, conf.int = TRUE, conf.level = 0.9)
  half <- qnorm(0.95) * estimates$std.err

  expect_equal(tidied$conf.low, estimates$estimate - half, tolerance = 1e-10)
  expect_equal(tidied$conf.high, estimates$estimate + half, tolerance = 1e-10)

  # a 0.9 request is not the stored 0.95 interval read back
  expect_true(all(tidied$conf.low > estimates$ci.lower))
  expect_true(all(tidied$conf.high < estimates$ci.upper))
})

test_that("tidy() recomputes when the fit stores a non-default conf.level", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    se_method = "mestimation",
    conf_level = 0.9
  )
  estimates <- res$estimates

  # asking for the level the fit stores returns those bounds exactly
  at_stored <- tidy(res, conf.int = TRUE, conf.level = 0.9)
  expect_identical(at_stored$conf.low, estimates$ci.lower)
  expect_identical(at_stored$conf.high, estimates$ci.upper)

  # the 0.95 default is the tidier's, not the fit's, so it is recomputed rather
  # than served from the stored 0.9 columns
  at_default <- tidy(res, conf.int = TRUE)
  half <- qnorm(0.975) * estimates$std.err
  expect_equal(
    at_default$conf.low,
    estimates$estimate - half,
    tolerance = 1e-10
  )
  expect_equal(
    at_default$conf.high,
    estimates$estimate + half,
    tolerance = 1e-10
  )
  expect_true(all(at_default$conf.low < estimates$ci.lower))
})

test_that("tidy() ignores conf.level when conf.int is FALSE", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # broom's convention: conf.level only governs bounds that were asked for, so
  # supplying it without conf.int must not add columns or change any value
  expect_identical(tidy(res, conf.level = 0.9), tidy(res))
})

test_that("tidy() rejects a conf.level that is not one number inside (0, 1)", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # The endpoints are excluded: a level of 0 or 1 has no normal-approximation
  # interval to report.
  expect_error(
    tidy(res, conf.int = TRUE, conf.level = 1),
    class = "propensity_conf_level_error"
  )
  expect_error(
    tidy(res, conf.int = TRUE, conf.level = 0),
    class = "propensity_conf_level_error"
  )
  expect_error(
    tidy(res, conf.int = TRUE, conf.level = c(0.9, 0.95)),
    class = "propensity_conf_level_error"
  )
  expect_error(
    tidy(res, conf.int = TRUE, conf.level = NA_real_),
    class = "propensity_conf_level_error"
  )
  expect_error(
    tidy(res, conf.int = TRUE, conf.level = NULL),
    class = "propensity_conf_level_error"
  )

  # The level is validated whether or not bounds were asked for, so a bad value
  # is reported rather than quietly going unused.
  expect_error(
    tidy(res, conf.level = "0.95"),
    class = "propensity_conf_level_error"
  )
})

test_that("tidy() rejects an argument that lands in the dots", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # The dots exist to match the generic. A misspelled argument absorbed there
  # would otherwise return the default interval as though it had been asked for.
  expect_error(
    tidy(res, conf.int = TRUE, conf_level = 0.9),
    class = "rlib_error_dots_nonempty"
  )
})

# ---- exponentiate ------------------------------------------------------------

test_that("tidy(exponentiate = TRUE) relabels and exponentiates ratio rows", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  plain <- tidy(res, conf.int = TRUE)
  expo <- tidy(res, conf.int = TRUE, exponentiate = TRUE)

  expect_named(expo, tidy_names(conf_int = TRUE))
  expect_identical(expo$term, c("rd", "rr", "or"))

  # the risk difference row is on its own scale already and is left alone
  ratio <- c(FALSE, TRUE, TRUE)
  expect_identical(expo$estimate[!ratio], plain$estimate[!ratio])
  expect_identical(expo$conf.low[!ratio], plain$conf.low[!ratio])
  expect_identical(expo$conf.high[!ratio], plain$conf.high[!ratio])

  expect_equal(
    expo$estimate[ratio],
    exp(plain$estimate[ratio]),
    tolerance = 1e-12
  )
  expect_equal(
    expo$conf.low[ratio],
    exp(plain$conf.low[ratio]),
    tolerance = 1e-12
  )
  expect_equal(
    expo$conf.high[ratio],
    exp(plain$conf.high[ratio]),
    tolerance = 1e-12
  )

  # the inference columns describe the log-scale estimate and stay there
  expect_identical(expo$std.error, plain$std.error)
  expect_identical(expo$statistic, plain$statistic)
  expect_identical(expo$p.value, plain$p.value)

  # the relabeling does not depend on bounds being requested
  bare <- tidy(res, exponentiate = TRUE)
  expect_named(bare, tidy_names())
  expect_identical(bare$term, c("rd", "rr", "or"))
  expect_identical(bare$estimate, expo$estimate)
})

test_that("tidy(exponentiate = TRUE) matches as.data.frame on a binary fit", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # `exponentiate` is documented as behaving exactly as it does in
  # `as.data.frame.ipw`, so the two surfaces must agree column for column.
  df <- as.data.frame(res, exponentiate = TRUE)
  tidied <- tidy(res, conf.int = TRUE, exponentiate = TRUE)

  expect_identical(tidied$term, df$effect)
  expect_identical(tidied$estimate, df$estimate)
  expect_identical(tidied$std.error, df$std.err)
  expect_identical(tidied$statistic, df$z)
  expect_identical(tidied$p.value, df$p.value)
  expect_identical(tidied$conf.low, df$ci.lower)
  expect_identical(tidied$conf.high, df$ci.upper)
})

test_that("tidy(exponentiate = TRUE) relabels ratio rows per comparison", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_tidy_categorical()
  mods <- fit_tidy_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  df <- as.data.frame(res, exponentiate = TRUE)
  tidied <- tidy(res, conf.int = TRUE, exponentiate = TRUE)

  expect_identical(tidied$term, rep(c("rd", "rr", "or"), times = 2))
  expect_identical(tidied$comparison, rep(c("b vs a", "c vs a"), each = 3))
  expect_identical(tidied$term, df$effect)
  expect_identical(tidied$comparison, df$comparison)
  expect_identical(tidied$estimate, df$estimate)
  expect_identical(tidied$conf.low, df$ci.lower)
  expect_identical(tidied$conf.high, df$ci.upper)
})

test_that("tidy(exponentiate = TRUE) relabels the continuous log odds ratio", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_continuous()
  mods <- fit_tidy_continuous_models(dat, outcome_family = "binomial")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  plain <- tidy(res, conf.int = TRUE)
  expo <- tidy(res, conf.int = TRUE, exponentiate = TRUE)

  expect_identical(expo$term, "or")
  expect_equal(expo$estimate, exp(plain$estimate), tolerance = 1e-12)
  expect_equal(expo$conf.low, exp(plain$conf.low), tolerance = 1e-12)
  expect_equal(expo$conf.high, exp(plain$conf.high), tolerance = 1e-12)
  expect_identical(expo$std.error, plain$std.error)
})

test_that("tidy() exponentiates the bounds it recomputed, not the stored ones", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
  estimates <- res$estimates

  # The two arguments compose: the interval is rebuilt at the requested level on
  # the log scale, and only then exponentiated on the ratio rows.
  tidied <- tidy(res, conf.int = TRUE, conf.level = 0.9, exponentiate = TRUE)
  half <- qnorm(0.95) * estimates$std.err
  ratio <- estimates$effect %in% c("log(rr)", "log(or)")

  expect_equal(
    tidied$conf.low,
    ifelse(
      ratio,
      exp(estimates$estimate - half),
      estimates$estimate - half
    ),
    tolerance = 1e-10
  )
  expect_equal(
    tidied$conf.high,
    ifelse(
      ratio,
      exp(estimates$estimate + half),
      estimates$estimate + half
    ),
    tolerance = 1e-10
  )
})

# ---- glance ------------------------------------------------------------------

test_that("glance() summarizes a binary mestimation fit in one row", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  glanced <- glance(res)

  # A logit propensity model describes a binary exposure, and the three estimate
  # rows the fit holds summarize to this one.
  expect_glance_contract(glanced, res, exposure_type = "binary", nobs = 600L)
  expect_identical(glanced$estimand, "ate")
  expect_identical(glanced$se_method, "mestimation")
  expect_identical(glanced$wt_model, "glm")
  expect_identical(glanced$outcome_model, "glm")
})

test_that("glance() reports the linearization standard error method", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # The linearization path stores no deli fit, so no column may be read from
  # one. The row is otherwise the M-estimation row.
  expect_null(res$fit)
  glanced <- glance(res)
  expect_glance_contract(glanced, res, exposure_type = "binary", nobs = 600L)
  expect_identical(glanced$se_method, "linearization")
})

test_that("glance() reports a categorical exposure for a multinom fit", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_tidy_categorical()
  mods <- fit_tidy_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  glanced <- glance(res)

  # Three exposure levels and six estimate rows still summarize to one row.
  expect_glance_contract(
    glanced,
    res,
    exposure_type = "categorical",
    nobs = 700L
  )
  expect_identical(glanced$wt_model, "multinom")
  expect_identical(glanced$outcome_model, "glm")
})

test_that("glance() reports a continuous exposure for an lm propensity model", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_continuous()
  mods <- fit_tidy_continuous_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  glanced <- glance(res)
  expect_glance_contract(
    glanced,
    res,
    exposure_type = "continuous",
    nobs = 600L
  )
  expect_identical(glanced$wt_model, "lm")
  expect_identical(glanced$outcome_model, "lm")
})

test_that("glance() reads the exposure type from the propensity model family", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_continuous()
  mods <- fit_tidy_continuous_models(dat, ps_model = "glm")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  glanced <- glance(res)

  # A glm propensity model is binary or continuous depending on its family, and
  # `ipw()` routes this gaussian one down the continuous path, so the class of
  # the model alone cannot decide the column.
  expect_glance_contract(
    glanced,
    res,
    exposure_type = "continuous",
    nobs = 600L
  )
  expect_identical(glanced$wt_model, "glm")
})

test_that("glance() reports the estimand the fit was built for", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat, estimand = "att")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # The estimand comes from the weights the outcome model was fit with, so the
  # column follows it rather than naming the ATE the other routes happen to use.
  expect_identical(res$estimand, "att")
  glanced <- glance(res)
  expect_glance_contract(glanced, res, exposure_type = "binary", nobs = 600L)
  expect_identical(glanced$estimand, "att")
})

test_that("glance() returns the same columns and types on every build route", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  binary <- fit_tidy_binary_models(sim_tidy_binary())
  categorical <- fit_tidy_categorical_models(sim_tidy_categorical())
  continuous <- fit_tidy_continuous_models(sim_tidy_continuous())

  glanced <- list(
    mestimation = glance(
      ipw(binary$ps_mod, binary$outcome_mod, se_method = "mestimation")
    ),
    linearization = glance(
      ipw(binary$ps_mod, binary$outcome_mod, se_method = "linearization")
    ),
    categorical = glance(ipw(categorical$ps_mod, categorical$outcome_mod)),
    continuous = glance(ipw(continuous$ps_mod, continuous$outcome_mod))
  )

  for (route in names(glanced)) {
    expect_identical(nrow(glanced[[route]]), 1L, info = route)
    expect_identical(
      vapply(glanced[[route]], typeof, character(1)),
      glance_types(),
      info = route
    )
  }

  # Type stability is what makes the rows of several fits stackable, so stack
  # them: a route disagreeing on a column name or type would fail here.
  stacked <- vctrs::vec_rbind(!!!glanced)
  expect_identical(nrow(stacked), 4L)
  expect_named(stacked, names(glance_types()))
})

test_that("glance() rejects an argument that lands in the dots", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # The dots exist to match the generic and carry nothing. An argument absorbed
  # there would otherwise be silently ignored.
  expect_error(
    glance(res, conf.int = TRUE),
    class = "rlib_error_dots_nonempty"
  )
})
