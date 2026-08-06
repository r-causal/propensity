# Tests for the broom tidiers on `ipw` results. `tidy()` reshapes the object's
# `estimates` table into the column names the tidymodels broom developer guide
# specifies, one row per estimate, adding interval bounds only when asked for
# them. The values themselves are the object's own: the tidier renames and
# selects, and recomputes only when a requested confidence level differs from
# the one the fit stored.
#
# `glance()` describes the fit rather than its estimates: a single row naming the
# estimand and counting the observations and the residual degrees of freedom of
# the system the standard errors came from. That row carries the same columns
# with the same types on every route.
#
# `augment()` works per observation instead of per estimate: the data the fit was
# produced from, carried through in full, with the propensity score, the weights,
# the fitted values, and the residuals added as dot-prefixed columns. It reads the
# two models the result holds, so the columns describe the fit and not the frame
# they are attached to. On the default path the source frame is the outcome
# model's own frame less the `(weights)` column `model.frame()` records, which
# leaves `.weights` as the single weight column; a frame supplied through `data`
# is carried exactly as it arrives, any weight column of its own included.
#
# A result reports its effects in one of two readings, and `tidy()` reports the
# one the result records unless the call names the other through `effects`. The
# marginal reading is the table of causal contrasts above; the conditional
# reading is the outcome model's coefficient surface, one row per coefficient,
# with the standard errors of the block of the joint sandwich that accounts for
# the weights having been estimated from the same data. `glance()` and
# `augment()` describe the fit and its observations rather than its estimates,
# and neither changes with the reading.
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
# through a second fixture. `outcome_link` picks the link of a binary outcome
# model, which is what the conditional reading exponentiates on: a log link puts
# the coefficients on the log risk scale, where the default logit link puts them
# on the log odds scale, and both are scales an exponential undoes.
fit_tidy_binary_models <- function(
  dat,
  outcome_family = "binomial",
  estimand = "ate",
  outcome_link = "logit"
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
      family = quasibinomial(link = outcome_link),
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
# exposure term the continuous path requires.
fit_tidy_continuous_models <- function(dat, outcome_family = "gaussian") {
  ps_mod <- lm(A ~ x1 + x2, data = dat)
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

# A propensity score model and an outcome model fit to different rows. Each
# argument names a row to blank out in a covariate one model uses and the other
# does not, so the model using that covariate drops the row and the other keeps
# it: a row blanked for the propensity score model alone leaves the two models
# with different counts, and a different row blanked for each leaves them with
# the same count and different rows. The weights come from a propensity score
# model of the complete data, which gives the outcome model a weight for every
# row it keeps whichever rows those are.
#
# The exposure is made to alternate rather than left as the simulated draw.
# Dropping one row from a frame shifts every row after it up by one, so a
# sequence that alternates disagrees with itself at every shifted position,
# which is what makes two frames of different rows provably different rather
# than different only if the draw happened not to repeat itself across the gap.
fit_tidy_misaligned_models <- function(
  dat,
  ps_row = integer(),
  outcome_row = integer()
) {
  withr::local_seed(825)
  dat$z <- rep(c(0L, 1L), length.out = nrow(dat))
  dat$xp <- rnorm(nrow(dat))
  dat$xo <- rnorm(nrow(dat))
  dat$xp[ps_row] <- NA
  dat$xo[outcome_row] <- NA

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(glm(z ~ x1 + x2, data = dat, family = binomial()))
  )

  list(
    ps_mod = glm(z ~ x1 + xp, data = dat, family = binomial()),
    outcome_mod = glm(
      y ~ z + xo,
      data = dat,
      family = quasibinomial(),
      weights = wts
    )
  )
}

# A pair of models fit to frames whose rows were renumbered from 1, which is what
# a frame arrives with when the rows missing a value were dropped before either
# model was fit rather than by the models themselves. `drop_row` names the row
# each model's frame is missing, so the two frames carry the same row labels
# whether or not they hold the same observations: equal labels for equal rows
# when the two arguments agree, and equal labels for different rows when they do
# not. The exposure alternates for the reason the fixture above sets it to.
#
# `ps_factor` holds the exposure of the propensity score model's data as a factor
# of the labels `"0"` and `"1"`, leaving the two frames recording one exposure in
# two encodings, as a fit of a factor exposure against an outcome model of the
# numbers behind it does.
fit_tidy_renumbered_models <- function(
  dat,
  ps_drop,
  outcome_drop,
  ps_factor = FALSE
) {
  dat$z <- rep(c(0L, 1L), length.out = nrow(dat))

  renumber <- function(rows) {
    frame <- dat[-rows, ]
    rownames(frame) <- NULL
    frame
  }
  ps_data <- renumber(ps_drop)
  outcome_data <- renumber(outcome_drop)

  if (ps_factor) {
    ps_data$z <- factor(ps_data$z)
  }

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(glm(z ~ x1 + x2, data = outcome_data, family = binomial()))
  )

  list(
    ps_mod = glm(z ~ x1 + x2, data = ps_data, family = binomial()),
    outcome_mod = glm(
      y ~ z,
      data = outcome_data,
      family = quasibinomial(),
      weights = wts
    )
  )
}

# A pair of models describing the same rows in the same order under different row
# labels. The propensity score model is fit to the frame holding the missing
# values and drops those rows itself, which leaves its frame labeled around the
# gaps; the outcome model is fit to the same rows with their labels renumbered
# from 1, as a frame filtered before fitting arrives. This is the shape of an
# analysis that dropped incomplete rows in the pipeline rather than in the
# models, and there is nothing wrong with it.
fit_tidy_relabeled_models <- function(dat, na_rows) {
  withr::local_seed(825)
  dat$xp <- rnorm(nrow(dat))
  dat$xp[na_rows] <- NA

  complete <- dat[!is.na(dat$xp), ]
  rownames(complete) <- NULL

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(glm(z ~ x1 + x2, data = complete, family = binomial()))
  )

  list(
    ps_mod = glm(z ~ x1 + xp, data = dat, family = binomial()),
    outcome_mod = glm(
      y ~ z,
      data = complete,
      family = quasibinomial(),
      weights = wts
    )
  )
}

# A result assembled around a pair of models without the checks `ipw()` makes on
# the pair it is given. The estimator refuses two models fit to different rows
# and refuses an outcome model fit without weights, so a result holding either
# reaches a method only through the class's own constructor. `augment()` reads
# the two models and never the estimates table, which is why the table here is
# the shape of one rather than a fit's own.
new_tidy_ipw <- function(ps_mod, outcome_mod) {
  causalgenerics::new_ipw(
    estimand = "ate",
    wt_mod = ps_mod,
    outcome_mod = outcome_mod,
    estimates = data.frame(
      effect = "rd",
      estimate = 0,
      std.err = 0,
      z = 0,
      ci.lower = 0,
      ci.upper = 0,
      conf.level = 0.95,
      p.value = 1
    ),
    se_method = "linearization",
    fit = NULL
  )
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

# The conditional reading of a result, read off the accessors rather than off the
# tidier: the outcome model's coefficients, the standard errors the joint block
# implies, and the normal statistics those two make. Reading them from the
# accessors is what keeps the tidied table and the printed one the same numbers,
# which is also why the p-value takes the form causalgenerics prints this reading
# with, `2 * pnorm(-abs(z))`, rather than the `2 * (1 - pnorm(abs(z)))` the
# marginal estimates table was built with, a form that loses its significant
# digits in the far tail.
conditional_expectations <- function(result, conf_level = 0.95) {
  estimate <- stats::coef(result, effects = "conditional")
  std_error <- sqrt(diag(stats::vcov(result, effects = "conditional")))
  statistic <- estimate / std_error
  limits <- stats::confint(result, level = conf_level, effects = "conditional")

  list(
    term = names(estimate),
    estimate = unname(estimate),
    std.error = unname(std_error),
    statistic = unname(statistic),
    p.value = unname(2 * stats::pnorm(-abs(statistic))),
    conf.low = unname(limits[, 1L]),
    conf.high = unname(limits[, 2L])
  )
}

# Assert the whole contract of the conditional reading against the result it was
# read from: the columns of the marginal reading in the same order, so that the
# two tables stack, holding one row per coefficient of the outcome model.
expect_conditional_tidy_contract <- function(
  tidied,
  result,
  conf_int = FALSE,
  conf_level = 0.95
) {
  expected <- conditional_expectations(result, conf_level)

  expect_s3_class(tidied, c("tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_named(tidied, tidy_names(conf_int))
  expect_identical(nrow(tidied), length(expected$term))
  expect_identical(tidied$term, expected$term)
  expect_identical(tidied$estimate, expected$estimate)
  expect_identical(tidied$std.error, expected$std.error)
  expect_equal(tidied$statistic, expected$statistic, tolerance = 1e-12)
  expect_equal(tidied$p.value, expected$p.value, tolerance = 1e-12)
  if (conf_int) {
    expect_identical(tidied$conf.low, expected$conf.low)
    expect_identical(tidied$conf.high, expected$conf.high)
  }

  # A column of a tibble is addressed by its position, so the names the
  # accessors label their vectors with belong to `term` and to nothing else.
  for (nm in setdiff(names(tidied), "term")) {
    expect_null(names(tidied[[nm]]))
  }

  invisible(tidied)
}

# The glance column contract, as a named vector of the storage types the columns
# hold. Every route reports these columns, in this order, at these types. They
# are statistics of the fit: what it targets and how much information went into
# it. The propensity model, the outcome model, and the standard error method
# describe how the result was built rather than what it found, and the result
# itself already carries them under those names.
glance_types <- function() {
  c(
    estimand = "character",
    nobs = "integer",
    df.residual = "integer"
  )
}

# Assert the whole contract against the result glanced at: one row of the
# documented columns. The two counts are the caller's expectations because the
# result records neither, both being read off whichever system produced the
# standard errors. Columns are addressed with `[[` so that a column the row does
# not hold is reported by the expectation that wanted it rather than by the
# warning `$` raises on a tibble.
expect_glance_contract <- function(glanced, result, nobs, df_residual) {
  expect_s3_class(glanced, c("tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_named(glanced, names(glance_types()))
  expect_identical(
    unname(vapply(glanced, typeof, character(1))),
    unname(glance_types())
  )
  expect_identical(nrow(glanced), 1L)

  expect_identical(glanced[["estimand"]], result$estimand)
  expect_identical(glanced[["nobs"]], nobs)
  expect_identical(glanced[["df.residual"]], df_residual)

  invisible(glanced)
}

# The counts the M-estimation routes report, asserted against the deli fit they
# are read from. The stacked estimating equations, not the outcome model alone,
# are what the standard errors of those routes come from, so they are what the
# two columns describe.
expect_glance_matches_fit <- function(glanced, result) {
  expect_identical(glanced[["nobs"]], as.integer(stats::nobs(result$fit)))
  expect_identical(
    glanced[["df.residual"]],
    as.integer(stats::df.residual(result$fit))
  )

  invisible(glanced)
}

# The source frame of the default `augment()` path: the outcome model's own
# model frame less the `(weights)` column `model.frame()` records. The weights
# are already reported as `.weights`, the `psw` vector the outcome model was fit
# with, which carries its estimand and its class; keeping the model frame's plain
# numeric copy beside it would name the same weight twice and leave the frame
# ambiguous about which column is an observation's weight. A frame supplied
# through `data` is the caller's own and is carried as it arrives, so a weight
# column a user keeps in their data stays where they put it.
augment_source_frame <- function(outcome_mod) {
  mf <- stats::model.frame(outcome_mod)
  mf[names(mf) != "(weights)"]
}

# The augment column contract: every column of the source frame, in its own
# order, followed by the dot-prefixed additions. A categorical exposure widens
# `.propensity` into one column per level, so the propensity columns are the
# caller's to name. `.resid` closes the frame and appears only when the source
# holds the outcome the fit modeled.
augment_names <- function(
  source,
  propensity_cols = ".propensity",
  resid = TRUE
) {
  c(names(source), propensity_cols, ".weights", ".fitted", if (resid) ".resid")
}

# Assert the whole contract against the frame the columns were attached to. The
# dot columns are checked against the two models the result holds rather than
# against the source, because that is what they describe: `augment()` reads `x`,
# and the source only says which frame the answers are carried on.
expect_augment_contract <- function(
  augmented,
  source,
  result,
  propensity_cols = ".propensity",
  resid = TRUE
) {
  expect_s3_class(augmented, c("tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_named(augmented, augment_names(source, propensity_cols, resid))

  # one row per observation the outcome model was fit on, in the source order,
  # with every source column arriving untouched
  expect_identical(nrow(augmented), nrow(source))
  expect_identical(nrow(augmented), as.integer(stats::nobs(result$outcome_mod)))
  for (nm in names(source)) {
    expect_identical(augmented[[nm]], source[[nm]])
  }

  # the weights the outcome model was actually fit with, their class and their
  # metadata intact, rather than a numeric copy of their values
  expect_identical(
    augmented$.weights,
    stats::model.weights(stats::model.frame(result$outcome_mod))
  )

  # the outcome model's own fitted values, on the response scale
  expect_identical(
    augmented$.fitted,
    unname(stats::predict(result$outcome_mod, type = "response"))
  )

  # A column of a tibble is addressed by its position, so the observation names
  # `predict()` carries over from the model frame's row names have no business
  # here.
  added <- c(propensity_cols, ".weights", ".fitted", if (resid) ".resid")
  for (nm in added) {
    expect_null(names(augmented[[nm]]))
  }

  invisible(augmented)
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

  # The three estimate rows the fit holds summarize to this one. Every
  # observation of the fixture is used, and the parameters the stacked system
  # solves for are what the observations are spent on.
  expect_glance_contract(
    glanced,
    res,
    nobs = nrow(dat),
    df_residual = nrow(dat) - res$fit@n_params
  )
  expect_glance_matches_fit(glanced, res)
  expect_identical(glanced[["estimand"]], "ate")
})

test_that("glance() reports NA degrees of freedom without a deli fit", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # The linearization path stores no deli fit, so there is no stacked system to
  # count parameters against and the observations are the outcome model's. The
  # column stays present at its type so that the row still stacks with the rest.
  expect_null(res$fit)
  glanced <- glance(res)
  expect_glance_contract(
    glanced,
    res,
    nobs = nrow(dat),
    df_residual = NA_integer_
  )
  expect_identical(glanced[["nobs"]], as.integer(stats::nobs(res$outcome_mod)))
})

test_that("glance() summarizes a categorical fit in one row", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_tidy_categorical()
  mods <- fit_tidy_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  glanced <- glance(res)

  # Three exposure levels and six estimate rows still summarize to one row. The
  # multinomial propensity model puts more parameters into the stacked system
  # than a binary exposure does, and the residual degrees of freedom follow.
  expect_glance_contract(
    glanced,
    res,
    nobs = nrow(dat),
    df_residual = nrow(dat) - res$fit@n_params
  )
  expect_glance_matches_fit(glanced, res)
})

test_that("glance() summarizes a continuous exposure fit in one row", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_continuous()
  mods <- fit_tidy_continuous_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  glanced <- glance(res)
  expect_glance_contract(
    glanced,
    res,
    nobs = nrow(dat),
    df_residual = nrow(dat) - res$fit@n_params
  )
  expect_glance_matches_fit(glanced, res)
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
  expect_glance_contract(
    glanced,
    res,
    nobs = nrow(dat),
    df_residual = nrow(dat) - res$fit@n_params
  )
  expect_identical(glanced[["estimand"]], "att")
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

# ---- augment, binary exposure ------------------------------------------------

test_that("augment() adds the fit's own columns to the outcome model frame", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  augmented <- augment(res)

  # With no `data`, the source is the frame the outcome model was fit from,
  # which is the default broom's own `augment.lm` uses. That frame holds the
  # response, the terms of the outcome formula, and the `(weights)` column
  # `model.frame()` names. The variables of the formula are carried through; the
  # `(weights)` column is not, because the weights are reported as `.weights`,
  # and a frame naming an observation's weight twice would leave the two columns
  # to be told apart by name alone.
  mf <- stats::model.frame(mods$outcome_mod)
  expect_named(mf, c("y", "z", "(weights)"))
  expect_augment_contract(
    augmented,
    augment_source_frame(mods$outcome_mod),
    res
  )
  expect_false("(weights)" %in% names(augmented))

  expect_identical(
    augmented$.propensity,
    unname(stats::predict(mods$ps_mod, type = "response"))
  )
  expect_identical(augmented$.resid, dat$y - augmented$.fitted)

  # nothing is summarized away: every observation the outcome model saw keeps
  # its row, and the exposure column arrives in the order the data holds it
  expect_identical(nrow(augmented), 600L)
  expect_identical(augmented$z, dat$z)
})

test_that("augment() reports the propensity score the weight functions take", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  augmented <- augment(res)

  # `.propensity` is the propensity score in the sense this package's weight
  # functions use the word, so feeding the column back with the exposure and the
  # estimand the weights record returns those weights exactly. That is the
  # column's definition, and it is what makes the two columns usable together.
  rebuilt <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(augmented$.propensity, dat$z)
  )
  expect_identical(rebuilt, augmented$.weights)
})

test_that("augment() carries the weights of a non-default estimand", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat, estimand = "att")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  augmented <- augment(res)
  expect_augment_contract(
    augmented,
    augment_source_frame(mods$outcome_mod),
    res
  )

  # The weights are the ones the outcome model was fit with, so they record the
  # estimand this fit targets rather than the ATE the other fixtures use, and
  # the round trip goes through that estimand's weight function.
  expect_identical(estimand(augmented$.weights), "att")
  rebuilt <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_att(augmented$.propensity, dat$z)
  )
  expect_identical(rebuilt, augmented$.weights)
})

test_that("augment() differences a factor outcome on the scale it was fit on", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)

  # The same weighted outcome model with its 0/1 response written as the factor
  # a user is as likely to hold it in. `glm()` reads the first level as the
  # failure, so this is the same fit, and `.resid` must be the same difference:
  # the outcome on the 0/1 scale its fitted values are on, not the level codes
  # a factor subtracts as.
  dat$yf <- factor(ifelse(dat$y == 1, "yes", "no"), levels = c("no", "yes"))
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(mods$ps_mod)
  )
  outcome_mod <- glm(
    yf ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  res <- ipw(mods$ps_mod, outcome_mod, se_method = "mestimation")

  augmented <- augment(res)
  mf <- augment_source_frame(outcome_mod)
  expect_augment_contract(augmented, mf, res)

  # the response arrives as the factor the frame holds, unconverted
  expect_identical(augmented$yf, mf$yf)
  expect_identical(augmented$.resid, dat$y - augmented$.fitted)
})

test_that("augment() returns the same frame on both standard error methods", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  mestimation <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    se_method = "mestimation"
  )
  linearization <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    se_method = "linearization"
  )

  # `augment()` reads the two models, and the standard error method changes
  # neither of them, so the per-observation frame is the same one on both
  # routes even though the linearization fit stores no deli object to read.
  expect_null(linearization$fit)
  expect_augment_contract(
    augment(linearization),
    augment_source_frame(mods$outcome_mod),
    linearization
  )
  expect_identical(augment(linearization), augment(mestimation))
})

# ---- augment, categorical exposure -------------------------------------------

test_that("augment() adds one propensity column per categorical level", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_tidy_categorical()
  mods <- fit_tidy_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  augmented <- augment(res)

  # A multinomial propensity model predicts a probability for every level, so
  # there is no single column to call `.propensity`: each level gets one, named
  # for it and ordered as the model orders its levels.
  lvls <- mods$ps_mod$lev
  expect_identical(lvls, c("a", "b", "c"))
  propensity_cols <- paste0(".propensity_", lvls)
  expect_augment_contract(
    augmented,
    augment_source_frame(mods$outcome_mod),
    res,
    propensity_cols = propensity_cols
  )

  probs <- stats::predict(mods$ps_mod, type = "probs")
  for (level in lvls) {
    expect_identical(
      augmented[[paste0(".propensity_", level)]],
      unname(probs[, level])
    )
  }

  # Those K columns are the propensity score matrix the categorical weight
  # functions take, so together they rebuild the weights the outcome model was
  # fit with, exactly as the single column does for a binary exposure.
  ps_matrix <- as.matrix(augmented[, propensity_cols])
  colnames(ps_matrix) <- lvls
  rebuilt <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_matrix, dat$a, exposure_type = "categorical")
  )
  expect_identical(rebuilt, augmented$.weights)
  expect_identical(augmented$.resid, dat$y - augmented$.fitted)
})

# ---- augment, continuous exposure --------------------------------------------

test_that("augment() reports the fitted exposure for a continuous exposure", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_continuous()
  mods <- fit_tidy_continuous_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  augmented <- augment(res)
  expect_augment_contract(
    augmented,
    augment_source_frame(mods$outcome_mod),
    res
  )

  # For a continuous exposure the propensity model predicts the conditional mean
  # of the exposure, and that is what this package calls the propensity score:
  # it is the `.propensity` argument `wt_ate()` centers the conditional density
  # on. The density itself is one step further on and is not reported here, so
  # the column means the same thing on this route as on the binary one.
  expect_identical(
    augmented$.propensity,
    unname(stats::predict(mods$ps_mod, type = "response"))
  )

  # The round trip is not exact here: the weights were built from `fitted()` and
  # the column comes from `predict()`, which an `lm` computes from the design
  # matrix rather than reading back, so the two agree only to rounding.
  rebuilt <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      augmented$.propensity,
      dat$A,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )
  expect_true(is_stabilized(augmented$.weights))
  expect_equal(
    as.double(rebuilt),
    as.double(augmented$.weights),
    tolerance = 1e-8
  )

  # a gaussian outcome model reports its residuals on the outcome's own scale
  expect_identical(augmented$.resid, dat$yc - augmented$.fitted)
})

# ---- augment, the data argument ----------------------------------------------

test_that("augment(data = ) returns that data with the fit's columns added", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  augmented <- augment(res, data = dat)

  # broom's convention: `data` is the data the fit was produced from, carried
  # through in full so that covariates the outcome formula left out sit beside
  # the fit's own columns. The outcome model frame holds `y`, `z`, and the
  # weights alone, so `x1` and `x2` arrive only this way.
  expect_augment_contract(augmented, dat, res)
  expect_identical(augmented$x1, dat$x1)
  expect_identical(augmented$x2, dat$x2)

  # The dot columns describe the fit, so which frame they were attached to does
  # not change one of them.
  default <- augment(res)
  for (nm in c(".propensity", ".weights", ".fitted", ".resid")) {
    expect_identical(augmented[[nm]], default[[nm]])
  }
})

test_that("augment(data = ) keeps a weight column the supplied data holds", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # The default path leaves `(weights)` out of the frame it builds for itself,
  # which is a decision about that frame and about no other. A frame supplied
  # through `data` is the caller's own, so a column of theirs survives whatever
  # it happens to be named, including the name the default path declines to
  # carry. Dropping it would delete a user's data for matching a name.
  user_weights <- 2 *
    as.double(stats::model.weights(stats::model.frame(mods$outcome_mod)))
  user_frame <- dat
  user_frame[["(weights)"]] <- user_weights

  augmented <- augment(res, data = user_frame)
  expect_augment_contract(augmented, user_frame, res)

  # The user's column arrives verbatim, as the scaled doubles they held rather
  # than as a copy of `.weights`, and the two sit side by side.
  expect_identical(augmented[["(weights)"]], user_weights)
  expect_s3_class(augmented$.weights, "psw")
  expect_identical(estimand(augmented$.weights), "ate")
})

test_that("augment(data = ) drops .resid when the data omits the outcome", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  covariates <- dat[c("x1", "x2", "z")]
  augmented <- augment(res, data = covariates)

  # A residual is an observed outcome minus a fitted value, and this frame holds
  # no outcome to subtract from, so the column is absent rather than missing.
  expect_augment_contract(augmented, covariates, res, resid = FALSE)
  expect_false(".resid" %in% names(augmented))
})

test_that("augment(data = ) rejects data that is not one row per observation", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # `data` names the frame the fit's own columns are carried on, so a frame of
  # another length has no row to carry them on. Recycling or truncating would
  # attach one observation's answers to another observation's covariates.
  expect_error(
    augment(res, data = dat[seq_len(100), ]),
    class = "propensity_augment_data_error"
  )
})

test_that("augment(data = ) rejects data holding a column augment adds", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # `.fitted` is one of the columns the fit's answers are reported in, so a
  # frame that already holds a column of that name has nowhere for the fit's own
  # to go. Adding it beside the caller's would return a frame naming two columns
  # `.fitted`, and writing over it would delete their data.
  collide <- dat
  collide$.fitted <- 0

  expect_error(
    augment(res, data = collide),
    class = "propensity_augment_column_error"
  )

  # the refusal names the column in the way and what to do about it
  expect_propensity_error(augment(res, data = collide))
})

test_that("augment(data = ) names every column of the data it clashes with", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # A frame can clash on more than one column, and a report of one of them
  # leaves the caller to rename that column and meet the next refusal. Every
  # column that clashes is named at once. `.propensity` is among them here: a
  # binary or continuous exposure reports its propensity score in the one column
  # of that name, which is as much in the way as the other three.
  collide <- dat
  collide$.propensity <- 0.5
  collide$.weights <- 1
  collide$.fitted <- 0

  err <- expect_error(
    augment(res, data = collide),
    class = "propensity_augment_column_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, ".propensity", fixed = TRUE)
  expect_match(msg, ".weights", fixed = TRUE)
  expect_match(msg, ".fitted", fixed = TRUE)

  # the sentence the three names are reported in
  expect_propensity_error(augment(res, data = collide))
})

test_that("augment(data = ) counts a .resid column only when it adds one", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # `.resid` is added only to a frame holding the outcome to difference, so
  # whether a `.resid` column of the caller's is in the way depends on the frame
  # it arrives in. A frame with no outcome column is given no residual, and the
  # caller's column of that name goes through as any other column of theirs
  # would.
  covariates <- dat[c("x1", "x2", "z")]
  covariates$.resid <- 0
  expect_augment_contract(
    augment(res, data = covariates),
    covariates,
    res,
    resid = FALSE
  )

  # The same column in a frame that does hold the outcome is a column the
  # residual would be written over.
  with_outcome <- dat
  with_outcome$.resid <- 0

  expect_error(
    augment(res, data = with_outcome),
    class = "propensity_augment_column_error"
  )
})

test_that("augment(data = ) clashes on the level columns of a categorical fit", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_tidy_categorical()
  mods <- fit_tidy_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  # A categorical exposure widens `.propensity` into one column per level, so
  # the propensity columns in the way of this fit are those, named for the
  # levels the propensity score model holds.
  collide <- dat
  collide$.propensity_b <- 0

  expect_error(
    augment(res, data = collide),
    class = "propensity_augment_column_error"
  )

  # the refusal names the level column rather than `.propensity`
  expect_propensity_error(augment(res, data = collide))

  # `.propensity` is not one of them on this route. No column of that name is
  # added, so a column of the caller's carrying it is theirs to keep, which is
  # what makes the clash a question about the columns this fit adds rather than
  # about a fixed list of names.
  kept <- dat
  kept$.propensity <- 0
  expect_augment_contract(
    augment(res, data = kept),
    kept,
    res,
    propensity_cols = paste0(".propensity_", mods$ps_mod$lev)
  )
})

test_that("augment() rejects an argument that lands in the dots", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # The dots exist to match the generic. `newdata` is a broom argument this
  # method does not take, and absorbing it there would return the fit's own rows
  # as though they were predictions for the frame supplied.
  expect_error(
    augment(res, newdata = dat),
    class = "rlib_error_dots_nonempty"
  )
})

# ---- augment, results ipw() would not have built ------------------------------

# `augment()` reads the two models a result holds and reports a column of each
# beside the other, which is an answer per observation only while the two models
# describe the same observations. `ipw()` checks that of the pair it is handed,
# so the pairs here are assembled into results directly: two fit to different
# rows, one fit to the same rows that the method has no reason to refuse, and
# one whose outcome model carries no weights at all.

test_that("augment() rejects a result whose models saw different counts", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_misaligned_models(dat, ps_row = 7)
  res <- new_tidy_ipw(mods$ps_mod, mods$outcome_mod)

  # The propensity score model dropped the row its covariate was missing on and
  # the outcome model kept every row, so from the seventh row on the two models
  # are describing different observations. Carrying the shorter column beside
  # the longer ones would report one observation's propensity score against
  # another's fitted value.
  expect_identical(as.integer(stats::nobs(mods$ps_mod)), 599L)
  expect_identical(nrow(stats::model.frame(mods$outcome_mod)), 600L)

  err <- expect_error(
    augment(res),
    class = "propensity_augment_alignment_error"
  )

  # Both counts are reported: which model dropped rows is the first thing a
  # reader needs, and the difference between the two is what says so.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "\\b599\\b")
  expect_match(msg, "\\b600\\b")

  # the sentence the two counts are reported in
  expect_propensity_error(augment(res))
})

test_that("augment() rejects a result whose models saw different rows", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_misaligned_models(dat, ps_row = 7, outcome_row = 12)
  res <- new_tidy_ipw(mods$ps_mod, mods$outcome_mod)

  # Each model dropped one row and they dropped different ones, so the counts
  # agree and the observations do not. A count is all that a column of the right
  # length proves; what says whether the observation in a row of one frame is
  # the observation in that row of the other is the data, and the exposure is
  # the variable both frames hold.
  ps_frame <- stats::model.frame(mods$ps_mod)
  outcome_frame <- stats::model.frame(mods$outcome_mod)
  expect_identical(nrow(ps_frame), nrow(outcome_frame))
  expect_false(identical(ps_frame$z, outcome_frame$z))

  expect_error(
    augment(res),
    class = "propensity_augment_alignment_error"
  )

  # the refusal reports observations that agree in number and not in identity
  expect_propensity_error(augment(res))
})

test_that("augment() rejects models of different rows carrying equal labels", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_renumbered_models(dat, ps_drop = 7, outcome_drop = 12)
  res <- new_tidy_ipw(mods$ps_mod, mods$outcome_mod)

  # Both frames were renumbered from 1 before either model was fit, so the two
  # carry exactly the same row labels while holding different observations. A
  # label records where a row came from and survives nothing that renumbers it,
  # so agreeing labels are no evidence of agreeing rows and this fit has to be
  # caught by the data itself.
  ps_frame <- stats::model.frame(mods$ps_mod)
  outcome_frame <- stats::model.frame(mods$outcome_mod)
  expect_identical(rownames(ps_frame), rownames(outcome_frame))
  expect_identical(rownames(outcome_frame), as.character(1:599))
  expect_false(identical(ps_frame$z, outcome_frame$z))

  expect_error(
    augment(res),
    class = "propensity_augment_alignment_error"
  )
})

test_that("augment() carries a factor exposure against its numeric encoding", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  factor_dat <- dat
  factor_dat$z <- factor(dat$z)

  # The propensity score model was fit to the exposure as a factor and the
  # outcome model to the 0/1 numbers those labels spell, which is a pair
  # `ipw()` builds a result from. The two frames record one exposure in two
  # encodings and describe the same rows, so the columns read from either belong
  # beside the columns read from the other, and a comparison that read `"0"`
  # against 0 as a disagreement would refuse a fit with nothing wrong with it.
  ps_mod <- glm(z ~ x1 + x2, data = factor_dat, family = binomial())
  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_mod))
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  res <- ipw(ps_mod, outcome_mod)

  augmented <- augment(res)
  expect_augment_contract(augmented, augment_source_frame(outcome_mod), res)

  # The encoding of the exposure the propensity score model was fit to is not
  # something the reported columns depend on, so the frame is the one the same
  # fit of the numeric exposure returns, column for column.
  numeric_mods <- fit_tidy_binary_models(dat)
  numeric_res <- ipw(numeric_mods$ps_mod, numeric_mods$outcome_mod)
  expect_identical(augmented, augment(numeric_res))
})

test_that("augment() carries an exposure recoded past what a number reads", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  labeled <- dat
  labeled$z <- factor(ifelse(dat$z == 1, "b", "a"), levels = c("a", "b"))

  # Here the two encodings are a factor of `"a"` and `"b"` against the 0 and 1 a
  # user recoded it to. The labels spell no numbers, so nothing lines the two up
  # position by position: their values differ at every position while their
  # observations agree at all of them. A comparison that cannot tell those apart
  # has proven nothing, and proving nothing is not a reason to refuse.
  ps_mod <- glm(z ~ x1 + x2, data = labeled, family = binomial())
  ps_labels <- as.character(stats::model.frame(ps_mod)$z)
  expect_true(anyNA(suppressWarnings(as.numeric(ps_labels))))

  wts <- withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_mod))
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  res <- ipw(ps_mod, outcome_mod)

  augmented <- augment(res)
  expect_augment_contract(augmented, augment_source_frame(outcome_mod), res)
  expect_identical(
    augmented$.propensity,
    unname(stats::predict(ps_mod, type = "response"))
  )
})

test_that("augment() rejects a factor exposure of rows the numbers deny", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_renumbered_models(
    dat,
    ps_drop = 7,
    outcome_drop = 12,
    ps_factor = TRUE
  )
  res <- new_tidy_ipw(mods$ps_mod, mods$outcome_mod)

  # The same two encodings as the fit above, over rows that are not the same
  # rows. Reading the labels as the numbers they spell is what makes the two
  # frames comparable, and once they are comparable the alternating exposure
  # shows the shift between them. An encoding the values can be read through is
  # no reason to stop looking at the values.
  ps_frame <- stats::model.frame(mods$ps_mod)
  outcome_frame <- stats::model.frame(mods$outcome_mod)
  expect_s3_class(ps_frame$z, "factor")
  expect_identical(rownames(ps_frame), rownames(outcome_frame))
  expect_false(identical(
    as.numeric(as.character(ps_frame$z)),
    as.numeric(outcome_frame$z)
  ))

  expect_error(
    augment(res),
    class = "propensity_augment_alignment_error"
  )
})

test_that("augment() carries models of the same rows carrying different labels", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_relabeled_models(dat, na_rows = c(7, 12, 400))
  res <- new_tidy_ipw(mods$ps_mod, mods$outcome_mod)

  # The propensity score model dropped the three incomplete rows itself and is
  # labeled around the gaps; the outcome model was fit to the same rows after
  # they had been renumbered from 1. The two describe the same 597 observations
  # in the same order, which is what makes the columns of one belong beside the
  # columns of the other, and their labels disagree throughout, which is nothing.
  ps_frame <- stats::model.frame(mods$ps_mod)
  outcome_frame <- stats::model.frame(mods$outcome_mod)
  expect_identical(nrow(ps_frame), 597L)
  expect_identical(rownames(outcome_frame), as.character(1:597))
  expect_false(identical(rownames(ps_frame), rownames(outcome_frame)))
  expect_identical(ps_frame$z, outcome_frame$z)

  augmented <- augment(res)
  expect_augment_contract(
    augmented,
    augment_source_frame(mods$outcome_mod),
    res
  )
  expect_identical(nrow(augmented), 597L)
  expect_identical(
    augmented$.propensity,
    unname(stats::predict(mods$ps_mod, type = "response"))
  )
})

test_that("augment() carries a fit whose two models dropped the same row", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_misaligned_models(dat, ps_row = 7, outcome_row = 7)
  res <- new_tidy_ipw(mods$ps_mod, mods$outcome_mod)

  # Both models were missing a covariate on the same row, so both dropped it and
  # both describe the same 599 observations. The rows they kept are numbered
  # around the gap that leaves, which is what a fit of incomplete data looks
  # like and is not a reason to refuse anything: what makes two models
  # comparable is that they kept the same rows, not that those rows run from one
  # to the number of them.
  outcome_rows <- rownames(stats::model.frame(mods$outcome_mod))
  expect_identical(outcome_rows, as.character(c(1:6, 8:600)))
  expect_identical(
    names(stats::predict(mods$ps_mod, type = "response")),
    outcome_rows
  )

  augmented <- augment(res)
  expect_augment_contract(
    augmented,
    augment_source_frame(mods$outcome_mod),
    res
  )
  expect_identical(nrow(augmented), 599L)
  expect_identical(
    augmented$.propensity,
    unname(stats::predict(mods$ps_mod, type = "response"))
  )
})

test_that("augment() rejects a result whose outcome model has no weights", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  unweighted <- glm(y ~ z, data = dat, family = quasibinomial())
  res <- new_tidy_ipw(mods$ps_mod, unweighted)

  # The outcome model of an IPW result is a weighted one by construction, so a
  # result reporting an outcome model with no weights to read is a result built
  # around something else. `.weights` has nothing to hold, and a frame of the
  # remaining columns would describe a fit the object does not name.
  expect_null(stats::model.weights(stats::model.frame(unweighted)))

  # The missing weights are the same fact `ipw()` refuses the pair for, so they
  # are reported under the class that names that fact. Which method noticed it
  # is not what a condition class is for.
  expect_error(
    augment(res),
    class = "propensity_ipw_weights_missing_error"
  )

  # the refusal says what an ipw() outcome model is fit with
  expect_propensity_error(augment(res))
})

# ---- the presentation mode ---------------------------------------------------

# `tidy()` reports the reading the result records unless `effects` names the
# other one for the call. The conditional reading is a coefficient table of the
# outcome model, which makes it the reading a `tidy()` method of any other model
# would return, so it takes the same columns in the same order as the marginal
# reading and the two stack. Its standard errors are the ones the joint
# estimation of the weights and the outcome implies, read through the accessors
# that report them, rather than the ones the outcome model computed for itself
# while treating the weights it was fit with as fixed.

test_that("tidy() reports the outcome model's coefficients conditionally", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  tidied <- tidy(res, effects = "conditional")

  # One row per coefficient of the outcome model in place of the three effect
  # measures the marginal reading of the same result reports, in the columns that
  # reading uses.
  expect_conditional_tidy_contract(tidied, res)
  expect_identical(tidied$term, names(coef(mods$outcome_mod)))
  expect_identical(tidied$term, c("(Intercept)", "z"))
  expect_identical(tidied$estimate, unname(coef(mods$outcome_mod)))
  expect_identical(tidy(res)$term, c("rd", "log(rr)", "log(or)"))

  expect_conditional_tidy_contract(
    tidy(res, conf.int = TRUE, effects = "conditional"),
    res,
    conf_int = TRUE
  )
})

test_that("tidy() reports the corrected standard errors conditionally", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  tidied <- tidy(res, conf.int = TRUE, effects = "conditional")

  # The standard errors are the diagonal of the block the stacked system leaves
  # for the outcome model, which is what the accessor reports and what the
  # printed reading of the same result tabulates.
  expect_identical(
    tidied$std.error,
    unname(sqrt(diag(vcov(res, effects = "conditional"))))
  )

  # That block is not the covariance the outcome model computed for itself: that
  # one treats the estimated weights as fixed and reports an uncertainty the
  # coefficients do not have. The two are far enough apart here that a tidier
  # reading the model instead of the result would be caught by the difference
  # rather than by rounding.
  naive <- unname(sqrt(diag(vcov(mods$outcome_mod))))
  expect_false(isTRUE(all.equal(tidied$std.error, naive, tolerance = 1e-3)))

  # Everything built from the standard error follows it there.
  expect_false(isTRUE(all.equal(
    tidied$statistic,
    tidied$estimate / naive,
    tolerance = 1e-3
  )))
  expect_false(isTRUE(all.equal(
    tidied$conf.low,
    tidied$estimate - qnorm(0.975) * naive,
    tolerance = 1e-3
  )))
})

test_that("tidy() follows the reading the result records", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
  built <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    se_method = "mestimation",
    effects = "conditional"
  )

  # `effects = NULL`, the default, is the reading the result itself records, so a
  # result built in the conditional one tidies to that table with nothing named
  # at the call site, and so does a result moved there afterwards.
  expect_conditional_tidy_contract(tidy(built), built)
  expect_identical(tidy(built), tidy(res, effects = "conditional"))
  expect_identical(
    tidy(as_conditional(res)),
    tidy(res, effects = "conditional")
  )

  # Naming a reading at the call site answers in it, in both directions, which is
  # what keeps the marginal table reachable from a result that presents the other
  # one.
  expect_identical(tidy(as_conditional(res), effects = "marginal"), tidy(res))
  expect_identical(
    tidy(as_conditional(res), conf.int = TRUE, effects = "marginal"),
    tidy(res, conf.int = TRUE)
  )
  expect_identical(tidy(res, effects = "marginal"), tidy(res))

  # Naming one also leaves the result where it was, so the next call with nothing
  # named answers in the stored reading.
  expect_identical(res$effects, "marginal")
})

test_that("tidy() rebuilds the conditional interval at the level asked for", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  at_default <- tidy(res, conf.int = TRUE, effects = "conditional")
  at_90 <- tidy(res, conf.int = TRUE, conf.level = 0.9, effects = "conditional")
  expect_conditional_tidy_contract(
    at_90,
    res,
    conf_int = TRUE,
    conf_level = 0.9
  )

  # The bounds the result stores describe the effects of the other reading, so
  # there are none to serve here and the limits are built at whatever level the
  # call asked for.
  half <- qnorm(0.95) * at_90$std.error
  expect_equal(at_90$conf.low, at_90$estimate - half, tolerance = 1e-10)
  expect_equal(at_90$conf.high, at_90$estimate + half, tolerance = 1e-10)
  expect_true(all(at_90$conf.low > at_default$conf.low))
  expect_true(all(at_90$conf.high < at_default$conf.high))

  # The level governs bounds that were asked for and nothing else, and it is
  # validated in this reading as in the other one.
  expect_identical(
    tidy(res, conf.level = 0.9, effects = "conditional"),
    tidy(res, effects = "conditional")
  )
  expect_error(
    tidy(res, conf.int = TRUE, conf.level = 1, effects = "conditional"),
    class = "propensity_conf_level_error"
  )
})

test_that("tidy() has no conditional reading of a linearization fit", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # Nothing was solved jointly here, so there is no corrected block for the
  # standard errors of a coefficient table to come from. A tidied row is its
  # estimate and its inference together, so the absence is reported rather than
  # filled in with the covariance the outcome model computed for itself, and it
  # is reported whether or not the call asked for bounds.
  expect_null(res$fit)
  expect_error(
    tidy(res, effects = "conditional"),
    class = "causalgenerics_no_conditional_vcov"
  )
  expect_error(
    tidy(res, conf.int = TRUE, effects = "conditional"),
    class = "causalgenerics_no_conditional_vcov"
  )
  expect_error(
    tidy(as_conditional(res)),
    class = "causalgenerics_no_conditional_vcov"
  )

  # The marginal reading of the same result is untouched by the absence, whether
  # it is the reading the result records or the one the call names.
  expect_tidy_contract(tidy(res), res)
  expect_tidy_contract(tidy(as_conditional(res), effects = "marginal"), res)
})

test_that("tidy() reports one conditional row per categorical coefficient", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- sim_tidy_categorical()
  mods <- fit_tidy_categorical_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  tidied <- tidy(res, conf.int = TRUE, effects = "conditional")
  expect_conditional_tidy_contract(tidied, res, conf_int = TRUE)
  expect_identical(tidied$term, c("(Intercept)", "ab", "ac"))

  # The comparison column belongs to the marginal reading, which reports every
  # effect measure once per contrast of levels. A coefficient is not such a
  # contrast, so the conditional table carries no such column and reports three
  # rows where the marginal reading of the same result reports six.
  expect_named(tidied, tidy_names(conf_int = TRUE))
  expect_identical(nrow(tidied), 3L)
  expect_identical(nrow(tidy(res, conf.int = TRUE)), 6L)
})

test_that("tidy() reports the conditional reading of a continuous fit", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_continuous()
  mods <- fit_tidy_continuous_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  tidied <- tidy(res, conf.int = TRUE, effects = "conditional")
  expect_conditional_tidy_contract(tidied, res, conf_int = TRUE)

  # The marginal reading of a continuous exposure is the one slope row; the
  # conditional one is the marginal structural model it came from, intercept
  # included.
  expect_identical(tidied$term, c("(Intercept)", "A"))
  expect_identical(tidy(res)$term, "slope")
  expect_false(isTRUE(all.equal(
    tidied$std.error,
    unname(sqrt(diag(vcov(mods$outcome_mod)))),
    tolerance = 1e-3
  )))

  # The outcome model of this fixture has an identity link, so its coefficients
  # are on the outcome's own scale and there is nothing for an exponential to
  # undo.
  expect_identical(family(mods$outcome_mod)$link, "identity")
  expect_error(
    tidy(res, exponentiate = TRUE, effects = "conditional"),
    class = "propensity_exponentiate_error"
  )
})

# The conditional reading exponentiates by the link of the outcome model rather
# than by the labels of the rows, which is the one place the two readings differ
# in what `exponentiate` means. The marginal reading exponentiates the rows it
# labels as ratios and relabels them; a coefficient table has no such labels, and
# every coefficient of a log odds or log risk model is on the scale an
# exponential undoes, the intercept with the rest.
#
# What is done to the other columns is broom's convention, which
# `broom:::exponentiate()` implements for `tidy.lm()`, `tidy.glm()`, and every
# other method that takes the argument: the estimate and, when they were asked
# for, the interval bounds are exponentiated, and the standard error, the
# statistic, and the p-value are left describing the link scale. There is no
# delta-method rescaling of the standard error and no renaming of a term.

test_that("tidy(exponentiate = TRUE) exponentiates every conditional row", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  expect_identical(family(mods$outcome_mod)$link, "logit")
  plain <- tidy(res, conf.int = TRUE, effects = "conditional")
  expo <- tidy(
    res,
    conf.int = TRUE,
    exponentiate = TRUE,
    effects = "conditional"
  )

  # The intercept is a log odds here and is exponentiated with the rest, where
  # the marginal reading of the same result leaves its risk difference row alone.
  expect_named(expo, tidy_names(conf_int = TRUE))
  expect_identical(expo$term, plain$term)
  expect_equal(expo$estimate, exp(plain$estimate), tolerance = 1e-12)
  expect_equal(expo$conf.low, exp(plain$conf.low), tolerance = 1e-12)
  expect_equal(expo$conf.high, exp(plain$conf.high), tolerance = 1e-12)

  # The inference columns describe the link scale estimate and stay there.
  expect_identical(expo$std.error, plain$std.error)
  expect_identical(expo$statistic, plain$statistic)
  expect_identical(expo$p.value, plain$p.value)
  expect_false(isTRUE(all.equal(
    expo$std.error,
    expo$estimate * plain$std.error,
    tolerance = 1e-6
  )))

  # The transformation does not depend on bounds having been asked for, so a call
  # that asked for none still returns the exponentiated estimates.
  bare <- tidy(res, exponentiate = TRUE, effects = "conditional")
  expect_named(bare, tidy_names())
  expect_identical(bare$estimate, expo$estimate)
  expect_identical(bare$term, plain$term)

  # The policy belongs to the reading rather than to the argument that named it,
  # so a result that records the conditional reading exponentiates the same way.
  expect_identical(
    tidy(as_conditional(res), conf.int = TRUE, exponentiate = TRUE),
    expo
  )

  # The marginal reading of the same result keeps the policy it had: the ratio
  # rows exponentiated and relabeled, the risk difference row untouched.
  expect_identical(
    tidy(as_conditional(res), exponentiate = TRUE, effects = "marginal"),
    tidy(res, exponentiate = TRUE)
  )
  expect_identical(tidy(res, exponentiate = TRUE)$term, c("rd", "rr", "or"))
})

test_that("tidy(exponentiate = TRUE) exponentiates a log link conditionally", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat, outcome_link = "log")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # A log link outcome model reports log risks rather than log odds, and the
  # policy follows the link rather than the effect measure a reading of the
  # result happens to name, so this table exponentiates exactly as the logit one
  # does.
  expect_identical(family(mods$outcome_mod)$link, "log")
  plain <- tidy(res, conf.int = TRUE, effects = "conditional")
  expo <- tidy(
    res,
    conf.int = TRUE,
    exponentiate = TRUE,
    effects = "conditional"
  )

  expect_identical(expo$term, plain$term)
  expect_equal(expo$estimate, exp(plain$estimate), tolerance = 1e-12)
  expect_equal(expo$conf.low, exp(plain$conf.low), tolerance = 1e-12)
  expect_equal(expo$conf.high, exp(plain$conf.high), tolerance = 1e-12)
  expect_identical(expo$std.error, plain$std.error)
  expect_identical(expo$statistic, plain$statistic)
  expect_identical(expo$p.value, plain$p.value)
})

test_that("tidy(exponentiate = TRUE) rejects an identity link conditionally", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # The coefficients of an identity link model are on the outcome's own scale,
  # and an exponential of one of them is a number that describes nothing. The
  # request is refused rather than answered with it, whichever way the reading
  # was arrived at and whether or not bounds were asked for.
  expect_identical(family(mods$outcome_mod)$link, "identity")
  expect_error(
    tidy(res, exponentiate = TRUE, effects = "conditional"),
    class = "propensity_exponentiate_error"
  )
  expect_error(
    tidy(res, conf.int = TRUE, exponentiate = TRUE, effects = "conditional"),
    class = "propensity_exponentiate_error"
  )
  expect_error(
    tidy(as_conditional(res), exponentiate = TRUE),
    class = "propensity_exponentiate_error"
  )

  # There is nothing to refuse without the argument, so the reading itself is
  # reported as it is on any other link.
  expect_conditional_tidy_contract(
    tidy(res, conf.int = TRUE, effects = "conditional"),
    res,
    conf_int = TRUE
  )

  # The refusal belongs to the conditional reading. The marginal reading of the
  # same result reports a difference in means, which is a row this argument has
  # always left alone rather than one it refuses.
  expect_identical(tidy(res, exponentiate = TRUE), tidy(res))
  expect_identical(
    tidy(as_conditional(res), exponentiate = TRUE, effects = "marginal"),
    tidy(res)
  )
})

test_that("tidy(exponentiate = TRUE) rejects the other binary links", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()

  # A probit and a cloglog link are the other two links a binary outcome model
  # of this estimator can be fit with. Both put their coefficients on a scale of
  # their own rather than on a log one, so both are refused for the same reason
  # the identity link is: the policy is an allowlist of the two links an
  # exponential undoes rather than a list of the ones it does not.
  for (link in c("probit", "cloglog")) {
    mods <- fit_tidy_binary_models(dat, outcome_link = link)
    res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")
    expect_identical(family(mods$outcome_mod)$link, link)

    expect_error(
      tidy(res, exponentiate = TRUE, effects = "conditional"),
      class = "propensity_exponentiate_error"
    )
    expect_error(
      tidy(res, conf.int = TRUE, exponentiate = TRUE, effects = "conditional"),
      class = "propensity_exponentiate_error"
    )
    expect_error(
      tidy(as_conditional(res), exponentiate = TRUE),
      class = "propensity_exponentiate_error"
    )

    # There is nothing to refuse without the argument, so the reading itself is
    # reported as it is on any other link.
    expect_conditional_tidy_contract(
      tidy(res, conf.int = TRUE, effects = "conditional"),
      res,
      conf_int = TRUE
    )

    # The refusal belongs to the conditional reading. The marginal reading of
    # the same result reports ratio measures whatever the link of the model they
    # were computed from, and exponentiates them as it always has.
    expect_identical(tidy(res, exponentiate = TRUE)$term, c("rd", "rr", "or"))
  }
})

test_that("tidy() rejects an effects value naming neither reading", {
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "linearization")

  # The argument belongs to the accessors, which resolve it and reject a value
  # that names neither reading. Validating it here as well would leave two
  # validators to keep in step, so the value goes to theirs and the condition
  # they raise is the one that reaches the caller.
  expect_error(
    tidy(res, effects = "banana"),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    tidy(res, effects = c("marginal", "conditional")),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    tidy(res, effects = 1),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    tidy(res, conf.int = TRUE, effects = NA_character_),
    class = "causalgenerics_invalid_argument_effects"
  )

  # The value is rejected before the reading it names is looked for, so a result
  # that presents one reading answers a nonsense request the same way.
  expect_error(
    tidy(as_conditional(res), effects = "banana"),
    class = "causalgenerics_invalid_argument_effects"
  )

  # `effects` follows the dots, where an abbreviation of it does not match. A
  # partial name is an argument left in the dots rather than a silent request for
  # a reading.
  expect_error(
    tidy(res, effect = "conditional"),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("glance() reports the same row in either reading", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # `glance()` describes the fit rather than its estimates: what it targets and
  # how much information went into it, neither of which is a reading of those
  # estimates. The row is the same row whichever reading the result presents, and
  # the reading is not an argument this method takes.
  expect_identical(glance(as_conditional(res)), glance(res))
  expect_identical(
    glance(ipw(
      mods$ps_mod,
      mods$outcome_mod,
      se_method = "mestimation",
      effects = "conditional"
    )),
    glance(res)
  )
  expect_error(
    glance(res, effects = "conditional"),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("augment() returns the same frame in either reading", {
  skip_if_not_installed("deli")
  dat <- sim_tidy_binary()
  mods <- fit_tidy_binary_models(dat)
  res <- ipw(mods$ps_mod, mods$outcome_mod, se_method = "mestimation")

  # `augment()` works per observation and reads the two models the result holds,
  # neither of which the reading changes, so the frame is the same on either
  # reading, whichever source frame it is carried on.
  expect_identical(augment(as_conditional(res)), augment(res))
  expect_identical(
    augment(as_conditional(res), data = dat),
    augment(res, data = dat)
  )
  expect_error(
    augment(res, effects = "conditional"),
    class = "rlib_error_dots_nonempty"
  )
})
