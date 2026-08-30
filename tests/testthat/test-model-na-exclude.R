# Model routes under `na.action = na.exclude`
#
# A model fit with `na.action = na.exclude` records only the rows it analyzed
# but pads everything it predicts back to the length of the data it was given,
# with a missing value where a row was dropped. The model routes read the scores
# from the padded side and the exposure from the recorded side, so the two have
# to be brought back to one length before anything is computed with them.
#
# These tests pin the reading `na.exclude` asks for: weights as long as the data
# the model was given, missing at exactly the rows it dropped and equal to the
# complete-case fit everywhere else. That is what the routes already do when the
# caller supplies the exposure themselves, and it is the reading `ipw()` gets
# without `.data`, so it is the one the whole family should agree on.

na_exclude_rows <- c(3L, 17L, 40L, 88L, 150L)

sim_na_exclude <- function(seed = 2026, n = 300) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  bin <- rbinom(n, 1, plogis(0.6 * x1 + 0.2 * x2))
  dose <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  eta_b <- -0.2 + 0.7 * x1
  eta_c <- 0.1 - 0.5 * x1 + 0.4 * x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  u <- runif(n)
  lab <- ifelse(
    u < 1 / denom,
    "a",
    ifelse(u < (1 + exp(eta_b)) / denom, "b", "c")
  )
  cat3 <- factor(lab, levels = c("a", "b", "c"))
  y <- rbinom(n, 1, plogis(-0.3 + 0.8 * bin + 0.5 * x1))
  dat <- data.frame(x1, x2, bin, dose, cat3, y)
  # Only a covariate goes missing: the exposure and the outcome are complete, so
  # every route drops the same rows and for the same reason.
  dat$x1[na_exclude_rows] <- NA
  dat
}

complete_rows <- function(dat) {
  dat[stats::complete.cases(dat), , drop = FALSE]
}

# ---- the weight routes ------------------------------------------------------

test_that("wt_ate() pads a binary glm fit with na.exclude", {
  dat <- sim_na_exclude()
  fit_padded <- glm(
    bin ~ x1 + x2,
    data = dat,
    family = binomial(),
    na.action = na.exclude
  )
  fit_complete <- glm(
    bin ~ x1 + x2,
    data = complete_rows(dat),
    family = binomial()
  )

  wts <- wt_ate(fit_padded)

  expect_length(wts, nrow(dat))
  expect_equal(which(is.na(as.numeric(wts))), na_exclude_rows)
  expect_equal(
    as.numeric(wts)[-na_exclude_rows],
    as.numeric(wt_ate(fit_complete))
  )
})

test_that("wt_att() pads a binary glm fit with na.exclude", {
  dat <- sim_na_exclude()
  fit_padded <- glm(
    bin ~ x1 + x2,
    data = dat,
    family = binomial(),
    na.action = na.exclude
  )
  fit_complete <- glm(
    bin ~ x1 + x2,
    data = complete_rows(dat),
    family = binomial()
  )

  wts <- wt_att(fit_padded)

  expect_length(wts, nrow(dat))
  expect_equal(which(is.na(as.numeric(wts))), na_exclude_rows)
  expect_equal(
    as.numeric(wts)[-na_exclude_rows],
    as.numeric(wt_att(fit_complete))
  )
})

test_that("wt_ate() pads a continuous lm fit with na.exclude", {
  dat <- sim_na_exclude()
  fit_padded <- lm(dose ~ x1 + x2, data = dat, na.action = na.exclude)
  fit_complete <- lm(dose ~ x1 + x2, data = complete_rows(dat))

  wts <- wt_ate(fit_padded, exposure_type = "continuous")

  expect_length(wts, nrow(dat))
  expect_equal(which(is.na(as.numeric(wts))), na_exclude_rows)
  expect_equal(
    as.numeric(wts)[-na_exclude_rows],
    as.numeric(wt_ate(fit_complete, exposure_type = "continuous"))
  )
})

test_that("wt_ate() pads a multinom fit with na.exclude", {
  skip_if_not_installed("nnet")
  dat <- sim_na_exclude()
  fit_padded <- nnet::multinom(
    cat3 ~ x1 + x2,
    data = dat,
    trace = FALSE,
    na.action = na.exclude
  )
  fit_complete <- nnet::multinom(
    cat3 ~ x1 + x2,
    data = complete_rows(dat),
    trace = FALSE
  )

  wts <- wt_ate(fit_padded, exposure_type = "categorical")

  expect_length(wts, nrow(dat))
  expect_equal(which(is.na(as.numeric(wts))), na_exclude_rows)
  expect_equal(
    as.numeric(wts)[-na_exclude_rows],
    as.numeric(wt_ate(fit_complete, exposure_type = "categorical")),
    tolerance = 1e-6
  )
})

# ---- the modifier routes ----------------------------------------------------

test_that("ps_trim() pads a fit with na.exclude on a method that reads the exposure", {
  dat <- sim_na_exclude()
  fit_padded <- glm(
    bin ~ x1 + x2,
    data = dat,
    family = binomial(),
    na.action = na.exclude
  )
  fit_complete <- glm(
    bin ~ x1 + x2,
    data = complete_rows(dat),
    family = binomial()
  )

  trimmed <- ps_trim(fit_padded, method = "pref")

  expect_length(trimmed, nrow(dat))
  expect_true(all(is.na(as.numeric(trimmed)[na_exclude_rows])))
  expect_equal(
    as.numeric(trimmed)[-na_exclude_rows],
    as.numeric(ps_trim(fit_complete, method = "pref"))
  )
})

# ---- the estimator ----------------------------------------------------------

test_that("ipw() accepts the full-length .data a na.exclude fit was given", {
  dat <- sim_na_exclude()
  complete <- complete_rows(dat)

  ps_padded <- glm(
    bin ~ x1 + x2,
    data = dat,
    family = binomial(),
    na.action = na.exclude
  )
  wts_padded <- wt_ate(as.numeric(fitted(ps_padded)), dat$bin)
  out_padded <- glm(
    y ~ bin + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts_padded,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  ps_complete <- glm(bin ~ x1 + x2, data = complete, family = binomial())
  wts_complete <- wt_ate(ps_complete)
  out_complete <- glm(
    y ~ bin + x1,
    data = complete,
    family = quasibinomial(),
    weights = wts_complete,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  padded <- ipw(ps_padded, out_padded, .data = dat)
  reference <- ipw(ps_complete, out_complete, .data = complete)

  expect_equal(padded$estimates, reference$estimates, tolerance = 1e-8)
})

test_that("ipw() still refuses a .data whose complete rows disagree with the fit", {
  # The relaxation is for the rows a fit dropped, not for any longer frame: a
  # `.data` holding more complete rows than the models were fit to is the wrong
  # data and keeps its report.
  dat <- sim_na_exclude()
  complete <- complete_rows(dat)

  ps_mod <- glm(bin ~ x1 + x2, data = complete, family = binomial())
  wts <- wt_ate(ps_mod)
  out_mod <- glm(
    y ~ bin + x1,
    data = complete,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  too_long <- rbind(complete, complete[1:5, , drop = FALSE])

  expect_error(
    ipw(ps_mod, out_mod, .data = too_long),
    class = "propensity_ipw_data_error"
  )
})
