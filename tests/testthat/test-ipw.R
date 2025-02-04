get_variance <- function(est, effect) {
  est$estimates$std.err[est$estimates$effect == effect]^2
}

test_that("variance works", {
  # these objects calculated from Kostouraki et al. See ?ipw
  l0_ATEW_cor <- readRDS(testthat::test_path("data", "l0_ATEW_cor.rds"))
  l0_ATTW_cor <- readRDS(testthat::test_path("data", "l0_ATTW_cor.rds"))
  l0_MW_cor   <- readRDS(testthat::test_path("data", "l0_MW_cor.rds"))
  l0_OW_cor   <- readRDS(testthat::test_path("data", "l0_OW_cor.rds"))
  l1_ATEW_cor <- readRDS(testthat::test_path("data", "l1_ATEW_cor.rds"))
  l1_ATTW_cor <- readRDS(testthat::test_path("data", "l1_ATTW_cor.rds"))
  l1_MW_cor   <- readRDS(testthat::test_path("data", "l1_MW_cor.rds"))
  l1_OW_cor   <- readRDS(testthat::test_path("data", "l1_OW_cor.rds"))
  .df <- readRDS(testthat::test_path("data", "df.rds"))

  ps_mod <- glm(
    Z ~ X1 + X2 + X3 + X4 + X5 + X6,
    family = binomial(),
    data = .df
  )

  ps <- predict(ps_mod, type = "response")

  outcome_mod_ate <- glm(Y ~ Z, weights = wt_ate(ps, .df$Z, exposure_type = "binary", .treated = 1), data = .df, family = quasibinomial())
  outcome_mod_att <- glm(Y ~ Z, weights = wt_att(ps, .df$Z, exposure_type = "binary", .treated = 1), data = .df, family = quasibinomial())
  outcome_mod_ato <- glm(Y ~ Z, weights = wt_ato(ps, .df$Z, exposure_type = "binary", .treated = 1), data = .df, family = quasibinomial())
  outcome_mod_atm <- glm(Y ~ Z, weights = wt_atm(ps, .df$Z, exposure_type = "binary", .treated = 1), data = .df, family = quasibinomial())

  est <- ipw(ps_mod, outcome_mod_ate)

  expect_equal(
    get_variance(est, "rd"),
    var(l1_ATEW_cor - l0_ATEW_cor) / nrow(.df)
  )

  est <- ipw(ps_mod, outcome_mod_att)

  expect_equal(
    get_variance(est, "rd"),
    var(l1_ATTW_cor - l0_ATTW_cor) / nrow(.df)
  )

  est <- ipw(ps_mod, outcome_mod_ato)

  expect_equal(
    get_variance(est, "rd"),
    var(l1_OW_cor - l0_OW_cor) / nrow(.df)
  )

  est <- ipw(ps_mod, outcome_mod_atm)

  expect_equal(
    get_variance(est, "rd"),
    var(l1_MW_cor - l0_MW_cor) / nrow(.df)
  )
})

test_that("ipw works for binary outcome with a confounder, using logistic ps, logistic outcome", {
  set.seed(101)
  n <- 100

  # Two confounders (continuous and binary):
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.3)

  # Exposure depends on both confounders:
  z <- rbinom(n, 1, plogis(0.2*x1 - 0.8*x2))

  # Binary outcome depends on exposure + confounders:
  y <- rbinom(n, 1, plogis(-1 + 1.5*z + 0.4*x1 - 0.5*x2))

  dat <- data.frame(x1, x2, z, y)

  # 1) Fit the PS model
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())

  # 2) Calculate ATE weights
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .treated = 1)

  # 3) Weighted outcome model
  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  # 4) ipw call
  res <- ipw(
    ps_mod = ps_mod,
    outcome_mod = outcome_mod,
    .df = dat,
    estimand = "ate"
  )

  expect_snapshot(res)

  # `ipw` checks
  expect_s3_class(res, "ipw")
  expect_true(is.list(res))
  expect_true("estimand" %in% names(res))
  expect_true("estimates" %in% names(res))

  # For binary outcomes, we should see 'rd', 'log(rr)', 'log(or)'
  est_df <- res$estimates
  expect_s3_class(est_df, "data.frame")
  expect_named(est_df, c("effect", "estimate", "std.err", "z",
                        "ci.lower", "ci.upper", "conf.level", "p.value"))

  expect_true("rd" %in% est_df$effect)
  expect_true("log(rr)" %in% est_df$effect)
  expect_true("log(or)" %in% est_df$effect)

  # No NAs in the main columns
  expect_false(anyNA(est_df$estimate))
})

test_that("ipw works for continuous outcome with a confounder, using logistic ps, linear outcome", {
  set.seed(102)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)

  # Exposure depends on confounders
  z <- rbinom(n, 1, plogis(0.2*x1 + 0.3*x2))

  # Continuous outcome depends on exposure and confounders
  y <- 5 + 2*z + 1*x1 - 0.5*x2 + rnorm(n)

  dat <- data.frame(x1, x2, z, y)

  # Propensity score model
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())

  # ATE weights
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .treated = 1)

  # Weighted outcome model (linear)
  outcome_mod1 <- lm(y ~ z, data = dat, weights = wts)
  outcome_mod2 <- glm(y ~ z, data = dat, weights = wts)

  # ipw call
  res <- ipw(ps_mod, outcome_mod1, .df = dat, estimand = "ate")

  expect_snapshot(res)

  # Should only have "diff" for continuous outcomes
  est_df <- res$estimates
  expect_s3_class(res, "ipw")
  expect_equal(unique(est_df$effect), "diff")
  expect_equal(nrow(est_df), 1)

  # Check columns
  expect_named(est_df, c("effect", "estimate", "std.err", "z",
                    "ci.lower", "ci.upper", "conf.level", "p.value"))

  expect_no_error(ipw(ps_mod, outcome_mod2, .df = dat, estimand = "ate"))
})

test_that("ps_mod must be glm, outcome_mod must be glm or lm", {
  set.seed(103)
  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.5)

  # valid ps_mod
  ps_mod <- glm(z ~ x, family = binomial())

  # invalid ps_mod
  bad_mod <- lm(z ~ x)

  # valid outcome mod (logistic)
  wts <- rep(1, n)
  outcome_mod <- glm(y ~ z, family = binomial(), weights = wts)

  expect_error(ipw(ps_mod = bad_mod, outcome_mod = outcome_mod),
               regexp = "inherits\\(ps_mod, \"glm\"\\) is not TRUE")

  # invalid outcome_mod
  bad_outcome <- list(call = quote(foo()), class = "list")

  expect_error(ipw(ps_mod = ps_mod, outcome_mod = bad_outcome),
               "inherits\\(outcome_mod, \"glm\"\\) \\|\\| inherits\\(outcome_mod, \"lm\"\\) is not TRUE")
})

test_that("ipw handles .df = NULL properly", {
  set.seed(104)
  n <- 200
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2*x))
  y <- rbinom(n, 1, plogis(-1 + 1.5*z + 0.5*x))

  data_no_df <- data.frame(x, z, y)

  ps_mod <- glm(z ~ x, data = data_no_df, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .treated = 1)

  outcome_mod <- glm(y ~ z, data = data_no_df, family = quasibinomial(), weights = wts)

  # .df = NULL => ipw should extract from model frames
  res <- ipw(ps_mod, outcome_mod, .df = NULL)
  expect_s3_class(res, "ipw")
})


test_that("ipw handles various errors correctly", {
  set.seed(104)
  n <- 200
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2*x))
  y <- rbinom(n, 1, plogis(-1 + 1.5*z + 0.5*x))

  df <- data.frame(x, z, y)

  ps_mod <- glm(z ~ x, data = df, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .treated = 1)

  outcome_mod_no_estimand <- glm(y ~ z, data = df, family = quasibinomial(), weights = as.double(wts))

  expect_error(
    ipw(ps_mod, outcome_mod_no_estimand),
    "Can't determine estimand from weights. Please specify `estimand`."
  )
})

test_that("exponentiate=TRUE in as.data.frame.ipw transforms log(rr), log(or)", {
  set.seed(105)
  n <- 500
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5*x))
  y <- rbinom(n, 1, plogis(-0.5 + 1.2*z + 0.3*x))

  dat <- data.frame(x, z, y)

  ps_mod <- glm(z ~ x, data = dat, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .treated = 1)

  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  ipw_res <- ipw(ps_mod, outcome_mod, .df = dat)

  df_log <- as.data.frame(ipw_res, exponentiate = FALSE)
  df_exp <- as.data.frame(ipw_res, exponentiate = TRUE)

  # The log scale has "log(rr)", "log(or)"
  expect_true(any(df_log$effect == "log(rr)"))
  expect_true(any(df_log$effect == "log(or)"))

  # Exponentiated scale => "rr", "or"
  expect_true(any(df_exp$effect == "rr"))
  expect_true(any(df_exp$effect == "or"))

  # "rd" should remain the same in both
  rd_log <- df_log[df_log$effect == "rd", "estimate"]
  rd_exp <- df_exp[df_exp$effect == "rd", "estimate"]
  expect_equal(rd_log, rd_exp)
})

test_that("Estimand mismatch triggers an error if outcome weights differ from user-specified", {
  # For example, outcome_mod has ATE weights but user tries to specify 'att'
  set.seed(106)
  n <- 300
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2*x))
  y <- rbinom(n, 1, plogis(1 + 0.8*z + 0.4*x))

  dat <- data.frame(x, z, y)

  ps_mod <- glm(z ~ x, data = dat, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts_ate <- wt_ate(ps, z, exposure_type = "binary", .treated = 1)

  # Weighted outcome model with ATE
  outcome_mod_ate <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts_ate)

  # If your code properly checks mismatch, this should raise an error
  expect_error(
    ipw(ps_mod = ps_mod, outcome_mod = outcome_mod_ate, .df = dat, estimand = "att"),
    "Estimand in weights different from `estimand`"
  )
})

test_that("ipw works for probit link in the propensity score model", {
  set.seed(2002)
  n <- 400
  x2 <- rnorm(n)
  z  <- rbinom(n, 1, pnorm(0.4*x2))
  y  <- rbinom(n, 1, plogis(-0.5 + 1.2*z + 0.3*x2))

  dat <- data.frame(x2, z, y)

  # Propensity score model with probit link
  ps_mod <- glm(z ~ x2, data = dat, family = binomial("probit"))
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .treated = 1)

  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  # ipw call
  res <- ipw(ps_mod, outcome_mod, .df = dat, estimand = "ate")

  expect_s3_class(res, "ipw")
  expect_equal(res$estimand, "ate")

  # Quick check: should have 'rd', 'log(rr)', 'log(or)'
  est_df <- res$estimates
  expect_true(all(c("rd", "log(rr)", "log(or)") %in% est_df$effect))
})

test_that("ipw works for cloglog link in the propensity score model", {
  set.seed(3003)
  n <- 400
  x3 <- rnorm(n)
  # Generating exposure from a cloglog perspective is trickier,
  # but we can just stick to logistic generation for simplicity
  # and fit cloglog anyway:
  z <- rbinom(n, 1, plogis(0.5*x3))
  y <- rbinom(n, 1, plogis(-1 + 1.5*z + 0.8*x3))

  dat <- data.frame(x3, z, y)

  # Fit the propensity score model with cloglog link
  ps_mod <- glm(z ~ x3, data = dat, family = binomial("cloglog"))
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .treated = 1)

  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  # ipw
  res <- ipw(ps_mod, outcome_mod, .df = dat, estimand = "ate")

  # `ipw` checks
  expect_s3_class(res, "ipw")
  expect_equal(res$estimand, "ate")
  expect_true("estimates" %in% names(res))
  expect_true(any(res$estimates$effect == "rd"))
  expect_true(any(res$estimates$effect == "log(rr)"))
  expect_true(any(res$estimates$effect == "log(or)"))
})

