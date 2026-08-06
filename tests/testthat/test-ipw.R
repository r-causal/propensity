get_variance <- function(est, effect) {
  est$estimates$std.err[est$estimates$effect == effect]^2
}

test_that("variance works", {
  # these objects calculated from Kostouraki et al. See ?ipw
  l0_ATEW_cor <- readRDS(testthat::test_path("data", "l0_ATEW_cor.rds"))
  l0_ATTW_cor <- readRDS(testthat::test_path("data", "l0_ATTW_cor.rds"))
  l0_MW_cor <- readRDS(testthat::test_path("data", "l0_MW_cor.rds"))
  l0_OW_cor <- readRDS(testthat::test_path("data", "l0_OW_cor.rds"))
  l1_ATEW_cor <- readRDS(testthat::test_path("data", "l1_ATEW_cor.rds"))
  l1_ATTW_cor <- readRDS(testthat::test_path("data", "l1_ATTW_cor.rds"))
  l1_MW_cor <- readRDS(testthat::test_path("data", "l1_MW_cor.rds"))
  l1_OW_cor <- readRDS(testthat::test_path("data", "l1_OW_cor.rds"))
  .df <- readRDS(testthat::test_path("data", "df.rds"))

  ps_mod <- glm(
    Z ~ X1 + X2 + X3 + X4 + X5 + X6,
    family = binomial(),
    data = .df
  )

  ps <- predict(ps_mod, type = "response")

  outcome_mod_ate <- glm(
    Y ~ Z,
    weights = wt_ate(ps, .df$Z, exposure_type = "binary", .focal_level = 1),
    data = .df,
    family = quasibinomial()
  )
  outcome_mod_att <- glm(
    Y ~ Z,
    weights = wt_att(ps, .df$Z, exposure_type = "binary", .focal_level = 1),
    data = .df,
    family = quasibinomial()
  )
  outcome_mod_ato <- glm(
    Y ~ Z,
    weights = wt_ato(ps, .df$Z, exposure_type = "binary", .focal_level = 1),
    data = .df,
    family = quasibinomial()
  )
  outcome_mod_atm <- glm(
    Y ~ Z,
    weights = wt_atm(ps, .df$Z, exposure_type = "binary", .focal_level = 1),
    data = .df,
    family = quasibinomial()
  )

  est <- ipw(ps_mod, outcome_mod_ate, se_method = "linearization")

  expect_equal(
    get_variance(est, "rd"),
    var(l1_ATEW_cor - l0_ATEW_cor) / nrow(.df)
  )

  est <- ipw(ps_mod, outcome_mod_att, se_method = "linearization")

  expect_equal(
    get_variance(est, "rd"),
    var(l1_ATTW_cor - l0_ATTW_cor) / nrow(.df)
  )

  est <- ipw(ps_mod, outcome_mod_ato, se_method = "linearization")

  expect_equal(
    get_variance(est, "rd"),
    var(l1_OW_cor - l0_OW_cor) / nrow(.df)
  )

  est <- ipw(ps_mod, outcome_mod_atm, se_method = "linearization")

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
  z <- rbinom(n, 1, plogis(0.2 * x1 - 0.8 * x2))

  # Binary outcome depends on exposure + confounders:
  y <- rbinom(n, 1, plogis(-1 + 1.5 * z + 0.4 * x1 - 0.5 * x2))

  dat <- data.frame(x1, x2, z, y)

  # 1) Fit the PS model
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())

  # 2) Calculate ATE weights
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)

  # 3) Weighted outcome model
  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  # 4) ipw call
  res <- ipw(
    wt_mod = ps_mod,
    outcome_mod = outcome_mod,
    .data = dat,
    estimand = "ate",
    se_method = "linearization"
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
  expect_named(
    est_df,
    c(
      "effect",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )

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
  z <- rbinom(n, 1, plogis(0.2 * x1 + 0.3 * x2))

  # Continuous outcome depends on exposure and confounders
  y <- 5 + 2 * z + 1 * x1 - 0.5 * x2 + rnorm(n)

  dat <- data.frame(x1, x2, z, y)

  # Propensity score model
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())

  # ATE weights
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)

  # Weighted outcome model (linear)
  outcome_mod1 <- lm(y ~ z, data = dat, weights = wts)
  outcome_mod2 <- glm(y ~ z, data = dat, weights = wts)

  # ipw call
  res <- ipw(
    ps_mod,
    outcome_mod1,
    .data = dat,
    estimand = "ate",
    se_method = "linearization"
  )

  expect_snapshot(res)

  # Should only have "diff" for continuous outcomes
  est_df <- res$estimates
  expect_s3_class(res, "ipw")
  expect_equal(unique(est_df$effect), "diff")
  expect_equal(nrow(est_df), 1)

  # Check columns
  expect_named(
    est_df,
    c(
      "effect",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )

  expect_no_error(ipw(ps_mod, outcome_mod2, .data = dat, estimand = "ate"))
})

test_that("wt_mod must be glm, outcome_mod must be glm or lm", {
  set.seed(103)
  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.5)

  # valid ps_mod
  ps_mod <- glm(z ~ x, family = binomial())

  # invalid ps_mod: a class ipw() has no method for
  bad_mod <- structure(list(), class = "not_a_model")

  # valid outcome mod (logistic)
  wts <- rep(1, n)
  outcome_mod <- glm(y ~ z, family = binomial(), weights = wts)

  expect_propensity_error(
    ipw(wt_mod = bad_mod, outcome_mod = outcome_mod)
  )

  expect_propensity_error(
    assert_class("a", "character", .length = 2)
  )

  expect_propensity_error(
    assert_class("a", c("numeric", "character"), .length = 2)
  )

  # invalid outcome_mod
  bad_outcome <- list(call = quote(foo()), class = "list")

  expect_propensity_error(
    ipw(wt_mod = ps_mod, outcome_mod = bad_outcome)
  )

  expect_propensity_error(
    ipw(wt_mod = ps_mod, outcome_mod = outcome_mod, .data = data.frame(x))
  )
})

test_that("ipw handles .data = NULL properly", {
  set.seed(104)
  n <- 200
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 * x))
  y <- rbinom(n, 1, plogis(-1 + 1.5 * z + 0.5 * x))

  data_no_df <- data.frame(x, z, y)

  ps_mod <- glm(z ~ x, data = data_no_df, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)

  outcome_mod <- glm(
    y ~ z,
    data = data_no_df,
    family = quasibinomial(),
    weights = wts
  )

  # .data = NULL => ipw should extract from model frames
  res <- ipw(ps_mod, outcome_mod, .data = NULL)
  expect_s3_class(res, "ipw")
})


test_that("ipw handles various errors correctly", {
  set.seed(104)
  n <- 200
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 * x))
  y <- rbinom(n, 1, plogis(-1 + 1.5 * z + 0.5 * x))

  df <- data.frame(x, z, y)

  ps_mod <- glm(z ~ x, data = df, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)

  outcome_mod_no_estimand <- glm(
    y ~ z,
    data = df,
    family = quasibinomial(),
    weights = as.double(wts)
  )

  expect_propensity_error(
    ipw(ps_mod, outcome_mod_no_estimand)
  )
})

test_that("exponentiate=TRUE in as.data.frame.ipw transforms log(rr), log(or)", {
  set.seed(105)
  n <- 500
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  y <- rbinom(n, 1, plogis(-0.5 + 1.2 * z + 0.3 * x))

  dat <- data.frame(x, z, y)

  ps_mod <- glm(z ~ x, data = dat, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)

  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  ipw_res <- ipw(ps_mod, outcome_mod, .data = dat)

  df_log <- as.data.frame(ipw_res, exponentiate = FALSE)
  df_exp <- as.data.frame(ipw_res, exponentiate = TRUE)

  # The log scale has "log(rr)", "log(or)"
  expect_true(any(df_log$term == "log(rr)"))
  expect_true(any(df_log$term == "log(or)"))

  # Exponentiated scale => "rr", "or"
  expect_true(any(df_exp$term == "rr"))
  expect_true(any(df_exp$term == "or"))

  # "rd" should remain the same in both
  rd_log <- df_log[df_log$term == "rd", "estimate"]
  rd_exp <- df_exp[df_exp$term == "rd", "estimate"]
  expect_equal(rd_log, rd_exp)
})

test_that("Estimand mismatch triggers an error if outcome weights differ from user-specified", {
  # For example, outcome_mod has ATE weights but user tries to specify 'att'
  set.seed(106)
  n <- 300
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 * x))
  y <- rbinom(n, 1, plogis(1 + 0.8 * z + 0.4 * x))

  dat <- data.frame(x, z, y)

  ps_mod <- glm(z ~ x, data = dat, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts_ate <- wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)

  # Weighted outcome model with ATE
  outcome_mod_ate <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts_ate
  )

  # If your code properly checks mismatch, this should raise an error
  expect_propensity_error(
    ipw(
      wt_mod = ps_mod,
      outcome_mod = outcome_mod_ate,
      .data = dat,
      estimand = "att"
    )
  )
})

test_that("ipw works for probit link in the propensity score model", {
  set.seed(2002)
  n <- 400
  x2 <- rnorm(n)
  z <- rbinom(n, 1, pnorm(0.4 * x2))
  y <- rbinom(n, 1, plogis(-0.5 + 1.2 * z + 0.3 * x2))

  dat <- data.frame(x2, z, y)

  # Propensity score model with probit link
  ps_mod <- glm(z ~ x2, data = dat, family = binomial("probit"))
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)

  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  # ipw call
  res <- ipw(ps_mod, outcome_mod, .data = dat, estimand = "ate")

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
  z <- rbinom(n, 1, plogis(0.5 * x3))
  y <- rbinom(n, 1, plogis(-1 + 1.5 * z + 0.8 * x3))

  dat <- data.frame(x3, z, y)

  # Fit the propensity score model with cloglog link
  ps_mod <- glm(z ~ x3, data = dat, family = binomial("cloglog"))
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)

  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  # ipw
  res <- ipw(ps_mod, outcome_mod, .data = dat, estimand = "ate")

  # `ipw` checks
  expect_s3_class(res, "ipw")
  expect_equal(res$estimand, "ate")
  expect_true("estimates" %in% names(res))
  expect_true(any(res$estimates$effect == "rd"))
  expect_true(any(res$estimates$effect == "log(rr)"))
  expect_true(any(res$estimates$effect == "log(or)"))
})

test_that("ipw errors on a length mismatch and a transformed exposure formula", {
  set.seed(3003)
  n <- 400
  x3 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x3))
  y <- rbinom(n, 1, plogis(-1 + 1.5 * z + 0.8 * x3))

  dat <- data.frame(x3, z, y)

  ps_mod <- glm(z ~ x3, data = dat, family = binomial())
  ps <- predict(ps_mod, type = "response")
  wts <- wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)

  outcome_mod_wrong <- glm(
    y ~ z,
    data = dat[1:100, ],
    family = quasibinomial(),
    weights = wts[1:100]
  )

  expect_propensity_error(
    ipw(ps_mod, outcome_mod_wrong, estimand = "ate")
  )

  outcome_mod_transformed <- glm(
    y ~ I(z^2),
    family = quasibinomial(),
    weights = wts
  )
  expect_propensity_error(
    ipw(ps_mod, outcome_mod_transformed, estimand = "ate")
  )
})

test_that("the cannot-determine-estimand error is attributed to ipw()", {
  # Plain numeric weights carry no estimand attribute, and no estimand argument
  # is supplied, so check_estimand cannot determine the estimand. The error must
  # be attributed to ipw(), not the internal check_estimand() helper. The
  # snapshot pins the "Error in `ipw()`" condition header. (The class of this
  # error is already covered by test-ipw-mestimation.R,
  # "ipw_spec_binary errors when the estimand cannot be determined".)
  set.seed(2024)
  n <- 300
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x1))
  y <- rbinom(n, 1, plogis(-0.4 + z))
  dat <- data.frame(x1, z, y)
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- as.double(wt_ate(
    predict(ps_mod, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  ))
  outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)

  # The stored call must be the ipw() frame, not the internal check_estimand()
  # helper. Without the call threaded through, the condition reports
  # check_estimand(...) and this assertion fails.
  cnd <- rlang::catch_cnd(ipw(ps_mod, outcome_mod, .data = dat))
  expect_identical(as.character(cnd$call[[1]]), "ipw")

  expect_snapshot(error = TRUE, ipw(ps_mod, outcome_mod, .data = dat))
})

# ---- arguments that fall into the dots --------------------------------------

# `ipw()` places `...` ahead of its named arguments, so anything the caller does
# not name lands in the dots: `.data` supplied positionally, which earlier
# releases accepted, and any misspelled argument name. Both must error rather
# than run to completion with the supplied value discarded and the default used.

dots_binary_fixture <- function() {
  withr::local_seed(2024)
  n <- 300
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x1))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * x1))
  dat <- data.frame(x1, z, y)
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- wt_ate(
    predict(ps_mod, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  outcome_mod <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  list(dat = dat, ps_mod = ps_mod, outcome_mod = outcome_mod)
}

test_that("ipw() glm rejects .data supplied positionally", {
  fx <- dots_binary_fixture()

  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, fx$dat),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("ipw() glm rejects a misspelled argument", {
  fx <- dots_binary_fixture()

  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, nonsense = 42),
    class = "rlib_error_dots_nonempty"
  )

  # A near miss of a real argument name is the case most likely to go unnoticed:
  # the call reads as though it selected linearization.
  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, se_methd = "linearization"),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("ipw() default rejects an argument it has no name for", {
  fx <- dots_binary_fixture()
  bad_mod <- structure(list(), class = "not_a_model")

  expect_error(
    ipw(bad_mod, fx$outcome_mod, nonsense = 42),
    class = "rlib_error_dots_nonempty"
  )

  # With the dots empty the unsupported class is still what the user hears
  # about.
  expect_error(
    ipw(bad_mod, fx$outcome_mod),
    class = "propensity_method_error"
  )
})

test_that("ipw() glm accepts every argument supplied by name", {
  fx <- dots_binary_fixture()

  # `ps_link` is deprecated but still accepted by name, which is all this pin
  # is about. The deprecation itself is pinned in test-ipw-se-method.R, and
  # silencing it here keeps this test independent of the order the suite runs
  # in, since the default verbosity warns only on the id's first use.
  withr::local_options(lifecycle_verbosity = "quiet")

  baseline <- ipw(fx$ps_mod, fx$outcome_mod)
  named <- ipw(
    fx$ps_mod,
    fx$outcome_mod,
    .data = fx$dat,
    estimand = "ate",
    ps_link = "logit",
    conf_level = 0.95,
    se_method = "mestimation"
  )

  expect_equal(named$estimates, baseline$estimates)
  expect_equal(named$estimates$effect, c("rd", "log(rr)", "log(or)"))
  expect_equal(
    named$estimates$estimate,
    c(0.293149648386106, 0.579445266649124, 1.21028245127107),
    tolerance = 1e-6
  )
  expect_equal(
    named$estimates$std.err,
    c(0.0545082497807004, 0.121118456441513, 0.238716341998527),
    tolerance = 1e-6
  )
})

test_that("ipw() glm accepts every argument supplied by name under linearization", {
  fx <- dots_binary_fixture()

  # As above: the deprecation of `ps_link` is pinned elsewhere, and this pin is
  # only about the argument still being accepted by name.
  withr::local_options(lifecycle_verbosity = "quiet")

  named <- ipw(
    fx$ps_mod,
    fx$outcome_mod,
    .data = fx$dat,
    estimand = "ate",
    ps_link = "logit",
    conf_level = 0.9,
    se_method = "linearization"
  )

  expect_equal(named$se_method, "linearization")
  expect_equal(unique(named$estimates$conf.level), 0.9)
  expect_equal(
    named$estimates$estimate,
    c(0.293149648386107, 0.579445266649124, 1.21028245127107),
    tolerance = 1e-6
  )
  expect_equal(
    named$estimates$std.err,
    c(0.054599327536842, 0.121320826160904, 0.239115226595756),
    tolerance = 1e-6
  )
})

# ---- the estimand has to be one ipw() knows ---------------------------------
#
# check_estimand() compared the argument against the estimand the weights
# recorded and never against the set of estimands that exist, so a typo was
# reported as whatever the comparison happened to notice. With psw weights it
# came back as a mismatch, which sent the user to reconcile two estimands when
# only one of them was real; with plain numeric weights nothing compared it at
# all and it reached the weighted means, where base R reported that `x` and `w`
# had different lengths.
#
# The estimand the weights themselves record is the other source and was checked
# against nothing at all. `psw()` is exported and records the estimand it is
# given, so weights built by hand carry whatever string they were built with
# straight past the argument check and into the same base R report.

estimand_fixture <- function() {
  set.seed(2024)
  n <- 300
  x1 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x1))
  y <- rbinom(n, 1, plogis(-0.4 + z))
  dat <- data.frame(x1, z, y)
  ps_mod <- glm(z ~ x1, data = dat, family = binomial())
  wts <- wt_ate(
    predict(ps_mod, type = "response"),
    dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  dat$plain_wts <- as.numeric(wts)
  outcome_psw <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = wts
  )
  outcome_plain <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = plain_wts
  )
  # weights whose recorded estimand is not one at all, which `psw()` accepts
  banana_wts <- psw(as.numeric(wts), estimand = "banana")
  outcome_banana <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = banana_wts
  )
  list(
    dat = dat,
    ps_mod = ps_mod,
    outcome_psw = outcome_psw,
    outcome_plain = outcome_plain,
    outcome_banana = outcome_banana
  )
}

test_that("an unknown estimand is rejected when the weights record one", {
  skip_if_not_installed("deli")
  fx <- estimand_fixture()

  err <- expect_error(
    ipw(fx$ps_mod, fx$outcome_psw, .data = fx$dat, estimand = "banana"),
    class = "propensity_ipw_estimand_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "banana", fixed = TRUE)
  expect_match(msg, "\"ate\"", fixed = TRUE)
  expect_match(msg, "\"entropy\"", fixed = TRUE)
  # not the mismatch redirect, which described a disagreement between two
  # estimands when only one of the two exists
  expect_false(grepl("Estimand in weights different", msg, fixed = TRUE))
})

test_that("an unknown estimand is rejected when the weights record none", {
  skip_if_not_installed("deli")
  fx <- estimand_fixture()

  for (se in c("mestimation", "linearization")) {
    err <- expect_error(
      ipw(
        fx$ps_mod,
        fx$outcome_plain,
        .data = fx$dat,
        estimand = "banana",
        se_method = se
      ),
      class = "propensity_ipw_estimand_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "banana", fixed = TRUE)
    expect_false(grepl("must have the same length", msg, fixed = TRUE))
  }
})

test_that("weights recording an unknown estimand are rejected on both paths", {
  skip_if_not_installed("deli")
  fx <- estimand_fixture()

  for (se in c("mestimation", "linearization")) {
    err <- expect_error(
      ipw(fx$ps_mod, fx$outcome_banana, .data = fx$dat, se_method = se),
      class = "propensity_ipw_estimand_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "banana", fixed = TRUE)
    expect_match(msg, "\"entropy\"", fixed = TRUE)
    # neither of the two reports this route used to reach: base R's, from the
    # weighted means the unmatched tilt returned nothing to, and the redirect
    # to the other standard error method, which named the value as an estimand
    # that method supports
    expect_false(grepl("must have the same length", msg, fixed = TRUE))
    expect_false(grepl("mestimation", msg, fixed = TRUE))
  }
})

test_that("the unknown-estimand report names the weights as the source", {
  skip_if_not_installed("deli")
  fx <- estimand_fixture()

  expect_snapshot(
    error = TRUE,
    ipw(fx$ps_mod, fx$outcome_banana, .data = fx$dat)
  )
})

# Membership is a set of names, and `%in%` reads its left side through
# `as.character()`, so a value of another type matched the name it prints as and
# passed. `list("ate")`, which is what single-bracket indexing of an options list
# returns, then reached the marginal means as a list and died in base R's terms;
# `factor("att")` reached the tilt, which is a `switch()`, where a factor is its
# integer level code, so a branch was selected by position, base R noted the
# coercion four times, and the fit either completed under an estimand nobody
# asked for or was reported as a weights mismatch.

test_that("an estimand argument of the wrong type is rejected on both paths", {
  skip_if_not_installed("deli")
  fx <- estimand_fixture()
  wrong_type <- list(list("ate"), factor("att"), factor("ate"))

  for (se in c("mestimation", "linearization")) {
    for (estimand in wrong_type) {
      err <- expect_error(
        ipw(
          fx$ps_mod,
          fx$outcome_plain,
          .data = fx$dat,
          estimand = estimand,
          se_method = se
        ),
        class = "propensity_ipw_estimand_error"
      )
      msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
      expect_match(msg, "single string", fixed = TRUE)
      expect_false(grepl("must have the same length", msg, fixed = TRUE))
    }
  }
})

test_that("an estimand argument of the wrong type leaks no base conditions", {
  skip_if_not_installed("deli")
  fx <- estimand_fixture()

  # `factor("ate")` spells an estimand that exists, so it passed membership and
  # ran to completion, reporting nothing but base R's note about the coercion
  for (se in c("mestimation", "linearization")) {
    seen <- character()
    withCallingHandlers(
      expect_error(
        ipw(
          fx$ps_mod,
          fx$outcome_plain,
          .data = fx$dat,
          estimand = factor("ate"),
          se_method = se
        ),
        class = "propensity_ipw_estimand_error"
      ),
      warning = function(w) {
        seen <<- c(seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_identical(seen, character())
  }
})

test_that("weights recording an estimand of the wrong type are rejected", {
  skip_if_not_installed("deli")
  fx <- estimand_fixture()
  # `psw()` rejects this at construction, so the attribute is set past it: the
  # weights-side check is what stands between a psw built some other way and the
  # tilt
  wts <- psw(as.numeric(fx$dat$plain_wts), estimand = "ate")
  attr(wts, "estimand") <- factor("att")
  outcome_mod <- glm(
    y ~ z,
    data = fx$dat,
    family = quasibinomial(),
    weights = wts
  )

  err <- expect_error(
    ipw(fx$ps_mod, outcome_mod, .data = fx$dat),
    class = "propensity_ipw_estimand_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "Their estimand has class", fixed = TRUE)
  expect_match(msg, "single string", fixed = TRUE)
})

# A small continuous fixture, local to this file.
sim_continuous_estimand <- function() {
  set.seed(2024)
  n <- 300
  x1 <- rnorm(n)
  A <- 0.5 + 0.8 * x1 + rnorm(n)
  yc <- 1 + 0.6 * A + 0.5 * x1 + rnorm(n)
  data.frame(x1, A, yc)
}

test_that("an unknown estimand is rejected on the continuous path", {
  skip_if_not_installed("deli")
  dat <- sim_continuous_estimand()
  ps_mod <- lm(A ~ x1, data = dat)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(fitted(ps_mod)),
      dat$A,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )
  outcome_mod <- lm(yc ~ A, data = dat, weights = wts)

  expect_error(
    ipw(ps_mod, outcome_mod, estimand = "banana"),
    class = "propensity_ipw_estimand_error"
  )
})

test_that("the continuous path reports an unknown estimand as unknown", {
  skip_if_not_installed("deli")
  dat <- sim_continuous_estimand()
  ps_mod <- lm(A ~ x1, data = dat)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(fitted(ps_mod)),
      dat$A,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )
  outcome_mod <- lm(yc ~ A, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, outcome_mod, estimand = "banana"),
    class = "propensity_ipw_estimand_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "must name an estimand", fixed = TRUE)
  # the ate-only restriction read as though the value were an estimand that
  # some other exposure type supports
  expect_false(grepl("not available for a continuous", msg, fixed = TRUE))

  # a real estimand this exposure type does not support still reads that way
  att <- expect_error(
    ipw(ps_mod, outcome_mod, estimand = "att"),
    class = "propensity_ipw_estimand_error"
  )
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(att)),
    "not available for a continuous",
    fixed = TRUE
  )
})

test_that("weights recording an unknown estimand are rejected on the continuous path", {
  skip_if_not_installed("deli")
  dat <- sim_continuous_estimand()
  ps_mod <- lm(A ~ x1, data = dat)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(fitted(ps_mod)),
      dat$A,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )
  banana_wts <- psw(as.numeric(wts), estimand = "banana")
  outcome_mod <- lm(yc ~ A, data = dat, weights = banana_wts)

  err <- expect_error(
    ipw(ps_mod, outcome_mod),
    class = "propensity_ipw_estimand_error"
  )
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(err)),
    "banana",
    fixed = TRUE
  )
})

test_that("a known estimand that disagrees with the weights still redirects", {
  skip_if_not_installed("deli")
  fx <- estimand_fixture()

  err <- expect_error(
    ipw(fx$ps_mod, fx$outcome_psw, .data = fx$dat, estimand = "att"),
    class = "propensity_ipw_estimand_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "Estimand in weights different", fixed = TRUE)
  expect_match(msg, "\"ate\" vs. \"att\"", fixed = TRUE)
})

test_that("the undeterminable estimand error carries the estimand class", {
  skip_if_not_installed("deli")
  fx <- estimand_fixture()

  expect_error(
    ipw(fx$ps_mod, fx$outcome_plain, .data = fx$dat),
    class = "propensity_ipw_estimand_error"
  )
})

test_that("every estimand ipw() supports passes the membership check", {
  for (estimand in ipw_estimands) {
    expect_identical(check_estimand(NULL, estimand), estimand)
  }
})

# ---- the exposure guards carry the exposure class ---------------------------
#
# Two guards report an exposure that is not binary, one per standard error
# path, and neither carried an error class: both arrived as the bare
# propensity_error every guard in the package shares, so a caller could not tell
# them from an offset, a link, or a response.

three_value_exposure_fixture <- function() {
  set.seed(9)
  n <- 300
  x1 <- rnorm(n)
  # a proportion response, which a binomial fit accepts and which takes three
  # distinct values rather than two
  p3 <- c(0, 0.5, 1)[cut(plogis(0.4 * x1), 3, labels = FALSE)]
  y <- rbinom(n, 1, 0.4)
  dat <- data.frame(x1, p3, y)
  ps_mod <- suppressWarnings(glm(p3 ~ x1, data = dat, family = binomial()))
  wts <- psw(rep(1, n), estimand = "ate")
  outcome_mod <- glm(
    y ~ p3,
    data = dat,
    family = quasibinomial(),
    weights = wts
  )
  list(dat = dat, ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

test_that("the mestimation binary-exposure guard carries the exposure class", {
  fx <- three_value_exposure_fixture()

  # the fixture is the shape it claims to be
  expect_length(
    unique(stats::model.response(stats::model.frame(fx$ps_mod))),
    3L
  )

  err <- expect_error(
    ipw_spec_binary(fx$ps_mod, fx$outcome_mod, .data = fx$dat),
    class = "propensity_ipw_exposure_error"
  )
  expect_match(conditionMessage(err), "binary exposures", fixed = TRUE)
})

test_that("the linearization binary-exposure guard carries the exposure class", {
  fx <- three_value_exposure_fixture()

  err <- expect_error(
    estimate_marginal_means(
      outcome_mod = fx$outcome_mod,
      wts = fx$wts,
      exposure = fx$dat$p3,
      exposure_name = "p3",
      .data = fx$dat
    ),
    class = "propensity_ipw_exposure_error"
  )
  expect_match(conditionMessage(err), "binary exposures", fixed = TRUE)
})

test_that("the exposure and outcome length guard carries the length class", {
  set.seed(3003)
  n <- 400
  x3 <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x3))
  y <- rbinom(n, 1, plogis(-1 + 1.5 * z + 0.8 * x3))
  dat <- data.frame(x3, z, y)
  ps_mod <- glm(z ~ x3, data = dat, family = binomial())
  wts <- wt_ate(
    predict(ps_mod, type = "response"),
    z,
    exposure_type = "binary",
    .focal_level = 1
  )
  outcome_mod <- glm(
    y ~ z,
    data = dat[1:100, ],
    family = quasibinomial(),
    weights = wts[1:100]
  )

  expect_error(
    ipw(ps_mod, outcome_mod, estimand = "ate", se_method = "linearization"),
    class = "propensity_length_error"
  )
})
