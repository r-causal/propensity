test_that("ATE works for binary cases", {
  expect_message(
    weights <- wt_ate(c(.1, .3, .4, .3), .exposure = c(0, 0, 1, 0)),
    "Treating `.exposure` as binary"
  )

  expect_silent(
    wt_ate(
      c(.1, .3, .4, .3),
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    )
  )

  expect_equal(
    weights,
    psw(c(1.11, 1.43, 2.50, 1.43), "ate"),
    tolerance = .01
  )
})

test_that("ATE works for continuous cases", {
  denom_model <- lm(mpg ~ gear + am + carb, data = mtcars)
  num <- dnorm(mtcars$mpg, mean(mtcars$mpg), sd(mtcars$mpg))
  denom <- dnorm(mtcars$mpg, predict(denom_model), mean(influence(denom_model)$sigma))
  wts <- 1 / denom
  stb_wts <- num / denom

  expect_message(
    weights <- wt_ate(
      predict(denom_model),
      .exposure = mtcars$mpg,
      .sigma = influence(denom_model)$sigma
    ),
    "Treating `.exposure` as continuous",
  )

  stablized_weights <- wt_ate(
    predict(denom_model),
    .exposure = mtcars$mpg,
    .sigma = influence(denom_model)$sigma,
    exposure_type = "continuous",
    stabilize = TRUE
  )

  expect_equal(weights, psw(wts, "ate"), tolerance = .01)
  expect_equal(stablized_weights, psw(stb_wts, "ate"), tolerance = .01)
})
test_that("wt_ate() with ps_trim issues refit warning if not refit, no warning if refit", {
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 1) Trim
  trimmed_ps <- ps_trim(ps, exposure = z, method = "ps", lower = 0.2, upper = 0.8)

  # not refit => expect a warning
  expect_warning(
    w_ate_unfit <- wt_ate(trimmed_ps, .exposure = z, exposure_type = "binary", .treated = 1),
    "did not refit"
  )
  expect_s3_class(w_ate_unfit, "psw")
  expect_true(grepl("; trimmed$", estimand(w_ate_unfit)))

  # 2) After refit => no warning
  trimmed_refit <- ps_refit(trimmed_ps, model = fit)
  expect_silent(
    w_ate_fit <- wt_ate(trimmed_refit, .exposure = z, exposure_type = "binary", .treated = 1)
  )
  expect_s3_class(w_ate_fit, "psw")
  expect_true(grepl("; trimmed$", estimand(w_ate_fit)))
})

test_that("wt_ate() with ps_trunc adds '; truncated' without refit warning", {
  set.seed(234)
  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.6 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # e.g. bounding at [0.2, 0.8]
  truncated_ps <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)

  # Should produce weighting with no refit warnings
  expect_silent(
    w_ate_trunc <- wt_ate(truncated_ps, .exposure = z, exposure_type = "binary", .treated = 1)
  )
  expect_s3_class(w_ate_trunc, "psw")
  # Estimand ends with "; truncated"
  expect_true(grepl("; truncated$", estimand(w_ate_trunc)))
})

test_that("Other estimands (att, atu, etc.) with ps_trim or ps_trunc", {
  set.seed(345)
  n <- 120
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 + 0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # Trim
  trimmed_ps <- ps_trim(ps, exposure = z, method = "ps")
  # No refit => warning
  expect_warning(
    w_att_trim <- wt_att(trimmed_ps, .exposure = z, exposure_type = "binary", .treated = 1),
    "did not refit"
  )
  # Check estimand
  expect_true(grepl("att; trimmed", estimand(w_att_trim)))

  # Trunc
  truncated_ps <- ps_trunc(ps, method = "pctl", lower = 0.2, upper = 0.8)
  # No warning
  expect_silent(
    w_att_trunc <- wt_att(truncated_ps, .exposure = z, exposure_type = "binary", .treated = 1)
  )
  expect_silent(
    w_atu_trunc <- wt_atu(truncated_ps, .exposure = z, exposure_type = "binary", .treated = 1)
  )
  expect_true(grepl("atu; truncated", estimand(w_atu_trunc)))
})

test_that("wt_ate() with ps_trunc sets truncated=TRUE in final psw", {
  set.seed(123)
  n <- 8
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  trunc_obj <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  w_ate <- wt_ate(trunc_obj, .exposure = z, exposure_type = "binary", .treated = 1)

  expect_true(is_truncated(w_ate))
  expect_false(is_trimmed(w_ate))
  expect_match(estimand(w_ate), "; truncated$")
})

test_that("wt_atu.ps_trim triggers refit check, sets 'atu; trimmed'", {
  set.seed(991)
  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(1.6 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 1) Trim the PS
  trimmed_obj <- ps_trim(ps, exposure = z, method = "ps", lower = 0.2, upper = 0.8)

  # Not refit => we get a warning
  expect_warning(
    w_atu_unfit <- wt_atu(trimmed_obj, .exposure = z, exposure_type = "binary", .treated = 1),
    "did not refit"
  )
  expect_s3_class(w_atu_unfit, "psw")
  expect_match(estimand(w_atu_unfit), "atu; trimmed")
  expect_true(attr(w_atu_unfit, "trimmed"))
  # ps_trim_meta copied
  expect_identical(attr(w_atu_unfit, "ps_trim_meta"), attr(trimmed_obj, "ps_trim_meta"))

  # 2) Now refit => no warning
  refit_obj <- ps_refit(trimmed_obj, model = fit)
  expect_silent(
    w_atu_fit <- wt_atu(refit_obj, .exposure = z, exposure_type = "binary", .treated = 1)
  )
  expect_s3_class(w_atu_fit, "psw")
  expect_match(estimand(w_atu_fit), "atu; trimmed")
  expect_true(attr(w_atu_fit, "trimmed"))
  # confirm ps_trim_meta matches
  expect_identical(attr(w_atu_fit, "ps_trim_meta"), attr(refit_obj, "ps_trim_meta"))
})

test_that("wt_atm.ps_trim triggers refit check, sets 'atm; trimmed'", {
  set.seed(992)
  n <- 50
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  trimmed_obj <- ps_trim(ps, exposure = z, method = "ps", lower = 0.2, upper = 0.8)

  # Not refit => warning
  expect_warning(
    w_atm_unfit <- wt_atm(trimmed_obj, .exposure = z, exposure_type = "binary", .treated = 1),
    "did not refit"
  )
  expect_s3_class(w_atm_unfit, "psw")
  expect_match(estimand(w_atm_unfit), "atm; trimmed")
  expect_true(attr(w_atm_unfit, "trimmed"))

  # Refit => no warning
  refit_obj <- ps_refit(trimmed_obj, model = fit)
  expect_silent(
    w_atm_fit <- wt_atm(refit_obj, .exposure = z, exposure_type = "binary", .treated = 1)
  )
  expect_s3_class(w_atm_fit, "psw")
  expect_match(estimand(w_atm_fit), "atm; trimmed")
  expect_true(attr(w_atm_fit, "trimmed"))
})

test_that("wt_ato.ps_trim triggers refit check, sets 'ato; trimmed'", {
  set.seed(993)
  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  trimmed_obj <- ps_trim(ps, exposure = z, method = "ps", lower = 0.1, upper = 0.9)

  # Not refit => warning
  expect_warning(
    w_ato_unfit <- wt_ato(trimmed_obj, .exposure = z, exposure_type = "binary", .treated = 1),
    "did not refit"
  )
  expect_s3_class(w_ato_unfit, "psw")
  expect_match(estimand(w_ato_unfit), "ato; trimmed")
  expect_true(attr(w_ato_unfit, "trimmed"))

  # Refit => no warning
  refit_obj <- ps_refit(trimmed_obj, model = fit)
  expect_silent(
    w_ato_fit <- wt_ato(refit_obj, .exposure = z, exposure_type = "binary", .treated = 1)
  )
  expect_s3_class(w_ato_fit, "psw")
  expect_match(estimand(w_ato_fit), "ato; trimmed")
  expect_true(attr(w_ato_fit, "trimmed"))
})
