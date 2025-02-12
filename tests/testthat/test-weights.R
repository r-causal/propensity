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
  n <- 10
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
  n <- 10
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
  n <- 12
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
