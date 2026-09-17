# ps_refit() refuses truncated scores ------------------------------------------

# Truncation keeps every unit, so refitting on the kept units is refitting the
# original model on every row, and would hand back the scores before they were
# bounded. The refusal says so rather than naming a class.

test_that("ps_refit() refuses a truncated score vector and says why", {
  set.seed(5)
  n <- 80
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(1.5 * x))
  fit <- glm(z ~ x, family = binomial())
  truncated <- ps_trunc(fitted(fit), lower = 0.2, upper = 0.8)

  expect_error(
    ps_refit(truncated, model = fit),
    class = "propensity_method_error"
  )
  expect_propensity_error(ps_refit(truncated, model = fit))
})

test_that("ps_refit() refuses a truncated score matrix and says why", {
  set.seed(6)
  n <- 90
  x <- rnorm(n)
  z <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  fit <- nnet::multinom(z ~ x, trace = FALSE)
  truncated <- ps_trunc(
    predict(fit, type = "probs"),
    method = "ps",
    lower = 0.3,
    .exposure = z
  )

  expect_error(
    ps_refit(truncated, model = fit),
    class = "propensity_method_error"
  )
  expect_propensity_error(ps_refit(truncated, model = fit))
})

test_that("ps_refit() keeps the generic refusal for other classes", {
  fit <- glm(c(0, 1, 0, 1) ~ c(1, 2, 3, 4), family = binomial())

  cnd <- expect_error(
    ps_refit(c(0.2, 0.4, 0.6, 0.8), model = fit),
    class = "propensity_class_error"
  )
  expect_false(inherits(cnd, "propensity_method_error"))
  expect_match(conditionMessage(cnd), "ps_trim")
})
