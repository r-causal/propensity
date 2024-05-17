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
    c(1.11, 1.43, 2.50, 1.43),
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

  expect_equal(weights, wts, tolerance = .01)
  expect_equal(stablized_weights, stb_wts, tolerance = .01)
})
