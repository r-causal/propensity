test_that("deli integration: M-estimation solves and returns sandwich variance", {
  skip_if_not_installed("deli")

  set.seed(1)
  y <- rnorm(50, mean = 3, sd = 2)

  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- deli::MEstimator(stacked_equations = psi, init = c(mu = 0))
  m <- deli::estimate(m)

  # the M-estimator for a mean solves to the sample mean
  expect_equal(unname(coef(m)), mean(y), tolerance = 1e-8)
  expect_named(coef(m), "mu")

  # empirical sandwich variance for a mean equals the finite-sample
  # variance of y divided by n
  se_sandwich <- sqrt(diag(vcov(m)))
  se_closed_form <- sqrt(mean((y - mean(y))^2) / length(y))
  expect_equal(unname(se_sandwich), se_closed_form, tolerance = 1e-6)
})
