test_that("errors for non-numeric ps", {
  expect_error(
    ps_calibrate("not numeric", c(0, 1)),
    "`ps` must be a numeric vector"
  )
})

test_that("errors for out-of-range ps", {
  expect_error(
    ps_calibrate(c(-0.1, 0.2), c(0, 1)),
    "`ps` values must be between 0 and 1"
  )
  expect_error(
    ps_calibrate(c(0.5, 1.1), c(0, 1)),
    "`ps` values must be between 0 and 1"
  )
})

test_that("errors when ps and treat have different lengths", {
  expect_error(
    ps_calibrate(runif(5), rep(0:1, length.out = 6)),
    "same length"
  )
})

test_that("returns a psw object of correct length and range", {
  set.seed(42)
  ps <- rep(0.5, 100)
  treat <- rbinom(100, 1, 0.3)

  out <- ps_calibrate(ps, treat)

  expect_s3_class(out, "psw")
  expect_length(out, 100)
  expect_true(all(out >= 0 & out <= 1))
})

test_that("constant ps yields calibrated = observed prevalence", {
  ps <- rep(0.5, 20)
  treat <- rep(c(0, 1), each = 10) # prevalence = 0.5

  out <- ps_calibrate(ps, treat)
  # all values should equal the 0.5 prevalence
  expect_equal(unique(as.numeric(out)), 0.5)
})

test_that("preserves psw attributes from an existing causal‐weights object", {
  # Create a dummy psw object
  ps_orig <- psw(
    x          = runif(10),
    estimand   = "ATT",
    stabilized = TRUE,
    trimmed    = TRUE,
    truncated  = FALSE
  )
  treat <- rbinom(10, 1, ps_orig)

  out <- ps_calibrate(ps_orig, treat)

  expect_equal(attr(out, "estimand"), "ATT")
  expect_true(attr(out, "stabilized"))
  expect_true(attr(out, "trimmed"))
  expect_false(attr(out, "truncated"))
})
