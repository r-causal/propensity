library(testthat)
library(vctrs)

test_that("new_psw creates valid psw objects", {
  x <- new_psw(c(0.1, 0.2, 0.3), estimand = "ate")
  expect_s3_class(
    x,
    c("psw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vec_data(x), c(0.1, 0.2, 0.3))
  expect_equal(estimand(x), "ate")
})

test_that("psw helper function works correctly", {
  x <- psw(c(0.1, 0.2, 0.3), estimand = "att")
  expect_s3_class(
    x,
    c("psw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vec_data(x), c(0.1, 0.2, 0.3))
  expect_equal(estimand(x), "att")
  expect_false(is_stabilized(x))
  expect_false(is_ps_trimmed(x))
  expect_false(is_ps_truncated(x))
  estimand(x) <- "ATT!"
  expect_equal(estimand(x), "ATT!")
})

test_that("is_psw identifies psw objects", {
  x <- psw(c(0.1, 0.2), estimand = "ato")
  expect_true(is_psw(x))
  expect_false(is_psw(c(0.1, 0.2)))
})

test_that("vec_ptype_abbr and vec_ptype_full return correct type labels", {
  x <- psw(c(0.1, 0.2), estimand = "atm")
  expect_equal(vec_ptype_abbr(x), "psw{atm}")
  expect_equal(vec_ptype_full(x), "psw{estimand = atm}")

  z <- psw(c(0.1, 0.2))
  expect_equal(vec_ptype_abbr(z), "psw")
  expect_equal(vec_ptype_full(z), "psw{estimand = unknown}")

  y <- psw(c(0.1, 0.2), estimand = "cens")
  x <- x * y
  expect_equal(vec_ptype_abbr(x), "psw{atm, cens}")
  expect_equal(vec_ptype_full(x), "psw{estimand = atm, cens}")
})

test_that("vec_cast works for psw and double", {
  x <- psw(c(0.1, 0.2, 0.3), estimand = "ate")

  # Cast to double
  double_x <- vec_cast(x, double())
  expect_equal(double_x, c(0.1, 0.2, 0.3))
  expect_true(is.numeric(double_x))

  # Cast back to psw
  psw_x <- as_psw(double_x, estimand = "ate")
  expect_s3_class(
    psw_x,
    c("psw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vec_data(psw_x), double_x)
  expect_equal(estimand(psw_x), "ate")
})

test_that("vec_cast works for psw and integer with precision checks", {
  x <- psw(c(1, 2, 3), estimand = "ate")

  # Cast to integer
  int_x <- vec_cast(x, integer())
  expect_equal(int_x, as.integer(c(1, 2, 3)))
  expect_true(is.integer(int_x))

  # Cast back to psw
  psw_x <- as_psw(int_x, estimand = "ate")
  expect_s3_class(
    psw_x,
    c("psw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vec_data(psw_x), as.numeric(int_x))
  expect_equal(estimand(psw_x), "ate")

  # Fail when precision is lost
  x_with_decimals <- psw(c(1.1, 2.2, 3.3), estimand = "ate")
  expect_propensity_error(
    vec_cast(x_with_decimals, integer())
  )
})

# A prototype carrying a distinct value in each of the six metadata fields, so a
# cast that drops any one of them is visible.
psw_cast_prototype <- function(stabilization_score = 0.42) {
  psw(
    double(),
    estimand = "ate",
    stabilized = TRUE,
    trimmed = TRUE,
    truncated = FALSE,
    calibrated = TRUE,
    stabilization_score = stabilization_score
  )
}

test_that("vec_cast from double carries every metadata field of the psw prototype", {
  to <- psw_cast_prototype()

  out <- vec_cast(c(1, 2, 3), to = to)

  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(1, 2, 3))
  expect_equal(estimand(out), "ate")
  expect_true(is_stabilized(out))
  expect_true(is_ps_trimmed(out))
  expect_false(is_ps_truncated(out))
  expect_true(is_ps_calibrated(out))
  expect_equal(stabilization_score(out), 0.42)
})

test_that("vec_cast from integer carries every metadata field of the psw prototype", {
  to <- psw_cast_prototype()

  out <- vec_cast(c(1L, 2L, 3L), to = to)

  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(1, 2, 3))
  expect_equal(estimand(out), "ate")
  expect_true(is_stabilized(out))
  expect_true(is_ps_trimmed(out))
  expect_false(is_ps_truncated(out))
  expect_true(is_ps_calibrated(out))
  expect_equal(stabilization_score(out), 0.42)
})

test_that("vec_cast from ps_trim carries every metadata field of the psw prototype", {
  to <- psw_cast_prototype()
  ps <- ps_trim(c(0.2, 0.5, 0.8), method = "ps", lower = 0.1, upper = 0.9)

  out <- vec_cast(ps, to = to)

  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(0.2, 0.5, 0.8))
  expect_equal(estimand(out), "ate")
  expect_true(is_stabilized(out))
  expect_true(is_ps_trimmed(out))
  expect_false(is_ps_truncated(out))
  expect_true(is_ps_calibrated(out))
  expect_equal(stabilization_score(out), 0.42)
})

test_that("vec_cast from ps_trunc carries every metadata field of the psw prototype", {
  to <- psw_cast_prototype()
  ps <- ps_trunc(c(0.2, 0.5, 0.8), method = "ps", lower = 0.1, upper = 0.9)

  out <- vec_cast(ps, to = to)

  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(0.2, 0.5, 0.8))
  expect_equal(estimand(out), "ate")
  expect_true(is_stabilized(out))
  expect_true(is_ps_trimmed(out))
  expect_false(is_ps_truncated(out))
  expect_true(is_ps_calibrated(out))
  expect_equal(stabilization_score(out), 0.42)
})

test_that("vec_cast keeps a per-observation stabilization score at a matching length", {
  score <- c(0.51, 0.52, 0.53)
  to <- psw_cast_prototype(stabilization_score = score)

  out <- expect_no_warning(vec_cast(c(1, 2, 3), to = to))

  expect_s3_class(out, "psw")
  expect_equal(stabilization_score(out), score)
  expect_true(is_stabilized(out))
})

test_that("vec_cast drops a per-observation stabilization score at a different length silently", {
  score <- c(0.51, 0.52, 0.53)
  to <- psw_cast_prototype(stabilization_score = score)

  # The score describes the observations behind `to`, not the data being cast,
  # so a length it does not match is nothing the caller can act on.
  out <- expect_silent(vec_cast(c(1, 2), to = to))

  expect_s3_class(out, "psw")
  expect_length(out, 2)
  expect_null(stabilization_score(out))
  expect_null(attr(out, "stabilization_score"))
  expect_true(is_stabilized(out))
  expect_equal(estimand(out), "ate")
})

test_that("vec_cast from psw to double returns the bare data", {
  x <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    trimmed = TRUE,
    calibrated = TRUE,
    stabilization_score = 0.42
  )

  out <- vec_cast(x, to = double())

  expect_identical(out, c(1, 2, 3))
  expect_null(attributes(out))
})

test_that("vec_ptype2 combines psw and other types correctly", {
  x <- psw(c(0.1, 0.2), estimand = "ate")
  y <- psw(c(0.3, 0.4), estimand = "ate")
  expect_equal(vec_ptype2(x, y), new_psw(estimand = "ate"))

  # Different estimands should warn and return numeric
  z <- psw(c(0.5, 0.6), estimand = "att")
  expect_propensity_warning(
    # jarl-ignore implicit_assignment: assignment keeps the return value out of the snapshot while the warning is captured
    result <- vec_ptype2(x, z)
  )
  expect_identical(result, double())

  # Combining with double
  expect_propensity_warning(
    expect_equal(vec_ptype2(x, double()), double())
  )
  expect_propensity_warning(
    expect_equal(vec_ptype2(double(), x), double())
  )

  # Combining with integer
  expect_propensity_warning(
    expect_equal(vec_ptype2(x, integer()), integer())
  )
  expect_propensity_warning(
    expect_equal(vec_ptype2(integer(), x), integer())
  )
})

test_that("vec_arith performs arithmetic on psw objects", {
  x <- psw(c(1, 2, 3), estimand = "ate")
  y <- psw(c(0.5, 1.5, 2.5), estimand = "ate")
  cens <- psw(c(3, 2, 1), estimand = "cens")

  classes <- c("psw", "causal_wts", "vctrs_vctr", "double")

  # same estimand
  result <- x + y
  expect_s3_class(result, classes, exact = TRUE)
  expect_equal(vec_data(result), c(1.5, 3.5, 5.5))
  expect_equal(estimand(result), "ate")

  # different estimand
  result <- x * cens
  expect_s3_class(result, classes, exact = TRUE)
  expect_equal(vec_data(result), c(3, 4, 3))
  expect_equal(estimand(result), "ate, cens")

  # Arithmetic with double
  result <- x * 2
  expect_s3_class(result, classes, exact = TRUE)
  expect_equal(vec_data(result), c(2, 4, 6))
  expect_equal(estimand(result), "ate")

  # Arithmetic with integer
  result <- x - 1L
  expect_s3_class(result, classes, exact = TRUE)
  expect_equal(vec_data(result), c(0, 1, 2))
  expect_equal(estimand(result), "ate")

  # various combos work in various orders
  expect_no_error(x * 2.1)
  expect_no_error(2.1 * x)
  expect_no_error(x / 2.1)
  expect_no_error(2.1 / x)
  expect_no_error(x + 2.1)
  expect_no_error(2.1 + x)
  expect_no_error(x * 2L)
  expect_no_error(2L * x)
  expect_no_error(x / 2L)
  expect_no_error(2L / x)
  expect_no_error(x + 2L)
  expect_no_error(2L + x)
  expect_no_error(x^2)

  # doesn't work with unsupported types
  thing <- new_vctr(runif(10), class = "thing")
  expect_propensity_error(x * thing)
  expect_propensity_error(x * list(runif(10)))
})

test_that("vec_math applies math functions to psw objects", {
  x <- psw(c(1, 4, 9), estimand = "ate")

  expect_false(is_psw(sqrt(x)))
  expect_false(is_psw(sum(x)))
})

test_that("Combination of arithmetic and math works correctly for psw", {
  # Example data
  y <- c(2, 4, 6, 8)
  x <- c(0, 1, 0, 1)
  wts <- psw(c(0.5, 1.0, 1.5, 2.0), estimand = "ate")

  # Compute expected results using double weights
  wts_double <- vec_data(wts) # Strip the psw class to use as raw weights

  expected_term1 <- sum(y * x * wts_double) / sum(x * wts_double)
  expected_term2 <- sum(y * (1 - x) * wts_double) / sum((1 - x) * wts_double)

  # Compute the actual results using psw weights
  term1 <- sum(y * x * wts) / sum(x * wts)
  term2 <- sum(y * (1 - x) * wts) / sum((1 - x) * wts)

  # Validate the results
  expect_equal(term1, expected_term1, tolerance = 1e-8)
  expect_equal(term2, expected_term2, tolerance = 1e-8)

  # Check for consistency in repeated calculations
  repeated_term2 <- sum(y * (1 - x) * wts) / sum((1 - x) * wts)
  expect_equal(repeated_term2, term2, tolerance = 1e-8)
})

test_that("Refit logic can be tracked if stored in attr", {
  # If you store a 'refit' attribute in the final psw object:
  w <- psw(c(0.1, 0.2), estimand = "ate", trimmed = TRUE)
  attr(w, "ps_trim_meta") <- list(refit = TRUE)
  # Suppose we define is_refit.psw() as checking attr(x, "refit")
  expect_true(is_refit(w))
})

test_that("psw objects can convert to character", {
  x <- as.character(new_psw(c(0.1, 0.2, 0.3), estimand = "ate"))
  expect_type(x, "character")
})

test_that("psw works with ggplot2", {
  skip_if_not_installed("ggplot2")
  type <- expect_silent(
    ggplot2::scale_type(psw(1))
  )

  expect_identical(type, "continuous")
})
