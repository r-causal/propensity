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
  expect_error(
    vec_cast(x_with_decimals, integer()),
    class = "vctrs_error_cast_lossy"
  )
})

test_that("vec_ptype2 combines psw and other types correctly", {
  x <- psw(c(0.1, 0.2), estimand = "ate")
  y <- psw(c(0.3, 0.4), estimand = "ate")
  expect_equal(vec_ptype2(x, y), new_psw(estimand = "ate"))

  # Different estimands should fail
  z <- psw(c(0.5, 0.6), estimand = "att")
  expect_error(
    vec_ptype2(x, z),
    "Can't combine weights with different estimands"
  )

  # Combining with double
  expect_equal(vec_ptype2(x, double()), double())
  expect_equal(vec_ptype2(double(), x), double())

  # Combining with integer
  expect_equal(vec_ptype2(x, integer()), integer())
  expect_equal(vec_ptype2(integer(), x), integer())
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
  expect_error(x * thing, class = "vctrs_error_incompatible_op")
  expect_error(x * list(runif(10)), class = "vctrs_error_incompatible_op")
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
  expect_silent(
    type <- ggplot2::scale_type(psw(1))
  )

  expect_identical(type, "continuous")
})
