# Test Base R Methods for ALL vctrs Classes
#
# This file tests common base R methods to discover which ones need implementation
# for psw, ps_trim, and ps_trunc classes.

library(testthat)

# =============================================================================
# Helper Functions
# =============================================================================

create_psw_example <- function(n = 10) {
  psw(runif(n, 0.5, 2.5), estimand = "ate", stabilized = TRUE)
}

create_ps_trim_example <- function(n = 10) {
  ps <- runif(n, 0.05, 0.95)
  ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
}

create_ps_trunc_example <- function(n = 10) {
  ps <- runif(n, 0.05, 0.95)
  ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)
}

# =============================================================================
# quantile() Tests
# =============================================================================

test_that("quantile works for psw", {
  x <- create_psw_example()
  expect_type(quantile(x), "double")
  expect_type(quantile(x, probs = c(0.25, 0.75)), "double")
  expect_length(quantile(x), 5) # Default quartiles
})

test_that("quantile works for ps_trim", {
  x <- create_ps_trim_example()
  expect_type(quantile(x, na.rm = TRUE), "double")
  expect_type(quantile(x, probs = c(0.25, 0.75), na.rm = TRUE), "double")
})

test_that("quantile works for ps_trunc", {
  x <- create_ps_trunc_example()
  expect_type(quantile(x), "double")
  expect_type(quantile(x, probs = c(0.25, 0.75)), "double")
})

# =============================================================================
# sort() Tests
# =============================================================================

test_that("sort works for psw", {
  set.seed(123)
  x <- create_psw_example()
  sorted <- sort(x)
  expect_s3_class(sorted, "psw")
  expect_equal(estimand(sorted), estimand(x))
  expect_true(all(diff(as.numeric(sorted)) >= 0))

  # Decreasing sort
  sorted_desc <- sort(x, decreasing = TRUE)
  expect_s3_class(sorted_desc, "psw")
  expect_true(all(diff(as.numeric(sorted_desc)) <= 0))
})

test_that("sort works for ps_trim", {
  set.seed(124)
  x <- create_ps_trim_example()
  sorted <- sort(x, na.last = TRUE)
  expect_s3_class(sorted, "ps_trim")

  # Check non-NA values are sorted
  non_na <- !is.na(sorted)
  if (any(non_na)) {
    sorted_vals <- as.numeric(sorted[non_na])
    expect_true(all(diff(sorted_vals) >= 0))
  }
})

test_that("sort works for ps_trunc", {
  set.seed(125)
  x <- create_ps_trunc_example()
  sorted <- sort(x)
  expect_s3_class(sorted, "ps_trunc")
  expect_true(all(diff(as.numeric(sorted)) >= 0))
})

# =============================================================================
# unique() Tests
# =============================================================================

test_that("unique works for psw", {
  x <- psw(c(1, 2, 2, 3, 3, 3), estimand = "att")
  uniq <- unique(x)
  expect_s3_class(uniq, "psw")
  expect_equal(estimand(uniq), "att")
  expect_equal(length(uniq), 3)
  expect_equal(as.numeric(uniq), c(1, 2, 3))
})

test_that("unique works for ps_trim", {
  # Create ps_trim with known duplicates
  ps <- c(0.3, 0.3, 0.5, 0.5, 0.1, 0.9) # 0.1 and 0.9 will be trimmed
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  uniq <- unique(x)
  expect_s3_class(uniq, "ps_trim")
  # Should have unique values plus NAs
  expect_true(length(uniq) <= length(x))
})

test_that("unique works for ps_trunc", {
  ps <- c(0.05, 0.15, 0.15, 0.85, 0.85, 0.95)
  x <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)
  uniq <- unique(x)
  expect_s3_class(uniq, "ps_trunc")
  expect_true(length(uniq) <= length(x))
})

# =============================================================================
# duplicated() and anyDuplicated() Tests
# =============================================================================

test_that("duplicated/anyDuplicated work for psw", {
  x <- psw(c(1, 2, 2, 3, 3, 3), estimand = "ato")

  dup <- duplicated(x)
  expect_type(dup, "logical")
  expect_length(dup, length(x))
  expect_equal(dup, c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE))

  any_dup <- anyDuplicated(x)
  expect_type(any_dup, "integer")
  expect_gt(any_dup, 0) # Should indicate position of first duplicate
})

test_that("duplicated/anyDuplicated work for ps_trim", {
  ps <- c(0.3, 0.3, 0.5, 0.5, 0.1, 0.9)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)

  dup <- duplicated(x)
  expect_type(dup, "logical")
  expect_length(dup, length(x))

  any_dup <- anyDuplicated(x)
  expect_type(any_dup, "integer")
})

test_that("duplicated/anyDuplicated work for ps_trunc", {
  ps <- c(0.15, 0.15, 0.5, 0.5, 0.85, 0.85)
  x <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)

  dup <- duplicated(x)
  expect_type(dup, "logical")
  expect_length(dup, length(x))

  any_dup <- anyDuplicated(x)
  expect_type(any_dup, "integer")
  expect_gt(any_dup, 0)
})

# =============================================================================
# diff() Tests
# =============================================================================

test_that("diff works for psw", {
  x <- psw(c(1, 3, 6, 10), estimand = "ate")
  d <- diff(x)
  expect_type(d, "double")
  expect_length(d, length(x) - 1)
  expect_equal(d, c(2, 3, 4))

  # lag > 1
  d2 <- diff(x, lag = 2)
  expect_type(d2, "double")
  expect_equal(d2, c(5, 7))
})

test_that("diff works for ps_trim", {
  x <- ps_trim(c(0.2, 0.4, 0.6, 0.8), method = "ps", lower = 0.1, upper = 0.9)
  d <- diff(x)
  expect_type(d, "double")
  expect_length(d, length(x) - 1)
})

test_that("diff works for ps_trunc", {
  x <- ps_trunc(c(0.2, 0.4, 0.6, 0.8), method = "ps", lower = 0.1, upper = 0.9)
  d <- diff(x)
  expect_type(d, "double")
  expect_length(d, length(x) - 1)
  expect_equal(d, rep(0.2, 3))
})

# =============================================================================
# rev() Tests
# =============================================================================

test_that("rev works for psw", {
  x <- psw(1:5, estimand = "att", stabilized = FALSE)
  reversed <- rev(x)
  expect_s3_class(reversed, "psw")
  expect_equal(as.numeric(reversed), 5:1)
  expect_equal(estimand(reversed), "att")
  expect_false(is_stabilized(reversed))
})

test_that("rev works for ps_trim", {
  ps <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  reversed <- rev(x)
  expect_s3_class(reversed, "ps_trim")
  expect_equal(length(reversed), length(x))
})

test_that("rev works for ps_trunc", {
  x <- ps_trunc(
    c(0.1, 0.3, 0.5, 0.7, 0.9),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  reversed <- rev(x)
  expect_s3_class(reversed, "ps_trunc")
  expect_equal(as.numeric(reversed), c(0.8, 0.7, 0.5, 0.3, 0.2))
})

# =============================================================================
# is.finite() and is.infinite() Tests
# =============================================================================

test_that("is.finite/is.infinite work for psw", {
  x <- psw(c(1, 2, Inf, -Inf, NaN), estimand = "ate")

  fin <- is.finite(x)
  expect_type(fin, "logical")
  expect_equal(fin, c(TRUE, TRUE, FALSE, FALSE, FALSE))

  inf <- is.infinite(x)
  expect_type(inf, "logical")
  expect_equal(inf, c(FALSE, FALSE, TRUE, TRUE, FALSE))
})

test_that("is.finite/is.infinite work for ps_trim", {
  # ps_trim shouldn't have infinite values normally, but test the method
  x <- create_ps_trim_example()

  fin <- is.finite(x)
  expect_type(fin, "logical")

  inf <- is.infinite(x)
  expect_type(inf, "logical")
  expect_true(all(!inf[!is.na(x)])) # Non-NA values should not be infinite
})

test_that("is.finite/is.infinite work for ps_trunc", {
  x <- create_ps_trunc_example()

  fin <- is.finite(x)
  expect_type(fin, "logical")
  expect_true(all(fin)) # All truncated values should be finite

  inf <- is.infinite(x)
  expect_type(inf, "logical")
  expect_true(all(!inf)) # No truncated values should be infinite
})

# =============================================================================
# NA-handling Methods Tests
# =============================================================================

test_that("is.na/anyNA work for psw", {
  x <- psw(c(1, 2, NA, 3), estimand = "ate")

  na_check <- is.na(x)
  expect_type(na_check, "logical")
  expect_equal(na_check, c(FALSE, FALSE, TRUE, FALSE))

  any_na <- anyNA(x)
  expect_type(any_na, "logical")
  expect_true(any_na)
})

test_that("is.na/anyNA work for ps_trim", {
  x <- ps_trim(c(0.1, 0.5, 0.9), method = "ps", lower = 0.2, upper = 0.8)

  na_check <- is.na(x)
  expect_type(na_check, "logical")
  expect_equal(na_check, c(TRUE, FALSE, TRUE)) # 0.1 and 0.9 are trimmed

  any_na <- anyNA(x)
  expect_type(any_na, "logical")
  expect_true(any_na)
})

test_that("is.na/anyNA work for ps_trunc", {
  x <- ps_trunc(c(0.1, 0.5, 0.9), method = "ps", lower = 0.2, upper = 0.8)

  na_check <- is.na(x)
  expect_type(na_check, "logical")
  expect_true(all(!na_check)) # Truncation doesn't create NAs

  any_na <- anyNA(x)
  expect_type(any_na, "logical")
  expect_false(any_na)
})

test_that("na.omit works for psw", {
  x <- psw(c(1, 2, NA, 3, NA), estimand = "atm")

  no_na <- na.omit(x)
  expect_s3_class(no_na, "psw")
  expect_false(anyNA(no_na))
  expect_equal(length(no_na), 3)
  expect_equal(estimand(no_na), "atm")
})

test_that("na.omit works for ps_trim", {
  x <- ps_trim(
    c(0.1, 0.3, 0.5, 0.7, 0.9),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  no_na <- na.omit(x)
  expect_s3_class(no_na, "ps_trim")
  expect_false(anyNA(no_na))
  expect_equal(as.numeric(no_na), c(0.3, 0.5, 0.7))
})

test_that("na.omit works for ps_trunc", {
  # Even though truncation doesn't create NAs, na.omit should still work
  x <- ps_trunc(c(0.2, 0.5, 0.8), method = "ps", lower = 0.1, upper = 0.9)

  no_na <- na.omit(x)
  expect_s3_class(no_na, "ps_trunc")
  expect_equal(length(no_na), length(x)) # No NAs removed
})

# =============================================================================
# Additional Useful Methods
# =============================================================================

test_that("which.min/which.max work for psw", {
  x <- psw(c(3, 1, 4, 1, 5), estimand = "ate")

  expect_equal(which.min(x), 2) # First occurrence of minimum
  expect_equal(which.max(x), 5)
})

test_that("head/tail work for psw", {
  x <- psw(1:10, estimand = "att")

  h <- head(x, 3)
  expect_s3_class(h, "psw")
  expect_equal(length(h), 3)
  expect_equal(as.numeric(h), 1:3)

  t <- tail(x, 3)
  expect_s3_class(t, "psw")
  expect_equal(length(t), 3)
  expect_equal(as.numeric(t), 8:10)
})

test_that("head/tail work for ps_trim", {
  x <- create_ps_trim_example(20)

  h <- head(x, 5)
  expect_s3_class(h, "ps_trim")
  expect_equal(length(h), 5)

  t <- tail(x, 5)
  expect_s3_class(t, "ps_trim")
  expect_equal(length(t), 5)
})

test_that("head/tail work for ps_trunc", {
  x <- create_ps_trunc_example(20)

  h <- head(x, 5)
  expect_s3_class(h, "ps_trunc")
  expect_equal(length(h), 5)

  t <- tail(x, 5)
  expect_s3_class(t, "ps_trunc")
  expect_equal(length(t), 5)
})

# =============================================================================
# complete.cases() for data.frame integration
# =============================================================================

test_that("complete.cases works with vctrs classes in data.frames", {
  df <- data.frame(
    id = 1:5,
    weight = psw(c(1, 2, NA, 3, 4), estimand = "ate"),
    ps_trimmed = ps_trim(
      c(0.1, 0.5, 0.6, 0.7, 0.9),
      method = "ps",
      lower = 0.2,
      upper = 0.8
    )
  )

  cc <- complete.cases(df)
  expect_type(cc, "logical")
  expect_equal(cc, c(FALSE, TRUE, FALSE, TRUE, FALSE)) # Rows 1,3,5 have NAs
})
