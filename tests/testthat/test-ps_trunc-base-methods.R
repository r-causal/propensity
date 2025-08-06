# Test base R methods for ps_trunc class
library(testthat)
library(vctrs)

test_that("ps_trunc subsetting with [ preserves class and updates indices", {
  ps <- runif(10, 0.05, 0.95)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  meta <- ps_trunc_meta(x)

  # Single element - not truncated
  non_trunc_idx <- setdiff(1:10, meta$truncated_idx)[1]
  sub1 <- x[non_trunc_idx]
  expect_s3_class(sub1, "ps_trunc")
  expect_equal(length(sub1), 1)

  # Single element - truncated
  if (length(meta$truncated_idx) > 0) {
    trunc_idx <- meta$truncated_idx[1]
    sub2 <- x[trunc_idx]
    expect_s3_class(sub2, "ps_trunc")
    expect_equal(length(sub2), 1)
    # Should be at boundary
    expect_true(
      as.numeric(sub2) == meta$lower_bound ||
        as.numeric(sub2) == meta$upper_bound
    )
  }

  # Multiple elements
  sub3 <- x[1:5]
  expect_s3_class(sub3, "ps_trunc")
  expect_equal(length(sub3), 5)
  sub3_meta <- ps_trunc_meta(sub3)
  expect_true(all(sub3_meta$truncated_idx <= 5))

  # Empty subsetting
  sub4 <- x[integer(0)]
  expect_s3_class(sub4, "ps_trunc")
  expect_equal(length(sub4), 0)
})

test_that("ps_trunc sort() preserves class", {
  ps <- c(0.05, 0.5, 0.95, 0.3)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.5, 0.8, 0.3

  sorted <- sort(x)
  expect_s3_class(sorted, "ps_trunc")
  expect_equal(as.numeric(sorted), c(0.2, 0.3, 0.5, 0.8))

  # Decreasing
  sorted_dec <- sort(x, decreasing = TRUE)
  expect_s3_class(sorted_dec, "ps_trunc")
  expect_equal(as.numeric(sorted_dec), c(0.8, 0.5, 0.3, 0.2))
})

test_that("Understanding current sort() behavior on ps_trunc with metadata", {
  ps <- c(0.05, 0.5, 0.95, 0.3, 0.1)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.5, 0.8, 0.3, 0.2

  # Check original metadata
  meta_original <- ps_trunc_meta(x)
  expect_equal(meta_original$truncated_idx, c(1, 3, 5)) # positions that were truncated
  expect_equal(meta_original$lower_bound, 0.2)
  expect_equal(meta_original$upper_bound, 0.8)

  # Values at truncated positions should be at bounds
  expect_true(all(x[meta_original$truncated_idx] %in% c(0.2, 0.8)))

  # Sort
  sorted <- sort(x)
  expect_s3_class(sorted, "ps_trunc")
  expect_equal(as.numeric(sorted), c(0.2, 0.2, 0.3, 0.5, 0.8))

  # Check metadata after sort
  meta_sorted <- ps_trunc_meta(sorted)
  # Our new sort method should have updated the indices correctly
  expect_equal(meta_sorted$truncated_idx, c(1, 2, 5)) # Correct positions after sorting

  # Verify the indices are correct
  # After sorting, positions 1 and 2 have value 0.2 (at lower bound)
  # Position 5 has value 0.8 (at upper bound)
  actual_truncated <- which(
    as.numeric(sorted) == meta_sorted$lower_bound |
      as.numeric(sorted) == meta_sorted$upper_bound
  )
  expect_equal(actual_truncated, c(1, 2, 5))
  expect_true(identical(actual_truncated, meta_sorted$truncated_idx))
})

test_that("ps_trunc unique() preserves class", {
  ps <- c(0.05, 0.05, 0.5, 0.95, 0.95, 0.5)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.2, 0.5, 0.8, 0.8, 0.5

  uniq <- unique(x)
  expect_s3_class(uniq, "ps_trunc")
  # Should have 0.2, 0.5, 0.8
  expect_equal(sort(as.numeric(uniq)), c(0.2, 0.5, 0.8))
})

test_that("ps_trunc duplicated() returns logical vector", {
  ps <- c(0.05, 0.05, 0.5, 0.95, 0.95, 0.5)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)

  dups <- duplicated(x)
  expect_type(dups, "logical")
  expect_equal(length(dups), 6)
  expect_equal(dups, c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE))

  # Check anyDuplicated
  expect_equal(anyDuplicated(x), 2)
})

test_that("ps_trunc rev() preserves class", {
  ps <- c(0.1, 0.3, 0.5, 0.9)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.3, 0.5, 0.8

  reversed <- rev(x)
  expect_s3_class(reversed, "ps_trunc")
  expect_equal(as.numeric(reversed), c(0.8, 0.5, 0.3, 0.2))
})

test_that("ps_trunc head() and tail() preserve class", {
  ps <- runif(20, 0.05, 0.95)
  x <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)

  # head
  h <- head(x, 5)
  expect_s3_class(h, "ps_trunc")
  expect_equal(length(h), 5)
  h_meta <- ps_trunc_meta(h)
  expect_true(all(h_meta$truncated_idx <= 5))

  # tail
  t <- tail(x, 5)
  expect_s3_class(t, "ps_trunc")
  expect_equal(length(t), 5)
})

test_that("ps_trunc rep() preserves class", {
  ps <- c(0.1, 0.5, 0.9)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.5, 0.8

  # times argument
  r1 <- rep(x, times = 2)
  expect_s3_class(r1, "ps_trunc")
  expect_equal(length(r1), 6)
  expect_equal(as.numeric(r1), c(0.2, 0.5, 0.8, 0.2, 0.5, 0.8))

  # each argument
  r2 <- rep(x, each = 2)
  expect_s3_class(r2, "ps_trunc")
  expect_equal(length(r2), 6)
  expect_equal(as.numeric(r2), c(0.2, 0.2, 0.5, 0.5, 0.8, 0.8))
})

test_that("ps_trunc is.na() and anyNA() work correctly", {
  ps <- c(0.1, 0.3, 0.5, 0.9)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)

  na_vec <- is.na(x)
  expect_type(na_vec, "logical")
  expect_equal(na_vec, c(FALSE, FALSE, FALSE, FALSE))

  expect_false(anyNA(x))
})

test_that("ps_trunc na.omit() preserves class", {
  # Create ps_trunc without NA values in original (ps_trunc doesn't create NAs)
  ps <- c(0.1, 0.5, 0.9)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.5, 0.8

  x_clean <- na.omit(x)
  expect_s3_class(x_clean, "ps_trunc")
  expect_false(anyNA(x_clean))
  expect_equal(length(x_clean), 3)
})

test_that("ps_trunc rejects infinite values", {
  # ps_trunc should reject Inf values
  ps <- c(0.1, 0.3, Inf, 0.5)
  expect_propensity_error(
    ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  )

  # But is.finite and is.infinite should work on valid ps_trunc objects
  ps_valid <- c(0.1, 0.3, 0.5, 0.9)
  x <- ps_trunc(ps_valid, method = "ps", lower = 0.2, upper = 0.8)

  finite_vec <- is.finite(x)
  expect_type(finite_vec, "logical")
  expect_equal(finite_vec, c(TRUE, TRUE, TRUE, TRUE))

  infinite_vec <- is.infinite(x)
  expect_type(infinite_vec, "logical")
  expect_equal(infinite_vec, c(FALSE, FALSE, FALSE, FALSE))
})

test_that("ps_trunc which() family functions work correctly", {
  ps <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.3, 0.5, 0.7, 0.8

  # which with logical condition
  idx <- which(as.numeric(x) > 0.4)
  expect_type(idx, "integer")
  expect_equal(idx, c(3, 4, 5))

  # which.min
  expect_equal(which.min(x), 1)

  # which.max
  expect_equal(which.max(x), 5)
})

test_that("ps_trunc order() and rank() work correctly", {
  ps <- c(0.9, 0.1, 0.7, 0.3)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.8, 0.2, 0.7, 0.3

  # order
  ord <- order(x)
  expect_type(ord, "integer")
  expect_equal(ord, c(2, 4, 3, 1))

  # rank
  rnk <- rank(x)
  expect_type(rnk, "double")
  expect_equal(rnk, c(4, 1, 3, 2))
})

test_that("ps_trunc match() and %in% work correctly", {
  ps1 <- c(0.1, 0.5, 0.9)
  x <- ps_trunc(ps1, method = "ps", lower = 0.2, upper = 0.8)
  # x is: 0.2, 0.5, 0.8

  ps2 <- c(0.5, 0.8)
  y <- ps_trunc(ps2, method = "ps", lower = 0.2, upper = 0.8)
  # y is: 0.5, 0.8

  # match
  m <- match(x, y)
  expect_type(m, "integer")
  expect_equal(m, c(NA, 1, 2))

  # %in%
  inn <- x %in% y
  expect_type(inn, "logical")
  expect_equal(inn, c(FALSE, TRUE, TRUE))
})

test_that("ps_trunc table() works correctly", {
  ps <- c(0.05, 0.05, 0.5, 0.5, 0.95, 0.95)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.2, 0.5, 0.5, 0.8, 0.8

  tbl <- table(x)
  expect_s3_class(tbl, "table")
  expect_equal(as.numeric(tbl), c(2, 2, 2))
  expect_equal(names(tbl), c("0.2", "0.5", "0.8"))
})

test_that("ps_trunc sample() preserves class", {
  set.seed(123)
  ps <- runif(10, 0.05, 0.95)
  x <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)

  # Sample without replacement
  s1 <- sample(x, 5)
  expect_s3_class(s1, "ps_trunc")
  expect_equal(length(s1), 5)

  # Sample with replacement
  s2 <- sample(x, 20, replace = TRUE)
  expect_s3_class(s2, "ps_trunc")
  expect_equal(length(s2), 20)
})

test_that("ps_trunc summary statistics methods work correctly", {
  ps <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.3, 0.5, 0.7, 0.8

  expect_equal(min(x), 0.2)
  expect_equal(max(x), 0.8)
  expect_equal(range(x), c(0.2, 0.8))
  expect_equal(median(x), 0.5)
  expect_equal(quantile(x, 0.5), c("50%" = 0.5))

  # These don't have specific methods but should work
  expect_equal(mean(x), 0.5)
  expect_equal(sum(x), 2.5)
})

test_that("ps_trunc diff() works correctly", {
  ps <- c(0.1, 0.3, 0.5, 0.9)
  x <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.2, 0.3, 0.5, 0.8

  d <- diff(x)
  expect_type(d, "double")
  expect_equal(length(d), 3)
  expect_equal(d, c(0.1, 0.2, 0.3))
})

test_that("ps_trunc works in data frames", {
  ps <- runif(5, 0.05, 0.95)
  df <- data.frame(
    id = 1:5,
    ps_scores = ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)
  )

  expect_s3_class(df$ps_scores, "ps_trunc")

  # Subsetting preserves class
  expect_s3_class(df$ps_scores[1:3], "ps_trunc")
  expect_s3_class(df[1:3, "ps_scores"], "ps_trunc")
})

# Additional edge cases ----
test_that("ps_trunc preserves relative ordering", {
  original <- c(0.01, 0.15, 0.5, 0.85, 0.99)
  truncated <- ps_trunc(original, method = "ps", lower = 0.2, upper = 0.8)

  # Check values
  expect_equal(as.numeric(truncated), c(0.2, 0.2, 0.5, 0.8, 0.8))

  # Check that relative ordering is preserved for non-truncated values
  expect_true(all(diff(truncated) >= 0))

  # More complex case
  original2 <- c(0.1, 0.3, 0.2, 0.7, 0.9, 0.6)
  truncated2 <- ps_trunc(original2, method = "ps", lower = 0.25, upper = 0.75)

  # Non-truncated values should maintain their relative order
  non_trunc_orig <- original2[c(2, 4, 6)] # 0.3, 0.7, 0.6
  non_trunc_result <- truncated2[c(2, 4, 6)]
  expect_equal(order(non_trunc_orig), order(non_trunc_result))
})

test_that("ps_trunc handles all values at same boundary", {
  # All below lower bound
  all_low <- ps_trunc(rep(0.05, 5), method = "ps", lower = 0.2, upper = 0.8)
  expect_true(all(as.numeric(all_low) == 0.2))
  expect_equal(ps_trunc_meta(all_low)$truncated_idx, 1:5)

  # All above upper bound
  all_high <- ps_trunc(rep(0.95, 5), method = "ps", lower = 0.2, upper = 0.8)
  expect_true(all(as.numeric(all_high) == 0.8))
  expect_equal(ps_trunc_meta(all_high)$truncated_idx, 1:5)

  # All exactly at lower bound
  at_lower <- ps_trunc(rep(0.2, 5), method = "ps", lower = 0.2, upper = 0.8)
  expect_true(all(as.numeric(at_lower) == 0.2))
  expect_equal(ps_trunc_meta(at_lower)$truncated_idx, integer(0)) # Not truncated, already at bound

  # All exactly at upper bound
  at_upper <- ps_trunc(rep(0.8, 5), method = "ps", lower = 0.2, upper = 0.8)
  expect_true(all(as.numeric(at_upper) == 0.8))
  expect_equal(ps_trunc_meta(at_upper)$truncated_idx, integer(0)) # Not truncated, already at bound
})
