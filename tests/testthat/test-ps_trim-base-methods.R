# Test base R methods for ps_trim class
library(testthat)
library(vctrs)

test_that("ps_trim subsetting with [ preserves class and updates indices", {
  ps <- runif(10, 0.05, 0.95)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  meta <- ps_trim_meta(x)

  # Single element from kept
  kept_idx <- meta$keep_idx[1]
  sub1 <- x[kept_idx]
  expect_s3_class(sub1, "ps_trim")
  expect_equal(length(sub1), 1)
  expect_false(is.na(sub1))

  # Single element from trimmed
  trim_idx <- meta$trimmed_idx[1]
  sub2 <- x[trim_idx]
  expect_s3_class(sub2, "ps_trim")
  expect_equal(length(sub2), 1)
  expect_true(is.na(sub2))

  # Multiple elements
  sub3 <- x[1:5]
  expect_s3_class(sub3, "ps_trim")
  expect_equal(length(sub3), 5)
  sub3_meta <- ps_trim_meta(sub3)
  expect_true(all(sub3_meta$trimmed_idx <= 5))
  expect_true(all(sub3_meta$keep_idx <= 5))

  # Logical subsetting - should error with wrong length
  expect_error(
    x[c(TRUE, FALSE)],
    "Logical subscript `i` must be size 1 or 10, not 2",
    class = "propensity_length_error"
  )

  # Logical subsetting - recycling with length 1
  sub4 <- x[TRUE] # All elements
  expect_s3_class(sub4, "ps_trim")
  expect_equal(length(sub4), 10)

  # Logical subsetting - full length
  logical_vec <- rep(c(TRUE, FALSE), 5)
  sub4b <- x[logical_vec]
  expect_s3_class(sub4b, "ps_trim")
  expect_equal(length(sub4b), 5)

  # Empty subsetting
  sub5 <- x[integer(0)]
  expect_s3_class(sub5, "ps_trim")
  expect_equal(length(sub5), 0)

  # Negative indices
  sub6 <- x[-c(1, 2)]
  expect_s3_class(sub6, "ps_trim")
  expect_equal(length(sub6), 8)
})

test_that("ps_trim sort() preserves class but changes indices", {
  ps <- c(0.3, 0.1, 0.9, 0.05, 0.5)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.3, NA, NA, NA, 0.5

  sorted <- sort(x)
  expect_s3_class(sorted, "ps_trim")
  expect_equal(as.numeric(sorted), c(0.3, 0.5)) # Default drops NAs

  # With na.last = FALSE
  sorted_na_first <- sort(x, na.last = FALSE)
  expect_s3_class(sorted_na_first, "ps_trim")
  expect_equal(as.numeric(sorted_na_first), c(NA, NA, NA, 0.3, 0.5))

  # Decreasing
  sorted_dec <- sort(x, decreasing = TRUE)
  expect_s3_class(sorted_dec, "ps_trim")
  expect_equal(as.numeric(sorted_dec), c(0.5, 0.3)) # Default drops NAs
})

test_that("Understanding current sort() behavior on ps_trim with metadata", {
  ps <- c(0.7, 0.1, 0.5, 0.9, 0.3)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.7, NA, 0.5, NA, 0.3

  # Check original metadata
  meta_original <- ps_trim_meta(x)
  expect_equal(meta_original$trimmed_idx, c(2, 4))
  expect_equal(meta_original$keep_idx, c(1, 3, 5))

  # Sort with default (drops NAs)
  sorted_default <- sort(x)
  expect_s3_class(sorted_default, "ps_trim")
  expect_equal(length(sorted_default), 3) # Only non-NA values
  expect_equal(as.numeric(sorted_default), c(0.3, 0.5, 0.7))

  # Check metadata after default sort
  meta_sorted_default <- ps_trim_meta(sorted_default)
  # After dropping NAs, no trimmed indices remain
  expect_equal(meta_sorted_default$trimmed_idx, integer(0)) # No trimmed values after dropping NAs

  # Sort with na.last = TRUE
  sorted_na_last <- sort(x, na.last = TRUE)
  expect_s3_class(sorted_na_last, "ps_trim")
  expect_equal(length(sorted_na_last), 5)
  expect_equal(as.numeric(sorted_na_last), c(0.3, 0.5, 0.7, NA, NA))

  # Check metadata after sort with na.last = TRUE
  meta_na_last <- ps_trim_meta(sorted_na_last)
  # Our new sort method should have updated the indices correctly
  expect_equal(meta_na_last$trimmed_idx, c(4, 5)) # NAs are now at positions 4 and 5

  # Verify the indices are correct
  expect_false(is.na(sorted_na_last[2])) # Position 2 is 0.5, not NA
  expect_true(is.na(sorted_na_last[4])) # Position 4 is NA
  expect_true(is.na(sorted_na_last[5])) # Position 5 is NA

  # Verify NA positions match metadata
  expect_equal(which(is.na(sorted_na_last)), meta_na_last$trimmed_idx)
})

test_that("ps_trim unique() preserves class", {
  ps <- c(0.3, 0.3, 0.1, 0.9, 0.5, 0.5)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.3, 0.3, NA, NA, 0.5, 0.5

  uniq <- unique(x)
  expect_s3_class(uniq, "ps_trim")
  # Should preserve one of each unique value including NA
  expect_true(0.3 %in% uniq)
  expect_true(0.5 %in% uniq)
  expect_true(anyNA(uniq))
  expect_true(length(uniq) <= 3) # At most 0.3, 0.5, and NA
})

test_that("ps_trim duplicated() returns logical vector", {
  ps <- c(0.3, 0.3, 0.1, 0.9, 0.5, 0.5)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)

  dups <- duplicated(x)
  expect_type(dups, "logical")
  expect_equal(length(dups), 6)

  # Check anyDuplicated
  expect_true(anyDuplicated(x) > 0)
})

test_that("ps_trim rev() preserves class and reverses indices", {
  ps <- c(0.1, 0.3, 0.5, 0.9)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: NA, 0.3, 0.5, NA

  reversed <- rev(x)
  expect_s3_class(reversed, "ps_trim")
  expect_equal(as.numeric(reversed), c(NA, 0.5, 0.3, NA))
})

test_that("ps_trim head() and tail() preserve class", {
  ps <- runif(20, 0.05, 0.95)
  x <- ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7)

  # head
  h <- head(x, 5)
  expect_s3_class(h, "ps_trim")
  expect_equal(length(h), 5)
  h_meta <- ps_trim_meta(h)
  expect_true(all(h_meta$trimmed_idx <= 5))
  expect_true(all(h_meta$keep_idx <= 5))

  # tail
  t <- tail(x, 5)
  expect_s3_class(t, "ps_trim")
  expect_equal(length(t), 5)

  # Default n = 6
  expect_equal(length(head(x)), 6)
  expect_equal(length(tail(x)), 6)
})

test_that("ps_trim rep() preserves class", {
  ps <- c(0.2, 0.5, 0.9)
  x <- ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7)
  # Original: NA, 0.5, NA

  # times argument
  r1 <- rep(x, times = 2)
  expect_s3_class(r1, "ps_trim")
  expect_equal(length(r1), 6)
  expect_equal(as.numeric(r1), c(NA, 0.5, NA, NA, 0.5, NA))

  # each argument
  r2 <- rep(x, each = 2)
  expect_s3_class(r2, "ps_trim")
  expect_equal(length(r2), 6)
  expect_equal(as.numeric(r2), c(NA, NA, 0.5, 0.5, NA, NA))
})

test_that("ps_trim is.na() and anyNA() work correctly", {
  ps <- c(0.1, 0.3, 0.5, 0.9)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)

  na_vec <- is.na(x)
  expect_type(na_vec, "logical")
  expect_equal(na_vec, c(TRUE, FALSE, FALSE, TRUE))

  expect_true(anyNA(x))
})

test_that("ps_trim na.omit() preserves class and updates indices", {
  ps <- c(0.1, 0.3, 0.5, 0.9)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: NA, 0.3, 0.5, NA

  x_clean <- na.omit(x)
  expect_s3_class(x_clean, "ps_trim")
  expect_equal(as.numeric(x_clean), c(0.3, 0.5))
  expect_equal(length(x_clean), 2)

  # Check that metadata is updated
  clean_meta <- ps_trim_meta(x_clean)
  expect_equal(length(clean_meta$trimmed_idx), 0)
  expect_equal(length(clean_meta$keep_idx), 2)
})

test_that("ps_trim rejects infinite values", {
  # ps_trim should reject Inf values
  ps <- c(0.1, 0.3, Inf, 0.5)
  expect_error(
    ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8),
    "The propensity score must be between 0 and 1",
    class = "propensity_range_error"
  )

  # But is.finite and is.infinite should work on valid ps_trim objects
  ps_valid <- c(0.1, 0.3, 0.5, 0.9)
  x <- ps_trim(ps_valid, method = "ps", lower = 0.2, upper = 0.8)

  finite_vec <- is.finite(x)
  expect_type(finite_vec, "logical")
  # NA values (trimmed) return FALSE for is.finite
  expect_equal(finite_vec, c(FALSE, TRUE, TRUE, FALSE))

  infinite_vec <- is.infinite(x)
  expect_type(infinite_vec, "logical")
  expect_equal(infinite_vec, c(FALSE, FALSE, FALSE, FALSE))
})

test_that("ps_trim which() family functions work correctly", {
  ps <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: NA, 0.3, 0.5, 0.7, NA

  # which with logical condition (ignores NAs)
  idx <- which(as.numeric(x) > 0.4)
  expect_type(idx, "integer")
  expect_equal(idx, c(3, 4))

  # which.min (ignores NAs)
  expect_equal(which.min(x), 2)

  # which.max (ignores NAs)
  expect_equal(which.max(x), 4)
})

test_that("ps_trim order() and rank() work correctly", {
  ps <- c(0.5, 0.2, 0.7, 0.1)
  x <- ps_trim(ps, method = "ps", lower = 0.3, upper = 0.8)
  # Original: 0.5, NA, 0.7, NA

  # order
  ord <- order(x, na.last = TRUE)
  expect_type(ord, "integer")
  expect_equal(as.numeric(x[ord]), c(0.5, 0.7, NA, NA))

  # rank
  rnk <- rank(x, na.last = TRUE)
  expect_type(rnk, "double")
  # 0.5 is rank 1, 0.7 is rank 2, NAs are ranks 3 and 4
  expect_equal(rnk[!is.na(x)], c(1, 2))
})

test_that("ps_trim match() and %in% work correctly", {
  ps1 <- c(0.3, 0.5, 0.7)
  x <- ps_trim(ps1, method = "ps", lower = 0.4, upper = 0.8)
  # x is: NA, 0.5, 0.7

  ps2 <- c(0.5, 0.9)
  y <- ps_trim(ps2, method = "ps", lower = 0.4, upper = 0.8)
  # y is: 0.5, NA

  # match
  m <- match(x, y)
  expect_type(m, "integer")
  # NA matches NA, 0.5 matches position 1, 0.7 has no match
  expect_equal(m[2], 1) # 0.5 matches position 1 in y

  # %in%
  inn <- x %in% y
  expect_type(inn, "logical")
  # NA is in y, 0.5 is in y, 0.7 is not in y
  expect_true(inn[1]) # NA is in y
  expect_true(inn[2]) # 0.5 is in y
})

test_that("ps_trim table() works correctly", {
  ps <- c(0.3, 0.3, 0.1, 0.5, 0.5, 0.9)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: 0.3, 0.3, NA, 0.5, 0.5, NA

  tbl <- table(x, useNA = "ifany")
  expect_s3_class(tbl, "table")
  # Check the actual table values
  expect_equal(as.vector(tbl), c(2, 2, 2)) # 2 of each: 0.3, 0.5, NA
  expect_equal(names(tbl), c("0.3", "0.5", NA))
})

test_that("ps_trim sample() preserves class", {
  set.seed(123)
  ps <- runif(10, 0.05, 0.95)
  x <- ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7)

  # Sample without replacement
  s1 <- sample(x, 5)
  expect_s3_class(s1, "ps_trim")
  expect_equal(length(s1), 5)

  # Sample with replacement
  s2 <- sample(x, 20, replace = TRUE)
  expect_s3_class(s2, "ps_trim")
  expect_equal(length(s2), 20)
})

test_that("ps_trim summary statistics methods work correctly", {
  ps <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  x <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Original: NA, 0.3, 0.5, 0.7, NA

  expect_equal(min(x, na.rm = TRUE), 0.3)
  expect_equal(max(x, na.rm = TRUE), 0.7)
  expect_equal(range(x, na.rm = TRUE), c(0.3, 0.7))
  expect_equal(median(x, na.rm = TRUE), 0.5)
  expect_equal(quantile(x, 0.5, na.rm = TRUE), c("50%" = 0.5))

  # These don't have specific methods but should work through defaults
  expect_equal(mean(x, na.rm = TRUE), 0.5)
  expect_equal(sum(x, na.rm = TRUE), 1.5)
})

test_that("ps_trim diff() works correctly", {
  ps <- c(0.2, 0.3, 0.5, 0.8)
  x <- ps_trim(ps, method = "ps", lower = 0.25, upper = 0.75)
  # Original: NA, 0.3, 0.5, NA

  d <- diff(x)
  expect_type(d, "double")
  expect_equal(length(d), 3)
  # diff of NA, 0.3, 0.5, NA is NA, 0.2, NA
  expect_true(is.na(d[1]))
  expect_equal(d[2], 0.2)
  expect_true(is.na(d[3]))
})

test_that("ps_trim works in data frames", {
  ps <- runif(5, 0.05, 0.95)
  df <- data.frame(
    id = 1:5,
    ps_scores = ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7)
  )

  expect_s3_class(df$ps_scores, "ps_trim")

  # Subsetting preserves class
  expect_s3_class(df$ps_scores[1:3], "ps_trim")
  expect_s3_class(df[1:3, "ps_scores"], "ps_trim")
})

# Additional edge cases ----
test_that("ps_trim handles all values outside bounds", {
  # All below lower bound
  all_low <- ps_trim(
    c(0.01, 0.05, 0.09),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  expect_true(all(is.na(all_low)))
  expect_equal(ps_trim_meta(all_low)$trimmed_idx, 1:3)
  expect_equal(ps_trim_meta(all_low)$keep_idx, integer(0))

  # All above upper bound
  all_high <- ps_trim(
    c(0.85, 0.9, 0.95),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  expect_true(all(is.na(all_high)))
  expect_equal(ps_trim_meta(all_high)$trimmed_idx, 1:3)
  expect_equal(ps_trim_meta(all_high)$keep_idx, integer(0))
})

test_that("ps_trim handles values exactly at boundaries", {
  # Values exactly at boundaries should be kept
  boundary_vals <- c(0.2, 0.5, 0.8)
  trimmed <- ps_trim(boundary_vals, method = "ps", lower = 0.2, upper = 0.8)
  expect_false(any(is.na(trimmed)))
  expect_equal(as.numeric(trimmed), boundary_vals)
  expect_equal(ps_trim_meta(trimmed)$keep_idx, 1:3)

  # Mix of boundary and out-of-bounds values
  mixed <- c(0.1, 0.2, 0.5, 0.8, 0.9)
  trimmed_mixed <- ps_trim(mixed, method = "ps", lower = 0.2, upper = 0.8)
  expect_equal(is.na(trimmed_mixed), c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(ps_trim_meta(trimmed_mixed)$trimmed_idx, c(1, 5))
  expect_equal(ps_trim_meta(trimmed_mixed)$keep_idx, c(2, 3, 4))
})

test_that("ps_trim handles edge cases for method = 'cr'", {
  # Test with no overlap in common range
  ps_vals <- c(0.1, 0.2, 0.8, 0.9)
  exposure <- c(0, 0, 1, 1)

  # Control: [0.1, 0.2], Treated: [0.8, 0.9]
  # Common range would be [0.8, 0.2] which is invalid - all should be trimmed
  trimmed_cr <- ps_trim(ps_vals, method = "cr", .exposure = exposure)

  expect_true(all(is.na(trimmed_cr)))
  expect_equal(ps_trim_meta(trimmed_cr)$keep_idx, integer(0))
  expect_equal(ps_trim_meta(trimmed_cr)$trimmed_idx, 1:4)

  # Test with valid common range
  ps_vals2 <- c(0.3, 0.4, 0.6, 0.5, 0.7)
  exposure2 <- c(0, 0, 0, 1, 1)

  # Control: [0.3, 0.4, 0.6], Treated: [0.5, 0.7]
  # Common range: [0.5, 0.6]
  trimmed_cr2 <- ps_trim(ps_vals2, method = "cr", .exposure = exposure2)

  # Values outside [0.5, 0.6] should be trimmed
  expect_true(is.na(trimmed_cr2[1])) # 0.3 < 0.5
  expect_true(is.na(trimmed_cr2[2])) # 0.4 < 0.5
  expect_false(is.na(trimmed_cr2[3])) # 0.6 in range
  expect_false(is.na(trimmed_cr2[4])) # 0.5 in range
  expect_true(is.na(trimmed_cr2[5])) # 0.7 > 0.6
})
