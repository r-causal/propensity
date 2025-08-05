# Test base R methods for psw class
library(testthat)
library(vctrs)

test_that("psw subsetting with [ preserves class and attributes", {
  x <- psw(c(0.1, 0.2, 0.3, 0.4, 0.5), estimand = "ate", stabilized = TRUE)

  # Single element
  expect_s3_class(x[1], "psw")
  expect_equal(estimand(x[1]), "ate")
  expect_true(is_stabilized(x[1]))
  expect_equal(length(x[1]), 1)
  expect_equal(as.numeric(x[1]), 0.1)

  # Multiple elements
  expect_s3_class(x[c(2, 4)], "psw")
  expect_equal(estimand(x[c(2, 4)]), "ate")
  expect_true(is_stabilized(x[c(2, 4)]))
  expect_equal(length(x[c(2, 4)]), 2)
  expect_equal(as.numeric(x[c(2, 4)]), c(0.2, 0.4))

  # Logical subsetting
  expect_s3_class(x[c(TRUE, FALSE, TRUE, FALSE, TRUE)], "psw")
  expect_equal(length(x[c(TRUE, FALSE, TRUE, FALSE, TRUE)]), 3)

  # Empty subsetting
  expect_s3_class(x[integer(0)], "psw")
  expect_equal(length(x[integer(0)]), 0)
  expect_equal(estimand(x[integer(0)]), "ate")

  # Negative indices
  expect_s3_class(x[-c(1, 5)], "psw")
  expect_equal(length(x[-c(1, 5)]), 3)
  expect_equal(as.numeric(x[-c(1, 5)]), c(0.2, 0.3, 0.4))
})

test_that("psw sort() preserves class and attributes", {
  x <- psw(c(0.3, 0.1, 0.4, 0.2), estimand = "att", stabilized = FALSE)

  sorted <- sort(x)
  expect_s3_class(sorted, "psw")
  expect_equal(estimand(sorted), "att")
  expect_false(is_stabilized(sorted))
  expect_equal(as.numeric(sorted), c(0.1, 0.2, 0.3, 0.4))

  # Decreasing order
  sorted_dec <- sort(x, decreasing = TRUE)
  expect_s3_class(sorted_dec, "psw")
  expect_equal(as.numeric(sorted_dec), c(0.4, 0.3, 0.2, 0.1))

  # With NAs
  x_na <- psw(c(0.3, NA, 0.1, 0.4, NA, 0.2), estimand = "att")
  sorted_na <- sort(x_na, na.last = TRUE)
  expect_s3_class(sorted_na, "psw")
  expect_equal(as.numeric(sorted_na), c(0.1, 0.2, 0.3, 0.4, NA, NA))
})

test_that("psw unique() preserves class and attributes", {
  x <- psw(c(0.1, 0.2, 0.1, 0.3, 0.2, 0.1), estimand = "ato")

  uniq <- unique(x)
  expect_s3_class(uniq, "psw")
  expect_equal(estimand(uniq), "ato")
  expect_equal(as.numeric(uniq), c(0.1, 0.2, 0.3))

  # With NAs
  x_na <- psw(c(0.1, NA, 0.2, NA, 0.1), estimand = "ato")
  uniq_na <- unique(x_na)
  expect_s3_class(uniq_na, "psw")
  expect_equal(as.numeric(uniq_na), c(0.1, NA, 0.2))
})

test_that("psw duplicated() returns logical vector", {
  x <- psw(c(0.1, 0.2, 0.1, 0.3, 0.2, 0.1), estimand = "ate")

  dups <- duplicated(x)
  expect_type(dups, "logical")
  expect_equal(dups, c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE))

  # anyDuplicated already tested in test-psw.R
  expect_equal(anyDuplicated(x), 3)
})

test_that("psw rev() preserves class and attributes", {
  x <- psw(c(0.1, 0.2, 0.3, 0.4), estimand = "att", stabilized = TRUE)

  reversed <- rev(x)
  expect_s3_class(reversed, "psw")
  expect_equal(estimand(reversed), "att")
  expect_true(is_stabilized(reversed))
  expect_equal(as.numeric(reversed), c(0.4, 0.3, 0.2, 0.1))
})

test_that("psw head() and tail() preserve class and attributes", {
  x <- psw(1:10 / 10, estimand = "ate")

  # head
  h <- head(x, 3)
  expect_s3_class(h, "psw")
  expect_equal(estimand(h), "ate")
  expect_equal(length(h), 3)
  expect_equal(as.numeric(h), c(0.1, 0.2, 0.3))

  # tail
  t <- tail(x, 3)
  expect_s3_class(t, "psw")
  expect_equal(estimand(t), "ate")
  expect_equal(length(t), 3)
  expect_equal(as.numeric(t), c(0.8, 0.9, 1.0))

  # Default n = 6
  expect_equal(length(head(x)), 6)
  expect_equal(length(tail(x)), 6)
})

test_that("psw rep() preserves class and attributes", {
  x <- psw(c(0.1, 0.2), estimand = "att")

  # times argument
  r1 <- rep(x, times = 3)
  expect_s3_class(r1, "psw")
  expect_equal(estimand(r1), "att")
  expect_equal(length(r1), 6)
  expect_equal(as.numeric(r1), c(0.1, 0.2, 0.1, 0.2, 0.1, 0.2))

  # each argument
  r2 <- rep(x, each = 2)
  expect_s3_class(r2, "psw")
  expect_equal(estimand(r2), "att")
  expect_equal(length(r2), 4)
  expect_equal(as.numeric(r2), c(0.1, 0.1, 0.2, 0.2))

  # length.out argument
  r3 <- rep(x, length.out = 5)
  expect_s3_class(r3, "psw")
  expect_equal(length(r3), 5)
  expect_equal(as.numeric(r3), c(0.1, 0.2, 0.1, 0.2, 0.1))
})

test_that("psw is.na() and anyNA() work correctly", {
  x <- psw(c(0.1, NA, 0.3, NA, 0.5), estimand = "ate")

  na_vec <- is.na(x)
  expect_type(na_vec, "logical")
  expect_equal(na_vec, c(FALSE, TRUE, FALSE, TRUE, FALSE))

  expect_true(anyNA(x))

  x_no_na <- psw(c(0.1, 0.2, 0.3), estimand = "ate")
  expect_false(anyNA(x_no_na))
})

test_that("psw na.omit() preserves class and attributes", {
  x <- psw(c(0.1, NA, 0.3, NA, 0.5), estimand = "ato", stabilized = TRUE)

  x_clean <- na.omit(x)
  expect_s3_class(x_clean, "psw")
  expect_equal(estimand(x_clean), "ato")
  expect_true(is_stabilized(x_clean))
  expect_equal(as.numeric(x_clean), c(0.1, 0.3, 0.5))
  expect_equal(length(x_clean), 3)
})

test_that("psw is.finite() and is.infinite() work correctly", {
  x <- psw(c(0.1, Inf, 0.3, -Inf, 0.5), estimand = "ate")

  finite_vec <- is.finite(x)
  expect_type(finite_vec, "logical")
  expect_equal(finite_vec, c(TRUE, FALSE, TRUE, FALSE, TRUE))

  infinite_vec <- is.infinite(x)
  expect_type(infinite_vec, "logical")
  expect_equal(infinite_vec, c(FALSE, TRUE, FALSE, TRUE, FALSE))
})

test_that("psw which() family functions work correctly", {
  x <- psw(c(0.1, 0.5, 0.2, 0.8, 0.3), estimand = "ate")

  # which with logical condition
  idx <- which(as.numeric(x) > 0.3)
  expect_type(idx, "integer")
  expect_equal(idx, c(2, 4))

  # which.min
  expect_equal(which.min(x), 1)

  # which.max
  expect_equal(which.max(x), 4)
})

test_that("psw order() and rank() work correctly", {
  x <- psw(c(0.3, 0.1, 0.4, 0.2), estimand = "ate")

  # order
  ord <- order(x)
  expect_type(ord, "integer")
  expect_equal(ord, c(2, 4, 1, 3))

  # rank
  rnk <- rank(x)
  expect_type(rnk, "double")
  expect_equal(rnk, c(3, 1, 4, 2))
})

test_that("psw match() and %in% work correctly", {
  x <- psw(c(0.1, 0.2, 0.3), estimand = "ate")
  y <- psw(c(0.2, 0.4), estimand = "ate")

  # match
  m <- match(x, y)
  expect_type(m, "integer")
  expect_equal(m, c(NA, 1, NA))

  # %in%
  inn <- x %in% y
  expect_type(inn, "logical")
  expect_equal(inn, c(FALSE, TRUE, FALSE))
})

test_that("psw table() works correctly", {
  x <- psw(c(0.1, 0.2, 0.1, 0.3, 0.2, 0.1), estimand = "ate")

  tbl <- table(x)
  expect_s3_class(tbl, "table")
  expect_equal(as.numeric(tbl), c(3, 2, 1))
  expect_equal(names(tbl), c("0.1", "0.2", "0.3"))
})

test_that("psw sample() preserves class and attributes", {
  set.seed(123)
  x <- psw(c(0.1, 0.2, 0.3, 0.4, 0.5), estimand = "att")

  # Sample without replacement
  s1 <- sample(x, 3)
  expect_s3_class(s1, "psw")
  expect_equal(estimand(s1), "att")
  expect_equal(length(s1), 3)
  expect_true(all(s1 %in% x))

  # Sample with replacement
  s2 <- sample(x, 10, replace = TRUE)
  expect_s3_class(s2, "psw")
  expect_equal(length(s2), 10)
})

test_that("psw summary statistics methods work correctly", {
  x <- psw(c(0.1, 0.2, 0.3, 0.4, 0.5), estimand = "ate")

  # These are already implemented and tested in test-psw.R
  expect_equal(min(x), 0.1)
  expect_equal(max(x), 0.5)
  expect_equal(range(x), c(0.1, 0.5))
  expect_equal(median(x), 0.3)
  expect_equal(quantile(x, 0.5), c("50%" = 0.3))

  # summary
  s <- summary(x)
  expect_type(s, "double")
  expect_true("Min." %in% names(s))
  expect_true("Max." %in% names(s))
  expect_true("Median" %in% names(s))
})

test_that("psw works in data frames", {
  df <- data.frame(
    id = 1:5,
    weights = psw(c(0.1, 0.2, 0.3, 0.4, 0.5), estimand = "ate")
  )

  expect_s3_class(df$weights, "psw")
  expect_equal(estimand(df$weights), "ate")

  # Subsetting preserves class
  expect_s3_class(df$weights[1:3], "psw")
  expect_s3_class(df[1:3, "weights"], "psw")
})

# Edge cases ----
test_that("psw handles empty vectors", {
  empty_psw <- psw(numeric(0), estimand = "ate")
  expect_s3_class(empty_psw, "psw")
  expect_equal(length(empty_psw), 0)
  expect_equal(estimand(empty_psw), "ate")

  # Operations on empty psw
  expect_equal(sum(empty_psw), 0)
  expect_true(is.na(mean(empty_psw))) # mean of empty numeric is NA
  expect_equal(length(summary(empty_psw)), 6)
})

test_that("psw handles single element vectors", {
  single_psw <- psw(0.5, estimand = "att")
  expect_s3_class(single_psw, "psw")
  expect_equal(length(single_psw), 1)
  expect_equal(as.numeric(single_psw), 0.5)

  # Subsetting single element
  expect_equal(single_psw[1], single_psw)
  expect_equal(length(single_psw[integer(0)]), 0)
})

test_that("psw handles NA values correctly", {
  na_psw <- psw(c(0.1, NA, 0.3, NA, 0.5), estimand = "ate")
  expect_s3_class(na_psw, "psw")
  expect_equal(sum(is.na(na_psw)), 2)

  # Operations with NAs
  expect_equal(sum(na_psw), NA_real_)
  expect_equal(sum(na_psw, na.rm = TRUE), 0.9)
  expect_equal(mean(na_psw, na.rm = TRUE), 0.3)

  # na.omit
  clean <- na.omit(na_psw)
  expect_equal(length(clean), 3)
  expect_false(anyNA(clean))
})

test_that("Combining psw objects with incompatible metadata warns", {
  psw1 <- psw(c(0.1, 0.2), estimand = "ate")
  psw2 <- psw(c(0.3, 0.4), estimand = "att")

  expect_warning(
    combined <- c(psw1, psw2),
    "incompatible estimands"
  )
  expect_type(combined, "double")
})
