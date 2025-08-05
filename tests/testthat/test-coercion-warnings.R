test_that("c() with psw objects of different estimands warns and returns correct numeric values", {
  x <- psw(c(0.5, 0.7), estimand = "ate")
  y <- psw(c(0.3, 0.8), estimand = "att")

  # Test c() function specifically
  expect_warning(
    result <- c(x, y),
    "incompatible estimands 'ate' and 'att'"
  )
  expect_type(result, "double")
  expect_equal(result, c(0.5, 0.7, 0.3, 0.8))
  expect_equal(length(result), 4)

  # Test with more values - note: when mixing with double, different warning
  z <- psw(c(0.1, 0.2, 0.9), estimand = "ato")
  expect_warning(
    result2 <- c(x, y, z),
    "Converting psw to numeric"
  )
  expect_equal(result2, c(0.5, 0.7, 0.3, 0.8, 0.1, 0.2, 0.9))
  expect_equal(length(result2), 7)
})

test_that("c() with psw objects of same estimand combines without warning", {
  x <- psw(c(0.5, 0.7), estimand = "ate")
  y <- psw(c(0.3, 0.8), estimand = "ate")
  z <- psw(c(0.1, 0.9), estimand = "ate")

  expect_silent(result <- c(x, y, z))
  expect_s3_class(result, "psw")
  expect_equal(estimand(result), "ate")
  expect_equal(as.numeric(result), c(0.5, 0.7, 0.3, 0.8, 0.1, 0.9))
  expect_equal(length(result), 6)
})

test_that("c() with psw and numeric values warns and returns numeric", {
  x <- psw(c(0.5, 0.7), estimand = "ate")

  # psw with single numeric
  expect_warning(
    result <- c(x, 0.9),
    "Converting psw to numeric"
  )
  expect_type(result, "double")
  expect_equal(result, c(0.5, 0.7, 0.9))

  # numeric first - no warning because numeric method is called
  result2 <- c(0.1, x)
  expect_type(result2, "double")
  expect_equal(result2, c(0.1, 0.5, 0.7))

  # psw with multiple numerics
  expect_warning(
    result3 <- c(x, c(0.2, 0.3)),
    "Converting psw to numeric"
  )
  expect_equal(result3, c(0.5, 0.7, 0.2, 0.3))

  # Mixed order - numeric first means no warning
  result4 <- c(0.1, x, 0.9, c(0.2, 0.3))
  expect_type(result4, "double")
  expect_equal(result4, c(0.1, 0.5, 0.7, 0.9, 0.2, 0.3))
})

test_that("c() with ps_trim objects of different parameters warns and returns numeric", {
  x <- ps_trim(c(0.1, 0.5, 0.9), lower = 0.1, upper = 0.9)
  y <- ps_trim(c(0.2, 0.6, 0.8), lower = 0.2, upper = 0.8)

  expect_warning(
    result <- c(x, y),
    "different trimming parameters"
  )
  expect_type(result, "double")
  # ps_trim might have NAs for trimmed values
  expect_equal(length(result), 6)
  expect_true(all(result[!is.na(result)] >= 0 & result[!is.na(result)] <= 1))
})

test_that("c() with ps_trim objects of same parameters combines correctly", {
  x <- ps_trim(c(0.2, 0.5, 0.8), method = "ps", lower = 0.1, upper = 0.9)
  y <- ps_trim(c(0.3, 0.6, 0.7), method = "ps", lower = 0.1, upper = 0.9)

  expect_silent(result <- c(x, y))
  expect_s3_class(result, "ps_trim")
  expect_equal(length(result), 6)
  # Values should be preserved (no NAs since all are within bounds)
  expect_equal(
    as.numeric(result[!is.na(result)]),
    c(0.2, 0.5, 0.8, 0.3, 0.6, 0.7)[!is.na(c(x, y))]
  )
})

test_that("c() with ps_trunc objects behaves correctly", {
  # Different parameters
  x <- ps_trunc(c(0.1, 0.5, 0.9), lower = 0.2, upper = 0.8)
  y <- ps_trunc(c(0.15, 0.6, 0.85), lower = 0.3, upper = 0.7)

  expect_warning(
    result <- c(x, y),
    "different truncation parameters"
  )
  expect_type(result, "double")
  expect_equal(result, c(0.2, 0.5, 0.8, 0.3, 0.6, 0.7)) # truncated values

  # Same parameters
  x2 <- ps_trunc(c(0.1, 0.5, 0.9), lower = 0.2, upper = 0.8)
  y2 <- ps_trunc(c(0.15, 0.6, 0.85), lower = 0.2, upper = 0.8)

  expect_silent(result2 <- c(x2, y2))
  expect_s3_class(result2, "ps_trunc")
  expect_equal(as.numeric(result2), c(0.2, 0.5, 0.8, 0.2, 0.6, 0.8))
})

test_that("c() with mixed propensity classes warns and returns numeric", {
  psw_obj <- psw(c(0.5, 0.7), estimand = "ate")
  trim_obj <- ps_trim(c(0.3, 0.8), lower = 0.1, upper = 0.9)
  trunc_obj <- ps_trunc(c(0.4, 0.6), lower = 0.1, upper = 0.9)

  # psw + ps_trim
  expect_warning(
    result1 <- c(psw_obj, trim_obj),
    "Converting psw and ps_trim to numeric"
  )
  expect_type(result1, "double")
  expect_equal(length(result1), 4)

  # psw + ps_trunc
  expect_warning(
    result2 <- c(psw_obj, trunc_obj),
    "Converting psw and ps_trunc to numeric"
  )
  expect_type(result2, "double")
  expect_equal(result2, c(0.5, 0.7, 0.4, 0.6))

  # ps_trim + ps_trunc
  expect_warning(
    result3 <- c(trim_obj, trunc_obj),
    "Converting ps_trim and ps_trunc to numeric"
  )
  expect_type(result3, "double")
  expect_equal(length(result3), 4)

  # All three - gets multiple warnings
  expect_warning(
    result4 <- c(psw_obj, trim_obj, trunc_obj),
    "Converting ps_trunc to numeric"
  )
  expect_type(result4, "double")
  expect_equal(length(result4), 6)
})

test_that("c() with empty vectors works correctly", {
  x <- psw(c(0.5, 0.7), estimand = "ate")
  empty_psw <- psw(double(), estimand = "ate")
  empty_numeric <- double()

  # psw with empty psw
  result1 <- c(x, empty_psw)
  expect_s3_class(result1, "psw")
  expect_equal(as.numeric(result1), c(0.5, 0.7))

  # psw with empty numeric - warns in current implementation
  expect_warning(
    result2 <- c(x, empty_numeric),
    "Converting psw to numeric"
  )
  expect_type(result2, "double")
  expect_equal(result2, c(0.5, 0.7))

  # empty psw with numeric
  expect_warning(
    result3 <- c(empty_psw, 0.5),
    "Converting psw to numeric"
  )
  expect_equal(result3, 0.5)
})

test_that("c() with single values works correctly", {
  x <- psw(0.5, estimand = "ate")
  y <- psw(0.7, estimand = "att")

  expect_warning(
    result <- c(x, y),
    "incompatible estimands"
  )
  expect_equal(result, c(0.5, 0.7))

  # Single with vector
  z <- psw(c(0.1, 0.2), estimand = "ato")
  expect_warning(
    result2 <- c(x, z),
    "incompatible estimands"
  )
  expect_equal(result2, c(0.5, 0.1, 0.2))
})

test_that("c() preserves attributes when metadata matches", {
  x <- new_psw(c(0.5, 0.7), estimand = "ate", stabilized = TRUE)
  y <- new_psw(c(0.3, 0.8), estimand = "ate", stabilized = TRUE)

  result <- c(x, y)
  expect_s3_class(result, "psw")
  expect_equal(estimand(result), "ate")
  expect_true(is_stabilized(result))
  expect_equal(as.numeric(result), c(0.5, 0.7, 0.3, 0.8))
})

test_that("c() with different stabilization status warns", {
  x <- new_psw(c(0.5, 0.7), estimand = "ate", stabilized = TRUE)
  y <- new_psw(c(0.3, 0.8), estimand = "ate", stabilized = FALSE)

  expect_warning(
    result <- c(x, y),
    "different stabilization status"
  )
  expect_type(result, "double")
  expect_equal(result, c(0.5, 0.7, 0.3, 0.8))
})

test_that("subsetting operations work correctly", {
  x <- psw(c(0.1, 0.5, 0.7, 0.9), estimand = "ate")
  y <- psw(c(0.2, 0.6), estimand = "att")

  # Subset then combine
  expect_warning(
    result <- c(x[1:2], y),
    "incompatible estimands"
  )
  expect_equal(result, c(0.1, 0.5, 0.2, 0.6))

  # Combine subsets
  expect_warning(
    result2 <- c(x[c(1, 3)], y[2]),
    "incompatible estimands"
  )
  expect_equal(result2, c(0.1, 0.7, 0.6))
})

test_that("append() works like c()", {
  x <- psw(c(0.5, 0.7), estimand = "ate")
  y <- psw(c(0.3, 0.8), estimand = "att")

  expect_warning(
    result <- append(x, y),
    "incompatible estimands"
  )
  expect_type(result, "double")
  expect_equal(result, c(0.5, 0.7, 0.3, 0.8))

  # With after argument - more complex due to subsetting
  expect_warning(
    result2 <- append(x, y, after = 1),
    "Converting psw to numeric"
  )
  expect_equal(result2, c(0.5, 0.3, 0.8, 0.7))
})

test_that("unlist() works correctly", {
  x <- psw(c(0.5, 0.7), estimand = "ate")
  y <- psw(c(0.3, 0.8), estimand = "att")

  lst <- list(a = x, b = y)

  # unlist doesn't go through vctrs so no warning
  result <- unlist(lst)
  expect_type(result, "double")
  expect_equal(as.numeric(result), c(0.5, 0.7, 0.3, 0.8))
  expect_equal(names(result), c("a1", "a2", "b1", "b2"))
})

test_that("data.frame operations work as expected", {
  df1 <- data.frame(
    id = 1:3,
    wt = psw(c(0.5, 0.7, 0.3), estimand = "ate")
  )

  df2 <- data.frame(
    id = 4:6,
    wt = psw(c(0.8, 0.2, 0.6), estimand = "att")
  )

  # rbind maintains psw class when estimands match for first data frame
  # but converts to numeric when they don't match
  result1 <- rbind(df1, df2)
  expect_equal(nrow(result1), 6)
  # rbind preserves the first object's class
  expect_s3_class(result1$wt, "psw")
  expect_equal(as.numeric(result1$wt), c(0.5, 0.7, 0.3, 0.8, 0.2, 0.6))

  # But vec_rbind does trigger the warning
  expect_warning(
    result2 <- vctrs::vec_rbind(df1, df2),
    "incompatible estimands"
  )
  expect_equal(nrow(result2), 6)
  expect_type(result2$wt, "double")
  expect_equal(result2$wt, c(0.5, 0.7, 0.3, 0.8, 0.2, 0.6))
})

test_that("vctrs vec_ptype2 returns appropriate prototypes", {
  x <- psw(c(0.5, 0.7), estimand = "ate")
  y <- psw(c(0.3, 0.8), estimand = "att")

  # This is what vctrs uses internally - verify it returns empty double
  expect_warning(
    proto <- vec_ptype2(x, y),
    class = "propensity_coercion_warning"
  )
  expect_identical(proto, double())
  expect_equal(length(proto), 0)

  # Compatible objects return psw prototype
  z <- psw(c(0.1, 0.2), estimand = "ate")
  expect_silent(proto2 <- vec_ptype2(x, z))
  expect_s3_class(proto2, "psw")
  expect_equal(length(proto2), 0)
})

test_that("arithmetic operations with different metadata work correctly", {
  x <- psw(c(0.5, 0.7), estimand = "ate")
  y <- psw(c(0.3, 0.8), estimand = "att")

  # These should create psw objects but with combined estimand info
  result <- x + y
  expect_s3_class(result, "psw")
  expect_equal(estimand(result), "ate, att")
  expect_equal(as.numeric(result), c(0.8, 1.5))

  result2 <- x * y
  expect_s3_class(result2, "psw")
  expect_equal(as.numeric(result2), c(0.15, 0.56))
})

test_that("comparison operations warn about class downgrade", {
  x <- psw(c(0.5, 0.7), estimand = "ate")

  expect_warning(
    result <- x > 0.6,
    "Converting psw to numeric"
  )
  expect_equal(result, c(FALSE, TRUE))

  expect_warning(
    result2 <- x == 0.5,
    "Converting psw to numeric"
  )
  expect_equal(result2, c(TRUE, FALSE))
})

test_that("interaction with base R functions work correctly", {
  x <- psw(c(0.5, 0.7, 0.3), estimand = "ate")

  # sum, mean, etc should return numeric
  expect_equal(sum(x), 1.5)
  expect_equal(mean(x), 0.5)
  expect_equal(min(x), 0.3)
  expect_equal(max(x), 0.7)

  # These return numeric as expected
  expect_type(sum(x), "double")
  expect_type(mean(x), "double")
})

test_that("c() ordering matters for warnings", {
  x <- psw(c(0.5, 0.7), estimand = "ate")
  num <- c(0.3, 0.4)

  # psw first triggers warning
  expect_warning(
    result1 <- c(x, num),
    "Converting psw to numeric"
  )
  expect_equal(result1, c(0.5, 0.7, 0.3, 0.4))

  # numeric first doesn't trigger our warning (uses base c())
  result2 <- c(num, x)
  expect_type(result2, "double")
  expect_equal(result2, c(0.3, 0.4, 0.5, 0.7))
})
