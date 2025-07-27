test_that("cut_quantile() - Basic functionality", {
  set.seed(123)
  x <- rnorm(100)
  
  # Default call should create deciles
  result <- cut_quantile(x)
  
  # Check basic properties
  expect_s3_class(result, "factor")
  expect_equal(length(result), length(x))
  
  # Should have 10 levels by default (deciles)
  expect_equal(nlevels(result), 10)
  
  # Check that it's equivalent to manual cut with quantiles
  manual_breaks <- quantile(x, probs = seq(0, 1, 0.1), na.rm = TRUE)
  manual_result <- cut(x, breaks = manual_breaks, include.lowest = TRUE)
  expect_equal(result, manual_result)
})

test_that("cut_quantile() - Custom quantiles", {
  set.seed(456)
  x <- rnorm(50)
  
  # Test quartiles
  result_quartiles <- cut_quantile(x, probs = seq(0, 1, 0.25))
  expect_equal(nlevels(result_quartiles), 4)
  
  # Test tertiles
  result_tertiles <- cut_quantile(x, probs = c(0, 1/3, 2/3, 1))
  expect_equal(nlevels(result_tertiles), 3)
  
  # Test with just two quantiles (median split)
  result_median <- cut_quantile(x, probs = c(0, 0.5, 1))
  expect_equal(nlevels(result_median), 2)
})

test_that("cut_quantile() - Error handling", {
  # Non-numeric input
  expect_error(
    cut_quantile(letters[1:10]), 
    "x.*must be numeric"
  )
  
  # Invalid probs - outside [0,1]
  expect_error(
    cut_quantile(1:10, probs = c(-0.1, 0.5, 1)), 
    "probs.*must be numeric values between 0 and 1"
  )
  
  expect_error(
    cut_quantile(1:10, probs = c(0, 0.5, 1.1)), 
    "probs.*must be numeric values between 0 and 1"
  )
  
  # Non-numeric probs
  expect_error(
    cut_quantile(1:10, probs = c("0", "0.5", "1")), 
    "probs.*must be numeric values between 0 and 1"
  )
  
  # Too few probs
  expect_error(
    cut_quantile(1:10, probs = 0.5), 
    "probs.*must have at least 2 values"
  )
  
  # All NA values
  expect_error(
    cut_quantile(rep(NA_real_, 10)), 
    "All values.*are NA"
  )
})

test_that("cut_quantile() - Handling NA values", {
  # Test with some NA values
  x_with_na <- c(1, 2, 3, NA, 5, 6, 7, NA, 9, 10)
  result <- cut_quantile(x_with_na, probs = c(0, 0.5, 1))
  
  # Should preserve NA positions
  expect_true(is.na(result[4]))
  expect_true(is.na(result[8]))
  
  # Non-NA values should be properly categorized
  expect_false(any(is.na(result[c(1, 2, 3, 5, 6, 7, 9, 10)])))
})

test_that("cut_quantile() - Additional arguments passed to cut()", {
  set.seed(789)
  x <- rnorm(20)
  
  # Test passing labels argument
  custom_labels <- c("Low", "High")
  result <- cut_quantile(x, probs = c(0, 0.5, 1), labels = custom_labels)
  expect_equal(levels(result), custom_labels)
  
  # Test passing right argument
  result_right_false <- cut_quantile(x, probs = c(0, 0.5, 1), right = FALSE)
  result_right_true <- cut_quantile(x, probs = c(0, 0.5, 1), right = TRUE)
  
  # Results should be different when right argument changes
  # (though specific difference depends on data)
  expect_s3_class(result_right_false, "factor")
  expect_s3_class(result_right_true, "factor")
})

test_that("cut_quantile() - Edge cases with identical values", {
  # Test with all identical values
  x_identical <- rep(5, 10)
  
  # This should warn and create a single interval
  expect_warning(
    result <- cut_quantile(x_identical, probs = c(0, 0.5, 1)),
    "All values are identical"
  )
  expect_s3_class(result, "factor")
  expect_equal(length(result), 10)
  expect_equal(nlevels(result), 1)
  
  # Test with mostly identical values that still have variation
  x_mostly_same <- c(rep(5, 8), 4, 6)
  result2 <- cut_quantile(x_mostly_same, probs = c(0, 0.5, 1))
  expect_s3_class(result2, "factor")
  expect_equal(length(result2), 10)
  # Should have at least 1 level, possibly 2 depending on data
  expect_gte(nlevels(result2), 1)
})