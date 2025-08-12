test_that("ps_trunc works with matrix propensity scores", {
  # Create test data
  set.seed(123)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.02, 0.8), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix) # Normalize to sum to 1
  colnames(ps_matrix) <- levels(exposure)

  # Test basic truncation
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )

  expect_s3_class(truncated, c("ps_trunc_matrix", "ps_trunc", "matrix"))
  expect_equal(dim(truncated), dim(ps_matrix))

  # Check metadata
  meta <- ps_trunc_meta(truncated)
  expect_equal(meta$method, "ps")
  expect_equal(meta$lower_bound, 0.05)
  expect_true(meta$is_matrix)

  # After truncation and renormalization, some values might still be below threshold
  # This is expected behavior - we truncate once and renormalize

  # Check that rows still sum to 1
  row_sums <- rowSums(truncated)
  expect_equal(row_sums, rep(1, n), tolerance = 1e-10)

  # Check which rows were truncated
  had_low_vals <- apply(ps_matrix, 1, function(x) any(x < 0.05))
  expect_equal(sort(meta$truncated_idx), sort(which(had_low_vals)))
})

test_that("ps_trunc works with data.frame propensity scores", {
  set.seed(456)
  n <- 50
  exposure <- factor(sample(c("Low", "Med", "High"), n, replace = TRUE))

  ps_df <- data.frame(
    Low = runif(n, 0.01, 0.6),
    Med = runif(n, 0.02, 0.5),
    High = runif(n, 0.01, 0.4)
  )
  # Normalize rows
  ps_df <- ps_df / rowSums(ps_df)

  truncated <- ps_trunc(ps_df, .exposure = exposure, method = "ps", lower = 0.1)

  expect_s3_class(truncated, c("ps_trunc_matrix", "ps_trunc", "matrix"))
  expect_equal(nrow(truncated), nrow(ps_df))
  expect_equal(ncol(truncated), ncol(ps_df))

  # After truncation and renormalization, some values might still be below threshold
  # This is expected behavior
})

test_that("ps_trunc preserves row sums equal to 1", {
  set.seed(789)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C", "D"), n, replace = TRUE))

  # Create propensity scores with some extreme values
  ps_matrix <- matrix(runif(n * 4, 0.001, 0.8), nrow = n, ncol = 4)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Apply truncation
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )

  # Check all rows sum to 1
  row_sums <- rowSums(truncated)
  expect_true(all(abs(row_sums - 1) < 1e-10))

  # Original and truncated should have same dimensions
  expect_equal(dim(truncated), dim(ps_matrix))
})

test_that("ps_trunc validates delta < 1/k", {
  n <- 50
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # delta >= 1/3 should error
  expect_propensity_error(
    ps_trunc(ps_matrix, .exposure = exposure, method = "ps", lower = 0.35)
  )
})

test_that("ps_trunc errors for unsupported methods with categorical", {
  n <- 20
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Only cr should error - pctl is now supported
  expect_propensity_error(
    ps_trunc(ps_matrix, .exposure = exposure, method = "cr")
  )
})

test_that("ps_trunc requires exposure for categorical", {
  n <- 20
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c("A", "B", "C")

  expect_propensity_error(
    ps_trunc(ps_matrix, method = "ps")
  )
})

test_that("is_ps_truncated works for matrix objects", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  expect_false(is_ps_truncated(ps_matrix))

  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )
  expect_true(is_ps_truncated(truncated))
})

test_that("is_unit_truncated works for matrix objects", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.1
  )
  unit_truncated <- is_unit_truncated(truncated)

  expect_equal(length(unit_truncated), n)
  expect_type(unit_truncated, "logical")

  meta <- ps_trunc_meta(truncated)
  expect_equal(which(unit_truncated), meta$truncated_idx)
})

test_that("ps_trunc handles extreme values correctly", {
  n <- 20
  exposure <- factor(c(rep("A", 10), rep("B", 5), rep("C", 5)))

  # Create matrix with very small values
  ps_matrix <- matrix(0, nrow = n, ncol = 3)
  ps_matrix[1:10, ] <- c(0.001, 0.499, 0.5) # Group A
  ps_matrix[11:15, ] <- c(0.5, 0.001, 0.499) # Group B
  ps_matrix[16:20, ] <- c(0.499, 0.5, 0.001) # Group C
  colnames(ps_matrix) <- c("A", "B", "C")

  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.01
  )

  # After truncation and renormalization, some values might still be below threshold
  # This is expected behavior

  # Check rows sum to 1
  expect_equal(rowSums(truncated), rep(1, n), tolerance = 1e-10)

  # All rows should have been truncated
  meta <- ps_trunc_meta(truncated)
  expect_equal(length(meta$truncated_idx), n)
})

test_that("ps_trunc handles parsnip-style column names", {
  n <- 40
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c(".pred_A", ".pred_B", ".pred_C")

  expect_no_error(
    truncated <- ps_trunc(
      ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.05
    )
  )

  expect_s3_class(truncated, "ps_trunc_matrix")
})

test_that("ps_trunc warns when no column names provided", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  # No column names

  expect_propensity_warning(
    truncated <- ps_trunc(
      ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.05
    )
  )
})

test_that("ps_trunc.ps_trunc warns about already truncated scores", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  truncated_once <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )

  expect_propensity_warning(
    truncated_twice <- ps_trunc(
      truncated_once,
      .exposure = exposure,
      method = "ps",
      lower = 0.1
    )
  )

  # Should return original truncated object
  expect_identical(truncated_twice, truncated_once)
})

# PSweight comparison tests
test_that("ps_trunc matches PSweight truncation behavior for entropy weights", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # PSweight automatically truncates for entropy weights
  # It uses fixed bounds [1e-6, 1-1e-6]

  set.seed(123)
  n <- 100

  # Create extreme propensity scores that sum to 1
  ps_matrix <- matrix(0, nrow = n, ncol = 3)
  # Use rep to fill rows properly
  ps_matrix[1:30, ] <- matrix(
    rep(c(1e-7, 0.5, 0.5 - 1e-7), each = 30),
    ncol = 3,
    byrow = TRUE
  )
  ps_matrix[31:60, ] <- matrix(
    rep(c(0.3, 0.3, 0.4), each = 30),
    ncol = 3,
    byrow = TRUE
  )
  ps_matrix[61:100, ] <- matrix(
    rep(c(0.99999, 5e-6, 5e-6), each = 40),
    ncol = 3,
    byrow = TRUE
  )
  # Normalize to ensure rows sum to 1
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c("A", "B", "C")

  exposure <- factor(c(rep("A", 30), rep("B", 30), rep("C", 40)))

  # Our truncation with PSweight's threshold
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 1e-6
  )

  # Check bounds match PSweight's approach (excluding NAs)
  expect_true(all(truncated >= 1e-6, na.rm = TRUE))
  expect_true(all(truncated <= 1 - 1e-6, na.rm = TRUE))

  # Rows should still sum to 1
  expect_equal(rowSums(truncated), rep(1, n), tolerance = 1e-10)
})

test_that("ps_trunc with method='pctl' works for categorical exposures", {
  set.seed(789)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  # Create propensity scores with some extreme values
  ps_matrix <- matrix(runif(n * 3, 0.001, 0.999), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Apply percentile truncation
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "pctl",
    lower = 0.05,
    upper = 0.95
  )

  expect_s3_class(truncated, c("ps_trunc_matrix", "ps_trunc", "matrix"))
  expect_equal(dim(truncated), dim(ps_matrix))

  # Check metadata
  meta <- ps_trunc_meta(truncated)
  expect_equal(meta$method, "pctl")
  expect_equal(meta$lower_pctl, 0.05)
  expect_equal(meta$upper_pctl, 0.95)
  expect_true(!is.na(meta$lower_bound))
  expect_true(!is.na(meta$upper_bound))

  # Check that rows still sum to 1
  row_sums <- rowSums(truncated)
  expect_equal(row_sums, rep(1, n), tolerance = 1e-10)

  # The thresholds should be based on the overall distribution
  all_ps_vals <- as.vector(ps_matrix)
  expected_lower <- quantile(all_ps_vals, 0.05)
  expected_upper <- quantile(all_ps_vals, 0.95)
  expect_equal(meta$lower_bound, expected_lower)
  expect_equal(meta$upper_bound, expected_upper)
})

test_that("categorical truncation preserves relative proportions", {
  # When truncating, the relative proportions of non-truncated values should be preserved
  set.seed(456)
  n <- 50

  ps_matrix <- matrix(
    c(
      0.001,
      0.599,
      0.4, # First value needs truncation
      0.3,
      0.001,
      0.699, # Second value needs truncation
      0.4,
      0.4,
      0.2
    ), # No truncation needed
    ncol = 3,
    byrow = TRUE
  )
  ps_matrix <- ps_matrix[rep(1:3, length.out = n), ]

  exposure <- factor(rep(c("A", "B", "C"), length.out = n))
  colnames(ps_matrix) <- levels(exposure)

  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.01
  )

  # For first row type: after truncating 0.001 to 0.01,
  # the ratio of second to third should be preserved
  row1_idx <- seq(1, n, by = 3)
  if (length(row1_idx) > 0) {
    original_ratio <- ps_matrix[row1_idx[1], 2] / ps_matrix[row1_idx[1], 3]
    new_ratio <- truncated[row1_idx[1], 2] / truncated[row1_idx[1], 3]
    expect_equal(unname(new_ratio), unname(original_ratio), tolerance = 1e-10)
  }
})

test_that("print.ps_trunc_matrix produces expected output", {
  set.seed(456)
  n <- 30
  exposure <- factor(sample(c("Low", "Med", "High"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.05, 0.95), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Test ps method
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.1
  )
  output <- capture.output(print(truncated, n = 3))

  expect_match(
    output[1],
    "<ps_trunc_matrix\\[30 x 3\\]; truncated \\d+ of 30; method=ps\\[0\\.1000,Inf\\]>"
  )
  expect_match(output[2], "High\\s+Low\\s+Med")
  expect_match(output[3], "\\[1,\\]")

  # Test pctl method
  truncated_pctl <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "pctl",
    lower = 0.05,
    upper = 0.95
  )
  output_pctl <- capture.output(print(truncated_pctl))

  expect_match(
    output_pctl[1],
    "<ps_trunc_matrix\\[30 x 3\\]; truncated \\d+ of 30; method=pctl\\[0\\.05,0\\.95\\]>"
  )

  # Test without column names
  ps_matrix_no_names <- ps_matrix
  colnames(ps_matrix_no_names) <- NULL
  suppressWarnings({
    truncated_no_names <- ps_trunc(
      ps_matrix_no_names,
      .exposure = exposure,
      method = "ps",
      lower = 0.1
    )
  })
  output_no_names <- capture.output(print(truncated_no_names))

  # Should not have column header
  expect_false(any(grepl("High\\s+Low\\s+Med", output_no_names)))
})
