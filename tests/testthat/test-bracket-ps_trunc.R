test_that("[.ps_trunc updates indices correctly", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)

  # Create ps_trunc object with known truncated positions
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)
  meta <- ps_trunc_meta(ps_truncated)

  # Record initial indices
  initial_truncated <- meta$truncated_idx

  # Test basic numeric subsetting
  subset1 <- ps_truncated[5:8]
  meta1 <- ps_trunc_meta(subset1)

  # Calculate expected indices after subsetting
  expected_truncated <- match(initial_truncated, 5:8)
  expected_truncated <- expected_truncated[!is.na(expected_truncated)]

  expect_equal(meta1$truncated_idx, expected_truncated)
  expect_length(subset1, 4)

  # Bounds should be preserved
  expect_equal(meta1$lower_bound, meta$lower_bound)
  expect_equal(meta1$upper_bound, meta$upper_bound)
})

test_that("[.ps_trunc handles logical indices", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)
  meta <- ps_trunc_meta(ps_truncated)

  # Logical subsetting
  keep_mask <- as.numeric(ps_truncated) > 0.5
  subset <- ps_truncated[keep_mask]
  meta_subset <- ps_trunc_meta(subset)

  # Calculate expected indices
  kept_positions <- which(keep_mask)
  expected_truncated <- match(meta$truncated_idx, kept_positions)
  expected_truncated <- expected_truncated[!is.na(expected_truncated)]

  expect_equal(meta_subset$truncated_idx, expected_truncated)
})

test_that("[.ps_trunc handles negative indices", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)
  meta <- ps_trunc_meta(ps_truncated)

  # Negative indexing
  subset <- ps_truncated[-(1:5)]
  meta_subset <- ps_trunc_meta(subset)

  # Calculate expected indices (positions 6:10)
  remaining_positions <- 6:10
  expected_truncated <- match(meta$truncated_idx, remaining_positions)
  expected_truncated <- expected_truncated[!is.na(expected_truncated)]

  expect_equal(meta_subset$truncated_idx, expected_truncated)
  expect_length(subset, 5)
})

test_that("[.ps_trunc handles empty subsets", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)

  # Empty subset
  empty <- ps_truncated[integer(0)]
  meta_empty <- ps_trunc_meta(empty)

  expect_length(empty, 0)
  expect_equal(meta_empty$truncated_idx, integer(0))
  expect_s3_class(empty, "ps_trunc")
})

test_that("[.ps_trunc_matrix updates indices correctly", {
  set.seed(123)
  n <- 20
  k <- 3

  # Create matrix with some extreme values
  ps_mat <- matrix(runif(n * k, 0.001, 0.999), ncol = k)
  ps_mat <- ps_mat / rowSums(ps_mat)
  # Add column names to avoid warnings
  colnames(ps_mat) <- LETTERS[1:k]

  # Create exposure
  exposure <- factor(sample(LETTERS[1:k], n, replace = TRUE))

  # Truncate the matrix
  ps_truncated <- ps_trunc(
    ps_mat,
    method = "ps",
    lower = 0.1,
    .exposure = exposure
  )
  meta <- ps_trunc_meta(ps_truncated)

  # Subset rows
  subset_rows <- ps_truncated[5:10, ]
  meta_subset <- ps_trunc_meta(subset_rows)

  # Calculate expected indices
  expected_truncated <- match(meta$truncated_idx, 5:10)
  expected_truncated <- expected_truncated[!is.na(expected_truncated)]

  expect_equal(meta_subset$truncated_idx, expected_truncated)
  expect_equal(nrow(subset_rows), 6)
  expect_equal(ncol(subset_rows), k)
})

test_that("is_unit_truncated works correctly after subsetting", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)

  # Get unit truncated status before subsetting
  unit_truncated_full <- is_unit_truncated(ps_truncated)

  # Subset
  indices <- 4:8
  subset <- ps_truncated[indices]
  unit_truncated_subset <- is_unit_truncated(subset)

  # The unit truncated status for the subset should match the original
  # for the selected indices
  expect_equal(unit_truncated_subset, unit_truncated_full[indices])
  expect_length(unit_truncated_subset, length(indices))
})

test_that("subsetting preserves trunc metadata fields", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)

  # Test percentile method
  ps_pctl <- ps_trunc(ps, method = "pctl", lower = 0.05, upper = 0.95)
  subset_pctl <- ps_pctl[1:5]
  meta_pctl <- ps_trunc_meta(subset_pctl)

  expect_equal(meta_pctl$method, "pctl")
  expect_true("lower_pctl" %in% names(meta_pctl))
  expect_true("upper_pctl" %in% names(meta_pctl))
  expect_equal(meta_pctl$lower_pctl, 0.05)
  expect_equal(meta_pctl$upper_pctl, 0.95)
})

test_that("subsetting maintains bounds after subsetting", {
  set.seed(123)
  ps <- runif(20, 0.1, 0.9)

  # Truncate with specific bounds
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)

  # Check that all values are within bounds
  expect_true(all(as.numeric(ps_truncated) >= 0.2))
  expect_true(all(as.numeric(ps_truncated) <= 0.8))

  # Subset and verify bounds are still respected
  subset <- ps_truncated[5:15]
  expect_true(all(as.numeric(subset) >= 0.2))
  expect_true(all(as.numeric(subset) <= 0.8))

  # Metadata should preserve bound information
  meta_sub <- ps_trunc_meta(subset)
  expect_equal(meta_sub$lower_bound, 0.2)
  expect_equal(meta_sub$upper_bound, 0.8)
})

test_that("subsetting handles complex index patterns for ps_trunc", {
  set.seed(123)
  ps <- runif(20, 0.1, 0.9)
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)
  meta <- ps_trunc_meta(ps_truncated)

  # Non-contiguous indices
  indices <- c(1, 3, 5, 7, 9, 11, 13, 15)
  subset <- ps_truncated[indices]
  meta_subset <- ps_trunc_meta(subset)

  # Calculate expected indices
  expected_truncated <- match(meta$truncated_idx, indices)
  expected_truncated <- expected_truncated[!is.na(expected_truncated)]

  expect_equal(meta_subset$truncated_idx, expected_truncated)

  # Repeated indices
  indices_rep <- c(1, 1, 2, 2, 3, 3)
  subset_rep <- ps_truncated[indices_rep]
  expect_length(subset_rep, length(indices_rep))
})
