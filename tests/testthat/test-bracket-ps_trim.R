test_that("[.ps_trim updates indices correctly", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)

  # Create ps_trim object with known trimmed positions
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.4, upper = 0.6)
  meta <- ps_trim_meta(ps_trimmed)

  # Record initial indices
  initial_trimmed <- meta$trimmed_idx
  initial_keep <- meta$keep_idx

  # Test basic numeric subsetting
  subset1 <- ps_trimmed[5:8]
  meta1 <- ps_trim_meta(subset1)

  # Calculate expected indices after subsetting
  expected_trimmed <- match(initial_trimmed, 5:8)
  expected_trimmed <- expected_trimmed[!is.na(expected_trimmed)]

  expected_keep <- match(initial_keep, 5:8)
  expected_keep <- expected_keep[!is.na(expected_keep)]

  expect_equal(meta1$trimmed_idx, expected_trimmed)
  expect_equal(meta1$keep_idx, expected_keep)
  expect_length(subset1, 4)
})

test_that("[.ps_trim handles logical indices", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.4, upper = 0.6)
  meta <- ps_trim_meta(ps_trimmed)

  # Logical subsetting
  keep_mask <- as.numeric(ps_trimmed) > 0.5 | is.na(ps_trimmed)
  subset <- ps_trimmed[keep_mask]
  meta_subset <- ps_trim_meta(subset)

  # Calculate expected indices
  kept_positions <- which(keep_mask)
  expected_trimmed <- match(meta$trimmed_idx, kept_positions)
  expected_trimmed <- expected_trimmed[!is.na(expected_trimmed)]

  expect_equal(meta_subset$trimmed_idx, expected_trimmed)
})

test_that("[.ps_trim handles negative indices", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.4, upper = 0.6)
  meta <- ps_trim_meta(ps_trimmed)

  # Negative indexing
  subset <- ps_trimmed[-(1:5)]
  meta_subset <- ps_trim_meta(subset)

  # Calculate expected indices (positions 6:10)
  remaining_positions <- 6:10
  expected_trimmed <- match(meta$trimmed_idx, remaining_positions)
  expected_trimmed <- expected_trimmed[!is.na(expected_trimmed)]

  expect_equal(meta_subset$trimmed_idx, expected_trimmed)
  expect_length(subset, 5)
})

test_that("[.ps_trim handles empty subsets", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.4, upper = 0.6)

  # Empty subset
  empty <- ps_trimmed[integer(0)]
  meta_empty <- ps_trim_meta(empty)

  expect_length(empty, 0)
  expect_equal(meta_empty$trimmed_idx, integer(0))
  expect_equal(meta_empty$keep_idx, integer(0))
  expect_s3_class(empty, "ps_trim")
})

test_that("[.ps_trim_matrix updates indices correctly", {
  set.seed(123)
  n <- 20
  k <- 3

  # Create matrix
  ps_mat <- matrix(runif(n * k), ncol = k)
  ps_mat <- ps_mat / rowSums(ps_mat)
  # Add column names to avoid warnings
  colnames(ps_mat) <- LETTERS[1:k]

  # Create exposure
  exposure <- factor(sample(LETTERS[1:k], n, replace = TRUE))

  # Trim the matrix
  ps_trimmed <- ps_trim(
    ps_mat,
    method = "ps",
    lower = 0.2,
    .exposure = exposure
  )
  meta <- ps_trim_meta(ps_trimmed)

  # Subset rows
  subset_rows <- ps_trimmed[5:10, ]
  meta_subset <- ps_trim_meta(subset_rows)

  # Calculate expected indices
  expected_trimmed <- match(meta$trimmed_idx, 5:10)
  expected_trimmed <- expected_trimmed[!is.na(expected_trimmed)]

  expect_equal(meta_subset$trimmed_idx, expected_trimmed)
  expect_equal(nrow(subset_rows), 6)
  expect_equal(ncol(subset_rows), k)
})

test_that("is_unit_trimmed works correctly after subsetting", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.4, upper = 0.6)

  # Get unit trimmed status before subsetting
  unit_trimmed_full <- is_unit_trimmed(ps_trimmed)

  # Subset
  indices <- 4:8
  subset <- ps_trimmed[indices]
  unit_trimmed_subset <- is_unit_trimmed(subset)

  # The unit trimmed status for the subset should match the original
  # for the selected indices
  expect_equal(unit_trimmed_subset, unit_trimmed_full[indices])
  expect_length(unit_trimmed_subset, length(indices))
})

test_that("subsetting preserves metadata fields", {
  set.seed(123)
  ps <- runif(10, 0.1, 0.9)

  # Test different trimming methods
  ps_adaptive <- ps_trim(ps, method = "adaptive")
  subset_adaptive <- ps_adaptive[1:5]
  meta_adaptive <- ps_trim_meta(subset_adaptive)

  expect_equal(meta_adaptive$method, "adaptive")
  expect_true("cutoff" %in% names(meta_adaptive))

  # Test percentile method
  ps_pctl <- ps_trim(ps, method = "pctl", lower = 0.1, upper = 0.9)
  subset_pctl <- ps_pctl[1:5]
  meta_pctl <- ps_trim_meta(subset_pctl)

  expect_equal(meta_pctl$method, "pctl")
  expect_true("q_lower" %in% names(meta_pctl))
  expect_true("q_upper" %in% names(meta_pctl))
})

test_that("subsetting handles complex index patterns", {
  set.seed(123)
  ps <- runif(20, 0.1, 0.9)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7)
  meta <- ps_trim_meta(ps_trimmed)

  # Non-contiguous indices
  indices <- c(1, 3, 5, 7, 9, 11, 13, 15)
  subset <- ps_trimmed[indices]
  meta_subset <- ps_trim_meta(subset)

  # Calculate expected indices
  expected_trimmed <- match(meta$trimmed_idx, indices)
  expected_trimmed <- expected_trimmed[!is.na(expected_trimmed)]

  expect_equal(meta_subset$trimmed_idx, expected_trimmed)

  # Repeated indices
  indices_rep <- c(1, 1, 2, 2, 3, 3)
  subset_rep <- ps_trimmed[indices_rep]
  expect_length(subset_rep, length(indices_rep))
})
