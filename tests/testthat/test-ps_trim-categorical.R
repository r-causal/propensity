test_that("ps_trim works with matrix propensity scores for symmetric trimming", {
  # Create test data
  set.seed(123)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)  # Normalize to sum to 1
  colnames(ps_matrix) <- levels(exposure)
  
  # Test basic symmetric trimming
  trimmed <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.1)
  
  expect_s3_class(trimmed, c("ps_trim_matrix", "ps_trim", "matrix"))
  expect_equal(dim(trimmed), dim(ps_matrix))
  
  # Check metadata
  meta <- ps_trim_meta(trimmed)
  expect_equal(meta$method, "ps")
  expect_equal(meta$delta, 0.1)
  expect_true(meta$is_matrix)
  
  # Check that rows with any value < 0.1 are trimmed
  min_vals <- apply(ps_matrix, 1, min)
  expected_trimmed <- which(min_vals <= 0.1)
  expect_equal(sort(meta$trimmed_idx), sort(expected_trimmed))
  
  # Check that trimmed rows are NA
  expect_true(all(is.na(trimmed[meta$trimmed_idx, ])))
  expect_false(any(is.na(trimmed[meta$keep_idx, ])))
})

test_that("ps_trim works with data.frame propensity scores", {
  set.seed(456)
  n <- 50
  exposure <- factor(sample(c("Low", "Med", "High"), n, replace = TRUE))
  
  ps_df <- data.frame(
    Low = runif(n, 0.1, 0.6),
    Med = runif(n, 0.2, 0.5),
    High = runif(n, 0.1, 0.4)
  )
  # Normalize rows
  ps_df <- ps_df / rowSums(ps_df)
  
  trimmed <- ps_trim(ps_df, .exposure = exposure, method = "ps", lower = 0.2)
  
  expect_s3_class(trimmed, c("ps_trim_matrix", "ps_trim", "matrix"))
  expect_equal(nrow(trimmed), nrow(ps_df))
  expect_equal(ncol(trimmed), ncol(ps_df))
})

test_that("optimal trimming works for categorical exposures", {
  set.seed(789)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C", "D"), n, replace = TRUE))
  
  # Create propensity scores with some extreme values
  ps_matrix <- matrix(runif(n * 4, 0.05, 0.8), nrow = n, ncol = 4)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)
  
  # Apply optimal trimming
  trimmed <- ps_trim(ps_matrix, .exposure = exposure, method = "optimal")
  
  expect_s3_class(trimmed, c("ps_trim_matrix", "ps_trim", "matrix"))
  
  # Check metadata
  meta <- ps_trim_meta(trimmed)
  expect_equal(meta$method, "optimal")
  expect_true(!is.null(meta$lambda))
  
  # Verify that high sum of inverse PS are trimmed
  sum_inv_ps <- rowSums(1/ps_matrix)
  if (!is.null(meta$lambda)) {
    expect_true(all(sum_inv_ps[meta$keep_idx] <= meta$lambda))
  }
})

test_that("ps_trim preserves all treatment groups", {
  set.seed(321)
  n <- 30
  exposure <- factor(c(rep("A", 10), rep("B", 10), rep("C", 10)))
  
  # Create propensity scores where group C has very low values
  ps_matrix <- matrix(0, nrow = n, ncol = 3)
  ps_matrix[1:10, ] <- c(0.5, 0.3, 0.2)  # Group A
  ps_matrix[11:20, ] <- c(0.3, 0.5, 0.2)  # Group B
  ps_matrix[21:30, ] <- c(0.05, 0.05, 0.9)  # Group C - will be trimmed
  colnames(ps_matrix) <- c("A", "B", "C")
  
  # Try trimming with high threshold
  expect_warning(
    trimmed <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.1),
    class = "propensity_no_data_warning"
  )
  
  # Should return original data
  meta <- ps_trim_meta(trimmed)
  expect_equal(length(meta$trimmed_idx), 0)
  expect_equal(length(meta$keep_idx), n)
})

test_that("ps_trim validates delta < 1/k", {
  n <- 50
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  
  # delta >= 1/3 should trigger warning
  expect_warning(
    trimmed <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.35),
    class = "propensity_range_warning"
  )
  
  # Should return original data
  meta <- ps_trim_meta(trimmed)
  expect_equal(length(meta$keep_idx), n)
})

test_that("ps_trim errors for unsupported methods with categorical", {
  n <- 20
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  
  # These methods should error with rlang argument matching error
  expect_error(
    ps_trim(ps_matrix, .exposure = exposure, method = "adaptive"),
    "must be one of"
  )
  
  expect_error(
    ps_trim(ps_matrix, .exposure = exposure, method = "pctl"),
    "must be one of"
  )
})

test_that("ps_trim requires exposure for categorical", {
  n <- 20
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  
  expect_error(
    ps_trim(ps_matrix, method = "ps"),
    class = "propensity_missing_arg_error"
  )
})

test_that("is_ps_trimmed works for matrix objects", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  
  expect_false(is_ps_trimmed(ps_matrix))
  
  trimmed <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.2)
  expect_true(is_ps_trimmed(trimmed))
})

test_that("is_unit_trimmed works for matrix objects", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  
  trimmed <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.25)
  unit_trimmed <- is_unit_trimmed(trimmed)
  
  expect_equal(length(unit_trimmed), n)
  expect_type(unit_trimmed, "logical")
  
  meta <- ps_trim_meta(trimmed)
  expect_equal(which(unit_trimmed), meta$trimmed_idx)
})

test_that("ps_trim handles parsnip-style column names", {
  n <- 40
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c(".pred_A", ".pred_B", ".pred_C")
  
  expect_no_error(
    trimmed <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.1)
  )
  
  expect_s3_class(trimmed, "ps_trim_matrix")
})

test_that("ps_trim warns when no column names provided", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  # No column names
  
  expect_warning(
    trimmed <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.1),
    class = "propensity_matrix_no_names_warning"
  )
})

test_that("ps_trim.ps_trim warns about already trimmed scores", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  
  trimmed_once <- ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.1)
  
  expect_warning(
    trimmed_twice <- ps_trim(trimmed_once, .exposure = exposure, method = "ps", lower = 0.2),
    class = "propensity_already_modified_warning"
  )
  
  # Should return original trimmed object
  expect_identical(trimmed_twice, trimmed_once)
})

# PSweight comparison tests
test_that("ps_trim symmetric trimming matches PSweight for categorical exposures", {
  skip_if_not_installed("PSweight")
  skip_on_cran()
  
  # Load PSweight data
  data(psdata, package = "PSweight")
  
  # Create a 3-category exposure from existing variables
  set.seed(123)
  psdata$trt_3cat <- cut(psdata$cov1, breaks = 3, labels = c("Low", "Med", "High"))
  
  # Estimate propensity scores using multinomial regression
  ps_formula <- trt_3cat ~ cov1 + cov2 + cov3 + cov4
  
  # PSweight approach
  psw_result <- PSweight::PStrim(
    data = psdata,
    ps.formula = ps_formula,
    zname = "trt_3cat",
    delta = 0.1
  )
  
  # Extract propensity scores from PSweight
  ps_matrix_psw <- psw_result$ps.estimate
  
  # Check if PSweight returned a matrix
  if (is.null(ps_matrix_psw) || !is.matrix(ps_matrix_psw)) {
    skip("PSweight didn't return a valid propensity score matrix")
  }
  
  # Our approach
  our_result <- ps_trim(
    ps = ps_matrix_psw,
    .exposure = psdata$trt_3cat,
    method = "ps",
    lower = 0.1
  )
  
  # Compare results
  meta <- ps_trim_meta(our_result)
  
  # PSweight returns the trimmed dataset
  # Count how many were retained
  n_retained_psw <- nrow(psw_result$data)
  n_retained_our <- length(meta$keep_idx)
  
  # PSweight warning indicates no trimming was applied due to group removal
  # In this case, both should retain all observations
  if (n_retained_psw == nrow(psdata)) {
    expect_equal(n_retained_our, nrow(psdata))
    expect_equal(length(meta$trimmed_idx), 0)
  } else {
    # Both should retain the same number of observations
    expect_equal(n_retained_our, n_retained_psw)
    
    # Check that both identified similar proportions to trim
    prop_trimmed_psw <- (nrow(psdata) - n_retained_psw) / nrow(psdata)
    prop_trimmed_our <- length(meta$trimmed_idx) / nrow(psdata)
    expect_equal(prop_trimmed_our, prop_trimmed_psw, tolerance = 0.01)
  }
})

test_that("ps_trim optimal trimming matches PSweight for multi-category", {
  skip_if_not_installed("PSweight")
  skip_on_cran()
  
  # Create test data with 4 categories
  set.seed(456)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  
  # Create 4-category treatment based on covariates
  # Use plogis to ensure positive probabilities
  trt_probs <- cbind(
    plogis(0.5 + 0.3 * x1),
    plogis(0.5 - 0.3 * x1 + 0.3 * x2),
    plogis(0.5 + 0.3 * x2),
    plogis(0.5 - 0.3 * x2)
  )
  trt_probs <- trt_probs / rowSums(trt_probs)
  
  trt <- factor(apply(trt_probs, 1, function(p) sample(1:4, 1, prob = p)))
  
  test_data <- data.frame(
    trt = trt,
    x1 = x1,
    x2 = x2
  )
  
  # PSweight optimal trimming
  psw_optimal <- PSweight::PStrim(
    data = test_data,
    ps.formula = trt ~ x1 + x2,
    zname = "trt",
    optimal = TRUE
  )
  
  # Our optimal trimming
  ps_matrix <- psw_optimal$ps.estimate
  
  # Check if PSweight returned valid propensity scores
  if (is.null(ps_matrix)) {
    skip("PSweight didn't return propensity scores for optimal trimming")
  }
  
  # Ensure ps_matrix is actually a matrix
  if (!is.matrix(ps_matrix)) {
    ps_matrix <- as.matrix(ps_matrix)
  }
  
  our_optimal <- ps_trim(
    ps = ps_matrix,
    .exposure = test_data$trt,
    method = "optimal"
  )
  
  # Compare lambda values if both found a solution
  meta <- ps_trim_meta(our_optimal)
  if (!is.null(meta$lambda) && !is.null(psw_optimal$lambda)) {
    expect_equal(meta$lambda, psw_optimal$lambda, tolerance = 1e-4)
  }
  
  # Compare which observations were kept
  psw_keep_idx <- seq_len(n) %in% seq_len(nrow(psw_optimal$data))
  our_keep_idx <- seq_len(n) %in% meta$keep_idx
  
  expect_equal(our_keep_idx, psw_keep_idx)
})

test_that("ps_trim handles edge cases consistently with PSweight", {
  skip_if_not_installed("PSweight")
  skip_on_cran()
  
  # Test with extreme propensity scores
  set.seed(789)
  n <- 50
  
  # Create data where one category has very low propensity scores
  ps_matrix <- matrix(
    c(rep(c(0.8, 0.15, 0.05), 20),  # First 20 obs
      rep(c(0.05, 0.8, 0.15), 15),   # Next 15 obs
      rep(c(0.15, 0.05, 0.8), 15)),  # Last 15 obs
    ncol = 3,
    byrow = TRUE
  )
  
  trt <- factor(c(rep("A", 20), rep("B", 15), rep("C", 15)))
  
  test_data <- data.frame(
    trt = trt,
    dummy = rnorm(n)  # PSweight needs at least one covariate
  )
  
  # With delta=0.06, groups with min propensity of 0.05 would be removed
  # causing both implementations to warn and return original data
  suppressWarnings({
    psw_trim <- PSweight::PStrim(
      data = test_data,
      ps.estimate = ps_matrix,
      zname = "trt",
      delta = 0.06
    )
  })
  
  expect_warning(
    our_trim <- ps_trim(
      ps = ps_matrix,
      .exposure = trt,
      method = "ps",
      lower = 0.06
    ),
    class = "propensity_no_data_warning"
  )
  
  # Both should return original data when group is removed
  expect_equal(nrow(psw_trim$data), n)
  expect_equal(length(ps_trim_meta(our_trim)$keep_idx), n)
})