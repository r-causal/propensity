test_that("ps_trim works with matrix propensity scores for symmetric trimming", {
  # Create test data
  set.seed(123)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix) # Normalize to sum to 1
  colnames(ps_matrix) <- levels(exposure)

  # Test basic symmetric trimming
  trimmed <- ps_trim(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.1
  )

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
  expect_false(anyNA(trimmed[meta$keep_idx, ]))
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
  sum_inv_ps <- rowSums(1 / ps_matrix)
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
  ps_matrix[1:10, ] <- c(0.5, 0.3, 0.2) # Group A
  ps_matrix[11:20, ] <- c(0.3, 0.5, 0.2) # Group B
  ps_matrix[21:30, ] <- c(0.05, 0.05, 0.9) # Group C - will be trimmed
  colnames(ps_matrix) <- c("A", "B", "C")

  # Try trimming with high threshold
  trimmed <- expect_propensity_warning(
    ps_trim(
      ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.1
    )
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
  colnames(ps_matrix) <- levels(exposure)

  # delta >= 1/3 should error
  expect_propensity_error(
    ps_trim(
      ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.35
    )
  )
})

test_that("ps_trim errors for unsupported methods with categorical", {
  n <- 20
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # These methods should error with rlang argument matching error
  expect_propensity_error(
    ps_trim(ps_matrix, .exposure = exposure, method = "adaptive")
  )

  expect_propensity_error(
    ps_trim(ps_matrix, .exposure = exposure, method = "pctl")
  )
})

test_that("ps_trim requires exposure for categorical", {
  n <- 20
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c("A", "B", "C")

  expect_propensity_error(
    ps_trim(ps_matrix, method = "ps")
  )
})

test_that("is_ps_trimmed works for matrix objects", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  expect_false(is_ps_trimmed(ps_matrix))

  trimmed <- ps_trim(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.2
  )
  expect_true(is_ps_trimmed(trimmed))
})

test_that("is_unit_trimmed works for matrix objects", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  trimmed <- expect_propensity_warning(
    ps_trim(
      ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.25
    )
  )
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

  trimmed <- expect_no_error(
    ps_trim(
      ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.1
    )
  )

  expect_s3_class(trimmed, "ps_trim_matrix")
})

test_that("ps_trim warns when no column names provided", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  # No column names

  expect_propensity_warning(
    ps_trim(
      ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.1
    )
  )
})

test_that("ps_trim.ps_trim warns about already trimmed scores", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  trimmed_once <- ps_trim(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.1
  )

  trimmed_twice <- expect_propensity_warning(
    ps_trim(
      trimmed_once,
      .exposure = exposure,
      method = "ps",
      lower = 0.2
    )
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
  psdata$trt_3cat <- cut(
    psdata$cov1,
    breaks = 3,
    labels = c("Low", "Med", "High")
  )

  # Estimate propensity scores using multinomial regression
  ps_formula <- trt_3cat ~ cov1 + cov2 + cov3 + cov4

  # PSweight approach
  expect_warning(
    psw_result <- PSweight::PStrim(
      data = psdata,
      ps.formula = ps_formula,
      zname = "trt_3cat",
      delta = 0.1
    ),
    "One or more groups removed after trimming"
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

  trt <- factor(apply(trt_probs, 1, function(p) sample.int(4, 1, prob = p)))

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
    c(
      rep(c(0.8, 0.15, 0.05), 20), # First 20 obs
      rep(c(0.05, 0.8, 0.15), 15), # Next 15 obs
      rep(c(0.15, 0.05, 0.8), 15)
    ), # Last 15 obs
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- c("A", "B", "C")

  trt <- factor(c(rep("A", 20), rep("B", 15), rep("C", 15)))

  test_data <- data.frame(
    trt = trt,
    dummy = rnorm(n) # PSweight needs at least one covariate
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

  our_trim <- expect_propensity_warning(
    ps_trim(
      ps = ps_matrix,
      .exposure = trt,
      method = "ps",
      lower = 0.06
    )
  )

  # Both should return original data when group is removed
  expect_equal(nrow(psw_trim$data), n)
  expect_equal(length(ps_trim_meta(our_trim)$keep_idx), n)
})

test_that("ps_refit works with categorical propensity score trimming", {
  skip_if_not_installed("nnet")

  # Create test data
  set.seed(234)
  n <- 150
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # Create 3-category exposure based on covariates
  trt_probs <- cbind(
    plogis(0.5 + 0.3 * x1),
    plogis(0.5 - 0.3 * x1 + 0.3 * x2),
    plogis(0.5 + 0.3 * x2)
  )
  trt_probs <- trt_probs / rowSums(trt_probs)

  trt <- factor(apply(
    trt_probs,
    1,
    function(p) sample(c("A", "B", "C"), 1, prob = p)
  ))

  test_data <- data.frame(
    trt = trt,
    x1 = x1,
    x2 = x2
  )

  # Fit multinomial model
  suppressWarnings({
    model <- nnet::multinom(trt ~ x1 + x2, data = test_data, trace = FALSE)
  })

  # Get propensity scores
  ps_matrix <- predict(model, type = "probs")
  if (!is.matrix(ps_matrix)) {
    ps_matrix <- matrix(ps_matrix, nrow = 1)
  }

  # Apply trimming
  trimmed_ps <- ps_trim(
    ps = ps_matrix,
    .exposure = test_data$trt,
    method = "ps",
    lower = 0.15
  )

  refitted_ps <- ps_refit(trimmed_ps, model, .data = test_data)

  # Check properties
  expect_s3_class(refitted_ps, c("ps_trim_matrix", "ps_trim", "matrix"))
  expect_true(is_refit(refitted_ps))
  expect_equal(dim(refitted_ps), dim(ps_matrix))

  # Check metadata is preserved
  orig_meta <- ps_trim_meta(trimmed_ps)
  refit_meta <- ps_trim_meta(refitted_ps)
  expect_equal(refit_meta$keep_idx, orig_meta$keep_idx)
  expect_equal(refit_meta$trimmed_idx, orig_meta$trimmed_idx)
  expect_true(refit_meta$refit)

  # Trimmed observations should still be NA
  expect_true(all(is.na(refitted_ps[refit_meta$trimmed_idx, ])))

  # Kept observations should have valid propensity scores
  kept_ps <- refitted_ps[refit_meta$keep_idx, ]
  expect_false(anyNA(kept_ps))
  expect_true(all(kept_ps >= 0))
  expect_true(all(kept_ps <= 1))

  # Rows should sum to 1
  row_sums <- rowSums(kept_ps)
  expect_equal(row_sums, rep(1, length(row_sums)), tolerance = 1e-10)
})

test_that("ps_refit errors when all observations are trimmed for categorical", {
  set.seed(345)
  n <- 30
  exposure <- factor(rep(c("A", "B", "C"), each = 10))

  # Create propensity scores where all will be trimmed
  ps_matrix <- matrix(
    rep(c(0.05, 0.05, 0.9), n),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- c("A", "B", "C")

  # Apply very strict trimming
  expect_propensity_warning(
    ps_trim(
      ps = ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.3 # This will trim everything
    )
  )

  # If warning occurred, all data was retained
  # So create a scenario where we manually trim everything
  ps_matrix_extreme <- matrix(
    rep(c(0.01, 0.01, 0.98), n),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix_extreme) <- c("A", "B", "C")

  trimmed_extreme <- expect_propensity_warning(
    ps_trim(
      ps = ps_matrix_extreme,
      .exposure = exposure,
      method = "ps",
      lower = 0.02
    )
  )

  # Manually set all to NA to simulate complete trimming
  trimmed_all_na <- trimmed_extreme
  trimmed_all_na[] <- NA_real_
  attr(trimmed_all_na, "ps_trim_meta")$keep_idx <- integer(0)
  attr(trimmed_all_na, "ps_trim_meta")$trimmed_idx <- seq_len(n)

  # Mock model
  model <- list(class = "multinom")
  test_data <- data.frame(x = rnorm(n))

  # Should error
  expect_propensity_error(
    ps_refit(trimmed_all_na, model, .data = test_data)
  )
})

test_that("ps_refit handles optimal trimming for categorical exposures", {
  skip_if_not_installed("nnet")

  # Create test data
  set.seed(456)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  # Create 4-category exposure
  trt <- factor(sample(c("A", "B", "C", "D"), n, replace = TRUE))

  test_data <- data.frame(
    trt = trt,
    x1 = x1,
    x2 = x2
  )

  # Fit multinomial model
  suppressWarnings({
    model <- nnet::multinom(trt ~ x1 + x2, data = test_data, trace = FALSE)
  })

  # Get propensity scores
  ps_matrix <- predict(model, type = "probs")
  if (!is.matrix(ps_matrix)) {
    ps_matrix <- matrix(ps_matrix, nrow = 1)
  }

  # Apply optimal trimming
  trimmed_ps <- ps_trim(
    ps = ps_matrix,
    .exposure = test_data$trt,
    method = "optimal"
  )

  refitted_ps <- ps_refit(trimmed_ps, model, .data = test_data)

  # Check that it worked
  expect_true(is_refit(refitted_ps))

  # Check that optimal trimming metadata is preserved
  meta <- ps_trim_meta(refitted_ps)
  expect_equal(meta$method, "optimal")
  expect_true(!is.null(meta$lambda))
})

test_that("ps_refit preserves column names and order for categorical", {
  skip_if_not_installed("nnet")

  set.seed(567)
  n <- 60
  x <- rnorm(n)

  # Create 3-category exposure
  trt <- factor(sample(c("Low", "Medium", "High"), n, replace = TRUE))

  test_data <- data.frame(
    trt = trt,
    x = x
  )

  # Fit model
  suppressWarnings({
    model <- nnet::multinom(trt ~ x, data = test_data, trace = FALSE)
  })

  # Get propensity scores
  ps_matrix <- predict(model, type = "probs")
  if (!is.matrix(ps_matrix)) {
    ps_matrix <- matrix(ps_matrix, nrow = 1)
  }

  # Ensure we have the right column names
  colnames(ps_matrix) <- levels(trt)

  # Apply trimming
  trimmed_ps <- ps_trim(
    ps = ps_matrix,
    .exposure = test_data$trt,
    method = "ps",
    lower = 0.1
  )

  # Refit
  refitted_ps <- ps_refit(trimmed_ps, model, .data = test_data)

  # Check column names are preserved
  expect_equal(colnames(refitted_ps), colnames(ps_matrix))
  expect_equal(colnames(refitted_ps), levels(trt))
})

test_that("ps_refit handles minimal data for categorical exposures", {
  skip_if_not_installed("nnet")

  # Create minimal data - one observation per category after trimming
  set.seed(678)
  n <- 30

  # Create data that will result in minimal observations after trimming
  trt <- factor(rep(c("A", "B", "C"), each = 10))
  x <- c(rnorm(10, -2), rnorm(10, 0), rnorm(10, 2)) # Separated groups

  test_data <- data.frame(
    trt = trt,
    x = x
  )

  # Fit model
  suppressWarnings({
    model <- nnet::multinom(trt ~ x, data = test_data, trace = FALSE)
  })

  # Get propensity scores
  ps_matrix <- predict(model, type = "probs")
  if (!is.matrix(ps_matrix)) {
    ps_matrix <- matrix(ps_matrix, nrow = 1)
  }

  # Apply aggressive trimming to leave minimal data
  trimmed_ps <- expect_propensity_warning(
    ps_trim(
      ps = ps_matrix,
      .exposure = test_data$trt,
      method = "ps",
      lower = 0.25 # This should trim many observations
    )
  )

  meta <- ps_trim_meta(trimmed_ps)

  # Only proceed if we have at least 3 observations kept (one per category)
  if (length(meta$keep_idx) >= 3) {
    refitted_ps <- ps_refit(trimmed_ps, model, .data = test_data)

    expect_true(is_refit(refitted_ps))
  }
})

test_that("wt_ate warns when using trimmed but not refitted categorical PS", {
  # Create test data
  set.seed(789)
  n <- 60
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.05, 0.95), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Apply trimming
  trimmed_ps <- ps_trim(
    ps = ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.2
  )

  # Using trimmed but not refitted PS should warn
  expect_propensity_warning(wt_ate(trimmed_ps, .exposure = exposure))
})

test_that("print.ps_trim_matrix produces expected output", {
  set.seed(123)
  n <- 20
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.05, 0.95), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Test basic print
  trimmed <- ps_trim(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.2
  )
  # Use smaller n to ensure consistent output
  output <- capture.output(print(trimmed, n = 3))

  expect_match(
    output[1],
    "<ps_trim_matrix\\[20 x 3\\]; trimmed \\d+ of 20; method=ps>"
  )
  expect_match(output[2], "A\\s+B\\s+C")
  expect_match(output[3], "\\[1,\\]")
  expect_true(any(grepl("# \\.\\.\\. with \\d+ more rows", output)))

  # Test without column names
  ps_matrix_no_names <- ps_matrix
  colnames(ps_matrix_no_names) <- NULL
  suppressWarnings({
    trimmed_no_names <- ps_trim(
      ps_matrix_no_names,
      .exposure = exposure,
      method = "ps",
      lower = 0.2
    )
  })
  output_no_names <- capture.output(print(trimmed_no_names))

  # Should not have column header
  expect_false(any(grepl("A\\s+B\\s+C", output_no_names)))
})

test_that("print methods handle large matrices correctly", {
  set.seed(789)
  n <- 100
  exposure <- factor(sample(LETTERS[1:5], n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 5, 0.05, 0.95), nrow = n, ncol = 5)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  trimmed <- ps_trim(ps_matrix, .exposure = exposure, method = "optimal")
  output <- capture.output(print(trimmed, n = 6))

  # Should only show 6 rows plus "..."
  expect_equal(sum(grepl("^\\[\\s*\\d+,\\]", output)), 6)
  expect_true(any(grepl("# \\.\\.\\. with 94 more rows", output)))
})

test_that("print methods respect n parameter", {
  set.seed(999)
  n <- 25
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  trimmed <- ps_trim(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.15
  )

  # Test different n values
  output_3 <- capture.output(print(trimmed, n = 3))
  output_10 <- capture.output(print(trimmed, n = 10))
  output_inf <- capture.output(print(trimmed, n = Inf))

  # Count rows shown (looking for [digit,] pattern)
  n_rows_3 <- sum(grepl("^\\s*\\[\\d+,\\]", output_3))
  n_rows_10 <- sum(grepl("^\\s*\\[\\d+,\\]", output_10))
  n_rows_inf <- sum(grepl("^\\s*\\[\\d+,\\]", output_inf))

  expect_equal(n_rows_3, 3)
  expect_equal(n_rows_10, 10)
  expect_equal(n_rows_inf, 25)

  # Check "more rows" message
  expect_true(any(grepl("# \\.\\.\\. with 22 more rows", output_3)))
  expect_true(any(grepl("# \\.\\.\\. with 15 more rows", output_10)))
  expect_false(any(grepl("# \\.\\.\\. with", output_inf)))
})

# What the matrix print method shows -----------------------------------------

# The trimming record is metadata, and `ps_trim_meta()` is how a caller reads
# it. The print method summarizes that record in the header and then shows the
# scores, and the scores are what it shows: the record is a set of index vectors
# as long as the data, so printed after the matrix it is a second and longer
# object the caller did not ask to see.

trim_print_fixture <- function() {
  set.seed(123)
  n <- 20
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.05, 0.95), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  ps_trim(ps_matrix, .exposure = exposure, method = "ps", lower = 0.2)
}

test_that("print.ps_trim_matrix() shows some rows without the trimming record", {
  testthat::skip("awaiting implementation")

  trimmed <- trim_print_fixture()
  output <- capture.output(print(trimmed, n = 3))

  expect_false(any(grepl('attr\\(,"ps_trim_meta"\\)', output)))

  # Header, column names, the rows asked for, and the count of the rest.
  expect_length(output, 6)
  expect_match(
    output[1],
    "<ps_trim_matrix\\[20 x 3\\]; trimmed \\d+ of 20; method=ps>"
  )
  expect_equal(sum(grepl("^\\s*\\[\\s*\\d+,\\]", output)), 3)
  expect_true(any(grepl("# \\.\\.\\. with 17 more rows", output)))
})

test_that("print.ps_trim_matrix() shows every row without the trimming record", {
  testthat::skip("awaiting implementation")

  trimmed <- trim_print_fixture()
  output <- capture.output(print(trimmed, n = Inf))

  expect_false(any(grepl('attr\\(,"ps_trim_meta"\\)', output)))

  # Header, column names, and every row.
  expect_length(output, 22)
  expect_equal(sum(grepl("^\\s*\\[\\s*\\d+,\\]", output)), 20)
})

# ---- the matrix method names ps_trim, whoever called it --------------------

# `ps_trim.matrix` validates the propensity score matrix from inside a
# dispatched method, so the frame it hands the validator decides what the error
# reports. A frame taken from the method's caller is the user's own function
# when `ps_trim()` is called from one, and nothing at all when it is called at
# the top level.

condition_call_name <- function(expr) {
  cnd <- rlang::catch_cnd(expr, classes = "error")
  if (is.null(cnd) || is.null(conditionCall(cnd))) {
    return(NA_character_)
  }

  paste(deparse(conditionCall(cnd)[[1]]), collapse = " ")
}

categorical_trim_fixture <- function() {
  list(
    exposure = factor(c("a", "b", "c", "a")),
    # the first row sums to 1.4, which the matrix method refuses
    bad_matrix = matrix(
      c(0.7, 0.3, 0.2, 0.5, 0.4, 0.4, 0.3, 0.2, 0.3, 0.3, 0.5, 0.3),
      nrow = 4,
      ncol = 3
    )
  )
}

test_that("a refused propensity score matrix names ps_trim", {
  fixture <- categorical_trim_fixture()

  cnd <- rlang::catch_cnd(
    ps_trim(fixture$bad_matrix, method = "ps", .exposure = fixture$exposure),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_matrix_sum_error")

  expect_identical(
    condition_call_name(
      ps_trim(fixture$bad_matrix, method = "ps", .exposure = fixture$exposure)
    ),
    "ps_trim"
  )
})

test_that("a refused propensity score matrix names ps_trim, not the caller", {
  fixture <- categorical_trim_fixture()

  wrapping_trim_fn <- function(ps, .exposure) {
    ps_trim(ps, method = "ps", .exposure = .exposure)
  }

  expect_identical(
    condition_call_name(
      wrapping_trim_fn(fixture$bad_matrix, fixture$exposure)
    ),
    "ps_trim"
  )
})

test_that("a matrix with the wrong number of columns names ps_trim", {
  fixture <- categorical_trim_fixture()

  wrapping_trim_fn <- function(ps, .exposure) {
    ps_trim(ps, method = "ps", .exposure = .exposure)
  }

  expect_identical(
    condition_call_name(
      wrapping_trim_fn(
        fixture$bad_matrix[, 1:2, drop = FALSE],
        fixture$exposure
      )
    ),
    "ps_trim"
  )
})

valid_trim_matrix_fixture <- function() {
  exposure <- factor(c("a", "b", "c", "a", "b", "c"))
  # every row sums to 1 and no cell falls below 0.1, so the default threshold
  # trims nothing and leaves all three groups in place
  ps_matrix <- rbind(
    c(0.50, 0.30, 0.20),
    c(0.30, 0.50, 0.20),
    c(0.20, 0.30, 0.50),
    c(0.40, 0.40, 0.20),
    c(0.25, 0.35, 0.40),
    c(0.34, 0.33, 0.33)
  )
  colnames(ps_matrix) <- levels(exposure)

  list(exposure = exposure, ps_matrix = ps_matrix)
}

test_that("ps_trim trims a matrix with the method left at its default", {
  fixture <- valid_trim_matrix_fixture()

  # The generic's default is the full six-value vector, which the matrix method
  # narrows to the two methods it supports. Omitting `method` has to resolve to
  # the first of those rather than failing to match the vector it was handed.
  trimmed <- ps_trim(fixture$ps_matrix, .exposure = fixture$exposure)

  expect_s3_class(trimmed, c("ps_trim_matrix", "ps_trim", "matrix"))
  expect_equal(ps_trim_meta(trimmed)$method, "ps")
  expect_equal(dim(trimmed), dim(fixture$ps_matrix))
})

test_that("ps_trim rejects an unsupported matrix method with a package error", {
  fixture <- valid_trim_matrix_fixture()

  cnd <- rlang::catch_cnd(
    ps_trim(
      fixture$ps_matrix,
      .exposure = fixture$exposure,
      method = "adaptive"
    ),
    classes = "error"
  )

  expect_s3_class(cnd, "propensity_method_error")
  # The message has to name the methods the matrix path does support.
  expect_match(conditionMessage(cnd), "optimal", fixed = TRUE)
})

test_that("ps_trim aborts when the categorical threshold reaches 1/k", {
  fixture <- valid_trim_matrix_fixture()

  # A threshold at or above 1/k cannot be met by every column of a row that
  # sums to one, so there is no trimming rule to apply. Returning the untrimmed
  # scores with the rejected threshold recorded in the metadata reports a
  # trimming that did not happen.
  cnd <- rlang::catch_cnd(
    ps_trim(
      fixture$ps_matrix,
      .exposure = fixture$exposure,
      method = "ps",
      lower = 0.35
    ),
    classes = "error"
  )

  expect_s3_class(cnd, "propensity_range_error")
})

test_that("ps_trim defaults the categorical threshold to 0.1", {
  fixture <- valid_trim_matrix_fixture()

  # Trimming and truncation deliberately default to different thresholds: 0.1
  # here against 0.01 in ps_trunc(), because trimming drops units and
  # truncation only bounds them.
  explicit <- ps_trim(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "ps"
  )
  expect_equal(ps_trim_meta(explicit)$delta, 0.1)

  defaulted <- ps_trim(fixture$ps_matrix, .exposure = fixture$exposure)
  expect_equal(ps_trim_meta(defaulted)$delta, 0.1)
})

test_that("ps_trim trims a data frame with the method left at its default", {
  fixture <- valid_trim_matrix_fixture()

  # The data frame method hands `method` to the matrix method untouched, so the
  # generic's default has to survive that hand-off as well.
  trimmed <- ps_trim(
    as.data.frame(fixture$ps_matrix),
    .exposure = fixture$exposure
  )

  expect_s3_class(trimmed, c("ps_trim_matrix", "ps_trim", "matrix"))
  expect_equal(ps_trim_meta(trimmed)$method, "ps")
  expect_equal(ps_trim_meta(trimmed)$delta, 0.1)
})

# Missing values ------------------------------------------------------------

# A row carrying a missing score is not one this package removed, so the record
# says nothing about it and the row comes back as it arrived. Trimming either
# discards a row or leaves it untouched, so unlike truncation it has no reason
# to rewrite the scores that row does hold.

missing_cell_trim_fixture <- function() {
  exposure <- factor(c("a", "a", "a", "b", "b", "b", "c", "c", "c"))
  # Three rows per level, so that dropping row 3 for its low score and row 7
  # for its missing one still leaves every group represented and the
  # group-preservation reset does not fire. Row 3 falls below the 0.1
  # threshold; row 7 has no score for level "a".
  ps_matrix <- rbind(
    c(0.60, 0.20, 0.20),
    c(0.50, 0.30, 0.20),
    c(0.05, 0.90, 0.05),
    c(0.20, 0.60, 0.20),
    c(0.30, 0.50, 0.20),
    c(0.40, 0.40, 0.20),
    c(NA, 0.50, 0.50),
    c(0.20, 0.20, 0.60),
    c(0.30, 0.20, 0.50)
  )
  colnames(ps_matrix) <- levels(exposure)

  list(exposure = exposure, ps_matrix = ps_matrix)
}

test_that("ps_trim() does not record a row with a missing score as trimmed", {
  fixture <- missing_cell_trim_fixture()

  # `min(x) > delta` is missing for row 7, which drops it from the retained
  # positions and lands it among the trimmed ones, where it is blanked out and
  # two observed scores are lost with it.
  with_na <- ps_trim(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "ps",
    lower = 0.1
  )
  without_na <- ps_trim(
    fixture$ps_matrix[-7, ],
    .exposure = fixture$exposure[-7],
    method = "ps",
    lower = 0.1
  )
  meta <- ps_trim_meta(with_na)

  expect_equal(meta$trimmed_idx, 3)
  expect_equal(meta$keep_idx, c(1, 2, 4, 5, 6, 8, 9))
  expect_equal(meta$n_obs, 9)
  expect_false(is_unit_trimmed(with_na)[7])

  # The row is neither trimmed nor rewritten, so it comes back as it arrived.
  expect_equal(unclass(with_na)[7, ], fixture$ps_matrix[7, ])

  # Every other row is decided by the same comparison against the same
  # threshold, so the missing row changes nothing about them.
  expect_equal(
    as.numeric(unclass(with_na)[-7, ]),
    as.numeric(unclass(without_na))
  )
})

test_that("ps_trim() takes its optimal threshold from the complete rows", {
  fixture <- missing_cell_trim_fixture()

  # The optimal threshold is a root of a function of the row sums of 1 / e,
  # which is missing for row 7, so the comparison that decides whether to solve
  # for it at all is an error rather than an answer.
  with_na <- ps_trim(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "optimal"
  )
  without_na <- ps_trim(
    fixture$ps_matrix[-7, ],
    .exposure = fixture$exposure[-7],
    method = "optimal"
  )
  meta <- ps_trim_meta(with_na)

  expect_equal(meta$lambda, ps_trim_meta(without_na)$lambda)
  expect_equal(meta$trimmed_idx, 3)
  expect_equal(meta$keep_idx, c(1, 2, 4, 5, 6, 8, 9))
  expect_false(is_unit_trimmed(with_na)[7])
  expect_equal(unclass(with_na)[7, ], fixture$ps_matrix[7, ])
  expect_equal(
    as.numeric(unclass(with_na)[-7, ]),
    as.numeric(unclass(without_na))
  )
})
