test_that("weight functions work with trimmed categorical propensity scores", {
  # Create test data
  set.seed(123)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.05, 0.95), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Apply trimming
  trimmed_ps <- ps_trim(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.1
  )

  # Calculate ATE weights - expect warning about refitting
  expect_warning(
    wt_ate_trimmed <- wt_ate(trimmed_ps, .exposure = exposure),
    class = "propensity_no_refit_warning"
  )

  expect_s3_class(wt_ate_trimmed, "psw")
  expect_equal(length(wt_ate_trimmed), n)
  expect_true(is_ps_trimmed(wt_ate_trimmed))
  expect_equal(estimand(wt_ate_trimmed), "ate; trimmed")

  # Check that weights are NA for trimmed units
  meta <- ps_trim_meta(trimmed_ps)
  expect_true(all(is.na(wt_ate_trimmed[meta$trimmed_idx])))
  expect_true(all(!is.na(wt_ate_trimmed[meta$keep_idx])))
})

test_that("weight functions work with truncated categorical propensity scores", {
  # Create test data
  set.seed(456)
  n <- 100
  exposure <- factor(sample(c("Low", "Med", "High"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Apply truncation
  truncated_ps <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )

  # Calculate ATE weights
  wt_ate_truncated <- wt_ate(truncated_ps, .exposure = exposure)

  expect_s3_class(wt_ate_truncated, "psw")
  expect_equal(length(wt_ate_truncated), n)
  expect_true(is_ps_truncated(wt_ate_truncated))
  expect_equal(estimand(wt_ate_truncated), "ate; truncated")

  # No weights should be NA for truncation
  expect_true(all(!is.na(wt_ate_truncated)))
})

test_that("weight functions work with data.frame propensity scores for categorical", {
  # Create test data
  set.seed(789)
  n <- 50
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_df <- data.frame(
    A = runif(n, 0.1, 0.6),
    B = runif(n, 0.2, 0.5),
    C = runif(n, 0.1, 0.4)
  )
  ps_df <- ps_df / rowSums(ps_df)

  # Trim
  trimmed_ps <- ps_trim(
    ps_df,
    .exposure = exposure,
    method = "ps",
    lower = 0.15
  )

  # Calculate weights - expect warning about refitting
  expect_warning(
    wt_ate_trimmed <- wt_ate(trimmed_ps, .exposure = exposure),
    class = "propensity_no_refit_warning"
  )

  expect_s3_class(wt_ate_trimmed, "psw")
  expect_equal(length(wt_ate_trimmed), n)
})

test_that("ATT weights work with categorical trimmed propensity scores", {
  set.seed(111)
  n <- 60
  exposure <- factor(sample(
    c("Control", "Treat1", "Treat2"),
    n,
    replace = TRUE
  ))

  ps_matrix <- matrix(runif(n * 3, 0.1, 0.8), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Trim
  trimmed_ps <- ps_trim(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.2
  )

  # Calculate ATT weights for Treat1 as focal - expect warning about refitting
  expect_warning(
    wt_att_trimmed <- wt_att(
      trimmed_ps,
      .exposure = exposure,
      focal = "Treat1"
    ),
    class = "propensity_no_refit_warning"
  )

  expect_s3_class(wt_att_trimmed, "psw")
  expect_equal(length(wt_att_trimmed), n)
  expect_equal(attr(wt_att_trimmed, "focal_category"), "Treat1")
  expect_equal(estimand(wt_att_trimmed), "att; trimmed")
})

test_that("Entropy weights work with truncated categorical propensity scores", {
  set.seed(222)
  n <- 80
  exposure <- factor(sample(LETTERS[1:4], n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 4, 0.02, 0.8), nrow = n, ncol = 4)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Truncate
  truncated_ps <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )

  # Calculate entropy weights - truncation doesn't need refitting so no warning
  wt_entropy_truncated <- wt_entropy(truncated_ps, .exposure = exposure)

  expect_s3_class(wt_entropy_truncated, "psw")
  expect_equal(length(wt_entropy_truncated), n)
  expect_equal(estimand(wt_entropy_truncated), "entropy; truncated")
})

test_that("weight functions work with refitted trimmed categorical propensity scores", {
  # Create test data
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  # Create exposure based on covariates
  # Generate probabilities for each category
  linear_pred <- cbind(x1, x1 + x2, x2)
  probs <- exp(linear_pred) / rowSums(exp(linear_pred))

  # Sample exposure for each observation
  exposure <- factor(apply(
    probs,
    1,
    function(p) sample(c("A", "B", "C"), 1, prob = p)
  ))

  # Create a simple data frame for modeling
  df <- data.frame(
    exposure = exposure,
    x1 = x1,
    x2 = x2
  )

  # Skip if nnet is not installed (needed for multinom)
  skip_if_not_installed("nnet")

  # Fit multinomial model
  fit <- nnet::multinom(exposure ~ x1 + x2, data = df, trace = FALSE)

  # Get propensity scores
  ps_matrix <- predict(fit, type = "probs")

  # Ensure it's a matrix
  if (!is.matrix(ps_matrix)) {
    ps_matrix <- as.matrix(ps_matrix)
  }

  # Apply trimming
  trimmed_ps <- ps_trim(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.1
  )

  # Refit the model
  refitted_ps <- ps_refit(trimmed_ps, fit, .df = df)

  # Calculate ATE weights - should NOT warn about refitting
  wt_ate_refitted <- wt_ate(refitted_ps, .exposure = exposure)

  expect_s3_class(wt_ate_refitted, "psw")
  expect_equal(length(wt_ate_refitted), n)
  expect_true(is_ps_trimmed(wt_ate_refitted))
  expect_true(is_refit(refitted_ps))
  expect_equal(estimand(wt_ate_refitted), "ate; trimmed")

  # Check that weights are NA for trimmed units
  meta <- ps_trim_meta(refitted_ps)
  expect_true(all(is.na(wt_ate_refitted[meta$trimmed_idx])))
  expect_true(all(!is.na(wt_ate_refitted[meta$keep_idx])))
})

test_that("truncated categorical propensity scores don't trigger refit warning", {
  # Create test data
  set.seed(456)
  n <- 100
  exposure <- factor(sample(c("Low", "Med", "High"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Apply truncation
  truncated_ps <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )

  # Calculate ATE weights - should NOT warn (truncation doesn't need refitting)
  wt_ate_truncated <- wt_ate(truncated_ps, .exposure = exposure)

  expect_s3_class(wt_ate_truncated, "psw")
  expect_equal(length(wt_ate_truncated), n)
  expect_true(is_ps_truncated(wt_ate_truncated))
  expect_equal(estimand(wt_ate_truncated), "ate; truncated")

  # No weights should be NA for truncation
  expect_true(all(!is.na(wt_ate_truncated)))
})
