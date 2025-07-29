test_that("wt_entropy works for binary cases", {
  expect_message(
    weights <- wt_entropy(c(.1, .3, .4, .3), .exposure = c(0, 0, 1, 0)),
    "Treating `.exposure` as binary"
  )

  expect_silent(
    weights2 <- wt_entropy(
      c(.1, .3, .4, .3),
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    )
  )

  expect_identical(weights, weights2)

  expect_message(
    weights3 <- wt_entropy(
      c(.1, .3, .4, .3),
      .exposure = as.logical(c(0, 0, 1, 0))
    ),
    "Treating `.exposure` as binary"
  )

  expect_identical(weights, weights3)

  expect_silent(
    weights4 <- wt_entropy(
      c(.1, .3, .4, .3),
      .exposure = c(2, 2, 1, 2),
      exposure_type = "binary",
      .untreated = 2
    )
  )

  expect_identical(weights, weights4)

  expect_error(
    wt_entropy(
      c(-.1, .3, .4, 3.3),
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    ),
    class = "propensity_range_error"
  )

  .exposure <- factor(
    c("untreated", "untreated", "treated", "untreated"),
    levels = c("untreated", "treated")
  )

  expect_message(
    weights5 <- wt_entropy(
      c(.1, .3, .4, .3),
      exposure_type = "binary",
      .exposure = .exposure
    ),
    "Setting treatment to `treated`"
  )

  expect_identical(weights, weights5)
})

test_that("entropy tilting function properties", {
  # Test symmetry: h(e) = h(1-e)
  ps1 <- c(0.2, 0.3, 0.4)
  ps2 <- 1 - ps1

  # Calculate tilting functions
  h1 <- -ps1 * log(ps1) - (1 - ps1) * log(1 - ps1)
  h2 <- -ps2 * log(ps2) - (1 - ps2) * log(1 - ps2)

  expect_equal(h1, h2, tolerance = 1e-10)

  # Test maximum at 0.5
  ps_seq <- seq(0.01, 0.99, by = 0.01)
  h_vals <- -ps_seq * log(ps_seq) - (1 - ps_seq) * log(1 - ps_seq)
  max_idx <- which.max(h_vals)
  expect_equal(ps_seq[max_idx], 0.5, tolerance = 0.01)

  # Test bounds
  # Maximum entropy is log(2) ≈ 0.693
  expect_true(all(h_vals <= log(2) + 1e-10))
  expect_true(all(h_vals >= 0))
})

test_that("entropy weights have expected properties", {
  # Generate random propensity scores
  set.seed(123)
  ps <- runif(100, 0.1, 0.9)
  treatment <- rbinom(100, 1, ps)

  weights <- wt_entropy(ps, .exposure = treatment, exposure_type = "binary")

  # Entropy weights should be positive and finite
  expect_true(all(weights > 0))
  expect_true(all(is.finite(weights)))

  # Weights at e=0.5 should be around log(2)/0.5 ≈ 1.386
  ps_near_half <- abs(ps - 0.5) < 0.01
  if (any(ps_near_half)) {
    expect_true(all(abs(weights[ps_near_half] - 1.386) < 0.1))
  }
})

test_that("entropy weights handle extreme propensity scores", {
  # Near 0 and 1 propensity scores
  ps_extreme <- c(0.001, 0.01, 0.99, 0.999)
  treatment_extreme <- c(0, 0, 1, 1)

  expect_silent(
    weights_extreme <- wt_entropy(
      ps_extreme,
      .exposure = treatment_extreme,
      exposure_type = "binary"
    )
  )

  expect_true(all(is.finite(weights_extreme)))
  expect_true(all(weights_extreme > 0))

  # Extreme weights can be large but should be finite
  expect_true(max(weights_extreme) < 10)
})

test_that("wt_entropy works with ps_trim objects", {
  ps <- c(.1, .3, .4, .3)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.15, upper = 0.85)

  expect_warning(
    weights <- wt_entropy(
      ps_trimmed,
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    ),
    "It appears you trimmed your propensity score but did not refit the model"
  )

  expect_s3_class(weights, "psw")
  expect_equal(estimand(weights), "entropy; trimmed")
  expect_true(attr(weights, "trimmed"))
})

test_that("wt_entropy works with ps_trunc objects", {
  ps <- c(.1, .3, .4, .3)
  ps_truncated <- ps_trunc(ps, lower = 0.15, upper = 0.85)

  weights <- wt_entropy(
    ps_truncated,
    .exposure = c(0, 0, 1, 0),
    exposure_type = "binary"
  )

  expect_s3_class(weights, "psw")
  expect_equal(estimand(weights), "entropy; truncated")
  expect_true(attr(weights, "truncated"))
})

test_that("entropy weights error on unsupported exposure types", {
  expect_error(
    wt_entropy(
      c(.1, .3, .4, .3),
      .exposure = c(1, 2, 3, 4),
      exposure_type = "categorical"
    ),
    class = "propensity_wt_not_supported_error"
  )

  expect_error(
    wt_entropy(
      rnorm(10),
      .exposure = rnorm(10),
      exposure_type = "continuous"
    ),
    class = "propensity_wt_not_supported_error"
  )
})

# Comparison with PSWeight package - weights
test_that("entropy weights match PSweight's raw weights", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Use a simple example where we can verify the calculations
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  ps_true <- plogis(0.5 * x)
  trt <- rbinom(n, 1, ps_true)
  
  # Use the true propensity scores for both implementations
  # This ensures we're comparing the weight calculation, not PS estimation
  
  # Our implementation
  our_weights <- wt_entropy(ps_true, .exposure = trt, exposure_type = "binary")
  
  # PSweight's implementation using SumStat
  # Create a data frame with the required structure
  test_data <- data.frame(
    trt = trt,
    ps = ps_true,
    x = x
  )
  
  # SumStat with provided propensity scores
  ps_sumstat <- PSweight::SumStat(
    ps.estimate = ps_true,
    zname = "trt",
    xname = "x",
    data = test_data,
    weight = "entropy"
  )
  
  # Extract PSweight's raw weights (before normalization)
  # PSweight stores weights in ps.weights$entropy, but these are normalized
  # We need to un-normalize them to compare
  psw_weights_norm <- ps_sumstat$ps.weights$entropy
  
  # Un-normalize by multiplying by the sum of our raw weights in each group
  # PSweight normalizes so weights sum to 1 within each treatment group
  our_sum1 <- sum(our_weights[trt == 1])
  our_sum0 <- sum(our_weights[trt == 0])
  psw_weights_raw <- numeric(n)
  psw_weights_raw[trt == 1] <- psw_weights_norm[trt == 1] * our_sum1
  psw_weights_raw[trt == 0] <- psw_weights_norm[trt == 0] * our_sum0
  
  # Compare raw weights
  expect_equal(as.numeric(our_weights), psw_weights_raw, tolerance = 1e-10)
})

# Comparison with PSWeight package - estimates
test_that("entropy weights give same treatment effect estimates as PSweight", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Use PSweight's example data
  data("psdata_cl", package = "PSweight")
  
  # Calculate treatment effect using PSweight
  ps_result <- PSweight::PSweight(
    ps.formula = trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6,
    yname = "Y",
    data = psdata_cl,
    weight = "entropy"
  )
  
  psweight_ate <- unname(ps_result$muhat[2] - ps_result$muhat[1])
  
  # Calculate using our implementation
  ps_fit <- glm(trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6, 
                data = psdata_cl, family = binomial)
  ps_scores <- fitted(ps_fit)
  
  our_weights <- wt_entropy(ps_scores, .exposure = psdata_cl$trt, 
                           exposure_type = "binary")
  
  # Calculate weighted means
  mu1 <- weighted.mean(psdata_cl$Y[psdata_cl$trt == 1], 
                      our_weights[psdata_cl$trt == 1])
  mu0 <- weighted.mean(psdata_cl$Y[psdata_cl$trt == 0], 
                      our_weights[psdata_cl$trt == 0])
  our_ate <- mu1 - mu0
  
  # Compare estimates - they should be very close
  expect_equal(our_ate, psweight_ate, tolerance = 1e-6)
})

test_that("entropy weighted estimates are reasonable", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Simulate data with known treatment effect
  set.seed(456)
  n <- 500
  x <- rnorm(n)
  ps <- plogis(0.5 * x)
  treatment <- rbinom(n, 1, ps)
  # True treatment effect = 2
  outcome <- 1 + 2 * treatment + 0.5 * x + rnorm(n)

  # Calculate weights
  weights <- wt_entropy(ps, .exposure = treatment, exposure_type = "binary")

  # Weighted means
  mu1 <- weighted.mean(outcome[treatment == 1], weights[treatment == 1])
  mu0 <- weighted.mean(outcome[treatment == 0], weights[treatment == 0])
  ate_est <- mu1 - mu0

  # Should be close to true value of 2
  expect_equal(ate_est, 2, tolerance = 0.5)
})


test_that("entropy weight calculation matches manual calculation", {
  # Test specific values
  ps <- c(0.2, 0.5, 0.8)
  treatment <- c(0, 1, 1)

  # Manual calculation
  ps_clipped <- pmax(pmin(ps, 1 - 1e-8), 1e-8)
  h_e <- -ps_clipped * log(ps_clipped) - (1 - ps_clipped) * log(1 - ps_clipped)

  expected <- numeric(3)
  expected[1] <- h_e[1] / (1 - ps[1]) # Control unit
  expected[2] <- h_e[2] / ps[2] # Treated unit
  expected[3] <- h_e[3] / ps[3] # Treated unit

  weights <- wt_entropy(ps, .exposure = treatment, exposure_type = "binary")

  expect_equal(as.numeric(weights), expected, tolerance = 1e-10)
})
