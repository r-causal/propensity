test_that("ATE works for binary cases", {
  expect_message(
    weights <- wt_ate(c(.1, .3, .4, .3), .exposure = c(0, 0, 1, 0)),
    "Treating `.exposure` as binary"
  )

  expect_silent(
    weights2 <- wt_ate(
      c(.1, .3, .4, .3),
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    )
  )

  expect_identical(weights, weights2)

  expect_message(
    weights3 <- wt_ate(
      c(.1, .3, .4, .3),
      .exposure = as.logical(c(0, 0, 1, 0))
    ),
    "Treating `.exposure` as binary"
  )

  expect_identical(weights, weights3)

  expect_silent(
    weights4 <- wt_ate(
      c(.1, .3, .4, .3),
      .exposure = c(2, 2, 1, 2),
      exposure_type = "binary",
      .untreated = 2
    )
  )

  expect_identical(weights, weights4)

  expect_error(
    wt_ate(
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
    weights5 <- wt_ate(
      c(.1, .3, .4, .3),
      exposure_type = "binary",
      .exposure = .exposure
    ),
    "Setting treatment to `treated`"
  )

  expect_identical(weights, weights5)

  expect_equal(
    weights,
    psw(c(1.11, 1.43, 2.50, 1.43), "ate"),
    tolerance = .01
  )
})

test_that("ATE works for continuous cases", {
  denom_model <- lm(mpg ~ gear + am + carb, data = mtcars)

  # Compute population variances
  un_mean <- mean(mtcars$mpg)
  un_var <- mean((mtcars$mpg - un_mean)^2)
  cond_var <- mean((mtcars$mpg - predict(denom_model))^2)

  # Compute z-scores and densities
  z_num <- (mtcars$mpg - un_mean) / sqrt(un_var)
  z_den <- (mtcars$mpg - predict(denom_model)) / sqrt(cond_var)
  f_num <- dnorm(z_num)
  f_den <- dnorm(z_den)

  # Expected weights
  wts <- 1 / f_den
  stb_wts <- f_num / f_den

  expect_message(
    weights <- wt_ate(
      predict(denom_model),
      .exposure = mtcars$mpg,
      .sigma = influence(denom_model)$sigma,
      exposure_type = "continuous"
    ),
    "Using unstabilized weights for continuous exposures is not recommended."
  )

  expect_equal(weights, psw(wts, "ate"), tolerance = 0.01)
  expect_message(
    stabilized_weights <- wt_ate(
      predict(denom_model),
      .exposure = mtcars$mpg,
      .sigma = influence(denom_model)$sigma,
      stabilize = TRUE,
    ),
    "Treating `.exposure` as continuous"
  )

  expect_equal(
    stabilized_weights,
    psw(stb_wts, "ate", stabilized = TRUE),
    tolerance = 0.01
  )
})

test_that("stabilized weights use P(A=1) and P(A=0) as numerators", {
  ps <- c(0.2, 0.5, 0.8, 0.4)
  A <- c(1, 0, 1, 0)

  p1 <- mean(A)
  p0 <- 1 - p1
  inv_ps <- 1 / ps
  inv_1m <- 1 / (1 - ps)
  expected <- A * inv_ps * p1 + (1 - A) * inv_1m * p0

  got <- ate_binary(ps, A, stabilize = TRUE)

  expect_equal(got, expected)
})


test_that("ATE works for categorical cases", {
  # we don't currently support this!
  expect_error(
    wt_ate(
      c(.1, .3, .4, .3),
      .exposure = c(0, 2, 1, 4),
      exposure_type = "categorical"
    ),
    class = "propensity_wt_not_supported_error"
  )
})

test_that("wt_ate() with ps_trim issues refit warning if not refit, no warning if refit", {
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 1) Trim
  trimmed_ps <- ps_trim(
    ps,
    .exposure = z,
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  # not refit => expect a warning
  expect_warning(
    w_ate_unfit <- wt_ate(
      trimmed_ps,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    ),
    class = "propensity_no_refit_warning"
  )
  expect_s3_class(w_ate_unfit, "psw")
  expect_true(grepl("; trimmed$", estimand(w_ate_unfit)))

  # 2) After refit => no warning
  trimmed_refit <- ps_refit(trimmed_ps, model = fit)
  expect_silent(
    w_ate_fit <- wt_ate(
      trimmed_refit,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    )
  )
  expect_s3_class(w_ate_fit, "psw")
  expect_true(grepl("; trimmed$", estimand(w_ate_fit)))
})

test_that("wt_ate() with ps_trunc adds '; truncated' without refit warning", {
  set.seed(234)
  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.6 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # e.g. bounding at [0.2, 0.8]
  truncated_ps <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)

  # Should produce weighting with no refit warnings
  expect_silent(
    w_ate_trunc <- wt_ate(
      truncated_ps,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    )
  )
  expect_s3_class(w_ate_trunc, "psw")
  # Estimand ends with "; truncated"
  expect_true(grepl("; truncated$", estimand(w_ate_trunc)))
})

test_that("Other estimands (att, atu, etc.) with ps_trim or ps_trunc", {
  set.seed(345)
  n <- 120
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 + 0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # Trim
  trimmed_ps <- ps_trim(ps, .exposure = z, method = "ps")
  # No refit => warning
  expect_warning(
    w_att_trim <- wt_att(
      trimmed_ps,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    ),
    class = "propensity_no_refit_warning"
  )
  # Check estimand
  expect_true(grepl("att; trimmed", estimand(w_att_trim)))

  # Trunc
  truncated_ps <- ps_trunc(ps, method = "pctl", lower = 0.2, upper = 0.8)
  # No warning
  expect_silent(
    w_att_trunc <- wt_att(
      truncated_ps,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    )
  )
  expect_silent(
    w_atu_trunc <- wt_atu(
      truncated_ps,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    )
  )
  expect_true(grepl("atu; truncated", estimand(w_atu_trunc)))
})

test_that("wt_ate() with ps_trunc sets truncated=TRUE in final psw", {
  set.seed(123)
  n <- 8
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  trunc_obj <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  w_ate <- wt_ate(
    trunc_obj,
    .exposure = z,
    exposure_type = "binary",
    .treated = 1
  )

  expect_true(is_ps_truncated(w_ate))
  expect_false(is_ps_trimmed(w_ate))
  expect_match(estimand(w_ate), "; truncated$")
})

test_that("wt_atu.ps_trim triggers refit check, sets 'atu; trimmed'", {
  set.seed(991)
  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(1.6 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 1) Trim the PS
  trimmed_obj <- ps_trim(
    ps,
    .exposure = z,
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  # Not refit => we get a warning
  expect_warning(
    w_atu_unfit <- wt_atu(
      trimmed_obj,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    ),
    class = "propensity_no_refit_warning"
  )
  expect_s3_class(w_atu_unfit, "psw")
  expect_match(estimand(w_atu_unfit), "atu; trimmed")
  expect_true(is_ps_trimmed(w_atu_unfit))
  # ps_trim_meta copied
  expect_identical(
    attr(w_atu_unfit, "ps_trim_meta"),
    attr(trimmed_obj, "ps_trim_meta")
  )

  # 2) Now refit => no warning
  refit_obj <- ps_refit(trimmed_obj, model = fit)
  expect_silent(
    w_atu_fit <- wt_atu(
      refit_obj,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    )
  )
  expect_s3_class(w_atu_fit, "psw")
  expect_match(estimand(w_atu_fit), "atu; trimmed")
  expect_true(is_ps_trimmed(w_atu_fit))
  # confirm ps_trim_meta matches
  expect_identical(
    attr(w_atu_fit, "ps_trim_meta"),
    attr(refit_obj, "ps_trim_meta")
  )
})

test_that("wt_atm.ps_trim triggers refit check, sets 'atm; trimmed'", {
  set.seed(992)
  n <- 50
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  trimmed_obj <- ps_trim(
    ps,
    .exposure = z,
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  # Not refit => warning
  expect_warning(
    w_atm_unfit <- wt_atm(
      trimmed_obj,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    ),
    class = "propensity_no_refit_warning"
  )
  expect_s3_class(w_atm_unfit, "psw")
  expect_match(estimand(w_atm_unfit), "atm; trimmed")
  expect_true(is_ps_trimmed(w_atm_unfit))

  # Refit => no warning
  refit_obj <- ps_refit(trimmed_obj, model = fit)
  expect_silent(
    w_atm_fit <- wt_atm(
      refit_obj,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    )
  )
  expect_s3_class(w_atm_fit, "psw")
  expect_match(estimand(w_atm_fit), "atm; trimmed")
  expect_true(is_ps_trimmed(w_atm_fit))
})

test_that("wt_ato.ps_trim triggers refit check, sets 'ato; trimmed'", {
  set.seed(993)
  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  trimmed_obj <- ps_trim(
    ps,
    .exposure = z,
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )

  # Not refit => warning
  expect_warning(
    w_ato_unfit <- wt_ato(
      trimmed_obj,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    ),
    class = "propensity_no_refit_warning"
  )
  expect_s3_class(w_ato_unfit, "psw")
  expect_match(estimand(w_ato_unfit), "ato; trimmed")
  expect_true(is_ps_trimmed(w_ato_unfit))

  # Refit => no warning
  refit_obj <- ps_refit(trimmed_obj, model = fit)
  expect_silent(
    w_ato_fit <- wt_ato(
      refit_obj,
      .exposure = z,
      exposure_type = "binary",
      .treated = 1
    )
  )
  expect_s3_class(w_ato_fit, "psw")
  expect_match(estimand(w_ato_fit), "ato; trimmed")
  expect_true(is_ps_trimmed(w_ato_fit))
})

# Entropy weight tests
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
    expect_true(all(abs(weights[ps_near_half] - log(2) / 0.5) < 0.1))
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
  # Theoretical upper bound for entropy weights based on extreme propensity scores
  # For extreme values near 0 or 1, weights can grow large but remain finite.
  # Here, we use a calculated bound derived from the entropy function properties.
  max_weight_bound <- log(2) / min(ps_extreme) # Example calculation
  expect_true(max(weights_extreme) < max_weight_bound)
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
  expect_true(is_ps_trimmed(weights))
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
  expect_true(is_ps_truncated(weights))
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

  # Now that continuous is not even an option for entropy,
  # the function will error during auto-detection if given continuous data
  expect_error(
    wt_entropy(
      rnorm(10),
      .exposure = rnorm(10)
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
  ps_fit <- glm(
    trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6,
    data = psdata_cl,
    family = binomial
  )
  ps_scores <- fitted(ps_fit)

  our_weights <- wt_entropy(
    ps_scores,
    .exposure = psdata_cl$trt,
    exposure_type = "binary"
  )

  # Calculate weighted means
  mu1 <- weighted.mean(
    psdata_cl$Y[psdata_cl$trt == 1],
    our_weights[psdata_cl$trt == 1]
  )
  mu0 <- weighted.mean(
    psdata_cl$Y[psdata_cl$trt == 0],
    our_weights[psdata_cl$trt == 0]
  )
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
  h_e <- -ps * log(ps) - (1 - ps) * log(1 - ps)

  expected <- numeric(3)
  expected[1] <- h_e[1] / (1 - ps[1]) # Control unit
  expected[2] <- h_e[2] / ps[2] # Treated unit
  expected[3] <- h_e[3] / ps[3] # Treated unit

  weights <- wt_entropy(ps, .exposure = treatment, exposure_type = "binary")

  expect_equal(as.numeric(weights), expected, tolerance = 1e-10)
})

# Tests for data.frame methods
test_that("wt_ate works with data frames", {
  # Create test data frame
  ps_df <- data.frame(
    control = c(0.9, 0.7, 0.3, 0.1),
    treated = c(0.1, 0.3, 0.7, 0.9)
  )
  exposure <- c(0, 0, 1, 1)

  # Test default behavior (uses second column)
  weights <- wt_ate(ps_df, exposure, exposure_type = "binary")
  expected <- wt_ate(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights, expected)

  # Test explicit column selection by name (quoted)
  weights_quoted <- wt_ate(
    ps_df,
    exposure,
    .propensity_col = "treated",
    exposure_type = "binary"
  )
  expect_equal(weights_quoted, expected)

  # Test unquoted column selection
  weights_unquoted <- wt_ate(
    ps_df,
    exposure,
    .propensity_col = treated,
    exposure_type = "binary"
  )
  expect_equal(weights_unquoted, expected)

  # Test column selection by index
  weights_idx <- wt_ate(
    ps_df,
    exposure,
    .propensity_col = 2,
    exposure_type = "binary"
  )
  expect_equal(weights_idx, expected)

  # Test single column data frame
  ps_single <- data.frame(prob = c(0.1, 0.3, 0.7, 0.9))
  weights_single <- wt_ate(ps_single, exposure, exposure_type = "binary")
  expected_single <- wt_ate(ps_single$prob, exposure, exposure_type = "binary")
  expect_equal(weights_single, expected_single)

  # Test error with empty data frame
  expect_error(
    wt_ate(data.frame(), exposure),
    class = "propensity_df_ncol_error"
  )

  # Test error with invalid column selection
  expect_error(
    wt_ate(ps_df, exposure, .propensity_col = "nonexistent"),
    class = "propensity_df_column_error"
  )
})

test_that("all wt_* functions work with data frames", {
  ps_df <- data.frame(
    control = c(0.9, 0.7, 0.3, 0.1),
    treated = c(0.1, 0.3, 0.7, 0.9)
  )
  exposure <- c(0, 0, 1, 1)

  # Test ATT
  weights_att <- wt_att(ps_df, exposure, exposure_type = "binary")
  expected_att <- wt_att(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights_att, expected_att)

  # Test ATU
  weights_atu <- wt_atu(ps_df, exposure, exposure_type = "binary")
  expected_atu <- wt_atu(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights_atu, expected_atu)

  # Test ATM
  weights_atm <- wt_atm(ps_df, exposure, exposure_type = "binary")
  expected_atm <- wt_atm(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights_atm, expected_atm)

  # Test ATO
  weights_ato <- wt_ato(ps_df, exposure, exposure_type = "binary")
  expected_ato <- wt_ato(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights_ato, expected_ato)

  # Test Entropy
  weights_entropy <- wt_entropy(ps_df, exposure, exposure_type = "binary")
  expected_entropy <- wt_entropy(
    ps_df$treated,
    exposure,
    exposure_type = "binary"
  )
  expect_equal(weights_entropy, expected_entropy)
})

test_that("wt_* functions work with parsnip output", {
  skip_if_not_installed("parsnip")
  skip_on_cran()

  # Simulate data
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment_numeric <- rbinom(n, 1, ps_true)
  treatment_factor <- factor(treatment_numeric, levels = c("0", "1"))

  df <- data.frame(
    treatment = treatment_factor,
    x1 = x1,
    x2 = x2
  )

  # Fit model with parsnip
  ps_spec <- parsnip::logistic_reg()
  ps_spec <- parsnip::set_engine(ps_spec, "glm")
  ps_model <- parsnip::fit(ps_spec, treatment ~ x1 + x2, data = df)

  # Get predictions
  ps_preds <- predict(ps_model, df, type = "prob")

  # Test that it works with default (second column)
  weights_ate <- wt_ate(ps_preds, treatment_numeric, exposure_type = "binary")
  expect_s3_class(weights_ate, "psw")
  expect_equal(estimand(weights_ate), "ate")

  # Test explicit column selection with parsnip column names
  weights_ate2 <- wt_ate(
    ps_preds,
    treatment_numeric,
    .propensity_col = ".pred_1",
    exposure_type = "binary"
  )
  expect_equal(weights_ate, weights_ate2)

  # Test with all estimands
  expect_s3_class(
    wt_att(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
  expect_s3_class(
    wt_atu(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
  expect_s3_class(
    wt_atm(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
  expect_s3_class(
    wt_ato(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
  expect_s3_class(
    wt_entropy(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
})

test_that("wt_* functions work with GLM objects", {
  # Simulate data
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, ps_true)

  # Fit GLM model
  ps_model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Test ATE
  weights_ate_glm <- wt_ate(ps_model, treatment, exposure_type = "binary")
  weights_ate_numeric <- wt_ate(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_ate_glm, weights_ate_numeric)
  expect_s3_class(weights_ate_glm, "psw")
  expect_equal(estimand(weights_ate_glm), "ate")

  # Test ATT
  weights_att_glm <- wt_att(ps_model, treatment, exposure_type = "binary")
  weights_att_numeric <- wt_att(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_att_glm, weights_att_numeric)
  expect_s3_class(weights_att_glm, "psw")
  expect_equal(estimand(weights_att_glm), "att")

  # Test ATU
  weights_atu_glm <- wt_atu(ps_model, treatment, exposure_type = "binary")
  weights_atu_numeric <- wt_atu(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_atu_glm, weights_atu_numeric)
  expect_s3_class(weights_atu_glm, "psw")
  expect_equal(estimand(weights_atu_glm), "atu")

  # Test ATM
  weights_atm_glm <- wt_atm(ps_model, treatment, exposure_type = "binary")
  weights_atm_numeric <- wt_atm(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_atm_glm, weights_atm_numeric)
  expect_s3_class(weights_atm_glm, "psw")
  expect_equal(estimand(weights_atm_glm), "atm")

  # Test ATO
  weights_ato_glm <- wt_ato(ps_model, treatment, exposure_type = "binary")
  weights_ato_numeric <- wt_ato(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_ato_glm, weights_ato_numeric)
  expect_s3_class(weights_ato_glm, "psw")
  expect_equal(estimand(weights_ato_glm), "ato")

  # Test Entropy
  weights_entropy_glm <- wt_entropy(
    ps_model,
    treatment,
    exposure_type = "binary"
  )
  weights_entropy_numeric <- wt_entropy(
    ps_fitted,
    treatment,
    exposure_type = "binary"
  )
  expect_equal(weights_entropy_glm, weights_entropy_numeric)
  expect_s3_class(weights_entropy_glm, "psw")
  expect_equal(estimand(weights_entropy_glm), "entropy")
})

test_that("GLM methods handle continuous exposures", {
  # Simulate continuous exposure data
  set.seed(456)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  exposure <- 2 + 0.5 * x1 + 0.3 * x2 + rnorm(n)

  # Fit linear model
  exposure_model <- glm(exposure ~ x1 + x2, family = gaussian)

  # Test ATE with continuous exposure
  weights_ate <- wt_ate(
    exposure_model,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  expect_s3_class(weights_ate, "psw")
  expect_equal(estimand(weights_ate), "ate")
  expect_true(is_stabilized(weights_ate))

  # Check that weights are reasonable
  expect_true(all(is.finite(weights_ate)))
  expect_true(all(weights_ate > 0))
})

test_that("GLM methods error on non-GLM objects", {
  # Try with a non-GLM object
  expect_error(
    wt_ate("not a glm", c(0, 1, 0, 1)),
    class = "propensity_method_error"
  )

  expect_error(
    wt_att(list(a = 1, b = 2), c(0, 1, 0, 1)),
    class = "propensity_method_error"
  )
})

# Edge case tests for all wt_* functions ----

test_that("wt_* functions handle edge case propensity scores", {
  # Edge case: propensity scores at boundaries
  ps_boundary <- c(1e-10, 0.001, 0.5, 0.999, 1 - 1e-10)
  exposure_boundary <- c(0, 0, 1, 1, 1)

  # Test all functions with boundary values
  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(ps_boundary, exposure_boundary, exposure_type = "binary")
    expect_s3_class(weights, "psw")
    expect_true(all(is.finite(weights)))
    expect_true(all(weights > 0))
  }

  # Edge case: all propensity scores identical
  ps_constant <- rep(0.5, 5)
  exposure_mixed <- c(0, 0, 1, 1, 0)

  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(ps_constant, exposure_mixed, exposure_type = "binary")
    expect_s3_class(weights, "psw")
    # For constant propensity scores, weights should be constant within treatment groups
    expect_equal(length(unique(weights[exposure_mixed == 0])), 1)
    expect_equal(length(unique(weights[exposure_mixed == 1])), 1)
  }

  # Edge case: single observation
  ps_single <- 0.3
  exposure_single <- 1

  # Single observation needs explicit .treated/.untreated
  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(
      ps_single,
      exposure_single,
      exposure_type = "binary",
      .treated = 1
    )
    expect_s3_class(weights, "psw")
    expect_length(weights, 1)
    expect_true(is.finite(weights))
  }

  # Edge case: all treated or all control
  ps_various <- c(0.2, 0.4, 0.6, 0.8)

  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    # All treated - need to specify .treated since there's only one level
    weights_all_1 <- fn(
      ps_various,
      rep(1, 4),
      exposure_type = "binary",
      .treated = 1
    )
    expect_s3_class(weights_all_1, "psw")
    expect_true(all(is.finite(weights_all_1)))

    # All control - need to specify .untreated since there's only one level
    weights_all_0 <- fn(
      ps_various,
      rep(0, 4),
      exposure_type = "binary",
      .untreated = 0
    )
    expect_s3_class(weights_all_0, "psw")
    expect_true(all(is.finite(weights_all_0)))
  }
})

test_that("wt_* functions handle extreme weight scenarios", {
  # Create scenario that produces very large weights
  ps_extreme <- c(0.001, 0.999, 0.5, 0.5)
  exposure_extreme <- c(1, 0, 0, 1) # Misaligned with propensity

  # ATE should produce extreme weights
  weights_ate <- wt_ate(ps_extreme, exposure_extreme, exposure_type = "binary")
  expect_true(max(weights_ate) > 100) # Very large weight expected

  # ATO and ATM should be bounded
  weights_ato <- wt_ato(ps_extreme, exposure_extreme, exposure_type = "binary")
  weights_atm <- wt_atm(ps_extreme, exposure_extreme, exposure_type = "binary")
  expect_true(all(weights_ato <= 1)) # ATO weights bounded by 1
  expect_true(all(weights_atm <= 2)) # ATM weights bounded
})

# Additional error handling tests ----

test_that("wt_* functions error appropriately on invalid inputs", {
  # Invalid propensity scores (outside [0,1])
  expect_error(
    wt_ate(c(-0.1, 0.5, 1.1), c(0, 1, 0), exposure_type = "binary"),
    class = "propensity_range_error"
  )

  expect_error(
    wt_att(c(0, 0.5, 1), c(0, 1, 0), exposure_type = "binary"),
    class = "propensity_range_error"
  )

  # Length mismatch should error
  expect_error(
    wt_ate(c(0.1, 0.5), c(0, 1, 0), exposure_type = "binary"),
    class = "propensity_length_error"
  )

  # Invalid exposure type
  expect_error(
    wt_ate(c(0.1, 0.5), c(0, 1), exposure_type = "invalid"),
    "must be one of"
  )

  # Non-numeric propensity scores
  expect_error(
    wt_ate(c("a", "b"), c(0, 1)),
    class = "propensity_method_error"
  )

  # Empty vectors
  expect_error(
    wt_ate(numeric(0), numeric(0))
  )

  # Categorical exposure (not supported for most estimands)
  expect_error(
    wt_att(c(0.3, 0.3, 0.4), c(1, 2, 3), exposure_type = "categorical"),
    class = "propensity_wt_not_supported_error"
  )
})

test_that("data frame methods error appropriately", {
  # Empty data frame
  expect_error(
    wt_ate(data.frame(), c(0, 1)),
    class = "propensity_df_ncol_error"
  )

  # Non-existent column
  df <- data.frame(a = c(0.1, 0.9), b = c(0.9, 0.1))
  expect_error(
    wt_ate(df, c(0, 1), .propensity_col = "nonexistent"),
    class = "propensity_df_column_error"
  )

  # Column index out of bounds
  expect_error(
    wt_ate(df, c(0, 1), .propensity_col = 5),
    class = "propensity_df_column_error"
  )

  # Non-numeric column
  df_char <- data.frame(a = c("high", "low"), b = c(0.9, 0.1))
  suppressWarnings({
    expect_error(
      wt_ate(df_char, c(0, 1), .propensity_col = 1)
    )
  })

  # Column with invalid values
  df_invalid <- data.frame(a = c(0.5, 1.5), b = c(0.9, 0.1))
  expect_error(
    wt_ate(df_invalid, c(0, 1), .propensity_col = 1),
    class = "propensity_range_error"
  )
})

test_that("GLM methods error appropriately", {
  # Non-GLM object
  expect_error(
    wt_ate(lm(mpg ~ wt, data = mtcars), rep(0:1, 16)),
    class = "propensity_method_error"
  )

  # GLM with wrong dimensions
  small_glm <- glm(c(0, 1) ~ c(1, 2), family = binomial)
  expect_error(
    wt_ate(small_glm, c(0, 1, 0, 1)), # Mismatch in length
    class = "propensity_length_error"
  )
})

# Default method dispatch tests ----

test_that("default methods provide informative errors", {
  # Custom class that doesn't have a method
  custom_obj <- structure(list(x = 1:5), class = "my_custom_class")

  expect_error(
    wt_ate(custom_obj, c(0, 1, 0, 1, 0)),
    class = "propensity_method_error"
  )

  expect_error(
    wt_att(custom_obj, c(0, 1, 0, 1, 0)),
    class = "propensity_method_error"
  )

  expect_error(
    wt_atu(custom_obj, c(0, 1, 0, 1, 0)),
    class = "propensity_method_error"
  )

  expect_error(
    wt_atm(custom_obj, c(0, 1, 0, 1, 0)),
    class = "propensity_method_error"
  )

  expect_error(
    wt_ato(custom_obj, c(0, 1, 0, 1, 0)),
    class = "propensity_method_error"
  )

  expect_error(
    wt_entropy(custom_obj, c(0, 1, 0, 1, 0)),
    class = "propensity_method_error"
  )
})

# Data frame method tests with various column types ----

test_that("data frame methods handle various column configurations", {
  # Data frame with many columns
  ps_multi <- data.frame(
    col1 = runif(10, 0.1, 0.9),
    col2 = runif(10, 0.1, 0.9),
    col3 = runif(10, 0.1, 0.9),
    col4 = runif(10, 0.1, 0.9),
    col5 = runif(10, 0.1, 0.9)
  )
  exposure <- rbinom(10, 1, 0.5)

  # Test column selection by position
  for (i in 1:5) {
    weights <- suppressWarnings(wt_ate(
      ps_multi,
      exposure,
      .propensity_col = i,
      exposure_type = "binary"
    ))
    expected <- wt_ate(ps_multi[[i]], exposure, exposure_type = "binary")
    expect_equal(weights, expected)
  }

  # Test with tibble
  if (requireNamespace("tibble", quietly = TRUE)) {
    ps_tibble <- tibble::as_tibble(ps_multi)
    weights_tibble <- wt_ate(
      ps_tibble,
      exposure,
      .propensity_col = 3,
      exposure_type = "binary"
    )
    weights_df <- wt_ate(
      ps_multi,
      exposure,
      .propensity_col = 3,
      exposure_type = "binary"
    )
    expect_equal(weights_tibble, weights_df)
  }

  # Test with column names containing spaces
  ps_spaces <- data.frame(
    `control probability` = runif(5, 0.1, 0.9),
    `treatment probability` = runif(5, 0.1, 0.9),
    check.names = FALSE
  )
  exposure_small <- c(0, 0, 1, 1, 0)

  weights_spaces <- wt_ate(
    ps_spaces,
    exposure_small,
    .propensity_col = "treatment probability",
    exposure_type = "binary"
  )
  expect_s3_class(weights_spaces, "psw")

  # Test tidyselect helpers
  weights_last <- wt_ate(
    ps_multi,
    exposure,
    .propensity_col = tidyselect::last_col(),
    exposure_type = "binary"
  )
  expected_last <- wt_ate(ps_multi$col5, exposure, exposure_type = "binary")
  expect_equal(weights_last, expected_last)
})

# GLM method tests with different families ----

test_that("GLM methods handle non-binomial families appropriately", {
  # Gaussian family (for continuous exposures)
  set.seed(123)
  n <- 50
  x <- rnorm(n)
  exposure_cont <- 2 + 0.5 * x + rnorm(n)

  glm_gaussian <- glm(exposure_cont ~ x, family = gaussian)

  # Should work for ATE with continuous exposure
  weights_gaussian <- wt_ate(
    glm_gaussian,
    exposure_cont,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  expect_s3_class(weights_gaussian, "psw")
  expect_true(all(is.finite(weights_gaussian)))

  # Should error for estimands that don't support continuous
  # ATT doesn't accept continuous as a valid exposure type
  expect_error(
    wt_att(glm_gaussian, exposure_cont, exposure_type = "continuous"),
    "must be one of"
  )

  # Poisson family (should extract fitted values)
  # For non-binomial GLMs, the fitted values might not be valid propensity scores
  # Skip this test as it's not a valid use case
})

# Attribute preservation tests ----

test_that("attributes are preserved across all weight methods", {
  # Create trimmed propensity scores with attributes
  ps <- runif(20, 0.05, 0.95)
  exposure <- rbinom(20, 1, ps)

  ps_trimmed <- ps_trim(
    ps,
    .exposure = exposure,
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )

  # Test that all weight functions preserve trim attributes
  for (fn_name in c(
    "wt_ate",
    "wt_att",
    "wt_atu",
    "wt_atm",
    "wt_ato",
    "wt_entropy"
  )) {
    fn <- get(fn_name)

    # Suppress refit warning
    suppressWarnings({
      weights <- fn(ps_trimmed, exposure, exposure_type = "binary")
    })

    expect_true(is_ps_trimmed(weights))
    expect_equal(
      attr(weights, "ps_trim_meta"),
      attr(ps_trimmed, "ps_trim_meta")
    )
    expect_match(estimand(weights), "; trimmed$")
  }

  # Test with truncated propensity scores
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)

  for (fn_name in c(
    "wt_ate",
    "wt_att",
    "wt_atu",
    "wt_atm",
    "wt_ato",
    "wt_entropy"
  )) {
    fn <- get(fn_name)
    weights <- fn(ps_truncated, exposure, exposure_type = "binary")

    expect_true(is_ps_truncated(weights))
    expect_equal(
      attr(weights, "ps_trunc_meta"),
      attr(ps_truncated, "ps_trunc_meta")
    )
    expect_match(estimand(weights), "; truncated$")
  }
})

test_that("stabilization attributes are set correctly", {
  ps <- runif(20, 0.1, 0.9)
  exposure <- rbinom(20, 1, ps)

  # Test ATE with stabilization
  weights_unstab <- wt_ate(
    ps,
    exposure,
    stabilize = FALSE,
    exposure_type = "binary"
  )
  weights_stab <- wt_ate(
    ps,
    exposure,
    stabilize = TRUE,
    exposure_type = "binary"
  )

  expect_false(is_stabilized(weights_unstab))
  expect_true(is_stabilized(weights_stab))

  # Test with custom stabilization score
  weights_custom <- wt_ate(
    ps,
    exposure,
    stabilize = TRUE,
    stabilization_score = 0.4,
    exposure_type = "binary"
  )
  expect_true(is_stabilized(weights_custom))
})

# Stabilization tests across methods ----

test_that("stabilization works correctly for all applicable methods", {
  set.seed(456)
  ps <- runif(30, 0.2, 0.8)
  exposure <- rbinom(30, 1, ps)

  # ATE supports stabilization
  weights_ate_unstab <- wt_ate(
    ps,
    exposure,
    stabilize = FALSE,
    exposure_type = "binary"
  )
  weights_ate_stab <- wt_ate(
    ps,
    exposure,
    stabilize = TRUE,
    exposure_type = "binary"
  )

  # Check that stabilization was applied
  expect_true(is_stabilized(weights_ate_stab))
})

# NA handling tests ----

test_that("all methods handle NAs appropriately", {
  # Propensity scores with NAs
  ps_na <- c(0.2, NA, 0.5, 0.8, NA)
  exposure_na <- c(0, 1, 1, 0, 1)

  # All methods should handle NAs by producing NA weights
  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(ps_na, exposure_na, exposure_type = "binary")
    expect_s3_class(weights, "psw")
    expect_true(any(is.na(weights)))
  }

  # Exposure with NAs
  ps_good <- c(0.2, 0.3, 0.5, 0.8, 0.7)
  exposure_with_na <- c(0, NA, 1, 0, 1)

  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(ps_good, exposure_with_na, exposure_type = "binary")
    expect_s3_class(weights, "psw")
    expect_true(any(is.na(weights)))
  }

  # Data frame with NAs
  df_na <- data.frame(
    col1 = c(0.1, NA, 0.5),
    col2 = c(0.9, 0.5, NA)
  )

  # Data frame with NAs produces NA weights
  weights_df_na <- wt_ate(df_na, c(0, 1, 0), exposure_type = "binary")
  expect_s3_class(weights_df_na, "psw")
  expect_true(any(is.na(weights_df_na)))

  # GLM predictions with NAs
  set.seed(789)
  n <- 20
  x <- c(rnorm(18), NA, NA)
  y <- c(rbinom(18, 1, plogis(0.5 * x[1:18])), 0, 1)

  # GLM will handle NAs in predictors
  glm_na <- glm(y ~ x, family = binomial)

  # GLM with NA predictions will have shorter output
  # The NA observations are dropped during model fitting
  suppressWarnings({
    # GLM drops NA observations, so output is shorter
    expect_error(
      wt_ate(glm_na, y, exposure_type = "binary"),
      "must have the same length"
    )
  })
})

# Integration tests ----

test_that("weight functions integrate correctly with ps_trim and ps_trunc", {
  set.seed(999)
  n <- 100
  ps <- runif(n, 0.05, 0.95)
  exposure <- rbinom(n, 1, ps)

  # Create a model for refitting
  x <- rnorm(n)
  model <- glm(exposure ~ x + I(x^2), family = binomial)

  # Trim and refit
  ps_trim_obj <- ps_trim(
    ps,
    .exposure = exposure,
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )
  ps_refit_obj <- ps_refit(ps_trim_obj, model = model)

  # Should not warn after refit
  suppressMessages({
    weights_refit <- wt_ate(ps_refit_obj, exposure, exposure_type = "binary")
  })
  expect_true(is_ps_trimmed(weights_refit))
  expect_true(is_refit(weights_refit))

  # Truncate
  ps_trunc_obj <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)

  # Should not warn for truncation
  suppressMessages({
    weights_trunc <- wt_ate(ps_trunc_obj, exposure, exposure_type = "binary")
  })
  expect_true(is_ps_truncated(weights_trunc))
  expect_false(is_ps_trimmed(weights_trunc))
})

test_that("data frame and GLM methods produce consistent results", {
  set.seed(111)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  true_ps <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, true_ps)

  # Fit GLM
  model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(model)

  # Create data frame
  ps_df <- data.frame(
    control_prob = 1 - ps_fitted,
    treatment_prob = ps_fitted
  )

  # All methods should give same results
  for (fn_name in c(
    "wt_ate",
    "wt_att",
    "wt_atu",
    "wt_atm",
    "wt_ato",
    "wt_entropy"
  )) {
    fn <- get(fn_name)

    # From GLM
    weights_glm <- fn(model, treatment, exposure_type = "binary")

    # From fitted values
    weights_numeric <- fn(ps_fitted, treatment, exposure_type = "binary")

    # From data frame (using treatment column)
    weights_df <- fn(
      ps_df,
      treatment,
      .propensity_col = 2,
      exposure_type = "binary"
    )

    expect_equal(weights_glm, weights_numeric)
    expect_equal(weights_glm, weights_df)
  }
})

# Additional mathematical property tests ----

test_that("weight functions satisfy mathematical properties", {
  set.seed(222)
  n <- 100
  ps <- runif(n, 0.1, 0.9)
  treatment <- rbinom(n, 1, ps)

  # ATT weights: treated units should have weight 1
  weights_att <- wt_att(ps, treatment, exposure_type = "binary")
  expect_true(all(weights_att[treatment == 1] == 1))

  # ATU weights: control units should have weight 1
  weights_atu <- wt_atu(ps, treatment, exposure_type = "binary")
  expect_true(all(weights_atu[treatment == 0] == 1))

  # ATO weights: A * (1-e) + (1-A) * e
  weights_ato <- wt_ato(ps, treatment, exposure_type = "binary")
  # ATO weights are bounded by 1
  expect_true(all(weights_ato <= 1 + 1e-10))

  # ATM weights: symmetric for e and 1-e
  # Create symmetric propensity scores
  ps_sym <- c(0.3, 0.7, 0.4, 0.6)
  treatment_sym <- c(0, 1, 1, 0)
  weights_atm <- wt_atm(ps_sym, treatment_sym, exposure_type = "binary")

  # Weights for symmetric PS pairs should be equal
  expect_equal(weights_atm[1], weights_atm[2], tolerance = 1e-10)
  expect_equal(weights_atm[3], weights_atm[4], tolerance = 1e-10)
})

test_that("continuous exposure weights have correct properties", {
  set.seed(333)
  n <- 100
  x <- rnorm(n)
  exposure <- 2 + 0.5 * x + rnorm(n, sd = 1.5)

  # Fit model
  model <- lm(exposure ~ x)
  mu_hat <- fitted(model)
  sigma <- summary(model)$sigma

  # Calculate weights
  weights_cont <- wt_ate(
    mu_hat,
    .exposure = exposure,
    .sigma = rep(sigma, n),
    exposure_type = "continuous",
    stabilize = TRUE
  )

  expect_s3_class(weights_cont, "psw")
  expect_true(all(is.finite(weights_cont)))
  expect_true(all(weights_cont > 0))

  # Stabilized continuous weights should have mean approximately 1
  expect_equal(mean(weights_cont), 1, tolerance = 0.2)
})
