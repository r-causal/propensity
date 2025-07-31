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
  expect_true(attr(w_atu_unfit, "trimmed"))
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
  expect_true(attr(w_atu_fit, "trimmed"))
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
  expect_true(attr(w_atm_unfit, "trimmed"))

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
  expect_true(attr(w_atm_fit, "trimmed"))
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
  expect_true(attr(w_ato_unfit, "trimmed"))

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
  expect_true(attr(w_ato_fit, "trimmed"))
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
  expect_true(attr(weights_ate, "stabilized"))

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
