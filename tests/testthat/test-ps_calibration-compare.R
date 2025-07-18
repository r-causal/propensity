# Cross-validation tests against WeightIt and probably packages

test_that("ps_calibrate matches WeightIt::calibrate for Platt scaling", {
  skip_if_not_installed("WeightIt")
  
  set.seed(789)
  n <- 500
  # Create some realistic propensity scores with miscalibration
  X <- rnorm(n)
  true_ps <- plogis(0.5 * X)
  # Add miscalibration
  obs_ps <- plogis(qlogis(true_ps) + 0.3 + 0.2 * X)
  treat <- rbinom(n, 1, true_ps)
  
  # Our calibration
  our_calib <- ps_calibrate(obs_ps, treat)
  
  # WeightIt calibration (Platt method is logistic calibration)
  weightit_calib <- WeightIt::calibrate(obs_ps, treat, method = "platt")
  
  # Should be very close (allowing for numerical differences)
  expect_equal(as.numeric(our_calib), as.numeric(weightit_calib), tolerance = 1e-10)
})

test_that("ps_calibrate handles different treatment encodings like WeightIt", {
  skip_if_not_installed("WeightIt")
  
  set.seed(321)
  ps <- runif(100, 0.2, 0.8)
  treat_num <- rbinom(100, 1, ps)
  treat_char <- ifelse(treat_num == 1, "T", "C")
  
  # Our calibration with character treatment
  our_calib <- ps_calibrate(ps, treat_char, .treated = "T", .untreated = "C")
  
  # WeightIt with numeric treatment
  weightit_calib <- WeightIt::calibrate(ps, treat_num, method = "platt")
  
  expect_equal(as.numeric(our_calib), as.numeric(weightit_calib), tolerance = 1e-10)
})

test_that("compare calibration performance metrics with WeightIt", {
  skip_if_not_installed("WeightIt")
  
  set.seed(456)
  n <- 1000
  # Generate data with known miscalibration
  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, 0.5)
  true_ps <- plogis(-1 + 0.5 * X1 + X2)
  # Observed PS with systematic bias
  obs_ps <- plogis(qlogis(true_ps) + 0.5)
  treat <- rbinom(n, 1, true_ps)
  
  # Calibrate with both methods
  our_calib <- ps_calibrate(obs_ps, treat)
  weightit_calib <- WeightIt::calibrate(obs_ps, treat, method = "platt")
  
  # Check they produce identical results
  expect_equal(as.numeric(our_calib), as.numeric(weightit_calib), tolerance = 1e-10)
  
  # Test that calibration changes the distribution
  # (calibration doesn't always improve slope toward 1)
  expect_false(identical(as.numeric(obs_ps), as.numeric(our_calib)))
})

test_that("ps_calibrate produces similar results to probably package", {
  skip_if_not_installed("probably")
  skip_if_not_installed("tidyselect")
  
  set.seed(654)
  n <- 200
  # Create miscalibrated probabilities
  true_ps <- runif(n, 0.2, 0.8)
  obs_ps <- plogis(qlogis(true_ps) + 0.3)  # Add miscalibration
  treat <- rbinom(n, 1, true_ps)  # Use true ps for treatment
  
  # Our calibration
  our_calib <- ps_calibrate(obs_ps, treat)
  
  # probably calibration
  df <- data.frame(
    treat = factor(treat, levels = c("0", "1")),
    .pred_0 = 1 - obs_ps,
    .pred_1 = obs_ps
  )
  
  suppressWarnings({
    cal_data <- probably::cal_estimate_logistic(
      df,
      truth = treat,
      estimate = tidyselect::starts_with(".pred_")
    )
    df_cal <- probably::cal_apply(df, cal_data)
  })
  
  # Compare calibrated probabilities
  prob_calib <- df_cal$.pred_1
  
  # Both should be in valid range
  expect_true(all(our_calib >= 0 & our_calib <= 1))
  expect_true(all(prob_calib >= 0 & prob_calib <= 1))
  
  # Both should improve calibration compared to original
  # Calculate calibration errors (binned)
  n_bins <- 5
  bins <- cut(obs_ps, breaks = n_bins, include.lowest = TRUE)
  
  # Original calibration error
  orig_error <- mean(abs(tapply(treat, bins, mean, na.rm = TRUE) - 
                        tapply(obs_ps, bins, mean, na.rm = TRUE)), na.rm = TRUE)
  
  # Our calibration error
  our_bins <- cut(as.numeric(our_calib), breaks = n_bins, include.lowest = TRUE)
  our_error <- mean(abs(tapply(treat, our_bins, mean, na.rm = TRUE) - 
                       tapply(as.numeric(our_calib), our_bins, mean, na.rm = TRUE)), na.rm = TRUE)
  
  # probably calibration error
  prob_bins <- cut(prob_calib, breaks = n_bins, include.lowest = TRUE)
  prob_error <- mean(abs(tapply(treat, prob_bins, mean, na.rm = TRUE) - 
                        tapply(prob_calib, prob_bins, mean, na.rm = TRUE)), na.rm = TRUE)
  
  # Both should reduce calibration error compared to original
  expect_true(our_error <= orig_error + 0.05)  # Allow small tolerance
  expect_true(prob_error <= orig_error + 0.05)
  
  # The correlation between our calibration and probably's should be high
  expect_true(cor(as.numeric(our_calib), prob_calib) > 0.8)
})

test_that("extreme values handled consistently with WeightIt", {
  skip_if_not_installed("WeightIt")
  
  # Test with extreme propensity scores
  ps <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
  treat <- c(0, 0, 0, 1, 1, 1, 1)
  
  suppressWarnings({
    our_calib <- ps_calibrate(ps, treat)
    weightit_calib <- WeightIt::calibrate(ps, treat, method = "platt")
  })
  
  # Both should handle extreme values similarly
  expect_equal(as.numeric(our_calib), as.numeric(weightit_calib), tolerance = 1e-10)
})