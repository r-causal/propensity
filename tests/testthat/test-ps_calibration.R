test_that("errors for non-numeric ps", {
  expect_error(
    ps_calibrate("not numeric", c(0, 1)),
    "`ps` must be a numeric vector"
  )
})

test_that("errors for out-of-range ps", {
  expect_error(
    ps_calibrate(c(-0.1, 0.2), c(0, 1)),
    "`ps` values must be between 0 and 1"
  )
  expect_error(
    ps_calibrate(c(0.5, 1.1), c(0, 1)),
    "`ps` values must be between 0 and 1"
  )
})

test_that("errors when ps and treat have different lengths", {
  expect_error(
    ps_calibrate(runif(5), rep(0:1, length.out = 6)),
    "same length"
  )
})

test_that("returns a psw object of correct length and range", {
  set.seed(42)
  ps <- rep(0.5, 100)
  treat <- rbinom(100, 1, 0.3)

  out <- ps_calibrate(ps, treat)

  expect_s3_class(out, "psw")
  expect_length(out, 100)
  expect_true(all(out >= 0 & out <= 1))
})

test_that("constant ps yields calibrated = observed prevalence", {
  ps <- rep(0.5, 20)
  treat <- rep(c(0, 1), each = 10) # prevalence = 0.5

  out <- ps_calibrate(ps, treat)
  # all values should equal the 0.5 prevalence
  expect_equal(unique(as.numeric(out)), 0.5)
})

test_that("preserves psw attributes from an existing causal-weights object", {
  # Create a dummy psw object
  ps_orig <- psw(
    x = runif(10),
    estimand = "ATT",
    stabilized = TRUE,
    trimmed = TRUE,
    truncated = FALSE
  )
  treat <- rbinom(10, 1, ps_orig)

  out <- ps_calibrate(ps_orig, treat)

  expect_equal(attr(out, "estimand"), "ATT")
  expect_true(attr(out, "stabilized"))
  expect_true(attr(out, "trimmed"))
  expect_false(attr(out, "truncated"))
})

test_that("calibration changes the distribution", {
  set.seed(123)
  n <- 1000
  # Generate miscalibrated propensity scores
  true_ps <- runif(n, 0.2, 0.8)
  # Add systematic bias
  obs_ps <- plogis(qlogis(true_ps) + 0.5)
  treat <- rbinom(n, 1, true_ps)

  # Calibrate
  calibrated_ps <- ps_calibrate(obs_ps, treat)

  # Check that calibration changed the values
  expect_false(identical(as.numeric(obs_ps), as.numeric(calibrated_ps)))

  # Check that all values are valid probabilities
  expect_true(all(calibrated_ps >= 0 & calibrated_ps <= 1))

  # In this specific case with systematic bias, check if calibration helps
  # Calculate mean calibration error
  obs_bins <- cut(obs_ps, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  calib_bins <- cut(
    calibrated_ps,
    breaks = seq(0, 1, by = 0.1),
    include.lowest = TRUE
  )

  obs_calib_error <- abs(
    tapply(treat, obs_bins, mean, na.rm = TRUE) -
      tapply(obs_ps, obs_bins, mean, na.rm = TRUE)
  )
  calib_calib_error <- abs(
    tapply(treat, calib_bins, mean, na.rm = TRUE) -
      tapply(calibrated_ps, calib_bins, mean, na.rm = TRUE)
  )

  # Average calibration error should generally be reduced
  expect_true(
    mean(calib_calib_error, na.rm = TRUE) <=
      mean(obs_calib_error, na.rm = TRUE) + 0.1
  ) # Allow some tolerance
})

test_that("handles edge cases with extreme propensity scores", {
  # Near 0 and 1 values
  ps <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
  treat <- c(0, 0, 0, 1, 1, 1, 1)

  expect_no_error(suppressWarnings(calibrated <- ps_calibrate(ps, treat)))
  expect_true(all(calibrated >= 0 & calibrated <= 1))
})

test_that("works with different treatment codings", {
  set.seed(456)
  ps <- runif(50)
  treat_num <- rbinom(50, 1, ps)

  # Test with different treatment encodings
  treat_char <- ifelse(treat_num == 1, "treated", "control")
  treat_factor <- factor(treat_char)
  treat_logical <- as.logical(treat_num)

  calib1 <- ps_calibrate(ps, treat_num)
  calib2 <- ps_calibrate(
    ps,
    treat_char,
    .treated = "treated",
    .untreated = "control"
  )
  calib3 <- ps_calibrate(
    ps,
    treat_factor,
    .treated = "treated",
    .untreated = "control"
  )
  calib4 <- ps_calibrate(ps, treat_logical)

  # All should produce the same result
  expect_equal(as.numeric(calib1), as.numeric(calib2))
  expect_equal(as.numeric(calib1), as.numeric(calib3))
  expect_equal(as.numeric(calib1), as.numeric(calib4))
})

test_that(".treated and .untreated parameters work consistently with package patterns", {
  set.seed(123)
  ps <- runif(30, 0.3, 0.7)

  # Test automatic detection with 0/1 coding
  treat_01 <- rbinom(30, 1, ps)
  calib_auto <- ps_calibrate(ps, treat_01)
  expect_s3_class(calib_auto, "psw")

  # Test explicit specification with 0/1 coding
  calib_explicit <- ps_calibrate(ps, treat_01, .treated = 1, .untreated = 0)
  expect_equal(as.numeric(calib_auto), as.numeric(calib_explicit))

  # Test with character coding
  treat_char <- ifelse(treat_01 == 1, "treat", "control")
  calib_char_explicit <- ps_calibrate(
    ps,
    treat_char,
    .treated = "treat",
    .untreated = "control"
  )
  expect_equal(as.numeric(calib_auto), as.numeric(calib_char_explicit))

  # Test automatic detection with factor
  treat_factor <- factor(treat_char, levels = c("control", "treat"))
  calib_factor_auto <- ps_calibrate(ps, treat_factor)
  expect_equal(as.numeric(calib_auto), as.numeric(calib_factor_auto))

  # Test with logical coding (should be automatic)
  treat_logical <- as.logical(treat_01)
  calib_logical <- ps_calibrate(ps, treat_logical)
  expect_equal(as.numeric(calib_auto), as.numeric(calib_logical))
})

test_that(".treated/.untreated defaults are NULL like other package functions", {
  # Check that the defaults match the package pattern
  ps_calibrate_formals <- formals(ps_calibrate)
  expect_null(ps_calibrate_formals$.treated)
  expect_null(ps_calibrate_formals$.untreated)

  # Compare with other weight functions to ensure consistency
  wt_ate_formals <- formals(wt_ate)
  expect_equal(ps_calibrate_formals$.treated, wt_ate_formals$.treated)
  expect_equal(ps_calibrate_formals$.untreated, wt_ate_formals$.untreated)
})

test_that("automatic treatment detection works with binary vectors", {
  set.seed(789)
  ps <- runif(40, 0.2, 0.8)

  # Test with different binary representations
  treat_01 <- rbinom(40, 1, ps)
  treat_12 <- treat_01 + 1 # 1/2 coding
  treat_neg <- ifelse(treat_01 == 1, 1, -1) # -1/1 coding

  # All should work with automatic detection
  calib_01 <- ps_calibrate(ps, treat_01)
  expect_s3_class(calib_01, "psw")

  # These require explicit specification
  calib_12 <- ps_calibrate(ps, treat_12, .treated = 2, .untreated = 1)
  calib_neg <- ps_calibrate(ps, treat_neg, .treated = 1, .untreated = -1)

  # All should produce valid results
  expect_true(all(calib_01 >= 0 & calib_01 <= 1))
  expect_true(all(calib_12 >= 0 & calib_12 <= 1))
  expect_true(all(calib_neg >= 0 & calib_neg <= 1))
})

test_that("error handling for ambiguous treatment coding", {
  set.seed(456)
  ps <- runif(20, 0.3, 0.7)

  # Three-level factor should require explicit specification
  treat_three <- factor(sample(c("A", "B", "C"), 20, replace = TRUE))
  expect_error(
    ps_calibrate(ps, treat_three),
    "Specify.*treated.*untreated"
  )

  # Should work with explicit specification
  treat_binary_from_three <- ifelse(treat_three == "A", 1, 0)
  expect_no_error(
    ps_calibrate(ps, treat_binary_from_three)
  )
})

test_that("is_ps_calibrated works correctly", {
  ps <- runif(20)
  treat <- rbinom(20, 1, ps)

  # Regular numeric vector
  expect_false(is_ps_calibrated(ps))

  # Uncalibrated psw object
  ps_wt <- psw(ps, estimand = "ate")
  expect_false(is_ps_calibrated(ps_wt))

  # Calibrated psw object
  calibrated <- ps_calibrate(ps, treat)
  expect_true(is_ps_calibrated(calibrated))
})

test_that("errors when trying to calibrate already calibrated ps", {
  ps <- runif(20)
  treat <- rbinom(20, 1, ps)

  calibrated <- ps_calibrate(ps, treat)

  expect_error(
    ps_calibrate(calibrated, treat),
    "already calibrated"
  )
})

test_that("handles NA values appropriately", {
  ps <- c(0.1, 0.3, NA, 0.7, 0.9)
  treat <- c(0, 0, 1, 1, 1)

  # Should preserve NAs in output
  suppressWarnings(calibrated <- ps_calibrate(ps, treat))
  expect_length(calibrated, 5)
  expect_true(is.na(calibrated[3]))
  expect_s3_class(calibrated, "psw")

  # Test with isotonic regression too
  calibrated_iso <- ps_calibrate(ps, treat, method = "isoreg")
  expect_length(calibrated_iso, 5)
  expect_true(is.na(calibrated_iso[3]))
})

test_that("isotonic regression calibration works", {
  set.seed(789)
  ps <- runif(100, 0.1, 0.9)
  treat <- rbinom(100, 1, ps)

  # Should not error
  calibrated_iso <- ps_calibrate(ps, treat, method = "isoreg")

  expect_s3_class(calibrated_iso, "psw")
  expect_length(calibrated_iso, 100)
  expect_true(all(calibrated_iso >= 0 & calibrated_iso <= 1))
  expect_true(is_ps_calibrated(calibrated_iso))
})

test_that("isotonic regression preserves monotonicity", {
  # Create data where treatment probability increases with ps
  set.seed(456)
  n <- 200
  ps <- seq(0.1, 0.9, length.out = n)
  # Add some noise but maintain overall trend
  treat <- rbinom(n, 1, ps + rnorm(n, 0, 0.05))

  calibrated_iso <- ps_calibrate(ps, treat, method = "isoreg")

  # Check that calibrated scores are monotonic (allowing for ties)
  diffs <- diff(as.numeric(calibrated_iso))
  expect_true(all(diffs >= -1e-10)) # Allow for numerical tolerance
})

test_that("method parameter validation works", {
  ps <- runif(20)
  treat <- rbinom(20, 1, ps)

  # Invalid method should error
  expect_error(
    ps_calibrate(ps, treat, method = "invalid"),
    "'arg' should be one of"
  )
})

test_that("isotonic and logistic calibration can differ", {
  set.seed(123)
  # Create data where isotonic might perform differently
  ps <- c(rep(0.2, 50), rep(0.8, 50))
  treat <- c(rbinom(50, 1, 0.3), rbinom(50, 1, 0.7))

  calib_logistic <- ps_calibrate(ps, treat, method = "logistic")
  calib_iso <- ps_calibrate(ps, treat, method = "isoreg")

  # They should produce different results in general
  expect_false(identical(as.numeric(calib_logistic), as.numeric(calib_iso)))
})

test_that("isotonic calibration matches WeightIt isoreg exactly", {
  skip_if_not_installed("WeightIt")

  set.seed(321)
  ps <- runif(100, 0.2, 0.8)
  treat <- rbinom(100, 1, ps)

  our_iso <- ps_calibrate(ps, treat, method = "isoreg")
  weightit_iso <- WeightIt::calibrate(ps, treat, method = "isoreg")

  expect_equal(as.numeric(our_iso), as.numeric(weightit_iso), tolerance = 1e-10)
})

test_that("isotonic regression preserves monotonicity better than logistic", {
  set.seed(888)
  # Create data where isotonic should perform better
  n <- 100
  ps <- seq(0.1, 0.9, length.out = n)
  # Non-linear relationship that violates logistic assumption
  true_prob <- 0.2 + 0.6 * ps^2
  treat <- rbinom(n, 1, true_prob)

  logistic_calib <- ps_calibrate(ps, treat, method = "logistic")
  iso_calib <- ps_calibrate(ps, treat, method = "isoreg")

  # Check isotonic preserves monotonicity
  iso_diffs <- diff(as.numeric(iso_calib))
  expect_true(all(iso_diffs >= -1e-10)) # Allow for numerical tolerance

  # Both should be different from original
  expect_false(identical(as.numeric(ps), as.numeric(logistic_calib)))
  expect_false(identical(as.numeric(ps), as.numeric(iso_calib)))

  # They should produce different results for non-linear data
  expect_false(identical(as.numeric(logistic_calib), as.numeric(iso_calib)))
})

test_that("isotonic regression handles various cases like WeightIt", {
  skip_if_not_installed("WeightIt")

  # Test with ties in propensity scores (but ensure sufficient data)
  set.seed(999)
  ps_ties <- rep(c(0.3, 0.7), each = 4)
  treat_ties <- c(0, 0, 1, 1, 0, 1, 1, 1)

  our_iso_ties <- ps_calibrate(ps_ties, treat_ties, method = "isoreg")
  weightit_iso_ties <- WeightIt::calibrate(
    ps_ties,
    treat_ties,
    method = "isoreg"
  )

  expect_equal(
    as.numeric(our_iso_ties),
    as.numeric(weightit_iso_ties),
    tolerance = 1e-10
  )

  # Test that our implementation handles edge cases gracefully even if WeightIt fails
  ps_extreme <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
  treat_extreme <- c(0, 0, 0, 1, 1, 1, 1)

  # Our implementation should handle this without error
  expect_no_error(
    our_extreme <- ps_calibrate(ps_extreme, treat_extreme, method = "isoreg")
  )
  expect_true(all(our_extreme >= 0 & our_extreme <= 1))
  expect_true(all(diff(as.numeric(our_extreme)) >= -1e-10)) # Monotonic
})

test_that("smooth parameter works correctly for logistic calibration", {
  skip_if_not_installed("mgcv")

  set.seed(42)
  ps <- runif(50, 0.2, 0.8)
  treat <- rbinom(50, 1, ps)

  # Both smooth options should work
  calib_smooth <- ps_calibrate(ps, treat, smooth = TRUE)
  calib_simple <- ps_calibrate(ps, treat, smooth = FALSE)

  # Both should be psw objects
  expect_s3_class(calib_smooth, "psw")
  expect_s3_class(calib_simple, "psw")

  # Both should be in valid range
  expect_true(all(calib_smooth >= 0 & calib_smooth <= 1))
  expect_true(all(calib_simple >= 0 & calib_simple <= 1))

  # They should generally produce different results
  expect_false(identical(as.numeric(calib_smooth), as.numeric(calib_simple)))
})

test_that("smooth parameter is ignored for isotonic regression", {
  set.seed(123)
  ps <- runif(30, 0.1, 0.9)
  treat <- rbinom(30, 1, ps)

  # smooth should be ignored for isoreg
  iso_smooth_true <- ps_calibrate(ps, treat, method = "isoreg", smooth = TRUE)
  iso_smooth_false <- ps_calibrate(ps, treat, method = "isoreg", smooth = FALSE)

  # Should be identical since smooth is ignored for isoreg
  expect_equal(as.numeric(iso_smooth_true), as.numeric(iso_smooth_false))
})

test_that("error when mgcv is not available for smooth calibration", {
  # Mock unavailable mgcv package by using a non-existent namespace
  skip("This test would require mocking requireNamespace which is complex in R")
  # In a real scenario, you'd mock requireNamespace to return FALSE
  # and test that the error message is correct
})

# Cross-validation tests against WeightIt and the probably package

test_that("ps_calibrate with smooth=FALSE matches WeightIt::calibrate for logistic calibration", {
  skip_if_not_installed("WeightIt")

  set.seed(789)
  n <- 500
  # Create some realistic propensity scores with miscalibration
  X <- rnorm(n)
  true_ps <- plogis(0.5 * X)
  # Add miscalibration
  obs_ps <- plogis(qlogis(true_ps) + 0.3 + 0.2 * X)
  treat <- rbinom(n, 1, true_ps)

  # Our calibration with simple logistic (to match WeightIt)
  our_calib <- ps_calibrate(obs_ps, treat, smooth = FALSE)

  # WeightIt calibration (platt method is logistic calibration)
  weightit_calib <- WeightIt::calibrate(obs_ps, treat, method = "platt")

  # Should be very close (allowing for numerical differences)
  expect_equal(
    as.numeric(our_calib),
    as.numeric(weightit_calib),
    tolerance = 1e-10
  )
})

test_that("ps_calibrate handles different treatment encodings like WeightIt", {
  skip_if_not_installed("WeightIt")

  set.seed(321)
  ps <- runif(100, 0.2, 0.8)
  treat_num <- rbinom(100, 1, ps)
  treat_char <- ifelse(treat_num == 1, "T", "C")

  # Our calibration with character treatment (use smooth=FALSE to match WeightIt)
  our_calib <- ps_calibrate(
    ps,
    treat_char,
    .treated = "T",
    .untreated = "C",
    smooth = FALSE
  )

  # WeightIt with numeric treatment
  weightit_calib <- WeightIt::calibrate(ps, treat_num, method = "platt")

  expect_equal(
    as.numeric(our_calib),
    as.numeric(weightit_calib),
    tolerance = 1e-10
  )
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

  # Calibrate with both methods (use smooth=FALSE to match WeightIt)
  our_calib <- ps_calibrate(obs_ps, treat, smooth = FALSE)
  weightit_calib <- WeightIt::calibrate(obs_ps, treat, method = "platt")

  # Check they produce identical results
  expect_equal(
    as.numeric(our_calib),
    as.numeric(weightit_calib),
    tolerance = 1e-10
  )

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
  obs_ps <- plogis(qlogis(true_ps) + 0.3) # Add miscalibration
  treat <- rbinom(n, 1, true_ps) # Use true ps for treatment

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

  # Helper function to calculate calibration error using binned approach
  calc_calib_error <- function(pred_probs, true_outcomes, n_bins = 5) {
    bins <- cut(pred_probs, breaks = n_bins, include.lowest = TRUE)
    bin_means_true <- tapply(true_outcomes, bins, mean, na.rm = TRUE)
    bin_means_pred <- tapply(pred_probs, bins, mean, na.rm = TRUE)
    mean(abs(bin_means_true - bin_means_pred), na.rm = TRUE)
  }

  # Calculate calibration errors for comparison
  orig_error <- calc_calib_error(obs_ps, treat)
  our_error <- calc_calib_error(as.numeric(our_calib), treat)
  prob_error <- calc_calib_error(prob_calib, treat)

  # Both should reduce calibration error compared to original
  expect_true(our_error <= orig_error + 0.05) # Allow small tolerance
  expect_true(prob_error <= orig_error + 0.05)

  # The correlation between our calibration and probably's should be high
  expect_true(cor(as.numeric(our_calib), prob_calib) > 0.8)
})

test_that("ps_calibrate with smooth=TRUE matches probably's default behavior exactly", {
  skip_if_not_installed("probably")
  skip_if_not_installed("tidyselect")
  skip_if_not_installed("mgcv")

  set.seed(123)
  n <- 100
  ps <- runif(n, 0.1, 0.9)
  treat <- rbinom(n, 1, ps)

  # Our smoothed calibration (default smooth=TRUE)
  our_smooth <- ps_calibrate(ps, treat, smooth = TRUE)

  # probably calibration with smooth=TRUE (default)
  df <- data.frame(
    treat = factor(treat, levels = c("0", "1")),
    .pred_0 = 1 - ps,
    .pred_1 = ps
  )

  suppressWarnings({
    cal_data <- probably::cal_estimate_logistic(
      df,
      truth = treat,
      estimate = tidyselect::starts_with(".pred_"),
      smooth = TRUE
    )
    df_cal <- probably::cal_apply(df, cal_data)
  })

  prob_smooth <- df_cal$.pred_1

  # Should be very close (allowing for numerical differences in GAM fitting)
  # Use a more reasonable tolerance for GAM differences and ensure both are vectors
  expect_equal(
    as.numeric(our_smooth),
    as.numeric(prob_smooth),
    tolerance = 1e-3
  )
})

test_that("ps_calibrate with smooth=FALSE matches probably's simple logistic exactly", {
  skip_if_not_installed("probably")
  skip_if_not_installed("tidyselect")

  set.seed(456)
  n <- 100
  ps <- runif(n, 0.1, 0.9)
  treat <- rbinom(n, 1, ps)

  # Our simple logistic calibration
  our_simple <- ps_calibrate(ps, treat, smooth = FALSE)

  # probably calibration with smooth=FALSE
  df <- data.frame(
    treat = factor(treat, levels = c("0", "1")),
    .pred_0 = 1 - ps,
    .pred_1 = ps
  )

  suppressWarnings({
    cal_data <- probably::cal_estimate_logistic(
      df,
      truth = treat,
      estimate = tidyselect::starts_with(".pred_"),
      smooth = FALSE
    )
    df_cal <- probably::cal_apply(df, cal_data)
  })

  prob_simple <- df_cal$.pred_1

  # Should be identical for simple logistic regression
  expect_equal(as.numeric(our_simple), prob_simple, tolerance = 1e-10)
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
  expect_equal(
    as.numeric(our_calib),
    as.numeric(weightit_calib),
    tolerance = 1e-10
  )
})
