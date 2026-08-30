test_that("errors for non-numeric ps", {
  expect_propensity_error(
    ps_calibrate("not numeric", c(0, 1))
  )
})

test_that("errors for out-of-range ps", {
  expect_propensity_error(
    ps_calibrate(c(-0.1, 0.2), c(0, 1))
  )
  expect_propensity_error(
    ps_calibrate(c(0.5, 1.1), c(0, 1))
  )
})

test_that("errors when ps and .exposure have different lengths", {
  expect_propensity_error(
    ps_calibrate(runif(5), rep(0:1, length.out = 6))
  )
})

test_that("returns a ps_calib object of correct length and range", {
  set.seed(42)
  ps <- rep(0.5, 100)
  treat <- rbinom(100, 1, 0.3)

  out <- ps_calibrate(ps, treat)

  expect_s3_class(out, "ps_calib")
  expect_length(out, 100)
  expect_true(all(as.numeric(out) >= 0 & as.numeric(out) <= 1))
})

test_that("constant ps yields calibrated = observed prevalence", {
  ps <- rep(0.5, 20)
  treat <- rep(c(0, 1), each = 10) # prevalence = 0.5

  out <- ps_calibrate(ps, treat)
  # all values should equal the 0.5 prevalence
  expect_equal(unique(as.numeric(out)), 0.5)
})

test_that("calibration metadata is properly stored", {
  ps <- runif(10)
  treat <- rbinom(10, 1, ps)

  # Test logistic calibration with smooth = TRUE
  out_smooth <- ps_calibrate(ps, treat, method = "logistic", smooth = TRUE)
  meta_smooth <- ps_calib_meta(out_smooth)
  expect_equal(meta_smooth$method, "logistic")
  expect_true(meta_smooth$smooth)

  # Test logistic calibration with smooth = FALSE
  out_simple <- ps_calibrate(ps, treat, method = "logistic", smooth = FALSE)
  meta_simple <- ps_calib_meta(out_simple)
  expect_equal(meta_simple$method, "logistic")
  expect_false(meta_simple$smooth)

  # Test isotonic regression
  out_iso <- ps_calibrate(ps, treat, method = "isoreg")
  meta_iso <- ps_calib_meta(out_iso)
  expect_equal(meta_iso$method, "isoreg")
  expect_false(meta_iso$smooth)
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
  expect_true(all(
    as.numeric(calibrated_ps) >= 0 & as.numeric(calibrated_ps) <= 1
  ))

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
  expect_true(all(as.numeric(calibrated) >= 0 & as.numeric(calibrated) <= 1))
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
    .focal_level = "treated",
    .reference_level = "control"
  )
  calib3 <- ps_calibrate(
    ps,
    treat_factor,
    .focal_level = "treated",
    .reference_level = "control"
  )
  calib4 <- ps_calibrate(ps, treat_logical)

  # All should produce the same result
  expect_equal(as.numeric(calib1), as.numeric(calib2))
  expect_equal(as.numeric(calib1), as.numeric(calib3))
  expect_equal(as.numeric(calib1), as.numeric(calib4))
})

test_that(".focal_level and .reference_level parameters work consistently with package patterns", {
  set.seed(123)
  ps <- runif(30, 0.3, 0.7)

  # Test automatic detection with 0/1 coding
  treat_01 <- rbinom(30, 1, ps)
  calib_auto <- ps_calibrate(ps, treat_01)
  expect_s3_class(calib_auto, "ps_calib")

  # Test explicit specification with 0/1 coding
  calib_explicit <- ps_calibrate(
    ps,
    treat_01,
    .focal_level = 1,
    .reference_level = 0
  )
  expect_equal(as.numeric(calib_auto), as.numeric(calib_explicit))

  # Test with character coding
  treat_char <- ifelse(treat_01 == 1, "treat", "control")
  calib_char_explicit <- ps_calibrate(
    ps,
    treat_char,
    .focal_level = "treat",
    .reference_level = "control"
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

test_that("ps_calibrate refuses a focal level the exposure never takes", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c("a", "b", "a", "b")

  # A focal level nobody holds leaves every unit in the reference group, so the
  # calibration model is fit against an outcome that never varies.
  expect_error(
    ps_calibrate(ps, exposure, .focal_level = "absent"),
    class = "propensity_focal_level_error"
  )
})

test_that(".focal_level/.reference_level defaults are NULL like other package functions", {
  # Check that the defaults match the package pattern
  ps_calibrate_formals <- formals(ps_calibrate)
  expect_null(ps_calibrate_formals$.focal_level)
  expect_null(ps_calibrate_formals$.reference_level)

  # Compare with other weight functions to ensure consistency
  wt_ate_formals <- formals(wt_ate)
  expect_equal(ps_calibrate_formals$.focal_level, wt_ate_formals$.focal_level)
  expect_equal(
    ps_calibrate_formals$.reference_level,
    wt_ate_formals$.reference_level
  )
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
  expect_s3_class(calib_01, "ps_calib")

  # These require explicit specification
  calib_12 <- ps_calibrate(ps, treat_12, .focal_level = 2, .reference_level = 1)
  calib_neg <- ps_calibrate(
    ps,
    treat_neg,
    .focal_level = 1,
    .reference_level = -1
  )

  # All should produce valid results
  expect_true(all(as.numeric(calib_01) >= 0 & as.numeric(calib_01) <= 1))
  expect_true(all(as.numeric(calib_12) >= 0 & as.numeric(calib_12) <= 1))
  expect_true(all(as.numeric(calib_neg) >= 0 & as.numeric(calib_neg) <= 1))
})

test_that("error handling for ambiguous treatment coding", {
  set.seed(456)
  ps <- runif(20, 0.3, 0.7)

  # Three-level factor should require explicit specification
  treat_three <- factor(sample(c("A", "B", "C"), 20, replace = TRUE))
  expect_propensity_error(
    ps_calibrate(ps, treat_three)
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

  expect_propensity_error(
    ps_calibrate(calibrated, treat)
  )
})

test_that("handles NA values appropriately", {
  ps <- c(0.1, 0.3, NA, 0.7, 0.9)
  treat <- c(0, 0, 1, 1, 1)

  # Should preserve NAs in output
  suppressWarnings(calibrated <- ps_calibrate(ps, treat))
  expect_length(calibrated, 5)
  expect_true(is.na(calibrated[3]))
  expect_s3_class(calibrated, "ps_calib")

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

  expect_s3_class(calibrated_iso, "ps_calib")
  expect_length(calibrated_iso, 100)
  expect_true(all(
    as.numeric(calibrated_iso) >= 0 & as.numeric(calibrated_iso) <= 1
  ))
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
  expect_propensity_error(
    ps_calibrate(ps, treat, method = "invalid")
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
  our_extreme <- expect_no_error(
    ps_calibrate(ps_extreme, treat_extreme, method = "isoreg")
  )
  expect_true(all(as.numeric(our_extreme) >= 0 & as.numeric(our_extreme) <= 1))
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
  expect_s3_class(calib_smooth, "ps_calib")
  expect_s3_class(calib_simple, "ps_calib")

  # Both should be in valid range
  expect_true(all(
    as.numeric(calib_smooth) >= 0 & as.numeric(calib_smooth) <= 1
  ))
  expect_true(all(
    as.numeric(calib_simple) >= 0 & as.numeric(calib_simple) <= 1
  ))

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

# Cross-validation tests against WeightIt and the probably package

# `probably::cal_estimate_logistic()` loads butcher the first time it is called,
# and butcher and generics both register `as.character.dev_topic`. propensity
# has already loaded generics by then, so the second registration prints an S3
# method overwrite note. The note is package-load output rather than a
# condition a test can catch, so the calls that trigger it are made through this
# helper, which quiets it along with the warnings the calibration fits raise.
quiet_calibration <- function(expr) {
  suppressPackageStartupMessages(suppressWarnings(expr))
}

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
    .focal_level = "T",
    .reference_level = "C",
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

  quiet_calibration({
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
  expect_true(all(as.numeric(our_calib) >= 0 & as.numeric(our_calib) <= 1))
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

  quiet_calibration({
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

  quiet_calibration({
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

# Standalone pava_weighted() tests ------------------------------------------

test_that("pava_weighted returns input unchanged when already non-decreasing", {
  y <- c(0.1, 0.3, 0.5, 0.8, 1.0)
  x <- seq_along(y)
  result <- pava_weighted(x, y)
  expect_equal(result, y)
})

test_that("pava_weighted merges violating pairs", {
  # y = c(0, 1, 0, 1): middle pair violates, should merge to 0.5
  result <- pava_weighted(1:4, c(0, 1, 0, 1))
  expect_equal(result, c(0, 0.5, 0.5, 1))
})

test_that("pava_weighted handles all-constant y", {
  result <- pava_weighted(1:5, rep(0.5, 5))
  expect_equal(result, rep(0.5, 5))
})

test_that("pava_weighted handles single observation", {
  result <- pava_weighted(1, 0.7)
  expect_equal(result, 0.7)
})

test_that("pava_weighted handles completely decreasing y", {
  result <- pava_weighted(1:4, c(1, 0.75, 0.5, 0.25))
  # All should merge to the grand mean
  expect_equal(result, rep(mean(c(1, 0.75, 0.5, 0.25)), 4))
})

test_that("pava_weighted respects observation weights", {
  # Two observations: y = c(1, 0) with equal weights -> mean = 0.5
  result_equal <- pava_weighted(1:2, c(1, 0), w = c(1, 1))
  expect_equal(result_equal, c(0.5, 0.5))

  # Same but with weight 3 on first obs: weighted mean = (1*3 + 0*1)/4 = 0.75
  result_weighted <- pava_weighted(1:2, c(1, 0), w = c(3, 1))
  expect_equal(result_weighted, c(0.75, 0.75))
})

test_that("pava_weighted handles tied x values", {
  # Tied x-values should NOT be grouped (unlike stats::isoreg)
  x <- c(1, 1, 2, 2)
  y <- c(0, 1, 0, 1)
  result <- pava_weighted(x, y)
  # Each observation is its own block initially; result must be monotonic
  expect_true(all(diff(result) >= -1e-10))
})

test_that("pava_weighted preserves original order", {
  x <- c(3, 1, 2)
  y <- c(0.9, 0.1, 0.5)
  result <- pava_weighted(x, y)
  # After ordering by x: (1, 0.1), (2, 0.5), (3, 0.9) - already non-decreasing
  # So result should be the same as input
  expect_equal(result, y)
})

# An exposure with dimensions reaches the same coding the weight functions use,
# where its cells were read in storage order rather than one value per
# observation.

test_that("ps_calibrate refuses an exposure with dimensions", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  dimensioned <- matrix(c(1, 0, 1, 0), nrow = 2, ncol = 2)

  expect_error(
    ps_calibrate(ps, dimensioned),
    class = "propensity_binary_transform_error"
  )
})

# What a comparison with a calibrated vector asks ----------------------------

# The scores this fixture calibrates span the unit interval and none of them
# land on a whole number, which is what the integer routes below turn on. The
# exposure alternates enough for the calibration model to fit without a seed.

calib_vctrs_fixture <- function(method = "logistic", smooth = FALSE) {
  ps_calibrate(
    c(0.05, 0.15, 0.3, 0.45, 0.55, 0.7, 0.85, 0.95),
    c(0, 0, 1, 0, 1, 1, 0, 1),
    method = method,
    smooth = smooth
  )
}

# A comparison asks about the values a vector holds, not about the type holding
# them, so it has no metadata to lose and nothing to announce. Left to the
# default vctrs route it goes through the common type of the two sides, which
# for a calibrated vector and a number is the numeric downgrade, and so every
# comparison reports dropping metadata that the answer never depended on. The
# sibling weight class answers comparisons from its own values and stays
# silent; a calibrated vector owes the same.

test_that("comparing a ps_calib with a number does not warn about class", {
  cal <- calib_vctrs_fixture()
  values <- as.numeric(cal)

  expect_no_warning(
    {
      expect_identical(cal > 0.5, values > 0.5)
      expect_identical(cal >= 0.5, values >= 0.5)
      expect_identical(cal < 0.5, values < 0.5)
      expect_identical(cal <= 0.5, values <= 0.5)
      expect_identical(cal == values[[1]], values == values[[1]])
      expect_identical(cal != values[[1]], values != values[[1]])
    },
    class = "propensity_class_downgrade_warning"
  )
})

test_that("comparing a ps_calib yields a plain logical", {
  cal <- calib_vctrs_fixture()

  expect_no_warning(
    {
      greater <- cal > 0.5
      equal <- cal == as.numeric(cal)[[1]]
    },
    class = "propensity_class_downgrade_warning"
  )

  expect_type(greater, "logical")
  expect_type(equal, "logical")
  expect_false(inherits(greater, "ps_calib"))
  expect_false(inherits(equal, "ps_calib"))
  expect_null(attributes(greater))
  expect_null(attributes(equal))
})

test_that("ps_calib comparisons enforce vctrs strict size semantics", {
  # The size half of this contract holds today and guards against a comparison
  # that answers directly falling back to base R recycling: anything other than
  # equal lengths, or one side of length 1, has no answer.
  cal <- calib_vctrs_fixture()
  short <- cal[1:3]

  expect_error(cal == short, class = "vctrs_error_incompatible_size")
  expect_error(cal > short, class = "vctrs_error_incompatible_size")
  expect_error(cal != short, class = "vctrs_error_incompatible_size")

  # Broadcasting a length-1 right side still works, and is the path a caller
  # asking which scores clear a threshold takes.
  expect_no_warning(
    {
      expect_identical(cal > 0.5, as.numeric(cal) > 0.5)
      expect_identical(cal[1] == cal, as.numeric(cal)[[1]] == as.numeric(cal))
    },
    class = "propensity_class_downgrade_warning"
  )
})

test_that("a refused comparison names the code that wrote it", {
  cal <- calib_vctrs_fixture()
  short <- cal[1:3]

  # The recycling the comparison operators enforce happens inside a helper the
  # caller has no way to reach, and vctrs raises the refusal from further down
  # still. Left to itself the refusal names that helper, which says nothing to
  # whoever wrote the comparison, so the frame it reports against travels down
  # with the sizes.
  compare_in_fn <- function(x, y) x == y
  cnd <- rlang::catch_cnd(compare_in_fn(cal, short), classes = "error")

  expect_s3_class(cnd, "vctrs_error_incompatible_size")
  expect_identical(
    paste(deparse(conditionCall(cnd)[[1]]), collapse = " "),
    "compare_in_fn"
  )
})

# What calibration reads from scores that carry a record ---------------------

# `ps_calibrate()` accepts propensity scores that have already been trimmed or
# truncated, and the first thing it does with them is compare them against 0
# and 1. Comparing a classed vector routes through the numeric downgrade, so
# the range check announces a conversion the caller never asked for, once per
# comparison, before any calibration happens. The values read are the same
# either way, and so is the calibration built from them.

test_that("calibrating trimmed scores does not warn about class", {
  ps <- c(0.05, 0.15, 0.3, 0.45, 0.55, 0.7, 0.85, 0.95)
  exposure <- c(0, 0, 1, 0, 1, 1, 0, 1)
  trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
  from_numeric <- ps_calibrate(as.numeric(trimmed), exposure, smooth = FALSE)

  expect_no_warning(
    {
      from_trimmed <- ps_calibrate(trimmed, exposure, smooth = FALSE)
      expect_s3_class(from_trimmed, "ps_calib")
      expect_equal(as.numeric(from_trimmed), as.numeric(from_numeric))
      expect_identical(ps_calib_meta(from_trimmed), ps_calib_meta(from_numeric))
    },
    class = "propensity_class_downgrade_warning"
  )
})

test_that("calibrating truncated scores does not warn about class", {
  ps <- c(0.05, 0.15, 0.3, 0.45, 0.55, 0.7, 0.85, 0.95)
  exposure <- c(0, 0, 1, 0, 1, 1, 0, 1)
  truncated <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)
  from_numeric <- ps_calibrate(as.numeric(truncated), exposure, smooth = FALSE)

  expect_no_warning(
    {
      from_truncated <- ps_calibrate(truncated, exposure, smooth = FALSE)
      expect_s3_class(from_truncated, "ps_calib")
      expect_equal(as.numeric(from_truncated), as.numeric(from_numeric))
      expect_identical(
        ps_calib_meta(from_truncated),
        ps_calib_meta(from_numeric)
      )
    },
    class = "propensity_class_downgrade_warning"
  )
})

test_that("calibrating trimmed scores by isotonic regression is also silent", {
  ps <- c(0.05, 0.15, 0.3, 0.45, 0.55, 0.7, 0.85, 0.95)
  exposure <- c(0, 0, 1, 0, 1, 1, 0, 1)
  trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
  from_numeric <- ps_calibrate(
    as.numeric(trimmed),
    exposure,
    method = "isoreg"
  )

  # Isotonic calibration reaches the scores through a different fit than the
  # logistic path, so the values it reads are worth pinning on their own.
  expect_no_warning(
    {
      from_trimmed <- ps_calibrate(trimmed, exposure, method = "isoreg")
      expect_equal(as.numeric(from_trimmed), as.numeric(from_numeric))
    },
    class = "propensity_class_downgrade_warning"
  )
})

# Where a calibrated vector and an integer meet ------------------------------

# Calibrated propensity scores lie strictly between 0 and 1, so the integers
# have no room for them: a combination that met an integer in the integers
# would be every score rounded away, and a cast to integer loses every value it
# was given. The sibling classes resolve an integer in the doubles, announcing
# the downgrade once, and refuse the cast to integer as lossy. A calibrated
# vector that has no answer at all for an integer is the odd one out.

test_that("combining a ps_calib with an integer keeps the calibrated scores", {
  cal <- calib_vctrs_fixture()
  combined <- expect_propensity_warning(vec_c(cal, 1L))

  expect_type(combined, "double")
  expect_equal(combined, c(as.numeric(cal), 1))
})

test_that("combining an integer with a ps_calib keeps the calibrated scores", {
  cal <- calib_vctrs_fixture()
  combined <- expect_propensity_warning(vec_c(1L, cal))

  expect_type(combined, "double")
  expect_equal(combined, c(1, as.numeric(cal)))
})

test_that("casting a ps_calib to integer refuses rather than rounds", {
  cal <- calib_vctrs_fixture()

  expect_error(
    vec_cast(cal, integer()),
    class = "vctrs_error_cast_lossy"
  )
})

test_that("casting an integer to ps_calib keeps the calibration of the target", {
  cal <- calib_vctrs_fixture()
  out <- vec_cast(c(0L, 1L), to = cal)

  expect_s3_class(out, "ps_calib")
  expect_equal(as.numeric(out), c(0, 1))
  expect_identical(ps_calib_meta(out), ps_calib_meta(cal))
})

test_that("comparing a ps_calib with an integer matches comparing a double", {
  cal <- calib_vctrs_fixture()

  expect_no_warning(
    {
      expect_identical(cal == 1L, cal == 1)
      expect_identical(cal > 0L, cal > 0)
    },
    class = "propensity_class_downgrade_warning"
  )
})

# What a cast to a calibrated vector owes its target -------------------------

# A cast returns the values it was handed in the type it was handed, and a
# calibrated vector's type is how the calibration was performed. A cast that
# writes its own description over the target's hands back a vector claiming a
# calibration nothing carried out, under a method name no argument accepts.

test_that("casting a double to ps_calib keeps the method of the target", {
  cal <- calib_vctrs_fixture()
  out <- vec_cast(c(0.3, 0.4), to = cal)

  expect_s3_class(out, "ps_calib")
  expect_equal(as.numeric(out), c(0.3, 0.4))
  expect_identical(ps_calib_meta(out), ps_calib_meta(cal))
  expect_identical(ps_calib_meta(out)$method, "logistic")
  expect_false(ps_calib_meta(out)$smooth)
})

test_that("casting a double to an isotonic ps_calib keeps that method", {
  cal <- calib_vctrs_fixture(method = "isoreg")
  out <- vec_cast(c(0.3, 0.4), to = cal)

  expect_identical(ps_calib_meta(out)$method, "isoreg")
})

test_that("casting a double to ps_calib keeps the smoothing of the target", {
  # Built directly so the spline fit, and with it mgcv, stays out of a test
  # about what the cast copies.
  to <- new_ps_calib(
    double(),
    ps_calib_meta = list(method = "logistic", smooth = TRUE)
  )
  out <- vec_cast(c(0.3, 0.4), to = to)

  expect_true(ps_calib_meta(out)$smooth)
})

# What the isotonic fit produces ---------------------------------------------

# Isotonic calibration pools adjacent violators separately within each exposure
# group and reads every unit's calibrated score from the fit for the group it
# belongs to, so each fitted value is the mean of a pooled block of zeros and
# ones rather than the output of an optimizer. The values below characterize
# that fit as it stands, inside the unit interval and at both of its ends, so
# any later change to the shape of the calibration curve, such as a guard
# against reading past the scores actually observed, has to be a deliberate
# one.

calib_isoreg_interior_fixture <- function() {
  ps_calibrate(
    c(0.05, 0.15, 0.3, 0.45, 0.55, 0.7, 0.85, 0.95),
    c(0, 0, 1, 0, 1, 1, 0, 1),
    method = "isoreg"
  )
}

# Scores hugging both ends of the interval, endpoints included, which
# calibration accepts.
calib_isoreg_boundary_fixture <- function() {
  ps_calibrate(
    c(0, 0.001, 0.02, 0.35, 0.5, 0.65, 0.98, 0.999, 1),
    c(0, 1, 0, 0, 1, 0, 1, 0, 1),
    method = "isoreg"
  )
}

# Regression guard: this fit holds today and has to survive any change to the
# isotonic route.
test_that("isotonic calibration fits pooled block means on interior scores", {
  expect_equal(
    as.numeric(calib_isoreg_interior_fixture()),
    c(
      0,
      0,
      0.5,
      0.5,
      0.666666666666667,
      0.666666666666667,
      0.666666666666667,
      1
    ),
    tolerance = 1e-12
  )
})

# Regression guard: the same fit read at the ends of the interval, where a
# change to the isotonic route would show up first.
test_that("isotonic calibration fits pooled block means on boundary scores", {
  calibrated <- calib_isoreg_boundary_fixture()

  expect_equal(
    as.numeric(calibrated),
    c(
      0,
      0.333333333333333,
      0.333333333333333,
      0.333333333333333,
      0.5,
      0.5,
      0.5,
      0.5,
      1
    ),
    tolerance = 1e-12
  )
  expect_true(all(diff(as.numeric(calibrated)) >= 0))
})

# What the fallback from a spline to a straight line announces ---------------

# `smooth = TRUE` asks for a spline, and a spline needs enough distinct scores
# to place its knots, so with fewer than ten the calibration falls back to a
# plain logistic regression. That changes which model produced the scores the
# caller is handed, and today it happens without a word: the returned metadata
# records the fallback, but nothing at the point of the call says the spline
# that was asked for was never fit.

calib_few_unique_fixture <- function() {
  list(
    ps = rep(c(0.2, 0.4, 0.6, 0.8), each = 5),
    exposure = c(
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      1,
      1,
      0,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      0
    )
  )
}

test_that("ps_calibrate announces falling back from a spline", {
  skip_if_not_installed("mgcv")
  withr::local_options(propensity.quiet = FALSE)
  fixture <- calib_few_unique_fixture()

  messages <- testthat::capture_messages(
    ps_calibrate(fixture$ps, fixture$exposure, smooth = TRUE)
  )

  # One message, naming the smoothing that was asked for and not carried out.
  # The wording itself is pinned by the snapshot below.
  expect_length(messages, 1)
  expect_match(messages, "smooth|spline", ignore.case = TRUE)
})

test_that("the fallback message is recorded", {
  skip_if_not_installed("mgcv")
  withr::local_options(propensity.quiet = FALSE)
  fixture <- calib_few_unique_fixture()

  expect_snapshot(
    calibrated <- ps_calibrate(fixture$ps, fixture$exposure, smooth = TRUE)
  )
})

test_that("the fallback is silent when informational output is turned off", {
  skip_if_not_installed("mgcv")
  # `alert_info()` respects `propensity.quiet`, which the suite sets, so the
  # announcement must not reach a caller who has asked for quiet.
  withr::local_options(propensity.quiet = TRUE)
  fixture <- calib_few_unique_fixture()

  expect_no_message(
    ps_calibrate(fixture$ps, fixture$exposure, smooth = TRUE)
  )
})

# Regression guard: the fallback fits the model it says it fits, which holds
# today and is what the announcement will describe.
test_that("the fallback fits the logistic model it records", {
  skip_if_not_installed("mgcv")
  fixture <- calib_few_unique_fixture()
  calibrated <- ps_calibrate(fixture$ps, fixture$exposure, smooth = TRUE)
  straight_line <- glm(
    exposure ~ ps,
    data = data.frame(exposure = fixture$exposure, ps = fixture$ps),
    family = binomial()
  )

  expect_false(ps_calib_meta(calibrated)$smooth)
  expect_equal(
    as.numeric(calibrated),
    unname(fitted(straight_line)),
    tolerance = 1e-12
  )
})

# Regression guard: with enough distinct scores the spline is fit, so there is
# no fallback to announce.
test_that("a spline fit with enough distinct scores announces nothing", {
  skip_if_not_installed("mgcv")
  withr::local_options(propensity.quiet = FALSE)
  ps <- seq(0.1, 0.9, length.out = 20)
  exposure <- rep(c(0, 1), 10)

  calibrated <- expect_no_message(ps_calibrate(ps, exposure, smooth = TRUE))
  expect_true(ps_calib_meta(calibrated)$smooth)
})

# Where calibration draws the ends of the unit interval ----------------------

# Trimming and weighting need scores strictly inside the unit interval, because
# a zero or a one divides into a weight nothing can use. Calibration is the one
# route that accepts them, since repairing scores a model pushed to the ends is
# part of what it is for. The logistic calibration curve maps them back inside
# the interval, which is what the fixture below is calibrated by. Isotonic
# calibration can return a score at an endpoint when its pooled block is pure,
# as the boundary characterization above records, and such scores are rejected
# by the weight functions.

calib_boundary_policy_fixture <- function() {
  list(
    ps = c(0, seq(0.05, 0.95, length.out = 20), 1),
    exposure = c(
      0,
      0,
      0,
      1,
      0,
      0,
      1,
      0,
      0,
      1,
      0,
      1,
      1,
      0,
      1,
      1,
      0,
      1,
      1,
      1,
      0,
      1
    )
  )
}

# Regression guard: the inclusive range is deliberate and has to keep holding.
test_that("calibration accepts scores at the ends of the unit interval", {
  fixture <- calib_boundary_policy_fixture()
  calibrated <- ps_calibrate(fixture$ps, fixture$exposure, smooth = FALSE)

  expect_s3_class(calibrated, "ps_calib")
  expect_length(calibrated, length(fixture$ps))
  expect_true(all(as.numeric(calibrated) > 0))
  expect_true(all(as.numeric(calibrated) < 1))
})

# Regression guard, and the contrast that makes the one above a choice rather
# than an oversight: every other route still refuses the endpoints.
test_that("trimming and weighting still refuse scores at the ends of the unit interval", {
  fixture <- calib_boundary_policy_fixture()

  expect_error(
    ps_trim(fixture$ps, method = "ps"),
    class = "propensity_range_error"
  )
  expect_error(
    wt_ate(fixture$ps, fixture$exposure),
    class = "propensity_range_error"
  )
})

# What calibration does with an exposure that has gaps -----------------------

# A unit with no recorded exposure tells the calibration model nothing, so the
# fit drops it, and that is the whole of its effect. The fitted curve maps a
# propensity score to a calibrated one, and the score of a unit with no
# exposure is as readable as any other, so the prediction covers every row.
# Calibrating introduces no gap the scores did not already have.

calib_na_exposure_fixture <- function() {
  list(
    ps = c(0.1, 0.25, 0.4, 0.5, 0.62, 0.75, 0.9),
    exposure = c(0, 1, NA, 0, 1, 1, 0)
  )
}

test_that("calibration predicts over every row when the exposure has gaps", {
  fixture <- calib_na_exposure_fixture()
  calibrated <- ps_calibrate(fixture$ps, fixture$exposure, smooth = FALSE)

  model_data <- data.frame(exposure = fixture$exposure, ps = fixture$ps)
  complete_case <- glm(
    exposure ~ ps,
    data = model_data[!is.na(fixture$exposure), ],
    family = binomial()
  )
  over_all_rows <- as.numeric(predict(
    complete_case,
    newdata = model_data,
    type = "response"
  ))

  expect_length(calibrated, length(fixture$ps))
  expect_false(anyNA(as.numeric(calibrated)))
  expect_equal(as.numeric(calibrated), over_all_rows, tolerance = 1e-12)
})

# Regression guard: the spline route already predicts over every row, which is
# the behavior the straight-line route above owes.
test_that("a spline calibration predicts over every row when the exposure has gaps", {
  skip_if_not_installed("mgcv")
  ps <- seq(0.1, 0.9, length.out = 20)
  exposure <- rep(c(0, 1), 10)
  exposure[7] <- NA

  calibrated <- ps_calibrate(ps, exposure, smooth = TRUE)

  expect_length(calibrated, length(ps))
  expect_false(anyNA(as.numeric(calibrated)))
})

# Regression guard and contrast: isotonic calibration reads a unit's score from
# the fit for the exposure group it belongs to, and a unit with no exposure
# belongs to neither, so it has no calibrated score to take.
test_that("isotonic calibration leaves a unit with no exposure uncalibrated", {
  fixture <- calib_na_exposure_fixture()
  calibrated <- as.numeric(
    ps_calibrate(fixture$ps, fixture$exposure, method = "isoreg")
  )

  expect_length(calibrated, length(fixture$ps))
  expect_true(is.na(calibrated[[3]]))
  expect_false(anyNA(calibrated[-3]))
})

# What calibration does with a matrix of scores ------------------------------

# The route through `ps_calibrate()` reads `ps` as one score per unit. A matrix
# of scores, the shape the categorical routes take, has no reading here: with
# one column it is flattened and its dimensions are dropped without a word, and
# with more than one the length check reports a count of cells against a count
# of units. Neither tells the caller that the shape is what is wrong. The same
# error the function already raises for a `ps` that is not a numeric vector is
# the answer a matrix deserves too.

test_that("ps_calibrate refuses a matrix of propensity scores", {
  exposure <- c(0, 1, 0, 1)
  one_column <- matrix(c(0.2, 0.4, 0.6, 0.8), ncol = 1)
  two_columns <- cbind(c(0.2, 0.4, 0.6, 0.8), c(0.8, 0.6, 0.4, 0.2))

  expect_error(
    ps_calibrate(one_column, exposure, smooth = FALSE),
    class = "propensity_type_error"
  )
  expect_error(
    ps_calibrate(two_columns, exposure, smooth = FALSE),
    class = "propensity_type_error"
  )
})

test_that("ps_calibrate reports a matrix of propensity scores informatively", {
  expect_propensity_error(
    ps_calibrate(
      cbind(c(0.2, 0.4, 0.6, 0.8), c(0.8, 0.6, 0.4, 0.2)),
      c(0, 1, 0, 1),
      smooth = FALSE
    )
  )
})

# Regression guard: a data frame of scores is already refused under the class
# a matrix is pinned to above.
test_that("ps_calibrate refuses a data frame of propensity scores", {
  expect_error(
    ps_calibrate(
      data.frame(a = c(0.2, 0.4, 0.6, 0.8), b = c(0.8, 0.6, 0.4, 0.2)),
      c(0, 1, 0, 1),
      smooth = FALSE
    ),
    class = "propensity_type_error"
  )
})

# One dimension is one score per unit, which is the shape calibration reads.
# The refusal above is about holding more than that, not about carrying a `dim`
# attribute, and `wt_ate()` reads a one-dimensional array too.
test_that("ps_calibrate accepts a one-dimensional array of propensity scores", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c(0, 1, 0, 1)
  one_d <- array(ps, dim = length(ps))

  calibrated <- ps_calibrate(one_d, exposure, smooth = FALSE)

  expect_s3_class(calibrated, "ps_calib")
  expect_equal(
    as.numeric(calibrated),
    as.numeric(ps_calibrate(ps, exposure, smooth = FALSE)),
    tolerance = 1e-12
  )
})

# Naming the propensity scores `.propensity` ----------------------------------

# The weight functions read the propensity scores from `.propensity` and this
# one reads them from `ps`, so a call written against one is refused by the
# other in both directions. The tests below pin the scores under the new name,
# the deprecated shim that keeps the old name working for a release, and the
# refusal to read both names at once. The positional pin comes first: whatever
# else the rename moves, it must not move what a call that names nothing
# returns.
#
# Unlike `ps_trim()` and `ps_trunc()`, `ps_calibrate()` is a plain function
# rather than a generic, so there is no dispatch to preserve. What stands in
# for the dispatch pins is the refusal of a matrix of scores, which names the
# argument and so has to follow the rename.

calib_rename_scores <- function() {
  c(0.10, 0.18, 0.27, 0.35, 0.44, 0.52, 0.61, 0.69, 0.78, 0.86, 0.92, 0.96)
}

calib_rename_exposure <- function() {
  c(0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1)
}

# Smoothing is turned off so the calibration curve is a logistic regression of
# the exposure on the scores, which the test can fit for itself as an oracle.
calib_rename_positional <- function() {
  ps_calibrate(
    calib_rename_scores(),
    calib_rename_exposure(),
    method = "logistic",
    smooth = FALSE
  )
}

test_that("ps_calibrate() calibrates a positional vector of scores", {
  out <- calib_rename_positional()

  expect_s3_class(out, "ps_calib")
  expect_length(out, 12)

  oracle_data <- data.frame(
    treat = calib_rename_exposure(),
    ps = calib_rename_scores()
  )
  oracle <- glm(treat ~ ps, data = oracle_data, family = binomial)
  expect_equal(
    as.numeric(out),
    unname(predict(oracle, type = "response"))
  )

  meta <- ps_calib_meta(out)
  expect_equal(meta$method, "logistic")
  expect_false(meta$smooth)
})

test_that("ps_calibrate() reads the propensity scores from .propensity", {
  expect_equal(
    ps_calibrate(
      .propensity = calib_rename_scores(),
      .exposure = calib_rename_exposure(),
      method = "logistic",
      smooth = FALSE
    ),
    calib_rename_positional()
  )
})

test_that("ps_calibrate() deprecates the propensity scores under ps", {
  with_always_deprecated({
    expect_warning(
      ps_calibrate(
        ps = calib_rename_scores(),
        .exposure = calib_rename_exposure(),
        method = "logistic",
        smooth = FALSE
      ),
      class = "lifecycle_warning_deprecated"
    )
  })

  # The old name still has to reach the same calibration, not merely warn. The
  # deprecation is pinned above, so it is silenced here rather than repeated.
  withr::local_options(lifecycle_verbosity = "quiet")
  expect_equal(
    ps_calibrate(
      ps = calib_rename_scores(),
      .exposure = calib_rename_exposure(),
      method = "logistic",
      smooth = FALSE
    ),
    calib_rename_positional()
  )
})

test_that("ps_calibrate() refuses the propensity scores under both names", {
  withr::local_options(lifecycle_verbosity = "quiet")

  # The condition subclass is the shim's to choose; what this pins is that the
  # refusal is one of the package's own errors and that it names both spellings,
  # so the caller can see which one to drop.
  err <- expect_error(
    ps_calibrate(
      .propensity = calib_rename_scores(),
      ps = calib_rename_scores(),
      .exposure = calib_rename_exposure(),
      method = "logistic",
      smooth = FALSE
    ),
    class = "propensity_error"
  )

  msg <- conditionMessage(err)
  expect_match(msg, "`.propensity`", fixed = TRUE)
  expect_match(msg, "`ps`", fixed = TRUE)
})

test_that("ps_calibrate() names .propensity when refusing a matrix of scores", {
  err <- expect_error(
    ps_calibrate(
      .propensity = matrix(calib_rename_scores(), nrow = 6),
      .exposure = calib_rename_exposure()[1:6],
      method = "logistic",
      smooth = FALSE
    ),
    class = "propensity_type_error"
  )

  msg <- conditionMessage(err)
  expect_match(msg, "`.propensity`", fixed = TRUE)
  expect_false(grepl("`ps`", msg, fixed = TRUE))
})

test_that("ps_calibrate() names .propensity in the length mismatch error", {
  err <- expect_error(
    ps_calibrate(
      .propensity = calib_rename_scores(),
      .exposure = calib_rename_exposure()[-1],
      method = "logistic",
      smooth = FALSE
    ),
    class = "propensity_length_error"
  )

  msg <- conditionMessage(err)
  expect_match(msg, "`.propensity`", fixed = TRUE)
  expect_false(grepl("`ps`", msg, fixed = TRUE))
})

# What a refused cast says, and what the fallback needs ------------------------

test_that("casting between differently calibrated ps_calib objects names the disagreement", {
  x <- calib_vctrs_fixture()

  # Built directly so the spline fit, and with it mgcv, stays out of a test
  # about what the cast compares.
  smoothed <- new_ps_calib(
    vec_data(x),
    ps_calib_meta = list(method = "logistic", smooth = TRUE)
  )

  # A `ps_calib` is printed with its method but not with whether the fit was
  # smoothed, so these two types render identically and a refusal that names
  # neither reads as a type that cannot be converted to itself.
  expect_identical(vec_ptype_full(x), vec_ptype_full(smoothed))
  expect_error(
    vec_cast(x, to = smoothed),
    regexp = "different calibration parameters",
    class = "vctrs_error_incompatible_type"
  )

  # The combine says what disagrees, and the cast makes the same comparison, so
  # it says the same thing.
  isotonic <- new_ps_calib(
    vec_data(x),
    ps_calib_meta = list(method = "isoreg", smooth = FALSE)
  )
  expect_error(
    vec_cast(x, to = isotonic),
    regexp = "different calibration parameters",
    class = "vctrs_error_incompatible_type"
  )
})

test_that("ps_calibrate() falls back to a straight line without reaching for mgcv", {
  # Five distinct scores are too few to place the knots of a spline, so the
  # fallback fits a logistic regression and mgcv takes no part in the
  # calibration. Checking for the package before working that out refuses the
  # call over a dependency the calibration it performs has no use for.
  scores <- rep(c(0.2, 0.3, 0.4, 0.5, 0.6), each = 4)
  exposure <- rep(c(0, 1), 10)

  local_mocked_bindings(
    check_installed = function(...) rlang::abort("mgcv was checked for"),
    .package = "rlang"
  )

  calibrated <- ps_calibrate(scores, exposure, smooth = TRUE)

  expect_s3_class(calibrated, "ps_calib")
  expect_false(ps_calib_meta(calibrated)$smooth)

  oracle <- glm(
    treat ~ ps,
    data = data.frame(treat = exposure, ps = scores),
    family = binomial()
  )
  expect_equal(
    as.numeric(calibrated),
    unname(predict(oracle, type = "response"))
  )
})

test_that("ps_calibrate() names the propensity scores it was not given", {
  err <- expect_error(
    ps_calibrate(.exposure = c(0, 1, 0, 1)),
    class = "propensity_missing_arg_error"
  )
  expect_match(conditionMessage(err), "`.propensity`", fixed = TRUE)
})

# ---- fitted-model methods ---------------------------------------------------

# `ps_calibrate()` reads a fitted propensity score model the way the weight
# functions read one: the scores come off the fit, and so does the exposure,
# which calibration always needs, unless the caller names an exposure of their
# own. Calibration reads one score per observation, so a binomial fit and a
# two-level multinomial fit are what it accepts; a multinomial fit of three or
# more levels reports one column per level and is refused as a model of
# something calibration has no reading for.

calib_model_data <- local({
  set.seed(20250930)

  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  odds_b <- exp(0.9 * x1 - 0.5 * x2)
  odds_c <- exp(-0.7 * x1 + 0.8 * x2)
  total <- 1 + odds_b + odds_c
  p_b <- odds_b / total
  p_c <- odds_c / total
  u <- runif(n)

  data.frame(
    x1 = x1,
    x2 = x2,
    trt = factor(ifelse(u < p_b, "b", ifelse(u < p_b + p_c, "c", "a"))),
    z = rbinom(n, 1, plogis(1.6 * x1 - 0.9 * x2)),
    a2 = factor(ifelse(
      runif(n) < plogis(1.2 * x1 - 0.6 * x2),
      "control",
      "treated"
    ))
  )
})

calib_binary_fit <- function() {
  glm(z ~ x1 + x2, data = calib_model_data, family = binomial())
}

calib_categorical_fit <- function() {
  nnet::multinom(trt ~ x1 + x2, data = calib_model_data, trace = FALSE)
}

calib_two_level_fit <- function() {
  nnet::multinom(a2 ~ x1 + x2, data = calib_model_data, trace = FALSE)
}

expect_same_calibration <- function(from_model, oracle) {
  testthat::expect_equal(
    as.numeric(from_model),
    as.numeric(oracle),
    tolerance = 1e-12
  )
  testthat::expect_identical(class(from_model), class(oracle))
  testthat::expect_identical(ps_calib_meta(from_model), ps_calib_meta(oracle))
}

test_that("ps_calibrate() calibrates the scores a binomial fit reports", {
  fit <- calib_binary_fit()
  scores <- predict(fit, type = "response")

  expect_same_calibration(
    ps_calibrate(fit, smooth = FALSE),
    ps_calibrate(scores, calib_model_data$z, smooth = FALSE)
  )
  expect_same_calibration(
    ps_calibrate(fit, method = "isoreg"),
    ps_calibrate(scores, calib_model_data$z, method = "isoreg")
  )
})

test_that("an explicit .exposure wins over the one a fit carries in ps_calibrate()", {
  fit <- calib_binary_fit()
  scores <- predict(fit, type = "response")
  reordered <- rev(calib_model_data$z)

  expect_same_calibration(
    ps_calibrate(fit, reordered, smooth = FALSE),
    ps_calibrate(scores, reordered, smooth = FALSE)
  )
  expect_false(isTRUE(all.equal(
    as.numeric(ps_calibrate(fit, reordered, smooth = FALSE)),
    as.numeric(ps_calibrate(fit, smooth = FALSE))
  )))
})

test_that("ps_calibrate() reads a two-level multinomial fit on the binary path", {
  skip_if_not_installed("nnet")

  fit <- calib_two_level_fit()
  scores <- as.numeric(fitted(fit))

  expect_same_calibration(
    ps_calibrate(fit, smooth = FALSE),
    ps_calibrate(scores, calib_model_data$a2, smooth = FALSE)
  )
})

test_that("ps_calibrate() refuses a multinomial fit of three or more levels", {
  skip_if_not_installed("nnet")

  expect_error(
    ps_calibrate(calib_categorical_fit(), smooth = FALSE),
    class = "propensity_model_family_error"
  )
})

test_that("ps_calibrate() refuses a fit it cannot read propensity scores from", {
  linear <- lm(z ~ x1 + x2, data = calib_model_data)

  expect_error(
    ps_calibrate(linear, calib_model_data$z, smooth = FALSE),
    class = "propensity_method_error"
  )
  expect_error(
    ps_calibrate(
      structure(list(), class = "not_a_model"),
      calib_model_data$z,
      smooth = FALSE
    ),
    class = "propensity_method_error"
  )
})

test_that("ps_calibrate() names the class of a fit it has no reading for", {
  expect_propensity_error(
    ps_calibrate(
      lm(z ~ x1 + x2, data = calib_model_data),
      calib_model_data$z,
      smooth = FALSE
    )
  )
})

test_that("ps_calibrate() reports a multinomial fit of too many levels", {
  skip_if_not_installed("nnet")

  expect_propensity_error(ps_calibrate(calib_categorical_fit(), smooth = FALSE))
})

# A missing exposure used to be reported by R against the argument it could not
# force, which now names the method dispatch reached rather than the function
# the caller wrote.
test_that("ps_calibrate() names the exposure it was not given", {
  expect_propensity_error(ps_calibrate(c(0.2, 0.4, 0.6, 0.8)))
})

# A fit reports the probability of the level the response's default coding
# treats as focal, so naming the other level means calibrating the complement of
# what the fit reports.
test_that("ps_calibrate() inverts a fit's scores for a named focal level", {
  fit <- calib_binary_fit()
  inverted <- 1 - predict(fit, type = "response")

  expect_same_calibration(
    ps_calibrate(fit, .focal_level = 0, smooth = FALSE),
    ps_calibrate(
      inverted,
      calib_model_data$z,
      .focal_level = 0,
      smooth = FALSE
    )
  )
  expect_false(isTRUE(all.equal(
    as.numeric(ps_calibrate(fit, .focal_level = 0, smooth = FALSE)),
    as.numeric(ps_calibrate(fit, smooth = FALSE))
  )))
})

test_that("ps_calibrate() inverts a two-level multinomial fit for a named level", {
  skip_if_not_installed("nnet")

  fit <- calib_two_level_fit()
  inverted <- 1 - as.numeric(fitted(fit))

  expect_same_calibration(
    ps_calibrate(fit, .focal_level = "control", smooth = FALSE),
    ps_calibrate(
      inverted,
      calib_model_data$a2,
      .focal_level = "control",
      smooth = FALSE
    )
  )
  expect_false(isTRUE(all.equal(
    as.numeric(ps_calibrate(fit, .focal_level = "control", smooth = FALSE)),
    as.numeric(ps_calibrate(fit, smooth = FALSE))
  )))
})

# The model route maps the deprecated spelling itself, before the scores reach
# the method that would map it again, so the reader is told about it once and
# under the name of the function they called.
test_that("the deprecated .treated on a fit is attributed to ps_calibrate()", {
  messages <- deprecation_warnings_from_user(
    quote(ps_calibrate(fit, .treated = 0, smooth = FALSE)),
    list(fit = calib_binary_fit())
  )

  expect_length(messages, 1)
  expect_match(messages[[1]], "ps_calibrate()", fixed = TRUE)
  expect_false(deprecation_misattributed(messages))
})
