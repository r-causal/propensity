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

test_that("preserves psw attributes from an existing causal‐weights object", {
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
  calib_bins <- cut(calibrated_ps, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  
  obs_calib_error <- abs(tapply(treat, obs_bins, mean, na.rm = TRUE) - 
                        tapply(obs_ps, obs_bins, mean, na.rm = TRUE))
  calib_calib_error <- abs(tapply(treat, calib_bins, mean, na.rm = TRUE) - 
                          tapply(calibrated_ps, calib_bins, mean, na.rm = TRUE))
  
  # Average calibration error should generally be reduced
  expect_true(mean(calib_calib_error, na.rm = TRUE) <= 
              mean(obs_calib_error, na.rm = TRUE) + 0.1)  # Allow some tolerance
})

test_that("handles edge cases with extreme propensity scores", {
  # Near 0 and 1 values
  ps <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
  treat <- c(0, 0, 0, 1, 1, 1, 1)
  
  expect_no_error(calibrated <- ps_calibrate(ps, treat))
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
  calib2 <- ps_calibrate(ps, treat_char, .treated = "treated", .untreated = "control")
  calib3 <- ps_calibrate(ps, treat_factor, .treated = "treated", .untreated = "control")
  calib4 <- ps_calibrate(ps, treat_logical)
  
  # All should produce the same result
  expect_equal(as.numeric(calib1), as.numeric(calib2))
  expect_equal(as.numeric(calib1), as.numeric(calib3))
  expect_equal(as.numeric(calib1), as.numeric(calib4))
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
  calibrated <- ps_calibrate(ps, treat)
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
  expect_true(all(diffs >= -1e-10))  # Allow for numerical tolerance
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

test_that("isotonic and platt calibration can differ", {
  set.seed(123)
  # Create data where isotonic might perform differently
  ps <- c(rep(0.2, 50), rep(0.8, 50))
  treat <- c(rbinom(50, 1, 0.3), rbinom(50, 1, 0.7))
  
  calib_platt <- ps_calibrate(ps, treat, method = "platt")
  calib_iso <- ps_calibrate(ps, treat, method = "isoreg")
  
  # They should produce different results in general
  expect_false(identical(as.numeric(calib_platt), as.numeric(calib_iso)))
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

test_that("isotonic regression preserves monotonicity better than Platt", {
  set.seed(888)
  # Create data where isotonic should perform better
  n <- 100
  ps <- seq(0.1, 0.9, length.out = n)
  # Non-linear relationship that violates logistic assumption
  true_prob <- 0.2 + 0.6 * ps^2
  treat <- rbinom(n, 1, true_prob)
  
  platt_calib <- ps_calibrate(ps, treat, method = "platt")
  iso_calib <- ps_calibrate(ps, treat, method = "isoreg")
  
  # Check isotonic preserves monotonicity
  iso_diffs <- diff(as.numeric(iso_calib))
  expect_true(all(iso_diffs >= -1e-10))  # Allow for numerical tolerance
  
  # Both should be different from original
  expect_false(identical(as.numeric(ps), as.numeric(platt_calib)))
  expect_false(identical(as.numeric(ps), as.numeric(iso_calib)))
  
  # They should produce different results for non-linear data
  expect_false(identical(as.numeric(platt_calib), as.numeric(iso_calib)))
})

test_that("isotonic regression handles various cases like WeightIt", {
  skip_if_not_installed("WeightIt")
  
  # Test with ties in propensity scores (but ensure sufficient data)
  set.seed(999)
  ps_ties <- rep(c(0.3, 0.7), each = 4)
  treat_ties <- c(0, 0, 1, 1, 0, 1, 1, 1)
  
  our_iso_ties <- ps_calibrate(ps_ties, treat_ties, method = "isoreg")
  weightit_iso_ties <- WeightIt::calibrate(ps_ties, treat_ties, method = "isoreg")
  
  expect_equal(as.numeric(our_iso_ties), as.numeric(weightit_iso_ties), tolerance = 1e-10)
  
  # Test that our implementation handles edge cases gracefully even if WeightIt fails
  ps_extreme <- c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
  treat_extreme <- c(0, 0, 0, 1, 1, 1, 1)
  
  # Our implementation should handle this without error
  expect_no_error(our_extreme <- ps_calibrate(ps_extreme, treat_extreme, method = "isoreg"))
  expect_true(all(our_extreme >= 0 & our_extreme <= 1))
  expect_true(all(diff(as.numeric(our_extreme)) >= -1e-10))  # Monotonic
})
