# Test Manual Weight Calculations with ALL vctrs Classes
#
# This file comprehensively tests manual calculations for ALL vctrs classes:
# - psw: propensity score weights
# - ps_trim: trimmed propensity scores
# - ps_trunc: truncated propensity scores

library(testthat)

# =============================================================================
# Manual Calculations with ps_trim Class
# =============================================================================

test_that("All arithmetic operations work with ps_trim", {
  set.seed(123)
  ps <- runif(100, 0.1, 0.9)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)

  # Basic arithmetic operations
  expect_type(ps_trimmed + 0.1, "double")
  expect_type(ps_trimmed - 0.1, "double")
  expect_type(ps_trimmed * 2, "double")
  expect_type(ps_trimmed / 2, "double")
  expect_type(ps_trimmed^2, "double")

  # Reciprocals
  expect_type(1 / ps_trimmed, "double")
  expect_type(1 / (1 - ps_trimmed), "double")

  # Complex arithmetic
  expect_type(ps_trimmed / (1 - ps_trimmed), "double") # Odds
  expect_type(ps_trimmed * (1 - ps_trimmed), "double") # Variance function
  expect_type((1 - ps_trimmed) / ps_trimmed, "double") # Inverse odds

  # Mathematical functions
  expect_type(sqrt(ps_trimmed), "double")
  expect_type(log(ps_trimmed), "double")
  expect_type(exp(ps_trimmed), "double")
  expect_type(abs(ps_trimmed - 0.5), "double")

  # Trigonometric functions
  expect_type(sin(ps_trimmed), "double")
  expect_type(cos(ps_trimmed), "double")
  expect_type(tan(ps_trimmed), "double")

  # Logical operations
  expect_type(as.numeric(ps_trimmed) > 0.5, "logical")
  expect_type(as.numeric(ps_trimmed) <= 0.3, "logical")
  expect_type(as.numeric(ps_trimmed) == 0.5, "logical")
})

test_that("Manual weight calculations work with ps_trim", {
  set.seed(456)
  ps <- runif(50, 0.15, 0.85)
  exposure <- rbinom(50, 1, ps)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)

  # ATE weights
  ate_weights <- ifelse(exposure == 1, 1 / ps_trimmed, 1 / (1 - ps_trimmed))
  expect_type(ate_weights, "double")
  expect_false(inherits(ate_weights, "ps_trim"))

  # ATT weights
  att_weights <- ifelse(exposure == 1, 1, ps_trimmed / (1 - ps_trimmed))
  expect_type(att_weights, "double")

  # ATC weights
  atc_weights <- ifelse(exposure == 1, (1 - ps_trimmed) / ps_trimmed, 1)
  expect_type(atc_weights, "double")

  # ATO weights
  ato_weights <- ps_trimmed * (1 - ps_trimmed)
  expect_type(ato_weights, "double")

  # ATM weights
  atm_weights <- pmin(as.numeric(ps_trimmed), 1 - as.numeric(ps_trimmed)) /
    ifelse(exposure == 1, as.numeric(ps_trimmed), 1 - as.numeric(ps_trimmed))
  expect_type(atm_weights, "double")

  # Entropy weights
  entropy_weights <- -ps_trimmed *
    log(ps_trimmed) -
    (1 - ps_trimmed) * log(1 - ps_trimmed)
  expect_type(entropy_weights, "double")
})

test_that("Complex calculations work with ps_trim", {
  set.seed(789)
  ps <- runif(40, 0.2, 0.8)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7)

  # Logit transformation
  logit_ps <- log(ps_trimmed / (1 - ps_trimmed))
  expect_type(logit_ps, "double")

  # Back transformation
  inv_logit <- exp(logit_ps) / (1 + exp(logit_ps))
  expect_type(inv_logit, "double")

  # Probit transformation
  probit_ps <- qnorm(ps_trimmed)
  expect_type(probit_ps, "double")

  # Stabilized weights calculation
  mean_ps <- mean(ps_trimmed, na.rm = TRUE)
  stab_factor <- mean_ps / ps_trimmed
  expect_type(stab_factor, "double")

  # Effective sample size components
  weights_squared <- (1 / ps_trimmed)^2
  expect_type(weights_squared, "double")
})

# =============================================================================
# Manual Calculations with ps_trunc Class
# =============================================================================

test_that("All arithmetic operations work with ps_trunc", {
  set.seed(234)
  ps <- runif(100, 0.05, 0.95)
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)

  # Basic arithmetic operations
  expect_type(ps_truncated + 0.1, "double")
  expect_type(ps_truncated - 0.1, "double")
  expect_type(ps_truncated * 2, "double")
  expect_type(ps_truncated / 2, "double")
  expect_type(ps_truncated^2, "double")

  # Reciprocals (no NAs since truncation doesn't create them)
  recip <- 1 / ps_truncated
  expect_type(recip, "double")
  expect_true(all(is.finite(recip)))

  comp_recip <- 1 / (1 - ps_truncated)
  expect_type(comp_recip, "double")
  expect_true(all(is.finite(comp_recip)))

  # Complex arithmetic
  expect_type(ps_truncated / (1 - ps_truncated), "double")
  expect_type(ps_truncated * (1 - ps_truncated), "double")
  expect_type((1 - ps_truncated) / ps_truncated, "double")

  # Mathematical functions
  expect_type(sqrt(ps_truncated), "double")
  expect_type(log(ps_truncated), "double")
  expect_type(exp(ps_truncated), "double")
  expect_type(abs(ps_truncated - 0.5), "double")
})

test_that("Manual weight calculations work with ps_trunc", {
  set.seed(567)
  ps <- runif(60, 0.05, 0.95)
  exposure <- rbinom(60, 1, ps)
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)

  # ATE weights - should have no NAs
  ate_weights <- ifelse(exposure == 1, 1 / ps_truncated, 1 / (1 - ps_truncated))
  expect_type(ate_weights, "double")
  expect_false(inherits(ate_weights, "ps_trunc"))
  expect_equal(sum(is.na(ate_weights)), 0)

  # ATT weights
  att_weights <- ifelse(exposure == 1, 1, ps_truncated / (1 - ps_truncated))
  expect_type(att_weights, "double")
  expect_equal(sum(is.na(att_weights)), 0)

  # ATC weights
  atc_weights <- ifelse(exposure == 1, (1 - ps_truncated) / ps_truncated, 1)
  expect_type(atc_weights, "double")
  expect_equal(sum(is.na(atc_weights)), 0)

  # ATO weights
  ato_weights <- ps_truncated * (1 - ps_truncated)
  expect_type(ato_weights, "double")
  expect_true(all(ato_weights >= 0))
  expect_true(all(ato_weights <= 0.25))

  # ATM weights
  atm_weights <- pmin(as.numeric(ps_truncated), 1 - as.numeric(ps_truncated)) /
    ifelse(exposure == 1, as.numeric(ps_truncated), 1 - as.numeric(ps_truncated))
  expect_type(atm_weights, "double")

  # All weights should be finite
  expect_true(all(is.finite(ate_weights)))
  expect_true(all(is.finite(att_weights)))
  expect_true(all(is.finite(atc_weights)))
})

test_that("Complex calculations work with ps_trunc", {
  set.seed(890)
  ps <- runif(45, 0.1, 0.9)
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.15, upper = 0.85)

  # Bounded transformations
  bounded_logit <- log((ps_truncated - 0.15) / (0.85 - ps_truncated))
  expect_type(bounded_logit, "double")
  # Some values may be infinite if ps_truncated is exactly at bounds
  expect_true(any(is.finite(bounded_logit)))

  # Variance stabilizing transformation
  arcsin_sqrt <- asin(sqrt(ps_truncated))
  expect_type(arcsin_sqrt, "double")

  # Beta distribution parameters
  alpha_param <- ps_truncated * 10
  beta_param <- (1 - ps_truncated) * 10
  expect_type(alpha_param, "double")
  expect_type(beta_param, "double")

  # Entropy calculation
  entropy <- -ps_truncated *
    log(ps_truncated) -
    (1 - ps_truncated) * log(1 - ps_truncated)
  expect_type(entropy, "double")
  expect_true(all(entropy >= 0))
})

# =============================================================================
# Manual Calculations with psw Class
# =============================================================================

test_that("All arithmetic operations work with psw", {
  set.seed(345)
  weights <- psw(runif(80, 0.5, 2.5), estimand = "ate", stabilized = TRUE)

  # Basic arithmetic operations (psw preserves class)
  expect_s3_class(weights + 0.5, "psw")
  expect_s3_class(weights - 0.5, "psw")
  expect_s3_class(weights * 2, "psw")
  expect_s3_class(weights / 2, "psw")
  expect_s3_class(weights^2, "psw")

  # Check attributes are preserved
  doubled <- weights * 2
  expect_equal(estimand(doubled), "ate")
  expect_true(is_stabilized(doubled))

  # Reciprocals preserve class
  expect_s3_class(1 / weights, "psw")
  expect_s3_class(2 / weights, "psw")

  # Mathematical functions return numeric
  expect_type(sqrt(weights), "double")
  expect_type(log(weights), "double")
  expect_type(exp(weights), "double")
  expect_type(abs(weights - 1), "double")

  # Trigonometric functions return numeric
  expect_type(sin(weights), "double")
  expect_type(cos(weights), "double")
  expect_type(tan(weights), "double")
})

test_that("psw weight manipulations work correctly", {
  set.seed(678)
  n <- 70
  ps <- runif(n, 0.2, 0.8)
  exposure <- rbinom(n, 1, ps)

  # Create initial weights
  ate_weights <- ifelse(exposure == 1, 1 / ps, 1 / (1 - ps))
  weights_obj <- psw(ate_weights, estimand = "ate")

  # Trimming weights by value
  trimmed_weights <- weights_obj
  trimmed_weights[as.numeric(weights_obj) > 10] <- 10 # Cap at 10
  expect_s3_class(trimmed_weights, "psw")

  # Normalizing weights
  normalized_weights <- weights_obj / sum(weights_obj)
  expect_s3_class(normalized_weights, "psw")
  expect_equal(sum(normalized_weights), 1, tolerance = 1e-10)

  # Stabilizing weights
  stab_factor <- mean(exposure) / weights_obj
  expect_s3_class(stab_factor, "psw")

  # Weight diagnostics
  cv_weights <- sd(weights_obj) / mean(weights_obj)
  expect_type(cv_weights, "double")

  # Effective sample size
  ess <- sum(weights_obj)^2 / sum(weights_obj^2)
  expect_type(ess, "double")
})

test_that("Complex calculations work with psw", {
  set.seed(901)
  weights <- psw(runif(50, 0.8, 1.2), estimand = "att", stabilized = FALSE)

  # Weight ratio calculations
  weight_ratio <- weights / mean(weights)
  expect_s3_class(weight_ratio, "psw")

  # Log weight calculations (returns numeric)
  log_weights <- log(weights)
  expect_type(log_weights, "double")

  # Exponential tilting
  tilted <- exp(0.5 * log_weights) # sqrt via exp-log
  expect_type(tilted, "double")

  # Weight percentiles via sorting
  sorted_weights <- sort(weights)
  expect_s3_class(sorted_weights, "psw")

  # Cumulative weights
  cumsum_weights <- cumsum(weights)
  expect_s3_class(cumsum_weights, "psw")

  # Running products
  cumprod_weights <- cumprod(weights)
  expect_s3_class(cumprod_weights, "psw")
})

# =============================================================================
# Summary Functions Return Numeric for ALL Classes
# =============================================================================

test_that("Summary functions return numeric for ps_trim", {
  set.seed(111)
  ps_trimmed <- ps_trim(
    runif(30, 0.2, 0.8),
    method = "ps",
    lower = 0.3,
    upper = 0.7
  )

  # All summary functions return numeric
  expect_type(min(ps_trimmed, na.rm = TRUE), "double")
  expect_type(max(ps_trimmed, na.rm = TRUE), "double")
  expect_type(range(ps_trimmed, na.rm = TRUE), "double")
  expect_type(sum(ps_trimmed, na.rm = TRUE), "double")
  expect_type(prod(ps_trimmed, na.rm = TRUE), "double")
  expect_type(mean(ps_trimmed, na.rm = TRUE), "double")
  expect_type(median(ps_trimmed, na.rm = TRUE), "double")

  # Should work in calculations
  expect_type(1 / max(ps_trimmed, na.rm = TRUE), "double")
  expect_type(min(ps_trimmed, na.rm = TRUE) * 100, "double")

  # Multiple arguments
  ps_trimmed2 <- ps_trim(
    runif(30, 0.2, 0.8),
    method = "ps",
    lower = 0.3,
    upper = 0.7
  )
  expect_type(min(ps_trimmed, ps_trimmed2, na.rm = TRUE), "double")
  expect_type(max(ps_trimmed, ps_trimmed2, na.rm = TRUE), "double")
})

test_that("Summary functions return numeric for ps_trunc", {
  set.seed(222)
  ps_truncated <- ps_trunc(
    runif(35, 0.1, 0.9),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  # All summary functions return numeric
  expect_type(min(ps_truncated), "double")
  expect_type(max(ps_truncated), "double")
  expect_type(range(ps_truncated), "double")
  expect_type(sum(ps_truncated), "double")
  expect_type(prod(ps_truncated), "double")
  expect_type(mean(ps_truncated), "double")
  expect_type(median(ps_truncated), "double")

  # No NA handling needed for truncation
  expect_equal(min(ps_truncated), min(ps_truncated, na.rm = TRUE))

  # Quantile functions
  expect_type(quantile(ps_truncated, 0.5), "double")
  expect_type(quantile(ps_truncated, c(0.25, 0.75)), "double")
})

test_that("Summary functions return numeric for psw", {
  set.seed(333)
  weights <- psw(runif(40, 0.5, 2.0), estimand = "ato")

  # All summary functions return numeric
  expect_type(min(weights), "double")
  expect_type(max(weights), "double")
  expect_type(range(weights), "double")
  expect_type(sum(weights), "double")
  expect_type(prod(weights), "double")
  expect_type(mean(weights), "double")
  expect_type(median(weights), "double")

  # Summary statistics for weight diagnostics
  expect_type(var(weights), "double")
  expect_type(sd(weights), "double")
  expect_type(IQR(weights), "double")
  expect_type(mad(weights), "double")
})

# =============================================================================
# vec_restore Works Correctly for ALL Classes
# =============================================================================

test_that("Subsetting preserves class and attributes for all classes", {
  set.seed(444)

  # ps_trim subsetting
  ps_trim_obj <- ps_trim(runif(50, 0.1, 0.9), method = "adaptive")
  subset1 <- ps_trim_obj[1:25]
  subset2 <- ps_trim_obj[seq(1, 50, by = 2)] # Every other element
  subset3 <- ps_trim_obj[as.numeric(ps_trim_obj) > 0.5]

  expect_s3_class(subset1, "ps_trim")
  expect_s3_class(subset2, "ps_trim")
  expect_s3_class(subset3, "ps_trim")
  expect_equal(ps_trim_meta(subset1)$method, "adaptive")

  # ps_trunc subsetting
  ps_trunc_obj <- ps_trunc(
    runif(50, 0.1, 0.9),
    method = "pctl",
    lower = 0.05,
    upper = 0.95
  )
  subset1 <- ps_trunc_obj[1:25]
  subset2 <- ps_trunc_obj[seq(1, 50, by = 2)]
  subset3 <- ps_trunc_obj[as.numeric(ps_trunc_obj) > 0.5]

  expect_s3_class(subset1, "ps_trunc")
  expect_s3_class(subset2, "ps_trunc")
  expect_s3_class(subset3, "ps_trunc")
  expect_equal(ps_trunc_meta(subset1)$method, "pctl")

  # psw subsetting
  psw_obj <- psw(
    runif(50, 0.5, 2.0),
    estimand = "atm",
    stabilized = TRUE,
    calibrated = TRUE
  )
  subset1 <- psw_obj[1:25]
  subset2 <- psw_obj[seq(1, 50, by = 2)]
  subset3 <- psw_obj[as.numeric(psw_obj) > 1]

  expect_s3_class(subset1, "psw")
  expect_s3_class(subset2, "psw")
  expect_s3_class(subset3, "psw")
  expect_equal(estimand(subset1), "atm")
  expect_true(is_stabilized(subset1))
  expect_true(attr(subset1, "calibrated"))
})

# =============================================================================
# Edge Cases and Special Values
# =============================================================================

test_that("Edge cases work correctly for all classes", {
  set.seed(555)

  # Empty vectors
  empty_trim <- ps_trim(numeric(0), method = "ps")
  empty_trunc <- ps_trunc(numeric(0), method = "ps")
  empty_psw <- psw(numeric(0), estimand = "ate")

  expect_length(empty_trim, 0)
  expect_length(empty_trunc, 0)
  expect_length(empty_psw, 0)
  expect_type(sum(empty_trim), "double")
  expect_type(sum(empty_trunc), "double")
  expect_type(sum(empty_psw), "double")

  # Single element
  single_trim <- ps_trim(0.5, method = "ps", lower = 0.1, upper = 0.9)
  single_trunc <- ps_trunc(0.5, method = "ps", lower = 0.1, upper = 0.9)
  single_psw <- psw(1.5, estimand = "ate")

  expect_type(1 / single_trim, "double")
  expect_type(1 / single_trunc, "double")
  expect_s3_class(1 / single_psw, "psw")

  # All trimmed case
  all_trimmed <- ps_trim(rep(0.05, 10), method = "ps", lower = 0.1, upper = 0.9)
  weights_from_trimmed <- 1 / all_trimmed
  expect_true(all(is.na(weights_from_trimmed)))

  # Boundary values for truncation
  boundary_ps <- c(0.001, 0.999)
  boundary_trunc <- ps_trunc(
    boundary_ps,
    method = "ps",
    lower = 0.01,
    upper = 0.99
  )
  boundary_weights <- 1 / boundary_trunc
  expect_true(all(is.finite(boundary_weights)))
  expect_true(all(boundary_weights > 0))
})
