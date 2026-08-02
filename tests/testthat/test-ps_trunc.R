test_that("ps_trunc() - PS method uses fixed bounds", {
  set.seed(123)
  ps <- c(0.01, 0.1, 0.3, 0.8, 0.95)
  out <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)

  expect_s3_class(out, "ps_trunc")
  out_data <- as.numeric(out)

  # Values <0.2 => 0.2, >0.8 => 0.8
  expect_equal(out_data[1], 0.2) # was 0.01
  expect_equal(out_data[2], 0.2) # was 0.1
  expect_equal(out_data[3], 0.3) # stays 0.3
  expect_equal(out_data[4], 0.8) # stays 0.8
  expect_equal(out_data[5], 0.8) # was 0.95 => truncated
})

test_that("ps_trunc() - pctl method uses quantiles", {
  set.seed(1)
  n <- 50
  ps <- plogis(rnorm(n, 0, 1.2))

  # default [0.05, 0.95]
  out1 <- ps_trunc(ps, method = "pctl")
  meta1 <- ps_trunc_meta(out1)
  expect_equal(meta1$lower_pctl, 0.05)
  expect_equal(meta1$upper_pctl, 0.95)
  out1_data <- as.numeric(out1)

  # Check boundary
  q_l <- quantile(ps, probs = 0.05)
  q_u <- quantile(ps, probs = 0.95)
  expect_true(all(out1_data >= q_l - 1e-8))
  expect_true(all(out1_data <= q_u + 1e-8))

  # custom [0.2, 0.8]
  out2 <- ps_trunc(ps, method = "pctl", lower = 0.2, upper = 0.8)
  meta2 <- ps_trunc_meta(out2)
  expect_equal(meta2$lower_pctl, 0.2)
  expect_equal(meta2$upper_pctl, 0.8)
  out2_data <- as.numeric(out2)
  q_l2 <- quantile(ps, probs = 0.2)
  q_u2 <- quantile(ps, probs = 0.8)

  # Everything below q_l2 => replaced with q_l2
  expect_true(all(out2_data >= q_l2 - 1e-8))
  expect_true(all(out2_data <= q_u2 + 1e-8))
})

test_that("ps_trunc() - cr method uses min(ps_treat)/max(ps_untrt)", {
  set.seed(2)
  n <- 30
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  out_cr <- ps_trunc(ps, .exposure = z, method = "cr", .focal_level = 1)
  expect_s3_class(out_cr, "ps_trunc")

  meta_cr <- ps_trunc_meta(out_cr)
  ps_treat <- ps[z == 1]
  ps_untrt <- ps[z == 0]
  cr_lower <- min(ps_treat)
  cr_upper <- max(ps_untrt)
  expect_equal(meta_cr$lower_bound, cr_lower)
  expect_equal(meta_cr$upper_bound, cr_upper)

  # check bounding
  out_data <- as.numeric(out_cr)
  expect_true(all(out_data >= cr_lower - 1e-8))
  expect_true(all(out_data <= cr_upper + 1e-8))
})

test_that("ps_trunc() errors on invalid usage or .exposure", {
  # if method="cr" but no .exposure => error
  expect_propensity_error(
    ps_trunc(runif(10), method = "cr")
  )

  # if .exposure not 0/1 => error
  expect_propensity_error(
    ps_trunc(runif(5), .exposure = 1:5, method = "cr")
  )

  # if lower >= upper => error for method="ps"
  expect_propensity_error(
    ps_trunc(runif(5), method = "ps", lower = 0.8, upper = 0.3)
  )
})

test_that("Truncation workflow yields truncated psw with no refit logic", {
  set.seed(888)
  n <- 10
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))

  # 1) Fit logistic model
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 2) Truncate (winsorize) the PS
  truncated_ps <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  expect_s3_class(truncated_ps, "ps_trunc")

  # 3) Compute ATE weights
  w_ate <- wt_ate(
    truncated_ps,
    .exposure = z,
    exposure_type = "binary",
    .focal_level = 1
  )
  expect_s3_class(w_ate, "psw")

  # 4) Verify truncated, not trimmed, not refit, estimand
  expect_true(is_ps_truncated(w_ate))
  expect_false(is_ps_trimmed(w_ate))
  expect_false(is_refit(w_ate))
  expect_match(estimand(w_ate), "; truncated$")
})

test_that("is_ps_truncated.default() -> FALSE, is_ps_truncated.ps_trunc() -> TRUE", {
  # 1) A plain numeric => default => FALSE
  expect_false(is_ps_truncated(runif(5)))

  # 2) A simple ps_trunc object => is_ps_truncated(...) => TRUE
  # Create via new_ps_trunc()
  my_trunc <- new_ps_trunc(
    x = c(0.2, 0.6, 0.8),
    meta = list(method = "ps", lower_bound = 0.2, upper_bound = 0.8)
  )
  expect_true(is_ps_truncated(my_trunc))
})

test_that("Arithmetic on ps_trunc returns numeric", {
  obj <- new_ps_trunc(
    x = c(0.2, 0.7, 0.9),
    meta = list(method = "ps", lower_bound = 0.2, upper_bound = 0.8)
  )

  # Arithmetic operations should return numeric
  expect_type(obj + 1, "double")
  expect_type(1 + obj, "double")
  expect_type(obj * 2, "double")
  expect_type(1 / obj, "double")

  # Verify values are correct
  expect_equal(obj + 1, c(1.2, 1.7, 1.9))
  expect_equal(1 / obj, c(5, 10 / 7, 10 / 9))

  # Combining two ps_trunc also returns numeric
  obj2 <- new_ps_trunc(
    x = c(0.1, 0.1, 0.3),
    meta = list(method = "ps", lower_bound = 0.1, upper_bound = 0.5)
  )
  expect_type(obj * obj2, "double")
  expect_equal(obj * obj2, c(0.02, 0.07, 0.27))
})

test_that("Combining & casting ps_trunc => correct ptype2, cast behavior", {
  obj <- new_ps_trunc(
    x = c(0.2, 0.6, 0.8),
    meta = list(method = "ps", lower_bound = 0.2, upper_bound = 0.8)
  )
  # 1) Combining two ps_trunc => error
  obj2 <- new_ps_trunc(
    x = c(0.4, 0.5, 0.7),
    meta = list(method = "ps", lower_bound = 0.3, upper_bound = 0.8)
  )

  # 3) Casting ps_trunc -> double => numeric data
  out_cast <- vctrs::vec_cast(obj, double())
  expect_type(out_cast, "double")
  expect_identical(out_cast, c(0.2, 0.6, 0.8))

  # 4) Casting double -> ps_trunc => new default meta
  new_vals <- runif(3)
  out_ps_trunc <- vctrs::vec_cast(new_vals, to = obj)
  expect_s3_class(out_ps_trunc, "ps_trunc")
  meta_new <- ps_trunc_meta(out_ps_trunc)
  expect_equal(meta_new$method, "unknown") # per your code
  expect_true(is.na(meta_new$lower_bound))
  expect_true(is.na(meta_new$upper_bound))
})

test_that("wt_atm.numeric calls atm_binary() for binary .exposure, returns psw", {
  set.seed(101)
  n <- 8
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x))

  # A numeric PS
  ps <- plogis(0.4 * x)

  # 1) Binary .exposure => calls atm_binary() => returns psw
  out_atm <- wt_atm.numeric(
    .propensity = ps,
    .exposure = z,
    exposure_type = "binary",
    .focal_level = 1
  )
  # Check it's a psw object with estimand "atm"
  expect_s3_class(out_atm, "psw")
  expect_equal(estimand(out_atm), "atm")
})

test_that("atm_binary() logic with transform_.exposure_binary() is triggered", {
  # atm_binary => pmin(ps, 1-ps) / (.exposure*ps + (1-.exposure)*(1-ps))
  ps_vec <- c(0.2, 0.8, 0.5)
  z_vec <- c(0, 1, 1)

  w <- atm_binary(
    .propensity = ps_vec,
    .exposure = z_vec,
    .focal_level = 1
  )
  # Just check dimension, no error
  expect_length(w, 3)

  # If .exposure isn't 0/1 or has different factor levels, transform_.exposure_binary
  # is tested. We'll do a quick check with factor( c("C","T","T") )
  w2 <- atm_binary(
    .propensity = ps_vec,
    .exposure = factor(c("C", "T", "T")),
    .focal_level = "T",
    .reference_level = "C"
  )
  expect_length(w2, 3)
})

test_that("wt_ato.numeric calls ato_binary() for binary .exposure, returns psw", {
  set.seed(202)
  n <- 6
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.6 * x))
  ps <- plogis(0.3 * x)

  # 1) Binary => calls ato_binary => returns psw
  out_ato <- wt_ato.numeric(
    .propensity = ps,
    .exposure = z,
    exposure_type = "binary",
    .focal_level = 1
  )
  expect_s3_class(out_ato, "psw")
  expect_equal(estimand(out_ato), "ato")
})

test_that("ato_binary() logic is triggered for p=0.3", {
  # (1 - p)*.exposure + p*(1-.exposure)
  ps_vec <- c(0.1, 0.9, 0.5)
  z_vec <- c(0, 1, 1)

  w <- ato_binary(
    .propensity = ps_vec,
    .exposure = z_vec,
    .focal_level = 1
  )
  expect_length(w, 3)
  # Just check no error, correct length
})

test_that("wt_atm.ps_trunc synergy with truncated object yields truncated psw", {
  set.seed(303)
  n <- 6
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 * x))

  ps <- plogis(0.7 * x)
  # Make a truncated object (like bounding ps in [0.2,0.8])
  trunc_obj <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)

  # Now call wt_atm() on ps_trunc => dispatches wt_atm.ps_trunc()
  w_atm <- wt_atm(
    trunc_obj,
    .exposure = z,
    exposure_type = "binary",
    .focal_level = 1
  )
  expect_s3_class(w_atm, "psw")
  # Estimand => "atm; truncated"
  expect_match(estimand(w_atm), "atm; truncated")
  # truncated=TRUE
  expect_true(is_ps_truncated(w_atm))
})

test_that("wt_ato.ps_trunc synergy with truncated object yields truncated psw", {
  set.seed(404)
  n <- 7
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.1 * x))
  ps <- plogis(0.5 * x)

  # bounding p in [0.1, 0.9], e.g.
  trunc_obj <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)
  w_ato <- wt_ato(
    trunc_obj,
    .exposure = z,
    exposure_type = "binary",
    .focal_level = 1
  )

  expect_s3_class(w_ato, "psw")
  expect_match(estimand(w_ato), "ato; truncated")
  expect_true(is_ps_truncated(w_ato))
})

test_that("is_unit_truncated.ps_trunc returns expected row-level booleans", {
  set.seed(101)
  ps_vec <- c(0.1, 0.2, 0.5, 0.85, 0.95)

  # Truncate outside [0.2, 0.8]
  truncated_obj <- ps_trunc(
    ps_vec,
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  expect_s3_class(truncated_obj, "ps_trunc")

  row_trunc <- is_unit_truncated(truncated_obj)
  expect_type(row_trunc, "logical")
  expect_length(row_trunc, length(ps_vec))

  truncated_data <- as.numeric(truncated_obj)
  expect_equal(which(row_trunc), c(1, 4, 5))
  expect_equal(truncated_data, c(0.2, 0.2, 0.5, 0.8, 0.8))
})

test_that("ps_trunc objects can convert to character", {
  ps <- c(0.01, 0.1, 0.3, 0.8, 0.95)
  out <- as.character(ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8))
  expect_type(out, "character")
})

test_that("ps_trunc works with summarize(mean = mean(ps))", {
  skip_if_not_installed("dplyr")
  library(dplyr, warn.conflicts = FALSE)

  set.seed(200)
  n <- 600
  x <- rnorm(n)
  z <- rbinom(n, size = 1, prob = plogis(x + rnorm(n)))
  fit <- glm(z ~ x, family = binomial)

  ps <- predict(fit, type = "response") |>
    ps_trunc(method = "ps", lower = 0.3, upper = 0.7)

  out <- tibble(x, z, ps) |>
    group_by(truncated = is_unit_truncated(ps)) |>
    summarize(mean = mean(ps), .groups = "drop")

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("truncated", "mean"))
  expect_type(out$mean, "double")
})

test_that("ps_trunc vec_ptype_full output matches expected format", {
  set.seed(123)
  ps <- runif(20, 0.05, 0.95)

  # Create ps_trunc with some values truncated
  ps_trunc_obj <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)

  # Test the vec_ptype_full output - should show bounds and method
  expect_equal(
    vctrs::vec_ptype_full(ps_trunc_obj),
    "ps_trunc{[0.2,0.8], method=ps}"
  )

  # Test with different bounds
  ps_trunc_narrow <- ps_trunc(ps, method = "ps", lower = 0.4, upper = 0.6)
  expect_equal(
    vctrs::vec_ptype_full(ps_trunc_narrow),
    "ps_trunc{[0.4,0.6], method=ps}"
  )

  # Test with very wide bounds (no actual truncation)
  ps_trunc_wide <- ps_trunc(ps, method = "ps", lower = 0.01, upper = 0.99)
  expect_equal(
    vctrs::vec_ptype_full(ps_trunc_wide),
    "ps_trunc{[0.01,0.99], method=ps}"
  )
})

test_that("ps_trunc index tracking works when combining objects", {
  set.seed(456)
  ps1 <- runif(10, 0.05, 0.95)
  ps2 <- runif(10, 0.05, 0.95)

  # Create ps_trunc objects with same parameters
  ps_trunc1 <- ps_trunc(ps1, method = "ps", lower = 0.2, upper = 0.8)
  ps_trunc2 <- ps_trunc(ps2, method = "ps", lower = 0.2, upper = 0.8)

  # Get original truncated indices
  meta1 <- ps_trunc_meta(ps_trunc1)
  meta2 <- ps_trunc_meta(ps_trunc2)
  n_truncated1 <- length(meta1$truncated_idx)
  n_truncated2 <- length(meta2$truncated_idx)

  # Combine the objects
  combined <- c(ps_trunc1, ps_trunc2)

  # Should be a ps_trunc object
  expect_s3_class(combined, "ps_trunc")

  # Check that indices are properly tracked
  combined_meta <- ps_trunc_meta(combined)
  expect_equal(length(combined), 20)

  # The total number of truncated should be the sum
  expect_equal(
    length(combined_meta$truncated_idx),
    n_truncated1 + n_truncated2
  )

  # Check that values at bounds are at the correct positions
  combined_data <- vec_data(combined)
  lower_bound <- combined_meta$lower_bound
  upper_bound <- combined_meta$upper_bound

  # All truncated indices should have values at the bounds
  truncated_values <- combined_data[combined_meta$truncated_idx]
  expect_true(all(
    truncated_values == lower_bound | truncated_values == upper_bound
  ))
})

test_that("ps_trunc warns when combining objects with different parameters", {
  ps1 <- runif(10, 0.05, 0.95)
  ps2 <- runif(10, 0.05, 0.95)

  # Create ps_trunc objects with different parameters
  ps_trunc1 <- ps_trunc(ps1, method = "ps", lower = 0.2, upper = 0.8)
  ps_trunc2 <- ps_trunc(ps2, method = "ps", lower = 0.3, upper = 0.7)

  # Should warn and return numeric
  expect_propensity_warning(
    # jarl-ignore implicit_assignment: assignment keeps the return value out of the snapshot while the warning is captured
    combined <- c(ps_trunc1, ps_trunc2)
  )

  expect_type(combined, "double")
  expect_false(inherits(combined, "ps_trunc"))
})

test_that("ps_trunc index tracking works with subsetting and combining", {
  set.seed(789)
  ps <- runif(20, 0.05, 0.95)

  # Create ps_trunc object
  ps_trunc_obj <- ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7)
  meta <- ps_trunc_meta(ps_trunc_obj)

  # Subset the object
  subset1 <- ps_trunc_obj[1:10]
  subset2 <- ps_trunc_obj[11:20]

  # Recombine
  recombined <- c(subset1, subset2)

  # Should maintain ps_trunc class
  expect_s3_class(recombined, "ps_trunc")

  # Check indices are properly tracked
  recombined_meta <- ps_trunc_meta(recombined)
  expect_equal(
    length(recombined_meta$truncated_idx),
    length(meta$truncated_idx)
  )

  # Check that truncated values are preserved at correct positions
  recombined_data <- vec_data(recombined)
  original_data <- vec_data(ps_trunc_obj)

  # Find which values were at the bounds
  lower_bound <- meta$lower_bound
  upper_bound <- meta$upper_bound
  original_at_bounds <- which(
    original_data == lower_bound | original_data == upper_bound
  )
  recombined_at_bounds <- which(
    recombined_data == lower_bound | recombined_data == upper_bound
  )

  expect_equal(recombined_at_bounds, original_at_bounds)
})

test_that("ps_trunc handles multiple combines correctly", {
  set.seed(321)

  # Create three ps_trunc objects
  ps1 <- runif(5, 0.05, 0.95)
  ps2 <- runif(5, 0.05, 0.95)
  ps3 <- runif(5, 0.05, 0.95)

  ps_trunc1 <- ps_trunc(ps1, method = "ps", lower = 0.25, upper = 0.75)
  ps_trunc2 <- ps_trunc(ps2, method = "ps", lower = 0.25, upper = 0.75)
  ps_trunc3 <- ps_trunc(ps3, method = "ps", lower = 0.25, upper = 0.75)

  # Combine all three
  combined <- c(ps_trunc1, ps_trunc2, ps_trunc3)

  # Should maintain ps_trunc class
  expect_s3_class(combined, "ps_trunc")
  expect_equal(length(combined), 15)

  # Check indices
  combined_meta <- ps_trunc_meta(combined)
  combined_data <- vec_data(combined)

  # All truncated indices should have values at bounds
  truncated_values <- combined_data[combined_meta$truncated_idx]
  expect_true(all(
    truncated_values == combined_meta$lower_bound |
      truncated_values == combined_meta$upper_bound
  ))
})

# Truncation record integrity ----------------------------------------------

# The record a `ps_trunc` carries names positions among the observations it was
# written for. An operation that changes how many observations there are is not
# handed the subscript, so it cannot re-index the record onto the result and
# drops it rather than leave positions describing units they were never about.

trunc_record_fixture <- function() {
  # Position 1 is pinned up to the lower bound. Position 4 already sits on the
  # upper bound and was left alone, so bound equality does not mark a unit.
  ps_trunc(
    c(0.05, 0.3, 0.5, 0.9, 0.6),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )
}

test_that("ps_trunc() records how many observations its positions describe", {
  meta <- ps_trunc_meta(trunc_record_fixture())

  expect_equal(meta$truncated_idx, 1L)
  expect_equal(meta$n_obs, 5L)
})

test_that("slicing a ps_trunc shorter drops the truncation record with a warning", {
  x <- trunc_record_fixture()

  cnd <- expect_warning(
    sliced <- vec_slice(x, 2:3),
    class = "propensity_trunc_record_warning"
  )
  expect_s3_class(cnd, "propensity_warning")

  expect_s3_class(sliced, "ps_trunc")
  expect_equal(as.numeric(sliced), c(0.3, 0.5))

  meta <- ps_trunc_meta(sliced)
  expect_null(meta$truncated_idx)
  expect_null(meta$n_obs)

  # Nothing that is not indexed by observation is touched by the drop.
  expect_equal(meta$method, "ps")
  expect_equal(meta$lower_bound, 0.1)
  expect_equal(meta$upper_bound, 0.9)
  expect_true(is_ps_truncated(sliced))

  expect_error(
    is_unit_truncated(sliced),
    class = "propensity_missing_meta_error"
  )
})

test_that("filtering a ps_trunc column drops the truncation record with a warning", {
  skip_if_not_installed("dplyr")

  df <- data.frame(id = 1:5)
  df$ps <- trunc_record_fixture()

  expect_warning(
    filtered <- dplyr::filter(df, id %in% 2:3),
    class = "propensity_trunc_record_warning"
  )

  expect_s3_class(filtered$ps, "ps_trunc")
  expect_equal(as.numeric(filtered$ps), c(0.3, 0.5))
  expect_null(ps_trunc_meta(filtered$ps)$truncated_idx)
  expect_error(
    is_unit_truncated(filtered$ps),
    class = "propensity_missing_meta_error"
  )
})

test_that("a length-preserving ps_trunc restore keeps the truncation record", {
  x <- trunc_record_fixture()
  meta <- ps_trunc_meta(x)
  truncated_units <- c(TRUE, FALSE, FALSE, FALSE, FALSE)

  whole <- expect_silent(vec_slice(x, seq_along(x)))
  expect_identical(ps_trunc_meta(whole), meta)
  expect_identical(is_unit_truncated(whole), truncated_units)

  empty_subscript <- expect_silent(x[])
  expect_identical(ps_trunc_meta(empty_subscript), meta)

  expect_silent({
    x[2] <- 0.4
  })
  expect_identical(ps_trunc_meta(x), meta)
  expect_identical(is_unit_truncated(x), truncated_units)
})

test_that("a zero-length ps_trunc slice keeps the record and answers nothing", {
  x <- trunc_record_fixture()

  proto <- expect_silent(vec_ptype(x))
  expect_length(proto, 0)
  expect_identical(ps_trunc_meta(proto), ps_trunc_meta(x))

  expect_identical(is_unit_truncated(proto), logical(0))
  expect_identical(is_unit_truncated(x[integer(0)]), logical(0))
})

test_that("a truncation record that no longer covers the scores refuses to answer", {
  x <- trunc_record_fixture()
  expect_silent({
    x[7] <- 0.5
  })

  expect_length(x, 7)
  expect_identical(ps_trunc_meta(x), ps_trunc_meta(trunc_record_fixture()))
  expect_true(is_ps_truncated(x))
  expect_error(
    is_unit_truncated(x),
    class = "propensity_missing_meta_error"
  )
})

test_that("combining ps_trunc objects drops the truncation record", {
  x <- trunc_record_fixture()

  combined <- expect_silent(c(x, x))
  expect_s3_class(combined, "ps_trunc")
  expect_length(combined, 10)
  expect_equal(as.numeric(combined), rep(as.numeric(x), 2))

  meta <- ps_trunc_meta(combined)
  expect_null(meta$truncated_idx)
  expect_null(meta$n_obs)
  expect_equal(meta$method, "ps")

  expect_true(is_ps_truncated(combined))
  expect_error(
    is_unit_truncated(combined),
    class = "propensity_missing_meta_error"
  )
})

test_that("combining ps_trunc objects does not read truncated units off the bounds", {
  # Position 4 holds a score that arrived equal to the upper bound and was left
  # alone. Reading membership back from bound equality would report it, and the
  # copy of it in the second half, as units this package winsorized.
  x <- trunc_record_fixture()
  expect_identical(
    is_unit_truncated(x),
    c(TRUE, FALSE, FALSE, FALSE, FALSE)
  )

  combined <- expect_silent(c(x, x))
  expect_null(ps_trunc_meta(combined)$truncated_idx)
  expect_error(
    is_unit_truncated(combined),
    class = "propensity_missing_meta_error"
  )
})

test_that("casting to ps_trunc records positions and a length", {
  x <- trunc_record_fixture()

  empty <- vec_cast(double(), to = x)
  expect_s3_class(empty, "ps_trunc")
  expect_true("truncated_idx" %in% names(ps_trunc_meta(empty)))
  expect_equal(ps_trunc_meta(empty)$truncated_idx, integer(0))
  expect_equal(ps_trunc_meta(empty)$n_obs, 0L)

  values <- vec_cast(c(0.3, 0.4), to = x)
  expect_s3_class(values, "ps_trunc")
  expect_equal(ps_trunc_meta(values)$truncated_idx, integer(0))
  expect_equal(ps_trunc_meta(values)$n_obs, 2L)
  expect_identical(is_unit_truncated(values), c(FALSE, FALSE))

  from_integer <- vec_cast(c(0L, 1L), to = x)
  expect_equal(ps_trunc_meta(from_integer)$truncated_idx, integer(0))
  expect_equal(ps_trunc_meta(from_integer)$n_obs, 2L)
})
