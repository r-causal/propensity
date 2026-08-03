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

  # 4) Casting double -> ps_trunc => the truncation of the target
  new_vals <- runif(3)
  out_ps_trunc <- vctrs::vec_cast(new_vals, to = obj)
  expect_s3_class(out_ps_trunc, "ps_trunc")
  meta_new <- ps_trunc_meta(out_ps_trunc)
  expect_equal(meta_new$method, "ps")
  expect_equal(meta_new$lower_bound, 0.2)
  expect_equal(meta_new$upper_bound, 0.8)
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

  # A grouped summary slices the column once per group, and each slice holds
  # scores the truncation record was not written for, so the record is dropped
  # and says so. The summary itself reads values rather than positions.
  summarized <- count_record_drops(
    tibble(x, z, ps) |>
      group_by(truncated = is_unit_truncated(ps)) |>
      summarize(mean = mean(ps), .groups = "drop")
  )
  expect_gt(summarized$drops, 0)

  out <- summarized$value
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

# A bound the method reads off the scores is a score, and carries every digit
# the score it came from was worked out to. The line a truncated vector is
# printed under reports the bound to a few significant digits instead.

test_that("vec_ptype_full() rounds a bound read off the scores", {
  set.seed(123)
  ps <- runif(20, 0.05, 0.95)

  truncated <- ps_trunc(ps, method = "pctl", lower = 0.05, upper = 0.95)
  meta <- ps_trunc_meta(truncated)

  # The quantiles the bounds came from, to the digits they were worked out to.
  expect_equal(unname(meta$lower_bound), 0.09084348598727956)
  expect_equal(unname(meta$upper_bound), 0.9091581205616238)

  expect_identical(
    vctrs::vec_ptype_full(truncated),
    "ps_trunc{[0.0908,0.909], method=pctl}"
  )
})

test_that("ps_trunc index tracking works when combining objects", {
  set.seed(456)
  ps1 <- runif(10, 0.05, 0.95)
  ps2 <- runif(10, 0.05, 0.95)

  # Create ps_trunc objects with same parameters
  ps_trunc1 <- ps_trunc(ps1, method = "ps", lower = 0.2, upper = 0.8)
  ps_trunc2 <- ps_trunc(ps2, method = "ps", lower = 0.2, upper = 0.8)

  # Combine the objects
  combined <- c(ps_trunc1, ps_trunc2)

  # Should be a ps_trunc object holding every value in order
  expect_s3_class(combined, "ps_trunc")
  expect_equal(length(combined), 20)
  expect_equal(
    vec_data(combined),
    c(vec_data(ps_trunc1), vec_data(ps_trunc2))
  )

  # Concatenation appends one set of observations to another, so the positions
  # either record names would describe units from the other input.
  combined_meta <- ps_trunc_meta(combined)
  expect_null(combined_meta$truncated_idx)
  expect_equal(combined_meta$method, "ps")
  expect_error(
    is_unit_truncated(combined),
    class = "propensity_missing_meta_error"
  )
})

test_that("ps_trunc warns when combining objects with different parameters", {
  ps1 <- runif(10, 0.05, 0.95)
  ps2 <- runif(10, 0.05, 0.95)

  # Create ps_trunc objects with different parameters
  ps_trunc1 <- ps_trunc(ps1, method = "ps", lower = 0.2, upper = 0.8)
  ps_trunc2 <- ps_trunc(ps2, method = "ps", lower = 0.3, upper = 0.7)

  # Should warn and return numeric
  combined <- expect_propensity_warning(c(ps_trunc1, ps_trunc2))

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

  # Each half kept a record written for itself, and the recombination has none:
  # the halves were subset with a subscript in hand and appended without one.
  expect_equal(
    length(ps_trunc_meta(subset1)$truncated_idx) +
      length(ps_trunc_meta(subset2)$truncated_idx),
    length(meta$truncated_idx)
  )
  expect_null(ps_trunc_meta(recombined)$truncated_idx)

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

  # Should maintain ps_trunc class and every value
  expect_s3_class(combined, "ps_trunc")
  expect_equal(length(combined), 15)
  expect_equal(
    vec_data(combined),
    c(vec_data(ps_trunc1), vec_data(ps_trunc2), vec_data(ps_trunc3))
  )

  # Folding three inputs together drops the record just as folding two does.
  expect_null(ps_trunc_meta(combined)$truncated_idx)
  expect_error(
    is_unit_truncated(combined),
    class = "propensity_missing_meta_error"
  )
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

test_that("a ps_trunc reordered through vctrs keeps the record for the old order", {
  # The documented limit of the coverage check, which counts observations and so
  # sees nothing in a reordering. No subscript reaches the restore, so the
  # record survives naming where the observations used to be.
  x <- trunc_record_fixture()

  reordered <- expect_silent(vec_slice(x, 5:1))
  expect_equal(as.numeric(reordered), c(0.6, 0.9, 0.5, 0.3, 0.1))
  expect_identical(ps_trunc_meta(reordered), ps_trunc_meta(x))

  # The winsorized unit now holds position 5, and the record still names 1, so
  # the answer is the one the record gives rather than the one the values show.
  expect_identical(
    is_unit_truncated(reordered),
    c(TRUE, FALSE, FALSE, FALSE, FALSE)
  )

  # `[` is handed the subscript and re-indexes, so the same reordering through
  # `[` reports the unit that holds that position.
  expect_identical(
    is_unit_truncated(x[5:1]),
    c(FALSE, FALSE, FALSE, FALSE, TRUE)
  )
})

# An exposure with dimensions reaches the same coding the weight functions use,
# where its cells were read in storage order rather than one value per
# observation.

test_that("ps_trunc refuses an exposure with dimensions", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  dimensioned <- matrix(c(1, 0, 1, 0), nrow = 2, ncol = 2)

  expect_error(
    ps_trunc(ps, method = "cr", .exposure = dimensioned),
    class = "propensity_binary_transform_error"
  )
})

test_that("ps_trunc refuses a focal level the exposure never takes", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c("a", "b", "a", "b")

  # A focal level nobody holds leaves every unit in the reference group, so the
  # bounds are computed over a split the caller did not ask for.
  expect_error(
    ps_trunc(
      ps,
      method = "cr",
      .exposure = exposure,
      .focal_level = "absent"
    ),
    class = "propensity_focal_level_error"
  )
})

test_that("ps_trunc names `.exposure` when the common range method requires one", {
  ps <- c(0.2, 0.4, 0.6, 0.8)

  # Without an exposure the common range method reaches the binary transform
  # with nothing to transform, and the caller is told that the exposure could
  # not be converted rather than that it was never supplied.
  cnd <- rlang::catch_cnd(ps_trunc(ps, method = "cr"), classes = "error")
  expect_s3_class(cnd, "propensity_missing_arg_error")
  expect_match(conditionMessage(cnd), "`.exposure`", fixed = TRUE)

  expect_propensity_error(ps_trunc(ps, method = "cr"))
})

# Missing values ------------------------------------------------------------

# Truncation keeps every unit, so a propensity score that arrives missing stays
# missing and is not one this package pinned to a bound. The methods that read
# their bounds off the scores read them off the scores they have.

test_that("ps_trunc() passes a score that arrived missing through unmarked", {
  # Regression guard: this is what the fixed-bound method already does, and the
  # policy the other methods are being brought to.
  truncated <- ps_trunc(
    c(0.2, 0.5, NA, 0.7),
    method = "ps",
    lower = 0.3,
    upper = 0.6
  )
  meta <- ps_trunc_meta(truncated)

  expect_equal(as.numeric(truncated), c(0.3, 0.5, NA, 0.6))
  expect_equal(meta$truncated_idx, c(1, 4))
  expect_equal(meta$n_obs, 4)
  expect_equal(is_unit_truncated(truncated), c(TRUE, FALSE, FALSE, TRUE))
})

test_that("ps_trunc() takes its percentile bounds from the complete scores", {
  ps <- c(0.05, 0.2, 0.5, NA, 0.7, 0.95)

  # `quantile()` refuses a missing value unless it is told to drop it, so the
  # percentile method is an error on any sample with one.
  with_na <- ps_trunc(ps, method = "pctl", lower = 0.2, upper = 0.8)
  without_na <- ps_trunc(ps[-4], method = "pctl", lower = 0.2, upper = 0.8)
  meta <- ps_trunc_meta(with_na)

  expect_equal(meta$lower_bound, ps_trunc_meta(without_na)$lower_bound)
  expect_equal(meta$upper_bound, ps_trunc_meta(without_na)$upper_bound)
  expect_equal(as.numeric(with_na)[-4], as.numeric(without_na))
  expect_true(is.na(as.numeric(with_na)[4]))
  expect_equal(meta$truncated_idx, c(1, 6))
  expect_false(is_unit_truncated(with_na)[4])
})

test_that("ps_trunc() takes its common range from the complete scores", {
  ps <- c(0.05, 0.2, 0.5, NA, 0.7, 0.95)
  z <- c(0, 0, 1, 1, 1, 0)

  # One missing score in the focal group makes the lower bound missing, and no
  # score compares below a missing bound, so nothing is pinned there.
  with_na <- ps_trunc(ps, method = "cr", .exposure = z, .focal_level = 1)
  without_na <- ps_trunc(
    ps[-4],
    method = "cr",
    .exposure = z[-4],
    .focal_level = 1
  )
  meta <- ps_trunc_meta(with_na)

  expect_equal(meta$lower_bound, ps_trunc_meta(without_na)$lower_bound)
  expect_equal(meta$upper_bound, ps_trunc_meta(without_na)$upper_bound)
  expect_equal(as.numeric(with_na)[-4], as.numeric(without_na))
  expect_true(is.na(as.numeric(with_na)[4]))
  expect_equal(meta$truncated_idx, c(1, 2))
  expect_false(is_unit_truncated(with_na)[4])
})

test_that("ps_trunc() refuses an exposure with missing values", {
  ps <- c(0.1, 0.2, 0.6, 0.7, 0.8, 0.9)
  z <- c(0, 0, 0, 1, 1, NA)
  z_factor <- factor(c("a", "a", "a", "b", "b", NA))

  # The common range is bounded by the extremes of each group, both of which
  # are missing once a unit belongs to neither. Nothing compares as outside a
  # missing bound, so every score is returned untouched and the object reports
  # bounds it never applied.
  cr_numeric <- rlang::catch_cnd(
    ps_trunc(ps, method = "cr", .exposure = z, .focal_level = 1),
    classes = "condition"
  )
  expect_s3_class(cr_numeric, "error")
  expect_s3_class(cr_numeric, "propensity_error")

  cr_factor <- rlang::catch_cnd(
    ps_trunc(ps, method = "cr", .exposure = z_factor, .focal_level = "b"),
    classes = "condition"
  )
  expect_s3_class(cr_factor, "error")
  expect_s3_class(cr_factor, "propensity_error")

  expect_propensity_error(
    ps_trunc(ps, method = "cr", .exposure = z, .focal_level = 1)
  )
})

# Bounds validation ---------------------------------------------------------

test_that("ps_trunc() requires lower below upper for the pctl method", {
  ps <- c(0.2, 0.4, 0.6, 0.8)

  # Bounds the wrong way around cross, and every score is then pinned to the
  # bound on the far side of the other one. Method "ps" already refuses this
  # and the percentile method owes the same refusal.
  cnd <- rlang::catch_cnd(
    ps_trunc(ps, method = "pctl", lower = 0.95, upper = 0.05),
    classes = "condition"
  )
  expect_s3_class(cnd, "error")
  expect_s3_class(cnd, "propensity_range_error")

  expect_propensity_error(
    ps_trunc(ps, method = "pctl", lower = 0.95, upper = 0.05)
  )
})

test_that("ps_trunc() refuses percentile bounds outside the unit interval", {
  ps <- c(0.2, 0.4, 0.6, 0.8)

  # For the percentile method the bounds are probabilities. `quantile()`
  # refuses one outside [0, 1] in a bare error naming `probs`, an argument
  # `ps_trunc()` does not have.
  low <- rlang::catch_cnd(
    ps_trunc(ps, method = "pctl", lower = -0.1),
    classes = "condition"
  )
  expect_s3_class(low, "error")
  expect_s3_class(low, "propensity_error")

  high <- rlang::catch_cnd(
    ps_trunc(ps, method = "pctl", upper = 1.5),
    classes = "condition"
  )
  expect_s3_class(high, "error")
  expect_s3_class(high, "propensity_error")

  expect_propensity_error(ps_trunc(ps, method = "pctl", lower = -0.1))
})

test_that("ps_trunc() refuses a bound that is missing", {
  ps <- c(0.2, 0.4, 0.6, 0.8)

  # A missing bound answers neither `TRUE` nor `FALSE` in the comparison that
  # decides which scores to pin, so the rule comes out as a bare base error
  # about a missing value rather than as an answer.
  fixed <- rlang::catch_cnd(
    ps_trunc(ps, method = "ps", lower = NA),
    classes = "condition"
  )
  expect_s3_class(fixed, "error")
  expect_s3_class(fixed, "propensity_error")

  pctl <- rlang::catch_cnd(
    ps_trunc(ps, method = "pctl", upper = NA),
    classes = "condition"
  )
  expect_s3_class(pctl, "error")
  expect_s3_class(pctl, "propensity_error")

  expect_propensity_error(ps_trunc(ps, method = "ps", lower = NA))
})

test_that("ps_trunc() refuses a common range the exposure groups do not share", {
  ps <- c(0.1, 0.2, 0.6, 0.7, 0.8, 0.9)
  z <- c(0, 0, 0, 1, 1, 1)

  # The lowest score among the exposed, 0.7, is above the highest among the
  # unexposed, 0.6, so the common range is empty. Every score is pinned to a
  # bound on the far side of the other one, which reports a common range where
  # the groups have none.
  cnd <- rlang::catch_cnd(
    ps_trunc(ps, method = "cr", .exposure = z, .focal_level = 1),
    classes = "condition"
  )
  expect_s3_class(cnd, "error")
  expect_s3_class(cnd, "propensity_error")
  expect_match(conditionMessage(cnd), "overlap")

  expect_propensity_error(
    ps_trunc(ps, method = "cr", .exposure = z, .focal_level = 1)
  )
})

test_that("ps_trunc() reports the crossed bounds at a readable precision", {
  ps <- c(0.1234567891, 0.2, 0.6, 0.7654321987, 0.8, 0.9)
  z <- c(0, 0, 0, 1, 1, 1)

  # The bounds are propensity scores read off the data, so they arrive at full
  # double precision and a message that interpolates them raw reads as noise.
  cnd <- rlang::catch_cnd(
    ps_trunc(ps, method = "cr", .exposure = z, .focal_level = 1),
    classes = "error"
  )
  msg <- conditionMessage(cnd)

  expect_match(msg, "0.765", fixed = TRUE)
  expect_no_match(msg, "0.7654321987", fixed = TRUE)
})

test_that("ps_trunc() keeps every score within the common range it records", {
  # Regression guard: where the groups do overlap, the bounds are a range and
  # every returned score lies inside it.
  truncated <- ps_trunc(
    c(0.1, 0.3, 0.6, 0.2, 0.5, 0.9),
    method = "cr",
    .exposure = c(0, 0, 0, 1, 1, 1),
    .focal_level = 1
  )
  meta <- ps_trunc_meta(truncated)

  expect_lt(meta$lower_bound, meta$upper_bound)
  expect_true(all(as.numeric(truncated) >= meta$lower_bound))
  expect_true(all(as.numeric(truncated) <= meta$upper_bound))
})

# Combining truncated propensity scores --------------------------------------

# Two `ps_trunc` objects are combined through the prototype they share. The
# prototype stands for the truncation that produced them, so it carries the
# method and the bounds that method settled on. A bound read off the scores
# cannot be worked out again from a prototype that holds none, so it is carried
# across rather than recomputed. The prototype describes no observations, so it
# names no positions and the combined result has no record.

trunc_combine_fixture <- function() {
  set.seed(4)
  n <- 40
  x <- rnorm(n, sd = 2)

  list(
    ps = plogis(x),
    exposure = rbinom(n, 1, plogis(x))
  )
}

test_that("combining ps_trunc objects keeps the bounds given to the truncation", {
  fixture <- trunc_combine_fixture()
  truncated <- ps_trunc(fixture$ps, method = "ps", lower = 0.1, upper = 0.9)

  combined <- expect_silent(c(truncated[1:20], truncated[21:40]))
  meta <- ps_trunc_meta(combined)

  expect_s3_class(combined, "ps_trunc")
  expect_length(combined, 40)
  expect_equal(as.numeric(combined), as.numeric(truncated))
  expect_equal(meta$method, "ps")
  expect_equal(meta$lower_bound, 0.1)
  expect_equal(meta$upper_bound, 0.9)
  expect_null(meta$truncated_idx)
  expect_null(meta$n_obs)
})

test_that("combining pctl ps_trunc objects keeps the bounds the truncation found", {
  fixture <- trunc_combine_fixture()
  truncated <- ps_trunc(fixture$ps, method = "pctl")
  original <- ps_trunc_meta(truncated)

  # The bounds come from the scores, and a prototype holds none, so a prototype
  # that works them out again reports missing bounds. Read back as the bounds
  # the truncation used, they no longer describe the scores being combined.
  combined <- expect_silent(c(truncated[1:20], truncated[21:40]))
  meta <- ps_trunc_meta(combined)

  expect_s3_class(combined, "ps_trunc")
  expect_length(combined, 40)
  expect_equal(as.numeric(combined), as.numeric(truncated))
  expect_equal(meta$method, "pctl")
  expect_equal(meta$lower_pctl, 0.05)
  expect_equal(meta$upper_pctl, 0.95)
  expect_false(is.na(meta$lower_bound))
  expect_false(is.na(meta$upper_bound))
  expect_equal(meta$lower_bound, original$lower_bound)
  expect_equal(meta$upper_bound, original$upper_bound)
  expect_null(meta$truncated_idx)
  expect_null(meta$n_obs)
})

test_that("combining cr ps_trunc objects keeps the common range", {
  fixture <- trunc_combine_fixture()
  truncated <- ps_trunc(
    fixture$ps,
    method = "cr",
    .exposure = fixture$exposure
  )
  original <- ps_trunc_meta(truncated)

  # The common range is defined against the exposure, which a prototype built by
  # truncating again would have to be handed and is not.
  combined <- expect_silent(c(truncated[1:20], truncated[21:40]))
  meta <- ps_trunc_meta(combined)

  expect_s3_class(combined, "ps_trunc")
  expect_length(combined, 40)
  expect_equal(as.numeric(combined), as.numeric(truncated))
  expect_equal(meta$method, "cr")
  expect_equal(meta$lower_bound, original$lower_bound)
  expect_equal(meta$upper_bound, original$upper_bound)
  expect_null(meta$truncated_idx)
  expect_null(meta$n_obs)
})

test_that("casting a double to a ps_trunc keeps the truncation of the target", {
  to <- ps_trunc(
    c(0.05, 0.3, 0.5, 0.95),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )

  # A cast returns the values it was given in the type it was given, and the
  # bounds are part of that type.
  out <- vec_cast(c(0.3, 0.4), to = to)
  meta <- ps_trunc_meta(out)

  expect_s3_class(out, "ps_trunc")
  expect_equal(as.numeric(out), c(0.3, 0.4))
  expect_equal(meta$method, "ps")
  expect_equal(meta$lower_bound, 0.1)
  expect_equal(meta$upper_bound, 0.9)
  expect_equal(meta$truncated_idx, integer(0))
  expect_equal(meta$n_obs, 2L)
})

test_that("combining a ps_trunc with an integer keeps the propensity scores", {
  x <- ps_trunc(c(0.2, 0.5, 0.85), method = "ps", lower = 0.1, upper = 0.9)

  # Propensity scores lie strictly between 0 and 1, so a combination that meets
  # an integer in the integers is every score rounded away.
  combined <- expect_propensity_warning(vec_c(x, 1L))

  expect_type(combined, "double")
  expect_equal(combined, c(0.2, 0.5, 0.85, 1))
})

test_that("casting a ps_trunc to integer refuses rather than rounds", {
  x <- ps_trunc(c(0.2, 0.5, 0.85), method = "ps", lower = 0.1, upper = 0.9)

  expect_error(
    vec_cast(x, integer()),
    class = "vctrs_error_cast_lossy"
  )
})

# Choosing a column from a data frame of propensity scores -------------------

# Predictions from a binary model arrive as a column per level, and truncation
# works on one column. The convention is the second column of a pair, which is
# the probability of the second level in the layout these predictions come in,
# and the only column otherwise. The caller did not make that choice, so it is
# announced rather than left to be inferred from the result.

trunc_frame_fixture <- function() {
  set.seed(2)
  n <- 40
  x <- rnorm(n)
  p <- plogis(x)

  list(
    ps = data.frame(.pred_0 = 1 - p, .pred_1 = p),
    exposure = rbinom(n, 1, p)
  )
}

test_that("ps_trunc() names the column it took from a data frame of two", {
  withr::local_options(propensity.quiet = FALSE)
  fixture <- trunc_frame_fixture()

  messages <- testthat::capture_messages(
    ps_trunc(
      fixture$ps,
      .exposure = fixture$exposure,
      method = "ps",
      lower = 0.3,
      upper = 0.7
    )
  )

  naming <- messages[grepl(".pred_1", messages, fixed = TRUE)]
  expect_length(naming, 1)
  expect_false(any(grepl(".pred_0", messages, fixed = TRUE)))
})

test_that("ps_trunc() names the only column of a one column data frame", {
  withr::local_options(propensity.quiet = FALSE)
  fixture <- trunc_frame_fixture()
  one_column <- fixture$ps[, ".pred_1", drop = FALSE]

  messages <- testthat::capture_messages(
    ps_trunc(one_column, method = "ps", lower = 0.3, upper = 0.7)
  )

  expect_length(messages, 1)
  expect_true(grepl(".pred_1", messages, fixed = TRUE))
})

test_that("ps_trunc() announces no column when the messages are quieted", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- trunc_frame_fixture()

  expect_silent(
    ps_trunc(
      fixture$ps,
      .exposure = fixture$exposure,
      method = "ps",
      lower = 0.3,
      upper = 0.7
    )
  )
})

test_that("the column ps_trunc() announces is named in a full sentence", {
  withr::local_options(propensity.quiet = FALSE)
  fixture <- trunc_frame_fixture()

  expect_snapshot(
    truncated <- ps_trunc(fixture$ps, method = "ps", lower = 0.3, upper = 0.7)
  )
})

test_that("ps_trunc() takes the second column of a data frame of two", {
  fixture <- trunc_frame_fixture()

  from_frame <- ps_trunc(fixture$ps, method = "ps", lower = 0.3, upper = 0.7)
  from_second <- ps_trunc(
    fixture$ps[[2]],
    method = "ps",
    lower = 0.3,
    upper = 0.7
  )
  from_first <- ps_trunc(
    fixture$ps[[1]],
    method = "ps",
    lower = 0.3,
    upper = 0.7
  )

  expect_gt(length(ps_trunc_meta(from_second)$truncated_idx), 0)
  expect_equal(as.numeric(from_frame), as.numeric(from_second))
  expect_equal(
    ps_trunc_meta(from_frame)$truncated_idx,
    ps_trunc_meta(from_second)$truncated_idx
  )
  expect_false(identical(as.numeric(from_frame), as.numeric(from_first)))
})
