test_that("ps_trim() - Basic structure and return types", {
  set.seed(42)

  n <- 50
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))

  # Fit a logistic regression
  fit <- glm(z ~ x, family = binomial)
  ps_vec <- predict(fit, type = "response")

  # 1) A default call with method="ps"
  out <- ps_trim(ps_vec, method = "ps")

  # Basic checks
  # Now 'out' is a ps_trim object
  expect_s3_class(out, "ps_trim")
  # The underlying data: same length
  expect_equal(length(out), length(ps_vec))

  # Inspect the internal meta
  meta <- ps_trim_meta(out)
  expect_true(is.list(meta))
  # e.g. method, lower, upper, keep_idx, trimmed_idx
  expect_true(all(
    c("method", "lower", "upper", "keep_idx", "trimmed_idx") %in% names(meta)
  ))

  # Check that the kept indices and trimmed indices are disjoint
  expect_length(intersect(meta$keep_idx, meta$trimmed_idx), 0)

  # Confirm that out-of-range entries are NA
  # By default, [0.1, 0.9]
  below_min <- ps_vec < 0.1
  above_max <- ps_vec > 0.9
  expect_true(all(is.na(out[below_min])))
  expect_true(all(is.na(out[above_max])))
  # The rest remain the same
})

test_that("ps method: default and custom cutoffs", {
  set.seed(1)

  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 1) Default cutoffs (0.1, 0.9)
  out1 <- ps_trim(ps, method = "ps")
  meta1 <- ps_trim_meta(out1)
  expect_equal(meta1$lower, 0.1)
  expect_equal(meta1$upper, 0.9)

  # 2) Custom cutoffs
  out2 <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  meta2 <- ps_trim_meta(out2)
  expect_equal(meta2$lower, 0.2)
  expect_equal(meta2$upper, 0.8)

  # Check that out-of-range entries are NA
  # i.e. everything <0.2 or >0.8 is NA
  out2_data <- as.numeric(out2)
  expect_true(all(out2_data[!is.na(out2_data)] >= 0.2))
  expect_true(all(out2_data[!is.na(out2_data)] <= 0.8))
})

test_that("adaptive method: ignores lower/upper, warns appropriately", {
  set.seed(2)

  n <- 80
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(-0.5 * x))
  ps <- predict(glm(z ~ x, family = binomial), type = "response")

  # 1) No user cutoffs
  out_adapt <- ps_trim(ps, method = "adaptive")
  meta_adapt <- ps_trim_meta(out_adapt)
  # The meta should have a 'cutoff' field
  expect_true("cutoff" %in% names(meta_adapt))

  # 2) If user sets lower/upper, we expect a warning
  expect_propensity_warning(
    ps_trim(
      ps,
      method = "adaptive",
      lower = 0.2,
      upper = 0.8
    )
  )
})

test_that("pctl method: percentile-based trimming", {
  set.seed(3)

  n <- 100
  x <- rnorm(n)
  ps <- plogis(0.8 * x)

  # 1) Default [0.05, 0.95]
  out1 <- ps_trim(ps, method = "pctl")
  meta1 <- ps_trim_meta(out1)
  expect_equal(meta1$lower, 0.05)
  expect_equal(meta1$upper, 0.95)

  q_l <- quantile(ps, probs = 0.05)
  q_u <- quantile(ps, probs = 0.95)
  out1_data <- as.numeric(out1)

  # Everything below q_l is NA, above q_u is NA
  expect_true(all(is.na(out1_data[ps < q_l])))
  expect_true(all(is.na(out1_data[ps > q_u])))

  # 2) Custom [0.2, 0.8]
  out2 <- ps_trim(ps, method = "pctl", lower = 0.2, upper = 0.8)
  meta2 <- ps_trim_meta(out2)
  q_l2 <- quantile(ps, probs = 0.2)
  q_u2 <- quantile(ps, probs = 0.8)
  out2_data <- as.numeric(out2)

  expect_true(all(is.na(out2_data[ps < q_l2])))
  expect_true(all(is.na(out2_data[ps > q_u2])))
})

test_that("pref method: requires exposure, fails with all 0 or all 1", {
  set.seed(4)

  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 1) If exposure = NULL, should fail
  expect_propensity_error(
    ps_trim(ps, method = "pref")
  )

  # 2) If exposure is all 0 or all 1 => fail
  expect_propensity_error(
    ps_trim(ps, .exposure = rep(0, n), method = "pref")
  )

  expect_propensity_error(
    ps_trim(ps, .exposure = rep(1, n), method = "pref")
  )

  # 3) Valid usage
  out_pref <- ps_trim(ps, .exposure = z, method = "pref", .focal_level = 1)
  meta_pref <- ps_trim_meta(out_pref)
  expect_equal(meta_pref$lower, 0.3)
  expect_equal(meta_pref$upper, 0.7)

  # Check final
  out_pref_data <- as.numeric(out_pref)
  # We know that we just set NA outside [0.3, 0.7] in preference-score space,
  # but the underlying values remain in [0,1].
  # So let's just confirm it is indeed a ps_trim object
  expect_s3_class(out_pref, "ps_trim")
})

test_that("cr method: uses min(ps_treat) / max(ps_untrt), warns if cutoffs given", {
  set.seed(5)

  n <- 50
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # Must have exposure
  expect_propensity_error(
    ps_trim(ps, method = "cr")
  )

  # If all 0 or all 1 => fail
  expect_propensity_error(
    ps_trim(ps, .exposure = rep(1, n), method = "cr")
  )

  # Valid usage
  out_cr <- ps_trim(ps, .exposure = z, method = "cr", .focal_level = 1)
  meta_cr <- ps_trim_meta(out_cr)
  ps_treat <- ps[z == 1]
  ps_untrt <- ps[z == 0]
  cr_l_exp <- min(ps_treat)
  cr_u_exp <- max(ps_untrt)
  expect_equal(meta_cr$cr_lower, cr_l_exp)
  expect_equal(meta_cr$cr_upper, cr_u_exp)

  # Check that user-specified lower/upper => warning
  expect_propensity_warning(
    ps_trim(
      ps,
      .exposure = z,
      method = "cr",
      lower = 0.2,
      upper = 0.8,
      .focal_level = 1
    ),
    print = TRUE
  )
})

test_that("Edge cases: ps near 0 or 1, empty trimming result", {
  ps_edge <- c(0.0001, 0.01, 0.5, 0.99, 0.9999)

  # If we do normal cut: e.g. [0.01, 0.99]
  out <- ps_trim(ps_edge, method = "ps", lower = 0.01, upper = 0.99)
  # out is same length as ps_edge, but the extremely small/large => NA
  out_data <- as.numeric(out)
  expect_true(is.na(out_data[1])) # 0.0001 <0.01 => NA
  expect_true(is.na(out_data[5])) # 0.9999>0.99 => NA
  expect_equal(sum(!is.na(out_data)), 3)

  # If we force a scenario with [1.1,1.2] => everything is out => all NA
  out_empty <- ps_trim(ps_edge, method = "ps", lower = 1.1, upper = 1.2)
  out_empty_data <- as.numeric(out_empty)
  expect_true(all(is.na(out_empty_data)))
  meta_e <- ps_trim_meta(out_empty)
  expect_length(meta_e$keep_idx, 0)
  expect_length(meta_e$trimmed_idx, length(ps_edge))
})

test_that("ps_refit() refits on keep_idx, warns if everything trimmed, etc.", {
  set.seed(123)
  n <- 20
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # Trim to [0.2, 0.8], then refit
  out <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  # Suppose we do a normal refit
  refit_out <- ps_refit(out, model = fit)
  expect_s3_class(refit_out, "ps_trim")
  meta_r <- ps_trim_meta(refit_out)
  expect_true(isTRUE(meta_r$refit))

  expect_propensity_error(
    ps_refit(out, model = fit, .data = data.frame(z, x)[1:10, ])
  )

  # If everything is trimmed => error
  ps_edge <- c(0.01, 0.01, 0.99, 0.99)
  z_edge <- c(0, 1, 1, 0)
  out_empty <- ps_trim(ps_edge, method = "ps", lower = 1.1, upper = 2)
  expect_propensity_error(
    ps_refit(out_empty, model = fit)
  )

  ps_trim(ps_edge, method = "ps", lower = 1.1, upper = 2)
})

test_that("Full workflow: trim -> refit -> weighting yields refit, trimmed psw", {
  set.seed(999)
  n <- 12
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))

  # 1) Fit initial logistic model, get ps
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 2) Trim the PS (e.g. method="ps" with [0.2, 0.8])
  trimmed_ps <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  expect_false(is_refit(trimmed_ps)) # not refit yet

  # 3) Refit on the subset
  trimmed_refit <- ps_refit(trimmed_ps, model = fit)
  expect_true(is_refit(trimmed_refit)) # now refit=TRUE in ps_trim_meta

  # 4) Create ATE weights with the refitted ps_trim object
  w_ate <- wt_ate(
    trimmed_refit,
    .exposure = z,
    exposure_type = "binary",
    .focal_level = 1
  )

  # 5) Check final psw object
  expect_s3_class(w_ate, "psw")

  # Should be trimmed, per the weighting method's logic
  expect_true(is_ps_trimmed(w_ate))
  # Should NOT be truncated or stabilized
  expect_false(is_ps_truncated(w_ate))
  expect_false(is_stabilized(w_ate))

  # Should preserve the refit info if you attach ps_trim_meta
  expect_true(is_refit(w_ate)) # e.g. if is_refit.psw() checks ps_trim_meta

  # The estimand should include "; trimmed"
  expect_match(estimand(w_ate), "; trimmed$")
})

test_that("adaptive method: triggers uniroot path (k < 0) coverage", {
  # We'll craft a scenario where k = 2*mean(sum_wt) - max(sum_wt) < 0
  # so that the else-branch is executed (lines 83-89).

  set.seed(1234)
  n <- 30
  # We create some extreme p near 0 or 1 so sum_wt = 1/(p*(1-p)) varies greatly
  # e.g. half near 0.01, half near 0.99
  p_vec <- c(
    runif(n / 2, min = 0.001, max = 0.01),
    runif(n / 2, min = 0.99, max = 0.999)
  )
  # We'll treat them as if we have a binary z, not relevant for "adaptive"
  z <- c(rep(0, n / 2), rep(1, n / 2))

  # Now call ps_trim with method="adaptive"
  # This should produce k < 0 => code path with uniroot
  out_adapt <- ps_trim(
    p_vec,
    .exposure = z,
    method = "adaptive",
    .focal_level = 1
  )

  # Check that the 'cutoff' field in the meta is present
  meta <- ps_trim_meta(out_adapt)
  expect_true("cutoff" %in% names(meta))

  # Also confirm the result is a ps_trim object
  expect_s3_class(out_adapt, "ps_trim")

  # Because of the extremes, we likely see a fairly small cutoff
  # Just check it's numeric and within (0, 0.5)
  expect_true(is.numeric(meta$cutoff))
  expect_gt(meta$cutoff, 0)
  expect_lt(meta$cutoff, 0.5)

  # Since many ps are out-of-range, we expect many NAs
  out_data <- as.numeric(out_adapt)
  # Just confirm there's at least some NA for the extreme values
  # e.g. near 0.001 or 0.999
  expect_true(anyNA(out_data))
})

test_that("Check defaults for helper functions", {
  # 1) Any random object that is not ps_trim => default method => FALSE
  not_trim_obj <- c(0.2, 0.4, 0.6)
  expect_false(is_ps_trimmed(not_trim_obj)) # triggers is_ps_trimmed.default()
  expect_false(is_refit(not_trim_obj))

  # 2) A mock ps_trim object => method => TRUE
  # For a real test, you'd create it via ps_trim(...). Here we simulate:
  fake_ps_trim <- structure(
    c(0.5, NA, 0.7),
    class = "ps_trim"
  )
  expect_true(is_ps_trimmed(fake_ps_trim)) # triggers is_ps_trimmed.ps_trim()
})

# tests/test-ps_trim-vctrs.R

library(testthat)
library(vctrs)

test_that("vec_ptype_abbr.ps_trim() and vec_ptype_full.ps_trim() coverage", {
  # Create a minimal ps_trim object
  # for demonstration, or use ps_trim() function if you like
  ps_obj <- new_trimmed_ps(
    c(0.1, NA, 0.7),
    ps_trim_meta = list(
      method = "ps",
      keep_idx = c(1, 3),
      trimmed_idx = 2
    )
  )

  # 1) Abbreviation
  abbr <- vec_ptype_abbr(ps_obj)
  expect_identical(abbr, "ps_trim")

  # 2) Full
  full <- vec_ptype_full(ps_obj)
  # E.g. "ps_trim; trimmed 1 of "
  # Just check it's a character containing "ps_trim"
  expect_true(is.character(full))
  expect_true(grepl("ps_trim;", full))
})

test_that("Arithmetic with ps_trim returns numeric", {
  # Create two ps_trim objects or combine with numeric
  x <- new_trimmed_ps(c(0.2, 0.3), ps_trim_meta = list())
  y <- new_trimmed_ps(c(0.4, 0.9), ps_trim_meta = list())

  # Arithmetic operations should return numeric
  expect_type(x + 1, "double")
  expect_type(x + 1L, "double")
  expect_type(1 + x, "double")
  expect_type(1L + x, "double")
  expect_type(x + y, "double")

  # Verify values are correct
  expect_equal(x + 1, c(1.2, 1.3))
  expect_equal(1 / x, c(5, 10 / 3))
  expect_equal(x * 2, c(0.4, 0.6))

  # List operations still fail
  expect_error(x + list(1))
})

test_that("Combining two ps_trim with different parameters triggers warning", {
  x <- ps_trim(
    c(0.2, 0.4, 0.8),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )
  y <- ps_trim(
    c(0.3, 0.5, 0.7),
    method = "ps",
    lower = 0.2, # Different lower bound
    upper = 0.8 # Different upper bound
  )

  # Attempt to combine with different parameters
  # This will warn about different trimming parameters and return numeric
  result <- expect_propensity_warning(vec_c(x, y))
  expect_type(result, "double")
})

test_that("Combining ps_trim with double => double", {
  x <- new_trimmed_ps(c(0.2, 0.5), ps_trim_meta = list())

  # vctrs logic => ptype2 => double
  combined <- expect_propensity_warning(vec_c(x, 0.7))
  expect_type(combined, "double")
  expect_false(inherits(combined, "ps_trim"))
})

test_that("Casting ps_trim -> double => underlying numeric data", {
  x <- new_trimmed_ps(
    c(0.2, NA, 0.9),
    ps_trim_meta = list(method = "ps", keep_idx = c(1, 3), trimmed_idx = 2)
  )
  casted <- vec_cast(x, to = double())
  expect_type(casted, "double")
  # Should match the underlying data
  expect_equal(casted, c(0.2, NA, 0.9))
})

test_that("Casting double -> ps_trim keeps the trimming of the target", {
  base_vec <- c(0.1, 0.7, NA, 0.4)
  to <- ps_trim(c(0.2, 0.5, 0.85), method = "ps", lower = 0.1, upper = 0.9)

  # A cast returns the values it was given in the type it was given, and the
  # trimming is part of that type.
  ps_t <- vec_cast(base_vec, to = to)
  expect_s3_class(ps_t, "ps_trim")
  meta <- attr(ps_t, "ps_trim_meta")
  expect_equal(meta$method, "ps")
  expect_equal(meta$lower, 0.1)
  expect_equal(meta$upper, 0.9)
  expect_equal(meta$keep_idx, seq_along(base_vec))
  expect_length(meta$trimmed_idx, 0)
  expect_equal(meta$n_obs, length(base_vec))
})

test_that("Casting integer->ps_trim likewise uses new_trimmed_ps", {
  base_int <- c(0L, 1L, 999L)
  ps_t <- vec_cast(base_int, to = structure(double(), class = "ps_trim"))
  expect_s3_class(ps_t, "ps_trim")
  # check the data is double
  expect_equal(as.numeric(ps_t), c(0, 1, 999))
})

test_that("is_unit_trimmed.ps_trim returns expected row-level booleans", {
  set.seed(100)
  ps_vec <- c(0.1, 0.2, 0.5, 0.85, 0.95)

  # Trim outside [0.2, 0.8]
  trimmed_obj <- ps_trim(
    ps_vec,
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  expect_s3_class(trimmed_obj, "ps_trim")

  row_trim <- is_unit_trimmed(trimmed_obj)
  expect_type(row_trim, "logical")
  expect_length(row_trim, length(ps_vec))

  expect_equal(which(row_trim), c(1, 4, 5))
})


test_that("ps_trim objects can convert to character", {
  ps <- c(0.01, 0.1, 0.3, 0.8, 0.95)
  out <- as.character(ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8))
  expect_type(out, "character")
})

test_that("ps_trim works with summarize(mean = mean(ps))", {
  skip_if_not_installed("dplyr")
  library(dplyr, warn.conflicts = FALSE)

  set.seed(200)
  n <- 600
  x <- rnorm(n)
  z <- rbinom(n, size = 1, prob = plogis(x + rnorm(n)))
  fit <- glm(z ~ x, family = binomial)

  ps <- predict(fit, type = "response") |>
    ps_trim(method = "ps", lower = 0.3, upper = 0.7) |>
    ps_refit(fit)

  # A grouped summary slices the column once per group, and each slice holds
  # scores the trimming record was not written for, so the record is dropped and
  # says so. The summary itself reads values rather than positions.
  summarized <- count_record_drops(
    tibble(x, z, ps) |>
      group_by(trimmed = is_unit_trimmed(ps)) |>
      summarize(mean = mean(ps), .groups = "drop")
  )
  expect_gt(summarized$drops, 0)

  out <- summarized$value
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("trimmed", "mean"))
  expect_type(out$mean, "double")
})

test_that("ps_trim errors when exposure is missing for methods that require it", {
  ps <- runif(20, 0.1, 0.9)

  # Test pref method without exposure
  expect_propensity_error(
    ps_trim(ps, method = "pref")
  )

  # Test cr method without exposure
  expect_propensity_error(
    ps_trim(ps, method = "cr")
  )

  # Should work fine with ps method (no exposure needed)
  expect_no_error(ps_trim(ps, method = "ps"))
})

test_that("ps_trim index tracking works when combining objects", {
  set.seed(456)
  ps1 <- runif(10, 0.05, 0.95)
  ps2 <- runif(10, 0.05, 0.95)

  # Create ps_trim objects with same parameters
  ps_trim1 <- ps_trim(ps1, method = "ps", lower = 0.2, upper = 0.8)
  ps_trim2 <- ps_trim(ps2, method = "ps", lower = 0.2, upper = 0.8)

  # Combine the objects
  combined <- c(ps_trim1, ps_trim2)

  # Should be a ps_trim object holding every value in order
  expect_s3_class(combined, "ps_trim")
  expect_equal(length(combined), 20)
  expect_equal(
    vec_data(combined),
    c(vec_data(ps_trim1), vec_data(ps_trim2))
  )

  # Concatenation appends one set of observations to another, so the positions
  # either record names would describe units from the other input.
  combined_meta <- ps_trim_meta(combined)
  expect_null(combined_meta$trimmed_idx)
  expect_null(combined_meta$keep_idx)
  expect_equal(combined_meta$method, "ps")
  expect_error(
    is_unit_trimmed(combined),
    class = "propensity_missing_meta_error"
  )
})

test_that("ps_trim warns when combining objects with different parameters", {
  ps1 <- runif(10, 0.05, 0.95)
  ps2 <- runif(10, 0.05, 0.95)

  # Create ps_trim objects with different parameters
  ps_trim1 <- ps_trim(ps1, method = "ps", lower = 0.2, upper = 0.8)
  ps_trim2 <- ps_trim(ps2, method = "ps", lower = 0.3, upper = 0.7)

  # Should warn and return numeric
  combined <- expect_propensity_warning(c(ps_trim1, ps_trim2))

  expect_type(combined, "double")
  expect_false(inherits(combined, "ps_trim"))
})

test_that("ps_trim index tracking works with subsetting and combining", {
  set.seed(789)
  ps <- runif(20, 0.05, 0.95)

  # Create ps_trim object
  ps_trim_obj <- ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7)
  meta <- ps_trim_meta(ps_trim_obj)

  # Subset the object
  subset1 <- ps_trim_obj[1:10]
  subset2 <- ps_trim_obj[11:20]

  # Recombine
  recombined <- c(subset1, subset2)

  # Should maintain ps_trim class
  expect_s3_class(recombined, "ps_trim")

  # Each half kept a record written for itself, and the recombination has none:
  # the halves were subset with a subscript in hand and appended without one.
  expect_equal(length(ps_trim_meta(subset1)$trimmed_idx), sum(is.na(subset1)))
  expect_null(ps_trim_meta(recombined)$trimmed_idx)

  # Check that NA values are preserved at correct positions
  recombined_data <- vec_data(recombined)
  original_data <- vec_data(ps_trim_obj)
  expect_equal(which(is.na(recombined_data)), which(is.na(original_data)))
})

test_that("ps_trim handles multiple combines correctly", {
  set.seed(321)

  # Create three ps_trim objects
  ps1 <- runif(5, 0.05, 0.95)
  ps2 <- runif(5, 0.05, 0.95)
  ps3 <- runif(5, 0.05, 0.95)

  ps_trim1 <- ps_trim(ps1, method = "ps", lower = 0.25, upper = 0.75)
  ps_trim2 <- ps_trim(ps2, method = "ps", lower = 0.25, upper = 0.75)
  ps_trim3 <- ps_trim(ps3, method = "ps", lower = 0.25, upper = 0.75)

  # Combine all three
  combined <- c(ps_trim1, ps_trim2, ps_trim3)

  # Should maintain ps_trim class and every value
  expect_s3_class(combined, "ps_trim")
  expect_equal(length(combined), 15)
  expect_equal(
    vec_data(combined),
    c(vec_data(ps_trim1), vec_data(ps_trim2), vec_data(ps_trim3))
  )

  # Folding three inputs together drops the record just as folding two does.
  expect_null(ps_trim_meta(combined)$trimmed_idx)
  expect_error(
    is_unit_trimmed(combined),
    class = "propensity_missing_meta_error"
  )
})

# Trimming record integrity ------------------------------------------------

# The record a `ps_trim` carries names positions among the observations it was
# written for. An operation that changes how many observations there are is not
# handed the subscript, so it cannot re-index the record onto the result and
# drops it rather than leave positions describing units they were never about.

trim_record_fixture <- function() {
  ps_trim(
    c(0.05, 0.3, 0.5, 0.95, 0.6),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )
}

test_that("ps_trim() records how many observations its positions describe", {
  meta <- ps_trim_meta(trim_record_fixture())

  expect_equal(meta$trimmed_idx, c(1L, 4L))
  expect_equal(meta$keep_idx, c(2L, 3L, 5L))
  expect_equal(meta$n_obs, 5L)
})

test_that("slicing a ps_trim shorter drops the trimming record with a warning", {
  x <- trim_record_fixture()

  cnd <- expect_warning(
    sliced <- vec_slice(x, 2:3),
    class = "propensity_trim_record_warning"
  )
  expect_s3_class(cnd, "propensity_warning")

  expect_s3_class(sliced, "ps_trim")
  expect_equal(as.numeric(sliced), c(0.3, 0.5))

  meta <- ps_trim_meta(sliced)
  expect_null(meta$trimmed_idx)
  expect_null(meta$keep_idx)
  expect_null(meta$n_obs)

  # Nothing that is not indexed by observation is touched by the drop.
  expect_equal(meta$method, "ps")
  expect_equal(meta$lower, 0.1)
  expect_equal(meta$upper, 0.9)
  expect_true(is_ps_trimmed(sliced))

  expect_error(
    is_unit_trimmed(sliced),
    class = "propensity_missing_meta_error"
  )
})

test_that("filtering a ps_trim column drops the trimming record with a warning", {
  skip_if_not_installed("dplyr")

  df <- data.frame(id = 1:5)
  df$ps <- trim_record_fixture()

  expect_warning(
    filtered <- dplyr::filter(df, id %in% 2:3),
    class = "propensity_trim_record_warning"
  )

  expect_s3_class(filtered$ps, "ps_trim")
  expect_equal(as.numeric(filtered$ps), c(0.3, 0.5))
  expect_null(ps_trim_meta(filtered$ps)$trimmed_idx)
  expect_error(
    is_unit_trimmed(filtered$ps),
    class = "propensity_missing_meta_error"
  )
})

test_that("a length-preserving ps_trim restore keeps the trimming record", {
  x <- trim_record_fixture()
  meta <- ps_trim_meta(x)
  trimmed_units <- c(TRUE, FALSE, FALSE, TRUE, FALSE)

  whole <- expect_silent(vec_slice(x, seq_along(x)))
  expect_identical(ps_trim_meta(whole), meta)
  expect_identical(is_unit_trimmed(whole), trimmed_units)

  empty_subscript <- expect_silent(x[])
  expect_identical(ps_trim_meta(empty_subscript), meta)

  # `[<-` casts the replacement and then leaves base R to assign into the
  # target, so the record never leaves the observations it was written for.
  expect_silent({
    x[2] <- 0.4
  })
  expect_identical(ps_trim_meta(x), meta)
  expect_identical(is_unit_trimmed(x), trimmed_units)
})

test_that("a zero-length ps_trim slice keeps the record and answers nothing", {
  x <- trim_record_fixture()
  meta <- ps_trim_meta(x)

  proto <- expect_silent(vec_ptype(x))
  expect_length(proto, 0)
  expect_identical(ps_trim_meta(proto), meta)

  # No observations, no answers. Indexing an empty logical by the positions the
  # record names would grow one as long as the original scores.
  expect_identical(is_unit_trimmed(proto), logical(0))
  expect_identical(is_unit_trimmed(x[integer(0)]), logical(0))
})

test_that("a trimming record that no longer covers the scores refuses to answer", {
  # Growing by subassignment crosses a length change no restore ever sees: the
  # record describes five observations and the scores now hold seven.
  x <- trim_record_fixture()
  expect_silent({
    x[7] <- 0.5
  })

  expect_length(x, 7)
  expect_identical(ps_trim_meta(x), ps_trim_meta(trim_record_fixture()))
  expect_true(is_ps_trimmed(x))
  expect_error(is_unit_trimmed(x), class = "propensity_missing_meta_error")

  # A refit is one fact about the propensity model rather than about any
  # position, so a record that no longer covers these scores still answers it.
  expect_false(expect_silent(is_refit(x)))
})

test_that("a record re-attached to shorter scores refuses rather than answers", {
  x <- trim_record_fixture()
  short <- new_trimmed_ps(c(0.3, 0.5), ps_trim_meta = ps_trim_meta(x))

  expect_length(short, 2)
  expect_error(is_unit_trimmed(short), class = "propensity_missing_meta_error")
})

test_that("ps_refit() refuses a ps_trim whose record was dropped", {
  set.seed(11)
  n <- 40
  covariate <- rnorm(n)
  exposure <- rbinom(n, 1, plogis(0.5 * covariate))
  model <- glm(exposure ~ covariate, family = binomial)
  trimmed <- ps_trim(
    predict(model, type = "response"),
    method = "ps",
    lower = 0.3,
    upper = 0.7
  )

  expect_warning(
    sliced <- vec_slice(trimmed, 1:20),
    class = "propensity_trim_record_warning"
  )

  expect_error(
    ps_refit(sliced, model),
    class = "propensity_missing_meta_error"
  )
})

test_that("combining ps_trim objects drops the trimming record", {
  x <- trim_record_fixture()

  combined <- expect_silent(c(x, x))
  expect_s3_class(combined, "ps_trim")
  expect_length(combined, 10)
  expect_equal(as.numeric(combined), rep(as.numeric(x), 2))

  # Concatenation appends one set of observations to another, so the positions
  # either record names would describe units from the other input.
  meta <- ps_trim_meta(combined)
  expect_null(meta$trimmed_idx)
  expect_null(meta$keep_idx)
  expect_null(meta$n_obs)
  expect_equal(meta$method, "ps")

  expect_true(is_ps_trimmed(combined))
  expect_error(
    is_unit_trimmed(combined),
    class = "propensity_missing_meta_error"
  )
})

test_that("combining ps_trim objects does not read trimmed units off the NAs", {
  # These scores record that nothing was trimmed, and one of them arrived
  # missing. Reading membership back from the `NA` pattern would report a score
  # that was never there as one this package removed.
  untrimmed <- vec_cast(
    c(0.3, NA, 0.5),
    to = structure(double(), class = "ps_trim")
  )
  expect_identical(is_unit_trimmed(untrimmed), c(FALSE, FALSE, FALSE))

  combined <- expect_silent(c(untrimmed, untrimmed))
  expect_length(combined, 6)
  expect_null(ps_trim_meta(combined)$trimmed_idx)
  expect_error(
    is_unit_trimmed(combined),
    class = "propensity_missing_meta_error"
  )
})

test_that("a ps_trim that lost its record says so instead of reporting none", {
  x <- trim_record_fixture()
  expect_warning(
    sliced <- vec_slice(x, 2:3),
    class = "propensity_trim_record_warning"
  )

  expect_match(vec_ptype_full(x), "trimmed 2 of", fixed = TRUE)
  expect_match(vec_ptype_full(sliced), "record dropped", fixed = TRUE)
})

test_that("a ps_trim reordered through vctrs keeps the record for the old order", {
  # The documented limit of the coverage check, which counts observations and so
  # sees nothing in a reordering. No subscript reaches the restore, so the
  # record survives naming where the observations used to be.
  x <- trim_record_fixture()

  reordered <- expect_silent(vec_slice(x, 5:1))
  expect_equal(as.numeric(reordered), c(0.6, NA, 0.5, 0.3, NA))
  expect_identical(ps_trim_meta(reordered), ps_trim_meta(x))

  # The trimmed units now hold positions 2 and 5, and the record still names 1
  # and 4, so the answer is the one the record gives rather than the one the
  # values show.
  expect_identical(
    is_unit_trimmed(reordered),
    c(TRUE, FALSE, FALSE, TRUE, FALSE)
  )

  # `[` is handed the subscript and re-indexes, so the same reordering through
  # `[` reports the units that hold those positions.
  expect_identical(is_unit_trimmed(x[5:1]), c(FALSE, TRUE, FALSE, FALSE, TRUE))
})

# An exposure with dimensions reaches the same coding the weight functions use,
# where its cells were read in storage order rather than one value per
# observation.

test_that("ps_trim refuses an exposure with dimensions", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  dimensioned <- matrix(c(1, 0, 1, 0), nrow = 2, ncol = 2)

  expect_error(
    ps_trim(ps, method = "pref", .exposure = dimensioned),
    class = "propensity_binary_transform_error"
  )
  expect_error(
    ps_trim(ps, method = "cr", .exposure = dimensioned),
    class = "propensity_binary_transform_error"
  )
})

test_that("ps_trim refuses a focal level the exposure never takes", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c("a", "b", "a", "b")

  # A focal level nobody holds leaves every unit in the reference group, so the
  # bounds are computed over a split the caller did not ask for.
  expect_error(
    ps_trim(ps, method = "pref", .exposure = exposure, .focal_level = "absent"),
    class = "propensity_focal_level_error"
  )
  expect_error(
    ps_trim(ps, method = "cr", .exposure = exposure, .focal_level = "absent"),
    class = "propensity_focal_level_error"
  )
})

test_that("ps_trim rejects the categorical-only optimal method on a vector", {
  set.seed(11)
  n <- 40
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  ps <- predict(glm(z ~ x, family = binomial), type = "response")

  # Optimal trimming is defined over the rows of a propensity score matrix. On
  # a vector the method falls through to common-range trimming and records
  # itself as "optimal", so the object misreports what was done to it.
  cnd <- rlang::catch_cnd(
    ps_trim(ps, method = "optimal", .exposure = z),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_wt_not_supported_error")

  # The rejection follows from the method and the type of `ps`, so it must not
  # depend on whether an exposure was supplied.
  cnd_no_exposure <- rlang::catch_cnd(
    ps_trim(ps, method = "optimal"),
    classes = "error"
  )
  expect_s3_class(cnd_no_exposure, "propensity_wt_not_supported_error")

  # The message has to point at the input the method is defined for.
  expect_propensity_error(ps_trim(ps, method = "optimal", .exposure = z))
})

test_that("ps_trim names `.exposure` when the method requires one", {
  ps <- c(0.2, 0.4, 0.6, 0.8)

  # `exposure` is not an argument of ps_trim(), so a message naming it sends
  # the caller after something they cannot supply.
  cnd_cr <- rlang::catch_cnd(ps_trim(ps, method = "cr"), classes = "error")
  expect_s3_class(cnd_cr, "propensity_missing_arg_error")
  expect_match(conditionMessage(cnd_cr), "`.exposure`", fixed = TRUE)

  cnd_pref <- rlang::catch_cnd(ps_trim(ps, method = "pref"), classes = "error")
  expect_s3_class(cnd_pref, "propensity_missing_arg_error")
  expect_match(conditionMessage(cnd_pref), "`.exposure`", fixed = TRUE)

  expect_propensity_error(ps_trim(ps, method = "cr"))
})

# Missing values ------------------------------------------------------------

# A propensity score that arrives missing is not one this package removed, so
# the record says nothing about it: it belongs to neither the retained nor the
# trimmed positions, and it propagates as `NA` into the result. The methods
# that read a cutoff off the scores read it off the scores they have.

test_that("ps_trim() does not record a score that arrived missing as trimmed", {
  trimmed <- ps_trim(
    c(0.2, 0.5, NA, 0.7),
    method = "ps",
    lower = 0.3,
    upper = 0.9
  )
  meta <- ps_trim_meta(trimmed)

  expect_equal(as.numeric(trimmed), c(NA, 0.5, NA, 0.7))
  expect_equal(meta$trimmed_idx, 1)
  expect_equal(meta$keep_idx, c(2, 4))
  expect_equal(meta$n_obs, 4)
  expect_equal(is_unit_trimmed(trimmed), c(TRUE, FALSE, FALSE, FALSE))
})

test_that("ps_trim() takes its adaptive cutoff from the complete scores", {
  ps <- c(0.01, 0.2, 0.5, NA, 0.7, 0.99)

  # The cutoff comes from the mean and the maximum of 1 / (e (1 - e)), both of
  # which are missing the moment one score is, and the comparison that chooses
  # between the two branches is then an error rather than an answer.
  with_na <- ps_trim(ps, method = "adaptive")
  without_na <- ps_trim(ps[-4], method = "adaptive")
  meta <- ps_trim_meta(with_na)

  expect_equal(meta$cutoff, ps_trim_meta(without_na)$cutoff)
  expect_equal(as.numeric(with_na)[-4], as.numeric(without_na))
  expect_true(is.na(as.numeric(with_na)[4]))
  expect_equal(meta$keep_idx, c(2, 3, 5))
  expect_equal(meta$trimmed_idx, c(1, 6))
  expect_false(is_unit_trimmed(with_na)[4])
})

test_that("ps_trim() takes its percentile cutoffs from the complete scores", {
  ps <- c(0.05, 0.2, 0.5, NA, 0.7, 0.95)

  # `quantile()` refuses a missing value unless it is told to drop it, so the
  # percentile method is an error on any sample with one.
  with_na <- ps_trim(ps, method = "pctl", lower = 0.2, upper = 0.8)
  without_na <- ps_trim(ps[-4], method = "pctl", lower = 0.2, upper = 0.8)
  meta <- ps_trim_meta(with_na)

  expect_equal(meta$q_lower, ps_trim_meta(without_na)$q_lower)
  expect_equal(meta$q_upper, ps_trim_meta(without_na)$q_upper)
  expect_equal(as.numeric(with_na)[-4], as.numeric(without_na))
  expect_true(is.na(as.numeric(with_na)[4]))
  expect_equal(meta$keep_idx, c(2, 3, 5))
  expect_equal(meta$trimmed_idx, c(1, 6))
  expect_false(is_unit_trimmed(with_na)[4])
})

test_that("ps_trim() takes its common range from the complete scores", {
  ps <- c(0.05, 0.2, 0.5, NA, 0.7, 0.95)
  z <- c(0, 0, 1, 1, 1, 0)

  # One missing score in the focal group makes the lower bound missing, which
  # puts every score outside the range and trims the whole sample.
  with_na <- ps_trim(ps, method = "cr", .exposure = z, .focal_level = 1)
  without_na <- ps_trim(
    ps[-4],
    method = "cr",
    .exposure = z[-4],
    .focal_level = 1
  )
  meta <- ps_trim_meta(with_na)

  expect_equal(meta$cr_lower, ps_trim_meta(without_na)$cr_lower)
  expect_equal(meta$cr_upper, ps_trim_meta(without_na)$cr_upper)
  expect_equal(as.numeric(with_na)[-4], as.numeric(without_na))
  expect_true(is.na(as.numeric(with_na)[4]))
  expect_equal(meta$keep_idx, c(3, 5, 6))
  expect_equal(meta$trimmed_idx, c(1, 2))
  expect_false(is_unit_trimmed(with_na)[4])
})

test_that("ps_trim() leaves a missing score out of the preference record", {
  ps <- c(0.05, 0.2, 0.5, NA, 0.7, 0.95)
  z <- c(0, 0, 0, 1, 1, 1)

  # The preference score of a missing propensity score is missing, so the unit
  # falls outside the bounds and is recorded as one this package trimmed.
  trimmed <- ps_trim(ps, method = "pref", .exposure = z, .focal_level = 1)
  meta <- ps_trim_meta(trimmed)

  expect_true(is.na(as.numeric(trimmed)[4]))
  expect_equal(meta$keep_idx, c(3, 5))
  expect_equal(meta$trimmed_idx, c(1, 2, 6))
  expect_false(is_unit_trimmed(trimmed)[4])
})

test_that("ps_trim() refuses an exposure with missing values", {
  ps <- c(0.1, 0.2, 0.6, 0.7, 0.8, 0.9)
  z <- c(0, 0, 0, 1, 1, NA)
  z_factor <- factor(c("a", "a", "a", "b", "b", NA))

  # The proportion exposed is missing with one exposure value missing, so every
  # preference score is missing and the whole sample is trimmed without a word.
  # A trimming rule that cannot be computed is reported rather than applied.
  pref_numeric <- rlang::catch_cnd(
    ps_trim(ps, method = "pref", .exposure = z, .focal_level = 1),
    classes = "condition"
  )
  expect_s3_class(pref_numeric, "error")
  expect_s3_class(pref_numeric, "propensity_error")

  pref_factor <- rlang::catch_cnd(
    ps_trim(ps, method = "pref", .exposure = z_factor, .focal_level = "b"),
    classes = "condition"
  )
  expect_s3_class(pref_factor, "error")
  expect_s3_class(pref_factor, "propensity_error")

  # The common range is bounded by the extremes of each group, both of which
  # are missing once a unit belongs to neither, and every score then falls
  # outside the range.
  cr_numeric <- rlang::catch_cnd(
    ps_trim(ps, method = "cr", .exposure = z, .focal_level = 1),
    classes = "condition"
  )
  expect_s3_class(cr_numeric, "error")
  expect_s3_class(cr_numeric, "propensity_error")

  cr_factor <- rlang::catch_cnd(
    ps_trim(ps, method = "cr", .exposure = z_factor, .focal_level = "b"),
    classes = "condition"
  )
  expect_s3_class(cr_factor, "error")
  expect_s3_class(cr_factor, "propensity_error")

  expect_propensity_error(
    ps_trim(ps, method = "pref", .exposure = z, .focal_level = 1)
  )
})

# Bounds validation ---------------------------------------------------------

test_that("ps_trim() requires lower below upper for the pctl and pref methods", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  z <- c(0, 0, 1, 1)

  # Bounds the wrong way around describe an empty interval, so every unit falls
  # outside it and the whole sample is trimmed. Method "ps" already refuses
  # this and the other bounded methods owe the same refusal.
  pctl <- rlang::catch_cnd(
    ps_trim(ps, method = "pctl", lower = 0.9, upper = 0.1),
    classes = "condition"
  )
  expect_s3_class(pctl, "error")
  expect_s3_class(pctl, "propensity_range_error")

  pref <- rlang::catch_cnd(
    ps_trim(
      ps,
      method = "pref",
      lower = 0.9,
      upper = 0.1,
      .exposure = z,
      .focal_level = 1
    ),
    classes = "condition"
  )
  expect_s3_class(pref, "error")
  expect_s3_class(pref, "propensity_range_error")

  expect_propensity_error(ps_trim(
    ps,
    method = "pctl",
    lower = 0.9,
    upper = 0.1
  ))
})

test_that("ps_trim() refuses percentile bounds outside the unit interval", {
  ps <- c(0.2, 0.4, 0.6, 0.8)

  # For the percentile method the bounds are probabilities. `quantile()`
  # refuses one outside [0, 1] in a bare error naming `probs`, an argument
  # `ps_trim()` does not have.
  low <- rlang::catch_cnd(
    ps_trim(ps, method = "pctl", lower = -0.1),
    classes = "condition"
  )
  expect_s3_class(low, "error")
  expect_s3_class(low, "propensity_error")

  high <- rlang::catch_cnd(
    ps_trim(ps, method = "pctl", upper = 1.5),
    classes = "condition"
  )
  expect_s3_class(high, "error")
  expect_s3_class(high, "propensity_error")

  expect_propensity_error(ps_trim(ps, method = "pctl", lower = -0.1))
})

test_that("ps_trim() refuses a bound that is missing", {
  ps <- c(0.2, 0.4, 0.6, 0.8)

  # A missing bound answers neither `TRUE` nor `FALSE` in the comparison that
  # decides which scores to keep, so the rule comes out as a bare base error
  # about a missing value rather than as an answer.
  fixed <- rlang::catch_cnd(
    ps_trim(ps, method = "ps", lower = NA),
    classes = "condition"
  )
  expect_s3_class(fixed, "error")
  expect_s3_class(fixed, "propensity_error")

  pctl <- rlang::catch_cnd(
    ps_trim(ps, method = "pctl", upper = NA),
    classes = "condition"
  )
  expect_s3_class(pctl, "error")
  expect_s3_class(pctl, "propensity_error")

  pref <- rlang::catch_cnd(
    ps_trim(
      ps,
      method = "pref",
      lower = NA,
      .exposure = c(0, 0, 1, 1),
      .focal_level = 1
    ),
    classes = "condition"
  )
  expect_s3_class(pref, "error")
  expect_s3_class(pref, "propensity_error")

  expect_propensity_error(ps_trim(ps, method = "ps", lower = NA))
})

# Combining trimmed propensity scores ----------------------------------------

# Two `ps_trim` objects are combined through the prototype they share. The
# prototype stands for the trimming that produced them, so it carries the method
# and the cutoffs that method settled on. A cutoff read off the scores cannot be
# worked out again from a prototype that holds none, so it is carried across
# rather than recomputed. The prototype describes no observations, so it names no
# positions and the combined result has no record.

trim_combine_fixture <- function() {
  set.seed(4)
  n <- 40
  x <- rnorm(n, sd = 2)

  list(
    ps = plogis(x),
    exposure = rbinom(n, 1, plogis(x))
  )
}

test_that("combining adaptive ps_trim objects keeps the cutoff the trimming found", {
  fixture <- trim_combine_fixture()
  trimmed <- ps_trim(fixture$ps, method = "adaptive")
  meta <- ps_trim_meta(trimmed)

  combined <- expect_silent(c(trimmed[1:20], trimmed[21:40]))
  combined_meta <- ps_trim_meta(combined)

  expect_s3_class(combined, "ps_trim")
  expect_length(combined, 40)
  expect_equal(as.numeric(combined), as.numeric(trimmed))
  expect_equal(combined_meta$method, "adaptive")
  expect_gt(meta$cutoff, 0)
  expect_equal(combined_meta$cutoff, meta$cutoff)
  expect_null(combined_meta$keep_idx)
  expect_null(combined_meta$trimmed_idx)
  expect_null(combined_meta$n_obs)
})

test_that("combining pctl ps_trim objects keeps the quantiles the trimming found", {
  fixture <- trim_combine_fixture()
  trimmed <- ps_trim(fixture$ps, method = "pctl")
  meta <- ps_trim_meta(trimmed)

  combined <- expect_silent(c(trimmed[1:20], trimmed[21:40]))
  combined_meta <- ps_trim_meta(combined)

  expect_s3_class(combined, "ps_trim")
  expect_length(combined, 40)
  expect_equal(as.numeric(combined), as.numeric(trimmed))
  expect_equal(combined_meta$method, "pctl")
  expect_equal(combined_meta$lower, 0.05)
  expect_equal(combined_meta$upper, 0.95)

  # The quantiles come from the scores, and a prototype holds none, so a
  # prototype that works them out again reports a missing cutoff instead.
  expect_false(is.na(combined_meta$q_lower))
  expect_false(is.na(combined_meta$q_upper))
  expect_equal(combined_meta$q_lower, meta$q_lower)
  expect_equal(combined_meta$q_upper, meta$q_upper)
  expect_null(combined_meta$keep_idx)
  expect_null(combined_meta$trimmed_idx)
  expect_null(combined_meta$n_obs)
})

test_that("combining pctl ps_trim objects trimmed at different quantiles warns", {
  # Same method and same percentiles, different scores, so the quantiles those
  # percentiles landed on differ. A shared prototype would report one object's
  # cutoffs over the other object's units, so there is no shared type to combine
  # in.
  x <- ps_trim(c(0.1, 0.3, 0.5, 0.7, 0.9), method = "pctl")
  y <- ps_trim(c(0.2, 0.4, 0.5, 0.6, 0.8), method = "pctl")

  expect_false(identical(
    ps_trim_meta(x)$q_lower,
    ps_trim_meta(y)$q_lower
  ))

  combined <- expect_propensity_warning(vec_c(x, y))

  expect_type(combined, "double")
  expect_false(inherits(combined, "ps_trim"))
  expect_equal(combined, c(as.numeric(x), as.numeric(y)))
})

test_that("combining pref ps_trim objects keeps the exposure prevalence", {
  fixture <- trim_combine_fixture()
  trimmed <- ps_trim(
    fixture$ps,
    method = "pref",
    .exposure = fixture$exposure
  )
  meta <- ps_trim_meta(trimmed)

  # The preference scale is defined against the exposure, which a prototype
  # built by trimming again would have to be handed and is not.
  combined <- expect_silent(c(trimmed[1:20], trimmed[21:40]))
  combined_meta <- ps_trim_meta(combined)

  expect_s3_class(combined, "ps_trim")
  expect_length(combined, 40)
  expect_equal(as.numeric(combined), as.numeric(trimmed))
  expect_equal(combined_meta$method, "pref")
  expect_equal(combined_meta$lower, 0.3)
  expect_equal(combined_meta$upper, 0.7)
  expect_equal(combined_meta$P, meta$P)
  expect_null(combined_meta$keep_idx)
  expect_null(combined_meta$trimmed_idx)
  expect_null(combined_meta$n_obs)
})

test_that("combining cr ps_trim objects keeps the common range", {
  fixture <- trim_combine_fixture()
  trimmed <- ps_trim(fixture$ps, method = "cr", .exposure = fixture$exposure)
  meta <- ps_trim_meta(trimmed)

  combined <- expect_silent(c(trimmed[1:20], trimmed[21:40]))
  combined_meta <- ps_trim_meta(combined)

  expect_s3_class(combined, "ps_trim")
  expect_length(combined, 40)
  expect_equal(as.numeric(combined), as.numeric(trimmed))
  expect_equal(combined_meta$method, "cr")
  expect_equal(combined_meta$cr_lower, meta$cr_lower)
  expect_equal(combined_meta$cr_upper, meta$cr_upper)
  expect_null(combined_meta$keep_idx)
  expect_null(combined_meta$trimmed_idx)
  expect_null(combined_meta$n_obs)
})

test_that("combining refit ps_trim objects keeps the refit flag", {
  set.seed(5)
  n <- 40
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  model <- glm(z ~ x, family = binomial)

  trimmed <- ps_trim(unname(fitted(model)), lower = 0.2, upper = 0.8)
  refit <- ps_refit(trimmed, model)
  expect_true(is_refit(refit))

  # Refitting is part of what produced these scores rather than a position among
  # them, so it means the same thing at any length.
  combined <- expect_silent(c(refit[1:20], refit[21:40]))

  expect_s3_class(combined, "ps_trim")
  expect_length(combined, 40)
  expect_true(is_refit(combined))
})

test_that("combining a ps_trim with an integer keeps the propensity scores", {
  x <- ps_trim(c(0.2, 0.5, 0.85), method = "ps", lower = 0.1, upper = 0.9)

  # Propensity scores lie strictly between 0 and 1, so a combination that meets
  # an integer in the integers is every score rounded away.
  combined <- expect_propensity_warning(vec_c(x, 1L))

  expect_type(combined, "double")
  expect_equal(combined, c(0.2, 0.5, 0.85, 1))
})

test_that("casting a ps_trim to integer refuses rather than rounds", {
  x <- ps_trim(c(0.2, 0.5, 0.85), method = "ps", lower = 0.1, upper = 0.9)

  expect_error(
    vec_cast(x, integer()),
    class = "vctrs_error_cast_lossy"
  )
})

# The description a trimmed vector is printed under --------------------------

# `vec_ptype_full()` writes the line a trimmed vector is printed under. It
# reports how many observations were trimmed, which means something to a reader
# only against the number of observations there were, so it names both.

trim_ptype_fixture <- function() {
  set.seed(123)
  ps_trim(runif(20, 0.05, 0.95), method = "ps", lower = 0.2, upper = 0.8)
}

test_that("vec_ptype_full() names the size the trimmed count is out of", {
  trimmed <- trim_ptype_fixture()
  n_trimmed <- length(ps_trim_meta(trimmed)$trimmed_idx)

  expect_gt(n_trimmed, 0)
  expect_identical(
    vec_ptype_full(trimmed),
    sprintf("ps_trim; trimmed %d of %d", n_trimmed, length(trimmed))
  )
})

test_that("vec_ptype_full() reports the size when nothing was trimmed", {
  set.seed(123)
  untrimmed <- ps_trim(
    runif(20, 0.05, 0.95),
    method = "ps",
    lower = 0,
    upper = 1
  )

  expect_length(ps_trim_meta(untrimmed)$trimmed_idx, 0)
  expect_identical(vec_ptype_full(untrimmed), "ps_trim; trimmed 0 of 20")
})

test_that("vec_ptype_full() counts against the subset, not the whole vector", {
  trimmed <- trim_ptype_fixture()
  sliced <- trimmed[1:8]
  n_sliced <- length(ps_trim_meta(sliced)$trimmed_idx)

  # The record is re-indexed to the subset, so both the count and the size the
  # description reports are the subset's own.
  expect_gt(n_sliced, 0)
  expect_lt(n_sliced, length(ps_trim_meta(trimmed)$trimmed_idx))
  expect_identical(
    vec_ptype_full(sliced),
    sprintf("ps_trim; trimmed %d of 8", n_sliced)
  )
})

# Choosing a column from a data frame of propensity scores -------------------

# Predictions from a binary model arrive as a column per level, and trimming
# works on one column. The convention is the second column of a pair, which is
# the probability of the second level in the layout these predictions come in,
# and the only column otherwise. The caller did not make that choice, so it is
# announced rather than left to be inferred from the result.

trim_frame_fixture <- function() {
  set.seed(2)
  n <- 40
  x <- rnorm(n)
  p <- plogis(x)

  list(
    ps = data.frame(.pred_0 = 1 - p, .pred_1 = p),
    exposure = rbinom(n, 1, p)
  )
}

test_that("ps_trim() names the column it took from a data frame of two", {
  withr::local_options(propensity.quiet = FALSE)
  fixture <- trim_frame_fixture()

  messages <- testthat::capture_messages(
    ps_trim(
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

test_that("ps_trim() names the only column of a one column data frame", {
  withr::local_options(propensity.quiet = FALSE)
  fixture <- trim_frame_fixture()
  one_column <- fixture$ps[, ".pred_1", drop = FALSE]

  messages <- testthat::capture_messages(
    ps_trim(one_column, method = "ps", lower = 0.3, upper = 0.7)
  )

  expect_length(messages, 1)
  expect_true(grepl(".pred_1", messages, fixed = TRUE))
})

test_that("ps_trim() announces no column when the messages are quieted", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- trim_frame_fixture()

  expect_silent(
    ps_trim(
      fixture$ps,
      .exposure = fixture$exposure,
      method = "ps",
      lower = 0.3,
      upper = 0.7
    )
  )
})

test_that("the column ps_trim() announces is named in a full sentence", {
  withr::local_options(propensity.quiet = FALSE)
  fixture <- trim_frame_fixture()

  expect_snapshot(
    trimmed <- ps_trim(fixture$ps, method = "ps", lower = 0.3, upper = 0.7)
  )
})

test_that("ps_trim() takes the second column of a data frame of two", {
  fixture <- trim_frame_fixture()

  from_frame <- ps_trim(fixture$ps, method = "ps", lower = 0.3, upper = 0.7)
  from_second <- ps_trim(
    fixture$ps[[2]],
    method = "ps",
    lower = 0.3,
    upper = 0.7
  )
  from_first <- ps_trim(
    fixture$ps[[1]],
    method = "ps",
    lower = 0.3,
    upper = 0.7
  )

  expect_gt(length(ps_trim_meta(from_second)$trimmed_idx), 0)
  expect_equal(as.numeric(from_frame), as.numeric(from_second))
  expect_equal(
    ps_trim_meta(from_frame)$trimmed_idx,
    ps_trim_meta(from_second)$trimmed_idx
  )
  expect_false(identical(as.numeric(from_frame), as.numeric(from_first)))
})

# What a cast between trimmed vectors owes its target ------------------------

# A cast returns the values it was handed in the type it was handed, and a
# `ps_trim`'s type is the whole description of the trimming: the method, what
# it was asked for, and the cutoffs it worked out from the scores it was given.
# Two objects that disagree on any of that are not each other's type, so the
# cast has no result to give and refuses, which is the comparison
# `vec_ptype2()` already makes when it refuses to find a common type. A cast
# that compares less than the combine does hands `x` back describing itself
# under the target's name.

test_that("casting between ps_trim objects trimmed at different cutoffs refuses", {
  # The same percentiles asked of different scores are different cutoffs, and
  # the record keeps both what was asked for and what came back.
  x <- ps_trim(
    c(0.1, 0.2, 0.5, 0.8, 0.9),
    method = "pctl",
    lower = 0.2,
    upper = 0.8
  )
  to <- ps_trim(
    c(0.05, 0.3, 0.5, 0.7, 0.95),
    method = "pctl",
    lower = 0.2,
    upper = 0.8
  )

  expect_false(identical(ps_trim_meta(x)$q_lower, ps_trim_meta(to)$q_lower))
  expect_warning(
    vec_ptype2(x, to),
    class = "propensity_coercion_warning"
  )
  expect_identical(suppressWarnings(vec_ptype2(x, to)), double())

  # Nothing a `ps_trim` is printed as names a cutoff, so the two types render
  # identically and the refusal names what disagrees, as the combine does.
  expect_identical(vec_ptype_full(x), vec_ptype_full(to))
  expect_error(
    vec_cast(x, to = to),
    regexp = "different trimming parameters",
    class = "vctrs_error_incompatible_type"
  )
})

test_that("casting between ps_trim objects with different refit status refuses", {
  # Refitting is part of what produced these scores, so scores the model was
  # refit on and scores it was not are different types at the same cutoffs.
  set.seed(11)
  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  model <- glm(z ~ x, family = binomial)

  trimmed <- ps_trim(unname(fitted(model)), lower = 0.2, upper = 0.8)
  refit <- ps_refit(trimmed, model)

  expect_false(is_refit(trimmed))
  expect_true(is_refit(refit))
  expect_warning(
    vec_ptype2(trimmed, refit),
    class = "propensity_coercion_warning"
  )
  expect_identical(suppressWarnings(vec_ptype2(trimmed, refit)), double())
  expect_error(
    vec_cast(trimmed, to = refit),
    regexp = "different refit status",
    class = "vctrs_error_incompatible_type"
  )
  expect_error(
    vec_cast(refit, to = trimmed),
    regexp = "different refit status",
    class = "vctrs_error_incompatible_type"
  )
})

test_that("casting between ps_trim objects describing the same trimming succeeds", {
  # The positional half of the record describes the values arriving rather than
  # the type they arrive in, so two objects trimmed the same way are each
  # other's type however many units either one kept.
  x <- ps_trim(
    c(0.05, 0.2, 0.5, 0.8, 0.95),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )
  to <- ps_trim(
    c(0.15, 0.25, 0.55, 0.85),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )

  out <- expect_silent(vec_cast(x, to = to))

  expect_s3_class(out, "ps_trim")
  expect_equal(as.numeric(out), as.numeric(x))
  expect_equal(ps_trim_meta(out)$lower, 0.1)
  expect_equal(ps_trim_meta(out)$upper, 0.9)
})

# Refitting a model whose formula transforms its terms ------------------------

# `ps_refit()` promises the propensity model re-estimated on the rows the
# trimming kept, so that is what these tests compare it against: the same
# formula fit by hand on those rows. A model frame names its columns after term
# expressions rather than after the variables the terms are built from, so it
# carries a column called `log(x)` and no column called `x`, and a refit that
# hands the model frame back as data leaves the formula with nothing to
# evaluate. The transformations below are the ones a propensity model reaches
# for most: a pointwise transformation, a basis whose knots depend on the rows
# it is fit on, and an interaction between two transformed terms.

refit_transform_fixture <- function() {
  set.seed(2024)
  n <- 120
  x <- runif(n, 0.2, 12)
  data.frame(
    z = rbinom(n, 1, plogis(-2.2 + 1.6 * log(x))),
    x = x,
    w = rnorm(n)
  )
}

# The scores `ps_refit()` owes the caller: the model fit on the retained rows,
# read back over the full sample with the trimmed positions left missing.
refit_by_hand <- function(model_formula, data, keep_idx, n_obs) {
  fit <- glm(
    model_formula,
    data = data[keep_idx, , drop = FALSE],
    family = binomial
  )

  expected <- rep(NA_real_, n_obs)
  expected[keep_idx] <- unname(predict(fit, type = "response"))
  expected
}

test_that("ps_refit() refits a log-transformed term with no data argument", {
  df <- refit_transform_fixture()
  z <- df$z
  x <- df$x

  fit <- glm(z ~ log(x), family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  expect_equal(
    as.numeric(refit),
    refit_by_hand(z ~ log(x), df, meta$keep_idx, nrow(df))
  )
})

test_that("ps_refit() refits a spline basis with no data argument", {
  df <- refit_transform_fixture()
  z <- df$z
  x <- df$x

  fit <- glm(z ~ splines::ns(x, 2), family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  # A spline places its knots from the rows it is fit on, so refitting on the
  # retained rows moves them. The refit owes the caller the basis the retained
  # rows imply, not the one the trimmed-away rows helped position.
  expect_equal(
    as.numeric(refit),
    refit_by_hand(z ~ splines::ns(x, 2), df, meta$keep_idx, nrow(df))
  )
})

test_that("ps_refit() refits an interaction between transformed terms", {
  df <- refit_transform_fixture()
  z <- df$z
  x <- df$x
  w <- df$w

  fit <- glm(z ~ log(x) * I(w^2), family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  expect_equal(
    as.numeric(refit),
    refit_by_hand(z ~ log(x) * I(w^2), df, meta$keep_idx, nrow(df))
  )
})

test_that("ps_refit() refits a transformed term fit with a data argument", {
  df <- refit_transform_fixture()

  fit <- glm(z ~ log(x), data = df, family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  expect_equal(
    as.numeric(refit),
    refit_by_hand(z ~ log(x), df, meta$keep_idx, nrow(df))
  )
})

test_that("ps_refit() records the refit and keeps the record through a transform", {
  df <- refit_transform_fixture()
  z <- df$z
  x <- df$x

  fit <- glm(z ~ log(x), family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)

  refit <- ps_refit(trimmed, fit)
  refit_meta <- ps_trim_meta(refit)

  expect_s3_class(refit, "ps_trim")
  expect_true(is_refit(refit))
  expect_false(is_refit(trimmed))
  expect_equal(refit_meta$method, meta$method)
  expect_equal(refit_meta$lower, meta$lower)
  expect_equal(refit_meta$upper, meta$upper)
  expect_equal(refit_meta$keep_idx, meta$keep_idx)
  expect_equal(refit_meta$trimmed_idx, meta$trimmed_idx)
  expect_equal(refit_meta$n_obs, meta$n_obs)

  # Refitting re-estimates the retained units and says nothing new about the
  # trimmed ones, so the missing positions are exactly the trimmed ones.
  expect_equal(which(is.na(as.numeric(refit))), meta$trimmed_idx)
  expect_length(refit, nrow(df))
})

test_that("ps_refit() still refits an untransformed formula from the model alone", {
  df <- refit_transform_fixture()
  z <- df$z
  x <- df$x

  fit <- glm(z ~ x, family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  expect_equal(
    as.numeric(refit),
    refit_by_hand(z ~ x, df, meta$keep_idx, nrow(df))
  )
})

test_that("ps_refit() still refits a transformed term from an explicit .data", {
  df <- refit_transform_fixture()

  fit <- glm(z ~ log(x), data = df, family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)

  refit <- ps_refit(trimmed, fit, .data = df)

  expect_equal(
    as.numeric(refit),
    refit_by_hand(z ~ log(x), df, meta$keep_idx, nrow(df))
  )
})

test_that("ps_refit() still refits a model whose rows were dropped as missing", {
  df <- refit_transform_fixture()
  z <- df$z
  covariate <- df$x
  covariate[c(3, 17)] <- NA

  fit <- glm(z ~ covariate, family = binomial)

  # Fitting drops the incomplete rows, so the scores, and with them the
  # trimming record, are about the rows the model kept rather than every row
  # the variables carry.
  analysis <- data.frame(z = z, covariate = covariate)
  analysis <- analysis[!is.na(covariate), , drop = FALSE]
  expect_equal(nrow(analysis), nrow(df) - 2L)

  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_length(trimmed, nrow(analysis))
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  expect_equal(
    as.numeric(refit),
    refit_by_hand(z ~ covariate, analysis, meta$keep_idx, nrow(analysis))
  )
})

test_that("ps_refit() refits on the retained rows of a model fit with a subset", {
  df <- refit_transform_fixture()
  z <- df$z
  x <- df$x

  fit <- glm(z ~ x, family = binomial, subset = -(1:40))

  # The original subset already chose the sample the scores are about, so the
  # trimming record indexes that sample and not the rows the variables carry.
  # Refitting narrows that sample to the retained rows; it must not put the
  # original subset to work a second time on rows it was never about. A
  # positional subset re-applied to the retained rows quietly names different
  # units rather than failing, so the wrong answer would arrive without a word.
  analysis <- df[-(1:40), , drop = FALSE]

  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_length(trimmed, nrow(analysis))
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  expect_equal(
    as.numeric(refit),
    refit_by_hand(z ~ x, analysis, meta$keep_idx, nrow(analysis))
  )
})

test_that("ps_refit() puts a subset to work once when given an explicit .data", {
  df <- refit_transform_fixture()

  fit <- glm(z ~ x, data = df, family = binomial, subset = -(1:40))

  # An explicit `.data` is already the sample the scores are about, so the
  # original subset has done its work and must not choose rows a second time.
  # A positional subset re-applied to the retained rows quietly names different
  # units rather than failing, so the wrong answer arrives without a word.
  analysis <- df[-(1:40), , drop = FALSE]

  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_length(trimmed, nrow(analysis))
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit, .data = analysis)

  expect_equal(
    as.numeric(refit),
    refit_by_hand(z ~ x, analysis, meta$keep_idx, nrow(analysis))
  )
})

test_that("ps_refit() says so when the model's data can no longer be reached", {
  df <- refit_transform_fixture()

  fit <- glm(z ~ log(x), data = df, family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  # A transformed term is stored already computed, so the raw variables have to
  # be read back from the data the model names. Nothing else can stand in for
  # them, and the data are gone.
  rm(df)

  expect_error(
    ps_refit(trimmed, fit),
    class = "propensity_no_data_error"
  )
})

test_that("ps_refit() honors a subset the caller passes through ...", {
  df <- refit_transform_fixture()

  fit <- glm(z ~ x, data = df, family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  narrowed <- seq_len(20)

  refit <- ps_refit(trimmed, fit, subset = narrowed)

  # The caller's subset chooses the rows the model is fit on; the scores are
  # still read back over every retained row.
  by_hand <- glm(
    z ~ x,
    data = df[meta$keep_idx, , drop = FALSE],
    family = binomial,
    subset = narrowed
  )
  expected <- rep(NA_real_, nrow(df))
  expected[meta$keep_idx] <- unname(predict(
    by_hand,
    newdata = df[meta$keep_idx, , drop = FALSE],
    type = "response"
  ))

  expect_equal(as.numeric(refit), expected)
})

test_that("ps_refit() says so when the model's data no longer hold its rows", {
  df <- refit_transform_fixture()

  fit <- glm(z ~ log(x), data = df, family = binomial)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  # The name the model gave its data now holds a narrower frame, so the rows the
  # model analyzed cannot all be found in it. Reading whatever rows are there
  # would refit on a sample the trimming record was never about.
  df <- df[41:nrow(df), , drop = FALSE]

  expect_error(
    ps_refit(trimmed, fit),
    class = "propensity_no_data_error"
  )
})

test_that("ps_refit() refits a weighted model whose weights live in the data", {
  df <- refit_transform_fixture()
  df$sampling_wt <- runif(nrow(df), 0.5, 2)

  fit <- glm(
    z ~ log(x),
    data = df,
    family = quasibinomial,
    weights = sampling_wt
  )
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  # `weights` names a column, so the recovered data have to carry it for the
  # weights to follow the retained rows.
  kept <- df[meta$keep_idx, , drop = FALSE]
  by_hand <- glm(
    z ~ log(x),
    data = kept,
    family = quasibinomial,
    weights = sampling_wt
  )
  expected <- rep(NA_real_, nrow(df))
  expected[meta$keep_idx] <- unname(predict(by_hand, type = "response"))

  expect_equal(as.numeric(refit), expected)
})

test_that("ps_refit() refits a model whose offset lives in the data", {
  df <- refit_transform_fixture()
  df$log_time <- rnorm(nrow(df), 0, 0.1)

  fit <- glm(z ~ log(x), data = df, family = binomial, offset = log_time)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  kept <- df[meta$keep_idx, , drop = FALSE]
  by_hand <- glm(
    z ~ log(x),
    data = kept,
    family = binomial,
    offset = log_time
  )
  expected <- rep(NA_real_, nrow(df))
  expected[meta$keep_idx] <- unname(predict(by_hand, type = "response"))

  expect_equal(as.numeric(refit), expected)
})

test_that("ps_refit() refits weights named in the data of an untransformed fit", {
  df <- refit_transform_fixture()
  df$sampling_wt <- runif(nrow(df), 0.5, 2)

  # Nothing in the formula is transformed, so the stored model frame holds every
  # variable the formula names; it still does not hold the weights.
  fit <- glm(z ~ x, data = df, family = quasibinomial, weights = sampling_wt)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_gt(length(meta$trimmed_idx), 0)

  refit <- ps_refit(trimmed, fit)

  kept <- df[meta$keep_idx, , drop = FALSE]
  by_hand <- glm(
    z ~ x,
    data = kept,
    family = quasibinomial,
    weights = sampling_wt
  )
  expected <- rep(NA_real_, nrow(df))
  expected[meta$keep_idx] <- unname(predict(by_hand, type = "response"))

  expect_equal(as.numeric(refit), expected)
})

test_that("ps_refit() names padded predictions when the row counts disagree", {
  df <- refit_transform_fixture()
  df$x[c(3, 17)] <- NA

  # `na.exclude` puts the dropped rows back as `NA` when predicting, so the
  # scores are about a longer sample than the fit ever read, and the record
  # built from them cannot index the rows the model analyzed.
  fit <- glm(z ~ x, data = df, family = binomial, na.action = na.exclude)
  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  expect_length(trimmed, nrow(df))

  expect_error(
    ps_refit(trimmed, fit),
    class = "propensity_length_error"
  )
  expect_propensity_error(ps_refit(trimmed, fit))
})

test_that("ps_refit() blames padding only when the scores outnumber the rows", {
  df <- refit_transform_fixture()

  fit <- glm(z ~ x, data = df, family = binomial)

  # Scores describing a narrower sample than the fit read. The counts disagree
  # the other way round from padding, which only ever lengthens them, so the
  # refusal has no business naming `na.exclude`.
  narrower <- ps_trim(
    unname(predict(fit, type = "response"))[1:100],
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  cnd <- expect_error(
    ps_refit(narrower, fit),
    class = "propensity_length_error"
  )
  msg <- conditionMessage(cnd)

  expect_match(msg, "has 120 rows", fixed = TRUE)
  expect_match(msg, "has 100 observations", fixed = TRUE)
  expect_no_match(msg, "na.exclude", fixed = TRUE)
})

test_that("ps_refit() replaces a stored subset with the one the caller passes", {
  df <- refit_transform_fixture()

  fit <- glm(z ~ x, data = df, family = binomial, subset = -(1:40))
  analysis <- df[-(1:40), , drop = FALSE]

  trimmed <- ps_trim(
    unname(predict(fit, type = "response")),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trim_meta(trimmed)
  expect_length(trimmed, nrow(analysis))

  narrowed <- seq_len(30)
  refit <- ps_refit(trimmed, fit, subset = narrowed)

  # The caller's subset stands in for the stored one rather than joining it:
  # it chooses among the retained rows once, and the rows the stored subset
  # already excluded are not excluded a second time.
  kept <- analysis[meta$keep_idx, , drop = FALSE]
  by_hand <- glm(z ~ x, data = kept, family = binomial, subset = narrowed)
  expected <- rep(NA_real_, nrow(analysis))
  expected[meta$keep_idx] <- unname(predict(
    by_hand,
    newdata = kept,
    type = "response"
  ))

  expect_equal(as.numeric(refit), expected)

  # Had the stored subset also been honored, the fit would have read a
  # different set of rows and said something different about them.
  applied_twice <- glm(
    z ~ x,
    data = kept[-(1:40), , drop = FALSE],
    family = binomial,
    subset = narrowed
  )
  not_expected <- rep(NA_real_, nrow(analysis))
  not_expected[meta$keep_idx] <- unname(predict(
    applied_twice,
    newdata = kept,
    type = "response"
  ))

  expect_false(isTRUE(all.equal(as.numeric(refit), not_expected)))
})

# Naming the propensity scores `.propensity` ----------------------------------

# The weight functions read the propensity scores from `.propensity` and these
# read them from `ps`, so a call written against one is refused by the other in
# both directions. The tests below pin the scores under the new name, the
# deprecated shim that keeps the old name working for a release, and the
# refusal to read both names at once. The positional pin comes first: whatever
# else the rename moves, it must not move what a call that names nothing
# returns.

trim_rename_scores <- function() {
  c(0.05, 0.15, 0.35, 0.5, 0.65, 0.85, 0.95)
}

trim_rename_positional <- function() {
  ps_trim(trim_rename_scores(), method = "ps", lower = 0.2, upper = 0.8)
}

# A categorical propensity score matrix whose rows sum to 1, small enough that
# the trimmed rows can be read off it: rows 1 and 3 hold a score of 0.1 and are
# the two the cutoff below discards.
trim_rename_matrix <- function() {
  m <- rbind(
    c(0.70, 0.20, 0.10),
    c(0.20, 0.60, 0.20),
    c(0.10, 0.30, 0.60),
    c(0.50, 0.30, 0.20),
    c(0.30, 0.40, 0.30),
    c(0.25, 0.35, 0.40)
  )
  colnames(m) <- c("a", "b", "c")
  m
}

trim_rename_exposure <- function() {
  factor(c("a", "b", "c", "a", "b", "c"), levels = c("a", "b", "c"))
}

test_that("ps_trim() trims a positional vector of scores at the cutoffs", {
  out <- trim_rename_positional()

  expect_s3_class(out, "ps_trim")
  expect_equal(as.numeric(out), c(NA, NA, 0.35, 0.5, 0.65, NA, NA))

  meta <- ps_trim_meta(out)
  expect_equal(meta$method, "ps")
  expect_equal(meta$lower, 0.2)
  expect_equal(meta$upper, 0.8)
  expect_equal(meta$keep_idx, 3:5)
  expect_equal(meta$trimmed_idx, c(1L, 2L, 6L, 7L))
  expect_equal(meta$n_obs, 7L)
})

test_that("ps_trim() reads the propensity scores from .propensity", {
  expect_equal(
    ps_trim(
      .propensity = trim_rename_scores(),
      method = "ps",
      lower = 0.2,
      upper = 0.8
    ),
    trim_rename_positional()
  )
})

test_that("ps_trim() deprecates the propensity scores under ps", {
  with_always_deprecated({
    expect_warning(
      ps_trim(
        ps = trim_rename_scores(),
        method = "ps",
        lower = 0.2,
        upper = 0.8
      ),
      class = "lifecycle_warning_deprecated"
    )
  })

  # The old name still has to reach the same trimming, not merely warn. The
  # deprecation is pinned above, so it is silenced here rather than repeated.
  withr::local_options(lifecycle_verbosity = "quiet")
  expect_equal(
    ps_trim(ps = trim_rename_scores(), method = "ps", lower = 0.2, upper = 0.8),
    trim_rename_positional()
  )
})

test_that("ps_trim() refuses the propensity scores under both names", {
  withr::local_options(lifecycle_verbosity = "quiet")

  # The condition subclass is the shim's to choose; what this pins is that the
  # refusal is one of the package's own errors and that it names both spellings,
  # so the caller can see which one to drop.
  err <- expect_error(
    ps_trim(
      .propensity = trim_rename_scores(),
      ps = trim_rename_scores(),
      method = "ps",
      lower = 0.2,
      upper = 0.8
    ),
    class = "propensity_error"
  )

  msg <- conditionMessage(err)
  expect_match(msg, "`.propensity`", fixed = TRUE)
  expect_match(msg, "`ps`", fixed = TRUE)
})

test_that("ps_trim() dispatches on a matrix supplied as .propensity", {
  out <- ps_trim(
    .propensity = trim_rename_matrix(),
    .exposure = trim_rename_exposure(),
    method = "ps",
    lower = 0.15
  )

  expect_s3_class(out, c("ps_trim_matrix", "ps_trim", "matrix"))
  expect_equal(ps_trim_meta(out)$trimmed_idx, c(1L, 3L))
  expect_equal(
    out,
    ps_trim(
      trim_rename_matrix(),
      .exposure = trim_rename_exposure(),
      method = "ps",
      lower = 0.15
    )
  )
})

test_that("ps_trim() dispatches on a data frame supplied as .propensity", {
  out <- ps_trim(
    .propensity = as.data.frame(trim_rename_matrix()),
    .exposure = trim_rename_exposure(),
    method = "ps",
    lower = 0.15
  )

  expect_s3_class(out, c("ps_trim_matrix", "ps_trim", "matrix"))
  expect_equal(
    out,
    ps_trim(
      as.data.frame(trim_rename_matrix()),
      .exposure = trim_rename_exposure(),
      method = "ps",
      lower = 0.15
    )
  )
})

test_that("ps_trim() names .propensity in the out-of-range error", {
  err <- expect_error(
    ps_trim(.propensity = c(-1, 0.5), method = "ps", lower = 0.2, upper = 0.8),
    class = "propensity_range_error"
  )

  msg <- conditionMessage(err)
  expect_match(msg, "`.propensity`", fixed = TRUE)
  expect_false(grepl("`ps`", msg, fixed = TRUE))
})

# Bounds, groups, and arguments the trimming cannot use ------------------------

test_that("ps_trim() names the propensity scores it was not given", {
  err <- expect_error(ps_trim(), class = "propensity_missing_arg_error")
  expect_match(conditionMessage(err), "`.propensity`", fixed = TRUE)
})

test_that("the ps deprecation is attributed to the caller", {
  messages <- deprecation_warnings_from_user(
    quote(ps_trim(ps = scores, method = "ps", lower = 0.2, upper = 0.8)),
    list(scores = trim_rename_scores())
  )

  expect_length(messages, 1)
  expect_false(deprecation_misattributed(messages))
})

# One exposure group's scores are all missing, so the common range has no bound
# to read off that group.
cr_missing_group_scores <- function() {
  c(NA, NA, 0.4, 0.55, 0.7)
}

test_that("ps_trim() refuses a common range read off no scores", {
  # The lower bound is the lowest score among the focal units. Over no scores
  # `min()` returns `Inf` under a base R warning, which reads as a range no unit
  # falls inside, so every unit is trimmed instead.
  expect_no_warning(
    expect_error(
      ps_trim(
        cr_missing_group_scores(),
        method = "cr",
        .exposure = c(1, 1, 0, 0, 0)
      ),
      class = "propensity_no_data_error"
    )
  )

  # The upper bound is the highest score among the reference units, and
  # `max()` returns `-Inf` over none of them.
  expect_no_warning(
    expect_error(
      ps_trim(
        cr_missing_group_scores(),
        method = "cr",
        .exposure = c(0, 0, 1, 1, 1)
      ),
      class = "propensity_no_data_error"
    )
  )
})

test_that("ps_trim() drops the bounds an adaptive trimming ignores", {
  set.seed(11)
  ps <- runif(60, 0.02, 0.98)

  # The bound is announced as ignored, which is what the announcement is for;
  # it is muffled below so that what the metadata records is what is read.
  expect_warning(
    ps_trim(ps, method = "adaptive", lower = 0.1),
    class = "propensity_warning"
  )
  ignored <- suppressWarnings(ps_trim(ps, method = "adaptive", lower = 0.1))
  meta <- ps_trim_meta(ignored)

  # An adaptive trimming reads its cutoff off the scores, so a bound it was
  # handed took no part in it and describes nothing about the result.
  expect_null(meta$lower)
  expect_null(meta$upper)

  plain <- ps_trim(ps, method = "adaptive")
  expect_equal(meta$cutoff, ps_trim_meta(plain)$cutoff)

  # Two trimmings that ran the same rule to the same cutoff describe the same
  # trimming, so combining them keeps the class rather than reporting different
  # parameters and falling back to numeric.
  combined <- expect_no_warning(c(plain, ignored))
  expect_s3_class(combined, "ps_trim")
  expect_equal(
    as.numeric(combined),
    c(as.numeric(plain), as.numeric(ignored))
  )
})

# The binary route names ps_trim, whichever method answered ------------------

# `ps_trim.data.frame` reaches the vector method by a plain call rather than by
# dispatch, so a condition the vector method raises reports the frame it was
# raised from, which is a method no caller wrote. The same condition on a bare
# vector names `ps_trim()`, because dispatch reports the generic. The scores a
# caller holds in one column of a data frame are the scores they would have
# passed as a vector, so the two routes owe the same report.

trim_binary_frame_fixture <- function() {
  list(
    scores = data.frame(.pred_1 = c(0.15, 0.4, 0.6, 0.85)),
    out_of_range = data.frame(.pred_1 = c(0, 0.5, 1)),
    exposure = c(0, 1, 0, 1)
  )
}

test_that("a refusal on the binary data frame route names ps_trim", {
  fixture <- trim_binary_frame_fixture()

  # The scores themselves.
  expect_identical(
    condition_call_name(ps_trim(fixture$out_of_range, method = "ps")),
    "ps_trim"
  )

  # The method asked for.
  expect_identical(
    condition_call_name(ps_trim(fixture$scores, method = "bogus")),
    "ps_trim"
  )

  # The cutoffs the method was handed.
  expect_identical(
    condition_call_name(
      ps_trim(fixture$scores, method = "ps", lower = 0.9, upper = 0.1)
    ),
    "ps_trim"
  )

  # The level the binary coding is read against.
  expect_identical(
    condition_call_name(
      ps_trim(fixture$scores, method = "ps", .focal_level = c(0, 1))
    ),
    "ps_trim"
  )

  # The exposure a preference score cannot be worked out without.
  expect_identical(
    condition_call_name(ps_trim(fixture$scores, method = "pref")),
    "ps_trim"
  )
})

test_that("a warning on the binary data frame route names ps_trim", {
  fixture <- trim_binary_frame_fixture()

  # Two of the methods read their cutoffs off the scores and say so when they
  # are handed cutoffs they will not use. The announcement is the caller's to
  # act on, so it names the call the caller wrote.
  expect_warning(
    ps_trim(fixture$scores, method = "adaptive", lower = 0.2),
    class = "propensity_warning"
  )
  expect_identical(
    condition_call_name(
      ps_trim(fixture$scores, method = "adaptive", lower = 0.2),
      classes = "warning"
    ),
    "ps_trim"
  )

  expect_warning(
    ps_trim(
      fixture$scores,
      method = "cr",
      upper = 0.8,
      .exposure = fixture$exposure
    ),
    class = "propensity_warning"
  )
  expect_identical(
    condition_call_name(
      ps_trim(
        fixture$scores,
        method = "cr",
        upper = 0.8,
        .exposure = fixture$exposure
      ),
      classes = "warning"
    ),
    "ps_trim"
  )
})

test_that("a refusal on the numeric route still names ps_trim", {
  fixture <- trim_binary_frame_fixture()
  scores <- fixture$scores[[1]]

  # Dispatch reports the generic without being handed a frame, so the vector
  # route already names `ps_trim()`. Threading a frame through the data frame
  # route has no business moving that.
  expect_identical(
    condition_call_name(ps_trim(fixture$out_of_range[[1]], method = "ps")),
    "ps_trim"
  )
  expect_identical(
    condition_call_name(ps_trim(scores, method = "bogus")),
    "ps_trim"
  )
  expect_identical(
    condition_call_name(
      ps_trim(scores, method = "ps", lower = 0.9, upper = 0.1)
    ),
    "ps_trim"
  )
  expect_identical(
    condition_call_name(
      ps_trim(scores, method = "ps", .focal_level = c(0, 1))
    ),
    "ps_trim"
  )
  expect_identical(
    condition_call_name(ps_trim(scores, method = "pref")),
    "ps_trim"
  )
})

test_that("ps_trim refuses a call argument on the binary route", {
  fixture <- trim_binary_frame_fixture()

  # The generic passes its dots to its methods, so the frame the binary path
  # reports against is reachable from user code, and a value the condition
  # system cannot read as one is refused where it arrives rather than left to
  # turn the next guard that fires into a report of rlang's internals.
  expect_error(
    ps_trim(fixture$scores[[1]], method = "ps", call = "bogus"),
    class = "propensity_call_arg_error"
  )
  expect_identical(
    condition_call_name(
      ps_trim(fixture$scores[[1]], method = "ps", call = "bogus")
    ),
    "ps_trim"
  )

  # The data frame method reads the value before it hands the frame on, so it
  # is the one the refusal names.
  expect_error(
    ps_trim(fixture$scores, method = "ps", call = "bogus"),
    class = "propensity_call_arg_error"
  )
  expect_identical(
    condition_call_name(ps_trim(fixture$scores, method = "ps", call = "bogus")),
    "ps_trim"
  )
})

# ---- fitted-model methods ---------------------------------------------------

# `ps_trim()` reads a fitted propensity score model the way the weight functions
# read one: the scores come off the fit, and the exposure comes off it too for
# the methods that need one, unless the caller names an exposure of their own.
# Each test below holds the model route to the route the same scores take when
# they are extracted by hand, which is the contract the weight functions are
# held to in tests/testthat/test-categorical-models.R.

trim_model_data <- local({
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

# The formulas are written out rather than passed in, because `ps_refit()`
# re-evaluates the fitting call and a formula held in a helper's argument is not
# in scope where that re-evaluation happens.
trim_binary_fit <- function() {
  glm(z ~ x1 + x2, data = trim_model_data, family = binomial())
}

trim_categorical_fit <- function() {
  nnet::multinom(trt ~ x1 + x2, data = trim_model_data, trace = FALSE)
}

trim_two_level_fit <- function() {
  nnet::multinom(a2 ~ x1 + x2, data = trim_model_data, trace = FALSE)
}

# Values, shape, class, and record together: a model route that agreed on the
# scores but described the trimming differently would still be a different
# trimming, and the record is what the rest of the package reads.
expect_same_trim <- function(from_model, oracle) {
  testthat::expect_equal(
    as.numeric(from_model),
    as.numeric(oracle),
    tolerance = 1e-12
  )
  testthat::expect_identical(dim(from_model), dim(oracle))
  testthat::expect_identical(dimnames(from_model), dimnames(oracle))
  testthat::expect_identical(class(from_model), class(oracle))
  testthat::expect_identical(ps_trim_meta(from_model), ps_trim_meta(oracle))
}

test_that("ps_trim() trims the scores a binomial fit reports", {
  fit <- trim_binary_fit()
  scores <- predict(fit, type = "response")

  expect_same_trim(ps_trim(fit, method = "ps"), ps_trim(scores, method = "ps"))
  expect_same_trim(
    ps_trim(fit, method = "adaptive"),
    ps_trim(scores, method = "adaptive")
  )
  expect_same_trim(
    ps_trim(fit, method = "pctl"),
    ps_trim(scores, method = "pctl")
  )
})

test_that("ps_trim() reads the exposure off a fit for the methods that need one", {
  fit <- trim_binary_fit()
  scores <- predict(fit, type = "response")

  expect_same_trim(
    ps_trim(fit, method = "pref"),
    ps_trim(scores, method = "pref", .exposure = trim_model_data$z)
  )
  expect_same_trim(
    ps_trim(fit, method = "cr"),
    ps_trim(scores, method = "cr", .exposure = trim_model_data$z)
  )
})

test_that("an explicit .exposure wins over the one a fit carries in ps_trim()", {
  fit <- trim_binary_fit()
  scores <- predict(fit, type = "response")
  flipped <- 1 - trim_model_data$z

  expect_same_trim(
    ps_trim(fit, method = "pref", .exposure = flipped),
    ps_trim(scores, method = "pref", .exposure = flipped)
  )
  expect_false(identical(
    ps_trim_meta(ps_trim(fit, method = "pref", .exposure = flipped))$keep_idx,
    ps_trim_meta(ps_trim(fit, method = "pref"))$keep_idx
  ))
})

test_that("ps_trim() trims the probabilities a multinomial fit reports", {
  skip_if_not_installed("nnet")

  fit <- trim_categorical_fit()
  probs <- fitted(fit)
  trt <- trim_model_data$trt

  expect_same_trim(
    ps_trim(fit, method = "ps"),
    ps_trim(probs, method = "ps", .exposure = trt)
  )
  expect_same_trim(
    ps_trim(fit, method = "optimal"),
    ps_trim(probs, method = "optimal", .exposure = trt)
  )
})

test_that("ps_trim() matches a multinomial fit's columns to the exposure given", {
  skip_if_not_installed("nnet")

  fit <- trim_categorical_fit()
  reordered <- relevel(trim_model_data$trt, "c")
  expect_false(identical(levels(reordered), fit$lev))

  expect_same_trim(
    ps_trim(fit, method = "ps", .exposure = reordered),
    ps_trim(fitted(fit), method = "ps", .exposure = reordered)
  )
  expect_identical(
    colnames(ps_trim(fit, method = "ps", .exposure = reordered)),
    levels(reordered)
  )
})

test_that("ps_trim() reads a two-level multinomial fit on the binary path", {
  skip_if_not_installed("nnet")

  fit <- trim_two_level_fit()
  scores <- as.numeric(fitted(fit))

  expect_same_trim(ps_trim(fit, method = "ps"), ps_trim(scores, method = "ps"))
  expect_same_trim(
    ps_trim(fit, method = "pref"),
    ps_trim(scores, method = "pref", .exposure = trim_model_data$a2)
  )
})

test_that("ps_trim() refuses a fit it cannot read propensity scores from", {
  linear <- lm(z ~ x1 + x2, data = trim_model_data)

  expect_error(
    ps_trim(linear, method = "ps"),
    class = "propensity_method_error"
  )
  expect_error(
    ps_trim(structure(list(), class = "not_a_model"), method = "ps"),
    class = "propensity_method_error"
  )
})

test_that("ps_refit() refits a binomial fit trimmed through the model route", {
  fit <- trim_binary_fit()
  refitted <- ps_refit(
    ps_trim(fit, method = "ps"),
    fit,
    .data = trim_model_data
  )

  expect_s3_class(refitted, "ps_trim")
  expect_true(ps_trim_meta(refitted)$refit)
  expect_same_trim(
    refitted,
    ps_refit(
      ps_trim(predict(fit, type = "response"), method = "ps"),
      fit,
      .data = trim_model_data
    )
  )
})

test_that("ps_refit() refits a multinomial fit trimmed through the model route", {
  skip_if_not_installed("nnet")

  fit <- trim_categorical_fit()
  refitted <- ps_refit(
    ps_trim(fit, method = "ps"),
    fit,
    .data = trim_model_data
  )

  expect_s3_class(refitted, "ps_trim")
  expect_true(ps_trim_meta(refitted)$refit)
  expect_same_trim(
    refitted,
    ps_refit(
      ps_trim(fitted(fit), method = "ps", .exposure = trim_model_data$trt),
      fit,
      .data = trim_model_data
    )
  )
})

test_that("ps_trim() names the class of a fit it has no reading for", {
  expect_propensity_error(
    ps_trim(lm(z ~ x1 + x2, data = trim_model_data), method = "ps")
  )
})
