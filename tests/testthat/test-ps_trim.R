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
  expect_true(all(c("method", "lower", "upper", "keep_idx", "trimmed_idx") %in% names(meta)))

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
  expect_warning(
    out_adapt_warn <- ps_trim(ps, method = "adaptive", lower = 0.2, upper = 0.8),
    "For `method = 'adaptive'`, `lower` and `upper` are ignored."
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
  expect_error(
    ps_trim(ps, method = "pref"),
    "For `method = 'pref'`, must supply `exposure`."
  )

  # 2) If exposure is all 0 or all 1 => fail
  expect_error(
    ps_trim(ps, .exposure = rep(0, n), method = "pref"),
    class = "propensity_binary_transform_error"
  )

  expect_error(
    ps_trim(ps, .exposure = rep(1, n), method = "pref"),
    class = "propensity_binary_transform_error"
  )

  # 3) Valid usage
  out_pref <- ps_trim(ps, .exposure = z, method = "pref", .treated = 1)
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
  expect_error(
    ps_trim(ps, method = "cr"),
    "For `method = 'cr'`, must supply `exposure`."
  )

  # If all 0 or all 1 => fail
  expect_error(
    ps_trim(ps, .exposure = rep(1, n), method = "cr"),
    class = "propensity_binary_transform_error"
  )

  # Valid usage
  out_cr <- ps_trim(ps, .exposure = z, method = "cr", .treated = 1)
  meta_cr <- ps_trim_meta(out_cr)
  ps_treat <- ps[z == 1]
  ps_untrt <- ps[z == 0]
  cr_l_exp <- min(ps_treat)
  cr_u_exp <- max(ps_untrt)
  expect_equal(meta_cr$cr_lower, cr_l_exp)
  expect_equal(meta_cr$cr_upper, cr_u_exp)

  # Check that user-specified lower/upper => warning
  expect_warning(
    ps_trim(ps, .exposure = z, method = "cr", lower = 0.2, upper = 0.8, .treated = 1),
    "For `method = 'cr'`, `lower` and `upper` are ignored."
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

  # If everything is trimmed => error
  ps_edge <- c(0.01, 0.01, 0.99, 0.99)
  z_edge <- c(0, 1, 1, 0)
  out_empty <- ps_trim(ps_edge, method = "ps", lower = 1.1, upper = 2)
  expect_error(
    ps_refit(out_empty, model = fit),
    "No retained rows to refit on"
  )
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
  w_ate <- wt_ate(trimmed_refit, .exposure = z, exposure_type = "binary", .treated = 1)

  # 5) Check final psw object
  expect_s3_class(w_ate, "psw")

  # Should be trimmed, per the weighting method's logic
  expect_true(is_trimmed(w_ate))
  # Should NOT be truncated or stabilized
  expect_false(is_truncated(w_ate))
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
  out_adapt <- ps_trim(p_vec, .exposure = z, method = "adaptive", .treated = 1)

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
  expect_true(any(is.na(out_data)))
})

test_that("Check defaults for helper functions", {
  # 1) Any random object that is not ps_trim => default method => FALSE
  not_trim_obj <- c(0.2, 0.4, 0.6)
  expect_false(is_trimmed(not_trim_obj)) # triggers is_trimmed.default()
  expect_false(is_refit(not_trim_obj))

  # 2) A mock ps_trim object => method => TRUE
  # For a real test, you'd create it via ps_trim(...). Here we simulate:
  fake_ps_trim <- structure(
    c(0.5, NA, 0.7),
    class = "ps_trim"
  )
  expect_true(is_trimmed(fake_ps_trim)) # triggers is_trimmed.ps_trim()
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
      method      = "ps",
      keep_idx    = c(1, 3),
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

test_that("Arithmetic with ps_trim fails (vec_arith.*)", {
  # Create two ps_trim objects or combine with numeric
  x <- new_trimmed_ps(c(0.2, 0.3), ps_trim_meta = list())
  y <- new_trimmed_ps(c(0.4, 0.9), ps_trim_meta = list())

  # Attempt numeric ops => should fail
  expect_error(x + 1, class = "vctrs_error_incompatible_op")
  expect_error(x + 1L, class = "vctrs_error_incompatible_op")
  expect_error(1 + x, class = "vctrs_error_incompatible_op")
  expect_error(1L + x, class = "vctrs_error_incompatible_op")
  expect_error(x + list(1), class = "vctrs_error_incompatible_op")
  # Combining two ps_trim
  expect_error(x + y, class = "vctrs_error_incompatible_op")
})

test_that("Combining two ps_trim => error in vec_ptype2.ps_trim.ps_trim", {
  x <- new_trimmed_ps(
    c(0.2, 0.4),
    ps_trim_meta = list(trimmed_idx = integer(0))
  )
  y <- new_trimmed_ps(
    c(0.3, 0.5),
    ps_trim_meta = list(trimmed_idx = integer(0))
  )

  # Attempt to combine, e.g. c(x,y) => triggers vctrs' combining rules
  # or vctrs::vec_c
  expect_error(vec_c(x, y), "Can't combine two ps_trim objects")
})

test_that("Combining ps_trim with double => double", {
  x <- new_trimmed_ps(c(0.2, 0.5), ps_trim_meta = list())

  # vctrs logic => ptype2 => double
  combined <- vec_c(x, 0.7)
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

test_that("Casting double -> ps_trim => minimal ps_trim object", {
  base_vec <- c(0.1, 0.7, NA, 0.4)
  # If we do vec_cast(base_vec, ps_trim())
  # => calls vec_cast.ps_trim.double
  ps_t <- vec_cast(base_vec, to = structure(double(), class = "ps_trim"))
  expect_s3_class(ps_t, "ps_trim")
  # The meta is "unknown" method or similar
  meta <- attr(ps_t, "ps_trim_meta")
  expect_equal(meta$method, "unknown")
  expect_equal(meta$keep_idx, seq_along(base_vec))
  expect_length(meta$trimmed_idx, 0)
})

test_that("Casting integer->ps_trim likewise uses new_trimmed_ps", {
  base_int <- c(0L, 1L, 999L)
  ps_t <- vec_cast(base_int, to = structure(double(), class = "ps_trim"))
  expect_s3_class(ps_t, "ps_trim")
  # check the data is double
  expect_equal(as.numeric(ps_t), c(0, 1, 999))
})
