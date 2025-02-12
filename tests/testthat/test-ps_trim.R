test_that("ps_trim() - Basic structure and return types", {
  set.seed(42)

  n <- 50
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))

  # Fit a logistic regression
  fit <- glm(z ~ x, family = binomial)
  ps_vec <- predict(fit, type = "response")

  # 1) A default call with method="ps"
  out <- ps_trim(ps_vec, exposure = z, method = "ps")

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
    "ignored"
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
    "must supply a binary `exposure`"
  )

  # 2) If exposure is all 0 or all 1 => fail
  expect_error(
    ps_trim(ps, exposure = rep(0, n), method = "pref"),
    "cannot compute preference score"
  )
  expect_error(
    ps_trim(ps, exposure = rep(1, n), method = "pref"),
    "cannot compute preference score"
  )

  # 3) Valid usage
  out_pref <- ps_trim(ps, exposure = z, method = "pref")
  meta_pref <- ps_trim_meta(out_pref)
  expect_equal(meta_pref$lower, 0.3)
  expect_equal(meta_pref$upper, 0.7)
  expect_true("pref_formula" %in% names(meta_pref))

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
  expect_error(ps_trim(ps, method = "cr"), "must supply a binary `exposure`")

  # If all 0 or all 1 => fail
  expect_error(
    ps_trim(ps, exposure = rep(1, n), method = "cr"),
    "cannot compute common range"
  )

  # Valid usage
  out_cr <- ps_trim(ps, exposure = z, method = "cr")
  meta_cr <- ps_trim_meta(out_cr)
  ps_treat <- ps[z == 1]
  ps_untrt <- ps[z == 0]
  cr_l_exp <- min(ps_treat)
  cr_u_exp <- max(ps_untrt)
  expect_equal(meta_cr$cr_lower, cr_l_exp)
  expect_equal(meta_cr$cr_upper, cr_u_exp)

  # Check that user-specified lower/upper => warning
  expect_warning(
    ps_trim(ps, exposure = z, method = "cr", lower = 0.2, upper = 0.8),
    "ignored"
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
  out <- ps_trim(ps, exposure = z, method = "ps", lower = 0.2, upper = 0.8)
  # Suppose we do a normal refit
  refit_out <- ps_refit(out, model = fit)
  expect_s3_class(refit_out, "ps_trim")
  meta_r <- ps_trim_meta(refit_out)
  expect_true(isTRUE(meta_r$refit))

  # If everything is trimmed => error
  ps_edge <- c(0.01, 0.01, 0.99, 0.99)
  z_edge <- c(0, 1, 1, 0)
  out_empty <- ps_trim(ps_edge, exposure = z_edge, method = "ps", lower = 1.1, upper = 2)
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
  ps  <- predict(fit, type = "response")

  # 2) Trim the PS (e.g. method="ps" with [0.2, 0.8])
  trimmed_ps <- ps_trim(ps, exposure = z, method = "ps", lower = 0.2, upper = 0.8)
  expect_false(is_refit(trimmed_ps))  # not refit yet

  # 3) Refit on the subset
  trimmed_refit <- ps_refit(trimmed_ps, model = fit)
  expect_true(is_refit(trimmed_refit))  # now refit=TRUE in ps_trim_meta

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
  expect_true(is_refit(w_ate))  # e.g. if is_refit.psw() checks ps_trim_meta

  # The estimand should include "; trimmed"
  expect_match(estimand(w_ate), "; trimmed$")
})

