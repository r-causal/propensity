test_that("ps_trim() - Basic structure and return types", {
  set.seed(42)

  n <- 50
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))

  # Fit ps
  fit <- glm(z ~ x, family = binomial)
  ps_vec <- predict(fit, type = "response")

  # Test a default call
  out <- ps_trim(ps_vec, exposure = z, method = "ps")

  # Basic checks
  expect_type(out, "list")
  expect_named(out, c("trimmed_ps", "method_info", "keep_idx", "trimmed_idx"))

  expect_true(all(out$trimmed_ps >= 0 & out$trimmed_ps <= 1))
  expect_type(out$method_info, "list")

  # Indices should not overlap
  expect_length(intersect(out$keep_idx, out$trimmed_idx), 0)
})

test_that("ps method: default and custom cutoffs", {
  set.seed(1)

  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x))

  fit  <- glm(z ~ x, family = binomial)
  ps   <- predict(fit, type = "response")

  # 1) Default cutoffs (0.1, 0.9)
  out1 <- ps_trim(ps, method = "ps")
  expect_equal(out1$method_info$lower, 0.1)
  expect_equal(out1$method_info$upper, 0.9)

  # 2) Custom cutoffs
  out2 <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
  expect_equal(out2$method_info$lower, 0.2)
  expect_equal(out2$method_info$upper, 0.8)

  # Check that everything is in the correct range
  expect_true(all(out2$trimmed_ps >= 0.2 & out2$trimmed_ps <= 0.8))
})

test_that("adaptive method: ignores lower/upper, warns appropriately", {
  set.seed(2)

  n <- 80
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(-0.5 * x))
  ps <- predict(glm(z ~ x, family = binomial), type = "response")

  # 1) No user cutoffs
  out_adapt <- ps_trim(ps, method = "adaptive")
  # The result should have a computed 'cutoff'
  expect_true("cutoff" %in% names(out_adapt$method_info))

  # 2) If user sets lower/upper, we expect a warning
  expect_warning(
     ps_trim(ps, method = "adaptive", lower = 0.2, upper = 0.8),
     "For method='adaptive', `lower`/`upper` are ignored."
  )
})

test_that("pctl method: percentile-based trimming", {
  set.seed(3)

  n <- 100
  x <- rnorm(n)
  ps <- plogis(0.8 * x)

  # 1) Default [0.05, 0.95]
  out1 <- ps_trim(ps, method = "pctl")
  expect_equal(out1$method_info$lower, 0.05)
  expect_equal(out1$method_info$upper, 0.95)

  # The actual numeric cutoffs
  q_l <- quantile(ps, probs = 0.05)
  q_u <- quantile(ps, probs = 0.95)

  expect_true(all(out1$trimmed_ps >= q_l & out1$trimmed_ps <= q_u))

  # 2) Custom [0.2, 0.8]
  out2 <- ps_trim(ps, method = "pctl", lower = 0.2, upper = 0.8)
  q_l2 <- quantile(ps, probs = 0.2)
  q_u2 <- quantile(ps, probs = 0.8)

  expect_true(all(out2$trimmed_ps >= q_l2 & out2$trimmed_ps <= q_u2))
})

test_that("pref method: requires exposure, fails with all 0 or all 1", {
  set.seed(4)

  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(x))
  fit <- glm(z ~ x, family = binomial)
  ps  <- predict(fit, type = "response")

  # 1) If exposure = NULL, should fail
  expect_error(
    ps_trim(ps, method = "pref"),
    "must supply a binary 'exposure'"
  )

  # 2) If exposure is all 0 or all 1, fail
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

  # Check default [0.3, 0.7]
  expect_equal(out_pref$method_info$lower, 0.3)
  expect_equal(out_pref$method_info$upper, 0.7)

  # Check that the returned trimmed_ps is still on [0,1] scale (these are original PS)
  expect_true(all(out_pref$trimmed_ps >= 0 & out_pref$trimmed_ps <= 1))
  expect_true("pref_formula" %in% names(out_pref$method_info))
})

test_that("cr method: uses min(ps_treat) / max(ps_untrt), warns if cutoffs given", {
  set.seed(5)

  n <- 50
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 * x))
  fit <- glm(z ~ x, family = binomial)
  ps  <- predict(fit, type = "response")

  # 1) Must have exposure
  expect_error(
    ps_trim(ps, method = "cr"),
    "must supply a binary 'exposure'"
  )

  # 2) If all 0 or all 1 => fail
  expect_error(
    ps_trim(ps, exposure = rep(1, n), method = "cr"),
    "cannot compute common range"
  )

  # 3) Valid usage
  out_cr <- ps_trim(ps, exposure = z, method = "cr")

  # Check that user-specified lower/upper is warned about
  expect_warning(
    out_cr_warn <- ps_trim(ps, exposure = z, method = "cr", lower = 0.2, upper = 0.8),
    "ignored"
  )

  # The derived cr_lower, cr_upper
  ps_treat <- ps[z == 1]
  ps_untrt <- ps[z == 0]
  cr_l_expected <- min(ps_treat)
  cr_u_expected <- max(ps_untrt)

  expect_equal(out_cr$method_info$cr_lower, cr_l_expected)
  expect_equal(out_cr$method_info$cr_upper, cr_u_expected)

  # Check keep_idx is consistent
  keep_ps <- ps[out_cr$keep_idx]
  expect_true(all(keep_ps >= cr_l_expected & keep_ps <= cr_u_expected))
})

test_that("Edge cases: ps near 0 or 1, empty trimming result", {
  ps_edge <- c(0.0001, 0.01, 0.5, 0.99, 0.9999)
  z_edge  <- c(0, 0, 1, 1, 1)

  # Check all in-bounds or out-of-bounds
  expect_error(
    ps_trim(ps_edge, method = "ps", lower = 0.01, upper = 0.001),
    "`lower` >= `upper`.`lower` must be smaller than `upper`, and they must not be the same."
  )

  expect_error(
    ps_trim(ps_edge, method = "ps", lower = 0.01, upper = 0.01),
    "`lower` >= `upper`.`lower` must be smaller than `upper`, and they must not be the same."
  )

  # If we do normal cut
  out <- ps_trim(ps_edge, method = "ps", lower = 0.01, upper = 0.99)
  # The 0.0001 and 0.9999 are out, so we keep middle 3
  expect_length(out$trimmed_ps, 3)
  expect_true(all(out$trimmed_ps >= 0.01 & out$trimmed_ps <= 0.99))
})

test_that("Return structure is always correct, even if everything is trimmed", {
  # Suppose we intentionally trim everything
  ps_all <- seq(0.01, 0.99, length.out = 10)
  out_empty <- ps_trim(ps_all, method = "ps", lower = 1.1, upper = 1.2)

  expect_type(out_empty, "list")
  expect_named(out_empty, c("trimmed_ps", "method_info", "keep_idx", "trimmed_idx"))
  expect_length(out_empty$trimmed_ps, 0)
  expect_length(out_empty$keep_idx, 0)
  expect_length(out_empty$trimmed_idx, length(ps_all))
  expect_type(out_empty$method_info, "list")
})
