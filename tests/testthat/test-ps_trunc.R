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

  out_cr <- ps_trunc(ps, exposure = z, method = "cr")
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

test_that("ps_trunc() errors on invalid usage or exposure", {
  # if method="cr" but no exposure => error
  expect_error(ps_trunc(runif(10), method = "cr"), "0/1")

  # if exposure not 0/1 => error
  expect_error(ps_trunc(runif(5), exposure = 1:5, method = "cr"), "0/1")

  # if lower >= upper => error for method="ps"
  expect_error(ps_trunc(runif(5), method = "ps", lower = 0.8, upper = 0.3), "need lower < upper")
})
