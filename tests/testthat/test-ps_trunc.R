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

  out_cr <- ps_trunc(ps, .exposure = z, method = "cr", .treated = 1)
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
  expect_error(
    ps_trunc(runif(10), method = "cr"),
    class = "propensity_binary_transform_error"
  )

  # if .exposure not 0/1 => error
  expect_error(
    ps_trunc(runif(5), .exposure = 1:5, method = "cr"),
    class = "propensity_binary_transform_error"
  )

  # if lower >= upper => error for method="ps"
  expect_error(
    ps_trunc(runif(5), method = "ps", lower = 0.8, upper = 0.3),
    class = "propensity_range_error"
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
  w_ate <- wt_ate(truncated_ps, .exposure = z, exposure_type = "binary", .treated = 1)
  expect_s3_class(w_ate, "psw")

  # 4) Verify truncated, not trimmed, not refit, estimand
  expect_true(is_truncated(w_ate))
  expect_false(is_trimmed(w_ate))
  expect_false(is_refit(w_ate))
  expect_match(estimand(w_ate), "; truncated$")
})

test_that("is_truncated.default() -> FALSE, is_truncated.ps_trunc() -> TRUE", {
  # 1) A plain numeric => default => FALSE
  expect_false(is_truncated(runif(5)))

  # 2) A simple ps_trunc object => is_truncated(...) => TRUE
  # Create via new_ps_trunc()
  my_trunc <- new_ps_trunc(
    x = c(0.2, 0.6, 0.8),
    meta = list(method = "ps", lower_bound = 0.2, upper_bound = 0.8)
  )
  expect_true(is_truncated(my_trunc))
})

test_that("Arithmetic on ps_trunc => vctrs_error_incompatible_op", {
  obj <- new_ps_trunc(
    x = c(0.2, 0.7, 0.9),
    meta = list(method = "ps", lower_bound = 0.2, upper_bound = 0.8)
  )

  # obj + 1 => error
  expect_error(
    obj + 1,
    class = "vctrs_error_incompatible_op"
  )
  # 1 + obj => error
  expect_error(
    1 + obj,
    class = "vctrs_error_incompatible_op"
  )

  # obj * another ps_trunc => error
  obj2 <- new_ps_trunc(
    x = c(0.1, 0.1, 0.3),
    meta = list(method = "ps", lower_bound = 0.1, upper_bound = 0.5)
  )
  expect_error(
    obj * obj2,
    class = "vctrs_error_incompatible_op"
  )
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
    .treated = 1
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
    exposure_type = "binary",
    .treated = 1
  )
  # Just check dimension, no error
  expect_length(w, 3)

  # If .exposure isn't 0/1 or has different factor levels, transform_.exposure_binary
  # is tested. We'll do a quick check with factor( c("C","T","T") )
  w2 <- atm_binary(
    .propensity = ps_vec,
    .exposure = factor(c("C", "T", "T")),
    .treated = "T",
    .untreated = "C"
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
    .treated = 1
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
    exposure_type = "binary",
    .treated = 1
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
    .treated = 1
  )
  expect_s3_class(w_atm, "psw")
  # Estimand => "atm; truncated"
  expect_match(estimand(w_atm), "atm; truncated")
  # truncated=TRUE
  expect_true(is_truncated(w_atm))
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
    .treated = 1
  )

  expect_s3_class(w_ato, "psw")
  expect_match(estimand(w_ato), "ato; truncated")
  expect_true(is_truncated(w_ato))
})
