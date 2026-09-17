# The Gruber et al. (2022) bound, c = sqrt(n) log(n) / 5, written out here so
# that every oracle below is computed independently of the package.
gruber_c <- function(n) {
  sqrt(n) * log(n) / 5
}

# Twenty scores, ten treated and ten untreated. No treated score lies above
# 1 - 1/c and no untreated score below 1/c, so the only bound either arm meets
# is the one that caps its own weights. Units 1 to 3 fall below 1/c = 0.373 and
# units 11 to 13 above 1 - 1/c = 0.627.
# fmt: skip
adaptive_ps <- c(
  0.05, 0.2, 0.3, 0.5, 0.6, 0.62, 0.45, 0.55, 0.4, 0.5,
  0.95, 0.8, 0.7, 0.5, 0.4, 0.38, 0.45, 0.55, 0.6, 0.39
)
adaptive_z <- rep(c(1, 0), each = 10)

# The same scores with a treated unit at 0.9 (unit 10) and an untreated one at
# 0.1 (unit 20). Each meets the bound that caps the other arm's weights.
nudged_ps <- replace(adaptive_ps, c(10, 20), c(0.9, 0.1))

# Winsorize by hand at [1/c, 1 - 1/c].
bound_by_hand <- function(ps, c) {
  pmin(pmax(ps, 1 / c), 1 - 1 / c)
}

# Unstabilized binary ATE weights, by hand.
ate_by_hand <- function(ps, z) {
  ifelse(z == 1, 1 / ps, 1 / (1 - ps))
}

adaptive_fit_data <- local({
  set.seed(4172)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  data.frame(
    x1 = x1,
    x2 = x2,
    z = rbinom(n, 1, plogis(2.2 * x1 - 1.4 * x2)),
    trt = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    a2 = factor(ifelse(
      runif(n) < plogis(2 * x1 - x2),
      "control",
      "treated"
    ))
  )
})

adaptive_binary_fit <- function() {
  glm(z ~ x1 + x2, data = adaptive_fit_data, family = binomial())
}

# ---- the bound ---------------------------------------------------------------

test_that("method = 'adaptive' bounds binary scores at [1/c, 1 - 1/c]", {
  c20 <- gruber_c(20)
  expect_equal(1 / c20, 0.3732089, tolerance = 1e-6)

  out <- expect_silent(ps_trunc(adaptive_ps, method = "adaptive"))

  expect_s3_class(out, "ps_trunc")
  expect_equal(
    as.numeric(out),
    bound_by_hand(adaptive_ps, c20),
    tolerance = 1e-12
  )
  expect_identical(which(is_unit_truncated(out)), c(1L, 2L, 3L, 11L, 12L, 13L))
  # The scores inside the bounds are handed back exactly as they arrived.
  untouched <- setdiff(seq_along(adaptive_ps), c(1:3, 11:13))
  expect_identical(as.numeric(out)[untouched], adaptive_ps[untouched])
})

test_that("the adaptive record names the method, both bounds, and the units", {
  c20 <- gruber_c(20)
  out <- ps_trunc(adaptive_ps, method = "adaptive")
  meta <- ps_trunc_meta(out)

  expect_identical(
    names(meta),
    c("method", "lower_bound", "upper_bound", "truncated_idx", "n_obs")
  )
  expect_identical(meta$method, "adaptive")
  expect_equal(meta$lower_bound, 1 / c20, tolerance = 1e-12)
  expect_equal(meta$upper_bound, 1 - 1 / c20, tolerance = 1e-12)
  expect_identical(meta$truncated_idx, c(1L, 2L, 3L, 11L, 12L, 13L))
  expect_identical(meta$n_obs, 20L)
  expect_true(is_ps_truncated(out))
})

test_that("the adaptive bound on a vector prints its method and bounds", {
  expect_snapshot(print(ps_trunc(adaptive_ps, method = "adaptive")))
})

test_that("the adaptive bound loosens as the number of scores grows", {
  # At n = 200 the bound is [0.0667, 0.9333], far wider than at n = 20.
  set.seed(52)
  ps <- plogis(rnorm(200, sd = 2.5))
  c200 <- gruber_c(200)

  out <- ps_trunc(ps, method = "adaptive")
  meta <- ps_trunc_meta(out)

  expect_equal(meta$lower_bound, 1 / c200, tolerance = 1e-12)
  expect_equal(meta$upper_bound, 1 - 1 / c200, tolerance = 1e-12)
  expect_equal(as.numeric(out), bound_by_hand(ps, c200), tolerance = 1e-12)
  expect_identical(
    meta$truncated_idx,
    which(ps < 1 / c200 | ps > 1 - 1 / c200)
  )
  expect_gt(length(meta$truncated_idx), 0)
})

test_that("the adaptive bound counts only the scores that are present", {
  # Twenty of 22 scores are present, so the floor is 1/c(20) = 0.3732. Counting
  # the two missing scores would put it at 1/c(22) = 0.3449 and leave the score
  # of 0.36 at unit 9 alone.
  ps <- c(replace(adaptive_ps, 9, 0.36), NA, NA)
  c20 <- gruber_c(20)
  expect_lt(1 / gruber_c(22), 0.36)
  expect_gt(1 / c20, 0.36)

  out <- ps_trunc(ps, method = "adaptive")
  meta <- ps_trunc_meta(out)

  expect_equal(meta$lower_bound, 1 / c20, tolerance = 1e-12)
  expect_equal(meta$upper_bound, 1 - 1 / c20, tolerance = 1e-12)
  expect_identical(meta$truncated_idx, c(1L, 2L, 3L, 9L, 11L, 12L, 13L))
  expect_identical(meta$n_obs, 22L)
  expect_equal(as.numeric(out), bound_by_hand(ps, c20), tolerance = 1e-12)
  # A missing score stays missing and is not a unit the bound moved.
  expect_true(all(is.na(out[21:22])))
  expect_false(any(is_unit_truncated(out)[21:22]))
})

# ---- the equivalence with a weight bound -------------------------------------

test_that("adaptive score bounds cap unstabilized ATE weights at c", {
  c20 <- gruber_c(20)
  w_raw <- ate_by_hand(adaptive_ps, adaptive_z)
  expect_gt(max(w_raw), c20)

  bounded <- ps_trunc(adaptive_ps, method = "adaptive")
  wts <- wt_ate(bounded, .exposure = adaptive_z, .focal_level = 1)

  expect_equal(as.numeric(wts), pmin(w_raw, c20), tolerance = 1e-12)
  expect_equal(max(as.numeric(wts)), c20, tolerance = 1e-12)
})

test_that("the adaptive score bound agrees with the adaptive weight bound", {
  # With no unit meeting the bound on the other arm's side, the two routes give
  # the same weights: the score bound and the weight bound share c.
  bounded <- ps_trunc(adaptive_ps, method = "adaptive")
  on_scores <- wt_ate(bounded, .exposure = adaptive_z, .focal_level = 1)
  on_weights <- wt_trunc(
    wt_ate(adaptive_ps, .exposure = adaptive_z, .focal_level = 1),
    method = "adaptive"
  )

  expect_equal(
    as.numeric(on_scores),
    as.numeric(on_weights),
    tolerance = 1e-12
  )
  expect_equal(
    attr(on_weights, "psw_trunc_meta")$upper_value,
    1 / ps_trunc_meta(bounded)$lower_bound,
    tolerance = 1e-12
  )
})

test_that("the score bound also nudges units near the smallest weight", {
  # Unit 10 is treated at 0.9 and unit 20 untreated at 0.1, so each has a
  # weight of 1.11. The score bound moves both to the bound on the other side,
  # giving each a weight of c / (c - 1) = 1.595, where the weight bound leaves
  # them alone. Every other unit gets the same weight either way.
  c20 <- gruber_c(20)
  nudge <- c20 / (c20 - 1)
  w_raw <- ate_by_hand(nudged_ps, adaptive_z)

  bounded <- ps_trunc(nudged_ps, method = "adaptive")
  on_scores <- wt_ate(bounded, .exposure = adaptive_z, .focal_level = 1)
  on_weights <- wt_trunc(
    wt_ate(nudged_ps, .exposure = adaptive_z, .focal_level = 1),
    method = "adaptive"
  )

  expect_identical(
    ps_trunc_meta(bounded)$truncated_idx,
    c(1L, 2L, 3L, 10L, 11L, 12L, 13L, 20L)
  )

  others <- setdiff(seq_along(nudged_ps), c(10, 20))
  expect_equal(
    as.numeric(on_scores)[others],
    pmin(w_raw, c20)[others],
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(on_scores)[others],
    as.numeric(on_weights)[others],
    tolerance = 1e-12
  )

  expect_equal(
    as.numeric(on_scores)[c(10, 20)],
    c(nudge, nudge),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(on_weights)[c(10, 20)],
    w_raw[c(10, 20)],
    tolerance = 1e-12
  )
  # The nudge raises a weight that was already below c / (c - 1) up to it, and
  # every weight above c / (c - 1) is the same under either bound.
  expect_true(all(w_raw[c(10, 20)] < nudge))
  large <- as.numeric(on_scores) > nudge + 1e-12
  expect_equal(
    as.numeric(on_scores)[large],
    as.numeric(on_weights)[large],
    tolerance = 1e-12
  )
})

test_that("the adaptive bound on a fitted model's scores caps ATE weights at c", {
  fit <- adaptive_binary_fit()
  ps <- unname(predict(fit, type = "response"))
  z <- adaptive_fit_data$z
  c200 <- gruber_c(200)
  w_raw <- ate_by_hand(ps, z)
  expect_gt(max(w_raw), c200)

  bounded <- ps_trunc(ps, method = "adaptive")
  wts <- as.numeric(wt_ate(bounded, .exposure = z, .focal_level = 1))

  # The units a bound on the other arm's side moves are the treated units above
  # 1 - 1/c and the untreated units below 1/c. Everywhere else the weights are
  # the capped weights, and there they are c / (c - 1).
  nudged <- (z == 1 & ps > 1 - 1 / c200) | (z == 0 & ps < 1 / c200)
  expect_true(any(nudged))
  expect_equal(wts[!nudged], pmin(w_raw, c200)[!nudged], tolerance = 1e-12)
  expect_equal(
    wts[nudged],
    rep(c200 / (c200 - 1), sum(nudged)),
    tolerance = 1e-12
  )
  expect_equal(
    wts,
    ate_by_hand(bound_by_hand(ps, c200), z),
    tolerance = 1e-12
  )
})

test_that("stabilized weights from the adaptive bound are capped per arm", {
  # Stabilizing multiplies each arm's weights by its prevalence, so the treated
  # weights are capped at p1 * c and the untreated at p0 * c. Unit 15 is moved
  # to the treated arm, which puts p1 at 0.55 and p0 at 0.45.
  z <- replace(adaptive_z, 15, 1)
  c20 <- gruber_c(20)
  p1 <- mean(z)
  expect_equal(p1, 0.55)

  bounded <- ps_trunc(adaptive_ps, method = "adaptive")
  wts <- as.numeric(wt_ate(bounded, .exposure = z, stabilize = TRUE))

  expect_equal(
    wts,
    ifelse(
      z == 1,
      p1 / bound_by_hand(adaptive_ps, c20),
      (1 - p1) / (1 - bound_by_hand(adaptive_ps, c20))
    ),
    tolerance = 1e-12
  )
  expect_equal(max(wts[z == 1]), p1 * c20, tolerance = 1e-12)
  expect_equal(max(wts[z == 0]), (1 - p1) * c20, tolerance = 1e-12)
})

# ---- weights built from the result -------------------------------------------

test_that("weights from adaptive scores are labelled as built from truncated scores", {
  bounded <- ps_trunc(adaptive_ps, method = "adaptive")

  ate <- wt_ate(bounded, .exposure = adaptive_z, .focal_level = 1)
  expect_identical(estimand(ate), "ate; truncated")
  expect_true(is_ps_truncated(ate))
  expect_false(is_wt_truncated(ate))
  expect_identical(ps_trunc_meta(ate), ps_trunc_meta(bounded))
  expect_identical(is_unit_truncated(ate), is_unit_truncated(bounded))

  att <- wt_att(bounded, .exposure = adaptive_z, .focal_level = 1)
  expect_identical(estimand(att), "att; truncated")
  expect_true(is_ps_truncated(att))
  expect_identical(ps_trunc_meta(att)$method, "adaptive")

  # The weight bound labels the weights, not the scores.
  on_weights <- wt_trunc(
    wt_ate(adaptive_ps, .exposure = adaptive_z, .focal_level = 1),
    method = "adaptive"
  )
  expect_identical(estimand(on_weights), "ate; weights truncated")
  expect_false(is_ps_truncated(on_weights))
  expect_true(is_wt_truncated(on_weights))
})

# ---- the record is told apart from "ps" --------------------------------------

test_that("an adaptive bound is not the same truncation as an equal 'ps' bound", {
  c20 <- gruber_c(20)
  adaptive <- ps_trunc(adaptive_ps, method = "adaptive")
  fixed <- ps_trunc(
    adaptive_ps,
    method = "ps",
    lower = 1 / c20,
    upper = 1 - 1 / c20
  )

  expect_equal(as.numeric(adaptive), as.numeric(fixed), tolerance = 1e-12)
  expect_equal(
    ps_trunc_meta(adaptive)$lower_bound,
    ps_trunc_meta(fixed)$lower_bound
  )

  params_adaptive <- trunc_parameters(ps_trunc_meta(adaptive))
  params_fixed <- trunc_parameters(ps_trunc_meta(fixed))
  expect_identical(params_adaptive$method, "adaptive")
  expect_identical(params_fixed$method, "ps")
  expect_false(identical(params_adaptive, params_fixed))

  expect_warning(
    combined <- c(adaptive, fixed),
    class = "propensity_coercion_warning"
  )
  expect_false(inherits(combined, "ps_trunc"))
})

test_that("two adaptive bounds over the same number of scores combine", {
  first <- ps_trunc(adaptive_ps, method = "adaptive")
  second <- ps_trunc(rev(adaptive_ps), method = "adaptive")

  expect_identical(
    trunc_parameters(ps_trunc_meta(first)),
    trunc_parameters(ps_trunc_meta(second))
  )
  combined <- expect_no_warning(c(first, second))
  expect_s3_class(combined, "ps_trunc")
  expect_identical(ps_trunc_meta(combined)$method, "adaptive")
})

# ---- supplied bounds ---------------------------------------------------------

test_that("bounds supplied to the adaptive method are ignored with a warning", {
  plain <- ps_trunc(adaptive_ps, method = "adaptive")

  expect_warning(
    with_lower <- ps_trunc(adaptive_ps, method = "adaptive", lower = 0.2),
    class = "propensity_warning",
    regexp = "ignored"
  )
  expect_warning(
    with_upper <- ps_trunc(adaptive_ps, method = "adaptive", upper = 0.8),
    class = "propensity_warning",
    regexp = "ignored"
  )
  expect_warning(
    with_both <- ps_trunc(
      adaptive_ps,
      method = "adaptive",
      lower = 0.2,
      upper = 0.8
    ),
    class = "propensity_warning",
    regexp = "ignored"
  )

  # Nothing the caller wrote reaches the bounds or the record.
  expect_identical(with_lower, plain)
  expect_identical(with_upper, plain)
  expect_identical(with_both, plain)

  expect_propensity_warning(
    ps_trunc(adaptive_ps, method = "adaptive", lower = 0.2, upper = 0.8)
  )
})

# ---- too few scores ----------------------------------------------------------

test_that("the adaptive bound needs at least 15 scores", {
  # Below c = 2 the floor 1/c lies above the ceiling 1 - 1/c, so there is no
  # interval to bound the scores to. c first exceeds 2 at n = 15.
  expect_lt(gruber_c(14), 2)
  expect_gt(gruber_c(15), 2)
  expect_lt(gruber_c(6), 1)

  ps15 <- seq(0.2, 0.8, length.out = 15)
  c15 <- gruber_c(15)
  out <- ps_trunc(ps15, method = "adaptive")
  expect_equal(as.numeric(out), bound_by_hand(ps15, c15), tolerance = 1e-12)
  expect_lt(ps_trunc_meta(out)$lower_bound, ps_trunc_meta(out)$upper_bound)

  expect_error(
    ps_trunc(ps15[-1], method = "adaptive"),
    class = "propensity_range_error"
  )
  expect_error(
    ps_trunc(ps15[1:6], method = "adaptive"),
    class = "propensity_range_error"
  )
  expect_error(
    ps_trunc(0.5, method = "adaptive"),
    class = "propensity_range_error"
  )
  # A missing score does not count toward the fifteen.
  expect_error(
    ps_trunc(replace(ps15, 3, NA), method = "adaptive"),
    class = "propensity_range_error"
  )

  expect_propensity_error(ps_trunc(ps15[1:14], method = "adaptive"))
  expect_propensity_error(ps_trunc(ps15[1:6], method = "adaptive"))
})

# ---- routes ------------------------------------------------------------------

test_that("the adaptive bound reads a binary data frame's score column", {
  df <- data.frame(control = 1 - adaptive_ps, treated = adaptive_ps)

  from_df <- ps_trunc(df, method = "adaptive")
  oracle <- ps_trunc(adaptive_ps, method = "adaptive")

  expect_equal(as.numeric(from_df), as.numeric(oracle), tolerance = 1e-12)
  expect_identical(class(from_df), class(oracle))
  expect_identical(ps_trunc_meta(from_df), ps_trunc_meta(oracle))

  from_df_z <- ps_trunc(df, method = "adaptive", .exposure = adaptive_z)
  expect_identical(ps_trunc_meta(from_df_z), ps_trunc_meta(oracle))
})

test_that("the adaptive bound on a binomial fit bounds the scores it reports", {
  fit <- adaptive_binary_fit()
  scores <- predict(fit, type = "response")

  from_fit <- expect_silent(ps_trunc(fit, method = "adaptive"))
  oracle <- ps_trunc(scores, method = "adaptive")

  expect_equal(as.numeric(from_fit), as.numeric(oracle), tolerance = 1e-12)
  expect_identical(class(from_fit), class(oracle))
  expect_identical(ps_trunc_meta(from_fit), ps_trunc_meta(oracle))
  expect_identical(ps_trunc_meta(from_fit)$method, "adaptive")
  expect_equal(
    ps_trunc_meta(from_fit)$lower_bound,
    1 / gruber_c(200),
    tolerance = 1e-12
  )
})

test_that("the adaptive bound reads a two-level multinomial fit as binary", {
  skip_if_not_installed("nnet")

  fit <- nnet::multinom(a2 ~ x1 + x2, data = adaptive_fit_data, trace = FALSE)
  scores <- as.numeric(fitted(fit))

  from_fit <- ps_trunc(fit, method = "adaptive")
  oracle <- ps_trunc(scores, method = "adaptive")

  expect_equal(as.numeric(from_fit), as.numeric(oracle), tolerance = 1e-12)
  expect_identical(ps_trunc_meta(from_fit), ps_trunc_meta(oracle))
})

# ---- categorical refusals ----------------------------------------------------

test_that("the adaptive bound refuses a categorical score matrix", {
  ps_mat <- matrix(
    c(0.02, 0.49, 0.49, 0.3, 0.4, 0.3, 0.6, 0.2, 0.2),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )
  trt <- factor(c("a", "b", "c"))

  expect_error(
    ps_trunc(ps_mat, method = "adaptive", .exposure = trt),
    class = "propensity_method_error"
  )
  expect_error(
    ps_trunc(as.data.frame(ps_mat), method = "adaptive", .exposure = trt),
    class = "propensity_method_error"
  )
  expect_propensity_error(
    ps_trunc(ps_mat, method = "adaptive", .exposure = trt)
  )
})

test_that("the adaptive bound refuses a multinomial fit of three levels", {
  skip_if_not_installed("nnet")

  fit <- nnet::multinom(trt ~ x1 + x2, data = adaptive_fit_data, trace = FALSE)

  expect_error(
    ps_trunc(fit, method = "adaptive"),
    class = "propensity_method_error"
  )
  expect_propensity_error(ps_trunc(fit, method = "adaptive"))
})
