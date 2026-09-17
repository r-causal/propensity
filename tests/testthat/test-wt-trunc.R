# Truncating the weights themselves ------------------------------------------

# `wt_trunc()` winsorizes a psw at a bound read on the scale of the weights and
# leaves a record of the bound and the units it moved. Every bound below is
# worked out from the data by hand, never read back from the function.

# The record the truncation should leave, written out.
expected_wt_trunc_record <- function(
  method,
  lower = NULL,
  upper = NULL,
  lower_value = NULL,
  upper_value = NULL,
  truncated_idx = integer(),
  n_obs
) {
  new_psw_trunc_meta(
    method = method,
    lower = lower,
    upper = upper,
    lower_value = lower_value,
    upper_value = upper_value,
    truncated_idx = truncated_idx,
    n_obs = n_obs
  )
}

# Winsorizing by hand. A missing weight is neither above nor below a bound, so
# it stays missing and is not among the units moved.
winsorize_by_hand <- function(x, lower_value = -Inf, upper_value = Inf) {
  x <- as.numeric(x)
  low <- which(x < lower_value)
  high <- which(x > upper_value)
  x[low] <- lower_value
  x[high] <- upper_value
  list(values = x, idx = sort(c(low, high)))
}

# The Gruber et al. bound, with `n` the number of weights present.
gruber_bound <- function(x) {
  n <- sum(!is.na(x))
  sqrt(n) * log(n) / 5
}

# Quantile type 7, the `quantile()` default that `ps_trunc(method = "pctl")`
# uses, written out: the order statistic at `h = (n - 1) p + 1`, interpolated
# linearly when `h` falls between two of them.
type7_quantile <- function(x, p) {
  x <- sort(x[!is.na(x)])
  h <- (length(x) - 1) * p + 1
  lo <- floor(h)
  hi <- ceiling(h)
  x[lo] + (h - lo) * (x[hi] - x[lo])
}

# The message a refusal carries, with the line wrapping cli adds squeezed out
# so that a phrase can be searched for whatever the console width.
refusal_message <- function(expr) {
  cnd <- rlang::catch_cnd(expr, classes = "error")
  if (is.null(cnd)) {
    return(NA_character_)
  }
  gsub("\\s+", " ", conditionMessage(cnd))
}

# Weights from real fits for each exposure type the verb accepts. The
# continuous weights are stabilized by a fitted numerator model, so they carry
# a density record and the model it names.
wt_trunc_data <- function(n = 60, seed = 20260916) {
  withr::local_seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  dat$z <- rbinom(n, 1, stats::plogis(1.2 * dat$x1))
  dat$k <- factor(
    sample(c("low", "mid", "high"), n, replace = TRUE),
    levels = c("low", "mid", "high")
  )
  dat$a <- 1 + 0.6 * dat$x1 - 0.4 * dat$x2 + rnorm(n)
  dat
}

binary_wt_trunc_fit <- function(dat = wt_trunc_data()) {
  stats::glm(z ~ x1 + x2, data = dat, family = stats::binomial())
}

binary_wt_trunc_weights <- function(dat = wt_trunc_data()) {
  wt_ate(binary_wt_trunc_fit(dat), .exposure = dat$z)
}

stabilized_binary_wt_trunc_weights <- function(dat = wt_trunc_data()) {
  numerator <- stats::glm(z ~ 1, data = dat, family = stats::binomial())
  wt_ate(binary_wt_trunc_fit(dat), .exposure = dat$z, stabilize = numerator)
}

categorical_wt_trunc_weights <- function(dat = wt_trunc_data()) {
  fit <- nnet::multinom(k ~ x1 + x2, data = dat, trace = FALSE)
  wt_ate(fit, .exposure = dat$k)
}

continuous_wt_trunc_weights <- function(dat = wt_trunc_data()) {
  dose <- stats::lm(a ~ x1 + x2, data = dat)
  numerator <- stats::lm(a ~ x1, data = dat)
  muffle_variance_warning(wt_ate(dose, stabilize = numerator))
}

# The attributes a truncation has no business changing: everything but the
# label, the flag, and the record it writes.
untouched_attrs <- function(x) {
  attrs <- attributes(x)
  attrs[setdiff(names(attrs), c("estimand", "wt_truncated", "psw_trunc_meta"))]
}

# ---- the adaptive bound ------------------------------------------------------

test_that("the adaptive bound is sqrt(n) log(n) / 5 on the upper tail only", {
  # Twenty weights put the bound at 2.679, so the two largest are pulled down to
  # it and the smallest, far below any sensible lower bound, is left alone.
  x <- c(rep(1.5, 16), 0.01, 2, 5, 9)
  w <- psw(x, estimand = "ate")
  bound <- sqrt(20) * log(20) / 5
  by_hand <- winsorize_by_hand(x, upper_value = bound)

  out <- expect_silent(wt_trunc(w))

  expect_s3_class(out, "psw")
  expect_equal(as.numeric(out), by_hand$values)
  expect_identical(by_hand$idx, c(19L, 20L))
  expect_true(is_wt_truncated(out))
  meta <- attr(out, "psw_trunc_meta")
  expect_s3_class(meta, "propensity_psw_trunc_meta")
  expect_equal(meta$upper_value, bound)
  expect_identical(
    unclass(meta)[c("method", "lower", "upper", "lower_value")],
    list(method = "adaptive", lower = NULL, upper = NULL, lower_value = NULL)
  )
  expect_identical(meta$truncated_idx, by_hand$idx)
  expect_identical(meta$n_obs, 20L)
  expect_identical(
    is_unit_wt_truncated(out),
    seq_along(x) %in% by_hand$idx
  )
})

test_that("adaptive is the default method", {
  x <- c(rep(1.5, 16), 0.01, 2, 5, 9)
  w <- psw(x, estimand = "ate")

  expect_identical(wt_trunc(w), wt_trunc(w, method = "adaptive"))
})

test_that("the adaptive bound counts only the weights that are present", {
  # With 20 of 22 weights present the bound is 2.679; counting the two missing
  # weights would put it at 2.720 and leave the weight of 2.7 alone.
  x <- c(rep(1.5, 16), NA, 2.7, 1, NA, 9, 2)
  w <- psw(x, estimand = "ate")
  bound <- sqrt(20) * log(20) / 5
  expect_lt(bound, 2.7)
  expect_gt(sqrt(22) * log(22) / 5, 2.7)

  out <- wt_trunc(w)
  meta <- attr(out, "psw_trunc_meta")

  expect_equal(meta$upper_value, bound)
  expect_identical(meta$truncated_idx, c(18L, 21L))
  expect_identical(meta$n_obs, 22L)
  expect_equal(
    as.numeric(out),
    winsorize_by_hand(x, upper_value = bound)$values
  )
  # A missing weight stays missing and is not a unit the bound moved.
  expect_true(all(is.na(out[c(17, 20)])))
  expect_false(any(is_unit_wt_truncated(out)[c(17, 20)]))
})

test_that("the adaptive bound on weights from a trimmed score counts the kept units", {
  # The trim keeps 77 of 80 units, which puts the bound at 7.62 rather than the
  # 7.84 all 80 would give. Two weights exceed it, and one of those lies
  # between the two bounds, so counting the trimmed units would move only one.
  dat <- wt_trunc_data(n = 80)
  fit <- binary_wt_trunc_fit(dat)
  trimmed <- ps_trim(fit, method = "ps", lower = 0.05, upper = 0.95)
  refit <- ps_refit(trimmed, fit)
  w <- wt_ate(refit, .exposure = dat$z)
  expect_identical(sum(!is.na(w)), 77L)

  bound <- gruber_bound(as.numeric(w))
  expect_equal(bound, sqrt(77) * log(77) / 5)
  by_hand <- winsorize_by_hand(w, upper_value = bound)
  expect_length(by_hand$idx, 2)
  all_units_bound <- sqrt(80) * log(80) / 5
  expect_identical(
    sum(as.numeric(w) > bound & as.numeric(w) <= all_units_bound, na.rm = TRUE),
    1L
  )

  out <- wt_trunc(w)

  expect_equal(attr(out, "psw_trunc_meta")$upper_value, bound)
  expect_identical(attr(out, "psw_trunc_meta")$truncated_idx, by_hand$idx)
  expect_equal(as.numeric(out), by_hand$values)
})

test_that("bounds supplied to the adaptive method are ignored with a warning", {
  x <- c(rep(1.5, 16), 0.01, 2, 5, 9)
  w <- psw(x, estimand = "ate")
  plain <- wt_trunc(w)

  expect_warning(
    with_upper <- wt_trunc(w, upper = 3),
    class = "propensity_warning",
    regexp = "ignored"
  )
  expect_warning(
    with_lower <- wt_trunc(w, lower = 0.5),
    class = "propensity_warning",
    regexp = "ignored"
  )
  expect_warning(
    with_both <- wt_trunc(w, method = "adaptive", lower = 0.5, upper = 3),
    class = "propensity_warning",
    regexp = "ignored"
  )

  # Nothing the caller wrote reaches the bound or the record.
  expect_identical(with_upper, plain)
  expect_identical(with_lower, plain)
  expect_identical(with_both, plain)

  expect_propensity_warning(wt_trunc(w, method = "adaptive", upper = 3))
})

# ---- an absolute bound -------------------------------------------------------

test_that("method = 'wt' winsorizes the upper tail at the value given", {
  x <- c(0.5, 1, 4, 2, 25, NA, 0.02)
  w <- psw(x, estimand = "ate")
  by_hand <- winsorize_by_hand(x, upper_value = 20)

  out <- wt_trunc(w, method = "wt", upper = 20)

  expect_equal(as.numeric(out), c(0.5, 1, 4, 2, 20, NA, 0.02))
  expect_equal(as.numeric(out), by_hand$values)
  expect_identical(
    attr(out, "psw_trunc_meta"),
    expected_wt_trunc_record(
      method = "wt",
      upper = 20,
      upper_value = 20,
      truncated_idx = 5L,
      n_obs = 7L
    )
  )
})

test_that("method = 'wt' bounds the lower tail only when asked to", {
  x <- c(0.5, 1, 4, 2, 25, NA, 0.02)
  w <- psw(x, estimand = "ate")

  out <- wt_trunc(w, method = "wt", lower = 0.05, upper = 20)

  expect_equal(as.numeric(out), c(0.5, 1, 4, 2, 20, NA, 0.05))
  expect_identical(
    attr(out, "psw_trunc_meta"),
    expected_wt_trunc_record(
      method = "wt",
      lower = 0.05,
      upper = 20,
      lower_value = 0.05,
      upper_value = 20,
      truncated_idx = c(5L, 7L),
      n_obs = 7L
    )
  )
  expect_identical(
    is_unit_wt_truncated(out),
    c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE)
  )
})

test_that("method = 'wt' accepts a floor of zero", {
  # Weights are never negative, so a floor of zero is a valid bound that
  # moves nothing.
  x <- c(0.5, 1, 4, 2, 25)
  w <- psw(x, estimand = "ate")

  out <- expect_silent(wt_trunc(w, method = "wt", lower = 0, upper = 20))

  expect_equal(as.numeric(out), c(0.5, 1, 4, 2, 20))
  meta <- attr(out, "psw_trunc_meta")
  expect_identical(meta$lower, 0)
  expect_identical(meta$lower_value, 0)
  expect_identical(meta$truncated_idx, 5L)
})

test_that("integer weights come back as doubles", {
  out <- wt_trunc(c(1L, 5L, 30L), method = "wt", upper = 2.5)

  expect_s3_class(out, "psw")
  expect_type(vctrs::vec_data(out), "double")
  expect_identical(as.numeric(out), c(1, 2.5, 2.5))
  expect_identical(attr(out, "psw_trunc_meta")$truncated_idx, c(2L, 3L))
})

test_that("a bound that moves nothing still leaves a record", {
  x <- c(0.5, 1, 4, 2)
  w <- psw(x, estimand = "ate")

  out <- wt_trunc(w, method = "wt", upper = 20)

  expect_equal(as.numeric(out), x)
  expect_true(is_wt_truncated(out))
  expect_identical(attr(out, "psw_trunc_meta")$truncated_idx, integer())
  expect_identical(is_unit_wt_truncated(out), rep(FALSE, 4))
})

test_that("an integer bound is recorded as the same double as its numeric twin", {
  w <- psw(c(1.2, 4, 1.5, 2.1, 3.5), estimand = "ate")

  from_integer <- wt_trunc(w, method = "wt", upper = 3L)
  from_double <- wt_trunc(w, method = "wt", upper = 3)
  meta <- attr(from_integer, "psw_trunc_meta")

  expect_identical(meta, attr(from_double, "psw_trunc_meta"))
  expect_identical(meta$upper, 3)
  expect_identical(meta$upper_value, 3)
  expect_identical(from_integer, from_double)

  # Two truncations at the same bound describe the same operation.
  combined <- expect_silent(c(from_integer, from_double))
  expect_true(is_wt_truncated(combined))

  lower_integer <- wt_trunc(w, method = "wt", lower = 2L, upper = 3L)
  lower_double <- wt_trunc(w, method = "wt", lower = 2, upper = 3)
  expect_identical(
    attr(lower_integer, "psw_trunc_meta"),
    attr(lower_double, "psw_trunc_meta")
  )
  expect_identical(attr(lower_integer, "psw_trunc_meta")$lower, 2)
  expect_silent(c(lower_integer, lower_double))
})

test_that("method = 'wt' refuses a bound it cannot apply", {
  w <- psw(c(1.2, 4, 1.5, 2.1, 3.5), estimand = "ate")

  expect_error(
    wt_trunc(w, method = "wt"),
    class = "propensity_missing_arg_error"
  )
  expect_error(
    wt_trunc(w, method = "wt", lower = 1),
    class = "propensity_missing_arg_error"
  )
  expect_error(
    wt_trunc(w, method = "wt", upper = NA),
    class = "propensity_missing_value_error"
  )
  for (upper in list(0, -1, Inf, c(2, 3))) {
    expect_error(
      wt_trunc(w, method = "wt", upper = upper),
      class = "propensity_range_error"
    )
  }
  expect_error(
    wt_trunc(w, method = "wt", upper = "3"),
    class = "propensity_type_error"
  )
  expect_error(
    wt_trunc(w, method = "wt", lower = 3, upper = 3),
    class = "propensity_range_error"
  )
  expect_error(
    wt_trunc(w, method = "wt", lower = 4, upper = 3),
    class = "propensity_range_error"
  )
  expect_error(
    wt_trunc(w, method = "wt", lower = -1, upper = 3),
    class = "propensity_range_error"
  )

  expect_propensity_error(wt_trunc(w, method = "wt"))
  expect_propensity_error(wt_trunc(w, method = "wt", upper = 0))
  expect_propensity_error(wt_trunc(w, method = "wt", lower = 4, upper = 3))
})

# ---- a percentile bound ------------------------------------------------------

# 101 distinct weights and two missing ones, in no particular order.
pctl_weights <- function() {
  withr::local_seed(99)
  x <- sample(round(0.5 + stats::rexp(101), 6))
  x <- c(x[1:40], NA, x[41:90], NA, x[91:101])
  psw(x, estimand = "ate")
}

test_that("method = 'pctl' reads the bound at a type 7 quantile of the weights present", {
  w <- pctl_weights()
  x <- as.numeric(w)
  present <- sort(x[!is.na(x)])
  expect_length(present, 101)
  expect_identical(anyDuplicated(present), 0L)

  # 0.99 of 101 falls on an order statistic, h = 100 * 0.99 + 1 = 100.
  exact <- wt_trunc(w, method = "pctl", upper = 0.99)
  expect_identical(attr(exact, "psw_trunc_meta")$upper_value, present[100])
  expect_identical(
    attr(exact, "psw_trunc_meta")$truncated_idx,
    which(x > present[100])
  )
  expect_equal(
    as.numeric(exact),
    winsorize_by_hand(x, upper_value = present[100])$values
  )

  # 0.975 of 101 falls between two, h = 98.5, so the bound is interpolated.
  between <- present[98] + 0.5 * (present[99] - present[98])
  expect_equal(type7_quantile(x, 0.975), between)
  interpolated <- wt_trunc(w, method = "pctl", upper = 0.975)
  meta <- attr(interpolated, "psw_trunc_meta")
  expect_equal(meta$upper_value, between)
  expect_identical(meta$truncated_idx, which(x > between))
  expect_length(meta$truncated_idx, 3)
  expect_identical(
    unclass(meta)[c("method", "lower", "upper", "lower_value", "n_obs")],
    list(
      method = "pctl",
      lower = NULL,
      upper = 0.975,
      lower_value = NULL,
      n_obs = 103L
    )
  )
  expect_null(names(meta$upper_value))
})

test_that("method = 'pctl' bounds the lower tail only when asked to", {
  w <- pctl_weights()
  x <- as.numeric(w)
  present <- sort(x[!is.na(x)])

  # h = 100 * 0.025 + 1 = 3.5 at the bottom, 98.5 at the top.
  lower_value <- present[3] + 0.5 * (present[4] - present[3])
  upper_value <- present[98] + 0.5 * (present[99] - present[98])
  by_hand <- winsorize_by_hand(x, lower_value, upper_value)

  out <- wt_trunc(w, method = "pctl", lower = 0.025, upper = 0.975)
  meta <- attr(out, "psw_trunc_meta")

  expect_identical(meta$method, "pctl")
  expect_identical(meta$lower, 0.025)
  expect_identical(meta$upper, 0.975)
  expect_equal(meta$lower_value, lower_value)
  expect_equal(meta$upper_value, upper_value)
  expect_identical(meta$truncated_idx, by_hand$idx)
  expect_length(meta$truncated_idx, 6)
  expect_equal(as.numeric(out), by_hand$values)
  expect_true(all(is.na(out[c(41, 92)])))
})

test_that("method = 'pctl' refuses an upper probability below one half and names both readings", {
  w <- pctl_weights()

  expect_error(
    wt_trunc(w, method = "pctl", upper = 0.1),
    class = "propensity_range_error"
  )
  message <- refusal_message(wt_trunc(w, method = "pctl", upper = 0.1))
  expect_match(message, "upper = 0.9", fixed = TRUE)
  expect_match(message, "lower = 0.1", fixed = TRUE)

  expect_propensity_error(wt_trunc(w, method = "pctl", upper = 0.1))
})

test_that("method = 'pctl' refuses probabilities outside their half of the unit interval", {
  w <- pctl_weights()

  expect_error(
    wt_trunc(w, method = "pctl"),
    class = "propensity_missing_arg_error"
  )
  expect_error(
    wt_trunc(w, method = "pctl", lower = 0.05),
    class = "propensity_missing_arg_error"
  )
  expect_error(
    wt_trunc(w, method = "pctl", upper = NA),
    class = "propensity_missing_value_error"
  )
  for (upper in list(1, 1.5, -0.2, c(0.9, 0.95))) {
    expect_error(
      wt_trunc(w, method = "pctl", upper = upper),
      class = "propensity_range_error"
    )
  }
  for (lower in list(0, 0.5, 0.6, -0.1)) {
    expect_error(
      wt_trunc(w, method = "pctl", lower = lower, upper = 0.95),
      class = "propensity_range_error"
    )
  }

  expect_propensity_error(wt_trunc(w, method = "pctl", upper = 1.5))
  expect_propensity_error(
    wt_trunc(w, method = "pctl", lower = 0.6, upper = 0.95)
  )
})

# ---- a count -----------------------------------------------------------------

test_that("method = 'count' bounds the largest weights at the next one down", {
  x <- c(3, 9, 1, 7, 2, 7.5, 4, NA, 0.5, 6)
  w <- psw(x, estimand = "ate")

  # The two largest are 9 and 7.5; the next one down is 7.
  out <- wt_trunc(w, method = "count", upper = 2)

  expect_equal(as.numeric(out), c(3, 7, 1, 7, 2, 7, 4, NA, 0.5, 6))
  expect_identical(
    attr(out, "psw_trunc_meta"),
    expected_wt_trunc_record(
      method = "count",
      upper = 2,
      upper_value = 7,
      truncated_idx = c(2L, 6L),
      n_obs = 10L
    )
  )
})

test_that("method = 'count' bounds the smallest weights at the next one up", {
  x <- c(3, 9, 1, 7, 2, 7.5, 4, NA, 0.5, 6)
  w <- psw(x, estimand = "ate")

  # The smallest is 0.5, the next one up is 1; the two largest go to 7.
  out <- wt_trunc(w, method = "count", lower = 1, upper = 2)

  expect_equal(as.numeric(out), c(3, 7, 1, 7, 2, 7, 4, NA, 1, 6))
  expect_identical(
    attr(out, "psw_trunc_meta"),
    expected_wt_trunc_record(
      method = "count",
      lower = 1,
      upper = 2,
      lower_value = 1,
      upper_value = 7,
      truncated_idx = c(2L, 6L, 9L),
      n_obs = 10L
    )
  )

  # Two smallest: 0.5 and 1 go to 2, the next one up.
  two_low <- wt_trunc(w, method = "count", lower = 2L, upper = 1L)
  expect_equal(as.numeric(two_low), c(3, 7.5, 2, 7, 2, 7.5, 4, NA, 2, 6))
  expect_identical(
    attr(two_low, "psw_trunc_meta")$truncated_idx,
    c(2L, 3L, 9L)
  )
  expect_identical(attr(two_low, "psw_trunc_meta")$lower, 2)
  expect_identical(attr(two_low, "psw_trunc_meta")$upper, 1)
})

test_that("method = 'count' moves only the weights above a tied bound", {
  # The largest weight ties with the next one down, so a count of one moves
  # nothing: the bound is 8 and no weight exceeds it.
  x <- c(4, 8, 8, 2, 5)
  w <- psw(x, estimand = "ate")

  one <- wt_trunc(w, method = "count", upper = 1)
  expect_equal(as.numeric(one), x)
  expect_identical(attr(one, "psw_trunc_meta")$upper_value, 8)
  expect_identical(attr(one, "psw_trunc_meta")$truncated_idx, integer())

  # A count of two reaches past the tie, to 5.
  two <- wt_trunc(w, method = "count", upper = 2)
  expect_equal(as.numeric(two), c(4, 5, 5, 2, 5))
  expect_identical(attr(two, "psw_trunc_meta")$upper_value, 5)
  expect_identical(attr(two, "psw_trunc_meta")$truncated_idx, c(2L, 3L))
})

test_that("method = 'count' refuses counts that meet or cross", {
  # Nine weights are present. Bounding the four smallest and the four largest
  # leaves one weight, which would have to be both the next one up and the
  # next one down.
  x <- c(3, 9, 1, 7, 2, 7.5, 4, NA, 0.5, 6)
  w <- psw(x, estimand = "ate")

  expect_error(
    wt_trunc(w, method = "count", lower = 4, upper = 4),
    class = "propensity_range_error"
  )
  expect_error(
    wt_trunc(w, method = "count", lower = 5, upper = 5),
    class = "propensity_range_error"
  )

  # Three smallest (0.5, 1, 2) go to 3; four largest (9, 7.5, 7, 6) go to 4.
  out <- expect_silent(wt_trunc(w, method = "count", lower = 3, upper = 4))
  expect_equal(as.numeric(out), c(3, 4, 3, 4, 3, 4, 4, NA, 3, 4))
  meta <- attr(out, "psw_trunc_meta")
  expect_identical(meta$lower_value, 3)
  expect_identical(meta$upper_value, 4)
  expect_identical(meta$truncated_idx, c(2L, 3L, 4L, 5L, 6L, 9L, 10L))

  expect_propensity_error(wt_trunc(w, method = "count", lower = 4, upper = 4))
})

test_that("method = 'count' refuses a count it cannot apply", {
  # Nine weights are present, so a count of nine leaves no next one down.
  w <- psw(c(3, 9, 1, 7, 2, 7.5, 4, NA, 0.5, 6), estimand = "ate")

  expect_error(
    wt_trunc(w, method = "count"),
    class = "propensity_missing_arg_error"
  )
  expect_error(
    wt_trunc(w, method = "count", upper = NA),
    class = "propensity_missing_value_error"
  )
  for (upper in list(0, -1, 1.5, 9, 10, c(1, 2), Inf)) {
    expect_error(
      wt_trunc(w, method = "count", upper = upper),
      class = "propensity_range_error"
    )
  }
  for (lower in list(0, 2.5, -1)) {
    expect_error(
      wt_trunc(w, method = "count", lower = lower, upper = 2),
      class = "propensity_range_error"
    )
  }
  expect_silent(wt_trunc(w, method = "count", upper = 8))

  expect_propensity_error(wt_trunc(w, method = "count", upper = 1.5))
  expect_propensity_error(wt_trunc(w, method = "count", upper = 9))
})

# ---- methods and arguments ---------------------------------------------------

test_that("wt_trunc() refuses a method it does not have", {
  w <- psw(c(1.2, 4, 1.5), estimand = "ate")

  # The score-scale methods belong to `ps_trunc()`.
  expect_error(wt_trunc(w, method = "ps", upper = 3), class = "rlang_error")
  expect_error(wt_trunc(w, method = "cr"), class = "rlang_error")
  expect_error(wt_trunc(w, method = "auto"), class = "rlang_error")
})

test_that("wt_trunc() refuses an argument it does not read", {
  # A misspelled bound would otherwise truncate at a bound the caller believed
  # they had replaced.
  w <- psw(c(1.2, 4, 1.5), estimand = "ate")

  expect_error(
    wt_trunc(w, method = "wt", uper = 3),
    class = "rlib_error_dots_nonempty"
  )
})

# ---- the input and its label -------------------------------------------------

test_that("the estimand gains a weights-truncated label", {
  dat <- wt_trunc_data()
  fit <- binary_wt_trunc_fit(dat)

  ate <- wt_trunc(wt_ate(fit, .exposure = dat$z), method = "wt", upper = 2)
  expect_identical(estimand(ate), "ate; weights truncated")

  att <- wt_trunc(wt_att(fit, .exposure = dat$z), method = "wt", upper = 2)
  expect_identical(estimand(att), "att; weights truncated")

  ato <- wt_trunc(wt_ato(fit, .exposure = dat$z), method = "wt", upper = 0.4)
  expect_identical(estimand(ato), "ato; weights truncated")
})

test_that("the label and the records stack on weights from a trimmed score", {
  dat <- wt_trunc_data(n = 200)
  fit <- binary_wt_trunc_fit(dat)
  trimmed <- ps_trim(fit, method = "ps", lower = 0.2, upper = 0.8)
  w <- wt_ate(ps_refit(trimmed, fit), .exposure = dat$z)
  expect_identical(estimand(w), "ate; trimmed")

  # Stacking a weight truncation on a trim is allowed, and says so.
  out <- expect_silent(wt_trunc(w, method = "wt", upper = 3))

  expect_identical(estimand(out), "ate; trimmed; weights truncated")
  expect_true(is_ps_trimmed(out))
  expect_true(is_wt_truncated(out))
  expect_false(is_ps_truncated(out))
  expect_identical(attr(out, "ps_trim_meta"), attr(w, "ps_trim_meta"))
  expect_identical(is_unit_trimmed(out), is_unit_trimmed(w))
  expect_identical(untouched_attrs(out), untouched_attrs(w))

  # The trimmed units stay missing and are not among those the bound moved.
  expect_identical(is.na(as.numeric(out)), is.na(as.numeric(w)))
  expect_false(any(is_unit_wt_truncated(out)[is.na(w)]))
  expect_equal(
    as.numeric(out),
    winsorize_by_hand(w, upper_value = 3)$values
  )
})

test_that("the label and the records stack on weights from a truncated score", {
  dat <- wt_trunc_data()
  fit <- binary_wt_trunc_fit(dat)
  truncated <- ps_trunc(fit, method = "ps", lower = 0.2, upper = 0.8)
  w <- wt_ate(truncated, .exposure = dat$z)
  expect_identical(estimand(w), "ate; truncated")

  out <- expect_silent(wt_trunc(w, method = "wt", upper = 3))

  expect_identical(estimand(out), "ate; truncated; weights truncated")
  expect_true(is_ps_truncated(out))
  expect_true(is_wt_truncated(out))
  expect_identical(is_unit_truncated(out), is_unit_truncated(w))
  expect_identical(untouched_attrs(out), untouched_attrs(w))
})

test_that("the label and the records stack on weights from a calibrated score", {
  dat <- wt_trunc_data(n = 200)
  ps <- stats::fitted(binary_wt_trunc_fit(dat))
  calibrated <- expect_silent(ps_calibrate(ps, dat$z))
  w <- wt_ate(calibrated, .exposure = dat$z)
  expect_identical(estimand(w), "ate; calibrated")

  out <- expect_silent(wt_trunc(w, method = "wt", upper = 3))

  expect_identical(estimand(out), "ate; calibrated; weights truncated")
  expect_true(is_ps_calibrated(out))
  expect_true(is_wt_truncated(out))
  expect_identical(untouched_attrs(out), untouched_attrs(w))
  expect_equal(
    as.numeric(out),
    winsorize_by_hand(w, upper_value = 3)$values
  )
})

test_that("a bare numeric vector becomes a truncated psw with no estimand", {
  x <- c(0.5, 1, 4, 2, 25)

  out <- expect_silent(wt_trunc(x, method = "wt", upper = 20))

  expect_s3_class(out, "psw")
  expect_null(estimand(out))
  expect_true(is_wt_truncated(out))
  expect_equal(as.numeric(out), c(0.5, 1, 4, 2, 20))
  expect_identical(
    attr(out, "psw_trunc_meta"),
    expected_wt_trunc_record(
      method = "wt",
      upper = 20,
      upper_value = 20,
      truncated_idx = 5L,
      n_obs = 5L
    )
  )

  adaptive <- wt_trunc(c(rep(1.5, 16), 0.01, 2, 5, 9))
  expect_s3_class(adaptive, "psw")
  expect_null(estimand(adaptive))
})

test_that("the names on the weights are kept", {
  unit_names <- paste0("u", 1:5)
  x <- stats::setNames(c(0.5, 1, 4, 2, 25), unit_names)
  w <- psw(x, estimand = "ate")
  names(w) <- unit_names

  from_psw <- list(
    wt = wt_trunc(w, method = "wt", upper = 3),
    count = wt_trunc(w, method = "count", upper = 1)
  )
  from_numeric <- list(
    wt = wt_trunc(x, method = "wt", upper = 3),
    count = wt_trunc(x, method = "count", upper = 1)
  )

  for (method in names(from_psw)) {
    expect_identical(names(from_psw[[method]]), unit_names, info = method)
    expect_identical(names(from_numeric[[method]]), unit_names, info = method)
  }
  expect_equal(unname(as.numeric(from_psw$wt)), c(0.5, 1, 3, 2, 3))
  expect_equal(unname(as.numeric(from_numeric$count)), c(0.5, 1, 4, 2, 4))

  # Integer input still becomes doubles when it is named.
  named_integer <- wt_trunc(
    c(a = 1L, b = 5L, c = 30L),
    method = "wt",
    upper = 2.5
  )
  expect_type(vctrs::vec_data(named_integer), "double")
  expect_identical(names(named_integer), c("a", "b", "c"))
})

test_that("a subset that drops the record prints a truncation line without counts", {
  w <- psw(c(0.5, 1, 4, 2, 25), estimand = "ate")

  expect_snapshot(vctrs::vec_slice(wt_trunc(w, method = "wt", upper = 3), 1:3))
})

test_that("a subset re-indexes the record and prints the truncation counts", {
  w <- psw(c(0.5, 1, 4, 2, 25), estimand = "ate")

  expect_snapshot(wt_trunc(w, method = "wt", upper = 3)[c(5, 1, 3)])
})

test_that("wt_trunc() refuses input that is not weights", {
  expect_error(wt_trunc(c("a", "b")), class = "propensity_type_error")
  expect_error(wt_trunc(list(1, 2)), class = "propensity_error")

  expect_propensity_error(wt_trunc(c("a", "b")))
})

test_that("wt_trunc() refuses a matrix of weights", {
  weights <- matrix(c(1, 2, 3, 40), nrow = 2)

  expect_error(
    wt_trunc(weights, method = "wt", upper = 5),
    class = "propensity_type_error"
  )
  expect_error(wt_trunc(weights), class = "propensity_type_error")

  expect_propensity_error(wt_trunc(weights, method = "wt", upper = 5))
})

test_that("wt_trunc() refuses propensity scores", {
  # A trimmed, truncated, or calibrated score is a score rather than a weight.
  # The refusal points to building weights from it first, or to bounding the
  # score itself with `ps_trunc()`.
  dat <- wt_trunc_data(n = 200)
  ps <- stats::fitted(binary_wt_trunc_fit(dat))
  scores <- list(
    ps_trim = ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9),
    ps_trunc = ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9),
    ps_calib = ps_calibrate(ps, dat$z)
  )

  for (label in names(scores)) {
    score <- scores[[label]]
    expect_s3_class(score, label)
    expect_error(
      wt_trunc(score, method = "wt", upper = 3),
      class = "propensity_type_error",
      info = label
    )
  }

  message <- refusal_message(
    wt_trunc(scores$ps_trim, method = "wt", upper = 3)
  )
  expect_match(message, "wt_ate()", fixed = TRUE)
  expect_match(message, "ps_trunc()", fixed = TRUE)

  trimmed <- scores$ps_trim
  expect_propensity_error(wt_trunc(trimmed, method = "wt", upper = 3))
})

test_that("the data-driven methods refuse fewer than two weights", {
  inputs <- list(
    empty = psw(double(), estimand = "ate"),
    all_missing = psw(c(NA_real_, NA_real_), estimand = "ate"),
    one_present = psw(c(NA, 2.5, NA), estimand = "ate")
  )

  for (label in names(inputs)) {
    w <- inputs[[label]]
    expect_error(
      wt_trunc(w),
      class = "propensity_range_error",
      info = label
    )
    expect_error(
      wt_trunc(w, method = "pctl", upper = 0.95),
      class = "propensity_range_error",
      info = label
    )
  }

  # A count must stay below the number of weights present, which none are.
  expect_error(
    wt_trunc(inputs$all_missing, method = "count", upper = 1),
    class = "propensity_range_error"
  )

  one_present <- inputs$one_present
  expect_propensity_error(wt_trunc(one_present))
})

test_that("an absolute bound on no weights leaves an empty record", {
  out <- expect_silent(
    wt_trunc(psw(double(), estimand = "ate"), method = "wt", upper = 3)
  )

  expect_s3_class(out, "psw")
  expect_length(out, 0)
  expect_true(is_wt_truncated(out))
  expect_identical(estimand(out), "ate; weights truncated")
  meta <- attr(out, "psw_trunc_meta")
  expect_identical(meta$truncated_idx, integer())
  expect_identical(meta$n_obs, 0L)
})

test_that("truncating truncated weights warns and returns them unchanged", {
  w <- psw(c(1.2, 4, 1.5, 2.1, 3.5), estimand = "ate")
  once <- wt_trunc(w, method = "wt", upper = 3)

  expect_warning(
    twice <- wt_trunc(once, method = "wt", upper = 2),
    class = "propensity_already_modified_warning"
  )
  expect_identical(twice, once)

  expect_warning(
    again <- wt_trunc(once),
    class = "propensity_already_modified_warning"
  )
  expect_identical(again, once)

  expect_propensity_warning(wt_trunc(once, method = "wt", upper = 2))
})

# ---- what the truncation keeps -----------------------------------------------

test_that("truncating binary weights keeps their records", {
  w <- stabilized_binary_wt_trunc_weights()
  expect_true(is_stabilized(w))
  expect_false(is.null(numerator_model(w)))
  bound <- stats::median(as.numeric(w))

  out <- wt_trunc(w, method = "wt", upper = bound)

  expect_identical(untouched_attrs(out), untouched_attrs(w))
  expect_identical(numerator_model(out), numerator_model(w))
  expect_true(is_stabilized(out))
  expect_identical(exposure_type(out), "binary")
  expect_identical(estimand(out), "ate; weights truncated")
  expect_equal(
    as.numeric(out),
    winsorize_by_hand(w, upper_value = bound)$values
  )
})

test_that("truncating categorical weights keeps their records", {
  w <- categorical_wt_trunc_weights()
  bound <- stats::quantile(as.numeric(w), 0.9, names = FALSE)
  by_hand <- winsorize_by_hand(w, upper_value = bound)
  expect_gt(length(by_hand$idx), 0)

  out <- wt_trunc(w, method = "wt", upper = bound)

  expect_identical(untouched_attrs(out), untouched_attrs(w))
  expect_identical(exposure_type(out), "categorical")
  expect_identical(attr(out, "category_names"), attr(w, "category_names"))
  expect_identical(estimand(out), "ate; weights truncated")
  expect_equal(as.numeric(out), by_hand$values)
  expect_identical(attr(out, "psw_trunc_meta")$truncated_idx, by_hand$idx)
})

test_that("truncating continuous weights keeps the density record and numerator model", {
  w <- continuous_wt_trunc_weights()
  expect_false(is.null(density_meta(w)))
  expect_s3_class(numerator_model(w), "lm")
  bound <- stats::quantile(as.numeric(w), 0.9, names = FALSE)
  by_hand <- winsorize_by_hand(w, upper_value = bound)

  out <- wt_trunc(w, method = "wt", upper = bound)

  expect_identical(density_meta(out), density_meta(w))
  expect_identical(numerator_model(out), numerator_model(w))
  expect_identical(untouched_attrs(out), untouched_attrs(w))
  expect_true(is_stabilized(out))
  expect_identical(exposure_type(out), "continuous")
  expect_identical(estimand(out), "ate; weights truncated")
  expect_equal(as.numeric(out), by_hand$values)

  # Every method keeps them, the default included.
  for (method in c("adaptive", "pctl", "count")) {
    upper <- switch(method, adaptive = NULL, pctl = 0.95, count = 3)
    kept <- wt_trunc(w, method = method, upper = upper)
    expect_identical(untouched_attrs(kept), untouched_attrs(w), info = method)
  }
})

test_that("the count method moves the largest weights from a real fit", {
  w <- binary_wt_trunc_weights()
  x <- as.numeric(w)
  ranked <- sort(x, decreasing = TRUE)
  expect_gt(ranked[3], ranked[4])

  out <- wt_trunc(w, method = "count", upper = 3)
  meta <- attr(out, "psw_trunc_meta")

  expect_identical(meta$upper_value, ranked[4])
  expect_identical(meta$truncated_idx, sort(order(x, decreasing = TRUE)[1:3]))
  expect_equal(
    as.numeric(out),
    winsorize_by_hand(x, upper_value = ranked[4])$values
  )
})

# ---- the printed record ------------------------------------------------------

# The footer prints each bound to three significant digits, as 43.69 prints
# as 43.7.
footer_number <- function(x) {
  format(x, digits = 3)
}

test_that("the footer reports the truncation", {
  # A thousand weights put the adaptive bound at 43.69, which two exceed.
  x <- c(rep(1, 997), 50, 60, 2)
  w <- psw(x, estimand = "ate")
  expect_equal(sqrt(1000) * log(1000) / 5, 43.69, tolerance = 1e-4)
  expect_identical(footer_number(sqrt(1000) * log(1000) / 5), "43.7")

  lines <- capture_output_lines(print(wt_trunc(w)))

  expect_identical(
    lines[[length(lines)]],
    paste0(
      "truncation: adaptive (upper ",
      footer_number(sqrt(1000) * log(1000) / 5),
      "), 2 of 1000 weights truncated"
    )
  )
})

test_that("the truncation line follows the density record", {
  w <- continuous_wt_trunc_weights()
  untruncated <- capture_output_lines(print(w))
  expect_false(any(grepl("^truncation:", untruncated)))

  upper_value <- sort(as.numeric(w), decreasing = TRUE)[4]
  lines <- capture_output_lines(print(wt_trunc(w, method = "count", upper = 3)))
  footer <- lines[seq(length(lines) - 4, length(lines))]

  expect_match(footer[[1]], "^density:")
  expect_match(footer[[4]], "^stabilize:")
  expect_identical(
    footer[[5]],
    paste0(
      "truncation: count 3 (upper ",
      footer_number(upper_value),
      "), 3 of 60 weights truncated"
    )
  )
})

test_that("the footer names each method's argument and its realized bounds", {
  last_line <- function(x) {
    lines <- capture_output_lines(print(x))
    lines[[length(lines)]]
  }

  w <- psw(c(0.5, 1, 4, 2, 25, NA, 0.02), estimand = "ate")
  expect_identical(
    last_line(wt_trunc(w, method = "wt", upper = 3)),
    "truncation: wt (upper 3), 2 of 7 weights truncated"
  )
  expect_identical(
    last_line(wt_trunc(w, method = "wt", lower = 0.05, upper = 3)),
    "truncation: wt (lower 0.05, upper 3), 3 of 7 weights truncated"
  )

  p <- pctl_weights()
  present <- sort(as.numeric(p)[!is.na(p)])
  low <- present[3] + 0.5 * (present[4] - present[3])
  high <- present[98] + 0.5 * (present[99] - present[98])
  expect_identical(
    last_line(wt_trunc(p, method = "pctl", upper = 0.99)),
    paste0(
      "truncation: pctl 0.99 (upper ",
      footer_number(present[100]),
      "), 1 of 103 weights truncated"
    )
  )
  expect_identical(
    last_line(wt_trunc(p, method = "pctl", lower = 0.025, upper = 0.975)),
    paste0(
      "truncation: pctl 0.025/0.975 (lower ",
      footer_number(low),
      ", upper ",
      footer_number(high),
      "), 6 of 103 weights truncated"
    )
  )

  counted <- psw(c(3, 9, 1, 7, 2, 7.5, 4, NA, 0.5, 6), estimand = "ate")
  expect_identical(
    last_line(wt_trunc(counted, method = "count", upper = 3)),
    "truncation: count 3 (upper 6), 3 of 10 weights truncated"
  )
  expect_identical(
    last_line(wt_trunc(counted, method = "count", lower = 2, upper = 3)),
    "truncation: count 2/3 (lower 2, upper 6), 5 of 10 weights truncated"
  )
})

test_that("truncated weights print their bound and count", {
  w <- psw(c(0.5, 1, 4, 2, 25, NA, 0.02), estimand = "ate")

  expect_snapshot({
    wt_trunc(w, method = "wt", upper = 20)
    wt_trunc(w, method = "wt", lower = 0.05, upper = 20)
    wt_trunc(c(0.5, 1, 4, 2, 25), method = "count", upper = 1)
  })

  dat <- wt_trunc_data(n = 12)
  dose <- stats::lm(a ~ x1, data = dat)
  weights <- muffle_variance_warning(wt_ate(dose))
  expect_snapshot(wt_trunc(weights, method = "pctl", upper = 0.9))
})

# ---- agreement with WeightIt -------------------------------------------------

# WeightIt's `trim()` winsorizes a weight vector as well, and its conventions
# match this verb's in two places only, which the two tests below use:
#
# * A count of the largest weights (`at` >= 1, `lower = FALSE`). Both sort the
#   weights and bound the `at` largest at the next one down, moving only those
#   strictly above it. WeightIt folds the count, `at <- min(at, n - at)`, where
#   this verb refuses a count that leaves no weight below it; WeightIt warns and
#   returns the weights unchanged for `at >= n`; its `lower = TRUE` trims the
#   same count at both ends, where this verb takes a separate `lower`; and its
#   `n` includes missing weights, where this verb counts the weights present.
# * An upper percentile where the two quantile rules pick the same order
#   statistic. WeightIt reads quantiles with `type = 3`; this verb follows
#   `ps_trunc(method = "pctl")`, which uses the `quantile()` default,
#   `type = 7`. At `n = 101` and `upper = 0.99` both select the 100th order
#   statistic. Where type 7 lands on an order statistic in the lower tail,
#   type 3 selects the one below it, so the two never agree on a lower bound
#   there. WeightIt also folds `at <- max(at, 1 - at)`, reading `at = 0.01` as
#   the upper tail, where this verb refuses an `upper` below one half; and in
#   WeightIt 2.0.0 its `quantile()` call errors on a missing weight.
#
# WeightIt has no absolute bound (an `at` of 1 or more is a count) and no
# Gruber bound, so `"wt"` and `"adaptive"` have no oracle there. WeightIt
# reports what it did with a message, which is muffled around the oracle only.

weightit_trim <- function(x, ...) {
  as.numeric(suppressMessages(WeightIt::trim(as.numeric(x), ...)))
}

test_that("a count of the largest weights agrees with WeightIt::trim()", {
  skip_if_not_installed("WeightIt")
  w <- binary_wt_trunc_weights()
  expect_false(anyNA(w))

  for (k in c(1, 3, 5)) {
    expect_equal(
      as.numeric(wt_trunc(w, method = "count", upper = k)),
      weightit_trim(w, at = k),
      info = k
    )
  }
})

test_that("an upper percentile agrees with WeightIt::trim() where the quantile rules coincide", {
  skip_if_not_installed("WeightIt")
  dat <- wt_trunc_data(n = 101)
  w <- continuous_wt_trunc_weights(dat)
  expect_false(anyNA(w))
  expect_length(w, 101)
  x <- as.numeric(w)
  expect_identical(
    stats::quantile(x, 0.99, type = 3, names = FALSE),
    sort(x)[100]
  )
  expect_identical(type7_quantile(x, 0.99), sort(x)[100])

  expect_equal(
    as.numeric(wt_trunc(w, method = "pctl", upper = 0.99)),
    weightit_trim(w, at = 0.99)
  )
})
