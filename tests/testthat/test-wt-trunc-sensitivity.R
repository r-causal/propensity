# The truncation sensitivity grid ----------------------------------------------

# `wt_trunc_sensitivity()` reports, for each bound in a grid, what `wt_trunc()`
# does to the weights at that bound: the bounds as given and as applied, the
# number of weights moved, and the range and mean of the weights that result.
# The first row is the untruncated reference, with missing bounds. There is no
# effective sample size column; that comes from `halfmoon::check_ess()`.

grid_columns <- c(
  "lower",
  "upper",
  "lower_value",
  "upper_value",
  "n_truncated",
  "min",
  "max",
  "range_ratio",
  "mean"
)

# One row of the grid, built from the weights and the bounds that produced
# them. A missing weight takes no part in the summaries.
grid_row <- function(
  values,
  lower = NA_real_,
  upper = NA_real_,
  lower_value = NA_real_,
  upper_value = NA_real_,
  n_truncated = 0L
) {
  values <- as.numeric(values)
  values <- values[!is.na(values)]
  tibble::tibble(
    lower = as.double(lower),
    upper = as.double(upper),
    lower_value = as.double(lower_value),
    upper_value = as.double(upper_value),
    n_truncated = as.integer(n_truncated),
    min = min(values),
    max = max(values),
    range_ratio = max(values) / min(values),
    mean = mean(values)
  )
}

null_to_na <- function(x) {
  if (is.null(x)) NA_real_ else x
}

# The row `wt_trunc()` implies at one bound, read from the weights it returns
# and the record it leaves.
grid_row_from_wt_trunc <- function(.weights, method, lower, upper) {
  out <- wt_trunc(.weights, method = method, lower = lower, upper = upper)
  meta <- attr(out, "psw_trunc_meta")
  grid_row(
    out,
    lower = null_to_na(lower),
    upper = upper,
    lower_value = null_to_na(meta$lower_value),
    upper_value = meta$upper_value,
    n_truncated = length(meta$truncated_idx)
  )
}

# The whole grid `wt_trunc()` implies: the untruncated row, then one row per
# bound in the order given, with a single `lower` paired with every `upper`.
grid_from_wt_trunc <- function(.weights, method, upper, lower = NULL) {
  if (length(lower) == 1L) {
    lower <- rep(lower, length(upper))
  }
  rows <- lapply(seq_along(upper), function(i) {
    grid_row_from_wt_trunc(.weights, method, lower[i], upper[[i]])
  })
  vctrs::vec_rbind(grid_row(.weights), !!!rows)
}

# Quantile type 7, the `quantile()` default that `wt_trunc(method = "pctl")`
# uses, written out.
type7_quantile <- function(x, p) {
  x <- sort(x[!is.na(x)])
  h <- (length(x) - 1) * p + 1
  lo <- floor(h)
  hi <- ceiling(h)
  x[lo] + (h - lo) * (x[hi] - x[lo])
}

# Ten weights whose truncated summaries are simple enough to work out by hand.
# Their sum is 30, so the untruncated mean is 3.
hand_weights <- function() {
  psw(c(0.5, 0.8, 1, 1.2, 1.5, 2, 3, 4, 6, 10), estimand = "ate")
}

sensitivity_data <- function(n = 80, seed = 20260917) {
  withr::local_seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  dat$z <- rbinom(n, 1, stats::plogis(1.4 * dat$x1))
  dat$a <- 1 + 0.6 * dat$x1 - 0.4 * dat$x2 + rnorm(n)
  dat
}

binary_sensitivity_fit <- function(dat = sensitivity_data()) {
  stats::glm(z ~ x1 + x2, data = dat, family = stats::binomial())
}

binary_sensitivity_weights <- function(dat = sensitivity_data()) {
  wt_ate(binary_sensitivity_fit(dat), .exposure = dat$z)
}

continuous_sensitivity_weights <- function(dat = sensitivity_data()) {
  dose <- stats::lm(a ~ x1 + x2, data = dat)
  numerator <- stats::lm(a ~ x1, data = dat)
  muffle_variance_warning(wt_ate(dose, stabilize = numerator))
}

# ---- the shape of the grid ---------------------------------------------------

test_that("the grid is a tibble with the documented columns in order", {
  out <- wt_trunc_sensitivity(hand_weights())

  expect_s3_class(out, "tbl_df")
  expect_identical(names(out), grid_columns)
  expect_false("ess" %in% tolower(names(out)))
  expect_identical(
    vapply(out, typeof, character(1)),
    c(
      lower = "double",
      upper = "double",
      lower_value = "double",
      upper_value = "double",
      n_truncated = "integer",
      min = "double",
      max = "double",
      range_ratio = "double",
      mean = "double"
    )
  )
  # The summaries are plain numbers, not weights.
  for (column in c("min", "max", "range_ratio", "mean")) {
    expect_false(inherits(out[[column]], "psw"), info = column)
  }
})

test_that("the first row is the untruncated reference with missing bounds", {
  # The reference comes first so that every truncated row below it reads as a
  # departure from it.
  w <- hand_weights()

  out <- wt_trunc_sensitivity(w, method = "wt", upper = c(5, 2.5))

  expect_identical(nrow(out), 3L)
  expect_identical(out$lower[[1]], NA_real_)
  expect_identical(out$upper[[1]], NA_real_)
  expect_identical(out$lower_value[[1]], NA_real_)
  expect_identical(out$upper_value[[1]], NA_real_)
  expect_identical(out$n_truncated[[1]], 0L)
  expect_equal(out$min[[1]], 0.5)
  expect_equal(out$max[[1]], 10)
  expect_equal(out$range_ratio[[1]], 20)
  expect_equal(out$mean[[1]], 3)
  # Only the reference row has missing bounds.
  expect_false(anyNA(out$upper[-1]))
})

test_that("the rows follow the grid in the order given", {
  w <- hand_weights()

  out <- wt_trunc_sensitivity(w, method = "wt", upper = c(2.5, 5, 4))

  expect_identical(out$upper, c(NA, 2.5, 5, 4))
})

# ---- the values in each row --------------------------------------------------

test_that("an absolute grid gives the summaries worked out by hand", {
  # At 5, the weights 6 and 10 are moved to 5: the sum becomes 24. At 2.5, the
  # weights 3, 4, 6, and 10 are moved to 2.5: the sum becomes 17.
  w <- hand_weights()

  out <- wt_trunc_sensitivity(w, method = "wt", upper = c(5, 2.5))

  expected <- tibble::tibble(
    lower = c(NA_real_, NA_real_, NA_real_),
    upper = c(NA, 5, 2.5),
    lower_value = c(NA_real_, NA_real_, NA_real_),
    upper_value = c(NA, 5, 2.5),
    n_truncated = c(0L, 2L, 4L),
    min = c(0.5, 0.5, 0.5),
    max = c(10, 5, 2.5),
    range_ratio = c(20, 10, 5),
    mean = c(3, 2.4, 1.7)
  )
  expect_equal(out, expected)
})

test_that("a two-sided absolute grid gives the summaries worked out by hand", {
  # A lower bound of 1 also moves 0.5 and 0.8 up to 1. With the upper bound at
  # 5 the sum is 24.7.
  w <- hand_weights()

  out <- wt_trunc_sensitivity(w, method = "wt", lower = 1, upper = 5)

  expect_equal(out$lower, c(NA, 1))
  expect_equal(out$lower_value, c(NA, 1))
  expect_identical(out$n_truncated, c(0L, 4L))
  expect_equal(out$min, c(0.5, 1))
  expect_equal(out$max, c(10, 5))
  expect_equal(out$range_ratio, c(20, 5))
  expect_equal(out$mean, c(3, 2.47))
})

test_that("a percentile grid reads its bounds at the type 7 quantiles", {
  # For ten weights the 0.9 quantile sits a tenth of the way from 6 to 10, at
  # 6.4, and the 0.99 quantile 91 hundredths of the way, at 9.64. Each moves
  # only the largest weight.
  w <- hand_weights()
  expect_equal(type7_quantile(as.numeric(w), 0.9), 6.4)
  expect_equal(type7_quantile(as.numeric(w), 0.99), 9.64)

  out <- wt_trunc_sensitivity(w, method = "pctl", upper = c(0.99, 0.9))

  expect_equal(out$upper, c(NA, 0.99, 0.9))
  expect_equal(out$upper_value, c(NA, 9.64, 6.4))
  expect_identical(out$n_truncated, c(0L, 1L, 1L))
  expect_equal(out$max, c(10, 9.64, 6.4))
  expect_equal(out$range_ratio, c(20, 19.28, 12.8))
  expect_equal(out$mean, c(3, 2.964, 2.64))
})

test_that("a count grid bounds the largest weights at the next one down", {
  # The two largest, 6 and 10, are bounded at 4, the third largest.
  w <- hand_weights()

  out <- wt_trunc_sensitivity(w, method = "count", upper = 2)

  expect_equal(out$upper, c(NA, 2))
  expect_equal(out$upper_value, c(NA, 4))
  expect_identical(out$n_truncated, c(0L, 2L))
  expect_equal(out$max, c(10, 4))
  expect_equal(out$range_ratio, c(20, 8))
  expect_equal(out$mean, c(3, 2.2))
})

test_that("every row agrees with wt_trunc() at its bound", {
  w <- binary_sensitivity_weights()
  grids <- list(
    wt = list(upper = c(8, 5, 3, 2), lower = NULL),
    wt_two_sided = list(upper = c(8, 5, 3), lower = 1.05),
    pctl = list(upper = c(0.99, 0.95, 0.8), lower = c(0.01, 0.05, 0.2)),
    count = list(upper = c(1, 3, 10), lower = 2)
  )

  for (label in names(grids)) {
    method <- sub("_two_sided", "", label, fixed = TRUE)
    grid <- grids[[label]]
    out <- wt_trunc_sensitivity(
      w,
      method = method,
      lower = grid$lower,
      upper = grid$upper
    )
    expected <- grid_from_wt_trunc(w, method, grid$upper, grid$lower)
    expect_equal(out, expected, info = label)
    # Each grid point moves at least one weight, so no row is the reference
    # repeated.
    expect_true(all(out$n_truncated[-1] > 0), info = label)
  }
})

test_that("an integer grid is recorded as doubles", {
  w <- hand_weights()

  out <- wt_trunc_sensitivity(w, method = "wt", lower = 1L, upper = c(5L, 3L))

  expect_identical(out$upper, c(NA, 5, 3))
  expect_identical(out$lower, c(NA, 1, 1))
})

# ---- the default grid --------------------------------------------------------

test_that("the default is a percentile grid at 0.99, 0.975, 0.95, and 0.90", {
  w <- binary_sensitivity_weights()
  default_grid <- c(0.99, 0.975, 0.95, 0.90)

  out <- wt_trunc_sensitivity(w)

  expect_identical(out$upper, c(NA, default_grid))
  expect_identical(out$lower, rep(NA_real_, 5))
  expect_identical(out, wt_trunc_sensitivity(w, method = "pctl"))
  expect_identical(
    out,
    wt_trunc_sensitivity(w, method = "pctl", upper = default_grid)
  )
  expect_equal(out, grid_from_wt_trunc(w, "pctl", default_grid))
})

# ---- pairing lower with upper ------------------------------------------------

test_that("a single lower bound is paired with every upper bound", {
  w <- binary_sensitivity_weights()

  out <- wt_trunc_sensitivity(w, lower = 0.05)

  expect_identical(out$lower, c(NA, rep(0.05, 4)))
  expect_equal(
    out,
    grid_from_wt_trunc(w, "pctl", c(0.99, 0.975, 0.95, 0.90), lower = 0.05)
  )
})

test_that("a lower grid of matching length is paired element by element", {
  w <- binary_sensitivity_weights()
  lower <- c(0.01, 0.025, 0.05, 0.10)
  upper <- c(0.99, 0.975, 0.95, 0.90)

  out <- wt_trunc_sensitivity(w, method = "pctl", lower = lower, upper = upper)

  expect_identical(out$lower, c(NA, lower))
  expect_identical(out$upper, c(NA, upper))
  for (i in seq_along(upper)) {
    expect_equal(
      out[i + 1, ],
      grid_row_from_wt_trunc(w, "pctl", lower[[i]], upper[[i]]),
      info = i
    )
  }
})

test_that("a lower grid of any other length is refused", {
  w <- hand_weights()

  expect_error(
    wt_trunc_sensitivity(w, lower = c(0.01, 0.05)),
    class = "propensity_length_error"
  )
  # `lower` is recycled to the length of `upper`, never the other way round.
  expect_error(
    wt_trunc_sensitivity(
      w,
      method = "wt",
      lower = c(0.6, 0.7),
      upper = 5
    ),
    class = "propensity_length_error"
  )

  expect_propensity_error(wt_trunc_sensitivity(w, lower = c(0.01, 0.05)))
})

# ---- the bounds refused ------------------------------------------------------

test_that("the absolute and count grids need upper, as wt_trunc() does", {
  w <- hand_weights()

  for (method in c("wt", "count")) {
    expect_error(
      wt_trunc_sensitivity(w, method = method),
      class = "propensity_missing_arg_error",
      info = method
    )
    expect_error(
      wt_trunc_sensitivity(w, method = method, lower = 1),
      class = "propensity_missing_arg_error",
      info = method
    )
  }

  expect_propensity_error(wt_trunc_sensitivity(w, method = "wt"))
})

test_that("one invalid bound refuses the whole grid with wt_trunc()'s class", {
  # A grid is refused as a whole rather than reported with a missing row, so a
  # typo in one bound cannot pass unnoticed in a table of many.
  w <- hand_weights()

  refused <- list(
    list(method = "pctl", lower = NULL, upper = c(0.99, 1.5)),
    list(method = "pctl", lower = NULL, upper = c(0.95, 0.1)),
    list(method = "pctl", lower = c(0.05, 0.6), upper = c(0.99, 0.95)),
    list(method = "wt", lower = NULL, upper = c(5, -1)),
    list(method = "wt", lower = NULL, upper = c(5, Inf)),
    list(method = "wt", lower = c(1, 6), upper = c(5, 5)),
    list(method = "count", lower = NULL, upper = c(2, 1.5)),
    list(method = "count", lower = NULL, upper = c(2, 10)),
    list(method = "count", lower = 5, upper = c(2, 4))
  )
  for (i in seq_along(refused)) {
    args <- refused[[i]]
    expect_error(
      wt_trunc_sensitivity(
        w,
        method = args$method,
        lower = args$lower,
        upper = args$upper
      ),
      class = "propensity_range_error",
      info = i
    )
  }

  expect_error(
    wt_trunc_sensitivity(w, method = "wt", upper = c(5, NA)),
    class = "propensity_missing_value_error"
  )
  expect_error(
    wt_trunc_sensitivity(w, method = "wt", lower = c(1, NA), upper = c(5, 5)),
    class = "propensity_missing_value_error"
  )
  expect_error(
    wt_trunc_sensitivity(w, method = "wt", upper = c("5", "3")),
    class = "propensity_type_error"
  )

  expect_propensity_error(
    wt_trunc_sensitivity(w, method = "pctl", upper = c(0.99, 1.5))
  )
  expect_propensity_error(
    wt_trunc_sensitivity(w, method = "wt", upper = c(5, NA))
  )
})

test_that("an empty grid is refused", {
  w <- hand_weights()

  expect_error(
    wt_trunc_sensitivity(w, method = "wt", upper = double()),
    class = "propensity_length_error"
  )
  expect_error(
    wt_trunc_sensitivity(w, method = "pctl", upper = double()),
    class = "propensity_length_error"
  )
  expect_error(
    wt_trunc_sensitivity(w, method = "count", lower = 1, upper = integer()),
    class = "propensity_length_error"
  )

  expect_propensity_error(
    wt_trunc_sensitivity(w, method = "wt", upper = double())
  )
})

test_that("a misspelled argument is an error", {
  # There is no `...`, so a misspelled bound cannot be silently ignored.
  w <- hand_weights()

  expect_error(wt_trunc_sensitivity(w, uper = 0.9))
  expect_error(wt_trunc_sensitivity(w, method = "wt", upper = 5, lowr = 1))
})

test_that("the adaptive method is not a grid method", {
  # The adaptive rule fixes a single bound from the number of weights, so there
  # is nothing to vary over.
  w <- hand_weights()

  expect_error(
    wt_trunc_sensitivity(w, method = "adaptive"),
    class = "rlang_error"
  )
  expect_error(
    wt_trunc_sensitivity(w, method = "ps", upper = 0.9),
    class = "rlang_error"
  )

  expect_propensity_error(wt_trunc_sensitivity(w, method = "adaptive"))
})

# ---- the weights accepted ----------------------------------------------------

test_that("a bare numeric vector gives the same grid as a psw", {
  x <- as.numeric(hand_weights())

  out <- wt_trunc_sensitivity(x, method = "wt", upper = c(5, 2.5))

  expect_identical(
    out,
    wt_trunc_sensitivity(hand_weights(), method = "wt", upper = c(5, 2.5))
  )
  expect_equal(out, grid_from_wt_trunc(x, "wt", c(5, 2.5)))
  expect_identical(wt_trunc_sensitivity(x), wt_trunc_sensitivity(psw(x)))
})

test_that("missing weights are left out of every row", {
  # The present weights are the hand weights, so every summary matches theirs,
  # and the missing weights are never counted as moved.
  x <- c(0.5, NA, 0.8, 1, 1.2, 1.5, NA, 2, 3, 4, 6, 10)
  w <- psw(x, estimand = "ate")

  out <- expect_silent(
    wt_trunc_sensitivity(w, method = "wt", upper = c(5, 2.5))
  )

  expect_equal(
    out,
    wt_trunc_sensitivity(hand_weights(), method = "wt", upper = c(5, 2.5))
  )
  expect_identical(out$n_truncated, c(0L, 2L, 4L))
  expect_false(anyNA(out[c("min", "max", "range_ratio", "mean")]))
  expect_equal(out, grid_from_wt_trunc(w, "wt", c(5, 2.5)))

  # A percentile grid reads its quantiles from the present weights alone.
  pctl <- wt_trunc_sensitivity(w, method = "pctl", upper = 0.9)
  expect_equal(pctl$upper_value, c(NA, 6.4))
  expect_equal(pctl, grid_from_wt_trunc(w, "pctl", 0.9))
})

test_that("weights from a trimmed score give a grid over the kept units", {
  dat <- sensitivity_data()
  fit <- binary_sensitivity_fit(dat)
  trimmed <- ps_refit(
    ps_trim(fit, method = "ps", lower = 0.05, upper = 0.95),
    fit
  )
  w <- wt_ate(trimmed, .exposure = dat$z)
  expect_true(anyNA(w))

  out <- expect_silent(wt_trunc_sensitivity(w))

  expect_equal(out, grid_from_wt_trunc(w, "pctl", c(0.99, 0.975, 0.95, 0.90)))
  expect_equal(out$mean[[1]], mean(as.numeric(w), na.rm = TRUE))
})

test_that("a pctl grid refuses fewer than two weights present, as wt_trunc() does", {
  w <- psw(c(NA, 2.5, NA), estimand = "ate")

  expect_error(wt_trunc_sensitivity(w), class = "propensity_range_error")
})

test_that("a grid over no weights is refused under every method", {
  # With no weight present there is no range or mean to report, so the grid is
  # refused before any summary is taken, whatever the method.
  inputs <- list(
    all_missing = psw(c(NA_real_, NA_real_), estimand = "ate"),
    empty = psw(double(), estimand = "ate"),
    empty_numeric = double()
  )
  grids <- list(
    pctl = NULL,
    wt = 5,
    count = 1
  )

  for (label in names(inputs)) {
    for (method in names(grids)) {
      expect_error(
        wt_trunc_sensitivity(
          inputs[[label]],
          method = method,
          upper = grids[[method]]
        ),
        class = "propensity_range_error",
        info = paste(label, method)
      )
    }
  }

  empty <- inputs$empty
  expect_propensity_error(wt_trunc_sensitivity(empty, method = "wt", upper = 5))
})

test_that("continuous weights with a density record give the grid wt_trunc() implies", {
  w <- continuous_sensitivity_weights()
  expect_false(is.null(density_meta(w)))
  expect_identical(exposure_type(w), "continuous")

  out <- expect_silent(wt_trunc_sensitivity(w))
  expect_equal(out, grid_from_wt_trunc(w, "pctl", c(0.99, 0.975, 0.95, 0.90)))

  by_weight <- wt_trunc_sensitivity(w, method = "wt", upper = c(3, 2, 1.5))
  expect_equal(by_weight, grid_from_wt_trunc(w, "wt", c(3, 2, 1.5)))
  expect_equal(by_weight$mean[[1]], mean(as.numeric(w)))

  by_count <- wt_trunc_sensitivity(w, method = "count", lower = 1, upper = 1:3)
  expect_equal(by_count, grid_from_wt_trunc(w, "count", 1:3, lower = 1))
})

test_that("a zero weight makes the range ratio infinite", {
  # Weights can be zero. The ratio is then the largest weight over zero, which
  # is Inf, and it stays Inf in every row that leaves the zero in place.
  x <- c(0, 1, 2, 3, 20)
  w <- psw(x, estimand = "ate")

  out <- expect_silent(wt_trunc_sensitivity(w, method = "wt", upper = 10))

  expect_identical(out$min, c(0, 0))
  expect_equal(out$max, c(20, 10))
  expect_identical(out$range_ratio, c(Inf, Inf))
  expect_equal(out$mean, c(26 / 5, 16 / 5))

  # A lower bound above zero moves the zero and gives a finite ratio.
  two_sided <- wt_trunc_sensitivity(w, method = "wt", lower = 0.5, upper = 10)
  expect_equal(two_sided$range_ratio, c(Inf, 20))
})

# ---- inputs refused ----------------------------------------------------------

test_that("propensity scores are refused, as wt_trunc() refuses them", {
  dat <- sensitivity_data(n = 200)
  ps <- stats::fitted(binary_sensitivity_fit(dat))
  scores <- list(
    ps_trim = ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9),
    ps_trunc = ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9),
    ps_calib = ps_calibrate(ps, dat$z)
  )

  for (label in names(scores)) {
    expect_error(
      wt_trunc_sensitivity(scores[[label]]),
      class = "propensity_type_error",
      info = label
    )
  }

  trimmed <- scores$ps_trim
  expect_propensity_error(wt_trunc_sensitivity(trimmed))
})

test_that("input that is not weights is refused", {
  expect_error(
    wt_trunc_sensitivity(c("a", "b")),
    class = "propensity_type_error"
  )
  expect_error(wt_trunc_sensitivity(list(1, 2)), class = "propensity_error")
})

test_that("weights that are already truncated are refused", {
  # A grid over weights that were already bounded would report a reference row
  # that is not the untruncated weights, so it is refused rather than
  # returned. The refusal points back to the original weights.
  once <- wt_trunc(hand_weights(), method = "wt", upper = 5)

  expect_error(
    wt_trunc_sensitivity(once),
    class = "propensity_already_modified_error"
  )
  expect_error(
    wt_trunc_sensitivity(once, method = "wt", upper = c(4, 3)),
    class = "propensity_already_modified_error"
  )

  message <- gsub(
    "\\s+",
    " ",
    conditionMessage(rlang::catch_cnd(wt_trunc_sensitivity(once), "error"))
  )
  expect_match(message, "already", fixed = TRUE)
  expect_match(message, "wt_trunc()", fixed = TRUE)

  expect_propensity_error(wt_trunc_sensitivity(once))
})
