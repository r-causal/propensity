# What a weight record keeps when its positions are dropped ------------------

# A trimming, truncation, or weight truncation record says two kinds of thing:
# how the modification was made (its method, bounds, and refit flag) and which
# units it touched. A route that cannot place the units drops only the second.
# The first still describes the weights, so a later combine compares bounds
# against it, `is_refit()` reads it, and the footer prints it.

count_truncated_psw <- function(upper) {
  wt_trunc(psw(c(1, 2, 50, 3, 4, 5, 6, 7, 8, 100)), "count", upper = upper)
}

# Units 1 and 4 fall outside the bounds; their weights are missing.
positions_trimmed_psw <- function() {
  w <- psw(c(NA, 2, 1.6, NA, 3.3), estimand = "ate; trimmed", trimmed = TRUE)
  attr(w, "ps_trim_meta") <- list(
    method = "ps",
    lower = 0.1,
    upper = 0.9,
    focal_inverted = FALSE,
    keep_idx = c(2L, 3L, 5L),
    trimmed_idx = c(1L, 4L),
    n_obs = 5L,
    refit = TRUE
  )
  w
}

positions_truncated_psw <- function() {
  w <- psw(
    c(10, 2, 1.6, 10, 3.3),
    estimand = "ate; truncated",
    truncated = TRUE
  )
  attr(w, "ps_trunc_meta") <- list(
    method = "ps",
    lower_bound = 0.1,
    upper_bound = 0.9,
    truncated_idx = c(1L, 4L),
    n_obs = 5L
  )
  w
}

positions_fixtures <- function() {
  list(
    trimmed = list(
      weights = positions_trimmed_psw(),
      attr = "ps_trim_meta",
      query = is_unit_trimmed
    ),
    truncated = list(
      weights = positions_truncated_psw(),
      attr = "ps_trunc_meta",
      query = is_unit_truncated
    ),
    wt_truncated = list(
      weights = count_truncated_psw(2),
      attr = "psw_trunc_meta",
      query = is_unit_wt_truncated
    )
  )
}

expect_fixture_positions_dropped <- function(out, fixture, label) {
  expect_s3_class(out, "psw")
  expect_positions_dropped(
    attr(out, fixture$attr),
    attr(fixture$weights, fixture$attr),
    info = label
  )
  expect_error(
    fixture$query(out),
    class = "propensity_missing_meta_error",
    info = label
  )
}

footer_lines <- function(x) {
  utils::capture.output(print(x))
}

test_that("a combine still compares bounds after vec_slice() reorders one input", {
  a <- count_truncated_psw(2)
  b <- count_truncated_psw(1)

  expect_warning(
    out <- vctrs::vec_c(vctrs::vec_slice(a, 10:1), b),
    class = "propensity_coercion_warning"
  )
  expect_false(is_psw(out))
  expect_type(out, "double")

  expect_warning(
    out <- vctrs::vec_c(b, vctrs::vec_slice(a, 1:5)),
    class = "propensity_coercion_warning"
  )
  expect_false(is_psw(out))
})

test_that("a combine still compares bounds after dplyr::arrange()", {
  skip_if_not_installed("dplyr")
  da <- data.frame(o = 10:1)
  da$w <- count_truncated_psw(2)
  db <- data.frame(o = 1:10)
  db$w <- count_truncated_psw(1)

  expect_warning(
    out <- dplyr::bind_rows(dplyr::arrange(da, o), db),
    class = "propensity_coercion_warning"
  )
  expect_false(is_psw(out$w))
})

test_that("weights bounded alike still combine after a reordering slice", {
  a <- count_truncated_psw(2)

  out <- expect_silent(vctrs::vec_c(vctrs::vec_slice(a, 10:1), a))
  expect_s3_class(out, "psw")
  expect_true(is_wt_truncated(out))
})

test_that("is_refit() answers for weights reordered by dplyr::arrange()", {
  skip_if_not_installed("dplyr")
  w <- positions_trimmed_psw()
  expect_true(is_refit(w))

  d <- data.frame(X = c(3, 1, 2, 5, 4))
  d$w <- w
  out <- expect_silent(dplyr::arrange(d, X))

  expect_true(is_refit(out$w))
  expect_true(is_refit(vctrs::vec_slice(w, 5:1)))
  expect_error(
    is_unit_trimmed(out$w),
    class = "propensity_missing_meta_error"
  )
})

test_that("is_refit() answers for weights from a refit model after a reordering", {
  skip_if_not_installed("dplyr")
  set.seed(42)
  n <- 60
  z <- rnorm(n)
  x <- rbinom(n, 1, plogis(2 * z))
  fit <- glm(x ~ z, family = binomial())
  trimmed <- ps_trim(
    predict(fit, type = "response"),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  wt <- wt_ate(
    ps_refit(trimmed, model = fit),
    .exposure = x,
    exposure_type = "binary",
    .focal_level = 1
  )
  d <- data.frame(X = z)
  d$w <- wt

  expect_true(is_refit(dplyr::arrange(d, X)$w))
})

test_that("the footer prints the bound of a record whose positions were dropped", {
  skip_if_not_installed("dplyr")
  w <- count_truncated_psw(2)
  expect_true(any(
    footer_lines(w) ==
      "truncation: count 2 (upper 8), 2 of 10 weights truncated"
  ))

  d <- data.frame(o = 10:1)
  d$w <- w
  arranged <- dplyr::arrange(d, o)$w
  expect_true(any(footer_lines(arranged) == "truncation: count 2 (upper 8)"))

  sliced <- vctrs::vec_slice(w, 10:1)
  expect_true(any(footer_lines(sliced) == "truncation: count 2 (upper 8)"))
})

test_that("the footer of flagged weights with no record names no bound", {
  w <- psw(c(1, 2, 3), estimand = "ate", wt_truncated = TRUE)

  expect_true(any(footer_lines(w) == "truncation: weights truncated"))
})

# Assignment -------------------------------------------------------------------

test_that("is.na<- keeps each record, since every unit stays in place", {
  for (label in names(positions_fixtures())) {
    fixture <- positions_fixtures()[[label]]
    w <- fixture$weights
    meta <- attr(w, fixture$attr)
    units <- fixture$query(w)

    expect_silent({
      is.na(w) <- 2
    })
    expect_s3_class(w, "psw")
    expect_true(is.na(w[[2]]), info = label)
    expect_identical(attr(w, fixture$attr), meta, info = label)
    expect_identical(fixture$query(w), units, info = label)

    expect_silent({
      is.na(w) <- c(FALSE, FALSE, TRUE, rep(FALSE, length(w) - 3))
    })
    expect_true(is.na(w[[3]]), info = label)
    expect_identical(attr(w, fixture$attr), meta, info = label)
  }
})

test_that("is.na<- keeps the other attributes of the weights", {
  w <- positions_trimmed_psw()
  before <- attributes(w)

  is.na(w) <- 3

  expect_identical(attributes(w)[names(before)], before)
  expect_identical(vec_data(w), c(NA, 2, NA, NA, 3.3))
})

test_that("vec_assign() drops each record's positions, since it cannot vouch for them", {
  for (label in names(positions_fixtures())) {
    fixture <- positions_fixtures()[[label]]
    w <- fixture$weights

    out <- expect_silent(vctrs::vec_assign(w, 1L, 3))
    expect_fixture_positions_dropped(out, fixture, label)
    expect_identical(vec_data(out)[[1]], 3, info = label)
  }
})

test_that("vec_assign() based helpers drop the record's positions", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  fixture <- positions_fixtures()$trimmed
  w <- fixture$weights

  replaced <- expect_silent(tidyr::replace_na(w, 0))
  expect_fixture_positions_dropped(replaced, fixture, "replace_na")
  expect_identical(vec_data(replaced), c(0, 2, 1.6, 0, 3.3))

  # These first combine their inputs into a common type, and a combined type
  # carries no positional record at all.
  combined <- list(
    coalesce = expect_silent(dplyr::coalesce(w, w * 0)),
    if_else = expect_silent(dplyr::if_else(is.na(w), w, w * 2)),
    case_when = expect_silent(dplyr::case_when(is.na(w) ~ w, .default = w))
  )
  for (label in names(combined)) {
    out <- combined[[label]]
    expect_s3_class(out, "psw")
    expect_true(is_ps_trimmed(out), info = label)
    expect_null(attr(out, "ps_trim_meta")$trimmed_idx, info = label)
    expect_error(
      is_unit_trimmed(out),
      class = "propensity_missing_meta_error",
      info = label
    )
  }
})

test_that("rep_len() and vec_c() of one input drop each record's positions", {
  for (label in names(positions_fixtures())) {
    fixture <- positions_fixtures()[[label]]
    w <- fixture$weights

    out <- expect_silent(rep_len(w, length(w)))
    expect_fixture_positions_dropped(out, fixture, label)

    out <- expect_silent(vctrs::vec_c(w))
    expect_fixture_positions_dropped(out, fixture, label)
  }
})

test_that("casting other data onto full-length weights drops each record's positions", {
  for (label in names(positions_fixtures())) {
    fixture <- positions_fixtures()[[label]]
    w <- fixture$weights

    out <- expect_silent(vctrs::vec_cast(rev(vec_data(w)), w))
    expect_fixture_positions_dropped(out, fixture, label)
  }
})
