# Slicing a modified propensity score without its subscript -------------------

# A `ps_trim` or `ps_trunc` record names units by position. `[`, `sort()`,
# `unique()`, and `c()` of one input place the record themselves. A slice that
# reaches the restore without a subscript, such as `vctrs::vec_slice()` or
# `dplyr::arrange()`, cannot tell a reordering from an identity, so it drops the
# positions at any length and keeps what the record says about the modification.

# Units 1 and 4 fall outside (0.1, 0.9).
slice_scores <- c(0.05, 0.5, 0.6, 0.95, 0.4)

slice_fixtures <- function() {
  list(
    trim = list(
      scores = ps_trim(slice_scores, lower = 0.1, upper = 0.9),
      meta = ps_trim_meta,
      query = is_unit_trimmed,
      flag = is_ps_trimmed
    ),
    trunc = list(
      scores = ps_trunc(slice_scores, lower = 0.1, upper = 0.9),
      meta = ps_trunc_meta,
      query = is_unit_truncated,
      flag = is_ps_truncated
    )
  )
}

expect_score_positions_dropped <- function(out, fixture, label) {
  expect_s3_class(out, class(fixture$scores)[[1]])
  expect_positions_dropped(
    fixture$meta(out),
    fixture$meta(fixture$scores),
    info = label
  )
  expect_true(fixture$flag(out), info = label)
  expect_error(
    fixture$query(out),
    class = "propensity_missing_meta_error",
    info = label
  )
}

test_that("vec_slice() reordering a modified score drops the positions", {
  for (label in names(slice_fixtures())) {
    fixture <- slice_fixtures()[[label]]
    x <- fixture$scores

    out <- expect_silent(vctrs::vec_slice(x, 5:1))
    expect_score_positions_dropped(out, fixture, label)
    expect_identical(vec_data(out), rev(vec_data(x)), info = label)
  }
})

test_that("the reported trim slice no longer names the wrong units", {
  x <- ps_trim(slice_scores, lower = 0.1, upper = 0.9)
  out <- vctrs::vec_slice(x, 5:1)

  expect_identical(which(is.na(out)), c(2L, 5L))
  expect_error(is_unit_trimmed(out), class = "propensity_missing_meta_error")
})

test_that("vec_slice() at the same length in the same order drops the positions", {
  for (label in names(slice_fixtures())) {
    fixture <- slice_fixtures()[[label]]
    x <- fixture$scores

    out <- expect_silent(vctrs::vec_slice(x, seq_along(x)))
    expect_score_positions_dropped(out, fixture, label)
  }
})

test_that("dplyr::arrange() on a modified score column drops the positions", {
  skip_if_not_installed("dplyr")

  for (label in names(slice_fixtures())) {
    fixture <- slice_fixtures()[[label]]
    df <- data.frame(o = 5:1)
    df$ps <- fixture$scores

    out <- expect_silent(dplyr::arrange(df, o))
    expect_score_positions_dropped(out$ps, fixture, label)
  }
})

test_that("weights built from a reordering slice carry no stale positions", {
  exposure <- c(1, 0, 1, 0, 1)

  for (label in names(slice_fixtures())) {
    fixture <- slice_fixtures()[[label]]
    sliced <- vctrs::vec_slice(fixture$scores, 5:1)

    w <- suppressWarnings(
      wt_ate(sliced, .exposure = rev(exposure), exposure_type = "binary"),
      classes = "propensity_no_refit_warning"
    )
    expect_error(
      fixture$query(w),
      class = "propensity_missing_meta_error",
      info = label
    )
  }
})

test_that("vec_assign() and vec_fill_missing() drop the positions", {
  for (label in names(slice_fixtures())) {
    fixture <- slice_fixtures()[[label]]
    x <- fixture$scores

    assigned <- expect_silent(vctrs::vec_assign(x, 2L, 0.55))
    expect_score_positions_dropped(assigned, fixture, label)

    filled <- expect_silent(vctrs::vec_fill_missing(x))
    expect_score_positions_dropped(filled, fixture, label)
  }
})

test_that("subassignment and is.na<- keep the record, since no unit moves", {
  for (label in names(slice_fixtures())) {
    fixture <- slice_fixtures()[[label]]
    x <- fixture$scores
    meta <- fixture$meta(x)
    units <- fixture$query(x)

    expect_silent({
      x[2] <- 0.55
    })
    expect_identical(fixture$meta(x), meta, info = label)
    expect_identical(fixture$query(x), units, info = label)

    expect_silent({
      is.na(x) <- 3
    })
    expect_s3_class(x, class(fixture$scores)[[1]])
    expect_true(is.na(x[[3]]), info = label)
    expect_identical(fixture$meta(x), meta, info = label)
    expect_identical(fixture$query(x), units, info = label)
  }
})

test_that("the routes that place the record still re-index it", {
  for (label in names(slice_fixtures())) {
    fixture <- slice_fixtures()[[label]]
    x <- fixture$scores
    units <- fixture$query(x)

    expect_identical(fixture$query(x[5:1]), rev(units), info = label)
    expect_identical(fixture$query(rev(x)), rev(units), info = label)
    expect_identical(fixture$query(c(x)), units, info = label)
    expect_identical(
      fixture$meta(vctrs::vec_slice(x, integer())),
      fixture$meta(x),
      info = label
    )
  }
})
