# A trimming or truncation record names units by position. Operations that
# regroup or collapse the units must hand back a record that describes the
# result exactly, or none at all.

# Units 1 and 4 fall outside the bounds.
record_trim_fixture <- function() {
  ps_trim(c(0.05, 0.3, 0.5, 0.95, 0.7, 0.2), lower = 0.1, upper = 0.9)
}

# Units 1 and 4 fall outside the bounds.
record_trunc_fixture <- function() {
  ps_trunc(c(0.05, 0.3, 0.5, 0.95, 0.7, 0.2), lower = 0.1, upper = 0.9)
}

record_frame <- function(scores) {
  df <- data.frame(g = c(1, 2, 1, 2, 1, 2))
  df$ps <- scores
  df
}

# Grouped mutate -------------------------------------------------------------

test_that("a grouped mutate copying a ps_trim column keeps its record silently", {
  skip_if_not_installed("dplyr")
  x <- record_trim_fixture()
  grouped <- dplyr::group_by(record_frame(x), g)

  out <- expect_silent(dplyr::mutate(grouped, copy = ps))

  expect_s3_class(out$copy, "ps_trim")
  expect_identical(ps_trim_meta(out$copy), ps_trim_meta(x))
  expect_identical(is_unit_trimmed(out$copy), is_unit_trimmed(x))
  expect_identical(ps_trim_meta(out$ps), ps_trim_meta(x))
})

test_that("a grouped mutate recombining a ps_trim column drops its positions silently", {
  skip_if_not_installed("dplyr")
  x <- record_trim_fixture()
  grouped <- dplyr::group_by(record_frame(x), g)

  out <- expect_silent(dplyr::mutate(grouped, copy = identity(ps)))

  expect_s3_class(out$copy, "ps_trim")
  expect_identical(vctrs::vec_data(out$copy), vctrs::vec_data(x))
  meta <- ps_trim_meta(out$copy)
  expect_null(meta$keep_idx)
  expect_null(meta$trimmed_idx)
  expect_null(meta$n_obs)
  expect_identical(meta$method, "ps")
  expect_true(is_ps_trimmed(out$copy))
  expect_error(
    is_unit_trimmed(out$copy),
    class = "propensity_missing_meta_error"
  )

  # The column the mutate read is untouched.
  expect_identical(ps_trim_meta(out$ps), ps_trim_meta(x))
})

test_that("a grouped mutate copying a ps_trunc column keeps its record silently", {
  skip_if_not_installed("dplyr")
  x <- record_trunc_fixture()
  grouped <- dplyr::group_by(record_frame(x), g)

  out <- expect_silent(dplyr::mutate(grouped, copy = ps))

  expect_s3_class(out$copy, "ps_trunc")
  expect_identical(ps_trunc_meta(out$copy), ps_trunc_meta(x))
  expect_identical(is_unit_truncated(out$copy), is_unit_truncated(x))
})

test_that("a grouped mutate recombining a ps_trunc column drops its positions silently", {
  skip_if_not_installed("dplyr")
  x <- record_trunc_fixture()
  grouped <- dplyr::group_by(record_frame(x), g)

  out <- expect_silent(dplyr::mutate(grouped, copy = identity(ps)))

  expect_s3_class(out$copy, "ps_trunc")
  expect_identical(vctrs::vec_data(out$copy), vctrs::vec_data(x))
  meta <- ps_trunc_meta(out$copy)
  expect_null(meta$truncated_idx)
  expect_null(meta$n_obs)
  expect_identical(meta$method, "ps")
  expect_true(is_ps_truncated(out$copy))
  expect_error(
    is_unit_truncated(out$copy),
    class = "propensity_missing_meta_error"
  )
})

# unique() -------------------------------------------------------------------

test_that("unique() of a ps_trim re-indexes a record every kept unit agrees with", {
  x <- ps_trim(c(0.05, 0.5, 0.95, 0.5, 0.3), lower = 0.1, upper = 0.9)

  out <- expect_silent(unique(x))

  expect_identical(vctrs::vec_data(out), c(NA, 0.5, 0.3))
  meta <- ps_trim_meta(out)
  expect_identical(meta$trimmed_idx, 1L)
  expect_identical(meta$keep_idx, c(2L, 3L))
  expect_identical(meta$n_obs, 3L)
  expect_identical(is_unit_trimmed(out), c(TRUE, FALSE, FALSE))
})

test_that("unique() of a ps_trim drops a record a missing score shares with a trimmed one", {
  # The one `NA` left stands for a score that arrived missing and for scores
  # the trimming removed, so no single status describes it.
  missing_first <- ps_trim(c(NA, 0.05, 0.5, 0.95), lower = 0.1, upper = 0.9)
  trimmed_first <- ps_trim(c(0.05, NA, 0.5, 0.95), lower = 0.1, upper = 0.9)

  for (x in list(missing_first, trimmed_first)) {
    expect_warning(
      out <- unique(x),
      class = "propensity_trim_record_warning"
    )

    expect_s3_class(out, "ps_trim")
    expect_length(out, 2)
    meta <- ps_trim_meta(out)
    expect_null(meta$keep_idx)
    expect_null(meta$trimmed_idx)
    expect_null(meta$n_obs)
    expect_identical(meta$method, "ps")
    expect_error(
      is_unit_trimmed(out),
      class = "propensity_missing_meta_error"
    )
  }
})

test_that("unique() of a ps_trunc re-indexes a record every kept unit agrees with", {
  x <- ps_trunc(c(0.05, 0.5, 0.95, 0.5, 0.01), lower = 0.1, upper = 0.9)

  out <- expect_silent(unique(x))

  expect_identical(vctrs::vec_data(out), c(0.1, 0.5, 0.9))
  meta <- ps_trunc_meta(out)
  expect_identical(meta$truncated_idx, c(1L, 3L))
  expect_identical(meta$n_obs, 3L)
  expect_identical(is_unit_truncated(out), c(TRUE, FALSE, TRUE))
})

test_that("unique() of a ps_trunc drops a record a bounded score shares with an untouched one", {
  # A score already at the bound and a score moved onto it collapse to one
  # value, and no single status describes it.
  untouched_first <- ps_trunc(c(0.1, 0.05, 0.5), lower = 0.1, upper = 0.9)
  bounded_first <- ps_trunc(c(0.05, 0.1, 0.5), lower = 0.1, upper = 0.9)

  for (x in list(untouched_first, bounded_first)) {
    expect_warning(
      out <- unique(x),
      class = "propensity_trunc_record_warning"
    )

    expect_s3_class(out, "ps_trunc")
    expect_identical(vctrs::vec_data(out), c(0.1, 0.5))
    meta <- ps_trunc_meta(out)
    expect_null(meta$truncated_idx)
    expect_null(meta$n_obs)
    expect_identical(meta$method, "ps")
    expect_error(
      is_unit_truncated(out),
      class = "propensity_missing_meta_error"
    )
  }
})
