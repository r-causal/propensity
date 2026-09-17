# Combining a single set of weights with base `c()` returns that set, so every
# record naming its observations still describes them. Every other combine goes
# through vctrs, which cannot tell the input's own prototype from one a caller
# supplied, and drops the positional records.

# fmt: skip
single_combine_ps <- c(
  0.05, 0.15, 0.3, 0.45, 0.6, 0.72, 0.85, 0.95, 0.4, 0.55,
  0.08, 0.22, 0.35, 0.5, 0.65, 0.78, 0.9, 0.12, 0.62, 0.33
)
single_combine_z <- rep(c(1, 0, 1, 0, 1), times = 4)

single_trimmed_psw <- function() {
  trimmed <- ps_trim(single_combine_ps, lower = 0.1, upper = 0.9)
  withr::local_options(propensity.quiet = TRUE)
  suppressWarnings(
    wt_ate(trimmed, single_combine_z),
    classes = "propensity_no_refit_warning"
  )
}

single_truncated_psw <- function() {
  withr::local_options(propensity.quiet = TRUE)
  wt_ate(
    ps_trunc(single_combine_ps, lower = 0.1, upper = 0.9),
    single_combine_z
  )
}

single_wt_truncated_psw <- function() {
  withr::local_options(propensity.quiet = TRUE)
  wt_trunc(
    wt_ate(single_combine_ps, single_combine_z),
    method = "pctl",
    upper = 0.9
  )
}

single_calibrated_psw <- function() {
  withr::local_options(propensity.quiet = TRUE)
  calibrated <- ps_calibrate(
    single_combine_ps,
    single_combine_z,
    smooth = FALSE
  )
  wt_ate(calibrated, single_combine_z)
}

# Weights built from a trimmed score -----------------------------------------

test_that("c() of one trimmed psw keeps its trimming record", {
  w <- single_trimmed_psw()
  expect_false(is.null(ps_trim_meta(w)$trimmed_idx))

  combined <- expect_silent(c(w))

  expect_identical(ps_trim_meta(combined), ps_trim_meta(w))
  expect_identical(is_unit_trimmed(combined), is_unit_trimmed(w))
  expect_identical(combined, w)
})

test_that("c() of one psw keeps its names, and a NULL alongside it", {
  w <- single_trimmed_psw()
  named <- w
  names(named) <- paste0("u", seq_along(w))

  expect_identical(c(named), named)

  # Base R discards a NULL argument before dispatch, so this is a combine of
  # one set of weights.
  expect_identical(c(w, NULL), w)
})

test_that("c() of one psw whose record no longer covers it returns it unchanged", {
  # A model frame shortens a weights column and re-attaches the attributes, so
  # the record describes more observations than the weights hold.
  w <- single_trimmed_psw()
  stale <- w
  attributes(stale) <- NULL
  stale <- stale[1:15]
  attributes(stale) <- attributes(w)

  combined <- expect_silent(c(stale))

  expect_identical(combined, stale)
  expect_identical(ps_trim_meta(combined), ps_trim_meta(w))
  expect_error(
    is_unit_trimmed(combined),
    class = "propensity_missing_meta_error"
  )
})

test_that("c() of one psw still refuses what vctrs refuses", {
  w <- single_trimmed_psw()

  expect_error(c(w, recursive = TRUE), "recursive")
  expect_error(c(w, use.names = FALSE), "use.names")
  expect_error(c(x = w), "outer name")
})

test_that("vec_c() and list_unchop() of one trimmed psw drop the record", {
  # vctrs restores the result against a zero-length prototype, which cannot be
  # told apart from one a caller supplied for data from elsewhere.
  w <- single_trimmed_psw()

  from_vec_c <- vctrs::vec_c(w)
  expect_null(attr(from_vec_c, "ps_trim_meta"))
  expect_true(is_ps_trimmed(from_vec_c))
  expect_error(
    is_unit_trimmed(from_vec_c),
    class = "propensity_missing_meta_error"
  )

  from_unchop <- vctrs::list_unchop(list(w))
  expect_null(attr(from_unchop, "ps_trim_meta"))
  expect_error(
    is_unit_trimmed(from_unchop),
    class = "propensity_missing_meta_error"
  )
})

test_that("c() of one truncated-score psw keeps its truncation record", {
  w <- single_truncated_psw()
  expect_false(is.null(ps_trunc_meta(w)$truncated_idx))

  combined <- expect_silent(c(w))

  expect_identical(ps_trunc_meta(combined), ps_trunc_meta(w))
  expect_identical(is_unit_truncated(combined), is_unit_truncated(w))
  expect_identical(combined, w)
})

test_that("c() of one weight-truncated psw keeps the weight truncation record", {
  w <- single_wt_truncated_psw()
  expect_true(any(is_unit_wt_truncated(w)))

  combined <- expect_silent(c(w))

  expect_identical(
    attr(combined, "psw_trunc_meta"),
    attr(w, "psw_trunc_meta")
  )
  expect_identical(is_unit_wt_truncated(combined), is_unit_wt_truncated(w))
  expect_identical(combined, w)
})

test_that("c() of one calibrated-score psw keeps its calibration record", {
  w <- single_calibrated_psw()
  expect_false(is.null(ps_calib_meta(w)))

  combined <- expect_silent(c(w))

  expect_identical(ps_calib_meta(combined), ps_calib_meta(w))
  expect_identical(combined, w)
})

# Two or more inputs ---------------------------------------------------------

test_that("c() of two trimmed psw drops the positional records", {
  w <- single_trimmed_psw()

  combined <- c(w, w)

  expect_length(combined, 2 * length(w))
  expect_null(attr(combined, "ps_trim_meta"))
  expect_true(is_ps_trimmed(combined))
})

test_that("c() of a psw split in two drops the records at the original length", {
  # The pieces add back up to the length the records were written for, but the
  # result is appended from two inputs, which nothing re-indexes.
  w <- single_trimmed_psw()
  first <- suppressWarnings(w[1:10])
  second <- suppressWarnings(w[11:20])

  combined <- c(first, second)

  expect_length(combined, length(w))
  expect_null(attr(combined, "ps_trim_meta"))
})

test_that("c() of two weight-truncated psw drops the weight truncation record", {
  w <- single_wt_truncated_psw()

  combined <- c(w, w)

  expect_null(attr(combined, "psw_trunc_meta"))
  expect_true(is_wt_truncated(combined))
})

test_that("c() of two truncated-score psw drops the truncation record", {
  w <- single_truncated_psw()

  combined <- c(w, w)

  expect_null(attr(combined, "ps_trunc_meta"))
})

# A prototype or cast target supplied by the caller -------------------------

# Each of these restores data that is not the weights' own, in their order,
# against a zero-length prototype carrying the weights' records. None of the
# records describes the result.

expect_trim_record_dropped <- function(x) {
  expect_null(attr(x, "ps_trim_meta"))
  expect_error(is_unit_trimmed(x), class = "propensity_missing_meta_error")
}

test_that("reordered pieces combined at a supplied prototype drop the record", {
  w <- single_trimmed_psw()
  first <- suppressWarnings(w[1:10])
  second <- suppressWarnings(w[11:20])

  expect_trim_record_dropped(
    vctrs::list_unchop(list(second, first), ptype = vctrs::vec_ptype(w))
  )
  expect_trim_record_dropped(vctrs::vec_c(second, first, .ptype = w))
})

test_that("other weights combined at a supplied prototype drop the record", {
  w <- single_trimmed_psw()
  other <- suppressWarnings(rev(w))

  expect_trim_record_dropped(vctrs::vec_c(other, .ptype = w))
})

test_that("a cast to or an initialization at the prototype drops the record", {
  w <- single_trimmed_psw()

  expect_trim_record_dropped(
    vctrs::vec_cast(seq(1, 2, length.out = 20), vctrs::vec_ptype(w))
  )
  expect_trim_record_dropped(vctrs::vec_init(w[0], 20))
})

test_that("a supplied prototype drops the weight truncation record", {
  w <- single_wt_truncated_psw()
  first <- suppressWarnings(w[1:10])
  second <- suppressWarnings(w[11:20])

  unchopped <- vctrs::list_unchop(
    list(second, first),
    ptype = vctrs::vec_ptype(w)
  )
  expect_null(attr(unchopped, "psw_trunc_meta"))
  expect_error(
    is_unit_wt_truncated(unchopped),
    class = "propensity_missing_meta_error"
  )

  cast <- vctrs::vec_cast(seq(1, 2, length.out = 20), vctrs::vec_ptype(w))
  expect_null(attr(cast, "psw_trunc_meta"))
  expect_error(
    is_unit_wt_truncated(cast),
    class = "propensity_missing_meta_error"
  )
})

# The propensity score classes -----------------------------------------------

test_that("c() of one ps_trim or ps_trunc keeps its record", {
  trimmed <- ps_trim(single_combine_ps, lower = 0.1, upper = 0.9)
  expect_identical(ps_trim_meta(c(trimmed)), ps_trim_meta(trimmed))
  expect_identical(is_unit_trimmed(c(trimmed)), is_unit_trimmed(trimmed))

  truncated <- ps_trunc(single_combine_ps, lower = 0.1, upper = 0.9)
  expect_identical(ps_trunc_meta(c(truncated)), ps_trunc_meta(truncated))
  expect_identical(
    is_unit_truncated(c(truncated)),
    is_unit_truncated(truncated)
  )
})

test_that("c() of two ps_trim or ps_trunc drops the positions", {
  trimmed <- ps_trim(single_combine_ps, lower = 0.1, upper = 0.9)
  expect_null(ps_trim_meta(c(trimmed, trimmed))$trimmed_idx)

  truncated <- ps_trunc(single_combine_ps, lower = 0.1, upper = 0.9)
  expect_null(ps_trunc_meta(c(truncated, truncated))$truncated_idx)
})
