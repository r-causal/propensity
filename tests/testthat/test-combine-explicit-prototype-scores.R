# A zero-length prototype carries the record of the scores it was taken from,
# and vctrs restores a combined or initialized result against it whether the
# prototype was sliced from the one input or supplied by the caller for data
# from elsewhere. The positions in that record describe none of the result, so
# the result keeps its class and the description of the modification, and the
# positional queries refuse. Base `c()` of a single score vector is the one
# combine that knows it was handed that vector alone, and it keeps the record.

# fmt: skip
prototype_ps <- c(
  0.05, 0.15, 0.3, 0.45, 0.6, 0.72, 0.85, 0.95, 0.4, 0.55,
  0.08, 0.22, 0.35, 0.5, 0.65, 0.78, 0.92, 0.12, 0.62, 0.33
)

prototype_trim <- function() {
  ps_trim(prototype_ps, lower = 0.1, upper = 0.9)
}

prototype_trunc <- function() {
  ps_trunc(prototype_ps, lower = 0.1, upper = 0.9)
}

expect_trim_positions_dropped <- function(x, like) {
  expect_s3_class(x, "ps_trim")
  meta <- ps_trim_meta(x)
  expect_null(meta$keep_idx)
  expect_null(meta$trimmed_idx)
  expect_null(meta$n_obs)
  expect_identical(meta, drop_trim_record(ps_trim_meta(like)))
  expect_error(is_unit_trimmed(x), class = "propensity_missing_meta_error")
}

expect_trunc_positions_dropped <- function(x, like) {
  expect_s3_class(x, "ps_trunc")
  meta <- ps_trunc_meta(x)
  expect_null(meta$truncated_idx)
  expect_null(meta$n_obs)
  expect_identical(meta, drop_trunc_record(ps_trunc_meta(like)))
  expect_error(is_unit_truncated(x), class = "propensity_missing_meta_error")
}

# The trimmed units are 1, 8, 11, and 17. Reordered as units 11 to 20 then 1
# to 10, they sit at 1, 7, 11, and 18, which no record written for the input
# order names.
reordered_pieces <- function(x) {
  list(x[11:20], x[1:10])
}

# ps_trim ---------------------------------------------------------------------

test_that("reordered trimmed pieces at a supplied prototype drop the positions", {
  t <- prototype_trim()
  pieces <- reordered_pieces(t)

  unchopped <- expect_silent(
    vctrs::list_unchop(pieces, ptype = vctrs::vec_ptype(t))
  )
  expect_identical(which(is.na(unchopped)), c(1L, 7L, 11L, 18L))
  expect_trim_positions_dropped(unchopped, t)

  combined <- expect_silent(vctrs::vec_c(!!!pieces, .ptype = t))
  expect_identical(vctrs::vec_data(combined), vctrs::vec_data(unchopped))
  expect_trim_positions_dropped(combined, t)
})

test_that("other trimmed scores at a supplied prototype drop the positions", {
  t <- prototype_trim()
  other <- rev(t)

  expect_trim_positions_dropped(vctrs::vec_c(other, .ptype = t), t)
})

test_that("vec_c() and list_unchop() of one ps_trim drop the positions", {
  t <- prototype_trim()

  expect_trim_positions_dropped(vctrs::vec_c(t), t)
  expect_trim_positions_dropped(vctrs::list_unchop(list(t)), t)
})

test_that("a double cast or initialized at a trimmed prototype names no unit", {
  t <- prototype_trim()

  # A cast writes a record for the values it was handed, none of them trimmed.
  cast <- vctrs::vec_cast(seq(0.2, 0.8, length.out = 20), vctrs::vec_ptype(t))
  expect_identical(ps_trim_meta(cast)$trimmed_idx, integer(0))
  expect_identical(ps_trim_meta(cast)$keep_idx, 1:20)
  expect_identical(ps_trim_meta(cast)$n_obs, 20L)
  expect_identical(is_unit_trimmed(cast), rep(FALSE, 20))

  expect_trim_positions_dropped(vctrs::vec_init(t[0], 20), t)
})

test_that("c() of one ps_trim returns it unchanged", {
  t <- prototype_trim()
  named <- t
  names(named) <- paste0("u", seq_along(t))

  expect_identical(expect_silent(c(t)), t)
  expect_identical(c(named), named)
  expect_identical(c(t, NULL), t)
  expect_identical(is_unit_trimmed(c(t)), is_unit_trimmed(t))

  expect_error(c(t, recursive = TRUE), "recursive")
  expect_error(c(t, use.names = FALSE), "use.names")
})

test_that("c() of two ps_trim still drops the positions", {
  t <- prototype_trim()

  combined <- c(t, t)

  expect_length(combined, 40)
  expect_trim_positions_dropped(combined, t)
})

test_that("slicing, subassignment, and arithmetic on a ps_trim are unchanged", {
  t <- prototype_trim()

  sliced <- t[c(1, 2, 8)]
  expect_identical(ps_trim_meta(sliced)$trimmed_idx, c(1L, 3L))
  expect_identical(ps_trim_meta(sliced)$n_obs, 3L)

  assigned <- t
  assigned[2] <- 0.5
  expect_identical(ps_trim_meta(assigned), ps_trim_meta(t))
  expect_identical(is_unit_trimmed(assigned), is_unit_trimmed(t))

  expect_identical(t + 0, vctrs::vec_data(t))
})

# ps_trunc --------------------------------------------------------------------

test_that("reordered truncated pieces at a supplied prototype drop the positions", {
  tr <- prototype_trunc()
  pieces <- reordered_pieces(tr)

  unchopped <- expect_silent(
    vctrs::list_unchop(pieces, ptype = vctrs::vec_ptype(tr))
  )
  expect_trunc_positions_dropped(unchopped, tr)

  combined <- expect_silent(vctrs::vec_c(!!!pieces, .ptype = tr))
  expect_trunc_positions_dropped(combined, tr)
})

test_that("other truncated scores at a supplied prototype drop the positions", {
  tr <- prototype_trunc()
  other <- rev(tr)

  expect_trunc_positions_dropped(vctrs::vec_c(other, .ptype = tr), tr)
})

test_that("vec_c() and list_unchop() of one ps_trunc drop the positions", {
  tr <- prototype_trunc()

  expect_trunc_positions_dropped(vctrs::vec_c(tr), tr)
  expect_trunc_positions_dropped(vctrs::list_unchop(list(tr)), tr)
})

test_that("a double cast or initialized at a truncated prototype names no unit", {
  tr <- prototype_trunc()

  cast <- vctrs::vec_cast(seq(0.2, 0.8, length.out = 20), vctrs::vec_ptype(tr))
  expect_identical(ps_trunc_meta(cast)$truncated_idx, integer(0))
  expect_identical(ps_trunc_meta(cast)$n_obs, 20L)
  expect_identical(is_unit_truncated(cast), rep(FALSE, 20))

  expect_trunc_positions_dropped(vctrs::vec_init(tr[0], 20), tr)
})

test_that("c() of one ps_trunc returns it unchanged", {
  tr <- prototype_trunc()
  named <- tr
  names(named) <- paste0("u", seq_along(tr))

  expect_identical(expect_silent(c(tr)), tr)
  expect_identical(c(named), named)
  expect_identical(c(tr, NULL), tr)
  expect_identical(is_unit_truncated(c(tr)), is_unit_truncated(tr))

  expect_error(c(tr, recursive = TRUE), "recursive")
  expect_error(c(tr, use.names = FALSE), "use.names")
})

test_that("c() of two ps_trunc still drops the positions", {
  tr <- prototype_trunc()

  expect_trunc_positions_dropped(c(tr, tr), tr)
})

test_that("slicing, subassignment, and arithmetic on a ps_trunc are unchanged", {
  tr <- prototype_trunc()

  sliced <- tr[c(1, 2, 8)]
  expect_identical(ps_trunc_meta(sliced)$truncated_idx, c(1L, 3L))
  expect_identical(ps_trunc_meta(sliced)$n_obs, 3L)

  assigned <- tr
  assigned[2] <- 0.5
  expect_identical(ps_trunc_meta(assigned), ps_trunc_meta(tr))

  expect_type(tr + 0, "double")
})
