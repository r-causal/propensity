# `unique()` of a matrix of scores keeps one row for each distinct row. The
# result is a matrix of the same class, and its record is re-indexed onto the
# kept rows only when every set of merged rows shares one status.

# fmt: skip
unique_matrix_scores <- matrix(
  c(
    0.8,  0.1,  0.1,
    0.4,  0.3,  0.3,
    0.02, 0.49, 0.49,
    0.4,  0.3,  0.3,
    0.98, 0.01, 0.01,
    0.2,  0.4,  0.4,
    0.3,  0.3,  0.4,
    NA,   NA,   NA
  ),
  ncol = 3,
  byrow = TRUE,
  dimnames = list(NULL, c("a", "b", "c"))
)
unique_matrix_exposure <- factor(c("a", "b", "c", "a", "b", "c", "a", "b"))

# Rows 3 and 5 fall below the threshold. Row 4 repeats row 2, and row 8 arrived
# missing.
unique_trim_matrix <- function(rows = 1:8) {
  ps_trim(
    unique_matrix_scores[rows, , drop = FALSE],
    .exposure = unique_matrix_exposure[rows],
    method = "ps",
    lower = 0.05
  )
}

unique_trunc_matrix <- function(rows = 1:8) {
  ps_trunc(
    unique_matrix_scores[rows, , drop = FALSE],
    .exposure = unique_matrix_exposure[rows],
    method = "ps",
    lower = 0.05
  )
}

# ps_trim_matrix --------------------------------------------------------------

test_that("unique() of a ps_trim_matrix with distinct rows returns it whole", {
  x <- unique_trim_matrix(c(1, 2, 3, 6, 7))

  out <- expect_silent(unique(x))

  expect_s3_class(out, "ps_trim_matrix")
  expect_identical(dim(out), c(5L, 3L))
  expect_identical(dimnames(out), dimnames(x))
  expect_identical(unclass(out)[,], unclass(x)[,])
  expect_identical(ps_trim_meta(out), ps_trim_meta(x))
  expect_identical(is_unit_trimmed(out), is_unit_trimmed(x))
})

test_that("unique() of a ps_trim_matrix re-indexes a record the merged rows agree on", {
  # Rows 2 and 4 are both retained, and rows 3 and 5 are both trimmed.
  x <- unique_trim_matrix(1:7)

  out <- expect_silent(unique(x))

  expect_s3_class(out, "ps_trim_matrix")
  expect_identical(dim(out), c(5L, 3L))
  expect_identical(colnames(out), c("a", "b", "c"))
  expect_identical(
    unclass(out)[,],
    unclass(x)[c(1, 2, 3, 6, 7), ]
  )

  meta <- ps_trim_meta(out)
  expect_identical(meta$trimmed_idx, 3L)
  expect_identical(meta$keep_idx, c(1L, 2L, 4L, 5L))
  expect_identical(meta$n_obs, 5L)
  expect_identical(meta$method, "ps")
  expect_identical(is_unit_trimmed(out), c(FALSE, FALSE, TRUE, FALSE, FALSE))
})

test_that("unique() of a ps_trim_matrix drops a record a missing row shares with a trimmed one", {
  # A trimmed row is all `NA`, so it merges with a row that arrived missing,
  # and neither status describes the row left.
  trimmed_first <- unique_trim_matrix(1:8)
  missing_first <- unique_trim_matrix(c(8, 1:7))

  for (x in list(trimmed_first, missing_first)) {
    out <- expect_silent(unique(x))

    expect_s3_class(out, "ps_trim_matrix")
    expect_identical(dim(out), c(5L, 3L))
    expect_identical(colnames(out), c("a", "b", "c"))

    meta <- ps_trim_meta(out)
    expect_null(meta$keep_idx)
    expect_null(meta$trimmed_idx)
    expect_null(meta$n_obs)
    expect_identical(meta$method, "ps")
    expect_true(is_ps_trimmed(out))
    expect_error(
      is_unit_trimmed(out),
      class = "propensity_missing_meta_error"
    )
  }
})

# ps_trunc_matrix -------------------------------------------------------------

test_that("unique() of a ps_trunc_matrix with distinct rows returns it whole", {
  x <- unique_trunc_matrix(c(1, 2, 3, 6, 7))

  out <- expect_silent(unique(x))

  expect_s3_class(out, "ps_trunc_matrix")
  expect_identical(dim(out), c(5L, 3L))
  expect_identical(dimnames(out), dimnames(x))
  expect_identical(unclass(out)[,], unclass(x)[,])
  expect_identical(ps_trunc_meta(out), ps_trunc_meta(x))
  expect_identical(is_unit_truncated(out), is_unit_truncated(x))
})

test_that("unique() of a ps_trunc_matrix re-indexes a record the merged rows agree on", {
  # Rows 2 and 4 are identical and both unmoved. Row 8 arrived missing and
  # merges with nothing.
  x <- unique_trunc_matrix(1:8)

  out <- expect_silent(unique(x))

  expect_s3_class(out, "ps_trunc_matrix")
  expect_identical(dim(out), c(7L, 3L))
  expect_identical(colnames(out), c("a", "b", "c"))
  expect_identical(
    unclass(out)[,],
    unclass(x)[c(1, 2, 3, 5, 6, 7, 8), ]
  )

  meta <- ps_trunc_meta(out)
  expect_identical(meta$truncated_idx, c(3L, 4L))
  expect_identical(meta$n_obs, 7L)
  expect_identical(meta$method, "ps")
  expect_identical(
    is_unit_truncated(out),
    c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
  )
})

test_that("unique() of a ps_trunc_matrix drops a record a moved row shares with an unmoved one", {
  # The record says row 1 was moved and row 2, which holds the same scores, was
  # not, so the one row left has no single status.
  template <- unique_trunc_matrix(1:3)
  # fmt: skip
  scores <- matrix(
    c(
      0.5, 0.25, 0.25,
      0.5, 0.25, 0.25,
      0.2, 0.4,  0.4
    ),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )
  meta <- ps_trunc_meta(template)
  meta$truncated_idx <- 1L
  meta$n_obs <- 3L
  moved_first <- new_ps_trunc(scores, meta)

  meta$truncated_idx <- 2L
  unmoved_first <- new_ps_trunc(scores, meta)

  for (x in list(moved_first, unmoved_first)) {
    out <- expect_silent(unique(x))

    expect_s3_class(out, "ps_trunc_matrix")
    expect_identical(dim(out), c(2L, 3L))
    expect_identical(colnames(out), c("a", "b", "c"))
    expect_identical(unclass(out)[,], scores[c(1, 3), ])

    out_meta <- ps_trunc_meta(out)
    expect_null(out_meta$truncated_idx)
    expect_null(out_meta$n_obs)
    expect_identical(out_meta$method, "ps")
    expect_true(is_ps_truncated(out))
    expect_error(
      is_unit_truncated(out),
      class = "propensity_missing_meta_error"
    )
  }
})

# One row and no rows ----------------------------------------------------------

test_that("unique() of a score matrix whose rows are all identical returns one row", {
  # Rows 2 and 4 hold the same scores, and both were retained and left unmoved.
  trimmed <- unique_trim_matrix(1:7)[c(2, 4), , drop = FALSE]
  truncated <- unique_trunc_matrix(1:7)[c(2, 4), , drop = FALSE]

  out_trim <- expect_silent(unique(trimmed))

  expect_s3_class(out_trim, "ps_trim_matrix")
  expect_identical(dim(out_trim), c(1L, 3L))
  expect_identical(colnames(out_trim), c("a", "b", "c"))
  expect_identical(
    unclass(out_trim)[,],
    unclass(trimmed)[1, , drop = FALSE][,]
  )
  trim_meta <- ps_trim_meta(out_trim)
  expect_identical(trim_meta$n_obs, 1L)
  expect_identical(trim_meta$keep_idx, 1L)
  expect_identical(trim_meta$trimmed_idx, integer(0))
  expect_true(is_ps_trimmed(out_trim))
  expect_identical(is_unit_trimmed(out_trim), FALSE)

  out_trunc <- expect_silent(unique(truncated))

  expect_s3_class(out_trunc, "ps_trunc_matrix")
  expect_identical(dim(out_trunc), c(1L, 3L))
  expect_identical(colnames(out_trunc), c("a", "b", "c"))
  expect_identical(
    unclass(out_trunc)[,],
    unclass(truncated)[1, , drop = FALSE][,]
  )
  trunc_meta <- ps_trunc_meta(out_trunc)
  expect_identical(trunc_meta$n_obs, 1L)
  expect_identical(trunc_meta$truncated_idx, integer(0))
  expect_true(is_ps_truncated(out_trunc))
  expect_identical(is_unit_truncated(out_trunc), FALSE)
})

test_that("unique() of a score matrix with no rows returns no rows", {
  trimmed <- unique_trim_matrix(1:7)[integer(0), , drop = FALSE]
  truncated <- unique_trunc_matrix(1:7)[integer(0), , drop = FALSE]

  out_trim <- expect_silent(unique(trimmed))

  expect_s3_class(out_trim, "ps_trim_matrix")
  expect_identical(dim(out_trim), c(0L, 3L))
  expect_identical(colnames(out_trim), c("a", "b", "c"))
  expect_identical(ps_trim_meta(out_trim)$n_obs, 0L)
  expect_identical(is_unit_trimmed(out_trim), logical(0))

  out_trunc <- expect_silent(unique(truncated))

  expect_s3_class(out_trunc, "ps_trunc_matrix")
  expect_identical(dim(out_trunc), c(0L, 3L))
  expect_identical(colnames(out_trunc), c("a", "b", "c"))
  expect_identical(ps_trunc_meta(out_trunc)$n_obs, 0L)
  expect_identical(is_unit_truncated(out_trunc), logical(0))
})
