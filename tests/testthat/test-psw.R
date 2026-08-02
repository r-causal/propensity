library(testthat)
library(vctrs)

test_that("new_psw creates valid psw objects", {
  x <- new_psw(c(0.1, 0.2, 0.3), estimand = "ate")
  expect_s3_class(
    x,
    c("psw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vec_data(x), c(0.1, 0.2, 0.3))
  expect_equal(estimand(x), "ate")
})

test_that("psw helper function works correctly", {
  x <- psw(c(0.1, 0.2, 0.3), estimand = "att")
  expect_s3_class(
    x,
    c("psw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vec_data(x), c(0.1, 0.2, 0.3))
  expect_equal(estimand(x), "att")
  expect_false(is_stabilized(x))
  expect_false(is_ps_trimmed(x))
  expect_false(is_ps_truncated(x))
  estimand(x) <- "ATT!"
  expect_equal(estimand(x), "ATT!")
})

test_that("is_psw identifies psw objects", {
  x <- psw(c(0.1, 0.2), estimand = "ato")
  expect_true(is_psw(x))
  expect_false(is_psw(c(0.1, 0.2)))
})

test_that("vec_ptype_abbr and vec_ptype_full return correct type labels", {
  x <- psw(c(0.1, 0.2), estimand = "atm")
  expect_equal(vec_ptype_abbr(x), "psw{atm}")
  expect_equal(vec_ptype_full(x), "psw{estimand = atm}")

  z <- psw(c(0.1, 0.2))
  expect_equal(vec_ptype_abbr(z), "psw")
  expect_equal(vec_ptype_full(z), "psw{estimand = unknown}")

  y <- psw(c(0.1, 0.2), estimand = "cens")
  x <- x * y
  expect_equal(vec_ptype_abbr(x), "psw{atm, cens}")
  expect_equal(vec_ptype_full(x), "psw{estimand = atm, cens}")
})

test_that("vec_cast works for psw and double", {
  x <- psw(c(0.1, 0.2, 0.3), estimand = "ate")

  # Cast to double
  double_x <- vec_cast(x, double())
  expect_equal(double_x, c(0.1, 0.2, 0.3))
  expect_true(is.numeric(double_x))

  # Cast back to psw
  psw_x <- as_psw(double_x, estimand = "ate")
  expect_s3_class(
    psw_x,
    c("psw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vec_data(psw_x), double_x)
  expect_equal(estimand(psw_x), "ate")
})

test_that("vec_cast works for psw and integer with precision checks", {
  x <- psw(c(1, 2, 3), estimand = "ate")

  # Cast to integer
  int_x <- vec_cast(x, integer())
  expect_equal(int_x, as.integer(c(1, 2, 3)))
  expect_true(is.integer(int_x))

  # Cast back to psw
  psw_x <- as_psw(int_x, estimand = "ate")
  expect_s3_class(
    psw_x,
    c("psw", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_equal(vec_data(psw_x), as.numeric(int_x))
  expect_equal(estimand(psw_x), "ate")

  # Fail when precision is lost
  x_with_decimals <- psw(c(1.1, 2.2, 3.3), estimand = "ate")
  expect_propensity_error(
    vec_cast(x_with_decimals, integer())
  )
})

# A prototype carrying a distinct value in each of the six metadata fields, so a
# cast that drops any one of them is visible.
psw_cast_prototype <- function(stabilization_score = 0.42) {
  psw(
    double(),
    estimand = "ate",
    stabilized = TRUE,
    trimmed = TRUE,
    truncated = FALSE,
    calibrated = TRUE,
    stabilization_score = stabilization_score
  )
}

test_that("vec_cast from double carries every metadata field of the psw prototype", {
  to <- psw_cast_prototype()

  out <- vec_cast(c(1, 2, 3), to = to)

  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(1, 2, 3))
  expect_equal(estimand(out), "ate")
  expect_true(is_stabilized(out))
  expect_true(is_ps_trimmed(out))
  expect_false(is_ps_truncated(out))
  expect_true(is_ps_calibrated(out))
  expect_equal(stabilization_score(out), 0.42)
})

test_that("vec_cast from integer carries every metadata field of the psw prototype", {
  to <- psw_cast_prototype()

  out <- vec_cast(c(1L, 2L, 3L), to = to)

  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(1, 2, 3))
  expect_equal(estimand(out), "ate")
  expect_true(is_stabilized(out))
  expect_true(is_ps_trimmed(out))
  expect_false(is_ps_truncated(out))
  expect_true(is_ps_calibrated(out))
  expect_equal(stabilization_score(out), 0.42)
})

test_that("vec_cast from ps_trim carries every metadata field of the psw prototype", {
  to <- psw_cast_prototype()
  ps <- ps_trim(c(0.2, 0.5, 0.8), method = "ps", lower = 0.1, upper = 0.9)

  out <- vec_cast(ps, to = to)

  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(0.2, 0.5, 0.8))
  expect_equal(estimand(out), "ate")
  expect_true(is_stabilized(out))
  expect_true(is_ps_trimmed(out))
  expect_false(is_ps_truncated(out))
  expect_true(is_ps_calibrated(out))
  expect_equal(stabilization_score(out), 0.42)
})

test_that("vec_cast from ps_trunc carries every metadata field of the psw prototype", {
  to <- psw_cast_prototype()
  ps <- ps_trunc(c(0.2, 0.5, 0.8), method = "ps", lower = 0.1, upper = 0.9)

  out <- vec_cast(ps, to = to)

  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(0.2, 0.5, 0.8))
  expect_equal(estimand(out), "ate")
  expect_true(is_stabilized(out))
  expect_true(is_ps_trimmed(out))
  expect_false(is_ps_truncated(out))
  expect_true(is_ps_calibrated(out))
  expect_equal(stabilization_score(out), 0.42)
})

test_that("vec_cast keeps a per-observation stabilization score at a matching length", {
  score <- c(0.51, 0.52, 0.53)
  to <- psw_cast_prototype(stabilization_score = score)

  out <- expect_no_warning(vec_cast(c(1, 2, 3), to = to))

  expect_s3_class(out, "psw")
  expect_equal(stabilization_score(out), score)
  expect_true(is_stabilized(out))
})

test_that("vec_cast drops a per-observation stabilization score at a different length silently", {
  score <- c(0.51, 0.52, 0.53)
  to <- psw_cast_prototype(stabilization_score = score)

  # The score describes the observations behind `to`, not the data being cast,
  # so a length it does not match is nothing the caller can act on.
  out <- expect_silent(vec_cast(c(1, 2), to = to))

  expect_s3_class(out, "psw")
  expect_length(out, 2)
  expect_null(stabilization_score(out))
  expect_null(attr(out, "stabilization_score"))
  expect_true(is_stabilized(out))
  expect_equal(estimand(out), "ate")
})

test_that("vec_cast from psw to double returns the bare data", {
  x <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    trimmed = TRUE,
    calibrated = TRUE,
    stabilization_score = 0.42
  )

  out <- vec_cast(x, to = double())

  expect_identical(out, c(1, 2, 3))
  expect_null(attributes(out))
})

test_that("vec_ptype2 combines psw and other types correctly", {
  x <- psw(c(0.1, 0.2), estimand = "ate")
  y <- psw(c(0.3, 0.4), estimand = "ate")
  expect_equal(vec_ptype2(x, y), new_psw(estimand = "ate"))

  # Different estimands should warn and return numeric
  z <- psw(c(0.5, 0.6), estimand = "att")
  expect_propensity_warning(
    # jarl-ignore implicit_assignment: assignment keeps the return value out of the snapshot while the warning is captured
    result <- vec_ptype2(x, z)
  )
  expect_identical(result, double())

  # Combining with double
  expect_propensity_warning(
    expect_equal(vec_ptype2(x, double()), double())
  )
  expect_propensity_warning(
    expect_equal(vec_ptype2(double(), x), double())
  )

  # Combining with integer
  expect_propensity_warning(
    expect_equal(vec_ptype2(x, integer()), integer())
  )
  expect_propensity_warning(
    expect_equal(vec_ptype2(integer(), x), integer())
  )
})

test_that("vec_arith performs arithmetic on psw objects", {
  x <- psw(c(1, 2, 3), estimand = "ate")
  y <- psw(c(0.5, 1.5, 2.5), estimand = "ate")
  cens <- psw(c(3, 2, 1), estimand = "cens")

  classes <- c("psw", "causal_wts", "vctrs_vctr", "double")

  # same estimand
  result <- x + y
  expect_s3_class(result, classes, exact = TRUE)
  expect_equal(vec_data(result), c(1.5, 3.5, 5.5))
  expect_equal(estimand(result), "ate")

  # different estimand
  result <- x * cens
  expect_s3_class(result, classes, exact = TRUE)
  expect_equal(vec_data(result), c(3, 4, 3))
  expect_equal(estimand(result), "ate, cens")

  # Arithmetic with double
  result <- x * 2
  expect_s3_class(result, classes, exact = TRUE)
  expect_equal(vec_data(result), c(2, 4, 6))
  expect_equal(estimand(result), "ate")

  # Arithmetic with integer
  result <- x - 1L
  expect_s3_class(result, classes, exact = TRUE)
  expect_equal(vec_data(result), c(0, 1, 2))
  expect_equal(estimand(result), "ate")

  # various combos work in various orders
  expect_no_error(x * 2.1)
  expect_no_error(2.1 * x)
  expect_no_error(x / 2.1)
  expect_no_error(2.1 / x)
  expect_no_error(x + 2.1)
  expect_no_error(2.1 + x)
  expect_no_error(x * 2L)
  expect_no_error(2L * x)
  expect_no_error(x / 2L)
  expect_no_error(2L / x)
  expect_no_error(x + 2L)
  expect_no_error(2L + x)
  expect_no_error(x^2)

  # doesn't work with unsupported types
  thing <- new_vctr(runif(10), class = "thing")
  expect_propensity_error(x * thing)
  expect_propensity_error(x * list(runif(10)))
})

test_that("vec_math applies math functions to psw objects", {
  x <- psw(c(1, 4, 9), estimand = "ate")

  expect_false(is_psw(sqrt(x)))
  expect_false(is_psw(sum(x)))
})

test_that("Combination of arithmetic and math works correctly for psw", {
  # Example data
  y <- c(2, 4, 6, 8)
  x <- c(0, 1, 0, 1)
  wts <- psw(c(0.5, 1.0, 1.5, 2.0), estimand = "ate")

  # Compute expected results using double weights
  wts_double <- vec_data(wts) # Strip the psw class to use as raw weights

  expected_term1 <- sum(y * x * wts_double) / sum(x * wts_double)
  expected_term2 <- sum(y * (1 - x) * wts_double) / sum((1 - x) * wts_double)

  # Compute the actual results using psw weights
  term1 <- sum(y * x * wts) / sum(x * wts)
  term2 <- sum(y * (1 - x) * wts) / sum((1 - x) * wts)

  # Validate the results
  expect_equal(term1, expected_term1, tolerance = 1e-8)
  expect_equal(term2, expected_term2, tolerance = 1e-8)

  # Check for consistency in repeated calculations
  repeated_term2 <- sum(y * (1 - x) * wts) / sum((1 - x) * wts)
  expect_equal(repeated_term2, term2, tolerance = 1e-8)
})

test_that("Refit logic can be tracked if stored in attr", {
  # If you store a 'refit' attribute in the final psw object:
  w <- psw(c(0.1, 0.2), estimand = "ate", trimmed = TRUE)
  attr(w, "ps_trim_meta") <- list(refit = TRUE)
  # Suppose we define is_refit.psw() as checking attr(x, "refit")
  expect_true(is_refit(w))
})

test_that("psw objects can convert to character", {
  x <- as.character(new_psw(c(0.1, 0.2, 0.3), estimand = "ate"))
  expect_type(x, "character")
})

test_that("psw works with ggplot2", {
  skip_if_not_installed("ggplot2")
  type <- expect_silent(
    ggplot2::scale_type(psw(1))
  )

  expect_identical(type, "continuous")
})

# Weights built from a modified propensity score carry the modification's
# record: `ps_trim_meta`, `ps_trunc_meta`, or `ps_calib_meta`. Those records are
# indexed by observation, so where an operation goes through vctrs they are kept
# whenever the weights come back at the length the record was written for,
# dropped when they do not, and left alone at zero length.
#
# The drop is silent, because vctrs is not the only route these records take. A
# `[<-` that grows the vector carries one across a length change, and
# `model.frame()` shortens a weights column in C and re-attaches the original
# variable's attributes to it, so weights built on a trimmed propensity score
# come back out of every outcome model fit on them still carrying a record for
# rows that were dropped. A warning on the one route vctrs controls would fire
# on ordinary work while saying nothing about the others.
#
# Honesty lives at query time instead. `is_unit_trimmed()` answers by position,
# so it refuses a record that does not cover the vector it is given, whichever
# route that vector arrived by. `is_refit()` reads one flag rather than a
# position, so it answers from any record present and refuses only when weights
# marked as trimmed have no record at all.
#
# The categorical attributes are not indexed by observation and describe the
# exposure rather than the units, so they follow no length rule and survive
# everything.

# `wt_*()` on a trimmed propensity score warns that the model was never refit.
# These fixtures are about metadata rather than refitting, so that one warning is
# muffled by class instead of silencing the whole call.
without_refit_warning <- function(expr) {
  withCallingHandlers(
    expr,
    propensity_no_refit_warning = function(cnd) invokeRestart("muffleWarning")
  )
}

# Every warning raised while evaluating `expr`, by class, so a case that must
# raise two can name both.
collect_warning_classes <- function(expr) {
  classes <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(cnd) {
      classes <<- c(classes, class(cnd)[[1]])
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, classes = classes)
}

# Units 1 and 5 fall outside the bounds, so `is_unit_trimmed()` has something to
# report and a result that lost the record reports something visibly different.
trimmed_psw <- function(stabilization_score = NULL) {
  ps <- ps_trim(
    c(0.05, 0.3, 0.5, 0.7, 0.95),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )

  without_refit_warning(wt_ate(
    ps,
    c(0, 0, 1, 1, 1),
    exposure_type = "binary",
    .focal_level = 1,
    stabilize = !is.null(stabilization_score),
    stabilization_score = stabilization_score
  ))
}

truncated_psw <- function() {
  ps <- ps_trunc(
    c(0.05, 0.3, 0.5, 0.7, 0.95),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )

  wt_ate(
    ps,
    c(0, 0, 1, 1, 1),
    exposure_type = "binary",
    .focal_level = 1
  )
}

calibrated_psw <- function() {
  exposure <- c(0, 1, 0, 0, 1, 0, 1, 1, 0, 1)
  ps <- ps_calibrate(
    c(0.14, 0.22, 0.31, 0.4, 0.48, 0.55, 0.62, 0.7, 0.78, 0.86),
    exposure
  )

  wt_ate(ps, exposure, exposure_type = "binary", .focal_level = 1)
}

categorical_psw <- function(.focal_level = "A") {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  wt_att(
    ps_matrix,
    exposure,
    .focal_level = .focal_level,
    exposure_type = "categorical"
  )
}

test_that("length-preserving psw arithmetic keeps the trimming record", {
  w <- trimmed_psw()
  meta <- ps_trim_meta(w)
  trimmed_units <- c(TRUE, FALSE, FALSE, FALSE, TRUE)

  expect_identical(is_unit_trimmed(w), trimmed_units)

  results <- list(
    `w * 1` = expect_silent(w * 1),
    `w / sum(w, na.rm = TRUE)` = expect_silent(w / sum(w, na.rm = TRUE)),
    `-w` = expect_silent(-w),
    `2 * w` = expect_silent(2 * w)
  )

  for (label in names(results)) {
    out <- results[[label]]
    expect_s3_class(out, "psw")
    expect_identical(ps_trim_meta(out), meta, info = label)
    expect_true(is_ps_trimmed(out))
    expect_identical(is_unit_trimmed(out), trimmed_units, info = label)
  }
})

test_that("a full-length psw subset keeps the trimming record", {
  w <- trimmed_psw()
  meta <- ps_trim_meta(w)

  whole <- expect_silent(w[seq_along(w)])
  expect_identical(ps_trim_meta(whole), meta)
  expect_identical(is_unit_trimmed(whole), c(TRUE, FALSE, FALSE, FALSE, TRUE))

  sliced <- expect_silent(vec_slice(w, seq_along(w)))
  expect_identical(ps_trim_meta(sliced), meta)
})

test_that("length-preserving psw arithmetic keeps truncation and calibration records", {
  truncated <- truncated_psw()
  trunc_meta <- ps_trunc_meta(truncated)
  out <- expect_silent(truncated * 1)
  expect_identical(ps_trunc_meta(out), trunc_meta)
  expect_true(is_ps_truncated(out))

  calibrated <- calibrated_psw()
  calib_meta <- ps_calib_meta(calibrated)
  out <- expect_silent(calibrated / sum(calibrated))
  expect_identical(ps_calib_meta(out), calib_meta)
  expect_true(is_ps_calibrated(out))
})

test_that("shortening a psw drops the trimming record silently", {
  w <- trimmed_psw()

  sub <- expect_silent(w[1:2])
  expect_s3_class(sub, "psw")
  expect_length(sub, 2)
  expect_null(ps_trim_meta(sub))
  expect_null(attr(sub, "ps_trim_meta"))

  # Everything that is not indexed by observation is untouched by the drop.
  expect_true(is_ps_trimmed(sub))
  expect_identical(estimand(sub), "ate; trimmed")

  sliced <- expect_silent(vec_slice(w, 1:2))
  expect_null(ps_trim_meta(sliced))
  expect_true(is_ps_trimmed(sliced))
})

test_that("shortening a psw drops truncation and calibration records silently", {
  truncated <- truncated_psw()
  sub <- expect_silent(truncated[1:2])
  expect_null(ps_trunc_meta(sub))
  expect_true(is_ps_truncated(sub))

  calibrated <- calibrated_psw()
  sub <- expect_silent(calibrated[1:2])
  expect_null(ps_calib_meta(sub))
  expect_true(is_ps_calibrated(sub))
})

test_that("a zero-length psw restore keeps the trimming record silently", {
  w <- trimmed_psw()
  meta <- ps_trim_meta(w)

  # A prototype or an empty subset holds no observations, so a record on it
  # lines up with nothing and contradicts nothing.
  empty <- expect_silent(w[integer(0)])
  expect_length(empty, 0)
  expect_identical(ps_trim_meta(empty), meta)

  proto <- expect_silent(vec_ptype(w))
  expect_length(proto, 0)
  expect_identical(ps_trim_meta(proto), meta)

  # No observations, no answers. Indexing an empty logical by the positions the
  # record names would grow one as long as the original weights, padded with
  # `NA`, which is neither an answer nor the documented return length.
  expect_identical(is_unit_trimmed(empty), logical(0))
  expect_identical(is_unit_trimmed(proto), logical(0))
})

test_that("shortening a psw with a score and a trimming record warns only for the score", {
  score <- c(0.51, 0.52, 0.53, 0.54, 0.55)
  w <- trimmed_psw(stabilization_score = score)
  expect_identical(stabilization_score(w), score)
  expect_false(is.null(ps_trim_meta(w)))

  out <- collect_warning_classes(w[1:2])

  # A user who recorded a score can recompute the weights on the subset, so that
  # drop is worth saying. The trimming record has a query-time guard instead, and
  # announcing it here would fire on every outcome model fit on these weights.
  expect_identical(out$classes, "propensity_stabilization_score_warning")
  expect_null(stabilization_score(out$value))
  expect_null(ps_trim_meta(out$value))
  expect_true(is_stabilized(out$value))
  expect_true(is_ps_trimmed(out$value))
})

test_that("the trimmed flag survives operations that drop the trimming record", {
  # The flag is vector-level truth about where the weights came from. Losing the
  # per-unit record does not make the propensity scores untrimmed.
  w <- trimmed_psw()

  expect_true(is_ps_trimmed(w * 1))
  expect_true(is_ps_trimmed(w[1:2]))
  expect_true(is_ps_trimmed(vec_slice(w, 1:2)))
  expect_true(is_ps_trimmed(w[integer(0)]))
})

test_that("casting to a psw keeps a length-matched trimming record and drops a shorter one silently", {
  w <- trimmed_psw()
  meta <- ps_trim_meta(w)

  matched <- expect_silent(vec_cast(c(1, 2, 3, 4, 5), to = w))
  expect_s3_class(matched, "psw")
  expect_identical(ps_trim_meta(matched), meta)

  # A cast takes its whole type from `to`, whose record describes `to`'s own
  # observations rather than the incoming data. A length it does not match is
  # nothing the caller can act on, so the record goes without comment.
  shorter <- expect_silent(vec_cast(c(1, 2), to = w))
  expect_s3_class(shorter, "psw")
  expect_length(shorter, 2)
  expect_null(ps_trim_meta(shorter))
})

test_that("subassigning into a decorated psw is silent and keeps its metadata", {
  score <- c(0.51, 0.52, 0.53, 0.54, 0.55)
  x <- trimmed_psw(stabilization_score = score)
  meta <- ps_trim_meta(x)

  # `[<-` casts the replacement to the target's type and then hands off to base
  # R, which assigns into the target and leaves its attributes alone; no restore
  # runs on the result. Only the cast passes through this package, and it sees a
  # length-2 value against a length-5 record, so anything that warned there
  # would warn on an operation that ends with both the score and the record
  # intact and correct.
  expect_silent({
    x[1:2] <- c(0.5, 0.6)
  })

  expect_s3_class(x, "psw")
  expect_length(x, 5)
  expect_identical(vec_data(x)[1:2], c(0.5, 0.6))
  expect_identical(stabilization_score(x), score)
  expect_identical(ps_trim_meta(x), meta)
  expect_identical(is_unit_trimmed(x), c(TRUE, FALSE, FALSE, FALSE, TRUE))
})

test_that("categorical psw attributes survive arithmetic and shortening slices", {
  w <- categorical_psw()

  categorical_attrs <- function(x) {
    list(
      n_categories = attr(x, "n_categories"),
      category_names = attr(x, "category_names"),
      focal_category = attr(x, "focal_category")
    )
  }
  expected <- list(
    n_categories = 3L,
    category_names = c("A", "B", "C"),
    focal_category = "A"
  )

  expect_identical(categorical_attrs(w), expected)
  expect_identical(categorical_attrs(expect_silent(w * 1)), expected)
  expect_identical(categorical_attrs(expect_silent(-w)), expected)

  # These attributes describe the exposure, not the observations, so a shorter
  # result carries them unchanged and has nothing to announce.
  expect_identical(categorical_attrs(expect_silent(w[1:3])), expected)
  expect_identical(
    categorical_attrs(expect_silent(vec_slice(w, 1:3))),
    expected
  )
  expect_identical(categorical_attrs(expect_silent(w[integer(0)])), expected)
})

test_that("a binary focal level survives psw arithmetic and shortening slices", {
  # `ipw()` reads `focal_category` to learn which level att and atu weights
  # target, and reaches the weights through the outcome model, which has already
  # subset and recycled them.
  w <- wt_att(
    c(0.2, 0.4, 0.6, 0.8),
    c("a", "b", "a", "b"),
    .focal_level = "a",
    exposure_type = "binary"
  )

  expect_identical(attr(w, "focal_category"), "a")
  expect_identical(attr(expect_silent(w * 1), "focal_category"), "a")
  expect_identical(attr(expect_silent(w[1:2]), "focal_category"), "a")
  expect_identical(
    attr(expect_silent(vec_slice(w, 1:2)), "focal_category"),
    "a"
  )
})

test_that("is_unit_trimmed() aborts on a trimmed psw with no trimming record", {
  # Without the record there is no way to tell a retained unit from a trimmed
  # one, and reporting every unit as retained is a wrong answer rather than a
  # missing one.
  w <- psw(c(1, 2, 3), estimand = "ate", trimmed = TRUE)

  cnd <- expect_error(
    is_unit_trimmed(w),
    class = "propensity_missing_meta_error"
  )
  expect_s3_class(cnd, "propensity_error")

  # Weights that were never trimmed have nothing missing.
  untrimmed <- psw(c(1, 2, 3), estimand = "ate")
  expect_identical(
    expect_silent(is_unit_trimmed(untrimmed)),
    c(FALSE, FALSE, FALSE)
  )
})

test_that("is_refit() aborts on a trimmed psw with no trimming record", {
  w <- psw(c(1, 2, 3), estimand = "ate", trimmed = TRUE)

  cnd <- expect_error(is_refit(w), class = "propensity_missing_meta_error")
  expect_s3_class(cnd, "propensity_error")

  untrimmed <- psw(c(1, 2, 3), estimand = "ate")
  expect_false(expect_silent(is_refit(untrimmed)))

  # A record that simply does not mention a refit is an answer, not a gap.
  expect_false(is_refit(trimmed_psw()))
})

test_that("weights read out of an outcome model's frame refuse the positional query", {
  # The route the silence design turns on, and the one no subscript in the
  # user's code goes near. `model.frame()` drops the `NA`-weighted rows in C and
  # re-attaches the original variable's attributes to the shortened column, so
  # the weights arrive carrying a record written for rows that are gone. Every
  # row that survived was retained, so the record's positions name trimmed units
  # among units that were all kept.
  set.seed(9)
  n <- 40
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  y <- rbinom(n, 1, 0.4)
  ps <- ps_trim(plogis(0.5 * x), method = "ps", lower = 0.35, upper = 0.65)
  w <- without_refit_warning(
    wt_ate(ps, z, exposure_type = "binary", .focal_level = 1)
  )
  expect_gt(sum(is.na(w)), 0)

  fit <- glm(
    y ~ x,
    data = data.frame(y = y, x = x),
    weights = w,
    family = quasibinomial()
  )
  model_wts <- model.frame(fit)[["(weights)"]]

  expect_s3_class(model_wts, "psw")
  expect_lt(length(model_wts), length(w))
  expect_identical(ps_trim_meta(model_wts), ps_trim_meta(w))

  expect_error(
    is_unit_trimmed(model_wts),
    class = "propensity_missing_meta_error"
  )

  # A refit is one fact about the propensity model rather than about any
  # position, so a record that no longer covers these weights still answers it.
  expect_false(expect_silent(is_refit(model_wts)))

  # Where the weights came from is vector-level provenance and is untouched.
  expect_true(is_ps_trimmed(model_wts))
})

test_that("a psw answers a positional query wherever its ps_trim does", {
  # A subscript holding `NA` takes a position that names no observation, so the
  # record it leaves holds one position fewer than the observations it
  # describes. Coverage is the number of observations the record was written
  # for rather than the total of its positions, which is the test the `ps_trim`
  # applies to itself, so both answer.
  ps <- ps_trim(
    c(0.05, 0.3, 0.5, 0.95),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )
  holed <- ps[c(1, NA, 2)]

  meta <- ps_trim_meta(holed)
  expect_equal(meta$n_obs, 3L)
  expect_equal(length(meta$keep_idx) + length(meta$trimmed_idx), 2L)
  expect_identical(is_unit_trimmed(holed), c(TRUE, FALSE, FALSE))

  w <- without_refit_warning(
    wt_ate(holed, c(1, 0, 1), exposure_type = "binary", .focal_level = 1)
  )
  expect_identical(ps_trim_meta(w), meta)
  expect_identical(is_unit_trimmed(w), c(TRUE, FALSE, FALSE))
})

test_that("growing a psw by subassignment leaves a record that no longer covers it", {
  # `[<-` casts the replacement and then leaves base R to preserve the target's
  # attributes, so the record crosses a length change no restore ever sees. It
  # describes five observations and the weights now hold seven.
  w <- trimmed_psw()
  expect_silent({
    w[7] <- 2
  })

  expect_length(w, 7)
  expect_identical(ps_trim_meta(w), ps_trim_meta(trimmed_psw()))
  expect_true(is_ps_trimmed(w))
  expect_error(is_unit_trimmed(w), class = "propensity_missing_meta_error")
  expect_false(expect_silent(is_refit(w)))
})

# An operation between two psw objects has two sets of metadata to answer for
# rather than one. The six fields describing the weights as a whole merge:
# estimands are pasted when they differ, the result is stabilized only when both
# operands are, and it is marked as trimmed, truncated, or calibrated when either
# operand is. Everything else is carried by agreement. An attribute only one
# operand records has nothing to disagree with, which is what a psw times a plain
# number does through the single-operand route; one both record identically
# describes the result either way; and one they record differently describes
# neither, so it is dropped and named in a warning.
#
# The per-observation length rules are unchanged. A carried record still has to
# cover the number of observations the result arrives at, and one dropped for
# length goes without comment, for the reasons the single-operand routes drop one
# without comment.
#
# Combining with `c()` keeps its own rule. It appends one set of observations to
# another, so the positions a record names would describe units from the other
# input, and nothing rebuilding the combined vector is handed the offsets to
# re-index by. Concatenation therefore drops the modification records whether or
# not the inputs agree on them. The categorical attributes name exposure levels
# rather than positions, so they mean the same thing at the combined length and
# carry when the inputs agree.

# A second trimming record over the same five observations, written at bounds
# that retain a different set, so a merge that kept either record where the two
# disagree would report visibly different units.
alt_trimmed_psw <- function() {
  ps <- ps_trim(
    c(0.05, 0.3, 0.5, 0.7, 0.95),
    method = "ps",
    lower = 0.35,
    upper = 0.65
  )

  without_refit_warning(wt_ate(
    ps,
    c(0, 0, 1, 1, 1),
    exposure_type = "binary",
    .focal_level = 1
  ))
}

categorical_attrs_of <- function(x) {
  list(
    n_categories = attr(x, "n_categories"),
    category_names = attr(x, "category_names"),
    focal_category = attr(x, "focal_category")
  )
}

test_that("a psw product carries a modification record only one operand has", {
  # The workflow the inconsistency showed up in: weights against confounding
  # built on a trimmed propensity score, multiplied by weights against
  # censoring, with the product handed to the outcome model.
  w <- trimmed_psw()
  meta <- ps_trim_meta(w)
  cens <- wt_cens(c(0.8, 0.7, 0.9, 0.6, 0.75), c(1, 1, 0, 1, 1))
  trimmed_units <- c(TRUE, FALSE, FALSE, FALSE, TRUE)

  out <- expect_silent(w * cens)
  expect_s3_class(out, "psw")
  expect_identical(ps_trim_meta(out), meta)
  expect_true(is_ps_trimmed(out))
  expect_identical(is_unit_trimmed(out), trimmed_units)

  # The record travels from whichever operand holds it.
  reversed <- expect_silent(cens * w)
  expect_identical(ps_trim_meta(reversed), meta)
  expect_identical(is_unit_trimmed(reversed), trimmed_units)

  # A psw of ones is the two-operand form of `w * 1`, which keeps the record.
  ones <- expect_silent(w * psw(rep(1, 5), estimand = "ate"))
  expect_identical(ps_trim_meta(ones), meta)
  expect_identical(is_unit_trimmed(ones), trimmed_units)
})

test_that("a psw product carries truncation and calibration records one operand has", {
  truncated <- truncated_psw()
  trunc_meta <- ps_trunc_meta(truncated)
  out <- expect_silent(truncated * psw(rep(1, 5), estimand = "cens"))
  expect_identical(ps_trunc_meta(out), trunc_meta)
  expect_true(is_ps_truncated(out))

  calibrated <- calibrated_psw()
  calib_meta <- ps_calib_meta(calibrated)
  out <- expect_silent(calibrated * psw(rep(1, 10), estimand = "cens"))
  expect_identical(ps_calib_meta(out), calib_meta)
  expect_true(is_ps_calibrated(out))
})

test_that("a psw product carries a modification record both operands share", {
  w <- trimmed_psw()
  meta <- ps_trim_meta(w)

  out <- expect_silent(w * w)
  expect_s3_class(out, "psw")
  expect_identical(ps_trim_meta(out), meta)
  expect_identical(is_unit_trimmed(out), c(TRUE, FALSE, FALSE, FALSE, TRUE))

  # Two separately built objects recording the same thing agree just as well.
  twin <- expect_silent(w * trimmed_psw())
  expect_identical(ps_trim_meta(twin), meta)
})

test_that("a psw product drops conflicting trimming records with one warning", {
  w <- trimmed_psw()
  alt <- alt_trimmed_psw()
  expect_false(identical(ps_trim_meta(w), ps_trim_meta(alt)))

  out <- collect_warning_classes(w * alt)
  expect_identical(out$classes, "propensity_metadata_conflict_warning")

  expect_s3_class(out$value, "psw")
  expect_null(ps_trim_meta(out$value))
  expect_true(is_ps_trimmed(out$value))

  # Neither record describes the product, so the positional query has nothing to
  # answer from and refuses rather than naming units from one of them.
  expect_error(
    is_unit_trimmed(out$value),
    class = "propensity_missing_meta_error"
  )

  cnd <- expect_warning(
    w * alt,
    class = "propensity_metadata_conflict_warning"
  )
  expect_true(grepl("ps_trim_meta", conditionMessage(cnd), fixed = TRUE))
})

test_that("a psw product carries a stabilization score both operands record", {
  x <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 0.4
  )
  y <- psw(
    c(2, 2, 2),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 0.4
  )

  # A result marked as stabilized with the score gone reads as one stabilized by
  # the default score rather than by the one the user supplied.
  out <- expect_silent(x * y)
  expect_identical(stabilization_score(out), 0.4)
  expect_true(is_stabilized(out))

  score <- c(0.51, 0.52, 0.53)
  a <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = score
  )
  b <- psw(
    c(2, 2, 2),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = score
  )
  per_observation <- expect_silent(a * b)
  expect_identical(stabilization_score(per_observation), score)

  # One operand recording a score has nothing to disagree with.
  alone <- expect_silent(
    x * psw(c(2, 2, 2), estimand = "cens", stabilized = TRUE)
  )
  expect_identical(stabilization_score(alone), 0.4)
})

test_that("a psw product carries a stabilization score only when it is stabilized", {
  x <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 0.4
  )

  # Stabilization takes both operands, and the score is the value the weights
  # were stabilized on, so a result that is not stabilized has nothing for it to
  # describe. Carrying it would report weights as stabilized on a score by the
  # side of a flag saying they are not stabilized at all.
  out <- expect_silent(x * psw(c(2, 2, 2), estimand = "cens"))
  expect_false(is_stabilized(out))
  expect_null(stabilization_score(out))

  reversed <- expect_silent(psw(c(2, 2, 2), estimand = "cens") * x)
  expect_false(is_stabilized(reversed))
  expect_null(stabilization_score(reversed))

  # Two scores that disagree are not reported as a disagreement when the result
  # drops the score for its status either way.
  other <- psw(c(2, 2, 2), estimand = "ate", stabilization_score = 0.5)
  both <- expect_silent(x * other)
  expect_false(is_stabilized(both))
  expect_null(stabilization_score(both))
})

test_that("a psw product drops conflicting stabilization scores with one warning", {
  x <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 0.4
  )
  y <- psw(
    c(2, 2, 2),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 0.5
  )

  out <- collect_warning_classes(x * y)
  expect_identical(out$classes, "propensity_metadata_conflict_warning")
  expect_null(stabilization_score(out$value))
  expect_true(is_stabilized(out$value))

  # A score recorded per observation is a different score from a single one.
  per_observation <- psw(
    c(2, 2, 2),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = c(0.4, 0.4, 0.4)
  )
  lengths <- collect_warning_classes(x * per_observation)
  expect_identical(lengths$classes, "propensity_metadata_conflict_warning")
  expect_null(stabilization_score(lengths$value))

  cnd <- expect_warning(x * y, class = "propensity_metadata_conflict_warning")
  expect_true(
    grepl("stabilization_score", conditionMessage(cnd), fixed = TRUE)
  )
})

test_that("a psw product carries categorical attributes the operands share", {
  w <- categorical_psw()
  expected <- list(
    n_categories = 3L,
    category_names = c("A", "B", "C"),
    focal_category = "A"
  )

  expect_identical(categorical_attrs_of(expect_silent(w * w)), expected)
  expect_identical(
    categorical_attrs_of(expect_silent(w * categorical_psw())),
    expected
  )

  # Only one operand describing the exposure is enough to describe the product.
  plain <- psw(rep(1, 6), estimand = "cens")
  expect_identical(categorical_attrs_of(expect_silent(w * plain)), expected)
  expect_identical(categorical_attrs_of(expect_silent(plain * w)), expected)
})

test_that("a psw product drops a conflicting focal category and keeps the rest", {
  a <- categorical_psw()
  b <- categorical_psw(.focal_level = "B")

  out <- collect_warning_classes(a * b)
  expect_identical(out$classes, "propensity_metadata_conflict_warning")

  # The two describe the same exposure and disagree only on which level the
  # weights target, so only that attribute is dropped.
  expect_identical(
    categorical_attrs_of(out$value),
    list(
      n_categories = 3L,
      category_names = c("A", "B", "C"),
      focal_category = NULL
    )
  )

  cnd <- expect_warning(a * b, class = "propensity_metadata_conflict_warning")
  message <- conditionMessage(cnd)
  expect_true(grepl("focal_category", message, fixed = TRUE))
  expect_false(grepl("category_names", message, fixed = TRUE))
})

test_that("the six core psw fields merge unchanged through a product", {
  x <- psw(c(1, 2), estimand = "ate", stabilized = TRUE, trimmed = TRUE)
  y <- psw(
    c(3, 4),
    estimand = "cens",
    stabilized = TRUE,
    truncated = TRUE,
    calibrated = TRUE
  )

  out <- expect_silent(x * y)
  expect_identical(estimand(out), "ate, cens")
  expect_true(is_stabilized(out))
  expect_true(is_ps_trimmed(out))
  expect_true(is_ps_truncated(out))
  expect_true(is_ps_calibrated(out))

  # A shared estimand is not pasted to itself.
  matched <- expect_silent(psw(c(1, 2), estimand = "ate") * y)
  expect_identical(estimand(matched), "ate, cens")
  same <- expect_silent(
    psw(c(1, 2), estimand = "ate") * psw(c(3, 4), estimand = "ate")
  )
  expect_identical(estimand(same), "ate")

  # Stabilization holds only when both operands are stabilized; where the
  # propensity scores came from holds when either operand records it.
  mixed <- expect_silent(x * psw(c(3, 4), estimand = "ate"))
  expect_false(is_stabilized(mixed))
  expect_true(is_ps_trimmed(mixed))
  expect_false(is_ps_truncated(mixed))
})

test_that("a psw product records an estimand only one operand names", {
  named <- psw(c(1, 2), estimand = "ate")
  unnamed <- psw(c(1, 1))

  # Only two estimands are pasted together. An operand naming none contributes
  # no half to paste, so the label the other supplies stands for the result
  # rather than trailing a separator with nothing after it.
  expect_identical(estimand(expect_silent(named * unnamed)), "ate")
  expect_identical(estimand(expect_silent(unnamed * named)), "ate")

  expect_null(estimand(expect_silent(unnamed * unnamed)))
})

test_that("combining psw objects drops the modification records", {
  w <- trimmed_psw()

  out <- expect_silent(c(w, w))
  expect_s3_class(out, "psw")
  expect_length(out, 10)
  expect_null(ps_trim_meta(out))

  # Everything that is not indexed by observation survives the concatenation.
  expect_true(is_ps_trimmed(out))
  expect_identical(estimand(out), "ate; trimmed")

  expect_error(is_unit_trimmed(out), class = "propensity_missing_meta_error")

  truncated <- truncated_psw()
  combined <- expect_silent(c(truncated, truncated))
  expect_null(ps_trunc_meta(combined))
  expect_true(is_ps_truncated(combined))

  calibrated <- calibrated_psw()
  combined <- expect_silent(c(calibrated, calibrated))
  expect_null(ps_calib_meta(combined))
  expect_true(is_ps_calibrated(combined))
})

test_that("combining psw objects drops modification records that disagree without comment", {
  w <- trimmed_psw()
  alt <- alt_trimmed_psw()
  expect_false(identical(ps_trim_meta(w), ps_trim_meta(alt)))

  # Concatenation drops every modification record for a reason that has nothing
  # to do with whether the inputs agree, so two that disagree leave nothing to
  # report. A warning here would name a record the result would have lost had
  # the two been identical.
  out <- expect_silent(c(w, alt))
  expect_s3_class(out, "psw")
  expect_length(out, 10)
  expect_null(ps_trim_meta(out))
  expect_true(is_ps_trimmed(out))
})

test_that("combining psw objects carries categorical attributes the inputs share", {
  w <- categorical_psw()
  expected <- list(
    n_categories = 3L,
    category_names = c("A", "B", "C"),
    focal_category = "A"
  )

  out <- expect_silent(c(w, w))
  expect_s3_class(out, "psw")
  expect_length(out, 12)
  expect_identical(categorical_attrs_of(out), expected)

  separate <- expect_silent(c(w, categorical_psw()))
  expect_identical(categorical_attrs_of(separate), expected)
})

test_that("combining psw objects drops a conflicting focal category", {
  a <- categorical_psw()
  b <- categorical_psw(.focal_level = "B")

  out <- collect_warning_classes(c(a, b))
  expect_identical(out$classes, "propensity_metadata_conflict_warning")

  expect_s3_class(out$value, "psw")
  expect_length(out$value, 12)
  expect_identical(
    categorical_attrs_of(out$value),
    list(
      n_categories = 3L,
      category_names = c("A", "B", "C"),
      focal_category = NULL
    )
  )
})

test_that("a conflicting attribute stays dropped however the inputs are ordered", {
  a <- categorical_psw()
  b <- categorical_psw(.focal_level = "B")
  dropped <- list(
    n_categories = 3L,
    category_names = c("A", "B", "C"),
    focal_category = NULL
  )

  # vctrs settles the common type two inputs at a time, so a third input meets a
  # prototype the earlier pair has already taken the attribute off. Carrying it
  # back from an input that has nothing to disagree with by then would make the
  # result depend on the order the inputs were written in and would leave the
  # warning naming an attribute the result still has.
  forward <- collect_warning_classes(c(a, b, b))
  backward <- collect_warning_classes(c(b, b, a))

  expect_identical(forward$classes, "propensity_metadata_conflict_warning")
  expect_identical(backward$classes, "propensity_metadata_conflict_warning")
  expect_identical(categorical_attrs_of(forward$value), dropped)
  expect_identical(categorical_attrs_of(backward$value), dropped)
  expect_length(forward$value, 18)
  expect_length(backward$value, 18)

  # A disagreement the operation has already reported is not reported again.
  repeated <- collect_warning_classes(c(a, b, a, b))
  expect_identical(repeated$classes, "propensity_metadata_conflict_warning")
  expect_identical(categorical_attrs_of(repeated$value), dropped)
})

test_that("the record of what an operation dropped stays off its result", {
  a <- categorical_psw()
  b <- categorical_psw(.focal_level = "B")

  # The record rides the prototype from one pair of inputs to the next. A result
  # carrying it would answer for a merge it took no part in, and would name an
  # attribute in nothing a reader of the weights can interpret.
  results <- list(
    suppressWarnings(c(a, b)),
    suppressWarnings(c(a, b, b)),
    suppressWarnings(c(b, b, a)),
    suppressWarnings(a * b),
    c(a, a),
    a * a
  )

  for (out in results) {
    expect_false(psw_conflicted_attr %in% names(attributes(out)))
  }

  expect_setequal(
    names(attributes(suppressWarnings(c(a, b)))),
    c(
      "estimand",
      "stabilized",
      "trimmed",
      "truncated",
      "calibrated",
      "class",
      "n_categories",
      "category_names"
    )
  )
})

test_that("a conflict warning is attributed to the operation the caller wrote", {
  a <- categorical_psw()
  b <- categorical_psw(.focal_level = "B")

  # Arithmetic reaches its method through vctrs' dispatch, whose call names
  # nothing the caller wrote, so the warning reports no call at all.
  product <- expect_warning(
    a * b,
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(conditionCall(product))

  # Combining reports the method, which is where the coercion warnings this
  # class already raises on the same route report from.
  combined <- expect_warning(
    c(a, b),
    class = "propensity_metadata_conflict_warning"
  )
  expect_identical(
    as.character(conditionCall(combined)[[1]]),
    "vec_ptype2.psw.psw"
  )
})

# ---- validating the estimand ------------------------------------------------
#
# The estimand is a name, and everything that reads it back reads it as one. A
# value of another type was recorded as it was given and compared against sets
# of names through `as.character()`, which matches the name it prints as, so it
# passed every check and then arrived at a `switch()` in its own type: a factor
# is its integer level code there, which selects a branch by position.

test_that("psw() rejects an estimand that is not a single string", {
  expect_error(
    psw(c(0.1, 0.2), estimand = list("ate")),
    class = "propensity_estimand_type_error"
  )
  expect_error(
    psw(c(0.1, 0.2), estimand = factor("ate")),
    class = "propensity_estimand_type_error"
  )
  expect_error(
    psw(c(0.1, 0.2), estimand = c("ate", "att")),
    class = "propensity_estimand_type_error"
  )
  expect_error(
    psw(c(0.1, 0.2), estimand = 1),
    class = "propensity_estimand_type_error"
  )

  cnd <- expect_error(psw(c(0.1, 0.2), estimand = list("ate")))
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "single string", fixed = TRUE)
  expect_match(message, "<list>", fixed = TRUE)
})

test_that("psw() accepts the estimands weights are built with", {
  for (estimand in ipw_estimands) {
    expect_identical(estimand(psw(c(0.1, 0.2), estimand = estimand)), estimand)
  }
  # the class does not own the set of names, so a name no weight function
  # produces is still recorded; what reads it decides whether it knows it
  expect_identical(estimand(psw(c(0.1, 0.2), estimand = "banana")), "banana")
  expect_null(estimand(psw(c(0.1, 0.2))))
})

# ---- validating a stabilization score ---------------------------------------

test_that("psw() rejects a stabilization score of the wrong length", {
  expect_error(
    psw(c(1, 2, 3, 4), estimand = "ate", stabilization_score = c(1, 2)),
    class = "propensity_stabilization_score_error"
  )
  expect_error(
    psw(c(1, 2, 3, 4), estimand = "ate", stabilization_score = c(1, 2, 3)),
    class = "propensity_stabilization_score_error"
  )
  expect_error(
    psw(c(1, 2, 3, 4), estimand = "ate", stabilization_score = numeric(0)),
    class = "propensity_stabilization_score_error"
  )

  cnd <- expect_error(
    psw(c(1, 2, 3, 4), estimand = "ate", stabilization_score = c(1, 2))
  )
  message <- conditionMessage(cnd)
  expect_true(grepl("stabilization_score", message, fixed = TRUE))
  expect_true(grepl("2 values", message, fixed = TRUE))
  expect_true(grepl("4 observations", message, fixed = TRUE))
})

test_that("psw() rejects a stabilization score that is not positive and finite", {
  scored <- function(score) {
    psw(c(1, 2, 3), estimand = "ate", stabilization_score = score)
  }

  expect_error(scored(0), class = "propensity_stabilization_score_error")
  expect_error(scored(-1), class = "propensity_stabilization_score_error")
  expect_error(
    scored(NA_real_),
    class = "propensity_stabilization_score_error"
  )
  expect_error(scored(NaN), class = "propensity_stabilization_score_error")
  expect_error(scored(Inf), class = "propensity_stabilization_score_error")
  expect_error(scored(-Inf), class = "propensity_stabilization_score_error")
  expect_error(
    scored(c(0.5, -0.5, 0.5)),
    class = "propensity_stabilization_score_error"
  )
  expect_error(
    scored(c(0.5, 0.5, NA_real_)),
    class = "propensity_stabilization_score_error"
  )

  cnd <- expect_error(scored(c(0.5, -0.5, 0.5)))
  message <- conditionMessage(cnd)
  expect_true(grepl("stabilization_score", message, fixed = TRUE))
  expect_true(grepl("position 2", message, fixed = TRUE))
})

test_that("psw() rejects a non-numeric stabilization score", {
  scored <- function(score) {
    psw(c(1, 2, 3), estimand = "ate", stabilization_score = score)
  }

  expect_error(scored("0.5"), class = "propensity_stabilization_score_error")
  expect_error(
    scored(list(0.5)),
    class = "propensity_stabilization_score_error"
  )
  expect_error(scored(TRUE), class = "propensity_stabilization_score_error")

  cnd <- expect_error(scored("0.5"))
  expect_true(
    grepl("stabilization_score", conditionMessage(cnd), fixed = TRUE)
  )
})

test_that("psw() accepts a scalar and a per-observation stabilization score", {
  scalar <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 0.4
  )
  expect_identical(stabilization_score(scalar), 0.4)

  score <- c(0.4, 0.5, 0.6)
  per_observation <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = score
  )
  expect_identical(stabilization_score(per_observation), score)
})

test_that("psw() stores an integer stabilization score as a double", {
  scalar <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 1L
  )
  expect_identical(stabilization_score(scalar), 1)
  expect_type(stabilization_score(scalar), "double")

  per_observation <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = c(1L, 2L, 3L)
  )
  expect_identical(stabilization_score(per_observation), c(1, 2, 3))

  # Names would reintroduce the difference the storage type no longer makes,
  # and nothing reads them, so a score is recorded as a plain double vector.
  named <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = c(a = 1, b = 2, c = 3)
  )
  expect_identical(stabilization_score(named), c(1, 2, 3))
})

test_that("an integer and a double stabilization score are the same score", {
  # The metadata merge compares scores with `identical()`, which reads a
  # difference in storage type as a difference in the score.
  from_integer <- psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 1L
  )
  from_double <- psw(
    c(2, 2, 2),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = 1
  )

  product <- expect_silent(from_integer * from_double)
  expect_s3_class(product, "psw")
  expect_identical(stabilization_score(product), 1)

  combined <- expect_silent(c(from_integer, from_double))
  expect_s3_class(combined, "psw")
  expect_identical(stabilization_score(combined), 1)
})

test_that("a psw prototype records a score for observations it does not hold", {
  # A zero-length psw carries metadata for observations that have not arrived,
  # so there is no length for a per-observation score to be checked against.
  # `vec_cast()` checks it against the length the data does arrive at.
  score <- c(0.51, 0.52, 0.53)
  proto <- expect_silent(psw(
    double(),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = score
  ))

  expect_identical(stabilization_score(proto), score)
  expect_identical(stabilization_score(vec_cast(c(1, 2, 3), to = proto)), score)
  expect_null(stabilization_score(vec_cast(c(1, 2), to = proto)))

  # The value checks apply to a prototype like anything else.
  expect_error(
    psw(double(), estimand = "ate", stabilization_score = c(0.5, -0.5)),
    class = "propensity_stabilization_score_error"
  )
})

test_that("the stabilization score rejections read clearly", {
  expect_propensity_error(
    psw(c(1, 2, 3, 4), estimand = "ate", stabilization_score = c(1, 2))
  )
  expect_propensity_error(
    psw(c(1, 2, 3), estimand = "ate", stabilization_score = c(0.5, -0.5, 0.5))
  )
  expect_propensity_error(
    psw(c(1, 2, 3), estimand = "ate", stabilization_score = 0)
  )
  expect_propensity_error(
    psw(c(1, 2, 3), estimand = "ate", stabilization_score = "0.5")
  )
})
