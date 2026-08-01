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
# indexed by observation, so they are kept whenever the weights come back at the
# length the record was written for, dropped when they do not, and left alone at
# zero length.
#
# The drop is silent. `model.frame()` slices `NA`-weighted rows out of the
# weights, so weights built on a trimmed propensity score are shortened inside
# every outcome model fit on them, with no subscript anywhere in the user's code
# and no indices available to a restore. A warning there would fire on ordinary
# work and say nothing actionable. Honesty lives at query time instead:
# `is_unit_trimmed()` and `is_refit()` refuse to answer once the record is gone
# rather than reporting every unit as retained.
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

categorical_psw <- function() {
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
    .focal_level = "A",
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

  # Subassignment casts the replacement to the target's type and then restores
  # the result at full length. The intermediate cast is shorter than the record,
  # so anything that warned there would warn on an operation that ends with both
  # the score and the record intact.
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
