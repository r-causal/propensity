valid_trunc_matrix_fixture <- function() {
  exposure <- factor(c("a", "b", "c", "a", "b", "c"))
  # every row sums to 1 and no cell falls below 0.01, so the default threshold
  # leaves the scores as they are
  ps_matrix <- rbind(
    c(0.50, 0.30, 0.20),
    c(0.30, 0.50, 0.20),
    c(0.20, 0.30, 0.50),
    c(0.40, 0.40, 0.20),
    c(0.25, 0.35, 0.40),
    c(0.34, 0.33, 0.33)
  )
  colnames(ps_matrix) <- levels(exposure)

  list(exposure = exposure, ps_matrix = ps_matrix)
}

test_that("ps_trunc works with matrix propensity scores", {
  # Create test data
  set.seed(123)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.02, 0.8), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix) # Normalize to sum to 1
  colnames(ps_matrix) <- levels(exposure)

  # Test basic truncation
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )

  expect_s3_class(
    truncated,
    c("ps_trunc_matrix", "ps_trunc", "matrix"),
    exact = TRUE
  )
  expect_equal(dim(truncated), dim(ps_matrix))

  # Check metadata
  meta <- ps_trunc_meta(truncated)
  expect_equal(meta$method, "ps")
  expect_equal(meta$lower_bound, 0.05)
  expect_true(meta$is_matrix)

  # After truncation and renormalization, some values might still be below threshold
  # This is expected behavior - we truncate once and renormalize

  # Check that rows still sum to 1
  row_sums <- rowSums(truncated)
  expect_equal(row_sums, rep(1, n), tolerance = 1e-10)

  # Check which rows were truncated
  had_low_vals <- apply(ps_matrix, 1, function(x) any(x < 0.05))
  expect_equal(sort(meta$truncated_idx), sort(which(had_low_vals)))
})

test_that("ps_trunc works with data.frame propensity scores", {
  set.seed(456)
  n <- 50
  exposure <- factor(sample(c("Low", "Med", "High"), n, replace = TRUE))

  ps_df <- data.frame(
    Low = runif(n, 0.01, 0.6),
    Med = runif(n, 0.02, 0.5),
    High = runif(n, 0.01, 0.4)
  )
  # Normalize rows
  ps_df <- ps_df / rowSums(ps_df)

  truncated <- ps_trunc(ps_df, .exposure = exposure, method = "ps", lower = 0.1)

  expect_s3_class(
    truncated,
    c("ps_trunc_matrix", "ps_trunc", "matrix"),
    exact = TRUE
  )
  expect_equal(nrow(truncated), nrow(ps_df))
  expect_equal(ncol(truncated), ncol(ps_df))

  # After truncation and renormalization, some values might still be below threshold
  # This is expected behavior
})

test_that("ps_trunc preserves row sums equal to 1", {
  set.seed(789)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C", "D"), n, replace = TRUE))

  # Create propensity scores with some extreme values
  ps_matrix <- matrix(runif(n * 4, 0.001, 0.8), nrow = n, ncol = 4)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Apply truncation
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )

  # Check all rows sum to 1
  row_sums <- rowSums(truncated)
  expect_true(all(abs(row_sums - 1) < 1e-10))

  # Original and truncated should have same dimensions
  expect_equal(dim(truncated), dim(ps_matrix))
})

test_that("ps_trunc validates delta < 1/k", {
  n <- 50
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # delta >= 1/3 should error
  expect_propensity_error(
    ps_trunc(ps_matrix, .exposure = exposure, method = "ps", lower = 0.35)
  )
})

test_that("ps_trunc errors for unsupported methods with categorical", {
  n <- 20
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Only cr should error - pctl is now supported
  expect_propensity_error(
    ps_trunc(ps_matrix, .exposure = exposure, method = "cr")
  )
})

test_that("ps_trunc requires exposure for categorical", {
  n <- 20
  ps_matrix <- matrix(runif(n * 3), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c("A", "B", "C")

  expect_propensity_error(
    ps_trunc(ps_matrix, method = "ps")
  )
})

test_that("is_ps_truncated works for matrix objects", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  expect_false(is_ps_truncated(ps_matrix))

  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )
  expect_true(is_ps_truncated(truncated))
})

test_that("is_unit_truncated works for matrix objects", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.1
  )
  unit_truncated <- is_unit_truncated(truncated)

  expect_equal(length(unit_truncated), n)
  expect_type(unit_truncated, "logical")

  meta <- ps_trunc_meta(truncated)
  expect_equal(which(unit_truncated), meta$truncated_idx)
})

test_that("ps_trunc handles extreme values correctly", {
  n <- 20
  exposure <- factor(c(rep("A", 10), rep("B", 5), rep("C", 5)))

  # Create matrix with very small values
  ps_matrix <- matrix(0, nrow = n, ncol = 3)
  ps_matrix[1:10, ] <- c(0.001, 0.499, 0.5) # Group A
  ps_matrix[11:15, ] <- c(0.5, 0.001, 0.499) # Group B
  ps_matrix[16:20, ] <- c(0.499, 0.5, 0.001) # Group C
  colnames(ps_matrix) <- c("A", "B", "C")

  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.01
  )

  # After truncation and renormalization, some values might still be below threshold
  # This is expected behavior

  # Check rows sum to 1
  expect_equal(rowSums(truncated), rep(1, n), tolerance = 1e-10)

  # All rows should have been truncated
  meta <- ps_trunc_meta(truncated)
  expect_equal(length(meta$truncated_idx), n)
})

test_that("ps_trunc handles parsnip-style column names", {
  n <- 40
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c(".pred_A", ".pred_B", ".pred_C")

  truncated <- expect_no_error(
    ps_trunc(
      ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.05
    )
  )

  expect_s3_class(
    truncated,
    c("ps_trunc_matrix", "ps_trunc", "matrix"),
    exact = TRUE
  )
})

test_that("ps_trunc warns when no column names provided", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  # No column names

  expect_propensity_warning(
    ps_trunc(
      ps_matrix,
      .exposure = exposure,
      method = "ps",
      lower = 0.05
    )
  )
})

test_that("ps_trunc.ps_trunc warns about already truncated scores", {
  n <- 30
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  ps_matrix <- matrix(runif(n * 3, 0.01, 0.9), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  truncated_once <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.05
  )

  truncated_twice <- expect_propensity_warning(
    ps_trunc(
      truncated_once,
      .exposure = exposure,
      method = "ps",
      lower = 0.1
    )
  )

  # Should return original truncated object
  expect_identical(truncated_twice, truncated_once)
})

# PSweight comparison tests
test_that("ps_trunc matches PSweight truncation behavior for entropy weights", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # PSweight automatically truncates for entropy weights
  # It uses fixed bounds [1e-6, 1-1e-6]

  set.seed(123)
  n <- 100

  # Create extreme propensity scores that sum to 1
  ps_matrix <- matrix(0, nrow = n, ncol = 3)
  # Use rep to fill rows properly
  ps_matrix[1:30, ] <- matrix(
    rep(c(1e-7, 0.5, 0.5 - 1e-7), each = 30),
    ncol = 3,
    byrow = TRUE
  )
  ps_matrix[31:60, ] <- matrix(
    rep(c(0.3, 0.3, 0.4), each = 30),
    ncol = 3,
    byrow = TRUE
  )
  ps_matrix[61:100, ] <- matrix(
    rep(c(0.99999, 5e-6, 5e-6), each = 40),
    ncol = 3,
    byrow = TRUE
  )
  # Normalize to ensure rows sum to 1
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- c("A", "B", "C")

  exposure <- factor(c(rep("A", 30), rep("B", 30), rep("C", 40)))

  # Our truncation with PSweight's threshold
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 1e-6
  )

  # Check bounds match PSweight's approach (excluding NAs)
  expect_true(all(truncated >= 1e-6, na.rm = TRUE))
  expect_true(all(truncated <= 1 - 1e-6, na.rm = TRUE))

  # Rows should still sum to 1
  expect_equal(rowSums(truncated), rep(1, n), tolerance = 1e-10)
})

test_that("ps_trunc with method='pctl' works for categorical exposures", {
  set.seed(789)
  n <- 100
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  # Create propensity scores with some extreme values
  ps_matrix <- matrix(runif(n * 3, 0.001, 0.999), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Apply percentile truncation
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "pctl",
    lower = 0.05,
    upper = 0.95
  )

  expect_s3_class(
    truncated,
    c("ps_trunc_matrix", "ps_trunc", "matrix"),
    exact = TRUE
  )
  expect_equal(dim(truncated), dim(ps_matrix))

  # Check metadata
  meta <- ps_trunc_meta(truncated)
  expect_equal(meta$method, "pctl")
  expect_equal(meta$lower_pctl, 0.05)
  expect_equal(meta$upper_pctl, 0.95)
  expect_true(!is.na(meta$lower_bound))
  expect_true(!is.na(meta$upper_bound))

  # Check that rows still sum to 1
  row_sums <- rowSums(truncated)
  expect_equal(row_sums, rep(1, n), tolerance = 1e-10)

  # The thresholds should be based on the overall distribution
  all_ps_vals <- as.vector(ps_matrix)
  # `quantile()` names its result for the probability it was asked for; the
  # bound is recorded without that name.
  expected_lower <- unname(quantile(all_ps_vals, 0.05))
  expected_upper <- unname(quantile(all_ps_vals, 0.95))
  expect_equal(meta$lower_bound, expected_lower)
  expect_equal(meta$upper_bound, expected_upper)
})

test_that("categorical truncation preserves relative proportions", {
  # When truncating, the relative proportions of non-truncated values should be preserved
  set.seed(456)
  n <- 50

  ps_matrix <- matrix(
    c(
      0.001,
      0.599,
      0.4, # First value needs truncation
      0.3,
      0.001,
      0.699, # Second value needs truncation
      0.4,
      0.4,
      0.2
    ), # No truncation needed
    ncol = 3,
    byrow = TRUE
  )
  ps_matrix <- ps_matrix[rep(1:3, length.out = n), ]

  exposure <- factor(rep(c("A", "B", "C"), length.out = n))
  colnames(ps_matrix) <- levels(exposure)

  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.01
  )

  # For first row type: after truncating 0.001 to 0.01,
  # the ratio of second to third should be preserved
  row1_idx <- seq(1, n, by = 3)
  if (length(row1_idx) > 0) {
    original_ratio <- ps_matrix[row1_idx[1], 2] / ps_matrix[row1_idx[1], 3]
    new_ratio <- truncated[row1_idx[1], 2] / truncated[row1_idx[1], 3]
    expect_equal(unname(new_ratio), unname(original_ratio), tolerance = 1e-10)
  }
})

test_that("print.ps_trunc_matrix produces expected output", {
  set.seed(456)
  n <- 30
  exposure <- factor(sample(c("Low", "Med", "High"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.05, 0.95), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  # Test ps method
  truncated <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "ps",
    lower = 0.1
  )
  output <- capture.output(print(truncated, n = 3))

  expect_match(
    output[1],
    "<ps_trunc_matrix\\[30 x 3\\]; truncated \\d+ of 30; method=ps\\[0\\.1000,Inf\\]>"
  )
  expect_match(output[2], "High\\s+Low\\s+Med")
  expect_match(output[3], "\\[1,\\]")

  # Test pctl method
  truncated_pctl <- ps_trunc(
    ps_matrix,
    .exposure = exposure,
    method = "pctl",
    lower = 0.05,
    upper = 0.95
  )
  output_pctl <- capture.output(print(truncated_pctl))

  expect_match(
    output_pctl[1],
    "<ps_trunc_matrix\\[30 x 3\\]; truncated \\d+ of 30; method=pctl\\[0\\.05,0\\.95\\]>"
  )

  # Test without column names
  ps_matrix_no_names <- ps_matrix
  colnames(ps_matrix_no_names) <- NULL
  suppressWarnings({
    truncated_no_names <- ps_trunc(
      ps_matrix_no_names,
      .exposure = exposure,
      method = "ps",
      lower = 0.1
    )
  })
  output_no_names <- capture.output(print(truncated_no_names))

  # Should not have column header
  expect_false(any(grepl("High\\s+Low\\s+Med", output_no_names)))
})

# What the matrix print method shows -----------------------------------------

# The truncation record is metadata, and `ps_trunc_meta()` is how a caller reads
# it. The print method summarizes that record in the header and then shows the
# scores, and the scores are what it shows: the record is a set of index vectors
# as long as the data, so printed after the matrix it is a second and longer
# object the caller did not ask to see.

trunc_print_fixture <- function() {
  set.seed(123)
  n <- 20
  exposure <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  ps_matrix <- matrix(runif(n * 3, 0.05, 0.95), nrow = n, ncol = 3)
  ps_matrix <- ps_matrix / rowSums(ps_matrix)
  colnames(ps_matrix) <- levels(exposure)

  ps_trunc(ps_matrix, .exposure = exposure, method = "ps", lower = 0.1)
}

test_that("print.ps_trunc_matrix() shows some rows without the record", {
  truncated <- trunc_print_fixture()
  output <- capture.output(print(truncated, n = 3))

  expect_false(any(grepl('attr\\(,"ps_trunc_meta"\\)', output)))

  # Header, column names, the rows asked for, and the count of the rest.
  expect_length(output, 6)
  expect_match(output[1], "<ps_trunc_matrix\\[20 x 3\\]; truncated \\d+ of 20;")
  expect_equal(sum(grepl("^\\s*\\[\\s*\\d+,\\]", output)), 3)
  expect_true(any(grepl("# \\.\\.\\. with 17 more rows", output)))
})

test_that("print.ps_trunc_matrix() shows every row without the record", {
  truncated <- trunc_print_fixture()
  output <- capture.output(print(truncated, n = Inf))

  expect_false(any(grepl('attr\\(,"ps_trunc_meta"\\)', output)))

  # Header, column names, and every row.
  expect_length(output, 22)
  expect_equal(sum(grepl("^\\s*\\[\\s*\\d+,\\]", output)), 20)
})

# `method = "ps"` clamps the categorical scores from below and renormalizes the
# row, so it settles on no upper bound and records none. The header reports the
# bounds the method worked to, and a bound it did not set is not one it can
# report as a number: what stands in for it has to read as its absence.

test_that("the ps_trunc_matrix header reports a bound the method never set", {
  truncated <- trunc_print_fixture()
  header <- capture.output(print(truncated, n = 3))[1]

  expect_true(is.na(ps_trunc_meta(truncated)$upper_bound))
  expect_false(grepl("\\bNA\\b", header))
  expect_false(grepl("NaN", header))
  expect_match(header, "method=ps")
  expect_match(header, "0.1000", fixed = TRUE)
})

# ---- the matrix method names ps_trunc, whoever called it -------------------

# `ps_trunc.matrix` validates the propensity score matrix from inside a
# dispatched method, so the frame it hands the validator decides what the error
# reports. A frame taken from the method's caller is the user's own function
# when `ps_trunc()` is called from one, and nothing at all when it is called at
# the top level.

condition_call_name <- function(expr) {
  cnd <- rlang::catch_cnd(expr, classes = "error")
  if (is.null(cnd) || is.null(conditionCall(cnd))) {
    return(NA_character_)
  }

  paste(deparse(conditionCall(cnd)[[1]]), collapse = " ")
}

categorical_trunc_fixture <- function() {
  list(
    exposure = factor(c("a", "b", "c", "a")),
    # the first row sums to 1.4, which the matrix method refuses
    bad_matrix = matrix(
      c(0.7, 0.3, 0.2, 0.5, 0.4, 0.4, 0.3, 0.2, 0.3, 0.3, 0.5, 0.3),
      nrow = 4,
      ncol = 3
    )
  )
}

test_that("a refused propensity score matrix names ps_trunc", {
  fixture <- categorical_trunc_fixture()

  cnd <- rlang::catch_cnd(
    ps_trunc(fixture$bad_matrix, method = "ps", .exposure = fixture$exposure),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_matrix_sum_error")

  expect_identical(
    condition_call_name(
      ps_trunc(fixture$bad_matrix, method = "ps", .exposure = fixture$exposure)
    ),
    "ps_trunc"
  )
})

test_that("a refused propensity score matrix names ps_trunc, not the caller", {
  fixture <- categorical_trunc_fixture()

  wrapping_trunc_fn <- function(ps, .exposure) {
    ps_trunc(ps, method = "ps", .exposure = .exposure)
  }

  expect_identical(
    condition_call_name(
      wrapping_trunc_fn(fixture$bad_matrix, fixture$exposure)
    ),
    "ps_trunc"
  )
})

test_that("a matrix with the wrong number of columns names ps_trunc", {
  fixture <- categorical_trunc_fixture()

  wrapping_trunc_fn <- function(ps, .exposure) {
    ps_trunc(ps, method = "ps", .exposure = .exposure)
  }

  expect_identical(
    condition_call_name(
      wrapping_trunc_fn(
        fixture$bad_matrix[, 1:2, drop = FALSE],
        fixture$exposure
      )
    ),
    "ps_trunc"
  )
})

test_that("a refusal on the delegated categorical route names ps_trunc", {
  fixture <- categorical_trunc_fixture()
  valid <- valid_trunc_matrix_fixture()

  # A data frame of categorical scores reaches the matrix method by a plain
  # call rather than by dispatch, so the frame the refusals report against has
  # to travel with them. Whichever guard fires, the caller wrote `ps_trunc()`.
  expect_identical(
    condition_call_name(
      ps_trunc(
        as.data.frame(fixture$bad_matrix),
        method = "ps",
        .exposure = fixture$exposure
      )
    ),
    "ps_trunc"
  )

  expect_identical(
    condition_call_name(
      ps_trunc(
        as.data.frame(valid$ps_matrix),
        .exposure = valid$exposure,
        method = "cr"
      )
    ),
    "ps_trunc"
  )

  expect_identical(
    condition_call_name(
      ps_trunc(
        as.data.frame(valid$ps_matrix),
        .exposure = valid$exposure,
        method = "pctl",
        lower = 1.2
      )
    ),
    "ps_trunc"
  )
})

test_that("ps_trunc refuses a call argument that names no frame", {
  fixture <- valid_trunc_matrix_fixture()

  # The generic passes its dots to its methods, so the frame the categorical
  # path reports against is reachable from user code, and a value the condition
  # system cannot read as one is refused where it arrives rather than left to
  # turn the next guard that fires into a report of rlang's internals.
  expect_error(
    ps_trunc(
      fixture$ps_matrix,
      .exposure = fixture$exposure,
      method = "ps",
      call = "bogus"
    ),
    class = "propensity_call_arg_error"
  )

  # The data frame method hands the frame on, so it is the one that reads the
  # value and the one the refusal names.
  delegated <- rlang::catch_cnd(
    ps_trunc(
      as.data.frame(fixture$ps_matrix),
      .exposure = fixture$exposure,
      method = "ps",
      call = "bogus"
    ),
    classes = "error"
  )

  expect_s3_class(delegated, "propensity_call_arg_error")
  expect_identical(
    paste(deparse(conditionCall(delegated)[[1]]), collapse = " "),
    "ps_trunc"
  )
})

test_that("ps_trunc truncates a matrix with the method left at its default", {
  fixture <- valid_trunc_matrix_fixture()

  # The generic's default is the full three-value vector, which the matrix
  # method narrows to the two methods it supports. Omitting `method` has to
  # resolve to the first of those rather than failing to match the vector it
  # was handed.
  truncated <- ps_trunc(fixture$ps_matrix, .exposure = fixture$exposure)

  # The whole chain, in order: a test that any one of these is inherited holds
  # for a matrix result that lost the class it was given.
  expect_s3_class(
    truncated,
    c("ps_trunc_matrix", "ps_trunc", "matrix"),
    exact = TRUE
  )
  expect_equal(ps_trunc_meta(truncated)$method, "ps")
  expect_equal(dim(truncated), dim(fixture$ps_matrix))
})

test_that("ps_trunc rejects an unsupported matrix method with a package error", {
  fixture <- valid_trunc_matrix_fixture()

  cnd <- rlang::catch_cnd(
    ps_trunc(fixture$ps_matrix, .exposure = fixture$exposure, method = "cr"),
    classes = "error"
  )

  expect_s3_class(cnd, "propensity_method_error")
  # The message has to name the methods the matrix path does support.
  expect_match(conditionMessage(cnd), "pctl", fixed = TRUE)
})

test_that("ps_trunc aborts when the categorical threshold reaches 1/k", {
  fixture <- valid_trunc_matrix_fixture()

  # The sibling of the condition ps_trim() raises for the same threshold,
  # pinned here so the two stay on one class.
  cnd <- rlang::catch_cnd(
    ps_trunc(
      fixture$ps_matrix,
      .exposure = fixture$exposure,
      method = "ps",
      lower = 0.35
    ),
    classes = "error"
  )

  expect_s3_class(cnd, "propensity_range_error")
})

test_that("ps_trunc defaults the categorical threshold to 0.01", {
  fixture <- valid_trunc_matrix_fixture()

  # Truncation and trimming deliberately default to different thresholds: 0.01
  # here against 0.1 in ps_trim(), because truncation only bounds units and
  # trimming drops them.
  explicit <- ps_trunc(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "ps"
  )
  expect_equal(ps_trunc_meta(explicit)$lower_bound, 0.01)

  defaulted <- ps_trunc(fixture$ps_matrix, .exposure = fixture$exposure)
  expect_equal(ps_trunc_meta(defaulted)$lower_bound, 0.01)
})

test_that("ps_trunc truncates a data frame with the method left at its default", {
  fixture <- valid_trunc_matrix_fixture()

  # The data frame method hands `method` to the matrix method untouched, so the
  # generic's default has to survive that hand-off as well.
  truncated <- ps_trunc(
    as.data.frame(fixture$ps_matrix),
    .exposure = fixture$exposure
  )

  expect_s3_class(
    truncated,
    c("ps_trunc_matrix", "ps_trunc", "matrix"),
    exact = TRUE
  )
  expect_equal(ps_trunc_meta(truncated)$method, "ps")
  expect_equal(ps_trunc_meta(truncated)$lower_bound, 0.01)
})

# Missing values ------------------------------------------------------------

# Truncation rewrites every cell of a row through that row's sum, which is
# unknown when one of its scores is missing, so the whole row comes back
# missing. It is still not a row this package pinned to a bound, so the record
# says nothing about it, and it contributes nothing to the distribution the
# bounds are read from.

missing_cell_trunc_fixture <- function() {
  exposure <- factor(c("a", "b", "c", "a"))
  # Row 3 has no score for level "a". Row 4 carries the extreme scores that
  # both methods pin back to a bound.
  ps_matrix <- rbind(
    c(0.60, 0.20, 0.20),
    c(0.20, 0.60, 0.20),
    c(NA, 0.50, 0.50),
    c(0.05, 0.90, 0.05)
  )
  colnames(ps_matrix) <- levels(exposure)

  list(exposure = exposure, ps_matrix = ps_matrix)
}

test_that("ps_trunc() leaves a row with a missing score out of the record", {
  # Regression guard: the fixed-threshold method already reports the row as one
  # it did not truncate, and this is the policy the percentile method is being
  # brought to.
  fixture <- missing_cell_trunc_fixture()

  with_na <- ps_trunc(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "ps",
    lower = 0.1
  )
  without_na <- ps_trunc(
    fixture$ps_matrix[-3, ],
    .exposure = fixture$exposure[-3],
    method = "ps",
    lower = 0.1
  )
  meta <- ps_trunc_meta(with_na)

  expect_true(all(is.na(unclass(with_na)[3, ])))
  expect_equal(meta$truncated_idx, 4)
  expect_equal(meta$n_obs, 4)
  expect_false(is_unit_truncated(with_na)[3])
  expect_equal(
    as.numeric(unclass(with_na)[-3, ]),
    as.numeric(unclass(without_na))
  )
})

test_that("ps_trunc() takes its percentile bounds from the complete rows", {
  fixture <- missing_cell_trunc_fixture()

  # The bounds are quantiles of every cell in the matrix, and `quantile()`
  # refuses a missing value unless it is told to drop it, so the percentile
  # method is an error on any matrix with one. A row this function cannot
  # renormalize contributes nothing to the distribution those bounds come from.
  with_na <- ps_trunc(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "pctl",
    lower = 0.2,
    upper = 0.8
  )
  without_na <- ps_trunc(
    fixture$ps_matrix[-3, ],
    .exposure = fixture$exposure[-3],
    method = "pctl",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trunc_meta(with_na)

  expect_equal(meta$lower_bound, ps_trunc_meta(without_na)$lower_bound)
  expect_equal(meta$upper_bound, ps_trunc_meta(without_na)$upper_bound)
  expect_true(all(is.na(unclass(with_na)[3, ])))
  expect_equal(meta$truncated_idx, 4)
  expect_false(is_unit_truncated(with_na)[3])
  expect_equal(
    as.numeric(unclass(with_na)[-3, ]),
    as.numeric(unclass(without_na))
  )
})

extreme_missing_cell_trunc_fixture <- function() {
  exposure <- factor(c("a", "b", "c", "a"))
  # Row 3 has no score for level "a" and carries observed scores that both
  # methods would otherwise pin back to a bound. Whether a row belongs in the
  # record is decided by a test that reads `any()` across it, and that test
  # answers `TRUE` off those observed scores whatever the missing one does. Row
  # 4 carries the extreme scores of a row this function really does truncate.
  ps_matrix <- rbind(
    c(0.60, 0.20, 0.20),
    c(0.20, 0.60, 0.20),
    c(NA, 0.05, 0.90),
    c(0.05, 0.90, 0.05)
  )
  colnames(ps_matrix) <- levels(exposure)

  list(exposure = exposure, ps_matrix = ps_matrix)
}

test_that("ps_trunc() keeps a missing row with an extreme score out of the record", {
  fixture <- extreme_missing_cell_trunc_fixture()

  truncated <- ps_trunc(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "ps",
    lower = 0.1
  )
  meta <- ps_trunc_meta(truncated)

  # The observed score that would put the row in the record if it were read on
  # its own.
  expect_lt(fixture$ps_matrix[3, 2], meta$lower_bound)

  # The row comes back missing throughout, so nothing in it was pinned to the
  # threshold and the record has nothing to say about it.
  expect_true(all(is.na(unclass(truncated)[3, ])))
  expect_equal(meta$truncated_idx, 4)
  expect_false(is_unit_truncated(truncated)[3])
  expect_true(is_unit_truncated(truncated)[4])
})

test_that("ps_trunc() keeps a missing row with an extreme score out of the percentile record", {
  fixture <- extreme_missing_cell_trunc_fixture()

  truncated <- ps_trunc(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "pctl",
    lower = 0.2,
    upper = 0.8
  )
  meta <- ps_trunc_meta(truncated)

  # The row holds a score below the lower bound and one above the upper bound,
  # so a test read across the row answers `TRUE` off its observed cells alone.
  expect_lt(fixture$ps_matrix[3, 2], meta$lower_bound)
  expect_gt(fixture$ps_matrix[3, 3], meta$upper_bound)

  expect_true(all(is.na(unclass(truncated)[3, ])))
  expect_equal(meta$truncated_idx, 4)
  expect_false(is_unit_truncated(truncated)[3])
  expect_true(is_unit_truncated(truncated)[4])
})

# What the bounds say, and what the categorical path cannot use ----------------

# Four exposure levels, so the threshold a row cannot meet is 1/4 and the
# message has an exact decimal to name.
four_level_trunc_fixture <- function() {
  exposure <- factor(c("a", "b", "c", "d"))
  ps_matrix <- rbind(
    c(0.40, 0.30, 0.20, 0.10),
    c(0.10, 0.40, 0.30, 0.20),
    c(0.20, 0.10, 0.40, 0.30),
    c(0.30, 0.20, 0.10, 0.40)
  )
  colnames(ps_matrix) <- levels(exposure)

  list(exposure = exposure, ps_matrix = ps_matrix)
}

test_that("ps_trunc() names the threshold and the limit it reached", {
  fixture <- four_level_trunc_fixture()

  cnd <- rlang::catch_cnd(
    ps_trunc(
      fixture$ps_matrix,
      .exposure = fixture$exposure,
      method = "ps",
      lower = 0.3
    ),
    classes = "error"
  )

  expect_s3_class(cnd, "propensity_range_error")

  # The sibling of the message ps_trim() gives for the same threshold: the one
  # the caller supplied, and the 1/k the k columns of the matrix impose.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(msg, "0.3", fixed = TRUE)
  expect_match(msg, "0.25", fixed = TRUE)
})

test_that("ps_trunc() validates the percentiles it bounds a matrix at", {
  fixture <- valid_trunc_matrix_fixture()

  # Bounds that cross pin every score to the bound on the far side of the
  # other, which leaves every row uniform and records a truncation nobody asked
  # for. The vector path refuses the same pair.
  crossed <- rlang::catch_cnd(
    ps_trunc(
      fixture$ps_matrix,
      .exposure = fixture$exposure,
      method = "pctl",
      lower = 0.99,
      upper = 0.01
    ),
    classes = "error"
  )
  expect_s3_class(crossed, "propensity_range_error")

  # A probability outside the unit interval reaches `quantile()` as one, which
  # refuses it in base R's words about an argument the caller never wrote.
  outside <- rlang::catch_cnd(
    ps_trunc(
      fixture$ps_matrix,
      .exposure = fixture$exposure,
      method = "pctl",
      lower = -0.1,
      upper = 0.9
    ),
    classes = "error"
  )
  expect_s3_class(outside, "propensity_range_error")
})

test_that("ps_trunc() records a matrix's percentile bounds without their quantile names", {
  fixture <- valid_trunc_matrix_fixture()
  truncated <- ps_trunc(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "pctl",
    lower = 0.1,
    upper = 0.9
  )
  meta <- ps_trunc_meta(truncated)

  # `quantile()` names its result for the probability it was asked for, which
  # says nothing about the bound and reappears wherever the bound is printed or
  # compared.
  expect_null(names(meta$lower_bound))
  expect_null(names(meta$upper_bound))
})

test_that("ps_trunc() refuses focal levels on the categorical path", {
  fixture <- valid_trunc_matrix_fixture()

  # The categorical path reads one column per exposure level and never resolves
  # a focal level, so an argument naming one describes a coding it does not
  # use. A matrix reaches the refusal by dispatch and a data frame by the
  # hand-off the data frame method makes, and the arguments used to be dropped
  # on both routes, which left the caller believing they were honored.
  expect_focal_refusal <- function(scores) {
    trunc_with <- function(...) {
      ps_trunc(scores, .exposure = fixture$exposure, method = "ps", ...)
    }

    focal <- expect_error(
      trunc_with(.focal_level = "a"),
      class = "propensity_unsupported_arg_error"
    )
    expect_match(conditionMessage(focal), "`.focal_level`", fixed = TRUE)

    reference <- expect_error(
      trunc_with(.reference_level = "a"),
      class = "propensity_unsupported_arg_error"
    )
    expect_match(
      conditionMessage(reference),
      "`.reference_level`",
      fixed = TRUE
    )

    # The deprecated spellings stand for the same two arguments, and are refused
    # on the same terms. The deprecation itself is quieted here so that what is
    # read is the refusal.
    withr::local_options(lifecycle_verbosity = "quiet")
    treated <- expect_error(
      trunc_with(.treated = "a"),
      class = "propensity_unsupported_arg_error"
    )
    expect_match(conditionMessage(treated), "`.treated`", fixed = TRUE)

    untreated <- expect_error(
      trunc_with(.untreated = "a"),
      class = "propensity_unsupported_arg_error"
    )
    expect_match(conditionMessage(untreated), "`.untreated`", fixed = TRUE)
  }

  expect_focal_refusal(fixture$ps_matrix)
  expect_focal_refusal(as.data.frame(fixture$ps_matrix))
})
