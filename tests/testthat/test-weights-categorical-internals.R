# The categorical weight route reads one probability per row out of the score
# matrix and validates that matrix before it does. These tests pin the values
# and the refusals so that a change to how either step is written has to keep
# the same answers, and they measure the allocation of both steps.

# A seeded problem with four levels, missing exposures, and a matrix whose
# columns arrive in a different order from `levels(.exposure)`.
categorical_case <- function() {
  set.seed(2024)
  n <- 200
  lvls <- c("a", "b", "c", "d")
  exposure <- factor(sample(lvls, n, replace = TRUE), levels = lvls)
  exposure[c(7, 42, 113)] <- NA

  ps <- matrix(
    runif(n * length(lvls), min = 0.1, max = 0.9),
    nrow = n,
    dimnames = list(NULL, lvls)
  )
  ps <- ps / rowSums(ps)

  list(
    n = n,
    levels = lvls,
    exposure = exposure,
    ps = ps,
    shuffled = ps[, c("c", "a", "d", "b"), drop = FALSE]
  )
}

# w_i = h(e_i) / e_{i, Z_i}, with the observed level's probability read by a
# plain row-column index rather than by an indicator matrix.
hand_categorical_weight <- function(ps, exposure, estimand, focal = NULL) {
  lvls <- levels(exposure)
  ps <- ps[, lvls, drop = FALSE]

  e_actual <- ps[cbind(seq_len(nrow(ps)), match(exposure, lvls))]

  h_e <- switch(
    estimand,
    ate = rep(1, nrow(ps)),
    att = ps[, focal],
    atu = 1 - ps[, focal],
    ato = 1 / rowSums(1 / ps),
    atm = do.call(pmin, lapply(seq_len(ncol(ps)), function(j) ps[, j])),
    entropy = -rowSums(ps * log(ps))
  )

  h_e / e_actual
}

# Each weight function as a call on a score matrix, with the exposure and the
# focal level the fixture uses held fixed.
categorical_weight_calls <- function(exposure, focal = "b") {
  list(
    ate = function(ps) wt_ate(ps, exposure, exposure_type = "categorical"),
    att = function(ps) wt_att(ps, exposure, .focal_level = focal),
    atu = function(ps) wt_atu(ps, exposure, .focal_level = focal),
    atm = function(ps) wt_atm(ps, exposure, exposure_type = "categorical"),
    ato = function(ps) wt_ato(ps, exposure, exposure_type = "categorical"),
    entropy = function(ps) {
      wt_entropy(ps, exposure, exposure_type = "categorical")
    }
  )
}

test_that("categorical weights match a hand computation for every estimand", {
  case <- categorical_case()

  weight_calls <- categorical_weight_calls(case$exposure)

  for (estimand in names(weight_calls)) {
    expected <- hand_categorical_weight(
      case$ps,
      case$exposure,
      estimand,
      focal = "b"
    )

    expect_equal(
      as.numeric(weight_calls[[estimand]](case$ps)),
      expected,
      ignore_attr = TRUE,
      info = estimand
    )

    # Columns are matched to levels by name, so the order they arrive in does
    # not change the answer.
    expect_equal(
      as.numeric(weight_calls[[estimand]](case$shuffled)),
      expected,
      ignore_attr = TRUE,
      info = estimand
    )
  }
})

test_that("units with a missing exposure get a missing weight", {
  case <- categorical_case()
  missing_rows <- c(7, 42, 113)

  for (estimand in c("ate", "atm", "ato", "entropy")) {
    weights <- switch(
      estimand,
      ate = wt_ate(case$ps, case$exposure, exposure_type = "categorical"),
      atm = wt_atm(case$ps, case$exposure, exposure_type = "categorical"),
      ato = wt_ato(case$ps, case$exposure, exposure_type = "categorical"),
      entropy = wt_entropy(
        case$ps,
        case$exposure,
        exposure_type = "categorical"
      )
    )

    expect_true(all(is.na(as.numeric(weights)[missing_rows])), info = estimand)
    expect_false(anyNA(as.numeric(weights)[-missing_rows]), info = estimand)
  }
})

test_that("a missing score anywhere in a row makes that row's weight missing", {
  case <- categorical_case()

  # Row 1 keeps the probability of its own observed level; the missing value
  # sits in a level the unit was not assigned to.
  observed <- match(case$exposure[1], case$levels)
  unobserved <- setdiff(seq_along(case$levels), observed)[1]

  ps <- case$ps
  ps[1, unobserved] <- NA_real_
  expect_false(is.na(ps[1, observed]))

  weights <- as.numeric(
    wt_ate(ps, case$exposure, exposure_type = "categorical")
  )

  # The whole row is missing, not just the tilts that read every column.
  expect_true(is.na(weights[1]))
  expect_false(is.na(weights[2]))
})

test_that("a missing score carries into its own row's weight and no other", {
  case <- categorical_case()

  # The step that carries a missing score into the weight runs on every call,
  # and what it has to produce is pinned here for both matrices it can be given.
  # Where nothing is missing the weights are the hand computation; where one
  # score is missing the row holding it is missing and every other row is what
  # it was. The missing score sits in a level the unit was not assigned to,
  # which is the case a check that read only the observed level would let
  # through.
  observed <- match(case$exposure[[1]], case$levels)
  unobserved <- setdiff(seq_along(case$levels), observed)[[1]]
  ps_missing <- case$ps
  ps_missing[1, unobserved] <- NA_real_

  weight_calls <- categorical_weight_calls(case$exposure)

  for (estimand in names(weight_calls)) {
    expected <- hand_categorical_weight(
      case$ps,
      case$exposure,
      estimand,
      focal = "b"
    )

    expect_equal(
      as.numeric(weight_calls[[estimand]](case$ps)),
      expected,
      ignore_attr = TRUE,
      info = estimand
    )

    with_missing <- as.numeric(weight_calls[[estimand]](ps_missing))
    expect_true(is.na(with_missing[[1]]), info = estimand)
    expect_equal(
      with_missing[-1],
      expected[-1],
      ignore_attr = TRUE,
      info = estimand
    )
  }
})

test_that("check_ps_matrix_range() refuses scores outside the open interval", {
  interval_cases <- list(
    below_zero = matrix(c(-0.1, 0.5, 1.1, 0.5), nrow = 2),
    above_one = matrix(c(0.5, 0.5, 1.5, 0.5), nrow = 2),
    exactly_zero = matrix(c(0, 0.5, 1, 0.5), nrow = 2),
    exactly_one = matrix(c(1, 0.5, 0, 0.5), nrow = 2),
    positive_infinite = matrix(c(Inf, 0.5, 0.5, 0.5), nrow = 2),
    negative_infinite = matrix(c(-Inf, 0.5, 0.5, 0.5), nrow = 2),
    logical_zero = matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2)
  )

  for (case_name in names(interval_cases)) {
    expect_error(
      check_ps_matrix_range(interval_cases[[case_name]]),
      class = "propensity_range_error",
      info = case_name
    )
  }

  expect_true(check_ps_matrix_range(matrix(c(0.2, 0.5, 0.8, 0.5), nrow = 2)))
})

test_that("check_ps_matrix_range() reads past missing values", {
  # Missing values are not scores, so they are neither refused nor allowed to
  # decide the range. A matrix with nothing to read is accepted here and left
  # to the rest of the pipeline.
  expect_true(check_ps_matrix_range(matrix(c(NA, 0.5, 0.5, 0.5), nrow = 2)))
  expect_true(check_ps_matrix_range(matrix(c(NaN, 0.5, 0.5, 0.5), nrow = 2)))
  expect_true(check_ps_matrix_range(matrix(NA_real_, nrow = 2, ncol = 2)))
  expect_true(check_ps_matrix_range(matrix(numeric(0), nrow = 0, ncol = 2)))

  # A missing value does not hide a score that is out of bounds.
  expect_error(
    check_ps_matrix_range(matrix(c(NA, 0.5, 1.5, 0.5), nrow = 2)),
    class = "propensity_range_error"
  )
})

test_that("the range refusal reports the range it read", {
  cnd <- expect_error(
    check_ps_matrix_range(matrix(c(0.2, 0.5, 1.5, 0.5), nrow = 2)),
    class = "propensity_range_error"
  )

  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(cnd)),
    "The range of values in `.propensity` is 0.2 and 1.5",
    fixed = TRUE
  )
})

test_that("rows that do not sum to 1 are refused", {
  expect_error(
    check_ps_matrix_rowsums(matrix(c(0.3, 0.5, 0.3, 0.5), nrow = 2)),
    class = "propensity_matrix_sum_error"
  )

  # A row containing a missing value has no sum to check, so it is passed over.
  expect_true(check_ps_matrix_rowsums(matrix(c(NA, 0.5, 0.3, 0.5), nrow = 2)))
  expect_true(check_ps_matrix_rowsums(matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2)))
})

# The whole of what `check_ps_matrix_rowsums()` promises, case by case, so that
# how it reads the row sums can change while what it accepts and refuses cannot.
# A row is refused when its sum is a number more than the tolerance away from 1,
# and passed over when its sum is not a number at all.
rowsum_cases_accepted <- function() {
  list(
    passing = matrix(c(0.3, 0.5, 0.7, 0.5), nrow = 2),
    missing_only = matrix(c(NA, 0.5, 0.3, 0.5), nrow = 2),
    all_missing = matrix(NA_real_, nrow = 3, ncol = 2),
    # `NaN` is a sum with no value to compare, which is the missing case again.
    not_a_number = matrix(c(NaN, 0.5, 0.3, 0.5), nrow = 2),
    # A row sum this far from 1 is the floating point error of a matrix that was
    # normalized, which is what the tolerance is there to admit.
    inside_tolerance = matrix(c(0.5, 0.5, 0.5 + 9e-7, 0.5), nrow = 2),
    single_row = matrix(c(0.25, 0.75), nrow = 1),
    # A logical matrix sums as the ones and zeros it is.
    logical = matrix(c(TRUE, FALSE, FALSE, TRUE), nrow = 2)
  )
}

rowsum_cases_refused <- function() {
  list(
    one_bad_row = matrix(c(0.3, 0.5, 0.3, 0.5), nrow = 2),
    several_bad_rows = matrix(rep(0.3, 14), nrow = 7),
    missing_with_a_bad_row = matrix(c(NA, 0.3, 0.3, 0.3), nrow = 2),
    # An infinite row sum is a number, and it is not 1.
    positive_infinite = matrix(c(Inf, 0.5, 0.3, 0.5), nrow = 2),
    negative_infinite = matrix(c(-Inf, 0.5, 0.3, 0.5), nrow = 2),
    outside_tolerance = matrix(c(0.5, 0.5, 0.5 + 2e-6, 0.5), nrow = 2)
  )
}

test_that("check_ps_matrix_rowsums() accepts every matrix it accepts today", {
  accepted <- rowsum_cases_accepted()

  for (case_name in names(accepted)) {
    # Silent as well as accepting: reading the extremes of a column of missing
    # values is where a fast path would report an empty minimum of its own.
    expect_silent(result <- check_ps_matrix_rowsums(accepted[[case_name]]))
    expect_true(result, info = case_name)
  }
})

test_that("check_ps_matrix_rowsums() refuses every matrix it refuses today", {
  refused <- rowsum_cases_refused()

  for (case_name in names(refused)) {
    expect_error(
      check_ps_matrix_rowsums(refused[[case_name]]),
      class = "propensity_matrix_sum_error",
      info = case_name
    )
  }
})

test_that("the row sum refusal names the rows it refused", {
  # The rows are named by their position in the matrix the caller handed over,
  # which is what the caller has to go and look at, so a row passed over for
  # holding a missing value still counts toward the position of the ones after
  # it. Reading the extremes of the row sums says only that some row is out of
  # bounds, so the comparison over every row is what this message needs and the
  # refusal is the one place worth building it.
  one_bad <- expect_error(
    check_ps_matrix_rowsums(matrix(c(0.3, 0.5, 0.3, 0.5), nrow = 2)),
    class = "propensity_matrix_sum_error"
  )
  expect_match(conditionMessage(one_bad), "Problem rows: 1", fixed = TRUE)

  after_missing <- expect_error(
    check_ps_matrix_rowsums(matrix(c(NA, 0.3, 0.3, 0.3), nrow = 2)),
    class = "propensity_matrix_sum_error"
  )
  expect_match(conditionMessage(after_missing), "Problem rows: 2", fixed = TRUE)

  # Five rows and an ellipsis, however many there are.
  several <- expect_error(
    check_ps_matrix_rowsums(matrix(rep(0.3, 14), nrow = 7)),
    class = "propensity_matrix_sum_error"
  )
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(several)),
    "Problem rows: 1, 2, 3, 4, and 5 ...",
    fixed = TRUE
  )
})

test_that("non-numeric scores are refused before the range is read", {
  exposure <- factor(c("a", "b", "c"))
  ps <- matrix(
    as.character(c(0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.3, 0.3, 0.3)),
    nrow = 3,
    dimnames = list(NULL, levels(exposure))
  )

  expect_error(
    wt_ate(ps, exposure, exposure_type = "categorical"),
    class = "propensity_error"
  )

  expect_error(
    ps_tilt(ps, "ate"),
    class = "propensity_matrix_type_error"
  )

  # The range check reads the bounds by comparison, which strings answer with
  # an ordering of their own, so it refuses a matrix it cannot read as numbers.
  expect_error(
    check_ps_matrix_range(ps),
    class = "propensity_matrix_type_error"
  )
})

test_that("check_ps_matrix() refuses a matrix that is not made of numbers", {
  exposure <- factor(c("a", "b", "c"))
  ps <- matrix(
    as.character(c(0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.3, 0.3, 0.3)),
    nrow = 3,
    dimnames = list(NULL, levels(exposure))
  )

  # The row sums are read before the range is, and `rowSums()` answers a matrix
  # of strings with base R's own refusal, which names neither the argument the
  # scores arrived in nor this package.
  expect_error(
    check_ps_matrix(ps, exposure),
    class = "propensity_matrix_type_error"
  )

  # The routes that read a matrix of scores against an exposure reach the same
  # check, so each of them reports the type rather than the row sums.
  expect_error(
    ps_trim(ps, .exposure = exposure),
    class = "propensity_matrix_type_error"
  )
  expect_error(
    ps_trunc(ps, .exposure = exposure),
    class = "propensity_matrix_type_error"
  )
  expect_error(
    wt_ate(as.data.frame(ps), exposure, exposure_type = "categorical"),
    class = "propensity_matrix_type_error"
  )
})

# ---- allocation guards ------------------------------------------------------
# These two expectations are the red assertions: they state the allocation the
# categorical route should need, and the current implementations do not meet
# them.

allocation_ratio <- function(fn, ps_matrix) {
  bench_result <- bench::mark(
    fn(),
    iterations = 1,
    check = FALSE,
    filter_gc = FALSE
  )

  as.numeric(bench_result$mem_alloc) / as.numeric(object.size(ps_matrix))
}

allocation_case <- function(n = 1e5, k = 5) {
  set.seed(9)
  lvls <- letters[seq_len(k)]
  exposure <- factor(sample(lvls, n, replace = TRUE), levels = lvls)
  ps <- matrix(
    runif(n * k, min = 0.1, max = 0.9),
    nrow = n,
    dimnames = list(NULL, lvls)
  )

  list(exposure = exposure, ps = ps / rowSums(ps))
}

test_that("calculate_categorical_weights() allocates a few copies of the matrix", {
  skip_on_cran()
  skip_if_not_installed("bench")

  case <- allocation_case()
  ratio <- allocation_ratio(
    function() calculate_categorical_weights(case$ps, case$exposure, "ate"),
    case$ps
  )

  skip_if(is.na(ratio), "memory profiling is unavailable")

  # RED: about 12 copies today, from the n by K indicator matrix and the
  # validation it runs first.
  expect_lt(ratio, 5)
})

test_that("check_ps_matrix_range() reads the matrix without copying it", {
  skip_on_cran()
  skip_if_not_installed("bench")

  case <- allocation_case()
  ratio <- allocation_ratio(function() check_ps_matrix_range(case$ps), case$ps)

  skip_if(is.na(ratio), "memory profiling is unavailable")

  # RED: about 7 copies today, one per temporary in the chained comparison.
  expect_lt(ratio, 1)
})

# ---- the reads each step makes ----------------------------------------------
# The two expectations below are the red assertions of this round: they state
# what the categorical route should read and allocate, and neither is met today.

# How many times a call scans the whole score matrix for missing values.
# `stats::complete.cases()` is that scan, and it is the only one the route
# makes; a call over a matrix holding no missing value has nothing to find.
count_missing_scans <- function(expr) {
  scans <- 0L
  bump <- function() scans <<- scans + 1L

  suppressMessages(trace(
    "complete.cases",
    tracer = bquote(.(bump)()),
    print = FALSE,
    where = asNamespace("stats")
  ))
  withr::defer(
    suppressMessages(untrace("complete.cases", where = asNamespace("stats")))
  )

  force(expr)
  scans
}

test_that("a matrix with no missing score is not scanned for one", {
  case <- allocation_case(n = 1e3)

  # RED: the scan runs once per call today, whether or not the matrix has a
  # missing value to find. Reading whether the matrix holds one at all answers
  # the same question and stops at the first one it finds.
  scans <- count_missing_scans(
    calculate_categorical_weights(case$ps, case$exposure, "ate")
  )
  expect_identical(scans, 0L)

  # The scan is still what settles the weights of a matrix that does hold one.
  ps_missing <- case$ps
  ps_missing[1, 2] <- NA_real_
  weights <- as.numeric(
    calculate_categorical_weights(ps_missing, case$exposure, "ate")
  )
  expect_true(is.na(weights[[1]]))
  expect_false(anyNA(weights[-1]))
})

test_that("check_ps_matrix_rowsums() reads the row sums without copying them", {
  skip_on_cran()
  skip_if_not_installed("bench")

  case <- allocation_case()
  ratio <- allocation_ratio(
    function() check_ps_matrix_rowsums(case$ps),
    case$ps
  )

  skip_if(is.na(ratio), "memory profiling is unavailable")

  # RED: 0.9 copies of the matrix today. The row sums themselves are a fifth of
  # it, and the comparison over every row is the rest; the extremes of the row
  # sums answer the same question and allocate nothing.
  expect_lt(ratio, 0.4)
})

test_that("the categorical weight call allocates a couple of copies", {
  skip_on_cran()
  skip_if_not_installed("bench")

  case <- allocation_case()
  ratio <- allocation_ratio(
    function() calculate_categorical_weights(case$ps, case$exposure, "ate"),
    case$ps
  )

  skip_if(is.na(ratio), "memory profiling is unavailable")

  # RED: 2.4 copies today. The row sum comparison is 0.7 of that and the
  # missing-value scan 0.3, and neither is a read the call has to make.
  expect_lt(ratio, 2)
})
