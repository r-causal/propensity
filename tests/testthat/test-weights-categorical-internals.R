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

test_that("categorical weights match a hand computation for every estimand", {
  case <- categorical_case()

  weight_calls <- list(
    ate = function(ps) wt_ate(ps, case$exposure, exposure_type = "categorical"),
    att = function(ps) wt_att(ps, case$exposure, .focal_level = "b"),
    atu = function(ps) wt_atu(ps, case$exposure, .focal_level = "b"),
    atm = function(ps) wt_atm(ps, case$exposure, exposure_type = "categorical"),
    ato = function(ps) wt_ato(ps, case$exposure, exposure_type = "categorical"),
    entropy = function(ps) {
      wt_entropy(ps, case$exposure, exposure_type = "categorical")
    }
  )

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
