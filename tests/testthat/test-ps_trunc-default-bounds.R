# A bound the caller supplies to `ps_trunc(method = "ps")` implies its mirror
# on the vector path, and both paths share one default floor of 0.1.

bounds_of <- function(x) {
  meta <- ps_trunc_meta(x)
  c(lower = meta$lower_bound, upper = meta$upper_bound)
}

default_bounds_ps <- c(0.01, 0.04, 0.07, 0.2, 0.5, 0.8, 0.93, 0.96, 0.99)

# Vector path ---------------------------------------------------------------

test_that("a lower bound alone is mirrored to 1 - lower", {
  out <- ps_trunc(default_bounds_ps, method = "ps", lower = 0.05)

  expect_equal(bounds_of(out), c(lower = 0.05, upper = 0.95))
  expect_equal(
    as.numeric(out),
    c(0.05, 0.05, 0.07, 0.2, 0.5, 0.8, 0.93, 0.95, 0.95)
  )
})

test_that("an upper bound alone is mirrored to 1 - upper", {
  out <- ps_trunc(default_bounds_ps, method = "ps", upper = 0.95)

  expect_equal(bounds_of(out), c(lower = 0.05, upper = 0.95))
  expect_equal(
    as.numeric(out),
    c(0.05, 0.05, 0.07, 0.2, 0.5, 0.8, 0.93, 0.95, 0.95)
  )
})

test_that("neither bound supplied keeps the defaults of 0.1 and 0.9", {
  out <- ps_trunc(default_bounds_ps, method = "ps")

  expect_equal(bounds_of(out), c(lower = 0.1, upper = 0.9))
  expect_equal(
    as.numeric(out),
    c(0.1, 0.1, 0.1, 0.2, 0.5, 0.8, 0.9, 0.9, 0.9)
  )
})

test_that("both bounds supplied are used as written, even when asymmetric", {
  out <- ps_trunc(default_bounds_ps, method = "ps", lower = 0.05, upper = 0.8)

  expect_equal(bounds_of(out), c(lower = 0.05, upper = 0.8))
  expect_equal(
    as.numeric(out),
    c(0.05, 0.05, 0.07, 0.2, 0.5, 0.8, 0.8, 0.8, 0.8)
  )
})

test_that("a lone bound whose mirror crosses it is refused", {
  expect_propensity_error(
    ps_trunc(default_bounds_ps, method = "ps", lower = 0.6)
  )
  expect_propensity_error(
    ps_trunc(default_bounds_ps, method = "ps", upper = 0.4)
  )

  # A bound of exactly one half is its own mirror, which leaves no interval.
  cnd <- rlang::catch_cnd(
    ps_trunc(default_bounds_ps, method = "ps", lower = 0.5),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_range_error")

  cnd <- rlang::catch_cnd(
    ps_trunc(default_bounds_ps, method = "ps", upper = 0.5),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_range_error")
})

test_that("a supplied bound must be a single score inside the unit interval", {
  expect_propensity_error(
    ps_trunc(default_bounds_ps, method = "ps", lower = -0.1)
  )
  expect_propensity_error(
    ps_trunc(default_bounds_ps, method = "ps", upper = 1.2)
  )
  expect_propensity_error(
    ps_trunc(default_bounds_ps, method = "ps", lower = c(0.1, 0.2))
  )

  # Both bounds supplied are checked as well, before they are compared.
  expect_error(
    ps_trunc(default_bounds_ps, method = "ps", lower = 0, upper = 0.9),
    class = "propensity_range_error"
  )
  expect_error(
    ps_trunc(default_bounds_ps, method = "ps", lower = 0.1, upper = 1),
    class = "propensity_range_error"
  )
  expect_error(
    ps_trunc(default_bounds_ps, method = "ps", lower = "0.1"),
    class = "propensity_range_error"
  )
  expect_error(
    ps_trunc(default_bounds_ps, method = "ps", upper = c(0.8, 0.9)),
    class = "propensity_length_error"
  )
  expect_error(
    ps_trunc(default_bounds_ps, method = "ps", upper = numeric()),
    class = "propensity_length_error"
  )

  # A missing bound keeps the refusal it has under every method, and names
  # only the bound that is missing.
  cnd <- rlang::catch_cnd(
    ps_trunc(default_bounds_ps, method = "ps", lower = NA_real_),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_missing_value_error")
  expect_no_match(conditionMessage(cnd), "upper", fixed = TRUE)
})

test_that("a lower bound of 1/c gives the bounds the adaptive method computes", {
  set.seed(4211)
  ps <- runif(200, 0.01, 0.99)
  n <- length(ps)
  c_bound <- sqrt(n) * log(n) / 5

  by_hand <- ps_trunc(ps, method = "ps", lower = 1 / c_bound)
  adaptive <- ps_trunc(ps, method = "adaptive")

  expect_equal(bounds_of(by_hand), bounds_of(adaptive), tolerance = 1e-12)
  expect_equal(as.numeric(by_hand), as.numeric(adaptive), tolerance = 1e-12)
  expect_equal(
    ps_trunc_meta(by_hand)$truncated_idx,
    ps_trunc_meta(adaptive)$truncated_idx
  )
})

test_that("the data frame route mirrors a lone bound", {
  ps_df <- data.frame(
    `0` = 1 - default_bounds_ps,
    `1` = default_bounds_ps,
    check.names = FALSE
  )

  withr::local_options(propensity.quiet = TRUE)
  from_lower <- ps_trunc(ps_df, method = "ps", lower = 0.05)
  from_upper <- ps_trunc(ps_df, method = "ps", upper = 0.95)
  from_neither <- ps_trunc(ps_df, method = "ps")

  expect_equal(bounds_of(from_lower), c(lower = 0.05, upper = 0.95))
  expect_equal(bounds_of(from_upper), c(lower = 0.05, upper = 0.95))
  expect_equal(bounds_of(from_neither), c(lower = 0.1, upper = 0.9))
  expect_equal(
    as.numeric(from_lower),
    as.numeric(ps_trunc(default_bounds_ps, method = "ps", lower = 0.05))
  )
})

test_that("the glm route mirrors a lone bound", {
  set.seed(52)
  x <- rnorm(200)
  z <- rbinom(200, 1, plogis(2 * x))
  fit <- glm(z ~ x, family = binomial)

  from_lower <- ps_trunc(fit, method = "ps", lower = 0.05)
  from_upper <- ps_trunc(fit, method = "ps", upper = 0.95)
  from_neither <- ps_trunc(fit, method = "ps")

  expect_equal(bounds_of(from_lower), c(lower = 0.05, upper = 0.95))
  expect_equal(bounds_of(from_upper), c(lower = 0.05, upper = 0.95))
  expect_equal(bounds_of(from_neither), c(lower = 0.1, upper = 0.9))
  expect_equal(
    unname(as.numeric(from_lower)),
    as.numeric(ps_trunc(unname(fitted(fit)), method = "ps", lower = 0.05))
  )
})

test_that("the glm route refuses a lone bound whose mirror crosses it", {
  set.seed(52)
  x <- rnorm(50)
  z <- rbinom(50, 1, plogis(x))
  fit <- glm(z ~ x, family = binomial)

  cnd <- rlang::catch_cnd(
    ps_trunc(fit, method = "ps", lower = 0.7),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_range_error")
})

# Matrix path ---------------------------------------------------------------

many_level_trunc_fixture <- function(k) {
  lvls <- letters[seq_len(k)]
  exposure <- factor(rep(lvls, length.out = 2 * k), levels = lvls)
  ps_matrix <- matrix(1 / k, nrow = 2 * k, ncol = k)
  # One row with a cell below every threshold these tests use.
  ps_matrix[1, ] <- c(0.001, rep((1 - 0.001) / (k - 1), k - 1))
  colnames(ps_matrix) <- lvls

  list(exposure = exposure, ps_matrix = ps_matrix)
}

test_that("the matrix path defaults the threshold to 0.1", {
  exposure <- factor(c("a", "b", "c", "a", "b", "c"))
  ps_matrix <- rbind(
    c(0.05, 0.45, 0.50),
    c(0.30, 0.50, 0.20),
    c(0.20, 0.30, 0.50),
    c(0.40, 0.40, 0.20),
    c(0.25, 0.35, 0.40),
    c(0.34, 0.33, 0.33)
  )
  colnames(ps_matrix) <- levels(exposure)

  out <- ps_trunc(ps_matrix, .exposure = exposure, method = "ps")
  meta <- ps_trunc_meta(out)

  expect_equal(meta$lower_bound, 0.1)
  expect_true(is.na(meta$upper_bound))
  expect_equal(meta$truncated_idx, 1L)

  expected_row <- c(0.1, 0.45, 0.50) / sum(c(0.1, 0.45, 0.50))
  expect_equal(unname(as.matrix(out)[1, ]), expected_row)

  # The same floor reaches a data frame with the method left at its default.
  from_df <- ps_trunc(as.data.frame(ps_matrix), .exposure = exposure)
  expect_equal(ps_trunc_meta(from_df)$lower_bound, 0.1)
})

test_that("the matrix path ignores a lone upper bound", {
  fixture <- many_level_trunc_fixture(3)

  out <- ps_trunc(
    fixture$ps_matrix,
    .exposure = fixture$exposure,
    method = "ps",
    upper = 0.8
  )
  meta <- ps_trunc_meta(out)

  expect_equal(meta$lower_bound, 0.1)
  expect_true(is.na(meta$upper_bound))
})

test_that("the default floor is accepted below ten levels", {
  fixture <- many_level_trunc_fixture(9)

  out <- ps_trunc(fixture$ps_matrix, .exposure = fixture$exposure)

  expect_equal(ps_trunc_meta(out)$lower_bound, 0.1)
  expect_equal(ps_trunc_meta(out)$truncated_idx, 1L)
})

test_that("the default floor meets the 1/k refusal at ten or more levels", {
  ten <- many_level_trunc_fixture(10)
  expect_propensity_error(
    ps_trunc(ten$ps_matrix, .exposure = ten$exposure, method = "ps")
  )

  twelve <- many_level_trunc_fixture(12)
  cnd <- rlang::catch_cnd(
    ps_trunc(twelve$ps_matrix, .exposure = twelve$exposure),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_range_error")

  # An explicit threshold below 1/k is the way through.
  out <- ps_trunc(
    ten$ps_matrix,
    .exposure = ten$exposure,
    method = "ps",
    lower = 0.05
  )
  expect_equal(ps_trunc_meta(out)$lower_bound, 0.05)
})
