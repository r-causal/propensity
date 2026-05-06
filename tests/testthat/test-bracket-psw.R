test_that("[.psw preserves class and metadata for vector subscripts", {
  w <- psw(c(1, 2, 3, 4, 5), estimand = "ate")

  out <- w[c(1, 3, 5)]
  expect_s3_class(out, "psw")
  expect_equal(estimand(out), "ate")
  expect_equal(vec_data(out), c(1, 3, 5))

  out <- w[c(TRUE, FALSE, TRUE, FALSE, TRUE)]
  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(1, 3, 5))

  out <- w[-1]
  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(2, 3, 4, 5))

  out <- w[c(1, NA, 3)]
  expect_s3_class(out, "psw")
  expect_equal(vec_data(out), c(1, NA, 3))

  out <- w[integer()]
  expect_s3_class(out, "psw")
  expect_length(out, 0)
})

test_that("[.psw with matrix subscript falls back to base linear indexing", {
  w <- psw(c(1, 2, 3, 4, 5), estimand = "ate")

  # A 2x2 logical matrix is unrolled column-major to c(T, F, T, F) and then
  # recycled to length(w) = 5, giving c(T, F, T, F, T). Base R selects
  # positions 1, 3, 5. The psw class is intentionally dropped since matrix
  # subscripts have no meaningful vctrs semantics.
  m <- matrix(c(TRUE, FALSE, TRUE, FALSE), nrow = 2)
  out <- w[m]

  expect_false(inherits(out, "psw"))
  expect_type(out, "double")
  expect_equal(out, c(1, 3, 5))
})

test_that("[.psw with no subscript returns the full psw unchanged", {
  w <- psw(c(1, 2, 3), estimand = "ate")
  out <- w[]
  expect_s3_class(out, "psw")
  expect_equal(estimand(out), "ate")
  expect_equal(vec_data(out), c(1, 2, 3))
})

test_that("psw comparison operators dispatch on either side and return logicals silently", {
  w <- psw(c(1, 2, 3), estimand = "ate")

  # All six ops, psw on the LHS.
  expect_no_warning(
    {
      expect_equal(w == 2, c(FALSE, TRUE, FALSE))
      expect_equal(w != 2, c(TRUE, FALSE, TRUE))
      expect_equal(w < 2, c(TRUE, FALSE, FALSE))
      expect_equal(w > 2, c(FALSE, FALSE, TRUE))
      expect_equal(w <= 2, c(TRUE, TRUE, FALSE))
      expect_equal(w >= 2, c(FALSE, TRUE, TRUE))
    },
    class = "propensity_class_downgrade_warning"
  )

  # psw on the RHS still dispatches to the psw method (S3 group dispatch picks
  # up the method from either operand). Behavior should mirror LHS.
  expect_no_warning(
    {
      expect_equal(2 == w, c(FALSE, TRUE, FALSE))
      expect_equal(2 != w, c(TRUE, FALSE, TRUE))
      expect_equal(2 > w, c(TRUE, FALSE, FALSE))
      expect_equal(2 < w, c(FALSE, FALSE, TRUE))
      expect_equal(2 >= w, c(TRUE, TRUE, FALSE))
      expect_equal(2 <= w, c(FALSE, TRUE, TRUE))
    },
    class = "propensity_class_downgrade_warning"
  )

  # psw vs psw: data-only comparison, no metadata check needed.
  w2 <- psw(c(1, 1, 3), estimand = "att")
  expect_no_warning(
    {
      expect_equal(w == w2, c(TRUE, FALSE, TRUE))
      expect_equal(w > w2, c(FALSE, TRUE, FALSE))
    },
    class = "propensity_class_downgrade_warning"
  )
})

test_that("psw comparisons enforce vctrs strict size semantics", {
  # vec_equal()/vec_compare() error on size-mismatched inputs (anything other
  # than equal length, or one side of length 1). The bypass through
  # psw_compare() must preserve that contract via vec_recycle_common(), not
  # silently fall back to base R recycling.
  a <- psw(c(1, 2, 3, 4), estimand = "ate")
  b <- psw(c(1, 2), estimand = "ate")

  expect_error(a == b, class = "vctrs_error_incompatible_size")
  expect_error(a > b, class = "vctrs_error_incompatible_size")
  expect_error(a != b, class = "vctrs_error_incompatible_size")

  # Scalar broadcasting (size-1 RHS) still works and remains silent — this is
  # the path glm.fit() exercises with `weights == 0` / `weights > 0`.
  expect_no_warning(
    {
      expect_equal(a == 2, c(FALSE, TRUE, FALSE, FALSE))
      expect_equal(a > 2, c(FALSE, FALSE, TRUE, TRUE))
    },
    class = "propensity_class_downgrade_warning"
  )
})

test_that("tidy(glm, conf.int = TRUE) works on glms weighted by psw vectors", {
  skip_if_not_installed("broom")

  set.seed(1)
  n <- 200
  d <- data.frame(
    y = rbinom(n, 1, 0.4),
    x = rbinom(n, 1, 0.5),
    z = rnorm(n)
  )
  ps <- glm(x ~ z, data = d, family = binomial())
  d$ps <- predict(ps, type = "response")
  d$w <- wt_att(d$ps, d$x)

  m <- glm(y ~ x, data = d, weights = w, family = quasibinomial())

  # The bug originally errored with `Subscript `i` must be a simple vector,
  # not a matrix.` Once the [.psw matrix-subscript fix is in place,
  # profile.glm() runs to completion. It also evaluates `weights == 0` and
  # `weights > 0` many times, which used to trigger
  # `propensity_class_downgrade_warning` from vec_ptype2.psw.double on every
  # comparison; the comparison-operator methods on psw silence that path.
  expect_no_warning(
    tidied <- broom::tidy(m, exponentiate = TRUE, conf.int = TRUE),
    class = "propensity_class_downgrade_warning"
  )
  expect_true(all(c("conf.low", "conf.high") %in% names(tidied)))
  expect_equal(nrow(tidied), 2)
})
