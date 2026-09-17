# Reordering weights stabilized on a per-observation score --------------------

# A per-observation stabilization score holds one value per unit, in the order
# of the units, like a position record. `[` knows where each unit went and
# reorders the score with the weights. A slice that is not handed its subscript
# cannot, so it drops the score with a warning, as a length change does.

reorder_score_data <- function(n = 300) {
  set.seed(2024)
  x <- rnorm(n)
  a <- 0.5 * x + rnorm(n)
  b <- rbinom(n, 1, plogis(0.4 * x))
  y <- a + b + x + rnorm(n)
  data.frame(x = x, a = a, b = b, y = y)
}

reorder_score <- function(dat) {
  stats::dnorm(dat$a, mean(dat$a), stats::sd(dat$a))
}

continuous_score_psw <- function(dat) {
  withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      lm(a ~ x, data = dat),
      stabilize = TRUE,
      stabilization_score = reorder_score(dat)
    )
  )
}

joint_score_psw <- function(dat) {
  withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(glm(b ~ x, data = dat, family = binomial())),
      wt_ate(
        lm(a ~ b + x, data = dat),
        stabilize = TRUE,
        stabilization_score = reorder_score(dat)
      ),
      exposure_type = c("binary", "continuous")
    )
  )
}

joint_scores <- function(w) {
  joint_wt_meta(w)$stabilization_score
}

expect_scores_follow <- function(out, w, i) {
  expect_identical(stabilization_score(out), stabilization_score(w)[i])
  expected <- joint_scores(w)
  if (!is.null(expected)) {
    expected <- lapply(expected, function(score) {
      if (length(score) > 1) score[i] else score
    })
  }
  expect_identical(joint_scores(out), expected)
}

test_that("the fixtures carry a per-observation score", {
  dat <- reorder_score_data()

  w <- continuous_score_psw(dat)
  expect_identical(stabilization_score(w), reorder_score(dat))

  joint <- joint_score_psw(dat)
  expect_true(any(lengths(joint_scores(joint)) == nrow(dat)))
})

test_that("`[` reorders a per-observation score with the weights", {
  dat <- reorder_score_data()

  for (w in list(continuous_score_psw(dat), joint_score_psw(dat))) {
    i <- c(5:1, 6:length(w))

    out <- expect_silent(w[i])
    expect_scores_follow(out, w, i)
    expect_identical(vec_data(out), vec_data(w)[i])

    shorter <- expect_silent(w[c(3, 1)])
    expect_scores_follow(shorter, w, c(3, 1))
  }
})

test_that("rev() reorders a per-observation score with the weights", {
  dat <- reorder_score_data()

  for (w in list(continuous_score_psw(dat), joint_score_psw(dat))) {
    out <- expect_silent(rev(w))
    expect_scores_follow(out, w, rev(seq_along(w)))
  }
})

test_that("reordering a data frame with order() reorders the score", {
  dat <- reorder_score_data()
  dat$w <- continuous_score_psw(dat)
  dat$j <- joint_score_psw(dat)
  ord <- order(dat$x)

  out <- expect_silent(dat[ord, ])

  expect_scores_follow(out$w, dat$w, ord)
  expect_scores_follow(out$j, dat$j, ord)
})

test_that("an identity subscript keeps the score unchanged", {
  dat <- reorder_score_data()
  w <- continuous_score_psw(dat)

  expect_identical(stabilization_score(w[seq_along(w)]), stabilization_score(w))
  expect_identical(stabilization_score(w[]), stabilization_score(w))
})

test_that("a scalar score is kept by any reordering", {
  w <- psw(c(1, 2, 3), estimand = "ate", stabilized = TRUE)
  attr(w, "stabilization_score") <- 0.4

  expect_identical(stabilization_score(w[3:1]), 0.4)
  expect_identical(stabilization_score(vctrs::vec_slice(w, 3:1)), 0.4)
})

test_that("vec_slice() at the same length drops the score with a warning", {
  dat <- reorder_score_data()

  for (w in list(continuous_score_psw(dat), joint_score_psw(dat))) {
    expect_warning(
      out <- vctrs::vec_slice(w, rev(seq_along(w))),
      class = "propensity_stabilization_score_warning"
    )
    expect_null(stabilization_score(out))
    expect_true(is_stabilized(out))
    expect_false(any(lengths(joint_scores(out)) > 1))
  }
})

test_that("dplyr::arrange() drops the score with a warning", {
  skip_if_not_installed("dplyr")
  dat <- reorder_score_data()
  dat$w <- continuous_score_psw(dat)

  expect_warning(
    out <- dplyr::arrange(dat, x),
    class = "propensity_stabilization_score_warning"
  )
  expect_null(stabilization_score(out$w))
})

test_that("arithmetic keeps a per-observation score in place", {
  dat <- reorder_score_data()
  w <- continuous_score_psw(dat)

  expect_identical(
    stabilization_score(expect_silent(w * 2)),
    reorder_score(dat)
  )
  expect_identical(
    stabilization_score(expect_silent(w * w)),
    reorder_score(dat)
  )
})

test_that("ipw() accepts weights reordered with their data", {
  dat <- reorder_score_data()
  dat$w <- continuous_score_psw(dat)
  sorted <- dat[order(dat$x), ]

  ps_mod <- lm(a ~ x, data = sorted)
  outcome_mod <- lm(y ~ a, data = sorted, weights = w)

  res <- expect_no_error(ipw(ps_mod, outcome_mod))
  expect_s3_class(res, "ipw")

  # The same weights built on the sorted data give the same estimate.
  rebuilt <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      ps_mod,
      stabilize = TRUE,
      stabilization_score = reorder_score(dat)[order(dat$x)]
    )
  )
  expect_equal(as.numeric(sorted$w), as.numeric(rebuilt))
  expected <- ipw(ps_mod, lm(y ~ a, data = sorted, weights = rebuilt))
  expect_equal(res$estimates, expected$estimates)
})

test_that("ipw() accepts joint weights reordered with their data", {
  dat <- reorder_score_data()
  dat$j <- joint_score_psw(dat)
  sorted <- dat[order(dat$x), ]

  models <- joint_wt_models(
    b = glm(b ~ x, data = sorted, family = binomial()),
    a = lm(a ~ b + x, data = sorted)
  )
  outcome_mod <- lm(y ~ a + b, data = sorted, weights = j)

  res <- expect_no_error(ipw(models, outcome_mod))
  expect_s3_class(res, "ipw")
})
