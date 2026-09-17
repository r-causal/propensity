# Scores and records through casts, subassignment, and is.na<- ----------------

score_stabilized_psw <- function() {
  psw(
    c(1, 2, 3),
    estimand = "ate",
    stabilized = TRUE,
    stabilization_score = c(0.5, 0.6, 0.7)
  )
}

muffle_score_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    propensity_stabilization_score_warning = function(cnd) {
      rlang::cnd_muffle(cnd)
    }
  )
}

test_that("a prototype's per-observation score is not cast onto other data", {
  w <- score_stabilized_psw()

  cast <- expect_silent(vctrs::vec_cast(runif(3), vctrs::vec_ptype(w)))
  expect_s3_class(cast, "psw")
  expect_null(stabilization_score(cast))
  expect_true(is_stabilized(cast))

  combined <- muffle_score_warnings(vctrs::vec_c(runif(3), .ptype = w))
  expect_null(stabilization_score(combined))

  # A single value scales every weight and still holds.
  scalar <- psw(c(1, 2), estimand = "ate", stabilized = TRUE)
  attr(scalar, "stabilization_score") <- 0.4
  cast <- vctrs::vec_cast(runif(3), vctrs::vec_ptype(scalar))
  expect_identical(stabilization_score(cast), 0.4)
})

test_that("assigning weights with a different score explains the refusal", {
  w <- score_stabilized_psw()
  y <- w

  expect_error(
    y[1:3] <- rev(w),
    class = "vctrs_error_cast"
  )
  expect_snapshot(error = TRUE, {
    y[1:3] <- rev(w)
  })

  # The values alone carry no score, so they are assigned under the target's.
  y[1:3] <- vec_data(rev(w))
  expect_identical(vec_data(y), c(3, 2, 1))
  expect_identical(stabilization_score(y), c(0.5, 0.6, 0.7))
})

test_that("is.na<- moves a retained score out of the trimming record", {
  set.seed(11)
  n <- 40
  z <- rnorm(n)
  exposure <- rbinom(n, 1, plogis(z))
  fit <- glm(exposure ~ z, family = binomial())
  x <- ps_trim(fitted(fit), lower = 0.2, upper = 0.8)
  meta <- ps_trim_meta(x)
  retained <- meta$keep_idx[[1]]
  trimmed <- meta$trimmed_idx

  expect_silent({
    is.na(x) <- retained
  })

  after <- ps_trim_meta(x)
  expect_identical(after$keep_idx, meta$keep_idx[meta$keep_idx != retained])
  expect_identical(after$trimmed_idx, trimmed)
  expect_identical(after$n_obs, meta$n_obs)
  expect_false(is_unit_trimmed(x)[[retained]])
  expect_true(is.na(x[[retained]]))

  # Marking a trimmed score missing changes neither set.
  expect_silent({
    is.na(x) <- trimmed[[1]]
  })
  expect_identical(ps_trim_meta(x)$trimmed_idx, trimmed)

  refit <- expect_no_error(ps_refit(x, model = fit))
  expect_true(is.na(refit[[retained]]))
  expect_identical(is_unit_trimmed(refit), is_unit_trimmed(x))
  expect_false(anyNA(vec_data(refit)[after$keep_idx]))
})

test_that("is.na<- moves a truncated score out of the truncation record", {
  x <- ps_trunc(c(0.05, 0.3, 0.5, 0.95, 0.4), lower = 0.1, upper = 0.9)
  expect_identical(ps_trunc_meta(x)$truncated_idx, c(1L, 4L))

  expect_silent({
    is.na(x) <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
  })

  expect_identical(ps_trunc_meta(x)$truncated_idx, 4L)
  expect_identical(ps_trunc_meta(x)$n_obs, 5L)
  expect_identical(is_unit_truncated(x), c(FALSE, FALSE, FALSE, TRUE, FALSE))
})
