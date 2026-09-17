# Assigning weights into weights with a modification record -------------------

# `[<-` keeps the target's attributes, so weights modified differently would be
# written in under a record that does not describe them. The cast behind `[<-`
# compares the records the way a combine does, by what they say about the
# modification rather than which units they name.

assign_wt_truncated <- function(upper) {
  wt_trunc(psw(c(1, 2, 50, 3, 4, 5, 6, 7, 8, 100)), "count", upper = upper)
}

assign_trimmed <- function(lower = 0.2, upper = 0.8, refit = FALSE) {
  w <- psw(c(NA, 2, 1.6, NA, 3.3), estimand = "ate; trimmed", trimmed = TRUE)
  meta <- list(
    method = "ps",
    lower = lower,
    upper = upper,
    focal_inverted = FALSE,
    keep_idx = c(2L, 3L, 5L),
    trimmed_idx = c(1L, 4L),
    n_obs = 5L
  )
  if (refit) {
    meta$refit <- TRUE
  }
  attr(w, "ps_trim_meta") <- meta
  w
}

assign_truncated <- function(lower = 0.1, upper = 0.9) {
  w <- psw(
    c(10, 2, 1.6, 10, 3.3),
    estimand = "ate; truncated",
    truncated = TRUE
  )
  attr(w, "ps_trunc_meta") <- list(
    method = "ps",
    lower_bound = lower,
    upper_bound = upper,
    truncated_idx = c(1L, 4L),
    n_obs = 5L
  )
  w
}

test_that("assigning weights truncated at another bound is refused", {
  target <- assign_wt_truncated(2)
  value <- assign_wt_truncated(1)
  meta <- attr(target, "psw_trunc_meta")

  cnd <- expect_error(
    {
      target[1:3] <- value[1:3]
    },
    class = "vctrs_error_cast"
  )
  expect_match(conditionMessage(cnd), "different weight truncation bounds")
  expect_identical(attr(target, "psw_trunc_meta"), meta)

  expect_error(
    vctrs::vec_cast(value, target),
    class = "vctrs_error_cast"
  )
})

test_that("assigning weights trimmed differently is refused", {
  target <- assign_trimmed(0.2, 0.8)

  cnd <- expect_error(
    {
      target[2] <- assign_trimmed(0.3, 0.7)[2]
    },
    class = "vctrs_error_cast"
  )
  expect_match(conditionMessage(cnd), "different trimming parameters")

  cnd <- expect_error(
    {
      target[2] <- assign_trimmed(refit = TRUE)[2]
    },
    class = "vctrs_error_cast"
  )
  expect_match(conditionMessage(cnd), "different trimming parameters")
})

test_that("assigning weights from differently truncated scores is refused", {
  target <- assign_truncated(0.1, 0.9)

  cnd <- expect_error(
    {
      target[2:3] <- assign_truncated(0.05, 0.95)[2:3]
    },
    class = "vctrs_error_cast"
  )
  expect_match(conditionMessage(cnd), "different truncation parameters")
})

test_that("assigning values that carry no record is allowed", {
  target <- assign_wt_truncated(2)
  meta <- attr(target, "psw_trunc_meta")

  expect_silent({
    target[1:2] <- c(1.5, 1.5)
  })
  expect_silent({
    target[3] <- vec_data(assign_wt_truncated(1))[3]
  })
  expect_identical(attr(target, "psw_trunc_meta"), meta)

  # Weights flagged as truncated with no record agree with any bound.
  bare <- psw(c(1.2, 1.3), wt_truncated = TRUE)
  expect_silent({
    target[4:5] <- bare
  })
  expect_identical(attr(target, "psw_trunc_meta"), meta)
})

test_that("assigning weights with the same settings keeps the target's record", {
  target <- assign_trimmed()
  meta <- attr(target, "ps_trim_meta")

  # A value with the same settings is allowed wherever its own units were
  # trimmed. The target's record describes the target's positions, which do
  # not move, so it is kept as is: here unit 2 is now missing but the record
  # still lists it as retained, as `is.na<-` leaves a psw's score records.
  value <- assign_trimmed()
  expect_silent({
    target[2] <- value[1]
  })
  expect_true(is.na(target[[2]]))
  expect_identical(attr(target, "ps_trim_meta"), meta)
  expect_identical(
    is_unit_trimmed(target),
    c(TRUE, FALSE, FALSE, TRUE, FALSE)
  )

  target[1] <- value[2]
  expect_identical(vec_data(target)[[1]], 2)
  expect_identical(attr(target, "ps_trim_meta"), meta)
})

test_that("combines still cast their inputs to the prototype", {
  a <- assign_wt_truncated(2)

  out <- expect_silent(vctrs::vec_c(a, a))
  expect_s3_class(out, "psw")
  expect_silent(vctrs::vec_cast(a, vctrs::vec_ptype(a)))
  expect_silent(vctrs::vec_cast(a, out))
})
