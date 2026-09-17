# What combined weights keep of their modification records --------------------

# Weights built from differently trimmed or truncated scores, or truncated at
# different bounds, target different estimands, so they have no common type.
# The combine compares what each record says about the modification, its method,
# bounds, and refit flag, and keeps that part on the result, so a later combine,
# `is_refit()`, and the footer can still read it. The positions are dropped,
# since nothing tells the combine where each input's units land.

combine_wt_truncated <- function(upper) {
  wt_trunc(psw(c(1, 2, 50, 3, 4, 5, 6, 7, 8, 100)), "count", upper = upper)
}

combine_trimmed <- function(lower = 0.2, upper = 0.8, refit = FALSE) {
  w <- psw(
    c(NA, 2, 1.6, NA, 3.3),
    estimand = "ate; trimmed",
    trimmed = TRUE
  )
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

combine_truncated <- function(lower = 0.1, upper = 0.9) {
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

# Every order of a list of inputs.
input_orders <- function(inputs) {
  orders <- list(integer())
  for (k in seq_along(inputs)) {
    orders <- unlist(
      lapply(orders, function(o) {
        lapply(setdiff(seq_along(inputs), o), function(j) c(o, j))
      }),
      recursive = FALSE
    )
  }
  lapply(orders, function(o) inputs[o])
}

# A nested combine can downgrade twice, once inside and once when the numeric
# result meets the remaining weights, so every downgrade warning is collected
# and at least one of them must name the incompatible records.
expect_combine_downgrades <- function(expr, label = NULL) {
  classes <- character()
  out <- withCallingHandlers(
    expr,
    propensity_warning = function(cnd) {
      classes <<- c(classes, class(cnd))
      rlang::cnd_muffle(cnd)
    }
  )
  expect_true("propensity_coercion_warning" %in% classes, info = label)
  expect_false(is_psw(out), info = label)
  expect_type(out, "double")
}

# Weight truncation bounds --------------------------------------------------------

test_that("nested combines compare weight truncation bounds", {
  a <- combine_wt_truncated(2)
  b <- combine_wt_truncated(1)

  expect_combine_downgrades(vctrs::vec_c(vctrs::vec_c(a, a), b))
  expect_combine_downgrades(vctrs::vec_c(b, vctrs::vec_c(a, a)))
  expect_combine_downgrades(c(c(a, a), b))
  expect_combine_downgrades(c(b, c(a, a)))
})

test_that("nested bind_rows() compares weight truncation bounds", {
  skip_if_not_installed("dplyr")
  da <- data.frame(id = 1:10)
  da$w <- combine_wt_truncated(2)
  db <- data.frame(id = 1:10)
  db$w <- combine_wt_truncated(1)

  expect_warning(
    out <- dplyr::bind_rows(dplyr::bind_rows(da, da), db),
    class = "propensity_coercion_warning"
  )
  expect_false(is_psw(out$w))

  expect_warning(
    out <- dplyr::bind_rows(db, dplyr::bind_rows(da, da)),
    class = "propensity_coercion_warning"
  )
  expect_false(is_psw(out$w))
})

test_that("agreeing weight truncation records survive a combine without positions", {
  a <- combine_wt_truncated(2)

  out <- expect_silent(vctrs::vec_c(a, a))
  expect_s3_class(out, "psw")
  expect_positions_dropped(
    attr(out, "psw_trunc_meta"),
    attr(a, "psw_trunc_meta")
  )
  expect_error(
    is_unit_wt_truncated(out),
    class = "propensity_missing_meta_error"
  )
  expect_true(any(
    utils::capture.output(print(out)) == "truncation: count 2 (upper 8)"
  ))

  again <- expect_silent(c(out, a))
  expect_positions_dropped(
    attr(again, "psw_trunc_meta"),
    attr(a, "psw_trunc_meta")
  )
})

# Trimming records ----------------------------------------------------------------

test_that("weights trimmed at different cutoffs do not combine", {
  t1 <- combine_trimmed(0.2, 0.8)
  t2 <- combine_trimmed(0.3, 0.7)

  expect_combine_downgrades(vctrs::vec_c(t1, t2))
  expect_combine_downgrades(vctrs::vec_c(t2, t1))
  expect_combine_downgrades(c(t1, t2))
  expect_combine_downgrades(vctrs::vec_c(vctrs::vec_c(t1, t1), t2))
  expect_combine_downgrades(vctrs::vec_c(t2, vctrs::vec_c(t1, t1)))
  expect_combine_downgrades(c(c(t1, t1), t2))
})

test_that("weights from a refit and an unrefit trim do not combine", {
  refit <- combine_trimmed(refit = TRUE)
  plain <- combine_trimmed(refit = FALSE)

  expect_combine_downgrades(vctrs::vec_c(refit, plain))
  expect_combine_downgrades(vctrs::vec_c(plain, refit))
  expect_combine_downgrades(vctrs::vec_c(vctrs::vec_c(refit, refit), plain))
  expect_combine_downgrades(c(plain, c(refit, refit)))
})

test_that("agreeing trimming records survive a combine without positions", {
  refit <- combine_trimmed(refit = TRUE)

  out <- expect_silent(vctrs::vec_c(refit, refit))
  expect_s3_class(out, "psw")
  expect_positions_dropped(ps_trim_meta(out), ps_trim_meta(refit))
  expect_true(is_refit(out))
  expect_error(is_unit_trimmed(out), class = "propensity_missing_meta_error")

  plain <- expect_silent(c(combine_trimmed(), combine_trimmed()))
  expect_false(is_refit(plain))
})

test_that("is_refit() answers for combined weights from a refit model", {
  set.seed(42)
  n <- 60
  z <- rnorm(n)
  x <- rbinom(n, 1, plogis(2 * z))
  fit <- glm(x ~ z, family = binomial())
  trimmed <- ps_trim(
    predict(fit, type = "response"),
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )
  w2 <- wt_ate(
    ps_refit(trimmed, model = fit),
    .exposure = x,
    exposure_type = "binary",
    .focal_level = 1
  )

  expect_true(is_refit(vctrs::vec_c(w2, w2)))
})

# Truncation records of the scores ------------------------------------------------

test_that("weights from scores truncated at different bounds do not combine", {
  u1 <- combine_truncated(0.1, 0.9)
  u2 <- combine_truncated(0.05, 0.95)

  expect_combine_downgrades(vctrs::vec_c(u1, u2))
  expect_combine_downgrades(vctrs::vec_c(vctrs::vec_c(u1, u1), u2))
  expect_combine_downgrades(c(u2, c(u1, u1)))

  out <- expect_silent(vctrs::vec_c(u1, u1))
  expect_positions_dropped(ps_trunc_meta(out), ps_trunc_meta(u1))
})

# Folds of three ------------------------------------------------------------------

test_that("a disagreement is caught in every order of a three-input fold", {
  sets <- list(
    wt = list(
      combine_wt_truncated(2),
      combine_wt_truncated(2),
      combine_wt_truncated(1)
    ),
    trim = list(
      combine_trimmed(0.2, 0.8),
      combine_trimmed(0.2, 0.8),
      combine_trimmed(0.3, 0.7)
    ),
    refit = list(
      combine_trimmed(refit = TRUE),
      combine_trimmed(refit = TRUE),
      combine_trimmed(refit = FALSE)
    ),
    trunc = list(
      combine_truncated(0.1, 0.9),
      combine_truncated(0.1, 0.9),
      combine_truncated(0.05, 0.95)
    )
  )

  for (label in names(sets)) {
    for (inputs in input_orders(sets[[label]])) {
      expect_combine_downgrades(vctrs::vec_c(!!!inputs), label)
      expect_combine_downgrades(
        vctrs::vec_c(vctrs::vec_c(!!!inputs[1:2]), inputs[[3]]),
        label
      )
      expect_combine_downgrades(
        vctrs::vec_c(inputs[[1]], vctrs::vec_c(!!!inputs[2:3])),
        label
      )
    }
  }
})

test_that("agreeing records survive every order of a three-input fold", {
  a <- combine_wt_truncated(2)
  # Flagged as truncated with no record, which agrees with any bound.
  bare <- psw(c(1, 2), wt_truncated = TRUE)

  for (inputs in input_orders(list(a, a, bare))) {
    out <- expect_silent(vctrs::vec_c(!!!inputs))
    expect_s3_class(out, "psw")
    expect_positions_dropped(
      attr(out, "psw_trunc_meta"),
      attr(a, "psw_trunc_meta")
    )

    nested <- expect_silent(
      vctrs::vec_c(vctrs::vec_c(!!!inputs[1:2]), inputs[[3]])
    )
    expect_positions_dropped(
      attr(nested, "psw_trunc_meta"),
      attr(a, "psw_trunc_meta")
    )
  }
})

test_that("a record without settings to compare still conflicts through a bare input", {
  a <- combine_wt_truncated(2)
  b <- combine_wt_truncated(1)
  bare <- psw(c(1, 2), wt_truncated = TRUE)

  for (inputs in input_orders(list(a, bare, b))) {
    expect_combine_downgrades(vctrs::vec_c(!!!inputs))
  }
})

test_that("a flagged psw with no trimming record takes the other record", {
  t1 <- combine_trimmed(refit = TRUE)
  bare <- psw(c(1, 2), estimand = "ate; trimmed", trimmed = TRUE)

  out <- expect_silent(vctrs::vec_c(bare, t1))
  expect_positions_dropped(ps_trim_meta(out), ps_trim_meta(t1))
  expect_true(is_refit(out))
})
