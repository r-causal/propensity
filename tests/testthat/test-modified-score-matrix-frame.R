# A modified score matrix as a data frame -------------------------------------

# A row of a score matrix is a unit, and so is a row of the data frame it
# becomes, so each column carries the matrix's trimming or truncation record and
# the weight functions read the frame as they read the matrix.

frame_scores_data <- function(n = 150) {
  set.seed(3)
  x <- rnorm(n)
  lp <- cbind(0, 2 * x, -2 * x)
  pr <- exp(lp) / rowSums(exp(lp))
  z <- factor(
    apply(pr, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)),
    levels = c("a", "b", "c")
  )
  dat <- data.frame(x = x, z = z, y = rnorm(n) + as.integer(z))
  fit <- nnet::multinom(z ~ x, data = dat, trace = FALSE)

  list(dat = dat, fit = fit, p = predict(fit, type = "probs"))
}

frame_trimmed <- function(fx) {
  ps_trim(fx$p, method = "ps", lower = 0.05, .exposure = fx$dat$z)
}

frame_truncated <- function(fx) {
  ps_trunc(fx$p, method = "ps", lower = 0.05, .exposure = fx$dat$z)
}

# Every condition an expression signals, by class, and its value.
frame_conditions <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(cnd) {
      seen <<- c(seen, paste("warning", class(cnd)[[1]]))
      invokeRestart("muffleWarning")
    },
    message = function(cnd) {
      seen <<- c(seen, paste("message", conditionMessage(cnd)))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, seen = seen)
}

test_that("the fixtures modify some units", {
  fx <- frame_scores_data()

  expect_true(any(is_unit_trimmed(frame_trimmed(fx))))
  expect_false(all(is_unit_trimmed(frame_trimmed(fx))))
  expect_true(any(is_unit_truncated(frame_truncated(fx))))
})

test_that("as.data.frame() keeps the trimming record on every column", {
  fx <- frame_scores_data()
  tr <- frame_trimmed(fx)

  out <- expect_silent(as.data.frame(tr))

  expect_s3_class(out, "data.frame")
  expect_identical(names(out), colnames(tr))
  expect_identical(rownames(out), rownames(tr))
  for (column in names(out)) {
    expect_s3_class(out[[column]], "ps_trim")
    expect_identical(ps_trim_meta(out[[column]]), ps_trim_meta(tr))
    expect_identical(
      vec_data(out[[column]]),
      unname(unclass(tr)[, column])
    )
    expect_identical(is_unit_trimmed(out[[column]]), is_unit_trimmed(tr))
  }
})

test_that("as.data.frame() keeps the truncation record on every column", {
  fx <- frame_scores_data()
  tc <- frame_truncated(fx)

  out <- expect_silent(as.data.frame(tc))

  expect_identical(names(out), colnames(tc))
  for (column in names(out)) {
    expect_s3_class(out[[column]], "ps_trunc")
    expect_identical(ps_trunc_meta(out[[column]]), ps_trunc_meta(tc))
    expect_identical(is_unit_truncated(out[[column]]), is_unit_truncated(tc))
  }
})

test_that("as.data.frame() honours row names and unnamed matrices", {
  fx <- frame_scores_data()
  tr <- frame_trimmed(fx)
  labels <- paste0("unit", seq_len(nrow(tr)))

  out <- as.data.frame(tr, row.names = labels)
  expect_identical(rownames(out), labels)

  bare <- fx$p
  dimnames(bare) <- NULL
  expect_warning(
    unnamed <- ps_trim(bare, method = "ps", lower = 0.05, .exposure = fx$dat$z),
    class = "propensity_matrix_no_names_warning"
  )
  out <- as.data.frame(unnamed)
  expect_identical(names(out), c("V1", "V2", "V3"))
  expect_s3_class(out$V1, "ps_trim")
})

test_that("tibble::as_tibble() keeps the record on every column", {
  fx <- frame_scores_data()
  tr <- frame_trimmed(fx)
  tc <- frame_truncated(fx)

  out <- expect_silent(tibble::as_tibble(tr))
  expect_s3_class(out, "tbl_df")
  expect_identical(names(out), colnames(tr))
  expect_s3_class(out$a, "ps_trim")
  expect_identical(ps_trim_meta(out$a), ps_trim_meta(tr))

  out <- expect_silent(tibble::as_tibble(tc))
  expect_s3_class(out$b, "ps_trunc")
  expect_identical(ps_trunc_meta(out$b), ps_trunc_meta(tc))
})

test_that("weights from the data frame match weights from the trimmed matrix", {
  fx <- frame_scores_data()
  tr <- frame_trimmed(fx)
  df <- as.data.frame(tr)
  z <- fx$dat$z

  calls <- list(
    ate = function(p) wt_ate(p, .exposure = z),
    att = function(p) wt_att(p, .exposure = z, .focal_level = "a"),
    atu = function(p) wt_atu(p, .exposure = z, .focal_level = "a"),
    atm = function(p) wt_atm(p, .exposure = z),
    ato = function(p) wt_ato(p, .exposure = z),
    entropy = function(p) wt_entropy(p, .exposure = z)
  )

  for (label in names(calls)) {
    from_matrix <- frame_conditions(calls[[label]](tr))
    from_frame <- frame_conditions(calls[[label]](df))

    expect_identical(from_frame$seen, from_matrix$seen, info = label)
    expect_true(
      "warning propensity_no_refit_warning" %in% from_frame$seen,
      info = label
    )
    expect_identical(from_frame$value, from_matrix$value, info = label)
    expect_true(is_ps_trimmed(from_frame$value), info = label)
    expect_match(estimand(from_frame$value), "; trimmed$", info = label)
    expect_identical(
      is_unit_trimmed(from_frame$value),
      is_unit_trimmed(tr),
      info = label
    )
  }
})

test_that("weights from the data frame match weights from the truncated matrix", {
  fx <- frame_scores_data()
  tc <- frame_truncated(fx)
  df <- as.data.frame(tc)
  z <- fx$dat$z

  from_matrix <- frame_conditions(wt_ate(tc, .exposure = z))
  from_frame <- frame_conditions(wt_ate(df, .exposure = z))

  expect_identical(from_frame$seen, from_matrix$seen)
  expect_identical(from_frame$value, from_matrix$value)
  expect_true(is_ps_truncated(from_frame$value))
  expect_identical(
    is_unit_truncated(from_frame$value),
    is_unit_truncated(tc)
  )

  from_matrix <- frame_conditions(wt_ato(tc, .exposure = z))
  from_frame <- frame_conditions(wt_ato(df, .exposure = z))
  expect_identical(from_frame$value, from_matrix$value)
})

test_that("weights from a refit trimmed matrix's data frame match its matrix", {
  fx <- frame_scores_data()
  refit <- ps_refit(frame_trimmed(fx), model = fx$fit)
  df <- as.data.frame(refit)
  z <- fx$dat$z

  expect_true(is_refit(df$a))

  from_matrix <- frame_conditions(wt_ate(refit, .exposure = z))
  from_frame <- frame_conditions(wt_ate(df, .exposure = z))

  expect_identical(from_frame$seen, from_matrix$seen)
  expect_false("warning propensity_no_refit_warning" %in% from_frame$seen)
  expect_identical(from_frame$value, from_matrix$value)
  expect_true(is_refit(from_frame$value))
})

test_that("a subset of the data frame's rows matches the same subset of the matrix", {
  fx <- frame_scores_data()
  tr <- frame_trimmed(fx)
  df <- as.data.frame(tr)
  rows <- c(40:1, 100:120)
  z <- fx$dat$z[rows]

  sub <- df[rows, ]
  expect_identical(ps_trim_meta(sub$a), ps_trim_meta(tr[rows, ]))
  expect_identical(is_unit_trimmed(sub$c), is_unit_trimmed(tr)[rows])

  from_matrix <- frame_conditions(wt_ate(tr[rows, ], .exposure = z))
  from_frame <- frame_conditions(wt_ate(sub, .exposure = z))
  expect_identical(from_frame$value, from_matrix$value)
})

test_that("ipw() refuses weights from the data frame as it refuses the matrix's", {
  fx <- frame_scores_data()
  dat <- fx$dat
  z <- dat$z

  trimmed <- frame_trimmed(fx)
  truncated <- frame_truncated(fx)
  cases <- list(
    trimmed = list(scores = trimmed, class = "propensity_ipw_trimmed_error"),
    truncated = list(
      scores = truncated,
      class = "propensity_ipw_truncated_error"
    )
  )

  for (label in names(cases)) {
    case <- cases[[label]]
    withr::local_options(propensity.quiet = TRUE)
    dat$w_matrix <- suppressWarnings(
      wt_ate(case$scores, .exposure = z),
      classes = "propensity_no_refit_warning"
    )
    dat$w_frame <- suppressWarnings(
      wt_ate(as.data.frame(case$scores), .exposure = z),
      classes = "propensity_no_refit_warning"
    )

    expect_error(
      ipw(fx$fit, lm(y ~ z, data = dat, weights = w_matrix)),
      class = case$class,
      info = label
    )
    expect_error(
      ipw(fx$fit, lm(y ~ z, data = dat, weights = w_frame)),
      class = case$class,
      info = label
    )
  }
})

test_that("a frame whose columns disagree about the record is refused", {
  fx <- frame_scores_data()
  tr <- frame_trimmed(fx)
  df <- as.data.frame(tr)
  df$b <- vctrs::vec_slice(df$b, seq_len(nrow(df)))

  expect_error(
    wt_ate(df, .exposure = fx$dat$z),
    class = "propensity_matrix_type_error"
  )

  mixed <- as.data.frame(tr)
  mixed$c <- vec_data(mixed$c)
  expect_error(
    wt_ate(mixed, .exposure = fx$dat$z),
    class = "propensity_matrix_type_error"
  )
})
