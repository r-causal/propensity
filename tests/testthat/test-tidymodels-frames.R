# Reading tidymodels prediction frames as propensity scores
#
# The data frame routes of the weight functions and of the propensity score
# modifiers exist for the tibbles `predict()` returns for a fitted parsnip
# model, and those tibbles come in four shapes:
#
#   * a regression fit returns one numeric column, `.pred`, which is a
#     conditional mean and so the continuous route's propensity score;
#   * a classification fit predicted with `type = "prob"` returns one column per
#     level, named `.pred_<level>`: two of them for a binary exposure, K for a
#     categorical one;
#   * a classification fit predicted with the default type returns a single
#     factor column, `.pred_class`, which holds predicted levels rather than
#     probabilities and so is no propensity score at all.
#
# These tests pin what each shape does on each route it can reach. The
# probability frames are compared against the same call on the numeric vector or
# the named matrix underneath them, which is the reading the package documents;
# the `.pred_class` frame is refused, and the refusal names the prediction type
# that would have given scores. A frame of several numeric columns is ambiguous
# on the continuous route, where nothing names the levels its columns could
# belong to, and is refused unless the caller names the column.

# ---- data and fits ----------------------------------------------------------

sim_tidymodels <- function(seed = 2026, n = 250) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  dose <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  bin <- factor(rbinom(n, 1, plogis(0.6 * x1 + 0.2 * x2)))
  eta_b <- -0.2 + 0.7 * x1
  eta_c <- 0.1 - 0.5 * x1 + 0.4 * x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  u <- runif(n)
  lab <- ifelse(
    u < 1 / denom,
    "a",
    ifelse(u < (1 + exp(eta_b)) / denom, "b", "c")
  )
  cat3 <- factor(lab, levels = c("a", "b", "c"))
  data.frame(x1, x2, dose, bin, cat3)
}

fit_parsnip_regression <- function(dat) {
  parsnip::fit(
    parsnip::set_engine(parsnip::linear_reg(), "lm"),
    dose ~ x1 + x2,
    data = dat
  )
}

fit_parsnip_binary <- function(dat) {
  parsnip::fit(
    parsnip::set_engine(parsnip::logistic_reg(), "glm"),
    bin ~ x1 + x2,
    data = dat
  )
}

fit_parsnip_categorical <- function(dat) {
  parsnip::fit(
    parsnip::set_engine(parsnip::multinom_reg(), "nnet"),
    cat3 ~ x1 + x2,
    data = dat
  )
}

# The matrix the probability frame stands for: the same numbers under the level
# names the package matches columns by, which is what a caller who bypassed the
# frame route would pass.
as_level_matrix <- function(probs, levels) {
  ps <- as.matrix(probs)
  colnames(ps) <- levels
  ps
}

# ---- the regression shape on the continuous route ---------------------------

test_that("a regression prediction frame is read as the continuous score", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  preds <- predict(fit_parsnip_regression(dat), new_data = dat)

  expect_named(preds, ".pred")

  expect_equal(
    as.numeric(wt_ate(preds, dat$dose, exposure_type = "continuous")),
    as.numeric(wt_ate(preds$.pred, dat$dose, exposure_type = "continuous"))
  )
})

test_that("a frame of several numeric columns is refused on the continuous route", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  preds <- predict(fit_parsnip_regression(dat), new_data = dat)
  ambiguous <- data.frame(mu1 = preds$.pred, mu2 = preds$.pred + 1)

  err <- expect_error(
    wt_ate(ambiguous, dat$dose, exposure_type = "continuous"),
    class = "propensity_df_ambiguous_column_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, ".propensity_col", fixed = TRUE)
})

test_that("naming the column reads the ambiguous frame on the continuous route", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  preds <- predict(fit_parsnip_regression(dat), new_data = dat)
  ambiguous <- data.frame(mu1 = preds$.pred, mu2 = preds$.pred + 1)

  expect_equal(
    as.numeric(wt_ate(
      ambiguous,
      dat$dose,
      exposure_type = "continuous",
      .propensity_col = mu1
    )),
    as.numeric(wt_ate(
      ambiguous$mu1,
      dat$dose,
      exposure_type = "continuous"
    ))
  )
})

# ---- the two-column probability shape on the binary route -------------------

test_that("a two-column probability frame is read as the binary score", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  probs <- predict(fit_parsnip_binary(dat), new_data = dat, type = "prob")

  expect_named(probs, c(".pred_0", ".pred_1"))

  expect_equal(
    as.numeric(wt_ate(probs, dat$bin)),
    as.numeric(wt_ate(probs$.pred_1, dat$bin))
  )
})

test_that("the binary probability column is resolved by the level it is named for", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  probs <- predict(fit_parsnip_binary(dat), new_data = dat, type = "prob")

  named <- expect_no_warning(wt_att(probs, dat$bin, .focal_level = "0"))

  expect_equal(
    as.numeric(named),
    as.numeric(wt_att(probs$.pred_0, dat$bin, .focal_level = "0"))
  )
})

test_that("a reordered binary probability frame is read by name, not by position", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  probs <- predict(fit_parsnip_binary(dat), new_data = dat, type = "prob")
  reversed <- probs[, c(".pred_1", ".pred_0")]

  expect_equal(
    as.numeric(wt_ate(reversed, dat$bin)),
    as.numeric(wt_ate(probs$.pred_1, dat$bin))
  )
})

# ---- the K-column probability shape on the categorical route ----------------

test_that("a K-column probability frame is read as the categorical score", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  probs <- predict(fit_parsnip_categorical(dat), new_data = dat, type = "prob")
  ps <- as_level_matrix(probs, levels(dat$cat3))

  expect_named(probs, c(".pred_a", ".pred_b", ".pred_c"))

  expect_equal(
    as.numeric(wt_ate(probs, dat$cat3, exposure_type = "categorical")),
    as.numeric(wt_ate(ps, dat$cat3, exposure_type = "categorical"))
  )
  expect_equal(
    as.numeric(wt_att(
      probs,
      dat$cat3,
      exposure_type = "categorical",
      .focal_level = "b"
    )),
    as.numeric(wt_att(
      ps,
      dat$cat3,
      exposure_type = "categorical",
      .focal_level = "b"
    ))
  )
})

test_that("a reordered categorical probability frame is read by name", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  probs <- predict(fit_parsnip_categorical(dat), new_data = dat, type = "prob")
  ps <- as_level_matrix(probs, levels(dat$cat3))
  shuffled <- probs[, c(".pred_c", ".pred_a", ".pred_b")]

  expect_equal(
    as.numeric(wt_ate(shuffled, dat$cat3, exposure_type = "categorical")),
    as.numeric(wt_ate(ps, dat$cat3, exposure_type = "categorical"))
  )
})

# ---- the predicted-class shape ----------------------------------------------

test_that("a predicted-class frame is refused on the binary route", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  classes <- predict(fit_parsnip_binary(dat), new_data = dat)

  expect_named(classes, ".pred_class")

  err <- expect_error(
    wt_ate(classes, dat$bin),
    class = "propensity_df_class_column_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, ".pred_class", fixed = TRUE)
  expect_match(msg, "prob", fixed = TRUE)
})

test_that("a predicted-class frame is refused on the categorical route", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  classes <- predict(fit_parsnip_categorical(dat), new_data = dat)

  err <- expect_error(
    wt_ate(classes, dat$cat3, exposure_type = "categorical"),
    class = "propensity_df_class_column_error"
  )
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(err)),
    "prob",
    fixed = TRUE
  )
})

test_that("a predicted-class frame is refused on the continuous route", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  classes <- predict(fit_parsnip_binary(dat), new_data = dat)

  expect_error(
    wt_ate(classes, dat$dose, exposure_type = "continuous"),
    class = "propensity_df_class_column_error"
  )
})

# ---- the modifiers ----------------------------------------------------------

test_that("ps_trim() and ps_trunc() read a binary probability frame", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  probs <- predict(fit_parsnip_binary(dat), new_data = dat, type = "prob")

  expect_equal(
    as.numeric(ps_trim(probs, lower = 0.1, upper = 0.9)),
    as.numeric(ps_trim(probs$.pred_1, lower = 0.1, upper = 0.9))
  )
  expect_equal(
    as.numeric(ps_trunc(probs, lower = 0.1, upper = 0.9)),
    as.numeric(ps_trunc(probs$.pred_1, lower = 0.1, upper = 0.9))
  )
})

test_that("the modifiers read a categorical probability frame", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  probs <- predict(fit_parsnip_categorical(dat), new_data = dat, type = "prob")
  ps <- as_level_matrix(probs, levels(dat$cat3))

  expect_equal(
    as.numeric(ps_trim(probs, .exposure = dat$cat3, method = "optimal")),
    as.numeric(ps_trim(ps, .exposure = dat$cat3, method = "optimal"))
  )
  expect_equal(
    as.numeric(ps_tilt(probs, "ato")),
    as.numeric(ps_tilt(ps, "ato"))
  )
  expect_equal(
    as.numeric(ps_tilt(probs, "att", .focal_level = "a")),
    as.numeric(ps_tilt(ps, "att", .focal_level = "a"))
  )
})

test_that("ps_calibrate() reads a binary probability frame", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  probs <- predict(fit_parsnip_binary(dat), new_data = dat, type = "prob")

  expect_equal(
    as.numeric(ps_calibrate(probs, dat$bin, smooth = FALSE)),
    as.numeric(ps_calibrate(probs$.pred_1, dat$bin, smooth = FALSE))
  )
})

test_that("the modifiers refuse a predicted-class frame", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  classes <- predict(fit_parsnip_binary(dat), new_data = dat)
  cat_classes <- predict(fit_parsnip_categorical(dat), new_data = dat)

  expect_error(
    ps_trim(classes, lower = 0.1),
    class = "propensity_df_class_column_error"
  )
  expect_error(
    ps_trunc(classes, lower = 0.1),
    class = "propensity_df_class_column_error"
  )
  expect_error(
    ps_tilt(cat_classes, "ato"),
    class = "propensity_df_class_column_error"
  )
  expect_error(
    ps_calibrate(classes, dat$bin, smooth = FALSE),
    class = "propensity_df_class_column_error"
  )
})

# ---- the wording of the new refusals ----------------------------------------

test_that("the frame refusals name the column to set and the type to predict", {
  skip_if_not_installed("parsnip")
  dat <- sim_tidymodels()
  preds <- predict(fit_parsnip_regression(dat), new_data = dat)
  ambiguous <- data.frame(mu1 = preds$.pred, mu2 = preds$.pred + 1)
  classes <- predict(fit_parsnip_binary(dat), new_data = dat)

  expect_propensity_error(
    wt_ate(ambiguous, dat$dose, exposure_type = "continuous")
  )
  expect_propensity_error(wt_ate(classes, dat$bin))
  expect_propensity_error(ps_trim(classes, lower = 0.1))
  expect_propensity_error(ps_calibrate(classes, dat$bin, smooth = FALSE))
})

test_that("a frame mixing a score column with a label column is still ambiguous", {
  # The refusal counts columns rather than candidates: a frame with one numeric
  # column beside a character one holds two things the caller might have meant
  # to select from, and reading the numeric one would be a choice made on the
  # types rather than on anything the caller expressed.
  dat <- sim_tidymodels()
  mixed <- data.frame(
    mu = dat$dose - dat$x1,
    label = as.character(dat$bin),
    stringsAsFactors = FALSE
  )

  expect_error(
    wt_ate(mixed, dat$dose, exposure_type = "continuous"),
    class = "propensity_df_ambiguous_column_error"
  )
  expect_equal(
    as.numeric(wt_ate(
      mixed,
      dat$dose,
      exposure_type = "continuous",
      .propensity_col = mu
    )),
    as.numeric(wt_ate(mixed$mu, dat$dose, exposure_type = "continuous"))
  )
})

# ---- levels that are themselves named like prediction columns ---------------

# The `.pred_` prefix is a convention of the frames parsnip returns, not a
# reserved string: an exposure is allowed levels that begin with it, and then
# the columns naming those levels begin with it twice over as far as the prefix
# rule is concerned. Reading such a frame under the stripped names matches the
# columns against halves of the level names, which either selects the wrong
# column or reports level names the caller never used. The exact names answer
# first, so a frame whose columns carry the levels verbatim is read as such.

sim_pred_like_levels <- function(seed = 2027, n = 200) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  p <- plogis(0.5 * x1)
  bin <- factor(
    ifelse(rbinom(n, 1, p) == 1, ".pred_b", ".pred_a"),
    levels = c(".pred_a", ".pred_b")
  )
  cat3 <- factor(
    sample(c(".pred_a", ".pred_b", ".pred_c"), n, replace = TRUE),
    levels = c(".pred_a", ".pred_b", ".pred_c")
  )
  list(x1 = x1, p = p, bin = bin, cat3 = cat3)
}

test_that("a binary frame whose columns name pred-like levels is read by name", {
  dat <- sim_pred_like_levels()
  probs <- data.frame(.pred_a = 1 - dat$p, .pred_b = dat$p)

  expect_equal(
    as.numeric(wt_att(probs, dat$bin)),
    as.numeric(wt_att(dat$p, dat$bin))
  )
  expect_equal(
    as.numeric(wt_att(probs, dat$bin, .focal_level = ".pred_a")),
    as.numeric(wt_att(1 - dat$p, dat$bin, .focal_level = ".pred_a"))
  )
})

test_that("a reordered binary frame of pred-like levels is read by name", {
  dat <- sim_pred_like_levels()
  probs <- data.frame(.pred_a = 1 - dat$p, .pred_b = dat$p)
  reversed <- probs[, c(".pred_b", ".pred_a")]

  expect_equal(
    as.numeric(wt_att(reversed, dat$bin)),
    as.numeric(wt_att(dat$p, dat$bin))
  )
  expect_equal(
    as.numeric(wt_att(reversed, dat$bin, .focal_level = ".pred_a")),
    as.numeric(wt_att(1 - dat$p, dat$bin, .focal_level = ".pred_a"))
  )
})

test_that("a categorical matrix of pred-like level names is read by name", {
  dat <- sim_pred_like_levels()
  lv <- levels(dat$cat3)
  ps <- matrix(runif(length(dat$cat3) * 3), ncol = 3)
  ps <- ps / rowSums(ps)
  colnames(ps) <- lv
  plain <- ps
  colnames(plain) <- sub("^\\.pred_", "", lv)
  plain_exposure <- factor(
    sub("^\\.pred_", "", as.character(dat$cat3)),
    levels = sub("^\\.pred_", "", lv)
  )

  expect_equal(
    as.numeric(wt_ate(ps, dat$cat3, exposure_type = "categorical")),
    as.numeric(wt_ate(plain, plain_exposure, exposure_type = "categorical"))
  )

  shuffled <- ps[, c(3L, 1L, 2L), drop = FALSE]
  expect_equal(
    as.numeric(wt_ate(shuffled, dat$cat3, exposure_type = "categorical")),
    as.numeric(wt_ate(plain, plain_exposure, exposure_type = "categorical"))
  )
})

test_that("a categorical frame of pred-like level names is read by name", {
  dat <- sim_pred_like_levels()
  lv <- levels(dat$cat3)
  ps <- matrix(runif(length(dat$cat3) * 3), ncol = 3)
  ps <- ps / rowSums(ps)
  colnames(ps) <- lv

  expect_equal(
    as.numeric(wt_att(
      as.data.frame(ps),
      dat$cat3,
      exposure_type = "categorical",
      .focal_level = ".pred_b"
    )),
    as.numeric(wt_att(
      ps,
      dat$cat3,
      exposure_type = "categorical",
      .focal_level = ".pred_b"
    ))
  )
})
