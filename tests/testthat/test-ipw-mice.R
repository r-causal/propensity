# Pooling `ipw()` results across multiply imputed datasets with mice.
#
# `mice::pool()` is the end of the standard multiple-imputation workflow: fit
# the analysis once per completed dataset, then combine the estimates by Rubin's
# rules. It reaches an analysis object through `tidy()`, so an `ipw` result pools
# by being tidied, and these tests run the whole path rather than the tidier
# alone: impute, fit once per imputation through `with()`, pool.
#
# What pooling needs from the table is a row per estimate, an `estimate` and a
# `std.error` column, and a column naming the rows so that the same estimate
# from different imputations is combined with itself rather than with its
# neighbour. mice keys rows on `term` and, where a table has one, `contrast`,
# both of which the tidier reports. A categorical exposure is where that matters:
# its rows repeat the same three effect measures once per contrast, so `term`
# alone would collapse six estimates into three groups of twice the imputations.
#
# The degrees of freedom mice reports are Barnard-Rubin, which needs the
# complete-data residual degrees of freedom. It takes those from
# `df.residual()` on the first fit. M-estimation records them, so the pooled
# degrees of freedom come back finite. Linearization stacks no system and
# records none, so they come back missing, and the remedy is `pool(dfcom = )`.
# Both are pinned below.

# ---- fixtures ---------------------------------------------------------------
#
# Small and seeded throughout: the point is the shape of the pooled table and
# the path that produces it, not the precision of any estimate, and each test
# imputes a dataset and fits the analysis once per imputation.

# A binary exposure with a covariate missing at random, the missingness
# depending on the exposure, which is observed. `mice()` is left on its defaults
# apart from the imputation count and the seed.
mice_binary_data <- function(seed = 2024, n = 250) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.5 * x2))
  y <- rbinom(n, 1, plogis(-0.4 + 0.9 * z + 0.5 * x1))
  dat <- data.frame(x1, x2, z, y)
  dat$x1[rbinom(n, 1, plogis(-1.2 + 0.8 * z)) == 1] <- NA
  dat
}

mice_categorical_data <- function(seed = 2025, n = 300) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  eta_b <- -0.2 + 0.5 * x1 + 0.3 * x2
  eta_c <- 0.1 - 0.4 * x1 + 0.6 * x2
  den <- 1 + exp(eta_b) + exp(eta_c)
  u <- runif(n)
  pa <- 1 / den
  pb <- exp(eta_b) / den
  a <- factor(
    ifelse(u < pa, "a", ifelse(u < pa + pb, "b", "c")),
    levels = c("a", "b", "c")
  )
  y <- rbinom(
    n,
    1,
    plogis(-0.3 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.5 * x1)
  )
  dat <- data.frame(x1, x2, a, y)
  dat$x1[rbinom(n, 1, plogis(-1.3 + 0.6 * x2)) == 1] <- NA
  dat
}

mice_continuous_data <- function(seed = 2026, n = 250) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  A <- 0.5 + 0.6 * x1 - 0.3 * x2 + rnorm(n)
  yc <- 1 + 0.4 * A + 0.5 * x1 + rnorm(n)
  dat <- data.frame(x1, x2, A, yc)
  dat$x1[rbinom(n, 1, plogis(-1.3 + 0.5 * x2)) == 1] <- NA
  dat
}

# `mice()` sets the global seed and does not put it back, so the call is wrapped
# rather than left to leak its state into whatever runs next, matching how the
# data fixtures above seed themselves.
mice_impute <- function(dat, seed) {
  withr::with_preserve_seed(mice::mice(dat, m = 3, print = FALSE, seed = seed))
}

# The analysis, once per completed dataset. Everything inside `with()` resolves
# against the completed data, so the propensity model, the weights, and the
# outcome model are all rebuilt on each imputation rather than carried in from
# outside it, which is what makes the pooled variance account for the imputation.
fit_mice_binary <- function(imp, se_method = "mestimation") {
  with(imp, {
    ps <- glm(z ~ x1 + x2, family = binomial())
    w <- wt_ate(ps)
    om <- glm(y ~ z, family = quasibinomial(), weights = w)
    ipw(ps, om, se_method = se_method)
  })
}

fit_mice_categorical <- function(imp) {
  with(imp, {
    ps <- nnet::multinom(a ~ x1 + x2, trace = FALSE)
    w <- wt_ate(predict(ps, type = "probs"), a, exposure_type = "categorical")
    om <- glm(y ~ a, family = quasibinomial(), weights = w)
    ipw(ps, om)
  })
}

fit_mice_continuous <- function(imp) {
  with(imp, {
    ps <- lm(A ~ x1 + x2)
    w <- wt_ate(
      as.double(fitted(ps)),
      A,
      exposure_type = "continuous",
      stabilize = TRUE
    )
    om <- lm(yc ~ A, weights = w)
    ipw(ps, om)
  })
}

# ---- binary, m-estimation ---------------------------------------------------

test_that("mice::pool() combines binary mestimation ipw results", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_binary_data(), seed = 1)
  fits <- fit_mice_binary(imp)

  expect_s3_class(fits, "mira")
  expect_length(fits$analyses, 3L)
  expect_s3_class(fits$analyses[[1]], "ipw")

  pooled <- mice::pool(fits)
  expect_s3_class(pooled, "mipo")
  expect_identical(pooled$m, 3L)

  # One row per effect measure, in the order the result reports them, and each
  # row combined over all three imputations.
  summarized <- summary(pooled)
  expect_identical(nrow(summarized), 3L)
  expect_identical(
    as.character(summarized$term),
    c("rd", "log(rr)", "log(or)")
  )
  expect_true(all(pooled$pooled$m == 3L))

  expect_true(all(is.finite(summarized$estimate)))
  expect_true(all(is.finite(summarized$std.error)))

  # The complete-data degrees of freedom mice reads off the first fit are the
  # ones M-estimation records, so the Barnard-Rubin degrees of freedom are
  # finite rather than the large-sample fallback.
  expect_identical(
    unique(pooled$pooled$dfcom),
    as.numeric(df.residual(fits$analyses[[1]]))
  )
  expect_true(all(is.finite(summarized$df)))

  # Rubin's first rule is the mean of the per-imputation estimates. Recomputing
  # it from the fits themselves pins that pooling read the tidier's estimates
  # rather than anything else the object carries.
  per_fit <- vapply(
    fits$analyses,
    function(fit) tidy(fit)$estimate[[1]],
    numeric(1)
  )
  expect_equal(pooled$pooled$estimate[[1]], mean(per_fit))
})

# ---- categorical: the grouping the contrast column carries ------------------

test_that("mice::pool() groups categorical ipw results by term and contrast", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_categorical_data(), seed = 2)
  fits <- fit_mice_categorical(imp)

  pooled <- mice::pool(fits)
  summarized <- summary(pooled)

  # The grouping is on the pooled table, which is where the column mice grouped
  # by survives. Three effect measures for each of two contrasts: grouping on
  # `term` alone would return three rows, each combining six estimates from
  # three imputations of two different contrasts.
  expect_identical(nrow(pooled$pooled), 6L)
  expect_true("contrast" %in% names(pooled$pooled))
  expect_identical(
    as.character(pooled$pooled$contrast),
    rep(c("b vs a", "c vs a"), times = 3)
  )
  expect_identical(
    as.character(pooled$pooled$term),
    rep(c("rd", "log(rr)", "log(or)"), each = 2)
  )

  # Each group holds one estimate per imputation, which is the half of the
  # contract a collapsed grouping would break while still returning numbers:
  # six groups of three rather than three groups of six.
  expect_true(all(pooled$pooled$m == 3L))

  # `summary()` reports the same six rows in the same order, keyed by `term`
  # alone, which repeats. What identifies a row there is its position beside the
  # pooled table, so the two are pinned against each other rather than the
  # summary being read on its own.
  expect_identical(nrow(summarized), 6L)
  expect_identical(
    as.character(summarized$term),
    rep(c("rd", "log(rr)", "log(or)"), each = 2)
  )
  expect_equal(summarized$estimate, pooled$pooled$estimate)

  expect_true(all(is.finite(summarized$estimate)))
  expect_true(all(is.finite(summarized$std.error)))
})

# ---- continuous -------------------------------------------------------------

test_that("mice::pool() combines continuous ipw results", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_continuous_data(), seed = 3)
  fits <- fit_mice_continuous(imp)

  pooled <- mice::pool(fits)
  summarized <- summary(pooled)

  # A continuous exposure reports the single coefficient of its marginal
  # structural model, named for the outcome link, so the pooled table holds one
  # row rather than one per effect measure.
  expect_identical(nrow(summarized), 1L)
  expect_identical(as.character(summarized$term), "slope")

  # Read off the pooled table, where a column mice grouped by would appear. The
  # summary drops such columns whether or not there were any, so the absence of
  # one has to be read where its presence would have been visible.
  expect_false("contrast" %in% names(pooled$pooled))

  expect_identical(pooled$pooled$m, 3L)
  expect_true(is.finite(summarized$estimate))
  expect_true(is.finite(summarized$std.error))
  expect_true(is.finite(summarized$df))
})

# ---- linearization: the degrees of freedom it does not record ---------------

test_that("mice::pool() combines linearization estimates but reports no df", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_binary_data(), seed = 1)
  fits <- fit_mice_binary(imp, se_method = "linearization")

  pooled <- mice::pool(fits)
  summarized <- summary(pooled)

  # The estimates pool: the point estimates and the within-imputation variances
  # are there to combine, and Rubin's first two rules need nothing else.
  expect_identical(nrow(summarized), 3L)
  expect_true(all(is.finite(summarized$estimate)))
  expect_true(all(is.finite(summarized$std.error)))

  # The degrees of freedom do not. Linearization stacks no system and records no
  # residual degrees of freedom, so `df.residual()` on the fit is missing and
  # the Barnard-Rubin adjustment mice builds from it is missing with it. This is
  # the current behavior rather than an aspiration: it is pinned so a change to
  # it is deliberate, and the documentation steers users to `dfcom`.
  expect_identical(df.residual(fits$analyses[[1]]), NA_integer_)
  expect_true(all(is.na(pooled$pooled$dfcom)))
  expect_true(all(is.na(summarized$df)))

  # The remedy, which is to name the complete-data degrees of freedom. The
  # outcome model has them, and supplying them restores the adjustment without
  # changing the estimates.
  dfcom <- df.residual(fits$analyses[[1]]$outcome_mod)
  expect_true(is.finite(dfcom))
  repaired <- summary(mice::pool(fits, dfcom = dfcom))
  expect_true(all(is.finite(repaired$df)))
  expect_equal(repaired$estimate, summarized$estimate)
})

# ---- the conditional reading ------------------------------------------------

test_that("mice::pool() combines a conditional reading by coefficient", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_binary_data(), seed = 1)
  fits <- with(imp, {
    ps <- glm(z ~ x1 + x2, family = binomial())
    w <- wt_ate(ps)
    om <- glm(y ~ z, family = quasibinomial(), weights = w)
    as_conditional(ipw(ps, om))
  })

  pooled <- mice::pool(fits)
  summarized <- summary(pooled)

  # A result recording the conditional reading tidies to the outcome model's
  # coefficients, so that is what pools: one row per coefficient rather than one
  # per effect measure.
  coefficients <- names(coef(fits$analyses[[1]]$outcome_mod))
  expect_identical(nrow(summarized), length(coefficients))
  expect_identical(as.character(summarized$term), coefficients)
  expect_true(all(pooled$pooled$m == 3L))
  expect_true(all(is.finite(summarized$estimate)))
  expect_true(all(is.finite(summarized$std.error)))
})

# ---- pool_ipw(), the verb this package recommends ---------------------------
#
# `mice::pool()` above combines the results through the tidier, by the rules it
# applies to any analysis it is given. `pool_ipw()` combines them knowing what
# they are: it reads the results themselves rather than a tidied table of them,
# so the effect labels survive, a categorical result keeps the contrast each row
# reports, and the complete-data degrees of freedom fall back to the outcome
# models when a result records none of its own. This is the path the
# documentation teaches; the one above is the interop.

test_that("pool_ipw() combines ipw results fitted to imputed data", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_binary_data(), seed = 1)
  fits <- fit_mice_binary(imp)

  pooled <- pool_ipw(fits)

  # `mice::with()` returns a `mira`, which the pooling takes directly rather
  # than asking the caller to unwrap it first.
  expect_s3_class(pooled, "ipw_pooled")
  expect_identical(pooled$m, 3L)
  expect_identical(pooled$estimand, "ate")

  tidied <- tidy(pooled)
  expect_identical(nrow(tidied), 3L)
  expect_identical(tidied$term, c("rd", "log(rr)", "log(or)"))
  expect_true(all(is.finite(tidied$estimate)))
  expect_true(all(is.finite(tidied$std.error)))
  expect_true(all(is.finite(tidied$df)))

  expect_identical(glance(pooled)$m, 3L)
})

test_that("pool_ipw() keeps the contrast labels of a categorical result", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_categorical_data(), seed = 2)
  fits <- fit_mice_categorical(imp)

  tidied <- tidy(pool_ipw(fits))

  # The labels reach the pooled table rather than being reconstructed from row
  # order, which is what reading the results themselves buys over reading a
  # tidied table of them.
  expect_identical(nrow(tidied), 6L)
  expect_identical(names(tidied)[seq(1, 2)], c("term", "contrast"))
  expect_identical(
    tidied$term,
    rep(c("rd", "log(rr)", "log(or)"), times = 2)
  )
  expect_identical(
    tidied$contrast,
    rep(c("b vs a", "c vs a"), each = 3)
  )
})

test_that("pool_ipw() reports degrees of freedom a linearization fit lacks", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_binary_data(), seed = 1)
  fits <- fit_mice_binary(imp, se_method = "linearization")

  # The results record no residual degrees of freedom, which is what leaves
  # `mice::pool()` above reporting none unless it is told them. The pooling that
  # reads the results falls back to the outcome models, so the caller is not
  # asked for a number the objects already hold.
  expect_identical(df.residual(fits$analyses[[1]]), NA_integer_)

  tidied <- tidy(pool_ipw(fits))
  expect_true(all(is.finite(tidied$df)))
  expect_true(is.finite(glance(pool_ipw(fits))$dfcom))
})

# ---- both readings of a pooled workflow -------------------------------------
#
# `pool_ipw()` pools both readings of the results it is given rather than the
# one it is asked for, so the reading a pooled result reports can be changed
# afterwards without pooling again. Where the results support only one reading
# the other is recorded as unavailable together with the reason. These run the
# whole imputation path rather than the pooling alone, which is where a reading
# that survives the pooling but not the imputation loop would show up.

test_that("pool_ipw() carries both readings of an mestimation workflow", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_binary_data(), seed = 1)
  fits <- fit_mice_binary(imp)

  pooled <- pool_ipw(fits)
  conditional <- as_conditional(pooled)

  # The conditional reading pools the outcome models' coefficients, so what it
  # reports is named for them rather than for the effect measures the marginal
  # reading contrasts. Read off the models the imputations were fitted with, so
  # the names are the ones the workflow produced rather than ones restated here.
  coefficients <- names(coef(fits$analyses[[1]]$outcome_mod))
  expect_identical(conditional$effects, "conditional")
  expect_identical(conditional$estimates$effect, coefficients)
  expect_identical(estimand(conditional), estimand(pooled))

  # Each reading is pooled once, so a flip changes which one the object reports
  # rather than pooling again: the round trip returns the object it started
  # from rather than one that agrees with it to within a re-estimate.
  expect_identical(as_marginal(conditional), pooled)
})

test_that("pool_ipw() records the reading a linearization workflow lacks", {
  skip_if_not_installed("mice", "3.18.0")
  skip_if_not_installed("deli")
  imp <- mice_impute(mice_binary_data(), seed = 1)
  fits <- fit_mice_binary(imp, se_method = "linearization")

  pooled <- pool_ipw(fits)

  # The conditional reading reports the covariance the joint estimation of the
  # weights and the outcome implies, and linearization stacks no such system, so
  # the outcome models carry none and only the marginal reading pools. The
  # unavailable reading is named and explained rather than left out, which is
  # what lets a flip toward it say why it cannot happen.
  expect_identical(pooled$effects, "marginal")
  expect_identical(names(pooled$alternate), c("effects", "reason"))
  expect_identical(pooled$alternate$effects, "conditional")
  expect_type(pooled$alternate$reason, "character")
  expect_length(pooled$alternate$reason, 1L)
  expect_true(nzchar(pooled$alternate$reason))
})
