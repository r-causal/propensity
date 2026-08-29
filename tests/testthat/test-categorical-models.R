# The weights a fitted `nnet::multinom()` model gives, and the fits the route
# refuses.
#
# A `multinom` reaches a method of its own rather than through inheritance: the
# class is `c("multinom", "nnet")`, so a fit that found no method here would
# reach `wt_ate.default()` and be refused as a class the weights cannot be read
# from.
#
# The contract each weighting test pins is the same one the model methods for a
# continuous exposure hold to in tests/testthat/test-continuous-models.R: a
# model route gives the weights the numeric route gives when it is handed the
# probabilities the model fitted, and nothing else.

# One categorical problem and one binary problem, fit by the same kind of model.
# The three-level exposure is what the categorical path is for; the two-level
# one is there because a `multinom` of two levels fits a single probability and
# so has to be read as a model of a binary exposure.
categorical_model_data <- local({
  set.seed(20250915)

  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)

  odds_b <- exp(0.6 * x1 - 0.3 * x2)
  odds_c <- exp(-0.4 * x1 + 0.5 * x2)
  total <- 1 + odds_b + odds_c
  p_a <- odds_b / total
  p_b <- odds_c / total

  u <- runif(n)

  data.frame(
    x1 = x1,
    x2 = x2,
    trt = factor(ifelse(u < p_a, "a", ifelse(u < p_a + p_b, "b", "c"))),
    a2 = factor(ifelse(
      runif(n) < plogis(0.7 * x1 - 0.4 * x2),
      "control",
      "treated"
    ))
  )
})

fit_multinom <- function(formula) {
  nnet::multinom(formula, data = categorical_model_data, trace = FALSE)
}

# ---- the categorical path ---------------------------------------------------

test_that("a multinomial fit gives the weights its fitted probabilities give", {
  skip_if_not_installed("nnet")

  fit <- fit_multinom(trt ~ x1 + x2)
  trt <- categorical_model_data$trt

  weights <- wt_ate(fit, trt, exposure_type = "categorical")
  oracle <- wt_ate(
    stats::predict(fit, type = "probs"),
    trt,
    exposure_type = "categorical"
  )

  expect_equal(
    as.numeric(weights),
    as.numeric(oracle),
    tolerance = 1e-12
  )
  expect_identical(estimand(weights), "ate")
  expect_identical(exposure_type(weights), "categorical")
})

test_that("the fitted probabilities are matched to the exposure's levels", {
  skip_if_not_installed("nnet")

  fit <- fit_multinom(trt ~ x1 + x2)
  # `fit$lev` is the order the response was coded in, which is the order the
  # fitted matrix carries its columns in. A caller may hand back the same
  # exposure under another level order, and the weights have to be read against
  # the levels the caller supplied rather than the ones the model was fit under.
  reordered <- relevel(categorical_model_data$trt, "c")

  expect_false(identical(levels(reordered), fit$lev))

  weights <- wt_ate(fit, reordered, exposure_type = "categorical")
  oracle <- wt_ate(
    stats::predict(fit, type = "probs")[, levels(reordered)],
    reordered,
    exposure_type = "categorical"
  )

  expect_equal(
    as.numeric(weights),
    as.numeric(oracle),
    tolerance = 1e-12
  )
})

# ---- the binary path --------------------------------------------------------

test_that("a two-level multinomial fit is weighted as a binary exposure", {
  skip_if_not_installed("nnet")

  fit <- fit_multinom(a2 ~ x1 + x2)
  a2 <- categorical_model_data$a2

  # A two-level `multinom` fits one probability rather than a column for each
  # level, and reports it as a single-column matrix of the probability of the
  # second level. That is what the binary path takes.
  weights <- wt_ate(fit, a2)
  oracle <- wt_ate(as.numeric(stats::fitted(fit)), a2)

  expect_equal(
    as.numeric(weights),
    as.numeric(oracle),
    tolerance = 1e-12
  )
  expect_identical(estimand(weights), "ate")
  expect_identical(exposure_type(weights), "binary")
})

test_that("naming the other level of a two-level fit inverts its probabilities", {
  skip_if_not_installed("nnet")

  fit <- fit_multinom(a2 ~ x1 + x2)
  a2 <- categorical_model_data$a2
  first_level <- levels(a2)[[1]]

  weights <- wt_ate(fit, a2, .focal_level = first_level)
  oracle <- wt_ate(
    1 - as.numeric(stats::fitted(fit)),
    a2,
    .focal_level = first_level
  )

  expect_equal(
    as.numeric(weights),
    as.numeric(oracle),
    tolerance = 1e-12
  )

  # The inversion is what separates these weights from the ones the same call
  # would give if the fitted probabilities were read as the probability of the
  # focal level rather than of the other one. The default-focal weights are no
  # such comparison: an ATE weight is the reciprocal of the probability of the
  # level each unit was observed at, which both focal levels give.
  expect_false(isTRUE(all.equal(
    as.numeric(weights),
    as.numeric(wt_ate(
      as.numeric(stats::fitted(fit)),
      a2,
      .focal_level = first_level
    ))
  )))
})

# ---- the exposure the model already holds -----------------------------------

test_that("a multinomial fit reads its own exposure and says which it read", {
  skip_if_not_installed("nnet")

  fit <- fit_multinom(trt ~ x1 + x2)
  trt <- categorical_model_data$trt

  # Read while the alerts are still silenced, so that only the call under test
  # has anything to say.
  supplied <- wt_ate(fit, trt, exposure_type = "categorical")

  withr::local_options(propensity.quiet = FALSE)

  # Both alerts are the route's own account of what it decided: which variable
  # it read as the exposure, and which type it resolved that variable to.
  expect_snapshot(from_model <- wt_ate(fit))

  expect_identical(as.numeric(from_model), as.numeric(supplied))
})

# ---- the fits the route refuses ---------------------------------------------

test_that("a multinomial fit refuses a continuous exposure", {
  skip_if_not_installed("nnet")

  fit <- fit_multinom(trt ~ x1 + x2)

  # A multinomial model fits a probability for each level, which is not a
  # conditional mean with a spread, so a continuous exposure is not one of the
  # types this route offers and the type is refused before anything is read off
  # the model.
  expect_error(
    wt_ate(
      fit,
      categorical_model_data$x1,
      exposure_type = "continuous"
    ),
    class = "causalgenerics_unsupported_exposure_type"
  )
})

test_that("a three-level fit refuses a two-level exposure", {
  skip_if_not_installed("nnet")

  fit <- fit_multinom(trt ~ x1 + x2)
  collapsed <- factor(ifelse(
    categorical_model_data$trt == "a",
    "a",
    "not a"
  ))

  # Two levels resolve as a binary exposure, which needs the probability of one
  # of them. A three-level fit has no such column to give.
  expect_error(
    wt_ate(fit, collapsed),
    class = "propensity_model_family_error"
  )
})

test_that("a multinomial fit of a response matrix is refused", {
  skip_if_not_installed("nnet")

  counts <- with(
    categorical_model_data,
    cbind(
      a = as.integer(trt == "a"),
      b = as.integer(trt == "b"),
      c = as.integer(trt == "c")
    )
  )
  fit <- nnet::multinom(
    counts ~ x1 + x2,
    data = categorical_model_data,
    trace = FALSE
  )

  # A matrix response is read as counts rather than as a factor, so the fit
  # records no levels and there is nothing to match the exposure's levels
  # against.
  expect_length(fit$lev, 0L)

  expect_error(
    wt_ate(
      fit,
      categorical_model_data$trt,
      exposure_type = "categorical"
    ),
    class = "propensity_model_family_error"
  )
})

test_that("censoring weights take no multinomial fit", {
  skip_if_not_installed("nnet")

  # Censoring is the one weight that is not a weight for an exposure, so the
  # multinomial route is not extended to it and a fit reaches the default
  # method. This stays a refusal after the other weight functions gain their
  # multinomial methods.
  fit <- fit_multinom(trt ~ x1 + x2)

  expect_error(
    wt_cens(fit, categorical_model_data$trt),
    class = "propensity_method_error"
  )
})

# ---- agreement with WeightIt ------------------------------------------------

test_that("multinomial weights match WeightIt for every estimand", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("WeightIt")
  skip_on_cran()

  # WeightIt fits the multinomial model with its own maximum likelihood routine
  # rather than with `nnet::multinom()`, so the two sets of fitted probabilities
  # agree to the accuracy of the optimizers rather than exactly.
  fit <- fit_multinom(trt ~ x1 + x2)
  trt <- categorical_model_data$trt

  weightit_weights <- function(estimand, focal = NULL) {
    args <- list(
      trt ~ x1 + x2,
      data = categorical_model_data,
      method = "glm",
      estimand = estimand
    )
    if (!is.null(focal)) {
      args$focal <- focal
    }

    do.call(WeightIt::weightit, args)$weights
  }

  expect_equal(
    as.numeric(wt_ate(fit, trt, exposure_type = "categorical")),
    weightit_weights("ATE"),
    tolerance = 1e-5,
    ignore_attr = "names"
  )

  # The tilting estimands have no multinomial method yet, so everything below
  # fails until they gain one.
  for (focal in levels(trt)) {
    expect_equal(
      as.numeric(wt_att(
        fit,
        trt,
        .focal_level = focal,
        exposure_type = "categorical"
      )),
      weightit_weights("ATT", focal = focal),
      tolerance = 1e-5,
      ignore_attr = "names",
      label = paste("ATT weights with focal =", focal)
    )
  }

  expect_equal(
    as.numeric(wt_ato(fit, trt, exposure_type = "categorical")),
    weightit_weights("ATO"),
    tolerance = 1e-5,
    ignore_attr = "names"
  )

  expect_equal(
    as.numeric(wt_atm(fit, trt, exposure_type = "categorical")),
    weightit_weights("ATM"),
    tolerance = 1e-5,
    ignore_attr = "names"
  )

  # WeightIt has no entropy estimand for a multinomial treatment, so the
  # entropy weights are held to the numeric route instead.
  expect_equal(
    as.numeric(wt_entropy(fit, trt, exposure_type = "categorical")),
    as.numeric(wt_entropy(
      stats::predict(fit, type = "probs"),
      trt,
      exposure_type = "categorical"
    )),
    tolerance = 1e-12
  )
})
