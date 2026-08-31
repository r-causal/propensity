# A fitted numerator model for a binary exposure. The continuous route already
# takes a model in `stabilize` and stabilizes the density ratio on the
# conditional density that model estimates; the same argument names the
# conditional probability of the exposure for a binary one. The numerator is
# P(Z = z_i | V_i) read at the level each unit took, so the weights are the
# unstabilized ones scaled by the fitted probability of the observed level.
#
# The reason to want one is an effect-modification marginal structural model:
# the numerator conditions on the modifier the reported model already reads,
# which leaves the estimator consistent for the same effects and cancels the
# part of the exposure's variation the modifier explains. `ipw()`'s `.by`
# documentation names that as the recommended numerator, and a fitted model is
# the form of it whose estimation the sandwich can account for.

# ---- fixtures ---------------------------------------------------------------

# A binary exposure, a modifier the numerator conditions on, and a covariate
# only the propensity score model reads, so the two models are genuinely
# different and the stabilized weights are not the unstabilized ones rescaled by
# a constant.
binary_numerator_data <- function(seed = 8102, n = 300) {
  withr::with_seed(seed, {
    x <- stats::rnorm(n)
    v <- stats::rbinom(n, 1, 0.5)
    z <- stats::rbinom(n, 1, stats::plogis(0.4 * x - 0.7 * v))
    y <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * z + 0.5 * x + 0.6 * v))
    data.frame(x, v, z, y)
  })
}

binary_numerator_dat <- binary_numerator_data()

# The propensity score model, its fitted scores, and the numerator model that
# conditions on the modifier alone.
binary_numerator_ps_mod <- stats::glm(
  z ~ x + v,
  data = binary_numerator_dat,
  family = stats::binomial()
)

binary_numerator_ps <- as.numeric(stats::fitted(binary_numerator_ps_mod))

binary_numerator_model <- stats::glm(
  z ~ v,
  data = binary_numerator_dat,
  family = stats::binomial()
)

# The numerator the model evaluates to: the fitted probability of the level each
# unit took, which is what multiplies the unstabilized weight.
binary_numerator_score <- function(model = binary_numerator_model) {
  z <- binary_numerator_dat$z
  p <- as.numeric(stats::fitted(model))
  z * p + (1 - z) * (1 - p)
}

binary_unstabilized_wt <- function() {
  z <- binary_numerator_dat$z
  e <- binary_numerator_ps
  (z / e) + ((1 - z) / (1 - e))
}

# ---- the weights a binary numerator model builds ----------------------------

test_that("a binary numerator model stabilizes on the probability of the level taken", {
  weights <- wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = binary_numerator_model
  )

  expect_s3_class(weights, "psw")
  expect_identical(estimand(weights), "ate")
  expect_true(is_stabilized(weights))
  expect_equal(
    as.numeric(weights),
    binary_unstabilized_wt() * binary_numerator_score(),
    tolerance = 1e-12
  )
})

test_that("a binary numerator model is not the marginal stabilizer", {
  # The default stabilizer is the marginal probability of the exposure, which is
  # one number rather than one per unit, so a conditional numerator builds a
  # different set of weights rather than the same ones under another name.
  modeled <- wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = binary_numerator_model
  )
  marginal <- wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = TRUE
  )

  expect_false(isTRUE(all.equal(
    as.numeric(modeled),
    as.numeric(marginal),
    tolerance = 1e-8
  )))

  # And it is the unstabilized weights scaled unit by unit, which is what makes
  # the comparison against the marginal one a comparison of numerators.
  expect_equal(
    as.numeric(modeled) / binary_unstabilized_wt(),
    binary_numerator_score(),
    tolerance = 1e-12
  )
})

test_that("the binary weights record the numerator model they were stabilized on", {
  weights <- wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = binary_numerator_model
  )

  # The model itself is recorded rather than the numerator it evaluates to, for
  # the reason the continuous route records it: `ipw()` rebuilds that numerator
  # at every value of theta, which takes the model's design and its
  # coefficients.
  expect_identical(numerator_model(weights), binary_numerator_model)
  expect_identical(exposure_type(weights), "binary")

  # A model is a numerator rather than a score, so nothing is recorded as one.
  expect_null(stabilization_score(weights))

  # A binary exposure's weights are a ratio of probabilities rather than of
  # densities, so there is still no density record to leave.
  expect_null(density_meta(weights))
})

test_that("weights stabilized on no model record no numerator model", {
  marginal <- wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = TRUE
  )
  unstabilized <- wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary"
  )
  scored <- wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = TRUE,
    stabilization_score = binary_numerator_score()
  )

  expect_null(numerator_model(marginal))
  expect_null(numerator_model(unstabilized))
  expect_null(numerator_model(scored))
})

test_that("the binary weights name the numerator model where the continuous ones do", {
  weights <- wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = binary_numerator_model
  )

  # A continuous set of weights prints the model's formula under its values,
  # labelled by the argument the model arrived in. A binary set records the same
  # model and names it the same way, since the reader's question is the same
  # one: which numerator are these weights divided by.
  printed <- capture.output(print(weights))
  expect_true(any(grepl("stabilize:", printed, fixed = TRUE)))
  expect_true(any(grepl("z ~ v", printed, fixed = TRUE)))
})

test_that("censoring weights reach the binary numerator model route too", {
  # `wt_cens()` builds its weights with the ATE formula, so a numerator model
  # arrives there through the same argument and is recorded the same way.
  weights <- wt_cens(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = binary_numerator_model
  )

  expect_identical(estimand(weights), "uncensored")
  expect_true(is_stabilized(weights))
  expect_equal(
    as.numeric(weights),
    binary_unstabilized_wt() * binary_numerator_score(),
    tolerance = 1e-12
  )
  expect_identical(numerator_model(weights), binary_numerator_model)
})

test_that("a binary numerator model reaches the weights through a fitted model", {
  # The model route forwards `stabilize` to the numeric one untouched, so a
  # propensity score model and a numerator model together build the same weights
  # the scores and the numerator model build.
  weights <- wt_ate(
    binary_numerator_ps_mod,
    stabilize = binary_numerator_model
  )

  expect_equal(
    as.numeric(weights),
    binary_unstabilized_wt() * binary_numerator_score(),
    tolerance = 1e-12
  )
  expect_identical(numerator_model(weights), binary_numerator_model)
})

# ---- what a binary numerator model has to be --------------------------------

test_that("a binary numerator model must model the exposure with a binomial family", {
  # The numerator is a probability of the level each unit took, which is what a
  # binomial fit reports and what a fit of any other family does not. The
  # numerator model is held to the family the propensity score model of a binary
  # exposure is held to, read from the other side of the ratio.
  gaussian_fit <- stats::lm(z ~ v, data = binary_numerator_dat)

  expect_error(
    wt_ate(
      binary_numerator_ps,
      binary_numerator_dat$z,
      exposure_type = "binary",
      stabilize = gaussian_fit
    ),
    class = "propensity_model_family_error"
  )
})

test_that("a binary numerator model fit to other observations is refused", {
  # A model fit to a different number of observations has one numerator for each
  # of its own units and none for the units here, so it is refused rather than
  # recycled against them.
  short <- stats::glm(
    z ~ v,
    data = utils::head(binary_numerator_dat, 100),
    family = stats::binomial()
  )

  expect_error(
    wt_ate(
      binary_numerator_ps,
      binary_numerator_dat$z,
      exposure_type = "binary",
      stabilize = short
    ),
    class = "propensity_length_error"
  )
})

test_that("a binary numerator model and a stabilization score are two numerators", {
  # A score the caller computed is itself the numerator of the weights, and a
  # model is a second one; the two together are two instructions about the same
  # quantity, which is the refusal the continuous route makes of the same pair.
  expect_error(
    wt_ate(
      binary_numerator_ps,
      binary_numerator_dat$z,
      exposure_type = "binary",
      stabilize = binary_numerator_model,
      stabilization_score = binary_numerator_score()
    ),
    class = "propensity_numerator_error"
  )
})

test_that("a numerator model of the wrong exposure type is refused either way", {
  # Which numerator a model estimates is a statement about the exposure it
  # models, so a model of a binary exposure names no density for a dose to
  # divide by and a model of a dose names no probability for a binary exposure
  # to divide by. Both are refused, and both for the same reason: the family the
  # model was fit with is not the family the ratio needs.
  dose <- withr::with_seed(
    8103,
    binary_numerator_dat$x + stats::rnorm(nrow(binary_numerator_dat))
  )
  dose_mu <- as.numeric(stats::fitted(
    stats::lm(dose ~ binary_numerator_dat$x)
  ))

  expect_error(
    wt_ate(
      dose_mu,
      dose,
      exposure_type = "continuous",
      stabilize = binary_numerator_model
    ),
    class = "propensity_model_family_error"
  )

  dose_numerator <- stats::lm(dose ~ binary_numerator_dat$v)
  expect_error(
    wt_ate(
      binary_numerator_ps,
      binary_numerator_dat$z,
      exposure_type = "binary",
      stabilize = dose_numerator
    ),
    class = "propensity_model_family_error"
  )
})

test_that("a categorical exposure still refuses a numerator model", {
  # A categorical exposure's weights are a ratio of probabilities over more than
  # two levels, and a fitted model of one is not what `stabilize` reads here.
  exposure <- factor(rep(c("a", "b", "c"), each = 4))
  categorical_ps <- matrix(
    rep(c(0.5, 0.3, 0.2), times = 12),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )
  fit <- stats::glm(
    rep(c(0, 1), length.out = 12) ~ seq_len(12),
    family = stats::binomial()
  )

  expect_error(
    wt_ate(
      categorical_ps,
      exposure,
      exposure_type = "categorical",
      stabilize = fit
    ),
    class = "propensity_numerator_error"
  )
})

# ---- combining weights that record a numerator model ------------------------

# The weights the combine tests are written about: stabilized on the numerator
# model, on a second model of other terms, and on the marginal probability of
# the level each unit took.
binary_modeled_wt <- function(model = binary_numerator_model) {
  wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = model
  )
}

binary_marginal_wt <- function() {
  wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary",
    stabilize = TRUE
  )
}

binary_other_numerator_model <- stats::glm(
  z ~ x,
  data = binary_numerator_dat,
  family = stats::binomial()
)

test_that("combining a modeled numerator with the marginal one drops the model", {
  # Half of the combined weights were divided by the model's fitted probability
  # and half by the marginal one, so neither numerator describes the result.
  # This is the disagreement the continuous route reports of the same pair,
  # where the numerator lives inside the density record: the record is dropped
  # and the drop is named, rather than the model being carried onto weights
  # that were not all stabilized on it.
  modeled <- binary_modeled_wt()
  marginal <- binary_marginal_wt()

  out <- NULL
  expect_warning(
    out <- vec_c(modeled, marginal),
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(numerator_model(out))

  # The result is still stabilized, and on nothing it can name. That is what
  # the continuous conflict leaves too: the density record goes and the
  # stabilization status stays.
  expect_true(is_stabilized(out))
  expect_false(any(grepl(
    "stabilize:",
    capture.output(print(out)),
    fixed = TRUE
  )))

  # Which operand was written first says nothing about which numerator the
  # result was divided by.
  reversed <- NULL
  expect_warning(
    reversed <- vec_c(marginal, modeled),
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(numerator_model(reversed))
  expect_true(is_stabilized(reversed))
})

test_that("combining weights stabilized on the same model keeps it", {
  # A numerator model is compared by what it says rather than by identity, the
  # way the model inside a density record is: a fit carries the frame it was
  # made in and the call that made it, so the same numerator fit in two calls
  # would otherwise read as a disagreement.
  again <- stats::glm(
    z ~ v,
    data = binary_numerator_dat,
    family = stats::binomial()
  )

  same <- expect_silent(vec_c(binary_modeled_wt(), binary_modeled_wt(again)))
  expect_identical(
    stats::coef(numerator_model(same)),
    stats::coef(binary_numerator_model)
  )
})

test_that("combining weights stabilized on different models drops both", {
  # Two models that read different terms estimate two different numerators, so
  # the combination was divided by neither.
  out <- NULL
  expect_warning(
    out <- vec_c(
      binary_modeled_wt(),
      binary_modeled_wt(binary_other_numerator_model)
    ),
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(numerator_model(out))
})

# ---- arithmetic on weights that record a numerator model --------------------

test_that("arithmetic that leaves the result unstabilized drops the model", {
  # A numerator model names what the weights were divided by, so a result only
  # half of which was divided by it was not. The merge already drops the
  # stabilization status and the stabilization score there, and the model is
  # the same record of the same thing: left standing, the result answers
  # `numerator_model()` with a model and prints a `stabilize:` line naming a
  # numerator that never multiplied it.
  modeled <- binary_modeled_wt()
  unstabilized <- wt_ate(
    binary_numerator_ps,
    binary_numerator_dat$z,
    exposure_type = "binary"
  )

  results <- list(
    `modeled * unstabilized` = expect_silent(modeled * unstabilized),
    `unstabilized * modeled` = expect_silent(unstabilized * modeled),
    `modeled + unstabilized` = expect_silent(modeled + unstabilized)
  )

  for (label in names(results)) {
    out <- results[[label]]
    expect_false(is_stabilized(out), label = label)
    expect_null(numerator_model(out), label = label)
    expect_false(
      any(grepl("stabilize:", capture.output(print(out)), fixed = TRUE)),
      label = label
    )
  }
})

test_that("arithmetic that stays stabilized keeps the numerator model", {
  # The drop above belongs to the unstabilized result rather than to every
  # product: two sets of weights divided by the same numerator give a product
  # that was divided by it too.
  same <- expect_silent(binary_modeled_wt() * binary_modeled_wt())
  expect_true(is_stabilized(same))
  expect_identical(
    stats::coef(numerator_model(same)),
    stats::coef(binary_numerator_model)
  )
})
