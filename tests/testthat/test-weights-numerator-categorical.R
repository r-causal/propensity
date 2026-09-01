# A fitted numerator model for a categorical exposure. The binary route already
# takes a model in `stabilize` and stabilizes the probability ratio on the
# conditional probability that model estimates; the same argument names the
# conditional probability of the level each unit took when the exposure has more
# than two of them. The numerator is P(Z = z_i | V_i) read off the column of the
# model's fitted matrix that belongs to the observed level, so the weights are
# the unstabilized ones scaled by the fitted probability of that level.
#
# The reason to want one is the reason the binary route has one: an
# effect-modification marginal structural model whose numerator conditions on
# the modifier the reported model already reads, which leaves the estimator
# consistent for the same effects and cancels the part of the exposure's
# variation the modifier explains. `nnet::multinom()` is the fit that reports a
# probability for every level, so it is the model the categorical route reads,
# the way it is the model the categorical propensity score is read from.

# ---- fixtures ---------------------------------------------------------------

# A three-level exposure, a modifier the numerator conditions on, and a
# covariate only the propensity score model reads, so the two models are
# genuinely different and the stabilized weights are not the unstabilized ones
# rescaled by one number per level.
categorical_numerator_data <- function(seed = 8104, n = 300) {
  withr::with_seed(seed, {
    x <- stats::rnorm(n)
    v <- stats::rbinom(n, 1, 0.5)
    eta_b <- -0.2 + 0.6 * x - 0.8 * v
    eta_c <- 0.1 - 0.5 * x + 0.7 * v
    denominator <- 1 + exp(eta_b) + exp(eta_c)
    p_a <- 1 / denominator
    p_b <- exp(eta_b) / denominator
    u <- stats::runif(n)
    level <- ifelse(u < p_a, "a", ifelse(u < p_a + p_b, "b", "c"))
    data.frame(x, v, z = factor(level, levels = c("a", "b", "c")))
  })
}

categorical_numerator_dat <- categorical_numerator_data()

# Every multinomial fixture is fit to a tighter convergence than
# `nnet::multinom()`'s default, the convention test-ipw-categorical.R states:
# at the default tolerance the reported coefficients sit a little off the
# multinomial score root, which is far too coarse for the exact equalities
# these tests are written at.
fit_categorical_numerator_model <- function(formula, data, ...) {
  nnet::multinom(
    formula,
    data = data,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000,
    ...
  )
}

# The fits the file is written about, built on first use rather than when the
# file is sourced, so that a session without nnet skips these tests rather than
# failing to source them. Each is built once, so the recorded model can be
# compared to the fixture by identity.
categorical_numerator_fit <- local({
  fits <- list()

  function(name) {
    if (is.null(fits[[name]])) {
      fits[[name]] <<- switch(
        name,
        ps_mod = fit_categorical_numerator_model(
          z ~ x + v,
          categorical_numerator_dat
        ),
        num_mod = fit_categorical_numerator_model(
          z ~ v,
          categorical_numerator_dat
        ),
        other_num_mod = fit_categorical_numerator_model(
          z ~ x,
          categorical_numerator_dat
        )
      )
    }

    fits[[name]]
  }
})

categorical_numerator_ps_mod <- function() {
  categorical_numerator_fit("ps_mod")
}

categorical_numerator_model <- function() {
  categorical_numerator_fit("num_mod")
}

# The propensity score matrix the weights are read off: the fit's own fitted
# values, whose columns are named for the levels they belong to.
categorical_numerator_ps <- function() {
  stats::fitted(categorical_numerator_ps_mod())
}

# The numerator the model evaluates to, gathered by hand: the fitted probability
# of the level each unit took, read out of the column named for that level
# rather than out of the column in that position.
categorical_numerator_gather <- function(
  model = categorical_numerator_model(),
  exposure = categorical_numerator_dat$z
) {
  fitted <- stats::fitted(model)
  column <- match(as.character(exposure), colnames(fitted))

  fitted[cbind(seq_along(exposure), column)]
}

# The unstabilized ATE weights, gathered the same way: the reciprocal of the
# propensity score of the level each unit took.
categorical_unstabilized_wt <- function(
  exposure = categorical_numerator_dat$z
) {
  ps <- categorical_numerator_ps()
  column <- match(as.character(exposure), colnames(ps))

  1 / ps[cbind(seq_along(exposure), column)]
}

# The weights the record, combine, and arithmetic tests are written about.
categorical_modeled_wt <- function(model = categorical_numerator_model()) {
  wt_ate(
    categorical_numerator_ps(),
    categorical_numerator_dat$z,
    exposure_type = "categorical",
    stabilize = model
  )
}

categorical_marginal_wt <- function() {
  wt_ate(
    categorical_numerator_ps(),
    categorical_numerator_dat$z,
    exposure_type = "categorical",
    stabilize = TRUE
  )
}

categorical_unstabilized_psw <- function() {
  wt_ate(
    categorical_numerator_ps(),
    categorical_numerator_dat$z,
    exposure_type = "categorical"
  )
}

# A condition's message with its line breaks and indentation flattened, so that
# a bullet written across several lines is matched as the sentence it reads as.
categorical_numerator_message <- function(cnd) {
  gsub("[[:space:]]+", " ", conditionMessage(cnd))
}

# ---- the weights a categorical numerator model builds ------------------------

test_that("a categorical numerator model stabilizes on the probability of the level taken", {
  skip_if_not_installed("nnet")

  weights <- wt_ate(
    categorical_numerator_ps(),
    categorical_numerator_dat$z,
    exposure_type = "categorical",
    stabilize = categorical_numerator_model()
  )

  expect_s3_class(weights, "psw")
  expect_identical(estimand(weights), "ate")
  expect_true(is_stabilized(weights))
  expect_equal(
    as.numeric(weights),
    categorical_unstabilized_wt() * categorical_numerator_gather(),
    tolerance = 1e-12
  )
})

test_that("a categorical numerator model is not the marginal stabilizer", {
  skip_if_not_installed("nnet")

  # The default stabilizer is the marginal probability of each level, which is
  # one number per level rather than one per unit, so a conditional numerator
  # builds a different set of weights rather than the same ones under another
  # name.
  modeled <- categorical_modeled_wt()
  marginal <- categorical_marginal_wt()

  expect_false(isTRUE(all.equal(
    as.numeric(modeled),
    as.numeric(marginal),
    tolerance = 1e-8
  )))

  # And it is the unstabilized weights scaled unit by unit, which is what makes
  # the comparison against the marginal one a comparison of numerators.
  expect_equal(
    as.numeric(modeled) / categorical_unstabilized_wt(),
    categorical_numerator_gather(),
    tolerance = 1e-12
  )
})

test_that("a unit with no observed level takes no numerator from the model", {
  skip_if_not_installed("nnet")

  # A unit whose level is missing has no column of the fitted matrix to read, so
  # the gather leaves its numerator missing, which is the answer the denominator
  # gather and the marginal stabilizer already give such a unit.
  missing_at <- c(3L, 17L)
  exposure <- categorical_numerator_dat$z
  exposure[missing_at] <- NA

  weights <- wt_ate(
    categorical_numerator_ps(),
    exposure,
    exposure_type = "categorical",
    stabilize = categorical_numerator_model()
  )

  expect_true(all(is.na(as.numeric(weights)[missing_at])))
  expect_equal(
    as.numeric(weights)[-missing_at],
    (categorical_unstabilized_wt() * categorical_numerator_gather())[
      -missing_at
    ],
    tolerance = 1e-12
  )
})

# ---- what the weights record about the numerator model -----------------------

test_that("the categorical weights record the numerator model they were stabilized on", {
  skip_if_not_installed("nnet")

  weights <- categorical_modeled_wt()

  # The model itself is recorded rather than the numerator it evaluates to, for
  # the reason the binary and continuous routes record it: `ipw()` rebuilds that
  # numerator at every value of theta, which takes the model's design and its
  # coefficients.
  expect_identical(numerator_model(weights), categorical_numerator_model())
  expect_identical(exposure_type(weights), "categorical")

  # A model is a numerator rather than a score, so nothing is recorded as one.
  expect_null(stabilization_score(weights))

  # A categorical exposure's weights are a ratio of probabilities rather than of
  # densities, so there is no density record to leave.
  expect_null(density_meta(weights))

  # The record rides alongside what the weights already recorded about the
  # levels rather than in place of it.
  expect_identical(attr(weights, "n_categories"), 3L)
  expect_identical(attr(weights, "category_names"), c("a", "b", "c"))
})

test_that("categorical weights stabilized on no model record no numerator model", {
  skip_if_not_installed("nnet")

  scored <- wt_ate(
    categorical_numerator_ps(),
    categorical_numerator_dat$z,
    exposure_type = "categorical",
    stabilize = TRUE,
    stabilization_score = categorical_numerator_gather()
  )

  expect_null(numerator_model(categorical_marginal_wt()))
  expect_null(numerator_model(categorical_unstabilized_psw()))
  expect_null(numerator_model(scored))
})

test_that("the categorical weights name the numerator model where the binary ones do", {
  skip_if_not_installed("nnet")

  # A binary set of weights prints the model's formula under its values,
  # labelled by the argument the model arrived in. A categorical set records the
  # same model and names it the same way, since the reader's question is the
  # same one: which numerator are these weights divided by.
  printed <- capture.output(print(categorical_modeled_wt()))

  expect_true(any(grepl("stabilize:", printed, fixed = TRUE)))
  expect_true(any(grepl("stabilize: z ~ v", printed, fixed = TRUE)))
})

test_that("a categorical numerator model reaches the weights through a fitted model", {
  skip_if_not_installed("nnet")

  # The multinomial route forwards `stabilize` to the numeric one untouched, so
  # a propensity score model and a numerator model together build the same
  # weights the fitted scores and the numerator model build.
  weights <- wt_ate(
    categorical_numerator_ps_mod(),
    categorical_numerator_dat$z,
    exposure_type = "categorical",
    stabilize = categorical_numerator_model()
  )

  expect_equal(
    as.numeric(weights),
    categorical_unstabilized_wt() * categorical_numerator_gather(),
    tolerance = 1e-12
  )
  expect_identical(numerator_model(weights), categorical_numerator_model())
})

# ---- what a categorical numerator model has to be ----------------------------

test_that("a categorical exposure takes a numerator model rather than refusing one", {
  skip_if_not_installed("nnet")

  # This route used to be closed: `stabilize` read a fitted model only for
  # binary and continuous exposures, and a categorical one refused every model
  # whatever it was, on the grounds that no single fit reports a probability for
  # every level. `nnet::multinom()` reports exactly that, so it is taken; what
  # is refused now is a fit that reports one value per unit instead.
  exposure <- categorical_numerator_dat$z
  taken <- wt_ate(
    categorical_numerator_ps(),
    exposure,
    exposure_type = "categorical",
    stabilize = categorical_numerator_model()
  )

  expect_true(is_stabilized(taken))
  expect_identical(numerator_model(taken), categorical_numerator_model())

  one_value_per_unit <- stats::glm(
    (exposure == "c") ~ v,
    data = categorical_numerator_dat,
    family = stats::binomial()
  )

  expect_error(
    wt_ate(
      categorical_numerator_ps(),
      exposure,
      exposure_type = "categorical",
      stabilize = one_value_per_unit
    ),
    class = "propensity_model_family_error"
  )
})

test_that("a categorical numerator model must fit a probability for every level", {
  skip_if_not_installed("nnet")

  # The numerator is a probability for the level each unit took, which is what a
  # multinomial fit reports and what a fit of one value per unit does not. The
  # numerator model is held to what the propensity score model of a categorical
  # exposure is held to, read from the other side of the ratio, so the refusal
  # names the fit that reports it.
  gaussian_fit <- stats::lm(x ~ v, data = categorical_numerator_dat)

  err <- expect_error(
    wt_ate(
      categorical_numerator_ps(),
      categorical_numerator_dat$z,
      exposure_type = "categorical",
      stabilize = gaussian_fit
    ),
    class = "propensity_model_family_error"
  )

  message <- categorical_numerator_message(err)
  expect_match(message, "stabilize", fixed = TRUE)
  expect_match(message, "nnet::multinom", fixed = TRUE)

  binomial_fit <- stats::glm(
    v ~ x,
    data = categorical_numerator_dat,
    family = stats::binomial()
  )

  expect_error(
    wt_ate(
      categorical_numerator_ps(),
      categorical_numerator_dat$z,
      exposure_type = "categorical",
      stabilize = binomial_fit
    ),
    class = "propensity_model_family_error"
  )
})

test_that("a multinomial numerator model is refused for the other exposure types", {
  skip_if_not_installed("nnet")

  # Which numerator a model estimates is a statement about the exposure it
  # models, so a fit reporting a probability for each of three levels names no
  # single probability for a binary exposure to divide by and no density for a
  # dose to divide by. Each is refused with the fit its own side of the ratio
  # takes.
  binary_exposure <- as.integer(categorical_numerator_dat$z == "c")
  binary_ps <- as.numeric(stats::fitted(stats::glm(
    binary_exposure ~ x + v,
    data = categorical_numerator_dat,
    family = stats::binomial()
  )))

  binary_err <- expect_error(
    wt_ate(
      binary_ps,
      binary_exposure,
      exposure_type = "binary",
      stabilize = categorical_numerator_model()
    ),
    class = "propensity_model_family_error"
  )
  expect_match(
    categorical_numerator_message(binary_err),
    "binomial",
    fixed = TRUE
  )

  dose <- withr::with_seed(
    8105,
    categorical_numerator_dat$x + stats::rnorm(nrow(categorical_numerator_dat))
  )
  dose_mu <- as.numeric(stats::fitted(
    stats::lm(dose ~ categorical_numerator_dat$x)
  ))

  dose_err <- expect_error(
    wt_ate(
      dose_mu,
      dose,
      exposure_type = "continuous",
      stabilize = categorical_numerator_model()
    ),
    class = "propensity_model_family_error"
  )
  expect_match(
    categorical_numerator_message(dose_err),
    "stats::lm",
    fixed = TRUE
  )
})

test_that("a categorical numerator model fit to other levels is refused", {
  skip_if_not_installed("nnet")

  # The columns the numerator is gathered from are named for the levels they
  # belong to, so a fit made to some other set of levels has no column to read
  # as the probability of the level a unit took. Order is no part of the
  # comparison here, the gather being by name.
  relabeled <- factor(
    as.character(categorical_numerator_dat$z),
    levels = c("a", "b", "c"),
    labels = c("d", "e", "f")
  )
  elsewhere <- fit_categorical_numerator_model(
    relabeled ~ v,
    transform(categorical_numerator_dat, relabeled = relabeled)
  )

  err <- expect_error(
    wt_ate(
      categorical_numerator_ps(),
      categorical_numerator_dat$z,
      exposure_type = "categorical",
      stabilize = elsewhere
    ),
    class = "propensity_model_family_error"
  )

  message <- categorical_numerator_message(err)
  expect_match(message, "stabilize", fixed = TRUE)

  # Dropping the exposure's unused levels is no remedy for a fit reporting
  # levels the exposure does not have, so that bullet belongs to the other
  # direction of the mismatch alone.
  expect_no_match(message, "droplevels", fixed = TRUE)
})

test_that("an exposure carrying an unused level sends the numerator model to droplevels()", {
  skip_if_not_installed("nnet")

  # `nnet::multinom()` drops a level no unit took, so a numerator fit to an
  # exposure with an empty level reports one column fewer than the exposure
  # declares levels. The two sets disagree for a reason the caller can fix
  # without refitting anything, and the refusal names it, the way the propensity
  # score model's refusal does.
  with_empty_level <- factor(
    as.character(categorical_numerator_dat$z),
    levels = c("a", "b", "c", "d")
  )

  # The empty level is what the fixture is about, so the fit's report of it is
  # expected rather than surfaced in the test run.
  numerator <- suppressWarnings(fit_categorical_numerator_model(
    with_empty_level ~ v,
    transform(categorical_numerator_dat, with_empty_level = with_empty_level)
  ))

  # A propensity score matrix wide enough for the declared levels, so the
  # disagreement under test is the numerator's rather than the denominator's.
  ps <- cbind(
    categorical_numerator_ps() * (1 - 1e-6),
    d = 1e-6
  )

  expect_length(numerator$lev, 3L)

  err <- expect_error(
    wt_ate(
      ps,
      with_empty_level,
      exposure_type = "categorical",
      stabilize = numerator
    ),
    class = "propensity_model_family_error"
  )

  expect_match(
    categorical_numerator_message(err),
    "droplevels",
    fixed = TRUE
  )
})

test_that("a numerator model fit to a response matrix is refused", {
  skip_if_not_installed("nnet")

  # A matrix response is read as counts rather than as a factor, so the fit
  # records no levels and its columns name no level the gather could match. It
  # is refused where the propensity score model fit the same way is.
  counts <- with(
    categorical_numerator_dat,
    cbind(
      a = as.integer(z == "a"),
      b = as.integer(z == "b"),
      c = as.integer(z == "c")
    )
  )
  matrix_fit <- fit_categorical_numerator_model(
    counts ~ v,
    transform(categorical_numerator_dat, counts = counts)
  )

  expect_length(matrix_fit$lev, 0L)

  expect_error(
    wt_ate(
      categorical_numerator_ps(),
      categorical_numerator_dat$z,
      exposure_type = "categorical",
      stabilize = matrix_fit
    ),
    class = "propensity_model_family_error"
  )
})

test_that("a categorical numerator model fit with case weights is refused", {
  skip_if_not_installed("nnet")

  # The refusal the binary and continuous numerator models take, taken here for
  # the same reasons: the numerator is a probability in the sample the weights
  # are being built for, and a fit made under case weights reports it for a
  # reweighted one; and the weights record the model for `ipw()` to rebuild the
  # numerator from, whose stacked equations are the model's unweighted score.
  dat <- categorical_numerator_dat
  dat$case_wt <- rep(c(1, 2), length.out = nrow(dat))
  # `nnet::multinom()` evaluates `weights` in the model frame it builds, so the
  # fit is written out rather than routed through the fixture helper's dots,
  # which would reach it as a promise it cannot substitute.
  weighted <- nnet::multinom(
    z ~ v,
    data = dat,
    weights = case_wt,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  expect_error(
    wt_ate(
      categorical_numerator_ps(),
      dat$z,
      exposure_type = "categorical",
      stabilize = weighted
    ),
    class = "propensity_numerator_error"
  )
})

test_that("a categorical numerator model and a stabilization score are two numerators", {
  skip_if_not_installed("nnet")

  # A score the caller computed is itself the numerator of the weights, and a
  # model is a second one; the two together are two instructions about the same
  # quantity, which is the refusal the other two exposure types make of the same
  # pair.
  expect_error(
    wt_ate(
      categorical_numerator_ps(),
      categorical_numerator_dat$z,
      exposure_type = "categorical",
      stabilize = categorical_numerator_model(),
      stabilization_score = categorical_numerator_gather()
    ),
    class = "propensity_numerator_error"
  )
})

test_that("a categorical numerator model fit to other observations is refused", {
  skip_if_not_installed("nnet")

  # A model fit to a different number of observations has one numerator row for
  # each of its own units and none for the units here, so it is refused rather
  # than recycled against them. A multinomial fit reports a row per unit and a
  # column per level, so it is the rows that are counted.
  short <- fit_categorical_numerator_model(
    z ~ v,
    utils::head(categorical_numerator_dat, 100)
  )

  err <- expect_error(
    wt_ate(
      categorical_numerator_ps(),
      categorical_numerator_dat$z,
      exposure_type = "categorical",
      stabilize = short
    ),
    class = "propensity_length_error"
  )

  expect_match(
    categorical_numerator_message(err),
    "100",
    fixed = TRUE
  )
})

test_that("a categorical numerator model stabilizes no estimand but the ATE", {
  skip_if_not_installed("nnet")

  # Only `wt_ate()` offers `stabilize` at all, so a model reaches no other
  # estimand through the exported surface. The refusal behind that is read
  # directly, since it is what keeps a tilt estimand from taking a numerator
  # whose cancellation it does not have.
  for (estimand in c("att", "atu", "atm", "ato", "entropy")) {
    expect_error(
      calculate_categorical_weights(
        categorical_numerator_ps(),
        categorical_numerator_dat$z,
        estimand,
        .focal_level = if (estimand %in% c("att", "atu")) "a",
        stabilize = TRUE,
        numerator_model = categorical_numerator_model()
      ),
      class = "propensity_stabilize_categorical_error"
    )
  }
})

# ---- combining weights that record a categorical numerator model -------------

test_that("combining a modeled categorical numerator with the marginal one drops the model", {
  skip_if_not_installed("nnet")

  # Half of the combined weights were divided by the model's fitted probability
  # and half by the marginal one, so neither numerator describes the result.
  # The record is dropped and the drop is named, rather than the model being
  # carried onto weights that were not all stabilized on it.
  modeled <- categorical_modeled_wt()
  marginal <- categorical_marginal_wt()

  out <- NULL
  expect_warning(
    out <- vec_c(modeled, marginal),
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(numerator_model(out))

  # The result is still stabilized, and on nothing it can name.
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

test_that("combining categorical weights stabilized on the same model keeps it", {
  skip_if_not_installed("nnet")

  # A numerator model is compared by what it says rather than by identity: a fit
  # carries the frame it was made in and the call that made it, so the same
  # numerator fit in two calls would otherwise read as a disagreement. A
  # multinomial fit says it in a matrix of coefficients, compared with its
  # dimnames.
  again <- fit_categorical_numerator_model(z ~ v, categorical_numerator_dat)

  same <- expect_silent(vec_c(
    categorical_modeled_wt(),
    categorical_modeled_wt(again)
  ))
  expect_identical(
    stats::coef(numerator_model(same)),
    stats::coef(categorical_numerator_model())
  )
})

test_that("combining categorical weights stabilized on different models drops both", {
  skip_if_not_installed("nnet")

  # Two models that read different terms estimate two different numerators, so
  # the combination was divided by neither.
  out <- NULL
  cnd <- expect_warning(
    out <- vec_c(
      categorical_modeled_wt(),
      categorical_modeled_wt(categorical_numerator_fit("other_num_mod"))
    ),
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(numerator_model(out))
  expect_true(
    grepl("numerator_model", conditionMessage(cnd), fixed = TRUE)
  )
})

# ---- arithmetic on weights that record a categorical numerator model ---------

test_that("categorical arithmetic that leaves the result unstabilized drops the model", {
  skip_if_not_installed("nnet")

  # A numerator model names what the weights were divided by, so a result only
  # half of which was divided by it was not. The merge already drops the
  # stabilization status, and the model is the same record of the same thing:
  # left standing, the result answers `numerator_model()` with a model and
  # prints a `stabilize:` line naming a numerator that never multiplied it.
  modeled <- categorical_modeled_wt()
  unstabilized <- categorical_unstabilized_psw()

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

test_that("categorical arithmetic that stays stabilized keeps the numerator model", {
  skip_if_not_installed("nnet")

  # The drop above belongs to the unstabilized result rather than to every
  # product: two sets of weights divided by the same numerator give a product
  # that was divided by it too.
  same <- expect_silent(categorical_modeled_wt() * categorical_modeled_wt())
  expect_true(is_stabilized(same))
  expect_identical(
    stats::coef(numerator_model(same)),
    stats::coef(categorical_numerator_model())
  )
})
