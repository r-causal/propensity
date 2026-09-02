# `ipw()` with a categorical exposure whose weights were stabilized on a fitted
# numerator model. The binary route already estimates such a model inside the
# stacked system, so the sandwich accounts for the numerator having been fitted
# rather than reading it as a constant; a categorical exposure's numerator is
# the conditional probability of the level each unit took, and the block that
# estimates it is the multinomial score of the fit that reported it. Everything
# else about the arrangement is the binary one: the numerator's coefficients sit
# in the stabilization block, in place of the k - 1 marginal proportions the
# default stabilizer seeds there.
#
# The oracle for every point estimate is the weighted marginal structural
# model's own coefficients, and the oracle for the standard error is that it
# differs from the one the same weights get when the numerator is handed over as
# a fixed `stabilization_score`: the two routes weight the units identically, so
# anything that separates them is the sandwich reading the numerator
# differently.

# ---- data simulator ---------------------------------------------------------

# A three-level exposure, a modifier that changes its effect, and a covariate
# the propensity score model reads and the numerator does not.
#
# The covariate is what the numerator has to be estimated against, for the
# reason the binary file records: it is the only thing the reported model reads
# that the numerator does not, so it is the only direction in which the
# numerator's own estimation can reach the reported standard errors. It also
# makes the numerator's design a proper subset of the propensity score model's,
# which is what keeps the coefficient block the numerator contributes from being
# the denominator's block under another name.
#
# `g` is a second three-level factor over the same levels the exposure carries.
# It is not the exposure and nothing models it, so it is what a numerator of the
# wrong response is fit to below.
sim_categorical_numerator <- function(seed = 4471, n = 500) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  v <- factor(rbinom(n, 1, 0.5))
  high <- as.numeric(v == "1")

  eta_b <- -0.2 + 1.0 * x1 - 0.7 * high
  eta_c <- 0.3 - 0.9 * x1 + 0.8 * high
  denom <- 1 + exp(eta_b) + exp(eta_c)
  pa <- 1 / denom
  pb <- exp(eta_b) / denom
  u <- runif(n)
  z <- factor(
    ifelse(u < pa, "a", ifelse(u < pa + pb, "b", "c")),
    levels = c("a", "b", "c")
  )

  y <- rbinom(
    n,
    1,
    plogis(
      -0.3 +
        0.5 * (z == "b") +
        0.9 * (z == "c") +
        0.5 * x1 -
        0.4 * high +
        0.8 * (z == "b") * high -
        0.6 * (z == "c") * high
    )
  )

  g <- factor(c("a", "b", "c")[1 + (rank(x1) %% 3)], levels = c("a", "b", "c"))

  data.frame(x1, v, z, y, g)
}

# ---- model fitting ----------------------------------------------------------

# Every multinomial fixture is fit to a tighter convergence than
# `nnet::multinom()`'s default, the convention test-ipw-categorical.R states. At
# the default tolerance the reported coefficients sit a little off the
# multinomial score root, far enough off that the solved stabilization block
# would not reproduce them and the two routes would not report the same
# estimates at the tolerances written here. The weight-consistency preflight is
# unaffected either way, `fitted()` being the softmax of the coefficients
# whatever they are, but the pins below are not.
#
# Each fit is written with a literal formula so that the frame behind it can be
# rebuilt from the fitting call, which is the only way `nnet::multinom()` keeps
# a design at all. The one fixture that deliberately cannot be rebuilt is
# written out in the test that is about it.
fit_categorical_ipw_ps <- function(dat) {
  nnet::multinom(
    z ~ x1 + v,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
}

fit_categorical_ipw_numerator <- function(dat) {
  nnet::multinom(
    z ~ v,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
}

fit_categorical_ipw_intercept <- function(dat) {
  nnet::multinom(
    z ~ 1,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
}

# A numerator of the wrong response: the same levels, fit to a column that is
# not the exposure.
fit_categorical_ipw_other <- function(dat) {
  nnet::multinom(
    g ~ v,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
}

# The numerator the model evaluates to: the fitted probability of the level each
# unit took, read out of the column named for that level rather than out of the
# column in that position.
categorical_numerator_probability <- function(model, exposure) {
  probabilities <- stats::fitted(model)
  column <- match(as.character(exposure), colnames(probabilities))

  probabilities[cbind(seq_along(exposure), column)]
}

# The names the stabilization block takes when a model fills it: the
# denominator's level-major, term-minor flattening of the coefficient matrix,
# under the prefix that says which side of the ratio the block belongs to. The
# denominator writes `<level>:<term>`, so the numerator writes
# `stab_<level>:<term>`, and the order is the order `as.vector(t(coef()))` puts
# the coefficients in.
categorical_numerator_theta_names <- function(model) {
  coefs <- stats::coef(model)
  terms <- colnames(coefs)

  as.vector(vapply(
    rownames(coefs),
    function(level) paste0("stab_", level, ":", terms),
    character(length(terms))
  ))
}

# The three fits every test below compares. All three weight the same
# multinomial propensity score model; what separates them is the numerator.
# `model` stabilizes on a fitted multinomial model of the exposure on the
# modifier, `score_fit` stabilizes on the very same numerator handed over as a
# vector of numbers, and `marginal` takes the default stabilizer, the marginal
# probability of each level. The first two build one and the same set of
# weights, so anything that differs between them is the sandwich rather than the
# weights.
categorical_numerator_fits <- function(
  dat,
  outcome_rhs = "z * v + x1",
  numerator = NULL
) {
  ps_mod <- fit_categorical_ipw_ps(dat)
  num_mod <- if (is.null(numerator)) {
    fit_categorical_ipw_numerator(dat)
  } else {
    numerator
  }
  ps <- stats::fitted(ps_mod)
  score <- categorical_numerator_probability(num_mod, dat$z)

  fit_one <- function(stabilize, stab_score) {
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_ate(
        ps,
        dat$z,
        exposure_type = "categorical",
        stabilize = stabilize,
        stabilization_score = stab_score
      )
    )
    outcome_mod <- glm(
      stats::reformulate(outcome_rhs, response = "y"),
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )

    list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
  }

  list(
    num_mod = num_mod,
    score = score,
    model = fit_one(num_mod, NULL),
    score_fit = fit_one(TRUE, score),
    marginal = fit_one(TRUE, NULL)
  )
}

# A condition's message with its line breaks and indentation flattened, so that
# a bullet written across several lines is matched as the sentence it reads as.
categorical_numerator_ipw_message <- function(cnd) {
  gsub("[[:space:]]+", " ", conditionMessage(cnd))
}

# ---- the estimates a stacked categorical numerator model reports -------------

test_that("ipw() reports the same estimate from a categorical numerator model and its score", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(dat)

  # The two routes weight the units identically, which is what makes the
  # standard errors below comparable at all.
  expect_equal(
    as.numeric(fits$model$wts),
    as.numeric(fits$score_fit$wts),
    tolerance = 1e-12
  )

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)

  expect_s3_class(res_model, "ipw")
  expect_equal(res_model$estimand, "ate")
  expect_equal(
    res_model$estimates$estimate,
    res_score$estimates$estimate,
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res_model$estimates$std.err)))
})

test_that("ipw() stacks the categorical numerator model it was given", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(dat)

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)
  res_marginal <- ipw(fits$marginal$ps_mod, fits$marginal$outcome_mod)

  # A supplied score is a constant of the stacked system and a supplied model is
  # not: the model's coefficients are parameters solved alongside everything
  # else. The default stabilizer estimates the k - 1 free marginal proportions;
  # a numerator model replaces them with its own coefficient block rather than
  # sitting beside them.
  theta_model <- coef(res_model$fit)
  theta_score <- coef(res_score$fit)
  theta_marginal <- coef(res_marginal$fit)

  stab_names <- categorical_numerator_theta_names(fits$num_mod)
  expect_identical(
    names(theta_model)[grepl("^stab_", names(theta_model))],
    stab_names
  )
  expect_false(any(grepl("^stab_", names(theta_score))))
  expect_identical(
    names(theta_marginal)[grepl("^stab_", names(theta_marginal))],
    c("stab_b", "stab_c")
  )
  expect_identical(
    length(theta_model),
    length(theta_score) + length(as.vector(t(coef(fits$num_mod))))
  )

  # The block is solved at the model it was given rather than carried at
  # whatever the seed was, and it is laid out the way the denominator's block
  # is: one run of terms per non-reference level, in the model's own level
  # order.
  expect_equal(
    unname(theta_model[stab_names]),
    as.vector(t(coef(fits$num_mod))),
    tolerance = 1e-6
  )
})

test_that("a stacked categorical numerator model reports a different standard error", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(dat)

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)

  se_model <- res_model$estimates$std.err
  se_score <- res_score$estimates$std.err

  expect_true(all(is.finite(se_model) & se_model > 0))

  # The weights are the same weights, so the difference is the numerator having
  # been estimated rather than known: the score route treats it as a constant
  # and the model route carries its estimation into the sandwich.
  expect_false(isTRUE(all.equal(se_model, se_score, tolerance = 1e-6)))
})

# ---- the weights the stacked system rebuilds --------------------------------

test_that("a categorical numerator model rebuilds the weights it was given", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(dat)

  # Everything the stacked system reports is a function of the weights it
  # rebuilds at each value of theta, so the preflight comparing the rebuild at
  # the seed against the weights the outcome model was fit with is what says the
  # two are the same weights. A numerator model seeds the stabilization block at
  # its own coefficients, and the softmax of those coefficients over its design
  # is the fitted matrix the weights were gathered from, so the rebuild is exact
  # rather than close.
  spec <- ipw_spec_categorical(fits$model$ps_mod, fits$model$outcome_mod)
  layout <- ipw_theta_layout(spec)

  expect_equal(
    ipw_weights_at_init(spec, layout),
    as.numeric(fits$model$wts),
    tolerance = 1e-10
  )
  expect_no_error(
    ipw_check_weight_consistency(spec, as.double(fits$model$wts))
  )
})

test_that("ipw_categorical_stab_probs() reads a marginal block and a modeled one", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  num_mod <- fit_categorical_ipw_numerator(dat)

  # The numerator the psi rows divide by and the numerator the preflight
  # rebuilds are the same quantity, read from the same block, so they are read
  # by one function. Without a model that is the marginal construction the two
  # wrote separately: the k - 1 free proportions with the reference level's
  # completing them.
  th_stab <- c(0.27, 0.51)
  expect_null(ipw_categorical_stab_probs(list(model = NULL), numeric(0), 3L))
  expect_equal(
    ipw_categorical_stab_probs(list(model = NULL), th_stab, 3L),
    c(1 - sum(th_stab), th_stab)
  )

  # With a model it is the n-by-K softmax over the model's own design, which at
  # the model's own coefficients is the matrix `fitted()` reports, columns in
  # the fit's level order.
  coefs <- as.vector(t(coef(num_mod)))
  block <- list(
    model = list(X = model.matrix(num_mod), coefs = coefs, k = 3L)
  )

  expect_equal(
    unname(ipw_categorical_stab_probs(block, coefs, 3L)),
    unname(stats::fitted(num_mod)),
    tolerance = 1e-10
  )
})

test_that("the categorical weight registry takes a numerator matrix as well as a row", {
  # A marginal numerator is one probability per level and a modeled one is a
  # probability per unit per level, so the registry that gathers the numerator
  # at the level each unit took has to read both shapes. A row is broadcast down
  # the units, and a matrix already has one row per unit and is read as it
  # stands. Broadcasting a matrix instead transposes it here, since the two are
  # square, and every unit is gathered at another unit's numerator.
  weight_fn <- ipw_weight_fn("categorical", "ate")

  ps <- rbind(
    c(0.5, 0.3, 0.2),
    c(0.2, 0.5, 0.3),
    c(0.6, 0.1, 0.3)
  )
  exposure <- rbind(
    c(0, 1, 0),
    c(0, 0, 1),
    c(1, 0, 0)
  )

  row_numerator <- c(0.25, 0.35, 0.40)
  expect_equal(
    weight_fn(ps, exposure, list(stab_probs = row_numerator)),
    c(0.35 / 0.3, 0.40 / 0.3, 0.25 / 0.6)
  )

  unit_numerator <- rbind(
    c(0.20, 0.50, 0.30),
    c(0.40, 0.25, 0.35),
    c(0.15, 0.45, 0.40)
  )
  expect_equal(
    weight_fn(ps, exposure, list(stab_probs = unit_numerator)),
    c(0.50 / 0.3, 0.35 / 0.3, 0.15 / 0.6)
  )
})

# ---- effect modification end to end -----------------------------------------

test_that("a categorical numerator model is the numerator `.by` asks for", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(dat)

  # A fit reporting effects within the strata of a modifier reports a default
  # stabilizer that conditions on nothing, because a numerator conditioning on
  # the modifier would be consistent for the same stratum effects and tighter.
  # This numerator model is a model of the exposure on that modifier, which is
  # the numerator the report asks for, so the model route is silent.
  expect_warning(
    ipw(fits$marginal$ps_mod, fits$marginal$outcome_mod, .by = v),
    class = "propensity_ipw_by_stabilizer_warning"
  )
  expect_no_warning(
    res <- ipw(fits$model$ps_mod, fits$model$outcome_mod, .by = v),
    class = "propensity_ipw_by_stabilizer_warning"
  )

  expect_s3_class(res, "ipw")
  expect_true("group" %in% names(res$estimates))
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("the stratum effects a categorical numerator model reports are the score's", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(dat)

  # The two routes build the same weights, so they fit the same outcome model
  # and report the same stratum effects and the same effect-modification rows.
  # Only the standard errors move, and they move because the numerator was
  # estimated rather than assumed.
  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod, .by = v)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod, .by = v)

  expect_identical(res_model$estimates$group, res_score$estimates$group)
  expect_equal(
    res_model$estimates$estimate,
    res_score$estimates$estimate,
    tolerance = 1e-8
  )
  expect_false(isTRUE(all.equal(
    res_model$estimates$std.err,
    res_score$estimates$std.err,
    tolerance = 1e-6
  )))
})

# ---- what a saturated model reads a numerator as ----------------------------

test_that("a marginal structural model saturated in the numerator reads it as none", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(dat, outcome_rhs = "z * v")

  # The binary file's companion, over three levels. A numerator of the exposure
  # on the modifier takes one value in each cell of the modifier and the
  # exposure, so a model saturated in those cells fits the same coefficients
  # with it and without it: a constant within a cell divides out of that cell's
  # weighted mean. The estimator does not move when the numerator does, so
  # nothing about having estimated the numerator can reach the standard error,
  # and the model route reports exactly what the score route reports.
  unstabilized <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      stats::fitted(fits$model$ps_mod),
      dat$z,
      exposure_type = "categorical",
      stabilize = FALSE
    )
  )
  saturated <- glm(
    y ~ z * v,
    data = dat,
    family = quasibinomial(),
    weights = unstabilized,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  expect_equal(
    unname(coef(fits$model$outcome_mod)),
    unname(coef(saturated)),
    tolerance = 1e-8
  )

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)

  # The tolerance is the one the block above uses to say the two routes report
  # different standard errors on an unsaturated model, so the same number
  # separates the two readings rather than each having its own.
  expect_equal(
    res_model$estimates$std.err,
    res_score$estimates$std.err,
    tolerance = 1e-6
  )
})

# ---- the numerator that conditions on nothing -------------------------------
#
# A numerator model of the exposure on an intercept alone estimates the marginal
# probability of each level, which is what the default stabilizer estimates. The
# two routes stack the same information in two parameterizations of it, the
# logits of the level probabilities against the free probabilities themselves,
# so they report the same estimates and the same standard errors. What separates
# them is what the fit is asked about: the model route records a model, names
# its coefficients in the solved theta, and reads its terms where the default
# stabilizer has none to read.

test_that("an intercept-only categorical numerator model reports the marginal stabilizer's answer", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(
    dat,
    numerator = fit_categorical_ipw_intercept(dat)
  )

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_marginal <- ipw(fits$marginal$ps_mod, fits$marginal$outcome_mod)

  # A smooth reparameterization of one block leaves the rest of the system
  # where it was, so this is an equality rather than a comparison. What it is
  # written at is where the two fixtures' weights part: the intercept-only fit
  # reaches the sample proportions to its own convergence rather than exactly,
  # which is 2.7e-8 relative on this fixture and reaches the outcome model's
  # coefficients and the standard errors built from them.
  expect_equal(
    res_model$estimates$estimate,
    res_marginal$estimates$estimate,
    tolerance = 1e-6
  )
  expect_equal(
    res_model$estimates$std.err,
    res_marginal$estimates$std.err,
    tolerance = 1e-6
  )
})

test_that("an intercept-only categorical numerator model is stacked as a model", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  num_mod <- fit_categorical_ipw_intercept(dat)
  fits <- categorical_numerator_fits(dat, numerator = num_mod)

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  theta <- coef(res_model$fit)

  # The block is the model's, named the way every numerator model's block is
  # named, so the intercept-only case is the general one with a design of one
  # column rather than a route back to the marginal proportions.
  expect_identical(
    names(theta)[grepl("^stab_", names(theta))],
    categorical_numerator_theta_names(num_mod)
  )
  expect_equal(
    unname(theta[categorical_numerator_theta_names(num_mod)]),
    as.vector(t(coef(num_mod))),
    tolerance = 1e-6
  )
})

test_that("an intercept-only categorical numerator model is the marginal one in a `.by` report", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(
    dat,
    numerator = fit_categorical_ipw_intercept(dat)
  )

  # The report reads the model's terms rather than its class, so a numerator
  # that reads nothing is reported as the numerator that conditions on nothing,
  # which is what it is. The sentence says so rather than naming an empty set of
  # variables.
  cnd <- expect_warning(
    ipw(fits$model$ps_mod, fits$model$outcome_mod, .by = v),
    class = "propensity_ipw_by_stabilizer_warning"
  )
  expect_match(
    categorical_numerator_ipw_message(cnd),
    "intercept alone",
    fixed = TRUE
  )
})

# ---- what a stacked categorical numerator model has to be -------------------

test_that("a categorical numerator model fit to the levels in another order is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  ps_mod <- fit_categorical_ipw_ps(dat)

  # Everything the block does is positional: the coefficients are flattened in
  # the fit's level order, the indicator matrix the score reads is laid out in
  # the denominator's, and the softmax the numerator is rebuilt from reads the
  # first level as the reference. A numerator fit to the same levels in another
  # order is therefore a different parameterization rather than a permutation of
  # the same one, and reordering its columns cannot repair it: releveling a
  # multinomial fit moves the reference level and with it every coefficient.
  #
  # The weight side gathers by name and takes this fit, so the refusal has to be
  # made here, and it has to name the cause: without it the mismatch surfaces
  # as the generic weight-consistency refusal, which is about how the weights
  # were built rather than about the order the numerator was fit in.
  dat_relevel <- dat
  dat_relevel$z <- factor(as.character(dat$z), levels = c("c", "a", "b"))
  relevel_mod <- fit_categorical_ipw_numerator(dat_relevel)

  expect_identical(relevel_mod$lev, c("c", "a", "b"))
  expect_identical(ps_mod$lev, c("a", "b", "c"))

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      stats::fitted(ps_mod),
      dat$z,
      exposure_type = "categorical",
      stabilize = relevel_mod
    )
  )
  outcome_mod <- glm(
    y ~ z * v + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(ps_mod, outcome_mod),
    class = "propensity_ipw_numerator_error"
  )

  message <- categorical_numerator_ipw_message(err)
  expect_match(message, "level order", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})

test_that("a categorical numerator model of another response is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  ps_mod <- fit_categorical_ipw_ps(dat)

  # The block's equations are the score of a model of the exposure, so a model
  # of something else sits away from the root of the rows seeded for it and the
  # solve would move it, reporting a numerator nobody fit. The levels agree, so
  # the weight side gathers a probability for every unit and has nothing to
  # object to; what is wrong is which variable those probabilities describe.
  other_mod <- fit_categorical_ipw_other(dat)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      stats::fitted(ps_mod),
      dat$z,
      exposure_type = "categorical",
      stabilize = other_mod
    )
  )
  outcome_mod <- glm(
    y ~ z * v + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(ps_mod, outcome_mod),
    class = "propensity_ipw_numerator_error"
  )

  message <- categorical_numerator_ipw_message(err)
  expect_match(message, "\"g\"", fixed = TRUE)
  expect_match(message, "\"z\"", fixed = TRUE)
})

test_that("a categorical numerator model fit with case weights is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  ps_mod <- fit_categorical_ipw_ps(dat)

  # The rows written for the numerator are its unweighted multinomial score, so
  # a fit made under case weights is not at their root. The weight layer already
  # refuses such a model, which is why the fit reaches the estimator here as a
  # record swapped onto weights the unweighted fit built: the guard is the one
  # `ipw()` owes a caller who assembled the weights themselves, and it is the
  # guard the denominator meets, read for the argument the numerator arrived in.
  dat$case_wt <- rep(c(1, 2), length.out = nrow(dat))
  # `nnet::multinom()` evaluates `weights` in the model frame it builds, so the
  # fit is written out rather than routed through a helper, which would reach it
  # as a promise it cannot substitute.
  weighted_mod <- nnet::multinom(
    z ~ v,
    data = dat,
    weights = case_wt,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  fits <- categorical_numerator_fits(dat)
  wts <- fits$model$wts
  attr(wts, "numerator_model") <- weighted_mod
  expect_identical(numerator_model(wts), weighted_mod)

  outcome_mod <- glm(
    y ~ z * v + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(ps_mod, outcome_mod),
    class = "propensity_ipw_ps_weights_error"
  )
  expect_match(
    categorical_numerator_ipw_message(err),
    "case weights",
    fixed = TRUE
  )
})

# ---- recovering the data behind a categorical numerator model ---------------
#
# `nnet::multinom()` stores no model frame, so the numerator's design is
# recovered by re-evaluating its fitting call, exactly as the denominator's is.
# A fit made inside a wrapper whose frame is gone cannot be re-evaluated, and
# what that owes the caller is the request the denominator's recovery already
# makes: name the model that cannot be rebuilt and ask for `.data`.

test_that("a categorical numerator model whose data is gone asks for .data", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  ps_mod <- fit_categorical_ipw_ps(dat)

  # The formula is written here, so its environment is this frame, and the
  # fitting call names a variable that lives only inside the function below.
  # Nothing can rebuild the numerator's design once that function has returned,
  # though its fitted probabilities are still readable and the weights were
  # built from them.
  fmla <- z ~ v
  fit_in_function <- function(fitting_data) {
    nnet::multinom(
      fmla,
      data = fitting_data,
      trace = FALSE,
      reltol = 1e-14,
      maxit = 2000
    )
  }
  gone <- fit_in_function(dat)

  expect_error(model.frame(gone))

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      stats::fitted(ps_mod),
      dat$z,
      exposure_type = "categorical",
      stabilize = gone
    )
  )
  outcome_mod <- glm(
    y ~ z * v + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(ps_mod, outcome_mod),
    class = "propensity_ipw_data_error"
  )

  # The propensity score model here is rebuildable, so a message naming it would
  # send the caller to the wrong fit.
  message <- categorical_numerator_ipw_message(err)
  expect_match(message, ".data", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})

test_that("a categorical numerator model whose data is gone is rebuilt from .data", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()

  # The other half of the contract: the `.data` the refusal above asks for
  # rebuilds the numerator's design from the fit's own terms and contrasts, and
  # what comes back is what the same numerator reports when its frame was never
  # lost.
  fmla <- z ~ v
  fit_in_function <- function(fitting_data) {
    nnet::multinom(
      fmla,
      data = fitting_data,
      trace = FALSE,
      reltol = 1e-14,
      maxit = 2000
    )
  }
  gone <- fit_in_function(dat)

  fits <- categorical_numerator_fits(dat)
  gone_fits <- categorical_numerator_fits(dat, numerator = gone)

  res_data <- ipw(
    gone_fits$model$ps_mod,
    gone_fits$model$outcome_mod,
    .data = dat
  )
  res_recovered <- ipw(fits$model$ps_mod, fits$model$outcome_mod)

  expect_s3_class(res_data, "ipw")
  expect_equal(
    res_data$estimates$estimate,
    res_recovered$estimates$estimate,
    tolerance = 1e-6
  )
  expect_equal(
    res_data$estimates$std.err,
    res_recovered$estimates$std.err,
    tolerance = 1e-6
  )
})

# ---- a numerator covariate `.data` supplies as another type ------------------
#
# A numerator model is rebuilt from `.data` the way the other models are, so a
# column it alone reads is a column the rebuild can be given as the wrong type.
# The design that comes back is then coded against different columns from the
# ones the coefficients were fit against, and the two are multiplied
# positionally. The guards the other models meet report the column and both
# types; without the numerator among them the same mistake surfaced as a raw
# error out of the seed or out of `model.matrix()`, naming neither `.data`, nor
# the column, nor the argument the model arrived in.

# A numerator of the exposure on one covariate, which nothing else reads. The
# fit is written out for the reason the fixtures above are: it is rebuilt from
# its own fitting call.
fit_categorical_ipw_numerator_x2 <- function(dat) {
  nnet::multinom(
    z ~ x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
}

test_that("a numerator covariate supplied as a factor where the fit read a number is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  dat$x2 <- round(dat$x1, 2)

  num_mod <- fit_categorical_ipw_numerator_x2(dat)
  fits <- categorical_numerator_fits(dat, numerator = num_mod)

  # A factor of three levels takes two design columns where the number it stands
  # in for took one, so the rebuilt design is wider than the coefficients it
  # would be multiplied against.
  supplied <- dat
  supplied$x2 <- dat$g

  err <- expect_error(
    muffle_coverage_warning(ipw(
      fits$model$ps_mod,
      fits$model$outcome_mod,
      .data = supplied
    )),
    class = "propensity_ipw_data_error"
  )

  message <- categorical_numerator_ipw_message(err)
  expect_match(message, "x2", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})

test_that("a numerator covariate supplied as a number where the fit read a factor is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  dat$x2 <- dat$g

  num_mod <- fit_categorical_ipw_numerator_x2(dat)
  fits <- categorical_numerator_fits(dat, numerator = num_mod)

  # The other direction, which never reaches a width to compare: the rebuild
  # asks for the fit's contrast coding of a column that arrives with no levels
  # to code.
  supplied <- dat
  supplied$x2 <- as.numeric(dat$g)

  err <- expect_error(
    muffle_coverage_warning(ipw(
      fits$model$ps_mod,
      fits$model$outcome_mod,
      .data = supplied
    )),
    class = "propensity_ipw_data_error"
  )

  message <- categorical_numerator_ipw_message(err)
  expect_match(message, "x2", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})

test_that("a numerator record that is not a multinomial fit is refused by class", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical_numerator()
  fits <- categorical_numerator_fits(dat)

  # The weight layer refuses everything but an `nnet::multinom()` for a
  # categorical exposure, so a record of another class arrives here only on
  # weights it was attached to. It has no levels to compare with the propensity
  # score model's, and a refusal reporting the levels it declares would report
  # an empty set; what is wrong is the class.
  other_mod <- glm(
    as.numeric(z == "b") ~ v,
    data = dat,
    family = binomial()
  )
  wts <- fits$model$wts
  attr(wts, "numerator_model") <- other_mod

  outcome_mod <- glm(
    y ~ z * v + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(fits$model$ps_mod, outcome_mod),
    class = "propensity_ipw_numerator_error"
  )

  message <- categorical_numerator_ipw_message(err)
  expect_match(message, "multinom", fixed = TRUE)
  expect_match(message, "glm", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})
