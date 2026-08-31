# `ipw()` with a binary exposure whose weights were stabilized on a fitted
# numerator model. The continuous route already estimates such a model in the
# stacked system, so the sandwich accounts for the numerator having been fitted
# rather than reading it as a constant. A binary exposure's numerator is a
# conditional probability rather than a conditional density, and the block that
# estimates it is the score of the binomial fit that reported it, but everything
# else about the arrangement is the same: the numerator parameters sit in the
# stabilization block, in place of the single marginal proportion the default
# stabilizer seeds there.
#
# The oracle for every point estimate is the weighted marginal structural model's
# own coefficients, and the oracle for the standard error is that it differs from
# the one the same weights get when the numerator is handed over as a fixed
# `stabilization_score`: the two routes weight the units identically, so anything
# that separates them is the sandwich reading the numerator differently.

# ---- data simulator ---------------------------------------------------------

# A binary exposure, a modifier that changes its effect, and a covariate the
# propensity score model reads and the numerator does not.
sim_binary_numerator <- function(seed = 6621, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  v <- factor(rbinom(n, 1, 0.5))
  high <- as.numeric(v == "1")
  z <- rbinom(n, 1, plogis(0.5 * x1 - 0.6 * high))
  y <- rbinom(
    n,
    1,
    plogis(-0.3 + 0.7 * z + 0.5 * x1 - 0.4 * high + 0.9 * z * high)
  )
  data.frame(x1, v, z, y)
}

# ---- model fitting ----------------------------------------------------------

# The three fits every test below compares. All three weight the same propensity
# score model; what separates them is the numerator. `model` stabilizes on a
# fitted binomial model of the exposure on the modifier, `score` stabilizes on
# the very same numerator handed over as a vector of numbers, and `marginal`
# takes the default stabilizer. The first two build one and the same set of
# weights, so anything that differs between them is the sandwich rather than the
# weights.
binary_numerator_fits <- function(dat, outcome_rhs = "z * v") {
  ps_mod <- glm(z ~ x1 + v, data = dat, family = binomial())
  num_mod <- glm(z ~ v, data = dat, family = binomial())
  p <- as.numeric(fitted(num_mod))
  score <- dat$z * p + (1 - dat$z) * (1 - p)

  fit_one <- function(stabilize, stab_score) {
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_ate(
        ps_mod,
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

# ---- the estimates a stacked binary numerator model reports -----------------

test_that("ipw() reports the same estimate from a binary numerator model and its score", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

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
  expect_equal(
    res_model$estimates$estimate,
    res_score$estimates$estimate,
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res_model$estimates$std.err)))
})

test_that("ipw() stacks the binary numerator model it was given", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

  res_model <- ipw(fits$model$ps_mod, fits$model$outcome_mod)
  res_score <- ipw(fits$score_fit$ps_mod, fits$score_fit$outcome_mod)
  res_marginal <- ipw(fits$marginal$ps_mod, fits$marginal$outcome_mod)

  # A supplied score is a constant of the stacked system and a supplied model is
  # not: the model's coefficients are parameters solved alongside everything
  # else. The default stabilizer estimates one parameter, the marginal
  # proportion; a numerator model replaces it with its own block rather than
  # sitting beside it.
  theta_model <- coef(res_model$fit)
  theta_score <- coef(res_score$fit)
  theta_marginal <- coef(res_marginal$fit)

  stab_names <- paste0("stab_", names(coef(fits$num_mod)))
  expect_true(all(stab_names %in% names(theta_model)))
  expect_false(any(grepl("^stab_", names(theta_score))))
  expect_identical(
    names(theta_marginal)[grepl("^stab_", names(theta_marginal))],
    "stab_pi"
  )
  expect_identical(
    length(theta_model),
    length(theta_score) + length(coef(fits$num_mod))
  )

  # The block is solved at the model it was given rather than carried at
  # whatever the seed was.
  expect_equal(
    unname(theta_model[stab_names]),
    unname(coef(fits$num_mod)),
    tolerance = 1e-6
  )
})

test_that("a stacked binary numerator model reports a different standard error", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

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

# ---- effect modification end to end -----------------------------------------

test_that("a binary numerator model is the numerator `.by` asks for", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

  # A fit reporting effects within the strata of a modifier reports a default
  # stabilizer that conditions on nothing, because a numerator conditioning on
  # the modifier would be consistent for the same stratum effects and tighter.
  # A numerator the caller supplied says nothing here, whether it arrived as a
  # score or as a fitted model.
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

test_that("the stratum effects a binary numerator model reports are the score's", {
  dat <- sim_binary_numerator()
  fits <- binary_numerator_fits(dat)

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
