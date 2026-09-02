# Weights stabilized on a fitted numerator model estimate an effect conditional
# on whatever that numerator reads. A numerator of the exposure on V divides the
# weights by a quantity that varies with V, so the weighted fit identifies the
# effect within levels of V rather than the marginal one, and only a marginal
# structural model carrying V reports what those weights were built for. The
# package says so in prose; these tests are what makes `ipw()` say it.
#
# The check compares raw variables on both sides. It asks which variables the
# numerator model reads, off the right-hand side of its terms, and which
# variables the outcome model reads anywhere in its formula, and reports the
# ones the first has and the second does not. Comparing variables rather than
# terms is what lets `factor(v)`, `splines::ns(v, 3)`, and `poly(v, 2)` count as
# carrying `v`: a marginal structural model reading a transformation of V is
# conditional on V, whatever shape the term takes. Only a truly absent variable
# is reported.
#
# The numerator model's own response is not among the variables asked about. It
# is the exposure, which the outcome model carries by construction, and a report
# naming it would fire on every stabilized fit there is.
#
# It is a warning rather than a refusal. The estimate a covered numerator would
# have produced is not the one this fit reports, but the fit is a real weighted
# estimator of something and the numbers are the caller's to keep, so `ipw()`
# reports what it was given and returns.

coverage_warning_class <- "propensity_ipw_stabilizer_coverage_warning"

# Runs `expr`, collecting every coverage report it raises and muffling each one,
# and returns the value beside the messages. Muffling is what keeps the reports
# out of the test run; collecting them is what lets a test assert how many
# arrived, which is the difference between one report per call and one report
# per missing variable.
#
# Each message is squashed to single spaces because cli wraps the sentence it
# formats at the console width, so a fragment that a message contains can be
# split across lines in the string the condition carries.
catch_coverage_warnings <- function(expr) {
  messages <- character()

  value <- withCallingHandlers(
    expr,
    propensity_ipw_stabilizer_coverage_warning = function(cnd) {
      messages <<- c(messages, gsub("\\s+", " ", conditionMessage(cnd)))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, messages = messages)
}

# ---- binary ------------------------------------------------------------------

# A binary exposure, a confounder both models read, a continuous baseline
# covariate the numerator can be fit to, and a three-valued one it can be fit to
# instead. `v` is continuous so that a marginal structural model can carry it
# through a spline or a polynomial; `k` takes three values so that one can carry
# it through `factor()` without fitting a coefficient per observation.
sim_coverage_binary <- function(seed = 8815, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  v <- rnorm(n)
  k <- sample(0:2, n, replace = TRUE)
  z <- rbinom(n, 1, plogis(1.0 * x1 - 0.7 * v + 0.3 * k))
  y <- rbinom(n, 1, plogis(-0.3 + 0.7 * z + 0.5 * x1 - 0.4 * v + 0.2 * k))

  data.frame(x1, v, k, z, y)
}

# The propensity score model, a numerator of the caller's choosing, and a
# weighted marginal structural model whose right-hand side is the thing every
# test below varies. The propensity score model reads everything, so the same
# denominator serves whichever numerator is asked for.
coverage_binary_fits <- function(
  dat,
  num_rhs = "v",
  outcome_rhs = "z + x1",
  stabilize = NULL,
  stab_score = NULL
) {
  ps_mod <- glm(z ~ x1 + v + k, data = dat, family = binomial())
  num_mod <- if (is.null(stabilize)) {
    glm(
      stats::reformulate(num_rhs, response = "z"),
      data = dat,
      family = binomial()
    )
  } else {
    stabilize
  }

  wts <- wt_ate(
    ps_mod,
    stabilize = num_mod,
    stabilization_score = stab_score
  )
  outcome_mod <- glm(
    stats::reformulate(outcome_rhs, response = "y"),
    data = dat,
    family = quasibinomial(),
    weights = wts
  )

  list(ps_mod = ps_mod, num_mod = num_mod, outcome_mod = outcome_mod, wts = wts)
}

test_that("ipw() reports a binary numerator variable the outcome model omits", {
  dat <- sim_coverage_binary()
  fits <- coverage_binary_fits(dat)

  caught <- catch_coverage_warnings(ipw(fits$ps_mod, fits$outcome_mod))

  # One report for the call, not one per model the call reads.
  expect_length(caught$messages, 1L)

  # The report names the variable that is missing, the argument the numerator
  # arrived in, and the model it is missing from. A reader who is told only that
  # a numerator is uncovered has to work out which of those three to change.
  expect_match(caught$messages, "\\bv\\b")
  expect_match(caught$messages, "stabilize")
  expect_match(caught$messages, "outcome_mod")

  # A warning rather than a refusal: the fit returns, whole.
  expect_s3_class(caught$value, "ipw")
  expect_true(all(is.finite(caught$value$estimates$estimate)))
  expect_true(all(is.finite(caught$value$estimates$std.err)))
})

test_that("the coverage report is what suppressing it leaves untouched", {
  dat <- sim_coverage_binary()
  fits <- coverage_binary_fits(dat)

  caught <- catch_coverage_warnings(ipw(fits$ps_mod, fits$outcome_mod))
  suppressed <- suppressWarnings(ipw(fits$ps_mod, fits$outcome_mod))

  # The report describes the fit rather than changing it, so the estimates it
  # arrives with are the estimates it does not.
  expect_equal(
    caught$value$estimates$estimate,
    suppressed$estimates$estimate,
    tolerance = 1e-12
  )
  expect_equal(
    caught$value$estimates$std.err,
    suppressed$estimates$std.err,
    tolerance = 1e-12
  )
})

test_that("one report names every numerator variable the outcome model omits", {
  dat <- sim_coverage_binary()
  fits <- coverage_binary_fits(dat, num_rhs = c("v", "k"))

  caught <- catch_coverage_warnings(ipw(fits$ps_mod, fits$outcome_mod))

  # Two variables are uncovered and one call was made, so one report arrives and
  # it carries both names. A reader told about `v` alone would add it, refit,
  # and be told about `k`.
  expect_length(caught$messages, 1L)
  expect_match(caught$messages, "\\bv\\b")
  expect_match(caught$messages, "\\bk\\b")
})

test_that("ipw() says nothing about a numerator variable the outcome model reads", {
  dat <- sim_coverage_binary()
  fits <- coverage_binary_fits(dat, outcome_rhs = "z + x1 + v")

  expect_no_warning(
    res <- ipw(fits$ps_mod, fits$outcome_mod),
    class = coverage_warning_class
  )
  expect_s3_class(res, "ipw")
})

test_that("a transformed term carries the variable inside it", {
  dat <- sim_coverage_binary()

  # A marginal structural model conditional on a spline of V is conditional on
  # V. What the check asks is whether the outcome model reads the variable at
  # all, and `all.vars()` reaches past the call to the names inside it, so every
  # one of these is covered.
  spline_fit <- coverage_binary_fits(
    dat,
    outcome_rhs = c("z", "x1", "splines::ns(v, 3)")
  )
  expect_no_warning(
    ipw(spline_fit$ps_mod, spline_fit$outcome_mod),
    class = coverage_warning_class
  )

  poly_fit <- coverage_binary_fits(
    dat,
    outcome_rhs = c("z", "x1", "poly(v, 2)")
  )
  expect_no_warning(
    ipw(poly_fit$ps_mod, poly_fit$outcome_mod),
    class = coverage_warning_class
  )

  factor_fit <- coverage_binary_fits(
    dat,
    num_rhs = "k",
    outcome_rhs = c("z", "x1", "factor(k)")
  )
  expect_no_warning(
    ipw(factor_fit$ps_mod, factor_fit$outcome_mod),
    class = coverage_warning_class
  )
})

test_that("the numerator's own response is not a variable it has to cover", {
  dat <- sim_coverage_binary()

  # The numerator is a model of the exposure, so the exposure is among the
  # variables its formula names. It is not among the ones it conditions on, and
  # the outcome model carries it by construction, so a check reading the whole
  # formula rather than its right-hand side would report nothing that could be
  # acted on.
  fits <- coverage_binary_fits(
    dat,
    num_rhs = "1",
    outcome_rhs = "z + x1 + v + k"
  )

  expect_no_warning(
    ipw(fits$ps_mod, fits$outcome_mod),
    class = coverage_warning_class
  )
})

test_that("only a fitted numerator model is asked about its coverage", {
  dat <- sim_coverage_binary()
  ps_mod <- glm(z ~ x1 + v + k, data = dat, family = binomial())

  # The default stabilizer is the marginal probability of the exposure. It
  # conditions on nothing, so there is nothing for a marginal structural model
  # to be missing, and an outcome model reading neither covariate is covered.
  marginal_wts <- wt_ate(ps_mod, stabilize = TRUE)
  marginal_out <- glm(
    y ~ z,
    data = dat,
    family = quasibinomial(),
    weights = marginal_wts
  )
  expect_no_warning(
    ipw(ps_mod, marginal_out),
    class = coverage_warning_class
  )

  # A numerator handed over as a vector of numbers names no variables at all.
  # What it conditions on is the caller's to know, and a report asserting that a
  # variable is absent from something that reads none would be a statement about
  # the vector rather than about the estimand.
  num_mod <- glm(z ~ v, data = dat, family = binomial())
  p <- as.numeric(fitted(num_mod))
  score_wts <- wt_ate(
    ps_mod,
    stabilize = TRUE,
    stabilization_score = dat$z * p + (1 - dat$z) * (1 - p)
  )
  score_out <- glm(
    y ~ z + x1,
    data = dat,
    family = quasibinomial(),
    weights = score_wts
  )
  expect_no_warning(
    ipw(ps_mod, score_out),
    class = coverage_warning_class
  )

  # Unstabilized weights have no numerator to be uncovered.
  plain_wts <- wt_ate(ps_mod)
  plain_out <- glm(
    y ~ z + x1,
    data = dat,
    family = quasibinomial(),
    weights = plain_wts
  )
  expect_no_warning(
    ipw(ps_mod, plain_out),
    class = coverage_warning_class
  )
})

# ---- categorical -------------------------------------------------------------

# The binary fixture with a three-level exposure. The numerator is a multinomial
# fit of the same shape as the denominator, which is the numerator this branch
# added to the stacked system, and the coverage question it raises is the binary
# one: a numerator reading V makes the weighted fit conditional on V whatever
# the exposure's arity.
sim_coverage_categorical <- function(seed = 8816, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  v <- rnorm(n)
  eta_b <- -0.2 + 0.9 * x1 - 0.6 * v
  eta_c <- 0.3 - 0.8 * x1 + 0.7 * v
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
    plogis(-0.3 + 0.5 * (z == "b") + 0.9 * (z == "c") + 0.5 * x1 - 0.4 * v)
  )

  data.frame(x1, v, z, y)
}

# Every multinomial fixture is fit to a tighter convergence than
# `nnet::multinom()`'s default, the convention test-ipw-categorical.R states.
coverage_categorical_fits <- function(
  dat,
  num_rhs = "v",
  outcome_rhs = "z + x1"
) {
  ps_mod <- nnet::multinom(
    z ~ x1 + v,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  num_mod <- nnet::multinom(
    stats::reformulate(num_rhs, response = "z"),
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  wts <- wt_ate(
    stats::fitted(ps_mod),
    dat$z,
    exposure_type = "categorical",
    stabilize = num_mod
  )
  outcome_mod <- glm(
    stats::reformulate(outcome_rhs, response = "y"),
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(ps_mod = ps_mod, num_mod = num_mod, outcome_mod = outcome_mod, wts = wts)
}

test_that("ipw() reports a categorical numerator variable the outcome model omits", {
  skip_if_not_installed("nnet")

  dat <- sim_coverage_categorical()
  fits <- coverage_categorical_fits(dat)

  caught <- catch_coverage_warnings(
    ipw(fits$ps_mod, fits$outcome_mod, .data = dat)
  )

  expect_length(caught$messages, 1L)
  expect_match(caught$messages, "\\bv\\b")
  expect_match(caught$messages, "stabilize")
  expect_match(caught$messages, "outcome_mod")

  expect_s3_class(caught$value, "ipw")
  expect_true(all(is.finite(caught$value$estimates$estimate)))
  expect_true(all(is.finite(caught$value$estimates$std.err)))
})

test_that("a categorical numerator the outcome model reads is not reported", {
  skip_if_not_installed("nnet")

  dat <- sim_coverage_categorical()
  fits <- coverage_categorical_fits(dat, outcome_rhs = "z + x1 + v")

  expect_no_warning(
    res <- ipw(fits$ps_mod, fits$outcome_mod, .data = dat),
    class = coverage_warning_class
  )
  expect_s3_class(res, "ipw")
})

# ---- continuous --------------------------------------------------------------

# A dose, a confounder the propensity score model reads, and a baseline
# covariate the numerator model reads. A continuous exposure is the case
# stabilization is all but required for, so it is also the case where a
# numerator conditioning on a covariate the marginal structural model omits is
# easiest to build without noticing.
sim_coverage_continuous <- function(seed = 8817, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  v <- rnorm(n)
  A <- 0.5 + 0.8 * x1 - 0.4 * v + rnorm(n)
  yc <- 1 + 0.6 * A + 0.5 * x1 - 0.3 * v + rnorm(n)

  data.frame(x1, v, A, yc)
}

coverage_continuous_fits <- function(
  dat,
  num_rhs = "v",
  msm_rhs = "A"
) {
  ps_mod <- lm(A ~ x1, data = dat)
  num_mod <- lm(stats::reformulate(num_rhs, response = "A"), data = dat)
  wts <- wt_ate(
    as.double(fitted(ps_mod)),
    dat$A,
    exposure_type = "continuous",
    stabilize = num_mod
  )
  outcome_mod <- lm(
    stats::reformulate(msm_rhs, response = "yc"),
    data = dat,
    weights = wts
  )

  list(ps_mod = ps_mod, num_mod = num_mod, outcome_mod = outcome_mod, wts = wts)
}

test_that("ipw() reports a continuous numerator variable the outcome model omits", {
  dat <- sim_coverage_continuous()
  fits <- coverage_continuous_fits(dat)

  caught <- catch_coverage_warnings(ipw(fits$ps_mod, fits$outcome_mod))

  expect_length(caught$messages, 1L)
  expect_match(caught$messages, "\\bv\\b")
  expect_match(caught$messages, "stabilize")
  expect_match(caught$messages, "outcome_mod")

  # The estimate a continuous marginal structural model reports under an
  # identity link is the weighted fit's own coefficient, so this is the oracle
  # that says estimation ran rather than merely returned.
  expect_s3_class(caught$value, "ipw")
  expect_equal(
    caught$value$estimates$estimate,
    unname(coef(fits$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(caught$value$estimates$std.err)))
})

test_that("a continuous numerator the marginal structural model reads is not reported", {
  dat <- sim_coverage_continuous()
  fits <- coverage_continuous_fits(dat, msm_rhs = c("A", "v"))

  expect_no_warning(
    res <- ipw(fits$ps_mod, fits$outcome_mod),
    class = coverage_warning_class
  )
  expect_s3_class(res, "ipw")
})

# ---- joint -------------------------------------------------------------------

# A binary treatment and a dose intervened on together. Either component's
# weights may carry a fitted numerator, and each numerator makes the product
# conditional on what it reads, so the coverage question is asked of the
# components rather than of the product.
sim_coverage_joint <- function(seed = 8818, n = 700) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  a <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * x2))
  e <- 0.4 + 0.5 * x1 - 0.3 * x2 - 0.8 * a + rnorm(n)
  y <- 1 + 0.7 * a + 0.5 * e + 0.6 * a * e + rnorm(n)

  data.frame(x1, x2, a, e, y)
}

coverage_joint_fits <- function(
  dat,
  a_stabilize = NULL,
  dose_stabilize = TRUE,
  outcome_rhs = "a * e"
) {
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- lm(e ~ a + x1 + x2, data = dat)

  wts <- wt_joint(
    wt_ate(ps_a, stabilize = a_stabilize),
    wt_ate(
      as.double(fitted(ps_e)),
      dat$e,
      exposure_type = "continuous",
      stabilize = dose_stabilize
    ),
    exposure_type = c("binary", "continuous")
  )
  outcome_mod <- lm(
    stats::reformulate(outcome_rhs, response = "y"),
    data = dat,
    weights = wts
  )

  list(
    models = joint_wt_models(a = ps_a, e = ps_e),
    outcome_mod = outcome_mod,
    wts = wts
  )
}

test_that("ipw() reports a joint dose numerator the outcome model omits", {
  dat <- sim_coverage_joint()
  fits <- coverage_joint_fits(dat, dose_stabilize = lm(e ~ x2, data = dat))

  caught <- catch_coverage_warnings(ipw(fits$models, fits$outcome_mod))

  expect_length(caught$messages, 1L)
  expect_match(caught$messages, "\\bx2\\b")
  expect_match(caught$messages, "stabilize")
  expect_match(caught$messages, "outcome_mod")

  expect_s3_class(caught$value, "ipw")
  expect_true(all(is.finite(caught$value$estimates$estimate)))
  expect_true(all(is.finite(caught$value$estimates$std.err)))
})

test_that("ipw() reports a joint discrete numerator the outcome model omits", {
  dat <- sim_coverage_joint()
  fits <- coverage_joint_fits(
    dat,
    a_stabilize = glm(a ~ x2, data = dat, family = binomial())
  )

  caught <- catch_coverage_warnings(ipw(fits$models, fits$outcome_mod))

  # A numerator on either component reaches the same product, so a component
  # whose numerator reads a variable the marginal structural model omits is
  # reported wherever it sits.
  expect_length(caught$messages, 1L)
  expect_match(caught$messages, "\\bx2\\b")

  expect_s3_class(caught$value, "ipw")
  expect_true(all(is.finite(caught$value$estimates$estimate)))
})

test_that("a joint numerator the outcome model reads is not reported", {
  dat <- sim_coverage_joint()
  fits <- coverage_joint_fits(
    dat,
    dose_stabilize = lm(e ~ x2, data = dat),
    outcome_rhs = c("a * e", "x2")
  )

  expect_no_warning(
    res <- ipw(fits$models, fits$outcome_mod),
    class = coverage_warning_class
  )
  expect_s3_class(res, "ipw")
})

test_that("the joint route's default numerators are not asked about coverage", {
  dat <- sim_coverage_joint()

  # The dose takes the default stabilizer and the binary component takes none,
  # so neither carries a model and the marginal structural model has nothing to
  # be missing.
  fits <- coverage_joint_fits(dat)

  expect_no_warning(
    res <- ipw(fits$models, fits$outcome_mod),
    class = coverage_warning_class
  )
  expect_s3_class(res, "ipw")
})
