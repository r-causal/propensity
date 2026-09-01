# The fitted models a continuous exposure can be weighted from: `lm()`, a
# gaussian `glm()` under any of its links, `mgcv::gam()`, and `MASS::rlm()`,
# together with the families each route refuses.
#
# Some of this is already true and is pinned here so that it stays true. The
# tests that pin behavior the package has today are the gaussian glm links, the
# gaussian gam that reaches the glm method by inheritance, the binomial gam and
# the quasibinomial glm on the binary path, and the refusal of a model by every
# estimand other than the ATE and the uncensored weights. Everything else is
# new: the `lm` method and the inheritance route it opens for `rlm`, the family
# checks on the continuous path, and the refusal of a linear model of a binary
# exposure.
#
# One test elsewhere changes when the `lm` method arrives, and is left alone
# here. "GLM methods error appropriately" in tests/testthat/test-weights.R asks
# for ATE weights from an `lm` and a binary exposure, and snapshots the answer
# as `propensity_method_error`. Under the contract these tests pin, that call
# becomes a refusal of a binary exposure, so the snapshot in
# tests/testthat/_snaps/weights.md has to be read and accepted along with the
# implementation rather than ahead of it.

# One continuous-exposure problem for the model methods to be fit against. The
# strictly positive exposure is there so that the log, inverse, and square root
# links have a conditional mean they can hold, and the binary one so that a
# linear model of a probability can be written down and refused.
continuous_model_data <- local({
  set.seed(20250901)

  n <- 150
  x <- rnorm(n)
  z <- rnorm(n)

  data.frame(
    x = x,
    z = z,
    dose = 1 + 0.7 * x - 0.4 * z + rnorm(n),
    positive = exp(0.35 + 0.3 * x + 0.2 * rnorm(n)),
    binary = rbinom(n, 1, plogis(0.5 * x))
  )
})

# The weights a model route is asked for, and the weights the numeric route
# gives for the same call when it is handed the model's fitted values rather
# than the model. Every model method has to agree with the second.
# Several of these fits explain more than half the variance of the exposure they
# are fit to, which puts the stabilized normal weights past the finite-variance
# boundary and raises the report that says so. These tests are written about
# which fitted values and which spread a model route reads, so the report is
# muffled on both sides of the comparison.
continuous_model_wt <- function(fit, exposure, ..., .weight_fn = wt_ate) {
  muffle_variance_warning(.weight_fn(
    fit,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    ...
  ))
}

continuous_model_oracle <- function(fit, exposure, ..., .weight_fn = wt_ate) {
  muffle_variance_warning(.weight_fn(
    as.numeric(stats::fitted(fit)),
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE,
    ...
  ))
}

# ---- linear models ----------------------------------------------------------

test_that("a linear model gives the weights its fitted values give", {
  fit <- lm(dose ~ x + z, data = continuous_model_data)

  expect_equal(
    as.numeric(continuous_model_wt(fit, continuous_model_data$dose)),
    as.numeric(continuous_model_oracle(fit, continuous_model_data$dose)),
    tolerance = 1e-12
  )
})

test_that("a linear model without an exposure reads its own response", {
  fit <- lm(dose ~ x + z, data = continuous_model_data)

  from_model <- wt_ate(fit, exposure_type = "continuous", stabilize = TRUE)

  expect_equal(
    as.numeric(from_model),
    as.numeric(continuous_model_oracle(fit, continuous_model_data$dose)),
    tolerance = 1e-12
  )
})

test_that("a linear model records a continuous exposure and a pooled spread", {
  fit <- lm(dose ~ x + z, data = continuous_model_data)

  weights <- continuous_model_wt(fit, continuous_model_data$dose)

  expect_identical(exposure_type(weights), "continuous")
  expect_identical(density_meta(weights)$sigma, "pooled")
  expect_identical(density_meta(weights)$numerator, "marginal")
})

test_that("a spread given to a linear model is recorded as supplied", {
  fit <- lm(dose ~ x + z, data = continuous_model_data)
  # A spread for each observation, laid out rather than drawn, so that the
  # weights it gives do not depend on the state of the generator.
  per_observation <- seq(0.7, 1.3, length.out = nrow(continuous_model_data))

  for (sigma in list(1.4, per_observation)) {
    weights <- continuous_model_wt(
      fit,
      continuous_model_data$dose,
      .sigma = sigma
    )

    expect_equal(
      as.numeric(weights),
      as.numeric(continuous_model_oracle(
        fit,
        continuous_model_data$dose,
        .sigma = sigma
      )),
      tolerance = 1e-12
    )
    expect_identical(density_meta(weights)$sigma, "supplied")
  }
})

test_that("a density reaches the weights through a linear model", {
  fit <- lm(dose ~ x + z, data = continuous_model_data)

  weights <- continuous_model_wt(
    fit,
    continuous_model_data$dose,
    .density = dens_t(df = 4)
  )

  expect_equal(
    as.numeric(weights),
    as.numeric(continuous_model_oracle(
      fit,
      continuous_model_data$dose,
      .density = dens_t(df = 4)
    )),
    tolerance = 1e-12
  )
  expect_identical(format(density_meta(weights)$density), "t(df = 4)")
})

test_that("an integrated numerator reaches the weights through a linear model", {
  fit <- lm(dose ~ x + z, data = continuous_model_data)

  weights <- continuous_model_wt(
    fit,
    continuous_model_data$dose,
    numerator = "integrated"
  )

  expect_equal(
    as.numeric(weights),
    as.numeric(continuous_model_oracle(
      fit,
      continuous_model_data$dose,
      numerator = "integrated"
    )),
    tolerance = 1e-12
  )
  expect_identical(density_meta(weights)$numerator, "integrated")
})

test_that("a linear model refuses a binary exposure", {
  fit <- lm(binary ~ x, data = continuous_model_data)

  expect_error(
    wt_ate(fit, continuous_model_data$binary, exposure_type = "binary"),
    class = "propensity_model_family_error"
  )

  # And when nothing names the type, so the zero-one response is what says the
  # exposure is binary.
  expect_error(
    wt_ate(fit),
    class = "propensity_model_family_error"
  )

  # A linear model of a probability is not held to the unit interval, and the
  # fitted values of this one leave it. The refusal has to read the model rather
  # than the values it fitted, so it has to come before the check on the range
  # of a propensity score, which would otherwise be all the caller hears about.
  escapes <- lm(binary ~ poly(x, 3), data = continuous_model_data)

  expect_error(
    wt_ate(escapes, continuous_model_data$binary, exposure_type = "binary"),
    class = "propensity_model_family_error"
  )
  expect_error(
    wt_ate(escapes),
    class = "propensity_model_family_error"
  )
})

test_that("a linear model refuses a categorical exposure", {
  fit <- lm(dose ~ x + z, data = continuous_model_data)
  levels <- factor(rep(
    c("a", "b", "c"),
    length.out = nrow(continuous_model_data)
  ))

  expect_error(
    wt_ate(fit, levels, exposure_type = "categorical"),
    class = "propensity_model_family_error"
  )
})

test_that("censoring weights take a linear model of their own", {
  fit <- lm(dose ~ x + z, data = continuous_model_data)

  weights <- continuous_model_wt(
    fit,
    continuous_model_data$dose,
    .weight_fn = wt_cens
  )

  expect_equal(
    as.numeric(weights),
    as.numeric(continuous_model_oracle(
      fit,
      continuous_model_data$dose,
      .weight_fn = wt_cens
    )),
    tolerance = 1e-12
  )
  expect_identical(estimand(weights), "uncensored")
  expect_identical(exposure_type(weights), "continuous")
  expect_identical(density_meta(weights)$sigma, "pooled")
})

test_that("censoring weights refuse a linear model of a binary exposure", {
  fit <- lm(binary ~ x, data = continuous_model_data)

  expect_error(
    wt_cens(fit, continuous_model_data$binary, exposure_type = "binary"),
    class = "propensity_model_family_error"
  )
  expect_error(
    wt_cens(fit),
    class = "propensity_model_family_error"
  )

  # A linear model of a probability is not held to the unit interval, and the
  # fitted values of this one leave it. The refusal has to read the model rather
  # than the values it fitted, so it has to come before the check on the range
  # of a propensity score, which would otherwise be all the caller hears about.
  escapes <- lm(binary ~ poly(x, 3), data = continuous_model_data)

  expect_error(
    wt_cens(escapes, continuous_model_data$binary, exposure_type = "binary"),
    class = "propensity_model_family_error"
  )
  expect_error(
    wt_cens(escapes),
    class = "propensity_model_family_error"
  )
})

# ---- gaussian links ---------------------------------------------------------

test_that("a gaussian model uses its fitted values under any link", {
  # Pins behavior the package already has: `fitted()` is on the scale of the
  # response whatever the link is, so the link never has to be undone.
  for (link in c("identity", "log", "inverse", "sqrt")) {
    fit <- glm(
      positive ~ x + z,
      data = continuous_model_data,
      family = gaussian(link = link)
    )

    expect_equal(
      as.numeric(continuous_model_wt(fit, continuous_model_data$positive)),
      as.numeric(continuous_model_oracle(fit, continuous_model_data$positive)),
      tolerance = 1e-12,
      info = link
    )
  }
})

test_that("a model of a continuous exposure refuses a family that is not gaussian", {
  counts <- round(continuous_model_data$positive * 5)

  poisson_fit <- glm(
    counts ~ x + z,
    data = continuous_model_data,
    family = poisson()
  )
  gamma_fit <- glm(
    positive ~ x + z,
    data = continuous_model_data,
    family = Gamma(link = "log")
  )
  binomial_fit <- glm(
    binary ~ x + z,
    data = continuous_model_data,
    family = binomial()
  )

  for (fit in list(poisson_fit, gamma_fit, binomial_fit)) {
    expect_error(
      continuous_model_wt(fit, continuous_model_data$dose),
      class = "propensity_model_family_error"
    )
  }
})

test_that("a quasi model with another variance refuses a continuous exposure", {
  # The constant variance is what makes a quasi model a model of a conditional
  # mean with one spread for every observation. Another variance function
  # describes a different spread at every fitted value, which the density these
  # weights are a ratio of has no way to take.
  fit <- glm(
    positive ~ x + z,
    data = continuous_model_data,
    family = quasi(link = "log", variance = "mu")
  )

  expect_error(
    continuous_model_wt(fit, continuous_model_data$dose),
    class = "propensity_model_family_error"
  )
})

test_that("a quasi model with the gaussian variance is accepted", {
  fit <- glm(
    dose ~ x + z,
    data = continuous_model_data,
    family = quasi(link = "identity", variance = "constant")
  )

  expect_equal(
    as.numeric(continuous_model_wt(fit, continuous_model_data$dose)),
    as.numeric(continuous_model_oracle(fit, continuous_model_data$dose)),
    tolerance = 1e-12
  )
})

# ---- additive models --------------------------------------------------------

test_that("a gaussian additive model gives the weights its fitted values give", {
  # Pins behavior the package already has: a `gam` inherits from `glm`, so it
  # reaches the same method.
  skip_if_not_installed("mgcv")

  fit <- mgcv::gam(
    dose ~ s(x) + s(z),
    data = continuous_model_data,
    family = gaussian()
  )

  weights <- continuous_model_wt(fit, continuous_model_data$dose)

  expect_equal(
    as.numeric(weights),
    as.numeric(continuous_model_oracle(fit, continuous_model_data$dose)),
    tolerance = 1e-12
  )
  expect_identical(exposure_type(weights), "continuous")
  expect_identical(density_meta(weights)$sigma, "pooled")
})

test_that("an additive model of a continuous exposure refuses a binomial family", {
  skip_if_not_installed("mgcv")

  fit <- mgcv::gam(
    binary ~ s(x),
    data = continuous_model_data,
    family = binomial()
  )

  expect_error(
    continuous_model_wt(fit, continuous_model_data$dose),
    class = "propensity_model_family_error"
  )
})

test_that("a binomial additive model weights a binary exposure as a glm does", {
  # Pins behavior the package already has.
  skip_if_not_installed("mgcv")

  fit <- mgcv::gam(
    binary ~ s(x),
    data = continuous_model_data,
    family = binomial()
  )

  expect_equal(
    as.numeric(wt_ate(
      fit,
      continuous_model_data$binary,
      exposure_type = "binary"
    )),
    as.numeric(wt_ate(
      as.numeric(stats::fitted(fit)),
      continuous_model_data$binary,
      exposure_type = "binary"
    )),
    tolerance = 1e-12
  )
})

# ---- robust linear models ---------------------------------------------------

test_that("a robust linear model is spread by the pooled residual root mean square", {
  skip_if_not_installed("MASS")

  fit <- MASS::rlm(dose ~ x + z, data = continuous_model_data)

  weights <- continuous_model_wt(fit, continuous_model_data$dose)

  expect_equal(
    as.numeric(weights),
    as.numeric(continuous_model_oracle(fit, continuous_model_data$dose)),
    tolerance = 1e-12
  )
  expect_identical(exposure_type(weights), "continuous")
  expect_identical(density_meta(weights)$sigma, "pooled")
})

test_that("the scale a robust linear model fit with is a spread of its own", {
  skip_if_not_installed("MASS")

  fit <- MASS::rlm(dose ~ x + z, data = continuous_model_data)

  # `rlm` reports a robust scale estimate, which is not the pooled root mean
  # square the weights are spread by unless it is asked for.
  from_robust_scale <- continuous_model_wt(
    fit,
    continuous_model_data$dose,
    .sigma = fit$s
  )

  expect_equal(
    as.numeric(from_robust_scale),
    as.numeric(continuous_model_oracle(
      fit,
      continuous_model_data$dose,
      .sigma = fit$s
    )),
    tolerance = 1e-12
  )
  expect_identical(density_meta(from_robust_scale)$sigma, "supplied")

  expect_false(isTRUE(all.equal(
    as.numeric(from_robust_scale),
    as.numeric(continuous_model_wt(fit, continuous_model_data$dose))
  )))
})

# ---- the binary path --------------------------------------------------------

test_that("a gaussian model refuses a binary exposure", {
  fit <- glm(binary ~ x, data = continuous_model_data, family = gaussian())

  expect_error(
    wt_ate(fit, continuous_model_data$binary, exposure_type = "binary"),
    class = "propensity_model_family_error"
  )

  # And when nothing names the type, so the zero-one response is what says the
  # exposure is binary.
  expect_error(
    wt_ate(fit),
    class = "propensity_model_family_error"
  )

  # A linear model of a probability is not held to the unit interval, and the
  # fitted values of this one leave it. The refusal has to read the model rather
  # than the values it fitted, so it has to come before the check on the range
  # of a propensity score, which would otherwise be all the caller hears about.
  escapes <- glm(
    binary ~ poly(x, 3),
    data = continuous_model_data,
    family = gaussian()
  )

  expect_error(
    wt_ate(escapes, continuous_model_data$binary, exposure_type = "binary"),
    class = "propensity_model_family_error"
  )
  expect_error(
    wt_ate(escapes),
    class = "propensity_model_family_error"
  )
})

test_that("the binomial families still weight a binary exposure", {
  # Pins behavior the package already has: quasibinomial is a model of a
  # probability, and is read as one.
  for (family in list(binomial(), quasibinomial())) {
    fit <- glm(binary ~ x, data = continuous_model_data, family = family)

    expect_equal(
      as.numeric(wt_ate(
        fit,
        continuous_model_data$binary,
        exposure_type = "binary"
      )),
      as.numeric(wt_ate(
        as.numeric(stats::fitted(fit)),
        continuous_model_data$binary,
        exposure_type = "binary"
      )),
      tolerance = 1e-12,
      info = family$family
    )
  }
})

# ---- the estimands that take no model ---------------------------------------

test_that("the other estimands take no linear model", {
  # Pins the scope of the model methods: only the ATE and the uncensored
  # weights are a ratio of densities, so only they read a model of a continuous
  # exposure.
  skip_if_not_installed("MASS")

  fits <- list(
    lm(dose ~ x + z, data = continuous_model_data),
    MASS::rlm(dose ~ x + z, data = continuous_model_data)
  )
  weight_fns <- list(wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)

  for (fit in fits) {
    for (weight_fn in weight_fns) {
      expect_error(
        weight_fn(fit, continuous_model_data$dose),
        class = "propensity_method_error"
      )
    }
  }
})

# ---- the refusals as the caller reads them ----------------------------------

test_that("a model of the wrong kind is refused by what it was fit with", {
  counts <- round(continuous_model_data$positive * 5)

  binary_lm <- lm(binary ~ x, data = continuous_model_data)
  binary_gaussian <- glm(
    binary ~ x,
    data = continuous_model_data,
    family = gaussian()
  )
  poisson_fit <- glm(
    counts ~ x + z,
    data = continuous_model_data,
    family = poisson()
  )
  quasi_fit <- glm(
    positive ~ x + z,
    data = continuous_model_data,
    family = quasi(link = "log", variance = "mu")
  )

  expect_propensity_error(
    wt_ate(binary_lm, continuous_model_data$binary, exposure_type = "binary")
  )
  expect_propensity_error(
    wt_ate(
      binary_gaussian,
      continuous_model_data$binary,
      exposure_type = "binary"
    )
  )
  expect_propensity_error(
    wt_ate(
      poisson_fit,
      continuous_model_data$dose,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )
  expect_propensity_error(
    wt_ate(
      quasi_fit,
      continuous_model_data$dose,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )

  # A family with a link and no name is described by what it is missing, rather
  # than by a pair of parentheses with nothing in front of them.
  expect_propensity_error(
    check_binary_model_family(
      structure(
        list(family = list(link = "logit")),
        class = "unnamed_family_fit"
      )
    )
  )
})

test_that("a family named with an empty string is described as unnamed", {
  # A name of no characters is no name: written into the message it leaves a
  # pair of backticks with nothing between them, which describes the family to
  # nobody. Such a family is described the way one carrying no name at all is.
  empty <- structure(
    list(family = list(family = "", link = "logit")),
    class = "empty_family_fit"
  )

  err <- expect_error(
    check_binary_model_family(empty),
    class = "propensity_model_family_error"
  )

  message <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(message, "unnamed family", fixed = TRUE)
  expect_no_match(message, "``", fixed = TRUE)

  expect_identical(model_family_label(list(family = "")), "an unnamed family")

  expect_propensity_error(check_binary_model_family(empty))
})

test_that("a family element that is not a family object is refused", {
  # `$` on an atomic vector is an error of base R's rather than one of this
  # package's, so a fit whose family element is a bare string reached the
  # family test and raised there instead of being refused by it. Neither check
  # can read such an element as a family, and both say so in their own terms.
  atomic_binary <- structure(
    list(family = "binomial"),
    class = "atomic_family_fit"
  )
  atomic_continuous <- structure(
    list(family = "gaussian"),
    class = "atomic_family_fit"
  )

  expect_error(
    check_binary_model_family(atomic_binary),
    class = "propensity_model_family_error"
  )
  expect_error(
    check_continuous_model_family(atomic_continuous),
    class = "propensity_model_family_error"
  )

  expect_propensity_error(check_binary_model_family(atomic_binary))
  expect_propensity_error(check_continuous_model_family(atomic_continuous))
})

test_that("model_family_label() writes a family the way it is called", {
  expect_identical(model_family_label(list(family = "gaussian")), "gaussian()")
  expect_identical(
    model_family_label(list(family = "quasi", varfun = "constant")),
    "quasi(variance = \"constant\")"
  )

  # An extended family names itself with the parameters it was fit at, so the
  # name is already a call and adding a pair of parentheses to it writes one
  # that could not be run.
  expect_identical(
    model_family_label(list(family = "Scaled t(Inf,1.012)")),
    "Scaled t(Inf,1.012)"
  )
})

test_that("a scaled t additive model is refused by the family it prints", {
  skip_if_not_installed("mgcv")

  fit <- mgcv::gam(
    dose ~ s(x),
    data = continuous_model_data,
    family = mgcv::scat()
  )

  expect_error(
    continuous_model_wt(fit, continuous_model_data$dose),
    class = "propensity_model_family_error"
  )

  expect_identical(model_family_label(fit$family), fit$family$family)
})
