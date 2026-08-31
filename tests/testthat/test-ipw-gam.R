# Tests for `ipw()` with an `mgcv::gam()` propensity score model of a continuous
# exposure. The stacked system reads an additive fit as penalized least squares
# with its smoothing parameters held at the values the fit chose: the design is
# `predict(type = "lpmatrix")`, the penalty is the sum of `sp[j] * S[[j]]` over
# the fitted smooths, and the equation the coefficients solve is
#
#   X'(a - X alpha) - S_lambda alpha = 0.
#
# The penalty is not a sum over observations, so it is spread evenly across the
# columns of the psi matrix; the column sums are what the equation says and the
# bread picks up `(X'X + S_lambda) / n`.
#
# Holding the smoothing parameters fixed is what makes the sandwich writable at
# all, and it is also the approximation. A simulation settled the terms before
# the route was opened: on identical data the coverage tracks the `glm` route
# the package already supports to within one to two points, and fixing the
# smoothing choice costs one to three percent of the standard error at n = 500.
# Four limitations came out of it, and the three that can be measured are pinned
# below.
#
#   1. The variance conditions on the smoothing parameters rather than
#      estimating them, so a fit handed those parameters reports what the fit
#      that chose them reports.
#   2. It is the frequentist variance and not mgcv's Bayesian `Vp`, so the
#      propensity score block reads as `vcov(fit, freq = TRUE)` rather than as
#      `vcov(fit)`.
#   3. The envelope is the one the `glm` route has: a gaussian family, an
#      identity or log link, and no prior weights.
#   4. Coverage degrades as the weights grow tails, whatever the dose model is.
#      That is a property of the weights rather than of this route, and no test
#      here can pin it.

# ---- data simulator ---------------------------------------------------------

# A dose whose conditional mean is a cubic in one covariate and a quadratic in
# the other, which is what an additive fit is reached for: smooth, nonlinear,
# and not a term any formula here writes. The scale keeps the marginal dose
# variance below twice the conditional variance, where stabilized density-ratio
# weights have a second moment at all. `Apos` is the same dose shifted to be
# positive everywhere, which a log-link fit needs to start its iteration.
sim_gam_dose <- function(seed = 11, n = 400) {
  withr::local_seed(seed)
  x1 <- runif(n, -2, 2)
  x2 <- runif(n, -2, 2)
  A <- 1 + 0.6 * (0.4 * x1^3 - 0.8 * x1 + 0.6 * x2^2 - 0.5 * x2) + rnorm(n)
  yc <- 2 + A + 0.8 * x1 + 0.5 * x2^2 + rnorm(n)
  data.frame(x1, x2, A, yc, Apos = A + 6)
}

# ---- model fitting ----------------------------------------------------------

gam_dose_weights <- function(fitted_ps, a) {
  withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(fitted_ps, a, exposure_type = "continuous", stabilize = TRUE)
  )
}

# Stabilized continuous ATE weights from a fitted dose model, and the weighted
# marginal structural model they fit. Every fixture in this file is one dose
# model and one marginal structural model of the dose alone, so the reported
# slope is the coefficient the outcome model already holds.
fit_gam_models <- function(dat, ps_mod, exposure = "A", outcome = "yc") {
  wts <- gam_dose_weights(as.double(fitted(ps_mod)), dat[[exposure]])
  outcome_mod <- lm(
    stats::reformulate(exposure, response = outcome),
    data = dat,
    weights = wts
  )

  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

# The stabilized normal ratio the weights are, written out by hand: the marginal
# density of the dose over its conditional density at the fit's own mean, both
# spread by the pooled residual root mean square.
gam_hand_weights <- function(mu, a) {
  stats::dnorm(a, mean(a), sqrt(mean((a - mean(a))^2))) /
    stats::dnorm(a, mu, sqrt(mean((a - mu)^2)))
}

# ---- assertion helpers ------------------------------------------------------

# The stacked system is seeded at the fits' own parameters, so every equation in
# it has to be at its root there. A propensity score block written without the
# penalty sits away from zero here, and the solve would walk the coefficients
# off the additive fit and report a dose model nobody fit.
expect_gam_root_seeded <- function(ps_mod, outcome_mod, wts) {
  spec <- ipw_spec_continuous(ps_mod, outcome_mod)
  layout <- ipw_theta_layout(spec)
  mat <- build_ipw_psi(spec, layout)(layout$init)

  testthat::expect_true(is.matrix(mat))
  testthat::expect_identical(dim(mat), c(length(layout$init), spec$n))
  testthat::expect_false(anyNA(mat))
  testthat::expect_true(all(abs(rowSums(mat)) / spec$n < 1e-8))
  testthat::expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(wts),
    tolerance = 1e-12
  )

  invisible(layout)
}

# ---- end to end -------------------------------------------------------------

test_that("ipw() stacks a gam dose model and reports the weighted slope", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()
  ps_mod <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, method = "REML")
  mods <- fit_gam_models(dat, ps_mod)

  # The weights are the ratio the record says they are, computed here from the
  # additive fit's conditional means without going through the package, and the
  # estimate is the slope the model those weights fit already reports.
  hand <- gam_hand_weights(as.double(fitted(ps_mod)), dat$A)
  expect_equal(as.double(mods$wts), hand, tolerance = 1e-12)
  hand_msm <- lm(yc ~ A, data = dat, weights = hand)

  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_s3_class(res, "ipw")
  expect_identical(res$se_method, "mestimation")
  expect_identical(res$estimand, "ate")
  expect_identical(nrow(res$estimates), 1L)
  expect_identical(res$estimates$effect, "slope")
  expect_equal(
    res$estimates$estimate,
    unname(coef(hand_msm)[["A"]]),
    tolerance = 1e-8
  )
  expect_true(is.finite(res$estimates$std.err))
  expect_gt(res$estimates$std.err, 0)
})

test_that("the solved system stays on the gam's own coefficients", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()
  ps_mod <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, method = "REML")
  mods <- fit_gam_models(dat, ps_mod)

  expect_gam_root_seeded(mods$ps_mod, mods$outcome_mod, mods$wts)

  res <- ipw(mods$ps_mod, mods$outcome_mod)
  theta <- coef(res$fit)
  ps_block <- theta[seq_along(coef(ps_mod))]
  expect_equal(unname(ps_block), unname(coef(ps_mod)), tolerance = 1e-8)

  # What the penalty is worth: the same basis fit by least squares alone sits
  # well away from where the smoothing put it, so a block written without the
  # penalty would be caught by the comparison above rather than agreeing with it
  # to eight places.
  basis <- stats::predict(ps_mod, type = "lpmatrix")
  unpenalized <- stats::lm.fit(x = basis, y = dat$A)$coefficients
  expect_gt(max(abs(coef(ps_mod) - unpenalized)), 0.05)
})

test_that("a gam dose model reaches the same fit through .data", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()
  ps_mod <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, method = "REML")
  mods <- fit_gam_models(dat, ps_mod)

  # An additive fit's design is a smooth basis rather than the columns its
  # formula names, so a `.data` rebuild has to ask the fit for the basis at the
  # rows it was handed rather than re-evaluating the formula. Both routes read
  # the same rows, so they are the same fit.
  expect_equal(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat)$estimates,
    ipw(mods$ps_mod, mods$outcome_mod)$estimates,
    tolerance = 1e-10
  )
})

test_that("a gam numerator model is estimated at its own penalized score", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()
  ps_mod <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, method = "REML")
  num_mod <- mgcv::gam(A ~ s(x1), data = dat, method = "REML")

  # A numerator model goes through the registry the propensity score model goes
  # through, so an additive fit is stacked there on the same terms rather than
  # refused in one argument and accepted in the other. Its coefficients are the
  # root of its own penalized score, and a block that dropped the penalty would
  # walk them off the fit and rebuild a numerator nobody estimated.
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(fitted(ps_mod)),
      dat$A,
      exposure_type = "continuous",
      stabilize = num_mod
    )
  )
  outcome_mod <- lm(yc ~ A, data = dat, weights = wts)

  layout <- expect_gam_root_seeded(ps_mod, outcome_mod, wts)

  res <- ipw(ps_mod, outcome_mod)
  theta <- coef(res$fit)
  numerator_block <- theta[layout$idx$stab][seq_along(coef(num_mod))]
  expect_equal(
    unname(numerator_block),
    unname(coef(num_mod)),
    tolerance = 1e-8
  )
  expect_equal(
    res$estimates$estimate,
    unname(coef(outcome_mod)[["A"]]),
    tolerance = 1e-8
  )
})

# ---- the least squares route, read as a gam ---------------------------------

test_that("a gam with no penalty reports what the lm route reports", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # Smoothing parameters of zero leave the penalty out of the score, so the
  # additive fit is ordinary least squares on its own basis and the stacked
  # system is the one the `lm` route already writes. That is the single case
  # where the two are the same equations, so it is where the new block can be
  # held to the old one's arithmetic rather than to a band.
  gam_mod <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, sp = c(0, 0))
  gam_fx <- fit_gam_models(dat, gam_mod)
  gam_res <- ipw(gam_fx$ps_mod, gam_fx$outcome_mod)

  basis <- stats::predict(gam_mod, type = "lpmatrix")[, -1, drop = FALSE]
  basis_dat <- as.data.frame(basis)
  names(basis_dat) <- paste0("b", seq_len(ncol(basis_dat)))
  basis_dat$A <- dat$A
  basis_dat$yc <- dat$yc
  lm_mod <- lm(
    stats::reformulate(paste0("b", seq_len(ncol(basis))), response = "A"),
    data = basis_dat
  )
  lm_fx <- fit_gam_models(basis_dat, lm_mod)
  lm_res <- ipw(lm_fx$ps_mod, lm_fx$outcome_mod)

  expect_equal(
    gam_res$estimates$estimate,
    lm_res$estimates$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    gam_res$estimates$std.err,
    lm_res$estimates$std.err,
    tolerance = 1e-5
  )
})

# ---- the classes the registry reads as additive fits ------------------------

test_that("a bam dose model is stacked as the gam it is", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # `bam()` is the same penalized fit computed for larger data, and it reports
  # its basis, its penalties, and its smoothing parameters exactly as `gam()`
  # does, so it is stacked through the same entry rather than refused for the
  # class it adds.
  ps_mod <- mgcv::bam(A ~ s(x1) + s(x2), data = dat, method = "fREML")
  mods <- fit_gam_models(dat, ps_mod)

  expect_gam_root_seeded(mods$ps_mod, mods$outcome_mod, mods$wts)

  res <- ipw(mods$ps_mod, mods$outcome_mod)
  theta <- coef(res$fit)
  expect_equal(
    unname(theta[seq_along(coef(ps_mod))]),
    unname(coef(ps_mod)),
    tolerance = 1e-8
  )
  expect_equal(
    res$estimates$estimate,
    unname(coef(mods$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )
})

test_that("the gam component of a gamm fit is refused for the class it is", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("nlme")
  dat <- sim_gam_dose()

  # `gamm()` returns its additive part as a bare `gam`, without the `glm` and
  # `lm` classes every route here reads a model through, so it is refused where
  # a class is first read rather than anywhere the additive machinery would
  # recognize it. The remedy is to refit with `gam()` or `bam()`.
  ps_mod <- mgcv::gamm(A ~ s(x1) + s(x2), data = dat)$gam
  expect_identical(class(ps_mod), "gam")
  mods <- fit_gam_models(dat, ps_mod)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_method_error"
  )
  expect_match(conditionMessage(err), "gam", fixed = TRUE)
})

# ---- limitation 1: the variance conditions on the smoothing parameters ------

test_that("a gam handed its own smoothing parameters reports the same fit", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  chosen <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, method = "REML")
  pinned <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, sp = chosen$sp)

  # The stacked system reads the smoothing parameters as constants, so a fit
  # given them and a fit that chose them are the same equations and report the
  # same standard error. That is the approximation stated as an equality: the
  # uncertainty of having chosen them is carried by neither.
  #
  # A fit given every smoothing parameter records them in `full.sp` and leaves
  # `sp` empty, so the two are read off the fit rather than off one field.
  expect_length(pinned$sp, 0L)

  chosen_fx <- fit_gam_models(dat, chosen)
  chosen_res <- ipw(chosen_fx$ps_mod, chosen_fx$outcome_mod)
  pinned_fx <- fit_gam_models(dat, pinned)
  pinned_res <- ipw(pinned_fx$ps_mod, pinned_fx$outcome_mod)

  expect_equal(
    pinned_res$estimates$estimate,
    chosen_res$estimates$estimate,
    tolerance = 1e-8
  )
  # The two fits are the same equations rather than the same arithmetic: their
  # coefficients agree to fourteen places and the sandwich's derivatives are
  # taken by finite differences, so the standard errors agree to six rather than
  # exactly. What the limitation would cost if the choice were carried is one to
  # three percent, four orders of magnitude above this.
  expect_equal(
    pinned_res$estimates$std.err,
    chosen_res$estimates$std.err,
    tolerance = 1e-5
  )
})

# ---- limitation 2: the frequentist variance rather than mgcv's Vp -----------

test_that("the propensity score block is the frequentist variance", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()
  ps_mod <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, method = "REML")
  mods <- fit_gam_models(dat, ps_mod)

  res <- ipw(mods$ps_mod, mods$outcome_mod)
  coefs <- seq_along(coef(ps_mod))
  block <- vcov(res$wt_mod)[coefs, coefs]

  # mgcv reports the Bayesian covariance of the smooth by default, which carries
  # the smoothing prior and is the wider of the two. What the sandwich
  # differentiates is the penalized score at fixed smoothing parameters, so the
  # block it reports is the frequentist reading: close to `freq = TRUE` and
  # nothing like the default.
  frequentist <- sqrt(diag(stats::vcov(ps_mod, freq = TRUE)))
  bayesian <- sqrt(diag(stats::vcov(ps_mod)))
  stacked <- sqrt(diag(block))

  expect_true(all(abs(stacked / frequentist - 1) < 0.15))
  expect_lt(min(stacked / bayesian), 0.7)
})

# ---- limitation 3: the envelope the glm route has ---------------------------

test_that("a log-link gam dose model is stacked at its own score", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # A gaussian additive fit through a log link reads its conditional mean as the
  # exponential of the smooth, and its penalized score carries the chain rule,
  # so the equation the coefficients solve is one this system can still write.
  ps_mod <- mgcv::gam(
    Apos ~ s(x1) + s(x2),
    data = dat,
    family = gaussian(link = "log"),
    method = "REML"
  )
  mods <- fit_gam_models(dat, ps_mod, exposure = "Apos")

  expect_gam_root_seeded(mods$ps_mod, mods$outcome_mod, mods$wts)

  res <- ipw(mods$ps_mod, mods$outcome_mod)
  expect_equal(
    unname(coef(res$fit)[seq_along(coef(ps_mod))]),
    unname(coef(ps_mod)),
    tolerance = 1e-8
  )
  expect_true(is.finite(res$estimates$std.err))
  expect_gt(res$estimates$std.err, 0)
})

test_that("the registry refuses an additive fit of another family", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()
  dat$counts <- rpois(nrow(dat), exp(0.3 + 0.2 * dat$x1))

  # The density the weights divide by has one spread, which is what a gaussian
  # family says and no other one does. A `gam` is read through the registry the
  # way every other class is, so the family it was fit with is refused there
  # with the class the rest of the package refuses a family with.
  ps_mod <- mgcv::gam(counts ~ s(x1), data = dat, family = poisson())

  expect_error(
    ipw_continuous_model(ps_mod),
    class = "propensity_model_family_error"
  )

  # The same fit handed to `ipw()` never reaches the additive stack either. A
  # model of another family is read as a model of a binary exposure there, one
  # step before the continuous registry, so what refuses it is that path's own
  # envelope rather than this one's.
  mods <- fit_gam_models(dat, ps_mod, exposure = "counts")

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_link_error"
  )
})

test_that("the registry refuses an additive fit through a link it cannot write", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # The links a gaussian fit can be stacked through are the registry's rather
  # than the additive machinery's: an identity link is least squares, a log link
  # is the score deli writes for a normal model of `exp(X alpha)`, and the
  # remaining gaussian links are refused by name because the coefficients an
  # IRLS iteration stops at under them are not a tight enough root to seed the
  # solve from. An additive fit is held to that envelope and not a wider one.
  ps_mod <- mgcv::gam(
    Apos ~ s(x1) + s(x2),
    data = dat,
    family = gaussian(link = "inverse"),
    method = "REML"
  )

  expect_error(
    ipw_continuous_model(ps_mod),
    class = "propensity_ipw_link_error"
  )

  mods <- fit_gam_models(dat, ps_mod, exposure = "Apos")

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_link_error"
  )
  expect_match(conditionMessage(err), "inverse", fixed = TRUE)
})

test_that("the registry refuses an additive fit whose penalty it cannot place", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # The penalty subtracted from the score is read off the fit's smooth terms,
  # one matrix at a time, at the coefficients each term owns. A penalty attached
  # to a parametric term belongs to no smooth and records a smoothing parameter
  # those terms do not account for, so the fit is refused rather than stacked at
  # a score missing a term of the equation its coefficients solve.
  ps_mod <- mgcv::gam(
    A ~ s(x1) + x2,
    data = dat,
    paraPen = list(x2 = list(diag(1))),
    method = "REML"
  )
  mods <- fit_gam_models(dat, ps_mod)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_se_method_unavailable_error"
  )
  expect_s3_class(err, "propensity_method_error")

  expect_propensity_error(ipw(mods$ps_mod, mods$outcome_mod))
})

test_that("the registry refuses a parametric penalty held at a fixed value", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # A fit handed some of its smoothing parameters records every penalty it holds
  # in `full.sp` and only the ones it chose in `sp`, so a parametric penalty
  # held fixed leaves `sp` exactly as long as the smooth terms are. Counting
  # against `sp` there would read the fit as one this path can place, and the
  # score it stacked would be missing the parametric penalty; the count is taken
  # against the full record whenever the fit keeps one.
  ps_mod <- mgcv::gam(
    A ~ s(x1) + x2,
    data = dat,
    paraPen = list(x2 = list(diag(1))),
    sp = c(2, -1),
    method = "REML"
  )
  expect_length(ps_mod$sp, 1L)
  expect_length(ps_mod$full.sp, 2L)

  mods <- fit_gam_models(dat, ps_mod)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_se_method_unavailable_error"
  )
  expect_s3_class(err, "propensity_method_error")
})

test_that("a gam handed some of its smoothing parameters is stacked", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # Every penalty in this fit belongs to a smooth term, so it is stacked even
  # though one smoothing parameter was handed to it and the other was chosen.
  # `full.sp` holds both, in the order the penalties are placed, and the penalty
  # built from it puts the fit's own coefficients at the root of the score.
  ps_mod <- mgcv::gam(
    A ~ s(x1) + s(x2),
    data = dat,
    sp = c(0.5, -1),
    method = "REML"
  )
  expect_length(ps_mod$sp, 1L)
  expect_length(ps_mod$full.sp, 2L)
  expect_equal(unname(ps_mod$full.sp[[1]]), 0.5)

  mods <- fit_gam_models(dat, ps_mod)

  expect_gam_root_seeded(mods$ps_mod, mods$outcome_mod, mods$wts)

  res <- ipw(mods$ps_mod, mods$outcome_mod)
  expect_equal(
    unname(coef(res$fit)[seq_along(coef(ps_mod))]),
    unname(coef(ps_mod)),
    tolerance = 1e-8
  )
})

test_that("the registry refuses a gam holding a smoothing floor", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # `min.sp` is a floor under the smoothing parameters, and the fit adds its
  # contribution to the penalty without recording it anywhere: `full.sp` is
  # documented as excluding that contribution, and the fitted object carries no
  # field holding it. The smoothing parameters read off such a fit therefore
  # rebuild a penalty the fit was never made under, its own coefficients are not
  # at the root of the score built from it, and a solve seeded there would walk
  # the propensity score block off the fit while reporting estimates that look
  # like any other. The count of the smoothing parameters says nothing about it,
  # so what catches it is the score itself.
  ps_mod <- mgcv::gam(
    A ~ s(x1) + s(x2),
    data = dat,
    min.sp = c(5, 5),
    method = "REML"
  )
  expect_equal(unname(ps_mod$full.sp), unname(ps_mod$sp))

  expect_false(ipw_continuous_model(ps_mod)$stackable)

  mods <- fit_gam_models(dat, ps_mod)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_se_method_unavailable_error"
  )
  expect_s3_class(err, "propensity_method_error")
  expect_match(conditionMessage(err), "min.sp", fixed = TRUE)

  expect_propensity_error(ipw(mods$ps_mod, mods$outcome_mod))
})

test_that("the registry refuses a smoothing floor too small to see", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # The unrecorded contribution is the floor itself, so a floor small enough to
  # leave the fit close to the one it would have been moves the score by little.
  # The refusal is on whether the record reproduces the score rather than on how
  # far off it is, so a floor four orders of magnitude below the smoothing
  # parameters the fit chose is refused the same way.
  ps_mod <- mgcv::gam(
    A ~ s(x1) + s(x2),
    data = dat,
    min.sp = c(1e-4, 1e-4),
    method = "REML"
  )
  mods <- fit_gam_models(dat, ps_mod)

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_se_method_unavailable_error"
  )
})

test_that("a bam handed a smoothing floor is stacked as the fit it made", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()

  # `bam()` does not implement a smoothing floor under its own criterion: it
  # says so and drops the argument, so what it returns is the fit its recorded
  # smoothing parameters describe. The check reads the score rather than the
  # arguments the call was written with, so this fit stacks.
  expect_warning(
    ps_mod <- mgcv::bam(
      A ~ s(x1) + s(x2),
      data = dat,
      min.sp = c(5, 5),
      method = "fREML"
    ),
    "min.sp"
  )

  mods <- fit_gam_models(dat, ps_mod)

  expect_gam_root_seeded(mods$ps_mod, mods$outcome_mod, mods$wts)

  res <- ipw(mods$ps_mod, mods$outcome_mod)
  expect_equal(
    unname(coef(res$fit)[seq_along(coef(ps_mod))]),
    unname(coef(ps_mod)),
    tolerance = 1e-8
  )
})

test_that("ipw() refuses a gam dose model fit with case weights", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()
  dat$w <- rep(c(1, 2), length.out = nrow(dat))

  # The stacked penalized score is unweighted, so a fit made under prior weights
  # is not at its root and the system would walk away from the coefficients the
  # user has. The refusal is the one the rest of the continuous path makes.
  ps_mod <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, weights = w)
  mods <- fit_gam_models(dat, ps_mod)

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_ps_weights_error"
  )
})

# ---- the weights side -------------------------------------------------------

test_that("wt_ate() reads a gam fit and records the ratio it built", {
  skip_if_not_installed("mgcv")
  dat <- sim_gam_dose()
  ps_mod <- mgcv::gam(A ~ s(x1) + s(x2), data = dat, method = "REML")

  # Pins behavior the package already has, which the stacked system now depends
  # on: the weights a `gam` gives are the ones its conditional means give, and
  # the record they carry is what the sandwich rebuilds them from.
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_mod, dat$A, exposure_type = "continuous", stabilize = TRUE)
  )

  expect_equal(
    as.double(wts),
    gam_hand_weights(as.double(fitted(ps_mod)), dat$A),
    tolerance = 1e-12
  )
  expect_identical(exposure_type(wts), "continuous")
  expect_true(is_stabilized(wts))
  expect_identical(density_meta(wts)$sigma, "pooled")
  expect_identical(density_meta(wts)$numerator, "marginal")
  expect_identical(format(density_meta(wts)$density), "normal")
})
