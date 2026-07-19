# Tests for ipw() with a continuous exposure: an lm (or gaussian-family glm)
# propensity score model of the exposure on covariates, stabilized continuous
# ATE weights, and a weighted marginal structural outcome model with exactly one
# exposure term. These exercise ipw() end to end the way a user would call it,
# pinning estimand detection, the estimates-table contract, point-estimate parity
# with the MSM coefficient, a bootstrap standard-error oracle, the stabilization
# variants, the estimand and MSM guards, the gaussian-glm routing, the fit
# object, and the print and as.data.frame surfaces.

# ---- data simulator ---------------------------------------------------------

sim_continuous <- function(seed = 2024, n = 800) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  yc <- 1 + 0.6 * A + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  yb <- rbinom(n, 1, plogis(-0.5 + 0.4 * A + 0.3 * x1 - 0.2 * x2))
  data.frame(x1, x2, A, yc, yb)
}

# ---- model fitting ----------------------------------------------------------

# Continuous ATE weights from a numeric fitted propensity, always computed
# silently. Stabilized weights are the recommended default for a continuous
# exposure; stabilize = FALSE emits an alert unless quieted.
continuous_weights <- function(
  fitted_ps,
  A,
  stabilize = TRUE,
  stab_score = NULL
) {
  withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      fitted_ps,
      A,
      exposure_type = "continuous",
      stabilize = stabilize,
      stabilization_score = stab_score
    )
  )
}

# Fit the propensity score model of A on the covariates, build continuous ATE
# weights from its fitted values, and fit a weighted MSM. `ps_type` selects the
# lm or the gaussian-family glm form of the propensity model; the two share
# fitted values and so produce identical weights. `msm_rhs` allows a
# multiple-term right-hand side for the MSM guard. The weights are kept as a psw
# object so the estimand survives into the outcome model frame for detection. The
# quasibinomial MSM tightens its IRLS tolerance so its coefficients sit at the
# weighted MLE to well below the point-estimate comparison tolerance.
fit_continuous_models <- function(
  dat,
  ps_type = c("lm", "glm"),
  outcome_family = c("gaussian", "binomial"),
  stabilize = TRUE,
  stab_score = NULL,
  msm_rhs = "A"
) {
  ps_type <- match.arg(ps_type)
  outcome_family <- match.arg(outcome_family)

  ps_mod <- if (ps_type == "lm") {
    lm(A ~ x1 + x2, data = dat)
  } else {
    glm(A ~ x1 + x2, data = dat, family = gaussian())
  }
  fitted_ps <- as.double(fitted(ps_mod))
  wts <- continuous_weights(
    fitted_ps,
    dat$A,
    stabilize = stabilize,
    stab_score = stab_score
  )

  outcome_var <- if (outcome_family == "binomial") "yb" else "yc"
  msm_fmla <- stats::reformulate(msm_rhs, response = outcome_var)
  if (outcome_family == "binomial") {
    outcome_mod <- glm(
      msm_fmla,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  } else {
    outcome_mod <- lm(msm_fmla, data = dat, weights = wts)
  }

  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts)
}

estimates_columns <- c(
  "effect",
  "estimate",
  "std.err",
  "z",
  "ci.lower",
  "ci.upper",
  "conf.level",
  "p.value"
)

# ---- end-to-end -------------------------------------------------------------

test_that("ipw() runs continuous ate end to end and auto-detects the estimand", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_s3_class(res, "ipw")
  expect_equal(res$se_method, "mestimation")
  expect_equal(res$estimand, "ate")
  expect_false(is.null(res$fit))
})

# ---- point-estimate parity with the MSM exposure coefficient ----------------

test_that("ipw() continuous point estimate equals the MSM exposure coefficient", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()

  # identity link (lm) is exact
  mods_g <- fit_continuous_models(dat, outcome_family = "gaussian")
  res_g <- ipw(mods_g$ps_mod, mods_g$outcome_mod)
  expect_equal(nrow(res_g$estimates), 1L)
  expect_equal(
    res_g$estimates$estimate,
    unname(coef(mods_g$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )

  # logit link (quasibinomial) tightened to its weighted MLE
  mods_b <- fit_continuous_models(dat, outcome_family = "binomial")
  res_b <- ipw(mods_b$ps_mod, mods_b$outcome_mod)
  expect_equal(nrow(res_b$estimates), 1L)
  expect_equal(
    res_b$estimates$estimate,
    unname(coef(mods_b$outcome_mod)[["A"]]),
    tolerance = 1e-8
  )
})

# ---- estimates table shape and effect labels --------------------------------

test_that("ipw() continuous labels the effect by the outcome link", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()

  # identity link: a single slope, the eight-column contract, no comparison
  mods_g <- fit_continuous_models(dat, outcome_family = "gaussian")
  res_g <- ipw(mods_g$ps_mod, mods_g$outcome_mod)
  expect_named(res_g$estimates, estimates_columns)
  expect_false("comparison" %in% names(res_g$estimates))
  expect_equal(nrow(res_g$estimates), 1L)
  expect_equal(res_g$estimates$effect, "slope")
  expect_true(all(res_g$estimates$conf.level == 0.95))

  # logit link: a single log odds ratio
  mods_b <- fit_continuous_models(dat, outcome_family = "binomial")
  res_b <- ipw(mods_b$ps_mod, mods_b$outcome_mod)
  expect_named(res_b$estimates, estimates_columns)
  expect_false("comparison" %in% names(res_b$estimates))
  expect_equal(nrow(res_b$estimates), 1L)
  expect_equal(res_b$estimates$effect, "log(or)")
})

# ---- standard-error oracle: nonparametric bootstrap -------------------------

# The M-estimation sandwich standard error, which accounts for estimating the
# propensity model, is validated against a seed-pinned nonparametric bootstrap
# that refits the propensity model, the weights, and the MSM on each resample.
# The two are asymptotically equivalent; a 15 percent relative band accommodates
# the finite-sample and Monte Carlo gap while keeping this a meaningful check.
test_that("ipw() continuous mestimation SE tracks a nonparametric bootstrap", {
  skip_if_not_installed("deli")
  skip_on_cran()

  dat <- sim_continuous(seed = 2024, n = 600)
  mods <- fit_continuous_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod)
  mest_se <- res$estimates$std.err

  boot_slope <- function(d) {
    ps <- lm(A ~ x1 + x2, data = d)
    fitted_ps <- as.double(fitted(ps))
    w <- continuous_weights(fitted_ps, d$A, stabilize = TRUE)
    msm <- lm(yc ~ A, data = d, weights = as.double(w))
    unname(coef(msm)[["A"]])
  }

  withr::local_seed(918)
  reps <- 400L
  n <- nrow(dat)
  boot <- vapply(
    seq_len(reps),
    function(i) {
      boot_slope(dat[sample.int(n, n, replace = TRUE), , drop = FALSE])
    },
    numeric(1)
  )
  boot_se <- stats::sd(boot)

  rel <- abs(mest_se - boot_se) / boot_se
  expect_lt(rel, 0.15)
})

# ---- stabilization variants -------------------------------------------------

test_that("ipw() continuous runs stabilized, unstabilized, and fixed-score weights", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()

  # stabilized: the recommended default for a continuous exposure
  mods_s <- fit_continuous_models(dat, stabilize = TRUE)
  res_s <- ipw(mods_s$ps_mod, mods_s$outcome_mod)
  expect_equal(res_s$estimand, "ate")
  expect_true(is.finite(res_s$estimates$estimate))
  expect_true(
    is.finite(res_s$estimates$std.err) && res_s$estimates$std.err > 0
  )

  # unstabilized: the ipw() call is silent and the SE is finite (typically
  # larger than the stabilized SE)
  mods_u <- fit_continuous_models(dat, stabilize = FALSE)
  res_u <- expect_silent(ipw(mods_u$ps_mod, mods_u$outcome_mod))
  expect_true(is.finite(res_u$estimates$estimate))
  expect_true(
    is.finite(res_u$estimates$std.err) && res_u$estimates$std.err > 0
  )

  # user-supplied stabilization score
  mods_f <- fit_continuous_models(dat, stabilize = TRUE, stab_score = 0.5)
  res_f <- ipw(mods_f$ps_mod, mods_f$outcome_mod)
  expect_true(is.finite(res_f$estimates$estimate))
  expect_true(
    is.finite(res_f$estimates$std.err) && res_f$estimates$std.err > 0
  )
})

# ---- estimand support -------------------------------------------------------

test_that("ipw() continuous rejects estimands other than ate", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  # a valid non-ate estimand is unsupported for a continuous exposure
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, estimand = "att"),
    class = "propensity_ipw_estimand_error"
  )

  # an entirely invalid estimand string errors with the same class
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, estimand = "nonsense"),
    class = "propensity_ipw_estimand_error"
  )
})

# ---- MSM guard --------------------------------------------------------------

test_that("ipw() continuous rejects an MSM with more than one exposure term", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, msm_rhs = c("A", "I(A^2)"))

  # more than one exposure term has no single reported effect; the error must
  # direct the user to the returned fit object for the full coefficient vector
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod),
    class = "propensity_ipw_msm_error",
    regexp = "fit"
  )
})

# ---- gaussian-glm routing ---------------------------------------------------

test_that("ipw() routes a gaussian-family glm ps model identically to lm", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()

  mods_lm <- fit_continuous_models(dat, ps_type = "lm")
  mods_glm <- fit_continuous_models(dat, ps_type = "glm")

  res_lm <- ipw(mods_lm$ps_mod, mods_lm$outcome_mod)
  res_glm <- ipw(mods_glm$ps_mod, mods_glm$outcome_mod)

  expect_equal(res_glm$estimand, res_lm$estimand)
  expect_equal(
    res_glm$estimates$estimate,
    res_lm$estimates$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    res_glm$estimates$std.err,
    res_lm$estimates$std.err,
    tolerance = 1e-8
  )
})

# ---- fit object -------------------------------------------------------------

test_that("ipw() continuous fit theta matches the reported estimate", {
  skip_if_not_installed("deli")
  skip_if_not_installed("generics")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  # the MSM exposure coefficient carries the outcome-model term name, unique
  # across theta since the propensity block uses the covariate names
  tidied <- generics::tidy(res$fit)
  a_rows <- tidied[tidied$term == "A", ]
  expect_equal(nrow(a_rows), 1L)
  expect_equal(a_rows$estimate, res$estimates$estimate, tolerance = 1e-8)
})

# ---- print and as.data.frame ------------------------------------------------

test_that("ipw() continuous print output is stable", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, outcome_family = "gaussian")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  expect_snapshot(print(res))
})

test_that("as.data.frame(exponentiate = TRUE) relabels the continuous log odds ratio", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, outcome_family = "binomial")
  res <- ipw(mods$ps_mod, mods$outcome_mod)

  df <- as.data.frame(res, exponentiate = TRUE)
  expect_equal(nrow(df), 1L)
  expect_equal(df$effect, "or")
})

# ---- offset guard -----------------------------------------------------------

test_that("ipw() continuous rejects an outcome model with an offset term", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  dat$off <- 0.3 * dat$x1

  ps_mod <- lm(A ~ x1 + x2, data = dat)
  fitted_ps <- as.double(fitted(ps_mod))
  wts <- continuous_weights(fitted_ps, dat$A, stabilize = TRUE)
  outcome_mod <- lm(yc ~ A + offset(off), data = dat, weights = wts)

  # The mestimation psi does not thread the offset into the MSM score block, so
  # the default path must reject an offset before any estimation and with no
  # message or other condition signaled first.
  expect_no_message(
    expect_error(
      ipw(ps_mod, outcome_mod),
      class = "propensity_ipw_offset_error",
      regexp = "offset"
    )
  )
})

# ---- linearization is unsupported for continuous ----------------------------

test_that("ipw() rejects linearization for a continuous exposure", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()

  # the lm propensity route
  mods_lm <- fit_continuous_models(dat, ps_type = "lm")
  expect_error(
    ipw(mods_lm$ps_mod, mods_lm$outcome_mod, se_method = "linearization"),
    class = "propensity_method_error"
  )

  # the gaussian-family glm propensity route rejects it identically
  mods_glm <- fit_continuous_models(dat, ps_type = "glm")
  expect_error(
    ipw(mods_glm$ps_mod, mods_glm$outcome_mod, se_method = "linearization"),
    class = "propensity_method_error"
  )
})

# ---- outcome-family validation ----------------------------------------------

test_that("ipw_spec_continuous rejects an unsupported outcome family", {
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)

  withr::local_seed(1)
  dat$ycount <- rpois(nrow(dat), exp(0.1 + 0.2 * dat$A))
  outcome_mod <- suppressWarnings(
    glm(ycount ~ A, data = dat, family = poisson(), weights = wts)
  )

  # A poisson marginal structural model is silently stacked as a binomial score
  # today.
  expect_error(
    ipw_spec_continuous(ps_mod, outcome_mod),
    class = "propensity_ipw_family_error"
  )
})

# ---- continuous-path link validation ----------------------------------------
#
# Two continuous-path inputs fail late or dishonestly today. A gaussian
# propensity model with a non-identity link reconstructs its linear predictor as
# the fitted mean, so the weights the user built from fitted() no longer match
# and the misleading weights-mismatch error blames the weights instead of the
# unsupported link. A probit marginal structural model is accepted by the family
# check but has no continuous effect label, so it errors late with a terse
# internal message rather than fast at entry listing the supported links. These
# tests pin an honest, fast, classed error; they cannot pass until the
# continuous-path link validation exists. (A gaussian-identity glm ps model,
# already covered by the gaussian-glm routing test above, and a logit msm,
# covered by the effect-label test above, must keep working; the log-link msm
# below adds the one supported msm link not otherwise exercised.)

test_that("ipw() rejects a non-identity link on the continuous propensity model", {
  dat <- sim_continuous()
  withr::local_seed(7)
  dat$Apos <- exp(0.3 + 0.4 * dat$x1 - 0.2 * dat$x2 + 0.3 * rnorm(nrow(dat)))
  dat$ypos <- 1 + 0.6 * dat$Apos + 0.5 * dat$x1 + rnorm(nrow(dat))
  ps_mod <- suppressWarnings(
    glm(Apos ~ x1 + x2, data = dat, family = gaussian(link = "log"))
  )
  wts <- continuous_weights(fitted(ps_mod), dat$Apos)
  msm <- lm(ypos ~ Apos, data = dat, weights = wts)

  # class propensity_ipw_link_error is the natural fit; the implementer may
  # adjust it, but the error must name the unsupported propensity model link
  # rather than blame the outcome model weights.
  expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_link_error"
  )
})

test_that("the continuous propensity-link error names the unsupported link", {
  dat <- sim_continuous()
  withr::local_seed(7)
  dat$Apos <- exp(0.3 + 0.4 * dat$x1 - 0.2 * dat$x2 + 0.3 * rnorm(nrow(dat)))
  dat$ypos <- 1 + 0.6 * dat$Apos + 0.5 * dat$x1 + rnorm(nrow(dat))
  ps_mod <- suppressWarnings(
    glm(Apos ~ x1 + x2, data = dat, family = gaussian(link = "log"))
  )
  wts <- continuous_weights(fitted(ps_mod), dat$Apos)
  msm <- lm(ypos ~ Apos, data = dat, weights = wts)

  expect_snapshot(error = TRUE, ipw(ps_mod, msm))
})

test_that("ipw() rejects a non-identity link on the continuous marginal structural model", {
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  msm <- suppressWarnings(
    glm(yb ~ A, data = dat, family = binomial(link = "probit"), weights = wts)
  )

  # The continuous effect label supports only identity, logit, and log. A probit
  # msm must fail fast at entry with a message listing the supported links, not
  # the terse late "Unsupported outcome link" message it produces today.
  expect_error(
    ipw(ps_mod, msm),
    class = "propensity_ipw_link_error",
    regexp = "identity"
  )
})

test_that("ipw() continuous works with a log-link marginal structural model", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  # Low baseline risk and a weak dose effect keep the log-binomial fitted
  # probabilities below 1 so the marginal structural model converges.
  withr::local_seed(11)
  dat$ylow <- rbinom(
    nrow(dat),
    1,
    exp(-1.6 + 0.15 * dat$A + 0.1 * dat$x1)
  )
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  msm <- suppressWarnings(glm(
    ylow ~ A,
    data = dat,
    family = binomial(link = "log"),
    weights = wts,
    start = c(-1.6, 0.15),
    control = glm.control(epsilon = 1e-14, maxit = 200)
  ))

  res <- ipw(ps_mod, msm)
  expect_equal(res$estimates$effect, "log(rr)")
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$std.err)))
})

# ---- factor outcome response ------------------------------------------------

test_that("continuous logistic MSM matches a factor outcome response to the numeric fit", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  dat$ybf <- factor(ifelse(dat$yb == 1, "yes", "no"), levels = c("no", "yes"))
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  ctrl <- glm.control(epsilon = 1e-14, maxit = 200)

  num <- glm(
    yb ~ A,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = ctrl
  )
  fac <- suppressWarnings(
    glm(
      ybf ~ A,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = ctrl
    )
  )
  res_num <- ipw(ps_mod, num)$estimates
  res_fac <- ipw(ps_mod, fac)$estimates
  expect_equal(res_fac$estimate, res_num$estimate, tolerance = 1e-8)
  expect_equal(res_fac$std.err, res_num$std.err, tolerance = 1e-8)
})
