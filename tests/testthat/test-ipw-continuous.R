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

  # Without the family guard a poisson marginal structural model is stacked as a
  # binomial score.
  expect_error(
    ipw_spec_continuous(ps_mod, outcome_mod),
    class = "propensity_ipw_family_error"
  )
})

# ---- continuous-path link validation ----------------------------------------
#
# Two continuous-path inputs are rejected at entry with a classed error naming
# the supported links. Left to their downstream failures both mislead: a
# gaussian propensity model with a non-identity link reconstructs its linear
# predictor as the fitted mean, so the weights the user built from fitted() no
# longer match and the weights-mismatch error blames the weights instead of the
# unsupported link; a probit marginal structural model passes the family check
# but has no continuous effect label, so it errors late with a terse internal
# message. (A gaussian-identity glm ps model, already covered by the
# gaussian-glm routing test above, and a logit msm, covered by the effect-label
# test above, must keep working; the log-link msm below adds the one supported
# msm link not otherwise exercised.)

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
  # msm fails at entry with a message listing the supported links, rather than
  # reaching the terse "Unsupported outcome link" message downstream.
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

# ---- ps_link on the continuous path -----------------------------------------
#
# ps_link is meaningful only for a binomial glm on the binary path, where it
# overrides the model's link. A continuous propensity model has no such link, so
# both continuous entries, the lm method and the gaussian-family branch of
# ipw.glm, must reject it rather than silently ignore it. The default
# ps_link = NULL is covered by every other test in this file, through both
# entries (see the end-to-end and gaussian-glm routing tests above).

test_that("ipw() rejects ps_link for an lm propensity model", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit"),
    class = "propensity_ipw_link_error",
    regexp = "ps_link"
  )
})

test_that("ipw() rejects ps_link for a gaussian glm propensity model", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, ps_type = "glm")

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit"),
    class = "propensity_ipw_link_error",
    regexp = "ps_link"
  )
})

test_that("the gaussian glm ps_link rejection deprecates nothing", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat, ps_type = "glm")

  # `ps_link` is deprecated on the binary path, and the gaussian branch of
  # ipw.glm returns before that warning is reached. A model with no link to
  # override must be rejected outright rather than told to drop the argument
  # and to correct it in the same breath.
  #
  # The signal is forced unconditionally, so the empty record below is evidence
  # that this path raises nothing rather than that something else reached the
  # `ipw(ps_link)` id first.
  deprecations <- character()
  err <- with_always_deprecated(
    withCallingHandlers(
      expect_error(
        ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit"),
        class = "propensity_ipw_link_error"
      ),
      lifecycle_warning_deprecated = function(cnd) {
        deprecations <<- c(deprecations, conditionMessage(cnd))
        rlang::cnd_muffle(cnd)
      }
    )
  )

  expect_length(deprecations, 0)

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "continuous propensity score model", fixed = TRUE)
})

test_that("the continuous ps_link error explains why the argument does not apply", {
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, mods$outcome_mod, ps_link = "logit")
  )
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

# ---- ps design reconstruction from .data ------------------------------------

test_that("continuous reconstructs the ps design from .data when the fit frame is gone", {
  skip_if_not_installed("deli")
  # The continuous path routes the ps design through ipw_extract_ps_design:
  # with the fitting data gone, no .data raises the guarded error and .data
  # rebuilds the design. The binary path matches this behavior; see
  # test-ipw-mestimation.R.
  dat <- sim_continuous()
  ps_gone <- lm(A ~ x1 + x2, data = dat, model = FALSE)
  wts <- continuous_weights(fitted(ps_gone), dat$A)
  om <- lm(yc ~ A, data = dat, weights = wts)
  ps_ref <- lm(A ~ x1 + x2, data = dat)
  dat_copy <- dat
  rm(dat)

  expect_error(ipw(ps_gone, om), class = "propensity_ipw_data_error")

  ref <- ipw(ps_ref, om)$estimates
  res <- ipw(ps_gone, om, .data = dat_copy)$estimates
  expect_equal(res, ref, tolerance = 1e-8)
})

test_that("continuous errors informatively when the outcome model frame is gone", {
  skip_if_not_installed("deli")
  # The companion to the binary and categorical pins in test-ipw-mestimation.R.
  # The guard lives in the weight extraction every ipw() method runs first, so
  # the continuous route gets it by sharing that helper rather than by having
  # its own check. This pin exists so a refactor that moves this path off the
  # shared helper is caught rather than silently losing the guided error.
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  outcome_gone <- local({
    d_local <- dat
    w_local <- wts
    m <- lm(yc ~ A, data = d_local, weights = w_local, model = FALSE)
    rm(d_local, w_local)
    m
  })

  # the fixture is genuinely unreconstructable
  expect_error(stats::model.frame(outcome_gone))

  for (supplied in list(NULL, dat)) {
    err <- expect_error(
      ipw(ps_mod, outcome_gone, .data = supplied),
      class = "propensity_ipw_data_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "outcome_mod", fixed = TRUE)
  }
})

test_that("continuous rejects a .data whose covariate types skew the design", {
  skip_if_not_installed("deli")
  # The continuous companion to the binary and categorical width pins. Same
  # rationale as the frame-gone pin above: the check lives in the shared
  # extraction, and this fixes that it covers this route.
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  dat_skew <- dat
  dat_skew$x1 <- cut(dat$x1, breaks = 3, labels = c("lo", "mid", "hi"))
  expect_equal(length(coef(mods$ps_mod)), 3)
  expect_equal(ncol(model.matrix(~ x1 + x2, data = dat_skew)), 4)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat_skew),
    class = "propensity_ipw_data_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, ".data", fixed = TRUE)
})

test_that("a log-transformed continuous outcome response through .data matches the model-frame route", {
  skip_if_not_installed("deli")
  # The continuous companion to the binary pins in test-ipw-se-method.R. Same
  # rationale as the pins above: the response is read out of `.data` by the
  # shared extraction, so a transformation of a column has to be computed rather
  # than looked up under the name of the function wrapping it. The exponentiated
  # outcome keeps `log()` defined on every observation.
  dat <- sim_continuous()
  dat$yp <- exp(dat$yc)
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(fitted(ps_mod), dat$A)
  outcome_mod <- lm(log(yp) ~ A, data = dat, weights = wts)

  expect_equal(
    ipw(ps_mod, outcome_mod, .data = dat)$estimates,
    ipw(ps_mod, outcome_mod)$estimates,
    tolerance = 1e-10
  )
})

# ---- a wrong-sized .data is a data problem -----------------------------------
#
# As on the binary and categorical paths, a `.data` whose row count disagrees
# with the fitted models is a data problem and reads as one only if it is
# reported as one. Without the guard it would reach the weight-consistency
# preflight, which talks about weights instead. The binary companion is "mestimation rejects a
# .data with too few rows" in test-ipw-mestimation.R.

test_that("ipw() continuous rejects a .data with too few rows", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  err <- expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, .data = dat[-1, ]),
    class = "propensity_ipw_data_error"
  )

  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "one row per observation", fixed = TRUE)
})

# ---- weight-consistency hint wording ----------------------------------------

test_that("the continuous weight-consistency error names no focal level", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)
  spec <- ipw_spec_continuous(mods$ps_mod, mods$outcome_mod)
  doubled <- 2 * as.double(model.weights(model.frame(mods$outcome_mod)))

  err <- expect_error(
    ipw_check_weight_consistency(spec, doubled),
    class = "propensity_ipw_weights_mismatch_error"
  )

  # A continuous exposure has no focal level, so the focal-level hint the binary
  # and categorical paths offer must be absent rather than reworded. The binary
  # wording is pinned by "the weight-consistency error names a focal level as a
  # possible cause" in test-ipw-mestimation.R.
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_false(grepl("focal", msg, fixed = TRUE))
})

# ---- propensity model class validation --------------------------------------
#
# The continuous path stacks an ordinary least-squares score block for the
# propensity model, so the M-estimator solves its coefficients to the
# least-squares root no matter how the supplied model was actually fit. An lm
# subclass whose coefficients are not that root would therefore yield estimates
# for a propensity model the user never fit, and no downstream guard catches it:
# MASS::rlm carries class c("rlm", "lm") and so reaches the lm method, and the
# weights built from its fitted values agree with the weights recomputed at the
# seeded init, so the weight-consistency preflight passes. The same holds
# through the gaussian-family branch of the glm method, which routes a gaussian
# mgcv::gam to the identical path. The continuous path therefore accepts only a
# plain lm or a gaussian glm and rejects any other subclass at entry, naming the
# class it was given.

# A fixture whose robust and least-squares propensity fits genuinely disagree:
# adding a block of large outliers to the exposure pulls the least-squares fit
# away from the robust fit, so the two imply different weights and different
# estimates. The outcome is regenerated from the contaminated exposure, under
# its own seed, so the marginal structural model remains coherent.
sim_continuous_outliers <- function(seed = 2024, n = 800, n_outliers = 40) {
  dat <- sim_continuous(seed = seed, n = n)
  dat$A[seq_len(n_outliers)] <- dat$A[seq_len(n_outliers)] + 12
  withr::local_seed(seed + 1L)
  dat$yc <- 1 + 0.6 * dat$A + 0.5 * dat$x1 - 0.3 * dat$x2 + rnorm(n)
  # yb came from the pre-contamination exposure and no caller reads it.
  dat$yb <- NULL
  dat
}

test_that("ipw() rejects a robust linear propensity model on the continuous path", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("deli")
  dat <- sim_continuous_outliers()

  ps_mod <- MASS::rlm(A ~ x1 + x2, data = dat)
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  # The test pins the general propensity_error class rather than a specific
  # subclass; what matters is that the error names the class it was handed
  # instead of accepting the model and reporting the least-squares analysis.
  expect_error(
    ipw(ps_mod, msm),
    class = "propensity_error",
    regexp = "rlm"
  )

  # Control: on the same fixture a plain lm propensity model still runs. Parity
  # between the lm and gaussian-glm routes is covered above by "ipw() routes a
  # gaussian-family glm ps model identically to lm".
  lm_mod <- lm(A ~ x1 + x2, data = dat)
  lm_wts <- continuous_weights(as.double(fitted(lm_mod)), dat$A)
  lm_msm <- lm(yc ~ A, data = dat, weights = lm_wts)
  expect_s3_class(ipw(lm_mod, lm_msm), "ipw")
})

test_that("ipw() rejects a gaussian gam propensity model on the continuous path", {
  skip_if_not_installed("mgcv")
  dat <- sim_continuous_outliers()

  ps_mod <- mgcv::gam(A ~ s(x1) + x2, data = dat, family = gaussian())
  wts <- continuous_weights(as.double(fitted(ps_mod)), dat$A)
  msm <- lm(yc ~ A, data = dat, weights = wts)

  # A gaussian gam reaches the continuous path through the glm method, so the
  # restriction has to hold at that entry point too and not only in the lm
  # method.
  expect_error(
    ipw(ps_mod, msm),
    class = "propensity_error",
    regexp = "gam"
  )
})

# ---- arguments that fall into the dots --------------------------------------

test_that("ipw() lm rejects arguments that fall into the dots", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  # `.data` was the third positional argument in earlier releases; it now lands
  # in the dots, where it would be discarded without a signal.
  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, dat),
    class = "rlib_error_dots_nonempty"
  )

  expect_error(
    ipw(mods$ps_mod, mods$outcome_mod, nonsense = 42),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("ipw() lm accepts every argument supplied by name", {
  skip_if_not_installed("deli")
  dat <- sim_continuous()
  mods <- fit_continuous_models(dat)

  baseline <- ipw(mods$ps_mod, mods$outcome_mod)
  named <- ipw(
    mods$ps_mod,
    mods$outcome_mod,
    .data = dat,
    estimand = "ate",
    conf_level = 0.95,
    se_method = "mestimation"
  )

  expect_equal(named$estimates, baseline$estimates)
})

# ---- the propensity model's response has to be a column ---------------------
#
# The continuous path derives the exposure's name by deparsing the propensity
# model's left-hand side, so a left-hand side that is not a bare symbol gives
# back a vector of several names and everything downstream is indexed by it. A
# matrix response and a transformed response both do this, and neither was
# guarded here: without `.data` the call died on an internal length assertion
# about `.exposure_name`, and with `.data` it asked for a column named after the
# transforming function, which cannot exist. The binary path already rejects
# both shapes, and the contract is the same one for a continuous exposure: the
# propensity model's response is the exposure column itself.

continuous_lhs_data <- function() {
  dat <- sim_continuous()
  # a strictly positive copy of the exposure, so log() of it is defined, and a
  # precomputed column holding that transformation
  dat$Apos <- exp(dat$A)
  dat$A1 <- round(pmax(dat$A, 0.1) * 10)
  dat$A2 <- round(pmax(-dat$A, 0.1) * 10)
  dat
}

# Continuous ATE weights and an MSM to carry them, built from whatever exposure
# vector the propensity model models. The models are only ever rejected here, so
# the weights need to be well formed rather than matched to the fit.
continuous_lhs_outcome <- function(dat, ps_mod, exposure) {
  wts <- continuous_weights(
    as.double(fitted(ps_mod))[seq_along(exposure)],
    exposure
  )
  lm(yc ~ A, data = dat, weights = wts)
}

test_that("a matrix-response continuous propensity model is rejected on its shape", {
  dat <- continuous_lhs_data()
  ps_mod <- lm(cbind(A1, A2) ~ x1 + x2, data = dat)
  ps_ok <- lm(A ~ x1 + x2, data = dat)
  outcome_mod <- continuous_lhs_outcome(dat, ps_ok, dat$A)

  # the fixture is the shape it claims to be
  expect_true(is.matrix(stats::model.response(stats::model.frame(ps_mod))))

  err <- expect_error(
    ipw(ps_mod, outcome_mod, .data = dat),
    class = "propensity_ipw_response_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "matrix response", fixed = TRUE)

  # and the spec constructor, which callers can reach directly, rejects it the
  # same way rather than on the internal length assertion it used to
  err_spec <- expect_error(
    ipw_spec_continuous(ps_mod, outcome_mod, .data = dat),
    class = "propensity_ipw_response_error"
  )
  expect_false(grepl(
    ".exposure_name",
    conditionMessage(err_spec),
    fixed = TRUE
  ))
})

test_that("a transformed continuous propensity response is rejected on both routes", {
  dat <- continuous_lhs_data()
  ps_mod <- lm(log(Apos) ~ x1 + x2, data = dat)
  outcome_mod <- continuous_lhs_outcome(dat, ps_mod, log(dat$Apos))

  for (data_arg in list(NULL, dat)) {
    err <- expect_error(
      ipw(ps_mod, outcome_mod, .data = data_arg),
      class = "propensity_ipw_response_error"
    )
    msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
    expect_match(msg, "log(Apos)", fixed = TRUE)
    # the two errors it used to raise, one per route
    expect_false(grepl(".exposure_name", msg, fixed = TRUE))
    expect_false(grepl("\"log\" column", msg, fixed = TRUE))
  }
})

test_that("a transformed continuous propensity response is rejected in the spec", {
  dat <- continuous_lhs_data()
  ps_mod <- lm(log(Apos) ~ x1 + x2, data = dat)
  outcome_mod <- continuous_lhs_outcome(dat, ps_mod, log(dat$Apos))

  expect_error(
    ipw_spec_continuous(ps_mod, outcome_mod, .data = dat),
    class = "propensity_ipw_response_error"
  )
})

test_that("the lm and gaussian glm propensity routes reject a transformed response alike", {
  dat <- continuous_lhs_data()
  ps_lm <- lm(log(Apos) ~ x1 + x2, data = dat)
  ps_glm <- glm(log(Apos) ~ x1 + x2, data = dat, family = gaussian())
  outcome_mod <- continuous_lhs_outcome(dat, ps_lm, log(dat$Apos))

  # The two propensity models share fitted values, so they have to reject the
  # same shape with the same message. The gaussian glm reached the binary
  # guard's cbind wording, which named a matrix this model does not have.
  messages <- vapply(
    list(ps_lm, ps_glm),
    function(ps_mod) {
      err <- expect_error(
        ipw(ps_mod, outcome_mod, .data = dat),
        class = "propensity_ipw_response_error"
      )
      gsub("[[:space:]]+", " ", conditionMessage(err))
    },
    character(1)
  )

  expect_equal(messages[[1]], messages[[2]])
  expect_false(grepl("cbind", messages[[1]], fixed = TRUE))
})

test_that("the transformed continuous propensity response error reads in the user's terms", {
  dat <- continuous_lhs_data()
  ps_mod <- lm(log(Apos) ~ x1 + x2, data = dat)
  outcome_mod <- continuous_lhs_outcome(dat, ps_mod, log(dat$Apos))

  expect_snapshot(
    error = TRUE,
    ipw(ps_mod, outcome_mod, .data = dat)
  )
})

test_that("a plain continuous propensity response is still accepted", {
  skip_if_not_installed("deli")
  dat <- continuous_lhs_data()
  mods <- fit_continuous_models(dat)

  expect_s3_class(ipw(mods$ps_mod, mods$outcome_mod, .data = dat), "ipw")
})
