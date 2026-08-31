# Tests for the internal psi-building layer that backs the M-estimation path of
# ipw(): the weight registry (ipw_weight_fn), the theta layout (ipw_theta_layout),
# and the stacked estimating-function builder (build_ipw_psi). Specs are
# constructed directly as plain lists following the M-estimation design contract.

# ---- data simulators --------------------------------------------------------

sim_binary <- function(seed = 2024, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.6 * x2))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * x1 - 0.3 * x2))
  yc <- 1.5 + 0.8 * z + 0.6 * x1 - 0.4 * x2 + rnorm(n)
  data.frame(x1, x2, z, y, yc)
}

sim_categorical <- function(seed = 512, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  eta_b <- -0.2 + 0.5 * x1 + 0.3 * x2
  eta_c <- 0.1 - 0.4 * x1 + 0.6 * x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  pa <- 1 / denom
  pb <- exp(eta_b) / denom
  u <- runif(n)
  lab <- ifelse(u < pa, "a", ifelse(u < pa + pb, "b", "c"))
  a <- factor(lab, levels = c("a", "b", "c"))
  y <- rbinom(
    n,
    1,
    plogis(-0.3 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.5 * x1)
  )
  yc <- 0.5 + 0.4 * (a == "b") + 0.9 * (a == "c") + 0.6 * x1 + rnorm(n)
  data.frame(x1, x2, a, y, yc)
}

sim_continuous <- function(seed = 77, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 0.5 + 0.7 * x1 - 0.3 * x2 + rnorm(n)
  y <- 1 + 0.6 * A + 0.4 * x1 + rnorm(n)
  yb <- rbinom(n, 1, plogis(-0.5 + 0.5 * A + 0.3 * x1))
  data.frame(x1, x2, A, y, yb)
}

# The same simulation with an exposure that is positive everywhere, which a
# log-link propensity model needs to start its iteration. The dose is the one
# above shifted up rather than exponentiated, so it stays conditionally normal
# with a single spread. The columns are named as they are in `sim_continuous()`,
# so `continuous_spec()` reads either one.
sim_continuous_positive <- function(seed = 77, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 5 + 0.7 * x1 - 0.3 * x2 + rnorm(n)
  y <- 1 + 0.6 * A + 0.4 * x1 + rnorm(n)
  yb <- rbinom(n, 1, plogis(-3 + 0.5 * A + 0.3 * x1))
  data.frame(x1, x2, A, y, yb)
}

# ---- shared builders --------------------------------------------------------

# Counterfactual design matrix: set the exposure column to a fixed value and
# rebuild the model matrix, mirroring how estimate_marginal_means() constructs
# counterfactual data in R/ipw.R.
counterfactual_mm <- function(mod, dat, exposure_name, value) {
  d <- dat
  d[[exposure_name]] <- value
  stats::model.matrix(stats::delete.response(stats::terms(mod)), data = d)
}

counterfactual_mm_factor <- function(mod, dat, exposure_name, level) {
  d <- dat
  d[[exposure_name]] <- factor(level, levels = levels(dat[[exposure_name]]))
  stats::model.matrix(stats::delete.response(stats::terms(mod)), data = d)
}

binary_weights <- function(
  ps,
  z,
  estimand,
  stabilize = FALSE,
  stab_score = NULL
) {
  wt_fun <- switch(
    estimand,
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )
  args <- list(ps, z, exposure_type = "binary", .focal_level = 1)
  if (estimand == "ate") {
    args$stabilize <- stabilize
    args$stabilization_score <- stab_score
  }
  withr::with_options(
    list(propensity.quiet = TRUE),
    do.call(wt_fun, args)
  )
}

categorical_weights <- function(
  ps_named,
  a,
  estimand,
  focal_level = NULL,
  stabilize = FALSE,
  stab_score = NULL
) {
  wt_fun <- switch(
    estimand,
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )
  args <- list(ps_named, a, exposure_type = "categorical")
  if (estimand %in% c("att", "atu")) {
    args$.focal_level <- focal_level
  }
  if (estimand == "ate") {
    args$stabilize <- stabilize
    args$stabilization_score <- stab_score
  }
  do.call(wt_fun, args)
}

binary_spec <- function(
  dat,
  estimand,
  outcome_family = "binomial",
  ps_link = "logit",
  stabilize = FALSE,
  stab_score = NULL
) {
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial(link = ps_link))
  ps <- as.double(predict(ps_mod, type = "response"))
  z <- dat$z

  wts <- binary_weights(ps, z, estimand, stabilize, stab_score)

  outcome_var <- if (outcome_family == "binomial") "y" else "yc"
  out_fmla <- stats::reformulate(c("z", "x1"), response = outcome_var)
  if (outcome_family == "binomial") {
    out_mod <- glm(
      out_fmla,
      data = dat,
      family = quasibinomial(),
      weights = as.double(wts)
    )
    out_link <- "logit"
    contrasts <- c("rd", "log(rr)", "log(or)")
  } else {
    out_mod <- lm(out_fmla, data = dat, weights = as.double(wts))
    out_link <- "identity"
    contrasts <- "diff"
  }

  list(
    exposure_type = "binary",
    estimand = estimand,
    n = nrow(dat),
    exposure = z,
    ps = list(
      X = model.matrix(ps_mod),
      link = ps_link,
      coefs = coef(ps_mod),
      k = 2
    ),
    stab = list(stabilized = isTRUE(stabilize), score = stab_score),
    outcome = list(
      X = model.matrix(out_mod),
      y = dat[[outcome_var]],
      family = outcome_family,
      link = out_link,
      coefs = coef(out_mod),
      X_counterfactual = list(
        X1 = counterfactual_mm(out_mod, dat, "z", 1),
        X0 = counterfactual_mm(out_mod, dat, "z", 0)
      )
    ),
    contrasts = contrasts,
    focal_level = if (estimand %in% c("att", "atu")) 1 else NULL,
    reference_level = NULL
  )
}

categorical_spec <- function(
  dat,
  estimand,
  focal_level = NULL,
  stabilize = FALSE,
  stab_score = NULL,
  outcome_family = "binomial"
) {
  ps_mod <- nnet::multinom(
    a ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  lev <- levels(dat$a)
  k <- length(lev)
  ps_matrix <- unname(predict(ps_mod, type = "probs"))
  ps_named <- ps_matrix
  colnames(ps_named) <- lev
  z_ind <- vapply(
    lev,
    function(l) as.integer(dat$a == l),
    integer(nrow(dat))
  )

  wts <- categorical_weights(
    ps_named,
    dat$a,
    estimand,
    focal_level,
    stabilize,
    stab_score
  )

  outcome_var <- if (outcome_family == "binomial") "y" else "yc"
  out_fmla <- stats::reformulate(c("a", "x1"), response = outcome_var)
  if (outcome_family == "binomial") {
    out_mod <- glm(
      out_fmla,
      data = dat,
      family = quasibinomial(),
      weights = as.double(wts)
    )
    out_link <- "logit"
    contrasts <- c("rd", "log(rr)", "log(or)")
  } else {
    out_mod <- lm(out_fmla, data = dat, weights = as.double(wts))
    out_link <- "identity"
    contrasts <- "diff"
  }
  x_cf <- lapply(lev, function(l) {
    counterfactual_mm_factor(out_mod, dat, "a", l)
  })
  names(x_cf) <- lev

  list(
    exposure_type = "categorical",
    estimand = estimand,
    n = nrow(dat),
    exposure = z_ind,
    ps = list(
      X = model.matrix(ps_mod),
      link = NULL,
      coefs = as.vector(t(coef(ps_mod))),
      k = k
    ),
    stab = list(stabilized = isTRUE(stabilize), score = stab_score),
    outcome = list(
      X = model.matrix(out_mod),
      y = dat[[outcome_var]],
      family = outcome_family,
      link = out_link,
      coefs = coef(out_mod),
      X_counterfactual = x_cf
    ),
    contrasts = contrasts,
    focal_level = focal_level,
    reference_level = lev[1]
  )
}

# A continuous spec, built by hand the way `ipw_spec_continuous()` builds one.
# `ps_type` selects the propensity model and, with it, the equation the ps score
# is written as: `kind` says which of the registry's entries wrote it and `link`
# says how the conditional mean is read back from the coefficients. A robust fit
# carries one thing more, the threshold its Huber score clips at, which is held
# fixed rather than estimated.
#
# That threshold is read off the psi function itself. `rlm` tunes the psi by
# rewriting its formals, so `formals(fit$psi)$k` is the constant the fit
# actually clipped the standardized residual at, and the raw residual is clipped
# at that constant times the scale the fit settled on. The `k2` element is a
# different constant: it tunes the scale estimator rather than the psi, it is
# unchanged when the caller passes `k`, and it changes when the caller passes
# `k2` without the psi's threshold moving at all. `rlm_args` exists to fit both
# of those cases.
#
# `.density`, `numerator`, and `.sigma` are the three choices the weights
# record, and the spec reads the first two straight off that record. The spread
# is the exception: the record says only where it came from, so a fixed one is
# carried here as the number it was.
#
# A fitted numerator model reaches `wt_ate()` through `stabilize` itself, so it
# arrives here the same way. The spec carries it as the block the stacked system
# estimates it in, described the way the propensity score block describes its
# own model: the design, the entry that writes its score, the link its mean is
# read back through, and the coefficients it was fit at.
continuous_spec <- function(
  dat,
  stabilize = TRUE,
  stab_score = NULL,
  outcome_family = "gaussian",
  ps_type = c("lm", "glm_log", "rlm"),
  rlm_args = list(),
  .density = "normal",
  numerator = "marginal",
  .sigma = NULL
) {
  ps_type <- match.arg(ps_type)
  num_mod <- if (inherits(stabilize, "lm")) stabilize
  ps_mod <- switch(
    ps_type,
    lm = lm(A ~ x1 + x2, data = dat),
    glm_log = glm(
      A ~ x1 + x2,
      data = dat,
      family = gaussian(link = "log"),
      control = glm.control(epsilon = 1e-14, maxit = 200)
    ),
    rlm = do.call(
      MASS::rlm,
      c(list(A ~ x1 + x2, data = dat, acc = 1e-10), rlm_args)
    )
  )
  ps_kind <- if (ps_type == "glm_log") "glm" else ps_type
  ps_link <- if (ps_type == "glm_log") "log" else "identity"
  huber_k <- if (ps_type == "rlm") {
    as.numeric(formals(ps_mod$psi)$k) * ps_mod$s
  }
  fitted_ps <- as.double(fitted(ps_mod))
  A <- dat$A
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      fitted_ps,
      A,
      .sigma = .sigma,
      exposure_type = "continuous",
      stabilize = stabilize,
      stabilization_score = stab_score,
      .density = .density,
      numerator = numerator
    )
  )
  meta <- density_meta(wts)

  outcome_var <- if (outcome_family == "binomial") "yb" else "y"
  msm_fmla <- stats::reformulate("A", response = outcome_var)
  if (outcome_family == "binomial") {
    msm <- glm(
      msm_fmla,
      data = dat,
      family = quasibinomial(),
      weights = as.double(wts)
    )
    out_link <- "logit"
  } else {
    msm <- lm(msm_fmla, data = dat, weights = as.double(wts))
    out_link <- "identity"
  }

  list(
    exposure_type = "continuous",
    estimand = "ate",
    n = nrow(dat),
    exposure = A,
    ps = list(
      X = model.matrix(ps_mod),
      kind = ps_kind,
      link = ps_link,
      huber_k = huber_k,
      coefs = coef(ps_mod),
      k = NULL
    ),
    stab = list(
      stabilized = isTRUE(stabilize) || !is.null(num_mod),
      score = stab_score,
      model = if (!is.null(num_mod)) {
        list(
          X = model.matrix(num_mod),
          kind = if (inherits(num_mod, "glm")) "glm" else "lm",
          link = if (inherits(num_mod, "glm")) num_mod$family$link else
            "identity",
          coefs = coef(num_mod)
        )
      }
    ),
    density = meta$density,
    numerator = meta$numerator,
    sigma = if (is.null(.sigma)) {
      list(kind = "pooled", value = NULL)
    } else {
      list(kind = "fixed", value = .sigma)
    },
    outcome = list(
      X = model.matrix(msm),
      y = dat[[outcome_var]],
      family = outcome_family,
      link = out_link,
      coefs = coef(msm),
      X_counterfactual = NULL
    ),
    contrasts = NULL,
    focal_level = NULL,
    reference_level = NULL
  )
}

# ---- assertion helpers ------------------------------------------------------

call_weight_fn <- function(
  exposure_type,
  estimand,
  ps,
  exposure,
  extras,
  silent = FALSE
) {
  fn <- ipw_weight_fn(exposure_type, estimand)
  expect_type(fn, "closure")
  if (silent) {
    return(expect_silent(fn(ps, exposure, extras)))
  }
  fn(ps, exposure, extras)
}

expect_weight_parity <- function(
  exposure_type,
  estimand,
  ps,
  exposure,
  extras,
  expected,
  silent = FALSE
) {
  got <- call_weight_fn(
    exposure_type,
    estimand,
    ps,
    exposure,
    extras,
    silent = silent
  )
  expect_equal(as.double(got), as.double(expected), tolerance = 1e-12)
}

expect_layout_partition <- function(layout, sizes) {
  block_order <- c(
    "ps",
    "stab",
    "out",
    "mu",
    "contrast",
    "by_mu",
    "by_contrast"
  )
  expect_named(layout$idx, block_order)
  expect_equal(
    lengths(layout$idx[block_order]),
    vapply(block_order, function(b) as.integer(sizes[[b]]), integer(1))
  )
  combined <- unlist(layout$idx[block_order], use.names = FALSE)
  expect_equal(combined, seq_along(layout$init))
  nm <- names(layout$init)
  expect_equal(length(nm), length(layout$init))
  expect_false(any(is.na(nm) | nm == ""))
}

expect_root_seeded <- function(spec) {
  layout <- ipw_theta_layout(spec)
  psi <- build_ipw_psi(spec, layout)
  mat <- psi(layout$init)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(length(layout$init), spec$n))
  expect_false(anyNA(mat))
  # init is the exact root, so the mean per-observation contribution of every
  # estimating equation is ~0. Normalize by n: the multinom ps block cannot
  # beat an unnormalized 1e-6 because nnet's optimizer stalls near 2e-6.
  expect_true(all(abs(rowSums(mat)) / spec$n < 1e-8))
}

# ---- weight registry: binary ------------------------------------------------

test_that("ipw_weight_fn reproduces binary weight functions at fitted params", {
  dat <- sim_binary()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  ps <- as.double(predict(ps_mod, type = "response"))
  z <- dat$z

  no_extras <- list(stab_prob = NULL, score = NULL)
  for (estimand in ipw_estimands) {
    expect_weight_parity(
      "binary",
      estimand,
      ps,
      z,
      no_extras,
      binary_weights(ps, z, estimand)
    )
  }
})

test_that("ipw_weight_fn binary ate matches across ps links", {
  dat <- sim_binary()
  z <- dat$z
  no_extras <- list(stab_prob = NULL, score = NULL)

  for (link in c("logit", "probit", "cloglog")) {
    ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial(link = link))
    ps <- as.double(predict(ps_mod, type = "response"))
    expect_weight_parity(
      "binary",
      "ate",
      ps,
      z,
      no_extras,
      binary_weights(ps, z, "ate")
    )
  }
})

test_that("ipw_weight_fn reproduces binary stabilized ate weights", {
  dat <- sim_binary()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  ps <- as.double(predict(ps_mod, type = "response"))
  z <- dat$z

  # default stabilizer pi = mean(z)
  expect_weight_parity(
    "binary",
    "ate",
    ps,
    z,
    list(stab_prob = mean(z), score = NULL),
    binary_weights(ps, z, "ate", stabilize = TRUE)
  )

  # fixed user-supplied stabilization score
  score <- 0.37
  expect_weight_parity(
    "binary",
    "ate",
    ps,
    z,
    list(stab_prob = NULL, score = score),
    binary_weights(ps, z, "ate", stabilize = TRUE, stab_score = score)
  )
})

test_that("binary entropy tilt stays finite at saturated propensity scores", {
  # The solver reaches saturated scores on its own: plogis() returns exactly 1
  # once the linear predictor exceeds about 36.7, and exactly 0 once it falls
  # below about -710, where exp(-eta) overflows to Inf.
  expect_identical(stats::plogis(40), 1)

  e <- c(0, 1e-3, 0.5, 1 - 1e-3, 1)
  tilt <- ps_tilt_binary(e, "entropy")

  expect_true(all(is.finite(tilt)))

  # The guard must leave every score strictly inside (0, 1) bit for bit.
  interior <- e[2:4]
  expect_identical(
    tilt[2:4],
    -interior * log(interior) - (1 - interior) * log(1 - interior)
  )

  # Saturated scores take the same clamp as the categorical tilt, which
  # substitutes .Machine$double.eps for a zero probability. The binary form
  # carries an extra -(1 - e) log(1 - e) term worth one machine epsilon at the
  # clamp, so the two agree on an absolute scale rather than exactly.
  degenerate <- ps_tilt_categorical(cbind(c(0, 1), c(1, 0)), "entropy")
  expect_lt(max(abs(tilt[c(1, 5)] - degenerate)), 1e-15)
  expect_identical(tilt[[1]], tilt[[5]])
})

# ---- tilt behavior at degenerate propensity scores --------------------------
#
# The tilts implement mathematical limits, and at a propensity score of exactly
# zero or one those limits are the values the arithmetic already produces. att's
# tilt at e = 0 is 0 because the weight on an unexposed unit's contribution to
# the treated population really is zero there; ato and atm vanish for the same
# reason. None of that needs correcting, and a clamp applied to them would move
# the answer away from the limit rather than toward it.
#
# entropy is the single exception, and the reason is specific: its -e log(e)
# form is 0 * log(0) at the boundary, an indeterminate whose limit is 0 but whose
# raw arithmetic is NaN. The clamp exists to resolve that one indeterminate. The
# tests below pin the others as unclamped so a future tidy-up cannot hoist the
# clamp out of the entropy branch and quietly change four estimands.
#
# The entropy branch's own behavior is pinned above by "binary entropy tilt stays
# finite at saturated propensity scores".

test_that("binary tilts other than entropy pass boundary scores through unclamped", {
  e <- c(0, 1)

  # Each is the limit, reached by the ordinary arithmetic.
  expect_identical(ps_tilt_binary(e, "att"), c(0, 1))
  expect_identical(ps_tilt_binary(e, "atu"), c(1, 0))
  expect_identical(ps_tilt_binary(e, "atm"), c(0, 0))
  expect_identical(ps_tilt_binary(e, "ato"), c(0, 0))
  expect_identical(ps_tilt_binary(e, "ate"), c(1, 1))

  # entropy is the exception: clamped, so strictly positive rather than the 0
  # the other tilts return, and finite rather than the NaN the raw form gives.
  entropy <- ps_tilt_binary(e, "entropy")
  expect_true(all(is.finite(entropy)))
  expect_true(all(entropy > 0))
  expect_true(is.nan(-0 * log(0) - 1 * log(1)))
})

test_that("categorical tilts other than entropy pass degenerate rows through unclamped", {
  # Rows one and three are degenerate: all the mass sits on a single level.
  ps <- rbind(c(0, 1, 0), c(0.2, 0.3, 0.5), c(1, 0, 0))

  # ato reaches exactly zero through an infinite row sum rather than by a guard.
  expect_identical(sum(1 / ps[1, ]), Inf)
  ato <- ps_tilt_categorical(ps, "ato")
  expect_identical(ato[c(1, 3)], c(0, 0))
  expect_gt(ato[[2]], 0)

  expect_identical(ps_tilt_categorical(ps, "atm")[c(1, 3)], c(0, 0))

  # att and atu read the focal column as it stands, boundary values included.
  expect_identical(
    ps_tilt_categorical(ps, "att", focal_idx = 2),
    c(1, 0.3, 0)
  )
  expect_identical(
    ps_tilt_categorical(ps, "atu", focal_idx = 2),
    c(0, 0.7, 1)
  )

  # entropy again the exception: clamped to a small positive number where the
  # unclamped tilts give exactly zero.
  entropy <- ps_tilt_categorical(ps, "entropy")
  expect_true(all(is.finite(entropy)))
  expect_true(all(entropy > 0))
})

test_that("categorical weights stay finite when an off-level score is zero", {
  # Only the score at the level a unit was actually assigned divides its weight,
  # so a zero somewhere else in the row is not a problem. This is what makes the
  # separation guard's condition the observed-level score rather than any zero in
  # the matrix; see "mestimation reports separation in the propensity model as
  # its own error" in test-ipw-mestimation.R.
  ps <- rbind(c(0.6, 0, 0.4), c(0.2, 0.3, 0.5))
  exposure <- rbind(c(1, 0, 0), c(1, 0, 0))
  expect_identical(rowSums(exposure * ps), c(0.6, 0.2))

  for (estimand in ipw_estimands) {
    weights <- ipw_weight_fn("categorical", estimand)(
      ps,
      exposure,
      list(focal_idx = 1)
    )
    expect_true(all(is.finite(weights)))
  }
})

test_that("a zero observed-level score gives a non-finite weight at the registry", {
  # Called directly, with no guard in front of it, the weight is the tilt divided
  # by the observed-level score, and dividing by zero is what it is. Which
  # non-finite value appears depends on the tilt: a nonzero numerator gives Inf,
  # and a tilt that vanishes at the same row gives NaN from 0/0. Pinned as the
  # arithmetic rather than as a contract anyone should rely on, since through
  # ipw() the separation guard rejects this before the registry is reached.
  ps <- rbind(c(0, 1, 0), c(0.2, 0.3, 0.5))
  exposure <- rbind(c(1, 0, 0), c(1, 0, 0))
  expect_identical(rowSums(exposure * ps)[[1]], 0)

  weight_at_row1 <- function(estimand) {
    ipw_weight_fn("categorical", estimand)(
      ps,
      exposure,
      list(focal_idx = 2)
    )[[1]]
  }

  # Tilts that are nonzero on this row: ate is flat, att reads the focal column
  # (one here), and entropy's clamp leaves it small but positive.
  expect_identical(weight_at_row1("ate"), Inf)
  expect_identical(weight_at_row1("att"), Inf)
  expect_identical(weight_at_row1("entropy"), Inf)

  # Tilts that vanish on this row divide zero by zero.
  expect_true(is.nan(weight_at_row1("atu")))
  expect_true(is.nan(weight_at_row1("atm")))
  expect_true(is.nan(weight_at_row1("ato")))
})

# ---- weight registry: categorical -------------------------------------------

test_that("ipw_weight_fn reproduces categorical weight functions at fitted params", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical()
  ps_mod <- nnet::multinom(a ~ x1 + x2, data = dat, trace = FALSE)
  ps_matrix <- unname(predict(ps_mod, type = "probs"))
  lev <- levels(dat$a)
  ps_named <- ps_matrix
  colnames(ps_named) <- lev
  z_ind <- vapply(lev, function(l) as.integer(dat$a == l), integer(nrow(dat)))

  focal <- "b"
  focal_idx <- which(lev == focal)

  configs <- list(
    ate = list(
      extras = list(focal_idx = NULL, stab_probs = NULL, score = NULL),
      oracle = categorical_weights(ps_named, dat$a, "ate")
    ),
    att = list(
      extras = list(focal_idx = focal_idx, stab_probs = NULL, score = NULL),
      oracle = categorical_weights(ps_named, dat$a, "att", focal_level = focal)
    ),
    atu = list(
      extras = list(focal_idx = focal_idx, stab_probs = NULL, score = NULL),
      oracle = categorical_weights(ps_named, dat$a, "atu", focal_level = focal)
    ),
    atm = list(
      extras = list(focal_idx = NULL, stab_probs = NULL, score = NULL),
      oracle = categorical_weights(ps_named, dat$a, "atm")
    ),
    ato = list(
      extras = list(focal_idx = NULL, stab_probs = NULL, score = NULL),
      oracle = categorical_weights(ps_named, dat$a, "ato")
    ),
    entropy = list(
      extras = list(focal_idx = NULL, stab_probs = NULL, score = NULL),
      oracle = categorical_weights(ps_named, dat$a, "entropy")
    )
  )

  for (estimand in names(configs)) {
    cfg <- configs[[estimand]]
    expect_weight_parity(
      "categorical",
      estimand,
      ps_matrix,
      z_ind,
      cfg$extras,
      cfg$oracle
    )
  }
})

test_that("ipw_weight_fn reproduces categorical stabilized ate weights", {
  skip_if_not_installed("nnet")
  dat <- sim_categorical()
  ps_mod <- nnet::multinom(a ~ x1 + x2, data = dat, trace = FALSE)
  ps_matrix <- unname(predict(ps_mod, type = "probs"))
  lev <- levels(dat$a)
  ps_named <- ps_matrix
  colnames(ps_named) <- lev
  z_ind <- vapply(lev, function(l) as.integer(dat$a == l), integer(nrow(dat)))

  marginal <- as.numeric(table(dat$a) / nrow(dat))

  # default marginal-probability stabilizer
  expect_weight_parity(
    "categorical",
    "ate",
    ps_matrix,
    z_ind,
    list(focal_idx = NULL, stab_probs = marginal, score = NULL),
    categorical_weights(ps_named, dat$a, "ate", stabilize = TRUE)
  )

  # fixed user-supplied stabilization score
  score <- 0.42
  expect_weight_parity(
    "categorical",
    "ate",
    ps_matrix,
    z_ind,
    list(focal_idx = NULL, stab_probs = NULL, score = score),
    categorical_weights(
      ps_named,
      dat$a,
      "ate",
      stabilize = TRUE,
      stab_score = score
    )
  )
})

# ---- categorical propensity score reconstruction ----------------------------

# Reference-first softmax written out without any shift, the reconstruction the
# psi layer performs. Used as an independent oracle rather than comparing
# ipw_categorical_ps() against itself. vapply() drops to a bare vector for a
# single-row design, so the linear predictors are reshaped explicitly.
# Assumes at least one design row: the vapply template is sized from nrow(x),
# so a zero-row design would not round-trip through the matrix() call below.
unshifted_softmax <- function(x, theta, k) {
  b <- ncol(x)
  eta <- matrix(
    vapply(
      seq_len(k - 1),
      function(j) as.vector(x %*% theta[((j - 1) * b + 1):(j * b)]),
      numeric(nrow(x))
    ),
    nrow = nrow(x)
  )
  num <- cbind(1, exp(eta))
  num / rowSums(num)
}

test_that("ipw_categorical_ps stays finite for extreme linear predictors", {
  # Level-major, term-minor theta for a three-level exposure: category b loads
  # only on xa, category c loads only on xb, so the two non-reference linear
  # predictors can be driven apart independently.
  x <- cbind(
    `(Intercept)` = 1,
    xa = c(100, 0, -100, 100, 0.3),
    xb = c(0, 100, -100, 99, -0.2)
  )
  theta <- c(0, 10, 0, 0, 0, 10)

  ps <- ipw_categorical_ps(x, theta, 3)

  expect_true(all(is.finite(ps)))
  expect_true(all(ps >= 0 & ps <= 1))
  expect_lt(max(abs(rowSums(ps) - 1)), 1e-12)

  # Rows 1 and 2 saturate on a single category, row 3 saturates on the
  # reference, and row 4 has two large predictors ten apart, so its answer is a
  # genuine softmax rather than a degenerate 0/1 split.
  expect_equal(ps[1, ], c(0, 1, 0))
  expect_equal(ps[2, ], c(0, 0, 1))
  expect_equal(ps[3, ], c(1, 0, 0))
  expect_equal(
    ps[4, ],
    c(0, 1, exp(-10)) / (1 + exp(-10)),
    tolerance = 1e-15
  )

  # Four levels, one category per design column, so the running maximum has to
  # accumulate over three non-reference categories rather than track any single
  # one of them. Row 1 has linear predictors (0, 0, 0, 1000), putting the row
  # maximum in the third non-reference category; row 2 has (0, 1000, 995, 0),
  # putting it in the first with a genuine softmax split behind it.
  x4 <- cbind(
    `(Intercept)` = 1,
    xa = c(0, 100),
    xb = c(0, 99.5),
    xc = c(100, 0)
  )
  theta4 <- c(0, 10, 0, 0, 0, 0, 10, 0, 0, 0, 0, 10)

  ps4 <- ipw_categorical_ps(x4, theta4, 4)

  expect_true(all(is.finite(ps4)))
  expect_true(all(ps4 >= 0 & ps4 <= 1))
  expect_lt(max(abs(rowSums(ps4) - 1)), 1e-12)
  expect_equal(ps4[1, ], c(0, 0, 0, 1))
  expect_equal(
    ps4[2, ],
    c(0, 1, exp(-5), 0) / (1 + exp(-5)),
    tolerance = 1e-15
  )
})

test_that("ipw_categorical_ps is unchanged for moderate linear predictors", {
  withr::local_seed(918)
  x <- cbind(1, rnorm(6), rnorm(6))
  theta <- c(-0.2, 0.5, 0.3, 0.1, -0.4, 0.6)

  expect_equal(
    ipw_categorical_ps(x, theta, 3),
    unshifted_softmax(x, theta, 3),
    tolerance = 1e-15
  )

  # When the reference category already holds the largest linear predictor the
  # shift is exactly zero, so the reconstruction is unchanged bit for bit.
  theta_ref <- c(-3, 0.2, 0.1, -4, -0.1, 0.05)
  expect_identical(
    ipw_categorical_ps(x, theta_ref, 3),
    unshifted_softmax(x, theta_ref, 3)
  )
})

# ---- weight registry: continuous --------------------------------------------

test_that("ipw_weight_fn reproduces continuous ate weights at fitted params", {
  dat <- sim_continuous()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  fitted_ps <- as.double(fitted(ps_mod))
  A <- dat$A

  sigma2_d <- mean((A - fitted_ps)^2)
  mu_a <- mean(A)
  sigma2_a <- mean((A - mu_a)^2)

  oracle_unstab <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(fitted_ps, A, exposure_type = "continuous", stabilize = FALSE)
  )
  oracle_stab <- wt_ate(
    fitted_ps,
    A,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  score <- 0.5
  oracle_fixed <- wt_ate(
    fitted_ps,
    A,
    exposure_type = "continuous",
    stabilize = TRUE,
    stabilization_score = score
  )

  # unstabilized: the registry reproduces the base weight and stays silent
  # (it must not emit the alert_info that wt_ate() does)
  expect_weight_parity(
    "continuous",
    "ate",
    fitted_ps,
    A,
    list(
      sigma2_d = sigma2_d,
      mu_a = NULL,
      sigma2_a = NULL,
      score = NULL,
      stabilized = FALSE
    ),
    oracle_unstab,
    silent = TRUE
  )

  expect_weight_parity(
    "continuous",
    "ate",
    fitted_ps,
    A,
    list(
      sigma2_d = sigma2_d,
      mu_a = mu_a,
      sigma2_a = sigma2_a,
      score = NULL,
      stabilized = TRUE
    ),
    oracle_stab
  )

  expect_weight_parity(
    "continuous",
    "ate",
    fitted_ps,
    A,
    list(
      sigma2_d = sigma2_d,
      mu_a = NULL,
      sigma2_a = NULL,
      score = score,
      stabilized = TRUE
    ),
    oracle_fixed
  )
})

test_that("ipw_weight_fn errors on unsupported estimand/exposure combinations", {
  expect_error(
    ipw_weight_fn("continuous", "att"),
    class = "propensity_ipw_estimand_error"
  )
  expect_error(
    ipw_weight_fn("continuous", "ato"),
    class = "propensity_ipw_estimand_error"
  )
})

# ---- theta layout -----------------------------------------------------------

test_that("ipw_theta_layout partitions theta for binary specs", {
  # stabilized ate: stab block present
  spec_stab <- binary_spec(sim_binary(), "ate", stabilize = TRUE)
  layout_stab <- ipw_theta_layout(spec_stab)
  expect_layout_partition(
    layout_stab,
    list(
      ps = ncol(spec_stab$ps$X),
      stab = 1,
      out = ncol(spec_stab$outcome$X),
      mu = 2,
      contrast = 3,
      by_mu = 0,
      by_contrast = 0
    )
  )

  # unstabilized ate: no stab block
  spec_unstab <- binary_spec(sim_binary(), "ate")
  layout_unstab <- ipw_theta_layout(spec_unstab)
  expect_layout_partition(
    layout_unstab,
    list(
      ps = ncol(spec_unstab$ps$X),
      stab = 0,
      out = ncol(spec_unstab$outcome$X),
      mu = 2,
      contrast = 3,
      by_mu = 0,
      by_contrast = 0
    )
  )

  # gaussian outcome: single linear contrast
  spec_lin <- binary_spec(sim_binary(), "att", outcome_family = "gaussian")
  layout_lin <- ipw_theta_layout(spec_lin)
  expect_layout_partition(
    layout_lin,
    list(
      ps = ncol(spec_lin$ps$X),
      stab = 0,
      out = ncol(spec_lin$outcome$X),
      mu = 2,
      contrast = 1,
      by_mu = 0,
      by_contrast = 0
    )
  )
})

test_that("ipw_theta_layout partitions theta for categorical specs", {
  skip_if_not_installed("nnet")

  spec <- categorical_spec(sim_categorical(), "ate")
  k <- spec$ps$k
  layout <- ipw_theta_layout(spec)
  expect_layout_partition(
    layout,
    list(
      ps = ncol(spec$ps$X) * (k - 1),
      stab = 0,
      out = ncol(spec$outcome$X),
      mu = k,
      contrast = (k - 1) * length(spec$contrasts),
      by_mu = 0,
      by_contrast = 0
    )
  )

  spec_stab <- categorical_spec(sim_categorical(), "ate", stabilize = TRUE)
  layout_stab <- ipw_theta_layout(spec_stab)
  expect_layout_partition(
    layout_stab,
    list(
      ps = ncol(spec_stab$ps$X) * (k - 1),
      stab = k - 1,
      out = ncol(spec_stab$outcome$X),
      mu = k,
      contrast = (k - 1) * length(spec_stab$contrasts),
      by_mu = 0,
      by_contrast = 0
    )
  )
})

test_that("ipw_theta_layout partitions theta for a continuous stabilized spec", {
  spec <- continuous_spec(sim_continuous(), stabilize = TRUE)
  layout <- ipw_theta_layout(spec)
  expect_layout_partition(
    layout,
    list(
      ps = ncol(spec$ps$X) + 1,
      stab = 2,
      out = ncol(spec$outcome$X),
      mu = 0,
      contrast = 0,
      by_mu = 0,
      by_contrast = 0
    )
  )
})

test_that("ipw_theta_layout gives a numerator model its own stabilization block", {
  dat <- sim_continuous()
  num_mod <- lm(A ~ x1, data = dat)
  spec <- continuous_spec(dat, stabilize = num_mod)
  layout <- ipw_theta_layout(spec)

  # The spec reads its numerator off the weights, so this is what the record
  # says as much as it is what the spec holds.
  expect_identical(spec$numerator, "model")

  # The numerator model is estimated where every other stabilizer is, in the
  # stab block, and it is the same shape as the propensity score block it
  # divides into: one parameter per coefficient and one for the spread its
  # density is read at.
  expect_layout_partition(
    layout,
    list(
      ps = ncol(spec$ps$X) + 1,
      stab = ncol(spec$stab$model$X) + 1,
      out = ncol(spec$outcome$X),
      mu = 0,
      contrast = 0,
      by_mu = 0,
      by_contrast = 0
    )
  )

  # The numerator's coefficients carry the block's prefix, which keeps them
  # apart from the propensity score coefficients over the same terms, and its
  # spread is named against `sigma2_d` the way it sits against it.
  expect_equal(
    names(layout$init)[layout$idx$stab],
    c(paste0("stab_", colnames(spec$stab$model$X)), "sigma2_n")
  )
  expect_equal(
    unname(layout$init[layout$idx$stab]),
    c(unname(coef(num_mod)), mean(residuals(num_mod)^2))
  )
})

test_that("ipw_theta_layout partitions theta for a continuous unstabilized spec", {
  # the ps block still carries the sigma2_d variance parameter (ncol(X) + 1),
  # but the stabilization block is empty
  spec <- continuous_spec(sim_continuous(), stabilize = FALSE)
  layout <- ipw_theta_layout(spec)
  expect_layout_partition(
    layout,
    list(
      ps = ncol(spec$ps$X) + 1,
      stab = 0,
      out = ncol(spec$outcome$X),
      mu = 0,
      contrast = 0,
      by_mu = 0,
      by_contrast = 0
    )
  )
})

# ---- build_ipw_psi: init is the exact root ----------------------------------

test_that("build_ipw_psi is root-seeded at init for binary ate binomial outcome", {
  expect_root_seeded(binary_spec(sim_binary(), "ate"))
})

test_that("build_ipw_psi is root-seeded at init for binary att gaussian outcome", {
  expect_root_seeded(binary_spec(
    sim_binary(),
    "att",
    outcome_family = "gaussian"
  ))
})

test_that("build_ipw_psi is root-seeded at init for binary stabilized ate", {
  expect_root_seeded(binary_spec(sim_binary(), "ate", stabilize = TRUE))
})

test_that("build_ipw_psi is root-seeded at init for categorical ate binomial outcome", {
  skip_if_not_installed("nnet")
  expect_root_seeded(categorical_spec(sim_categorical(), "ate"))
})

test_that("build_ipw_psi is root-seeded at init for categorical ate gaussian outcome", {
  skip_if_not_installed("nnet")
  expect_root_seeded(categorical_spec(
    sim_categorical(),
    "ate",
    outcome_family = "gaussian"
  ))
})

test_that("build_ipw_psi is root-seeded at init for continuous stabilized ate MSM", {
  expect_root_seeded(continuous_spec(sim_continuous(), stabilize = TRUE))
})

test_that("build_ipw_psi is root-seeded at init for continuous unstabilized ate MSM", {
  expect_root_seeded(continuous_spec(sim_continuous(), stabilize = FALSE))
})

test_that("build_ipw_psi is root-seeded at init for a numerator model", {
  # The numerator model's own equations sit at their root at the coefficients it
  # was fit at and at the second moment of its residuals, which is what makes
  # the stacked system carry the uncertainty of having fit it.
  dat <- sim_continuous()
  spec <- continuous_spec(dat, stabilize = lm(A ~ x1, data = dat))

  expect_identical(spec$numerator, "model")
  expect_root_seeded(spec)
})

test_that("build_ipw_psi is root-seeded at init for continuous stabilized logistic MSM", {
  expect_root_seeded(continuous_spec(
    sim_continuous(),
    stabilize = TRUE,
    outcome_family = "binomial"
  ))
})

# ---- the densities, numerators, and spreads the stacked system carries ------
#
# A continuous exposure's weights are a ratio of densities, and the stacked
# system has to rebuild that ratio from theta at every evaluation. Which
# parameters it needs depends on how the weights were built: the family and the
# spread decide the propensity block, and the numerator decides whether there is
# a stabilization block at all.

test_that("build_ipw_psi is root-seeded for a log-link propensity model", {
  expect_root_seeded(continuous_spec(
    sim_continuous_positive(),
    ps_type = "glm_log"
  ))
})

test_that("build_ipw_psi is root-seeded for a robust Huber propensity model", {
  skip_if_not_installed("MASS")
  expect_root_seeded(continuous_spec(sim_continuous(), ps_type = "rlm"))
})

test_that("build_ipw_psi is root-seeded for a robust fit tuned by its psi", {
  skip_if_not_installed("MASS")

  # A caller who passes `k` moves the psi's own threshold, and the fit settles
  # at the root of the score that threshold defines. Reading `k2` here would
  # write a different equation and leave the seed away from its root.
  expect_root_seeded(continuous_spec(
    sim_continuous(),
    ps_type = "rlm",
    rlm_args = list(k = 2)
  ))
})

test_that("build_ipw_psi is root-seeded for a robust fit tuned by its scale", {
  skip_if_not_installed("MASS")

  # A caller who passes `k2` moves the constant of the scale estimator and
  # leaves the psi's threshold where it was, so the score is the default one
  # read at whatever scale the fit settled on.
  expect_root_seeded(continuous_spec(
    sim_continuous(),
    ps_type = "rlm",
    rlm_args = list(k2 = 2)
  ))
})

test_that("build_ipw_psi is root-seeded for a t-density marginal ratio", {
  expect_root_seeded(continuous_spec(sim_continuous(), .density = dens_t(3)))
})

test_that("build_ipw_psi is root-seeded for a laplace-density marginal ratio", {
  expect_root_seeded(continuous_spec(sim_continuous(), .density = "laplace"))
})

test_that("build_ipw_psi is root-seeded for a normal integrated ratio", {
  expect_root_seeded(continuous_spec(
    sim_continuous(),
    numerator = "integrated"
  ))
})

test_that("build_ipw_psi is root-seeded for a t-density integrated ratio", {
  expect_root_seeded(continuous_spec(
    sim_continuous(),
    .density = dens_t(3),
    numerator = "integrated"
  ))
})

test_that("build_ipw_psi is root-seeded for a fixed scalar spread", {
  expect_root_seeded(continuous_spec(sim_continuous(), .sigma = 1.25))
})

test_that("build_ipw_psi is root-seeded for a density the user wrote", {
  expect_root_seeded(continuous_spec(
    sim_continuous(),
    .density = function(z) stats::dt(z, df = 5)
  ))
})

# ---- theta partitions for the continuous variants ---------------------------

test_that("ipw_theta_layout partitions theta for a log-link continuous spec", {
  # A link changes what the propensity score block means, not how much of it
  # there is: the coefficients and the one conditional variance.
  spec <- continuous_spec(sim_continuous_positive(), ps_type = "glm_log")
  layout <- ipw_theta_layout(spec)
  expect_layout_partition(
    layout,
    list(
      ps = ncol(spec$ps$X) + 1,
      stab = 2,
      out = ncol(spec$outcome$X),
      mu = 0,
      contrast = 0,
      by_mu = 0,
      by_contrast = 0
    )
  )
})

test_that("ipw_theta_layout gives an integrated numerator no stabilization block", {
  # The marginalized conditional density is a function of the propensity block
  # and the data, so there is nothing left for the stabilization block to
  # estimate; only a marginal numerator carries the exposure's own moments.
  spec <- continuous_spec(sim_continuous(), numerator = "integrated")
  layout <- ipw_theta_layout(spec)
  expect_layout_partition(
    layout,
    list(
      ps = ncol(spec$ps$X) + 1,
      stab = 0,
      out = ncol(spec$outcome$X),
      mu = 0,
      contrast = 0,
      by_mu = 0,
      by_contrast = 0
    )
  )
})

test_that("ipw_theta_layout stacks no variance parameter for a fixed spread", {
  # A spread the caller fixed is a constant the weights were built with rather
  # than something the data estimate, so the propensity block is the
  # coefficients alone.
  spec <- continuous_spec(sim_continuous(), .sigma = 1.25)
  layout <- ipw_theta_layout(spec)
  expect_layout_partition(
    layout,
    list(
      ps = ncol(spec$ps$X),
      stab = 2,
      out = ncol(spec$outcome$X),
      mu = 0,
      contrast = 0,
      by_mu = 0,
      by_contrast = 0
    )
  )
})

# ---- the rows the continuous blocks are ------------------------------------

test_that("the continuous stabilization block is deli's mean and variance", {
  # The marginal numerator is read at the exposure's own mean and variance, and
  # those two moments are exactly the estimating equations deli writes for them.
  spec <- continuous_spec(sim_continuous())
  layout <- ipw_theta_layout(spec)
  psi <- build_ipw_psi(spec, layout)

  rows <- psi(layout$init)[layout$idx$stab, , drop = FALSE]
  expect_equal(
    unname(rows),
    unname(deli::ee_mean_variance(
      layout$init[layout$idx$stab],
      y = spec$exposure
    ))
  )
})

test_that("the continuous ps block for a log link is deli's glm score", {
  # The conditional mean is exp(X alpha) under a log link, so the score is the
  # one deli writes for a normal glm with that link and the variance row reads
  # its residuals rather than the least-squares ones.
  spec <- continuous_spec(sim_continuous_positive(), ps_type = "glm_log")
  layout <- ipw_theta_layout(spec)
  psi <- build_ipw_psi(spec, layout)

  theta <- layout$init
  rows <- psi(theta)[layout$idx$ps, , drop = FALSE]

  n_alpha <- ncol(spec$ps$X)
  alpha <- theta[layout$idx$ps][seq_len(n_alpha)]
  sigma2_d <- theta[layout$idx$ps][[n_alpha + 1]]
  mu <- exp(as.vector(spec$ps$X %*% alpha))
  expected <- rbind(
    deli::ee_glm(
      alpha,
      X = spec$ps$X,
      y = spec$exposure,
      distribution = "normal",
      link = "log"
    ),
    matrix((spec$exposure - mu)^2 - sigma2_d, nrow = 1)
  )

  expect_equal(unname(rows), unname(expected))
})

test_that("the continuous ps block for a robust fit is deli's Huber score", {
  skip_if_not_installed("MASS")

  # `rlm` clips the residual divided by its scale estimate at the psi's own
  # threshold, and deli's robust regression clips the raw residual at `k`, so
  # the two are the same equation when `k` is the product of the two. The scale
  # is the one the fit settled on, held fixed: nothing in theta reproduces it,
  # and the variance row below it is the pooled residual moment the density
  # ratio divides by, which is a different quantity entirely.
  spec <- continuous_spec(sim_continuous(), ps_type = "rlm")
  layout <- ipw_theta_layout(spec)
  psi <- build_ipw_psi(spec, layout)

  theta <- layout$init
  rows <- psi(theta)[layout$idx$ps, , drop = FALSE]

  n_alpha <- ncol(spec$ps$X)
  alpha <- theta[layout$idx$ps][seq_len(n_alpha)]
  sigma2_d <- theta[layout$idx$ps][[n_alpha + 1]]
  mu <- as.vector(spec$ps$X %*% alpha)
  expected <- rbind(
    deli::ee_robust_regression(
      alpha,
      X = spec$ps$X,
      y = spec$exposure,
      model = "linear",
      k = spec$ps$huber_k,
      loss = "huber"
    ),
    matrix((spec$exposure - mu)^2 - sigma2_d, nrow = 1)
  )

  expect_equal(unname(rows), unname(expected))
})

test_that("the robust ps block reads the psi's threshold, not the scale constant", {
  skip_if_not_installed("MASS")

  # The two constants `rlm` reports are not the same number and do not tune the
  # same thing. `formals(fit$psi)$k` is where the psi clips the standardized
  # residual, and `fit$k2` is the tuning constant of the scale estimator the fit
  # iterates alongside it. They agree only at the defaults, so a spec built from
  # either one looks correct until a caller changes one of them.
  expect_ps_block_huber <- function(spec) {
    layout <- ipw_theta_layout(spec)
    psi <- build_ipw_psi(spec, layout)

    theta <- layout$init
    rows <- psi(theta)[layout$idx$ps, , drop = FALSE]

    n_alpha <- ncol(spec$ps$X)
    alpha <- theta[layout$idx$ps][seq_len(n_alpha)]
    sigma2_d <- theta[layout$idx$ps][[n_alpha + 1]]
    mu <- as.vector(spec$ps$X %*% alpha)
    expected <- rbind(
      deli::ee_robust_regression(
        alpha,
        X = spec$ps$X,
        y = spec$exposure,
        model = "linear",
        k = spec$ps$huber_k,
        loss = "huber"
      ),
      matrix((spec$exposure - mu)^2 - sigma2_d, nrow = 1)
    )

    expect_equal(unname(rows), unname(expected))
  }

  dat <- sim_continuous()

  # A rewritten psi threshold. `rlm` records it by rewriting the psi's formals
  # rather than by storing it anywhere else, so `fit$k2` still reads its default
  # here and a block built from that constant would be the score of a fit nobody
  # made.
  tuned_psi_fit <- MASS::rlm(A ~ x1 + x2, data = dat, acc = 1e-10, k = 2)
  expect_equal(as.numeric(formals(tuned_psi_fit$psi)$k), 2)
  expect_equal(tuned_psi_fit$k2, 1.345)
  expect_ps_block_huber(continuous_spec(
    dat,
    ps_type = "rlm",
    rlm_args = list(k = 2)
  ))

  # A rewritten scale constant. The psi clips where it always did, so the block
  # is the default score read at the scale this fit settled on, and a block
  # built from `k2` would clip somewhere the fit never did.
  tuned_scale_fit <- MASS::rlm(A ~ x1 + x2, data = dat, acc = 1e-10, k2 = 2)
  expect_equal(as.numeric(formals(tuned_scale_fit$psi)$k), 1.345)
  expect_equal(tuned_scale_fit$k2, 2)
  expect_ps_block_huber(continuous_spec(
    dat,
    ps_type = "rlm",
    rlm_args = list(k2 = 2)
  ))
})

test_that("a degenerate maximum likelihood scale is refused in the caller's name", {
  # The seed of the variance parameter is the scale the weights were built at,
  # so the estimator that refuses residuals it has no maximum for refuses here
  # too, and the user reads the refusal against the function they called rather
  # than against the seed.
  cnd <- rlang::catch_cnd(
    ipw_continuous_sigma2_seed(
      list(kind = "mle"),
      c(rep(0, 9), 1),
      dens_t(4, sigma_method = "mle"),
      call = rlang::call2("ipw")
    ),
    classes = "propensity_density_error"
  )

  expect_identical(conditionCall(cnd), quote(ipw()))
})
