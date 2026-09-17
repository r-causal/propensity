# Refitting a trimmed dose model, and the weights its record describes
#
# A density trim keeps the conditional means of the units whose observed dose
# the model finds plausible, together with the family and the spread the
# conditional density was read at. `ps_refit()` refits the mean on the retained
# rows and re-estimates the spread from their residuals, and the weight
# functions read the family and the spread from the record rather than from
# their own arguments. Every expectation below is held against a refit and a
# density ratio worked out by hand in helper-dose-trim.R.

# The dose model, its density trim, and the same model refit by hand on the
# rows the trim kept.
dose_refit_fixture <- function(
  density = dens_normal(),
  family = "normal",
  sigma_method = "rms",
  df = NULL,
  lower = 0.05
) {
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)
  trimmed <- ps_trim(fit, method = "density", lower = lower, .density = density)
  oracle <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "density",
    family = family,
    sigma_method = sigma_method,
    df = df,
    lower = lower
  )
  keep <- oracle$keep_idx

  refit_fit <- lm(a ~ x1 + x2, data = dat[keep, , drop = FALSE])
  refit_mu <- rep(NA_real_, nrow(dat))
  refit_mu[keep] <- fitted(refit_fit)

  full_mu <- as.numeric(fitted(fit))
  full_mu[oracle$trimmed_idx] <- NA_real_

  list(
    dat = dat,
    n = nrow(dat),
    fit = fit,
    trimmed = trimmed,
    oracle = oracle,
    keep = keep,
    refit_fit = refit_fit,
    refit_mu = refit_mu,
    full_mu = full_mu,
    refit_sigma = if (family %in% c("normal", "t", "laplace")) {
      dose_trim_sigma(
        dat$a[keep] - fitted(refit_fit),
        family,
        sigma_method,
        df
      )
    }
  )
}

# ---- ps_refit(): the mean and the spread on the retained rows ----------------

test_that("ps_refit() refits a density trim's mean and spread on the retained rows", {
  cases <- list(
    normal = list(density = dens_normal(), family = "normal"),
    t_mle = list(
      density = dens_t(4),
      family = "t",
      sigma_method = "mle",
      df = 4
    ),
    t_rms = list(
      density = dens_t(4, sigma_method = "rms"),
      family = "t",
      df = 4
    ),
    laplace = list(
      density = dens_laplace(),
      family = "laplace",
      sigma_method = "mle"
    )
  )

  for (case in names(cases)) {
    fx <- do.call(dose_refit_fixture, cases[[case]])
    expect_dose_trim(fx$trimmed, fx$oracle, fx$n)

    refit <- ps_refit(fx$trimmed, fx$fit)
    before <- ps_trim_meta(fx$trimmed)
    after <- ps_trim_meta(refit)

    expect_s3_class(refit, "ps_trim")
    expect_true(is_refit(refit), info = case)

    # The retained values are the refit model's conditional means, and the
    # trimmed units stay set aside.
    values <- as.numeric(refit)
    expect_true(all(is.na(values[fx$oracle$trimmed_idx])), info = case)
    expect_equal(
      values[fx$keep],
      fx$refit_mu[fx$keep],
      tolerance = 1e-10,
      info = case
    )

    # The spread is the family's own estimate on the retained residuals, which
    # is not the full-sample spread the trim was decided at.
    expect_equal(after$sigma, fx$refit_sigma, tolerance = 1e-8, info = case)
    expect_gt(abs(after$sigma - before$sigma), 1e-3)

    # What describes the cut that was made is left as it was.
    expect_identical(after$threshold, before$threshold, info = case)
    expect_identical(after$keep_idx, before$keep_idx, info = case)
    expect_identical(after$trimmed_idx, before$trimmed_idx, info = case)
    expect_identical(after$sigma_kind, before$sigma_kind, info = case)
    expect_identical(after$method, "density", info = case)
    expect_identical(after$lower, before$lower, info = case)
    expect_true(
      density_specs_agree(after$density, before$density),
      info = case
    )
  }
})

test_that("a refit density trim reads the density the refit model gives", {
  fx <- dose_refit_fixture(
    density = dens_t(4),
    family = "t",
    sigma_method = "mle",
    df = 4
  )
  refit <- ps_refit(fx$trimmed, fx$fit)
  meta <- ps_trim_meta(refit)

  # The density the refit record describes, read from its means, its spread,
  # and its family, is the inverse of the unstabilized weight the refit model
  # gives on the retained rows, which is how the simulation read it.
  a_kept <- fx$dat$a[fx$keep]
  from_record <- stats::dt(
    (a_kept - as.numeric(refit)[fx$keep]) / meta$sigma,
    df = 4
  ) /
    meta$sigma
  from_model <- 1 /
    as.numeric(wt_ate(
      fx$refit_fit,
      stabilize = FALSE,
      .density = dens_t(4)
    ))

  expect_equal(from_record, from_model, tolerance = 1e-10)
})

test_that("ps_refit() refits a gaussian glm trim on the retained rows", {
  dat <- sim_dose_trim()
  n <- nrow(dat)

  fit_glm <- glm(a ~ x1 + x2, data = dat, family = gaussian())
  trimmed <- ps_trim(
    fit_glm,
    method = "resid",
    upper = 2,
    .density = dens_laplace()
  )
  oracle <- dose_trim_oracle(
    dat$a,
    fitted(fit_glm),
    method = "resid",
    family = "laplace",
    sigma_method = "mle",
    upper = 2
  )
  expect_dose_trim(trimmed, oracle, n)

  refit <- ps_refit(trimmed, fit_glm)
  by_hand <- glm(
    a ~ x1 + x2,
    data = dat[oracle$keep_idx, ],
    family = gaussian()
  )
  expect_equal(
    as.numeric(refit)[oracle$keep_idx],
    as.numeric(fitted(by_hand)),
    tolerance = 1e-10
  )
  expect_equal(
    ps_trim_meta(refit)$sigma,
    mean(abs(residuals(by_hand, type = "response"))),
    tolerance = 1e-10
  )
  expect_identical(ps_trim_meta(refit)$upper, 2)
  expect_identical(
    ps_trim_meta(refit)$threshold,
    ps_trim_meta(trimmed)$threshold
  )
})

test_that("ps_refit() refits a robust dose model's trim", {
  skip_if_not_installed("MASS")
  dat <- sim_dose_trim()
  n <- nrow(dat)

  fit_rlm <- MASS::rlm(a ~ x1 + x2, data = dat)
  trimmed_rlm <- ps_trim(fit_rlm, method = "density", lower = 0.05)
  oracle_rlm <- dose_trim_oracle(
    dat$a,
    fitted(fit_rlm),
    method = "density",
    lower = 0.05
  )
  expect_dose_trim(trimmed_rlm, oracle_rlm, n)

  refit_rlm <- ps_refit(trimmed_rlm, fit_rlm)
  by_hand_rlm <- MASS::rlm(a ~ x1 + x2, data = dat[oracle_rlm$keep_idx, ])
  expect_equal(
    as.numeric(refit_rlm)[oracle_rlm$keep_idx],
    as.numeric(fitted(by_hand_rlm)),
    tolerance = 1e-8
  )
  expect_equal(
    ps_trim_meta(refit_rlm)$sigma,
    dose_trim_sigma(dat$a[oracle_rlm$keep_idx] - fitted(by_hand_rlm)),
    tolerance = 1e-8
  )
})

test_that("ps_refit() refits an additive dose model's trim", {
  skip_if_not_installed("mgcv")
  dat <- sim_dose_trim()
  n <- nrow(dat)

  fit_gam <- mgcv::gam(a ~ s(x1) + x2, data = dat, family = gaussian())
  trimmed_gam <- ps_trim(fit_gam, method = "density", lower = 0.05)
  oracle_gam <- dose_trim_oracle(
    dat$a,
    fitted(fit_gam),
    method = "density",
    lower = 0.05
  )
  expect_dose_trim(trimmed_gam, oracle_gam, n)

  refit_gam <- ps_refit(trimmed_gam, fit_gam)
  by_hand_gam <- mgcv::gam(
    a ~ s(x1) + x2,
    data = dat[oracle_gam$keep_idx, ],
    family = gaussian()
  )
  expect_equal(
    as.numeric(refit_gam)[oracle_gam$keep_idx],
    as.numeric(fitted(by_hand_gam)),
    tolerance = 1e-8
  )
  expect_equal(
    ps_trim_meta(refit_gam)$sigma,
    dose_trim_sigma(dat$a[oracle_gam$keep_idx] - fitted(by_hand_gam)),
    tolerance = 1e-8
  )
})

test_that("ps_refit() keeps a spread the caller supplied", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)

  trimmed <- ps_trim(fit, method = "resid", upper = 2, .sigma = 1.5)
  oracle <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "resid",
    sigma = 1.5,
    upper = 2
  )
  expect_dose_trim(trimmed, oracle, n)

  refit <- ps_refit(trimmed, fit)
  meta <- ps_trim_meta(refit)
  by_hand <- lm(a ~ x1 + x2, data = dat[oracle$keep_idx, ])

  expect_identical(meta$sigma, 1.5)
  expect_identical(meta$sigma_kind, "supplied")
  expect_true(is_refit(refit))
  expect_equal(
    as.numeric(refit)[oracle$keep_idx],
    as.numeric(fitted(by_hand)),
    tolerance = 1e-10
  )
})

test_that("ps_refit() reads a density trim's spread over the rows the refit analyzed", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)
  trimmed <- ps_trim(fit, method = "density", lower = 0.05)
  keep <- ps_trim_meta(trimmed)$keep_idx
  expect_gt(length(ps_trim_meta(trimmed)$trimmed_idx), 0)

  # A `subset` passed through narrows the rows the refit analyzes below the
  # rows the trim kept, and the spread describes the analyzed rows alone. It
  # indexes the kept rows, which are the data the refit is handed.
  positive <- dat$x1[keep] > 0
  expect_no_warning(refit <- ps_refit(trimmed, fit, subset = positive))
  by_hand <- lm(a ~ x1 + x2, data = dat[keep, ], subset = positive)
  expect_lt(length(fitted(by_hand)), length(keep))
  expect_equal(
    ps_trim_meta(refit)$sigma,
    sqrt(mean(residuals(by_hand)^2)),
    tolerance = 1e-10
  )

  # A covariate the refit formula adds, missing on some kept rows, leaves those
  # rows out of the refit, and so out of the spread.
  withr::local_seed(53)
  dat$x3 <- rnorm(n)
  dat$x3[keep[c(2, 5, 9)]] <- NA_real_
  expect_no_warning(
    refit_x3 <- ps_refit(trimmed, fit, .data = dat, formula. = ~ . + x3)
  )
  by_hand_x3 <- lm(a ~ x1 + x2 + x3, data = dat[keep, ])
  expect_equal(
    ps_trim_meta(refit_x3)$sigma,
    sqrt(mean(residuals(by_hand_x3)^2)),
    tolerance = 1e-10
  )
})

# ---- ps_refit(): refusals ----------------------------------------------------

test_that("ps_refit() refuses a model of level probabilities for a density trim", {
  skip_if_not_installed("nnet")
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)
  trimmed <- ps_trim(fit, method = "density", lower = 0.05)
  expect_gt(length(ps_trim_meta(trimmed)$trimmed_idx), 0)

  withr::local_seed(31)
  dat$z2 <- factor(sample(c("a", "b"), nrow(dat), replace = TRUE))
  dat$z3 <- factor(sample(c("a", "b", "c"), nrow(dat), replace = TRUE))
  fits <- list(
    two = nnet::multinom(z2 ~ x1 + x2, data = dat, trace = FALSE),
    three = nnet::multinom(z3 ~ x1 + x2, data = dat, trace = FALSE)
  )

  # A two-level `multinom` answers as a single probability rather than one per
  # level, and is still refused for a record of conditional means.
  for (kind in names(fits)) {
    expect_true(
      inherits(fits[[kind]], "multinom") || model_fits_levels(fits[[kind]]),
      info = kind
    )
    expect_error(
      ps_refit(trimmed, fits[[kind]]),
      class = "propensity_model_family_error",
      info = kind
    )
  }

  expect_propensity_error(ps_refit(trimmed, fits$two))
  expect_propensity_error(ps_refit(trimmed, fits$three))
})

test_that("ps_refit() refuses a model of a probability for a density trim and names every dose model", {
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)
  trimmed <- ps_trim(fit, method = "density", lower = 0.05)
  expect_gt(length(ps_trim_meta(trimmed)$trimmed_idx), 0)

  withr::local_seed(37)
  dat$z <- rbinom(nrow(dat), 1, 0.5)
  dat$k <- rpois(nrow(dat), 2)
  fits <- list(
    binomial = glm(z ~ x1 + x2, data = dat, family = binomial()),
    poisson = glm(k ~ x1 + x2, data = dat, family = poisson())
  )

  # The remedy names the dose models the density trim and the weights accept,
  # which include the robust and additive fits read by inheritance.
  for (kind in names(fits)) {
    cnd <- expect_error(
      ps_refit(trimmed, fits[[kind]]),
      class = "propensity_model_family_error",
      info = kind
    )
    message <- cli::ansi_strip(conditionMessage(cnd))
    expect_match(message, "mgcv::gam", fixed = TRUE, info = kind)
    expect_match(message, "MASS::rlm", fixed = TRUE, info = kind)
  }

  expect_propensity_error(ps_refit(trimmed, fits$binomial))
})

# ---- weights read from the record --------------------------------------------

test_that("weights from a refit density trim are the refit model's weights", {
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)
  t4 <- function(z) stats::dt(z, df = 4)

  densities <- list(
    normal = dens_normal(),
    t_mle = dens_t(4),
    t_rms = dens_t(4, sigma_method = "rms"),
    laplace = dens_laplace(),
    kernel = dens_kernel(),
    written = dens_fn(t4)
  )
  settings <- list(
    none = list(stabilize = FALSE, numerator = "marginal"),
    marginal = list(stabilize = TRUE, numerator = "marginal"),
    integrated = list(stabilize = TRUE, numerator = "integrated")
  )
  weight_fns <- list(wt_ate = wt_ate, wt_cens = wt_cens)

  for (family in names(densities)) {
    density <- densities[[family]]
    trimmed <- ps_trim(fit, method = "density", lower = 0.1, .density = density)
    meta <- ps_trim_meta(trimmed)
    keep <- meta$keep_idx
    expect_gt(length(meta$trimmed_idx), 0)

    refit <- ps_refit(trimmed, fit)
    refit_fit <- lm(a ~ x1 + x2, data = dat[keep, ])

    for (setting in names(settings)) {
      args <- settings[[setting]]
      for (fn_name in names(weight_fns)) {
        weight_fn <- weight_fns[[fn_name]]
        info <- paste(family, setting, fn_name)

        from_record <- weight_fn(
          refit,
          dat$a,
          exposure_type = "continuous",
          stabilize = args$stabilize,
          numerator = args$numerator
        )
        from_model <- weight_fn(
          refit_fit,
          exposure_type = "continuous",
          stabilize = args$stabilize,
          numerator = args$numerator,
          .density = density
        )

        expect_length(from_record, nrow(dat))
        expect_true(all(is.na(as.numeric(from_record)[-keep])), info = info)
        expect_equal(
          as.numeric(from_record)[keep],
          as.numeric(from_model),
          tolerance = 1e-12,
          info = info
        )

        recorded <- density_meta(from_record)
        expect_true(
          density_specs_agree(recorded$density, meta$density),
          info = info
        )
        expect_identical(recorded$sigma, meta$sigma_kind, info = info)
        expect_identical(
          recorded$numerator,
          if (args$stabilize) args$numerator else "none",
          info = info
        )
      }
    }
  }
})

test_that("weights from a refit density trim are the density ratio worked out by hand", {
  cases <- list(
    normal = list(
      density = dens_normal(),
      family = "normal",
      sigma_method = "rms"
    ),
    t_mle = list(
      density = dens_t(4),
      family = "t",
      sigma_method = "mle",
      df = 4
    ),
    laplace = list(
      density = dens_laplace(),
      family = "laplace",
      sigma_method = "mle"
    )
  )

  for (case in names(cases)) {
    spec <- cases[[case]]
    fx <- do.call(dose_refit_fixture, spec)
    expect_gt(length(fx$oracle$trimmed_idx), 0)
    refit <- ps_refit(fx$trimmed, fx$fit)
    sigma_kind <- if (identical(spec$sigma_method, "mle")) "mle" else "pooled"

    for (numerator in c("none", "marginal", "integrated")) {
      info <- paste(case, numerator)
      weights <- wt_ate(
        refit,
        fx$dat$a,
        exposure_type = "continuous",
        stabilize = numerator != "none",
        numerator = if (numerator == "none") "marginal" else numerator
      )
      by_hand <- dose_trim_weights(
        fx$dat$a,
        fx$refit_mu,
        sigma = fx$refit_sigma,
        numerator = numerator,
        family = spec$family,
        sigma_method = spec$sigma_method,
        df = spec$df
      )

      expect_equal(
        as.numeric(weights),
        by_hand,
        tolerance = 1e-8,
        info = info
      )

      recorded <- density_meta(weights)
      expect_identical(recorded$density$family, spec$family, info = info)
      expect_identical(recorded$sigma, sigma_kind, info = info)
      expect_equal(
        recorded$sigma_value,
        ps_trim_meta(refit)$sigma,
        tolerance = 1e-12,
        info = info
      )
    }
  }

  # The heavy-tailed record is the one a normal default would misread, so its
  # record is spelled out.
  fx <- do.call(dose_refit_fixture, cases$t_mle)
  weights <- wt_ate(
    ps_refit(fx$trimmed, fx$fit),
    fx$dat$a,
    exposure_type = "continuous"
  )
  recorded <- density_meta(weights)
  expect_identical(recorded$density$family, "t")
  expect_identical(recorded$density$params$df, 4)
  expect_identical(recorded$density$sigma_method, "mle")
  expect_identical(recorded$sigma, "mle")
  expect_identical(recorded$numerator, "marginal")
})

test_that("weights from a refit density trim reproduce the simulation's trimming arm", {
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)

  # The arm read each unit's conditional density as the inverse of its
  # unstabilized weight, kept the units at or above a quantile of it, refit the
  # dose model on them, and weighted by the refit model's stabilized weights.
  density <- 1 / as.numeric(wt_ate(fit, stabilize = FALSE))
  keep <- density >= stats::quantile(density, 0.05, names = FALSE)
  expect_gt(sum(!keep), 0)
  arm_fit <- lm(a ~ x1 + x2, data = dat[keep, ])
  arm_weights <- as.numeric(wt_ate(arm_fit))

  trimmed <- ps_trim(fit, method = "density", lower = 0.05)
  expect_identical(is_unit_trimmed(trimmed), !keep)

  weights <- wt_ate(
    ps_refit(trimmed, fit),
    dat$a,
    exposure_type = "continuous"
  )

  expect_equal(as.numeric(weights)[keep], arm_weights, tolerance = 1e-12)
  expect_true(all(is.na(as.numeric(weights)[!keep])))

  mu <- rep(NA_real_, nrow(dat))
  mu[keep] <- fitted(arm_fit)
  expect_equal(
    as.numeric(weights),
    dose_trim_weights(
      dat$a,
      mu,
      sigma = sqrt(mean(residuals(arm_fit)^2)),
      numerator = "marginal"
    ),
    tolerance = 1e-10
  )

  # The weights say what was done to them: trimmed, refit, and which units.
  expect_identical(estimand(weights), "ate; trimmed")
  expect_true(is_ps_trimmed(weights))
  expect_true(isTRUE(attr(weights, "trimmed")))
  expect_true(is_refit(weights))
  expect_true(is_stabilized(weights))
  expect_identical(is_unit_trimmed(weights), !keep)
  expect_identical(ps_trim_meta(weights)$method, "density")
  expect_identical(
    ps_trim_meta(weights)$sigma,
    ps_trim_meta(ps_refit(trimmed, fit))$sigma
  )
})

test_that("weights from a density trim that was not refit read the record's full-sample spread", {
  cases <- list(
    normal = list(
      density = dens_normal(),
      family = "normal",
      sigma_method = "rms"
    ),
    t_mle = list(
      density = dens_t(4),
      family = "t",
      sigma_method = "mle",
      df = 4
    )
  )

  for (case in names(cases)) {
    spec <- cases[[case]]
    fx <- do.call(dose_refit_fixture, spec)
    expect_gt(length(fx$oracle$trimmed_idx), 0)
    sigma_kind <- if (identical(spec$sigma_method, "mle")) "mle" else "pooled"

    for (numerator in c("none", "marginal", "integrated")) {
      info <- paste(case, numerator)
      expect_warning(
        weights <- wt_ate(
          fx$trimmed,
          fx$dat$a,
          exposure_type = "continuous",
          stabilize = numerator != "none",
          numerator = if (numerator == "none") "marginal" else numerator
        ),
        class = "propensity_no_refit_warning"
      )

      # The same object that made the cut makes the weights: the full-sample
      # means of the retained units, read at the full-sample spread.
      by_hand <- dose_trim_weights(
        fx$dat$a,
        fx$full_mu,
        sigma = fx$oracle$sigma,
        numerator = numerator,
        family = spec$family,
        sigma_method = spec$sigma_method,
        df = spec$df
      )
      expect_equal(
        as.numeric(weights),
        by_hand,
        tolerance = 1e-8,
        info = info
      )

      # Re-estimating the spread from the retained residuals alone gives other
      # weights, so the comparison above tells the two readings apart.
      reestimated <- wt_ate(
        fx$full_mu,
        fx$dat$a,
        exposure_type = "continuous",
        stabilize = numerator != "none",
        numerator = if (numerator == "none") "marginal" else numerator,
        .density = spec$density
      )
      expect_gt(
        max(abs(as.numeric(reestimated) - by_hand), na.rm = TRUE),
        1e-4
      )

      recorded <- density_meta(weights)
      expect_identical(recorded$density$family, spec$family, info = info)
      expect_identical(recorded$sigma, sigma_kind, info = info)
      expect_equal(
        recorded$sigma_value,
        fx$oracle$sigma,
        tolerance = 1e-8,
        info = info
      )

      expect_identical(estimand(weights), "ate; trimmed", info = info)
      expect_true(is_ps_trimmed(weights), info = info)
      expect_false(is_refit(weights), info = info)
      expect_identical(
        is_unit_trimmed(weights),
        seq_len(fx$n) %in% fx$oracle$trimmed_idx,
        info = info
      )
    }
  }
})

test_that("an unrefit density trim warns that it was not refit", {
  fx <- dose_refit_fixture()
  expect_gt(length(fx$oracle$trimmed_idx), 0)

  expect_propensity_warning(wt_ate(
    fx$trimmed,
    fx$dat$a,
    exposure_type = "continuous"
  ))
  expect_warning(
    wt_cens(fx$trimmed, fx$dat$a, exposure_type = "continuous"),
    class = "propensity_no_refit_warning"
  )
  expect_no_warning(wt_ate(
    ps_refit(fx$trimmed, fx$fit),
    fx$dat$a,
    exposure_type = "continuous"
  ))
})

test_that("censoring weights from a density trim carry the trimmed label", {
  fx <- dose_refit_fixture(
    density = dens_t(4),
    family = "t",
    sigma_method = "mle",
    df = 4
  )
  refit <- ps_refit(fx$trimmed, fx$fit)

  weights <- wt_cens(refit, fx$dat$a, exposure_type = "continuous")
  expect_identical(estimand(weights), "uncensored; trimmed")
  expect_true(is_ps_trimmed(weights))
  expect_true(is_refit(weights))
  expect_identical(
    is_unit_trimmed(weights),
    seq_len(fx$n) %in% fx$oracle$trimmed_idx
  )
  expect_identical(density_meta(weights)$density$family, "t")
  expect_identical(density_meta(weights)$sigma, "mle")
  expect_equal(
    as.numeric(weights),
    dose_trim_weights(
      fx$dat$a,
      fx$refit_mu,
      sigma = fx$refit_sigma,
      numerator = "marginal",
      family = "t",
      sigma_method = "mle",
      df = 4
    ),
    tolerance = 1e-8
  )
})

test_that("the integrated numerator reads a record's spread and still refuses a caller's", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)

  trimmed <- ps_trim(fit, method = "resid", upper = 2, .sigma = 1.5)
  oracle <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "resid",
    sigma = 1.5,
    upper = 2
  )
  expect_dose_trim(trimmed, oracle, n)
  refit <- ps_refit(trimmed, fit)
  keep <- oracle$keep_idx

  mu <- rep(NA_real_, n)
  mu[keep] <- fitted(lm(a ~ x1 + x2, data = dat[keep, ]))

  # A spread the caller fixed when trimming is the record's, so the integrated
  # numerator marginalizes over it rather than refusing it.
  for (fn in list(wt_ate, wt_cens)) {
    for (numerator in c("marginal", "integrated")) {
      weights <- fn(
        refit,
        dat$a,
        exposure_type = "continuous",
        numerator = numerator
      )
      expect_equal(
        as.numeric(weights),
        dose_trim_weights(dat$a, mu, sigma = 1.5, numerator = numerator),
        tolerance = 1e-8,
        info = numerator
      )
      expect_identical(density_meta(weights)$sigma, "supplied")
      expect_identical(density_meta(weights)$sigma_value, 1.5)
      expect_identical(density_meta(weights)$numerator, numerator)
    }
  }

  # The same spread handed to the numeric route by the caller is still refused
  # under the integrated numerator, so only a record opens it.
  expect_error(
    wt_ate(
      mu,
      dat$a,
      exposure_type = "continuous",
      .sigma = 1.5,
      numerator = "integrated"
    ),
    class = "propensity_numerator_error"
  )
  expect_error(
    wt_cens(
      mu,
      dat$a,
      exposure_type = "continuous",
      .sigma = 1.5,
      numerator = "integrated"
    ),
    class = "propensity_numerator_error"
  )
})

test_that("a density trim in a data frame is weighted as the trim on its own", {
  fx <- dose_refit_fixture(
    density = dens_t(4),
    family = "t",
    sigma_method = "mle",
    df = 4
  )
  refit <- ps_refit(fx$trimmed, fx$fit)
  by_hand <- dose_trim_weights(
    fx$dat$a,
    fx$refit_mu,
    sigma = fx$refit_sigma,
    numerator = "marginal",
    family = "t",
    sigma_method = "mle",
    df = 4
  )

  frames <- list(
    alone = data.frame(ps = refit),
    selected = data.frame(other = fx$refit_mu, ps = refit)
  )
  for (frame in names(frames)) {
    weights <- wt_ate(
      frames[[frame]],
      fx$dat$a,
      exposure_type = "continuous",
      .propensity_col = "ps"
    )
    expect_equal(as.numeric(weights), by_hand, tolerance = 1e-8, info = frame)
    expect_identical(estimand(weights), "ate; trimmed", info = frame)
    expect_identical(density_meta(weights)$density$family, "t", info = frame)
    expect_identical(density_meta(weights)$sigma, "mle", info = frame)
  }
})

# ---- weights from the record: refusals ---------------------------------------

test_that("weights from a density trim refuse a density the record does not hold", {
  normal <- dose_refit_fixture()
  normal_refit <- ps_refit(normal$trimmed, normal$fit)
  heavy <- dose_refit_fixture(
    density = dens_t(4),
    family = "t",
    sigma_method = "mle",
    df = 4
  )
  heavy_refit <- ps_refit(heavy$trimmed, heavy$fit)

  conflicts <- list(
    list(refit = normal_refit, density = "laplace"),
    list(refit = normal_refit, density = dens_t(4)),
    list(refit = heavy_refit, density = "normal"),
    list(refit = heavy_refit, density = dens_t(5)),
    list(refit = heavy_refit, density = dens_t(4, sigma_method = "rms"))
  )

  for (i in seq_along(conflicts)) {
    conflict <- conflicts[[i]]
    for (fn in list(wt_ate, wt_cens)) {
      expect_error(
        fn(
          conflict$refit,
          normal$dat$a,
          exposure_type = "continuous",
          .density = conflict$density
        ),
        class = "propensity_density_error",
        info = paste("conflict", i)
      )
    }
  }

  # A density that agrees with the record is the record's, and changes
  # nothing.
  expect_equal(
    wt_ate(
      normal_refit,
      normal$dat$a,
      exposure_type = "continuous",
      .density = "normal"
    ),
    wt_ate(normal_refit, normal$dat$a, exposure_type = "continuous")
  )
  agreeing <- wt_ate(
    heavy_refit,
    heavy$dat$a,
    exposure_type = "continuous",
    .density = dens_t(4)
  )
  expect_equal(
    as.numeric(agreeing),
    dose_trim_weights(
      heavy$dat$a,
      heavy$refit_mu,
      sigma = heavy$refit_sigma,
      numerator = "marginal",
      family = "t",
      sigma_method = "mle",
      df = 4
    ),
    tolerance = 1e-8
  )

  expect_propensity_error(wt_ate(
    normal_refit,
    normal$dat$a,
    exposure_type = "continuous",
    .density = "laplace"
  ))
})

test_that("weights from a density trim refuse a spread of the caller's", {
  estimated <- dose_refit_fixture()
  estimated_refit <- ps_refit(estimated$trimmed, estimated$fit)

  dat <- estimated$dat
  fit <- estimated$fit
  supplied <- ps_trim(fit, method = "resid", upper = 2, .sigma = 1.5)
  expect_gt(length(ps_trim_meta(supplied)$trimmed_idx), 0)
  supplied_refit <- ps_refit(supplied, fit)

  # The record holds the spread the trim was made at, so a second one is
  # refused even when it repeats the recorded number.
  for (fn in list(wt_ate, wt_cens)) {
    expect_error(
      fn(estimated_refit, dat$a, exposure_type = "continuous", .sigma = 1),
      class = "propensity_sigma_error"
    )
    expect_error(
      fn(supplied_refit, dat$a, exposure_type = "continuous", .sigma = 1.5),
      class = "propensity_sigma_error"
    )
    expect_error(
      fn(
        supplied_refit,
        dat$a,
        exposure_type = "continuous",
        .sigma = rep(1.5, nrow(dat))
      ),
      class = "propensity_sigma_error"
    )
  }

  expect_propensity_error(wt_ate(
    estimated_refit,
    dat$a,
    exposure_type = "continuous",
    .sigma = 1
  ))
})

test_that("weights from a density trim are refused for an exposure that is not a dose", {
  fx <- dose_refit_fixture()
  refit <- ps_refit(fx$trimmed, fx$fit)
  expect_gt(length(fx$oracle$trimmed_idx), 0)

  withr::local_seed(41)
  z <- rbinom(fx$n, 1, 0.5)

  for (fn in list(wt_ate, wt_cens)) {
    for (exposure_type in c("binary", "auto")) {
      expect_error(
        fn(refit, z, exposure_type = exposure_type),
        class = "propensity_modified_continuous_error",
        info = exposure_type
      )
      expect_error(
        fn(data.frame(ps = refit), z, exposure_type = exposure_type),
        class = "propensity_modified_continuous_error",
        info = exposure_type
      )
      expect_error(
        fn(
          data.frame(other = fx$refit_mu, ps = refit),
          z,
          exposure_type = exposure_type,
          .propensity_col = "ps"
        ),
        class = "propensity_modified_continuous_error",
        info = exposure_type
      )
    }
  }

  expect_propensity_error(wt_ate(refit, z, exposure_type = "binary"))
  expect_propensity_error(wt_ate(
    data.frame(ps = refit),
    z,
    exposure_type = "binary"
  ))
})

test_that("ipw() refuses weights built from a density trim", {
  fx <- dose_refit_fixture(
    density = dens_t(4),
    family = "t",
    sigma_method = "mle",
    df = 4
  )
  expect_gt(length(fx$oracle$trimmed_idx), 0)
  dat <- fx$dat
  withr::local_seed(43)
  dat$y <- 0.5 * dat$a + dat$x1 + rnorm(nrow(dat))

  weights <- wt_ate(
    ps_refit(fx$trimmed, fx$fit),
    dat$a,
    exposure_type = "continuous"
  )
  outcome_mod <- lm(y ~ a, data = dat, weights = weights)

  expect_error(
    ipw(fx$fit, outcome_mod, .data = dat),
    class = "propensity_ipw_trimmed_error"
  )
})
