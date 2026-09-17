# Trimming a dose model on the scale of its conditional density
#
# For a continuous exposure the generalized propensity score is the conditional
# density f(a | x), and a unit whose observed dose the model finds implausible
# has a small one. `ps_trim()` on a dose model sets those units aside, by a
# quantile floor on the density (`"density"`) or by a bound on the absolute
# standardized residual (`"resid"`), and keeps the conditional means of the rest.
# Every expectation below is held against a density worked out by hand in
# helper-dose-trim.R from the model's residuals, the family's density, and the
# family's own scale.

# ---- "density": a quantile floor --------------------------------------------

test_that("ps_trim() floors the conditional density of a lm fit at a quantile", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)

  trimmed <- ps_trim(fit, method = "density")
  oracle <- dose_trim_oracle(dat$a, fitted(fit), method = "density")
  expect_dose_trim(trimmed, oracle, n)

  # The default floor is the first percentile, which on sixty units sets one
  # aside.
  meta <- ps_trim_meta(trimmed)
  expect_identical(meta$method, "density")
  expect_identical(meta$lower, 0.01)
  expect_null(meta$upper)
  expect_length(meta$trimmed_idx, 1)

  # The spread is the pooled root mean square, the one a normal density is
  # read at, and the family is the default normal, resolved to a specification.
  expect_identical(meta$sigma_kind, "pooled")
  expect_s3_class(meta$density, "propensity_density")
  expect_identical(meta$density$family, "normal")

  # A higher floor sets more aside, and at the tenth percentile it is six.
  trimmed_10 <- ps_trim(fit, method = "density", lower = 0.1)
  oracle_10 <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "density",
    lower = 0.1
  )
  expect_dose_trim(trimmed_10, oracle_10, n)
  expect_length(ps_trim_meta(trimmed_10)$trimmed_idx, 6)
  expect_identical(ps_trim_meta(trimmed_10)$lower, 0.1)

  expect_true(is_ps_trimmed(trimmed))
  expect_identical(
    is_unit_trimmed(trimmed_10),
    seq_len(n) %in% oracle_10$trimmed_idx
  )
})

test_that("ps_trim() reads a t density at its own scale unless asked for the root mean square", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)

  mle <- ps_trim(fit, method = "density", lower = 0.05, .density = dens_t(4))
  oracle_mle <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "density",
    family = "t",
    df = 4,
    sigma_method = "mle",
    lower = 0.05
  )
  expect_dose_trim(mle, oracle_mle, n)
  expect_identical(ps_trim_meta(mle)$sigma_kind, "mle")
  expect_identical(ps_trim_meta(mle)$density$family, "t")
  expect_identical(ps_trim_meta(mle)$density$params$df, 4)
  expect_identical(ps_trim_meta(mle)$density$sigma_method, "mle")

  rms <- ps_trim(
    fit,
    method = "density",
    lower = 0.05,
    .density = dens_t(4, sigma_method = "rms")
  )
  oracle_rms <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "density",
    family = "t",
    df = 4,
    sigma_method = "rms",
    lower = 0.05
  )
  expect_dose_trim(rms, oracle_rms, n)
  expect_identical(ps_trim_meta(rms)$sigma_kind, "pooled")
  expect_identical(ps_trim_meta(rms)$density$sigma_method, "rms")

  # The two scales are different numbers, so the floors are too. A quantile
  # floor on a symmetric unimodal density ranks the units by their absolute
  # residual whatever the scale, so the two agree on which units go; the fixed
  # residual bound below is where the scale decides membership.
  expect_gt(abs(oracle_mle$sigma - oracle_rms$sigma), 0.1)
  expect_false(isTRUE(all.equal(
    ps_trim_meta(mle)$threshold,
    ps_trim_meta(rms)$threshold
  )))
})

test_that("ps_trim() reads a Laplace density at its mean absolute residual", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)

  trimmed <- ps_trim(
    fit,
    method = "density",
    lower = 0.05,
    .density = dens_laplace()
  )
  oracle <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "density",
    family = "laplace",
    sigma_method = "mle",
    lower = 0.05
  )
  expect_dose_trim(trimmed, oracle, n)
  expect_equal(ps_trim_meta(trimmed)$sigma, mean(abs(residuals(fit))))
  expect_identical(ps_trim_meta(trimmed)$sigma_kind, "mle")

  # The name of a family resolves to the specification the constructor builds,
  # which for the Laplace carries its own scale estimator.
  by_name <- ps_trim(
    fit,
    method = "density",
    lower = 0.05,
    .density = "laplace"
  )
  expect_dose_trim(by_name, oracle, n)
  expect_s3_class(ps_trim_meta(by_name)$density, "propensity_density")
  expect_identical(ps_trim_meta(by_name)$density$family, "laplace")
  expect_identical(ps_trim_meta(by_name)$density$sigma_method, "mle")
})

test_that("ps_trim() floors a kernel or user-written density on the same rule", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)

  # A kernel is fit to the standardized residuals, so the quantile floor
  # composes with a density no family describes.
  kernel <- ps_trim(
    fit,
    method = "density",
    lower = 0.1,
    .density = dens_kernel()
  )
  oracle_kernel <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "density",
    family = "kernel",
    lower = 0.1
  )
  expect_dose_trim(kernel, oracle_kernel, n)
  expect_identical(ps_trim_meta(kernel)$sigma_kind, "pooled")
  expect_identical(ps_trim_meta(kernel)$density$family, "kernel")

  t4 <- function(z) stats::dt(z, df = 4)
  written <- ps_trim(
    fit,
    method = "density",
    lower = 0.1,
    .density = dens_fn(t4)
  )
  oracle_written <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "density",
    family = "function",
    fn = t4,
    lower = 0.1
  )
  expect_dose_trim(written, oracle_written, n)
  expect_identical(ps_trim_meta(written)$sigma_kind, "pooled")
  expect_identical(ps_trim_meta(written)$density$family, "function")

  # A bare function is wrapped the way the weight functions wrap one.
  bare <- ps_trim(fit, method = "density", lower = 0.1, .density = t4)
  expect_dose_trim(bare, oracle_written, n)
  expect_identical(ps_trim_meta(bare)$density$family, "function")
})

test_that("ps_trim() trims a gaussian glm as it trims the lm it matches", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit_lm <- lm(a ~ x1 + x2, data = dat)
  fit_glm <- glm(a ~ x1 + x2, data = dat, family = gaussian())

  from_glm <- ps_trim(fit_glm, method = "density", lower = 0.05)
  oracle <- dose_trim_oracle(
    dat$a,
    fitted(fit_glm),
    method = "density",
    lower = 0.05
  )
  expect_dose_trim(from_glm, oracle, n)
  expect_identical(
    ps_trim_meta(from_glm)$keep_idx,
    ps_trim_meta(ps_trim(fit_lm, method = "density", lower = 0.05))$keep_idx
  )

  resid <- ps_trim(fit_glm, method = "resid", upper = 2)
  expect_dose_trim(
    resid,
    dose_trim_oracle(dat$a, fitted(fit_glm), method = "resid", upper = 2),
    n
  )

  # A family estimating its own scale reads it off the glm's residuals too.
  t_glm <- ps_trim(
    fit_glm,
    method = "density",
    lower = 0.05,
    .density = dens_t(4)
  )
  expect_dose_trim(
    t_glm,
    dose_trim_oracle(
      dat$a,
      fitted(fit_glm),
      method = "density",
      family = "t",
      df = 4,
      sigma_method = "mle",
      lower = 0.05
    ),
    n
  )
  expect_identical(ps_trim_meta(t_glm)$sigma_kind, "mle")
})

test_that("ps_trim() reads a gaussian glm under a log link on the dose's own scale", {
  dat <- sim_dose_trim()
  dat$a <- dat$a + 10
  n <- nrow(dat)
  fit <- glm(a ~ x1 + x2, data = dat, family = gaussian(link = "log"))

  # `fitted()` is on the scale of the dose whatever the link, so the residuals
  # are the dose less its fitted mean, not anything on the linear predictor.
  trimmed <- ps_trim(fit, method = "density", lower = 0.05)
  oracle <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "density",
    lower = 0.05
  )
  expect_dose_trim(trimmed, oracle, n)

  resid <- ps_trim(fit, method = "resid", upper = 2)
  expect_dose_trim(
    resid,
    dose_trim_oracle(dat$a, fitted(fit), method = "resid", upper = 2),
    n
  )
})

test_that("ps_trim() trims the dose models wt_ate() reads by inheritance", {
  dat <- sim_dose_trim()
  n <- nrow(dat)

  skip_if_not_installed("MASS")
  fit_rlm <- MASS::rlm(a ~ x1 + x2, data = dat)
  # The spread is the family's estimator on the residuals, not the robust
  # scale `rlm` reports, which a caller passes as `.sigma` if they want it.
  expect_dose_trim(
    ps_trim(fit_rlm, method = "density", lower = 0.05),
    dose_trim_oracle(dat$a, fitted(fit_rlm), method = "density", lower = 0.05),
    n
  )
  expect_dose_trim(
    ps_trim(fit_rlm, method = "resid", upper = 2, .density = dens_t(4)),
    dose_trim_oracle(
      dat$a,
      fitted(fit_rlm),
      method = "resid",
      family = "t",
      df = 4,
      sigma_method = "mle",
      upper = 2
    ),
    n
  )

  skip_if_not_installed("mgcv")
  fit_gam <- mgcv::gam(a ~ s(x1) + x2, data = dat, family = gaussian())
  expect_dose_trim(
    ps_trim(fit_gam, method = "density", lower = 0.05),
    dose_trim_oracle(dat$a, fitted(fit_gam), method = "density", lower = 0.05),
    n
  )
  expect_dose_trim(
    ps_trim(fit_gam, method = "resid", upper = 2, .density = dens_laplace()),
    dose_trim_oracle(
      dat$a,
      fitted(fit_gam),
      method = "resid",
      family = "laplace",
      sigma_method = "mle",
      upper = 2
    ),
    n
  )
})

# ---- "resid": a standardized-residual bound ----------------------------------

test_that("ps_trim() bounds the absolute standardized residual under each symmetric family", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)

  cases <- list(
    normal = list(
      density = dens_normal(),
      family = "normal",
      sigma_method = "rms",
      kind = "pooled"
    ),
    t_mle = list(
      density = dens_t(4),
      family = "t",
      sigma_method = "mle",
      kind = "mle"
    ),
    t_rms = list(
      density = dens_t(4, sigma_method = "rms"),
      family = "t",
      sigma_method = "rms",
      kind = "pooled"
    ),
    laplace_mle = list(
      density = dens_laplace(),
      family = "laplace",
      sigma_method = "mle",
      kind = "mle"
    ),
    laplace_rms = list(
      density = dens_laplace(sigma_method = "rms"),
      family = "laplace",
      sigma_method = "rms",
      kind = "pooled"
    )
  )

  kept <- list()
  for (case in names(cases)) {
    spec <- cases[[case]]
    trimmed <- ps_trim(
      fit,
      method = "resid",
      upper = 2,
      .density = spec$density
    )
    oracle <- dose_trim_oracle(
      dat$a,
      fitted(fit),
      method = "resid",
      family = spec$family,
      df = 4,
      sigma_method = spec$sigma_method,
      upper = 2
    )
    expect_dose_trim(trimmed, oracle, n)

    meta <- ps_trim_meta(trimmed)
    expect_identical(meta$method, "resid", info = case)
    expect_identical(meta$upper, 2, info = case)
    expect_null(meta$lower)
    expect_identical(meta$sigma_kind, spec$kind, info = case)

    # The bound is a density floor written in residual units: the record keeps
    # that floor, g(k) / sigma, and it keeps exactly the units whose density
    # reaches it.
    expect_equal(
      meta$threshold,
      oracle$threshold,
      tolerance = 1e-12,
      info = case
    )
    expect_identical(
      as.integer(meta$keep_idx),
      which(oracle$f >= meta$threshold),
      info = case
    )

    kept[[case]] <- as.integer(meta$keep_idx)
  }

  # At a fixed bound the scale decides membership: the t and Laplace scales a
  # family estimates for itself are smaller than the root mean square, so the
  # same bound reaches a unit the root mean square leaves in.
  expect_gt(length(setdiff(kept$t_rms, kept$t_mle)), 0)
  expect_gt(length(setdiff(kept$laplace_rms, kept$laplace_mle)), 0)
  expect_identical(kept$normal, kept$t_rms)
})

test_that("ps_trim() reads a supplied spread and records it as supplied", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)

  resid <- ps_trim(fit, method = "resid", upper = 2, .sigma = 1.5)
  oracle_resid <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "resid",
    sigma = 1.5,
    upper = 2
  )
  expect_dose_trim(resid, oracle_resid, n)
  expect_identical(ps_trim_meta(resid)$sigma, 1.5)
  expect_identical(ps_trim_meta(resid)$sigma_kind, "supplied")

  # A family asked for the root mean square takes a spread of the caller's in
  # its place; only a family estimating its own scale refuses one.
  density <- ps_trim(
    fit,
    method = "density",
    lower = 0.05,
    .sigma = 0.8,
    .density = dens_t(4, sigma_method = "rms")
  )
  oracle_density <- dose_trim_oracle(
    dat$a,
    fitted(fit),
    method = "density",
    family = "t",
    df = 4,
    sigma = 0.8,
    lower = 0.05
  )
  expect_dose_trim(density, oracle_density, n)
  expect_identical(ps_trim_meta(density)$sigma_kind, "supplied")
})

# ---- missing values ----------------------------------------------------------

test_that("ps_trim() leaves a unit with a missing mean or dose out of a density trim", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  missing_x <- c(4L, 30L)
  dat$x1[missing_x] <- NA
  fit <- lm(a ~ x1 + x2, data = dat, na.action = na.exclude)

  # A covariate the model dropped leaves a missing mean, padded back to the
  # length of the data.
  mu <- fitted(fit)
  expect_true(all(is.na(mu[missing_x])))

  trimmed <- ps_trim(fit, method = "density", lower = 0.05)
  oracle <- dose_trim_oracle(dat$a, mu, method = "density", lower = 0.05)
  expect_dose_trim(trimmed, oracle, n)
  meta <- ps_trim_meta(trimmed)
  expect_false(any(missing_x %in% meta$keep_idx))
  expect_false(any(missing_x %in% meta$trimmed_idx))
  expect_false(any(is_unit_trimmed(trimmed)[missing_x]))

  # The floor is the one the complete units alone give.
  complete <- setdiff(seq_len(n), missing_x)
  oracle_complete <- dose_trim_oracle(
    dat$a[complete],
    mu[complete],
    method = "density",
    lower = 0.05
  )
  expect_equal(meta$threshold, oracle_complete$threshold, tolerance = 1e-12)

  # A dose missing from the exposure supplied beside a complete fit leaves a
  # missing residual in the same way.
  fit_full <- lm(a ~ x2, data = sim_dose_trim())
  dose <- sim_dose_trim()$a
  missing_a <- c(9L, 51L)
  dose[missing_a] <- NA
  resid <- ps_trim(fit_full, method = "resid", upper = 1.5, .exposure = dose)
  oracle_resid <- dose_trim_oracle(
    dose,
    fitted(fit_full),
    method = "resid",
    upper = 1.5
  )
  expect_dose_trim(resid, oracle_resid, n)
  meta_resid <- ps_trim_meta(resid)
  expect_false(any(missing_a %in% meta_resid$keep_idx))
  expect_false(any(missing_a %in% meta_resid$trimmed_idx))
})

# ---- the record --------------------------------------------------------------

# The record names all ten fields whichever method made it, holding `NULL` for
# the bound the method does not read, as `list()` keeps a `NULL` entry by name.
test_that("the density trim record carries the fields that describe the cut", {
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)

  for (method in c("density", "resid")) {
    trimmed <- if (method == "density") {
      ps_trim(fit, method = method)
    } else {
      ps_trim(fit, method = method, upper = 2)
    }
    meta <- ps_trim_meta(trimmed)
    expect_true(
      all(
        c(
          "method",
          "lower",
          "upper",
          "threshold",
          "sigma",
          "sigma_kind",
          "density",
          "keep_idx",
          "trimmed_idx",
          "n_obs"
        ) %in%
          names(meta)
      ),
      info = method
    )
    expect_true(is.numeric(meta$threshold) && length(meta$threshold) == 1)
    expect_true(meta$threshold > 0)
    expect_true(is.numeric(meta$sigma) && length(meta$sigma) == 1)
    expect_null(meta$refit)
  }
})

test_that("trim_parameters() tells apart density trims cut at different places", {
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)
  meta <- ps_trim_meta(
    ps_trim(fit, method = "resid", upper = 2, .density = dens_t(4))
  )
  params <- trim_parameters(meta)

  for (field in c("threshold", "sigma", "sigma_kind", "density")) {
    expect_true(field %in% names(params), info = field)
  }

  differs <- function(changed) {
    !identical(params, trim_parameters(changed))
  }

  other_spread <- meta
  other_spread$sigma <- meta$sigma * 1.1
  expect_true(differs(other_spread))

  other_threshold <- meta
  other_threshold$threshold <- meta$threshold * 1.1
  expect_true(differs(other_threshold))

  other_kind <- meta
  other_kind$sigma_kind <- "supplied"
  expect_true(differs(other_kind))

  other_family <- meta
  other_family$density <- dens_laplace()
  expect_true(differs(other_family))

  other_df <- meta
  other_df$density <- dens_t(5)
  expect_true(differs(other_df))

  other_method <- meta
  other_method$density <- dens_t(4, sigma_method = "rms")
  expect_true(differs(other_method))

  # Each constructor call writes a fresh function, so a specification built
  # again is compared on what identifies it rather than on its closure.
  same_family <- meta
  same_family$density <- dens_t(4)
  expect_false(differs(same_family))

  # A density written by the caller is identified by the function itself.
  meta_fn <- ps_trim_meta(
    ps_trim(
      fit,
      method = "density",
      .density = dens_fn(function(z) stats::dt(z, df = 4))
    )
  )
  other_fn <- meta_fn
  other_fn$density <- dens_fn(function(z) stats::dt(z, df = 6))
  expect_false(identical(trim_parameters(meta_fn), trim_parameters(other_fn)))
})

test_that("c() keeps a density trim only when the two cuts agree", {
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)

  base <- ps_trim(fit, method = "resid", upper = 2, .density = dens_t(4))
  again <- ps_trim(fit, method = "resid", upper = 2, .density = dens_t(4))

  combined <- c(base, again)
  expect_s3_class(combined, "ps_trim")
  meta <- ps_trim_meta(combined)
  expect_identical(meta$method, "resid")
  expect_identical(meta$sigma, ps_trim_meta(base)$sigma)
  expect_identical(meta$threshold, ps_trim_meta(base)$threshold)
  expect_identical(meta$sigma_kind, "mle")
  expect_identical(meta$density$family, "t")

  # Trims of the same model under the same bound that differ in how the spread
  # was arrived at, in the spread itself, or in the family are different cuts.
  others <- list(
    sigma_method = ps_trim(
      fit,
      method = "resid",
      upper = 2,
      .density = dens_t(4, sigma_method = "rms")
    ),
    spread = ps_trim(
      fit,
      method = "resid",
      upper = 2,
      .sigma = 1,
      .density = dens_t(4, sigma_method = "rms")
    ),
    family = ps_trim(
      fit,
      method = "resid",
      upper = 2,
      .density = dens_laplace()
    )
  )

  for (kind in names(others)) {
    expect_warning(
      mixed <- c(base, others[[kind]]),
      class = "propensity_coercion_warning"
    )
    expect_false(inherits(mixed, "ps_trim"), info = kind)
    expect_type(mixed, "double")
  }

  # Two floors at the same quantile under different spreads trim the same
  # units, and still describe different cuts.
  spread_a <- ps_trim(fit, method = "density", .sigma = 1)
  spread_b <- ps_trim(fit, method = "density", .sigma = 2)
  expect_identical(
    ps_trim_meta(spread_a)$keep_idx,
    ps_trim_meta(spread_b)$keep_idx
  )
  expect_warning(
    mixed <- c(spread_a, spread_b),
    class = "propensity_coercion_warning"
  )
  expect_false(inherits(mixed, "ps_trim"))
})

test_that("a density trim keeps its description through subsetting and loses its class to arithmetic", {
  dat <- sim_dose_trim()
  n <- nrow(dat)
  fit <- lm(a ~ x1 + x2, data = dat)
  trimmed <- ps_trim(fit, method = "density", lower = 0.1, .density = dens_t(4))
  meta <- ps_trim_meta(trimmed)
  trimmed_units <- is_unit_trimmed(trimmed)
  expect_identical(which(trimmed_units), as.integer(meta$trimmed_idx))

  sub <- trimmed[10:n]
  sub_meta <- ps_trim_meta(sub)
  expect_s3_class(sub, "ps_trim")
  expect_identical(is_unit_trimmed(sub), trimmed_units[10:n])
  expect_identical(sub_meta$n_obs, length(10:n))
  for (field in c("method", "lower", "threshold", "sigma", "sigma_kind")) {
    expect_identical(sub_meta[[field]], meta[[field]], info = field)
  }
  expect_identical(sub_meta$density$family, "t")
  expect_identical(sub_meta$density$params, meta$density$params)
  expect_identical(sub_meta$density$sigma_method, meta$density$sigma_method)

  # A transformed conditional mean is not a dose model.
  inverse <- 1 / trimmed
  expect_false(inherits(inverse, "ps_trim"))
  expect_type(inverse, "double")
  expect_identical(which(is.na(inverse)), which(trimmed_units))
  expect_false(inherits(trimmed + 1, "ps_trim"))

  expect_identical(
    vctrs::vec_ptype_full(trimmed),
    paste0("ps_trim; trimmed ", length(meta$trimmed_idx), " of ", n)
  )
  expect_output(print(trimmed), "ps_trim; trimmed 6 of 60", fixed = TRUE)
})

# ---- refusals ----------------------------------------------------------------

test_that("ps_trim() refuses an upper bound on a density floor", {
  fit <- lm(a ~ x1 + x2, data = sim_dose_trim())

  # A large conditional density is a plausible dose, so the floor has one side.
  expect_error(
    ps_trim(fit, method = "density", upper = 0.99),
    class = "propensity_unsupported_arg_error"
  )
  expect_propensity_error(ps_trim(fit, method = "density", upper = 0.99))
})

test_that("ps_trim() refuses a lower bound on a residual bound, and needs an upper one", {
  fit <- lm(a ~ x1 + x2, data = sim_dose_trim())

  # The bound is on the absolute residual, so it has no lower end.
  expect_error(
    ps_trim(fit, method = "resid", lower = 1, upper = 3),
    class = "propensity_unsupported_arg_error"
  )
  expect_propensity_error(ps_trim(fit, method = "resid", lower = 1, upper = 3))

  # There is no default bound: it is the analyst's to state.
  expect_error(
    ps_trim(fit, method = "resid"),
    class = "propensity_missing_arg_error"
  )
  expect_propensity_error(ps_trim(fit, method = "resid"))
})

test_that("ps_trim() refuses a density floor outside (0, 0.5) and a residual bound that is not positive", {
  fit <- lm(a ~ x1 + x2, data = sim_dose_trim())

  for (lower in list(0, 0.5, 0.7, -0.1, NA_real_, c(0.01, 0.02))) {
    expect_error(
      ps_trim(fit, method = "density", lower = lower),
      class = "propensity_range_error",
      info = deparse(lower)
    )
  }
  expect_propensity_error(ps_trim(fit, method = "density", lower = 0.5))

  for (upper in list(0, -2, NA_real_, Inf, c(2, 3))) {
    expect_error(
      ps_trim(fit, method = "resid", upper = upper),
      class = "propensity_range_error",
      info = deparse(upper)
    )
  }
  expect_propensity_error(ps_trim(fit, method = "resid", upper = -2))
})

test_that("ps_trim() refuses a residual bound under a density that need not be symmetric", {
  fit <- lm(a ~ x1 + x2, data = sim_dose_trim())

  # Only for a symmetric unimodal family is a residual bound a density floor,
  # so a kernel or a function of the caller's is pointed at "density".
  densities <- list(
    kernel = dens_kernel(),
    kernel_name = "kernel",
    written = dens_fn(function(z) stats::dt(z, df = 4))
  )
  for (kind in names(densities)) {
    cnd <- expect_error(
      ps_trim(fit, method = "resid", upper = 2, .density = densities[[kind]]),
      class = "propensity_density_error",
      info = kind
    )
    message <- gsub("[[:space:]]+", " ", cli::ansi_strip(conditionMessage(cnd)))
    expect_match(message, "\"density\"", fixed = TRUE, info = kind)
  }

  expect_propensity_error(
    ps_trim(fit, method = "resid", upper = 2, .density = dens_kernel())
  )
})

test_that("ps_trim() refuses a spread it cannot record on a dose model", {
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)

  # One spread is a constant a refit can keep; one per unit is not.
  expect_error(
    ps_trim(fit, method = "density", .sigma = rep(1, nrow(dat))),
    class = "propensity_sigma_error"
  )
  expect_propensity_error(
    ps_trim(fit, method = "density", .sigma = rep(1, nrow(dat)))
  )

  # A spread of the caller's and a family estimating its own scale are two
  # instructions about the same quantity, refused as `wt_ate()` refuses them.
  expect_error(
    ps_trim(fit, method = "resid", upper = 2, .sigma = 1, .density = dens_t(4)),
    class = "propensity_density_error"
  )
  expect_propensity_error(
    ps_trim(fit, method = "resid", upper = 2, .sigma = 1, .density = dens_t(4))
  )
})

test_that("ps_trim() refuses a dose that does not line up with the fitted means", {
  dat <- sim_dose_trim()
  n <- nrow(dat)

  # A fit that omits the rows it dropped reports fewer means than the data has
  # doses, so the full dose cannot be read against them.
  dat_missing <- dat
  dat_missing$x1[c(4L, 30L)] <- NA
  fit_omit <- lm(a ~ x1 + x2, data = dat_missing, na.action = na.omit)
  expect_length(fitted(fit_omit), n - 2L)
  expect_error(
    ps_trim(fit_omit, method = "density", .exposure = dat_missing$a),
    class = "propensity_length_error"
  )
  expect_propensity_error(
    ps_trim(fit_omit, method = "density", .exposure = dat_missing$a)
  )

  # The dose is read as a number for each unit, so neither text nor a matrix
  # will do.
  fit <- lm(a ~ x1 + x2, data = dat)
  expect_error(
    ps_trim(fit, method = "density", .exposure = as.character(dat$a)),
    class = "propensity_type_error"
  )
  expect_error(
    ps_trim(fit, method = "resid", upper = 2, .exposure = matrix(dat$a)),
    class = "propensity_type_error"
  )
  expect_propensity_error(
    ps_trim(fit, method = "density", .exposure = as.character(dat$a))
  )
})

test_that("ps_trim() refuses a dose model of several responses", {
  dat <- sim_dose_trim()
  dat$b <- dat$a + dat$x1
  fit <- lm(cbind(a, b) ~ x1 + x2, data = dat)

  expect_error(
    ps_trim(fit, method = "density"),
    class = "propensity_ps_shape_error"
  )
  expect_propensity_error(ps_trim(fit, method = "density"))
})

test_that("ps_trim() refuses a dose model fit with a spread that moves with its mean", {
  withr::local_seed(5)
  dat <- sim_dose_trim()
  dat$count <- rpois(nrow(dat), lambda = exp(0.3 * dat$x1 + 1))
  fits <- list(
    poisson = glm(count ~ x1, data = dat, family = poisson()),
    quasi_mu = glm(
      count ~ x1,
      data = dat,
      family = quasi(link = "log", variance = "mu")
    )
  )
  fit <- fits$poisson

  for (kind in names(fits)) {
    for (method in c("density", "resid")) {
      expect_error(
        ps_trim(
          fits[[kind]],
          method = method,
          upper = if (method == "resid") 2
        ),
        class = "propensity_model_family_error",
        info = paste(kind, method)
      )
    }
  }
  expect_propensity_error(ps_trim(fit, method = "density"))
})

test_that("ps_trim() refuses a score method on a dose model and names the dose methods", {
  dat <- sim_dose_trim()
  fits <- list(
    lm = lm(a ~ x1 + x2, data = dat),
    gaussian = glm(a ~ x1 + x2, data = dat, family = gaussian())
  )

  for (kind in names(fits)) {
    for (method in c("ps", "adaptive", "pctl", "pref", "cr", "optimal")) {
      cnd <- expect_error(
        ps_trim(fits[[kind]], method = method),
        class = "propensity_method_error",
        info = paste(kind, method)
      )
      message <- cli::ansi_strip(conditionMessage(cnd))
      expect_match(message, "\"density\"", fixed = TRUE, info = kind)
      expect_match(message, "\"resid\"", fixed = TRUE, info = kind)
    }

    # Left unnamed, the method is the generic's first choice, a score method, so
    # a dose model is refused rather than trimmed on a scale nobody asked for.
    implicit <- expect_error(
      ps_trim(fits[[kind]]),
      class = "propensity_method_error",
      info = kind
    )
    implicit_message <- cli::ansi_strip(conditionMessage(implicit))
    expect_match(implicit_message, "\"density\"", fixed = TRUE, info = kind)
    expect_match(implicit_message, "\"resid\"", fixed = TRUE, info = kind)
    explicit <- rlang::catch_cnd(ps_trim(fits[[kind]], method = "ps"))
    expect_identical(
      implicit_message,
      cli::ansi_strip(conditionMessage(explicit)),
      info = kind
    )
  }

  expect_propensity_error(ps_trim(fits$lm, method = "ps"))
  expect_propensity_error(ps_trim(fits$gaussian, method = "adaptive"))
})

test_that("ps_trim() refuses levels on a dose model", {
  fit <- lm(a ~ x1 + x2, data = sim_dose_trim())

  # A dose has no levels to name.
  expect_error(
    ps_trim(fit, method = "density", .focal_level = 1),
    class = "propensity_focal_level_error"
  )
  expect_error(
    ps_trim(fit, method = "resid", upper = 2, .reference_level = 0),
    class = "propensity_focal_level_error"
  )
  expect_propensity_error(ps_trim(fit, method = "density", .focal_level = 1))
})

test_that("ps_trim() refuses the deprecated level arguments on a dose model", {
  dat <- sim_dose_trim()
  fit <- lm(a ~ x1 + x2, data = dat)
  fit_glm <- glm(a ~ x1 + x2, data = dat, family = gaussian())

  # The deprecated pair names levels too, and is refused the same way once its
  # deprecation has been reported. A gaussian glm reaches the dose route by its
  # own method, so it is held to the same refusal.
  deprecated <- list(
    treated = quote(ps_trim(fit, method = "density", .treated = 1)),
    untreated = quote(
      ps_trim(fit, method = "resid", upper = 2, .untreated = 0)
    ),
    glm_treated = quote(ps_trim(fit_glm, method = "density", .treated = 1))
  )
  for (kind in names(deprecated)) {
    warned <- 0L
    with_always_deprecated({
      expect_error(
        withCallingHandlers(
          eval(deprecated[[kind]]),
          lifecycle_warning_deprecated = function(cnd) {
            warned <<- warned + 1L
            invokeRestart("muffleWarning")
          }
        ),
        class = "propensity_focal_level_error",
        info = kind
      )
    })
    expect_identical(warned, 1L, info = kind)
  }
})

test_that("ps_trim() refuses a density or a spread on the score routes", {
  withr::local_seed(8)
  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  fit <- glm(z ~ x, family = binomial())
  ps <- as.numeric(fitted(fit))

  # The score routes read no density, so a family or a spread given to them
  # would be an instruction silently ignored.
  expect_error(
    ps_trim(ps, method = "ps", .density = dens_t(4)),
    class = "propensity_density_error"
  )
  expect_error(
    ps_trim(ps, method = "pctl", .density = "laplace"),
    class = "propensity_density_error"
  )
  expect_error(
    ps_trim(fit, method = "ps", .density = dens_t(4)),
    class = "propensity_density_error"
  )
  expect_error(
    ps_trim(ps, method = "ps", .sigma = 1),
    class = "propensity_sigma_error"
  )
  expect_error(
    ps_trim(fit, method = "adaptive", .sigma = 1),
    class = "propensity_sigma_error"
  )

  # A matrix of scores for a categorical exposure, and a data frame of them,
  # are score routes too.
  trt <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  ps_matrix <- matrix(
    c(rep(0.3, n), rep(0.3, n), rep(0.4, n)),
    ncol = 3,
    dimnames = list(NULL, levels(trt))
  )
  expect_error(
    ps_trim(ps_matrix, method = "ps", .exposure = trt, .sigma = 1),
    class = "propensity_sigma_error"
  )
  expect_error(
    ps_trim(ps_matrix, method = "ps", .exposure = trt, .density = "laplace"),
    class = "propensity_density_error"
  )
  expect_error(
    ps_trim(
      as.data.frame(ps_matrix),
      method = "ps",
      .exposure = trt,
      .sigma = 1
    ),
    class = "propensity_sigma_error"
  )
  expect_error(
    ps_trim(
      as.data.frame(ps_matrix),
      method = "ps",
      .exposure = trt,
      .density = "laplace"
    ),
    class = "propensity_density_error"
  )
  expect_error(
    ps_trim(data.frame(ps = ps), method = "ps", .density = dens_t(4)),
    class = "propensity_density_error"
  )

  # The default family is no instruction at all.
  expect_s3_class(ps_trim(ps, method = "ps", .density = "normal"), "ps_trim")
  expect_s3_class(
    ps_trim(ps_matrix, method = "ps", .exposure = trt, .density = "normal"),
    "ps_trim"
  )

  expect_propensity_error(ps_trim(ps, method = "ps", .density = dens_t(4)))
  expect_propensity_error(ps_trim(ps, method = "ps", .sigma = 1))
})

test_that("ps_trim() refuses the density methods on a model of a probability", {
  withr::local_seed(9)
  n <- 90
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  score_data <- data.frame(z = z, x = x)

  fits <- list(
    binomial = glm(z ~ x, family = binomial(), data = score_data),
    quasibinomial = glm(z ~ x, family = quasibinomial(), data = score_data)
  )

  if (rlang::is_installed("nnet")) {
    score_data$z2 <- factor(z)
    score_data$z3 <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
    fits$multinom_two <- nnet::multinom(
      z2 ~ x,
      data = score_data,
      trace = FALSE
    )
    fits$multinom_three <- nnet::multinom(
      z3 ~ x,
      data = score_data,
      trace = FALSE
    )
  }

  # A model of a probability has no conditional density to trim, and the
  # refusal says so rather than listing the methods a score takes.
  for (kind in names(fits)) {
    for (method in c("density", "resid")) {
      cnd <- expect_error(
        ps_trim(
          fits[[kind]],
          method = method,
          upper = if (method == "resid") 2
        ),
        class = "propensity_method_error",
        info = paste(kind, method)
      )
      message <- cli::ansi_strip(conditionMessage(cnd))
      expect_no_match(message, "must be one of", fixed = TRUE)
    }
  }

  expect_propensity_error(ps_trim(fits$binomial, method = "density"))
  expect_propensity_error(ps_trim(fits$binomial, method = "resid", upper = 2))

  skip_if_not_installed("nnet")
  expect_propensity_error(ps_trim(fits$multinom_two, method = "density"))
  expect_propensity_error(ps_trim(
    fits$multinom_three,
    method = "resid",
    upper = 2
  ))
})
