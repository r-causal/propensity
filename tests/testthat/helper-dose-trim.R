# Fixtures and oracles for trimming a dose model on the scale of its conditional
# density.
#
# The dose is drawn with t(3) residuals so that a few units sit far enough from
# their predicted dose for a trim to find, and so that the scale estimators the
# density families take disagree by enough to move a fixed residual bound.
sim_dose_trim <- function(seed = 1234, n = 60) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  a <- 1 + 0.6 * x1 - 0.4 * x2 + stats::rt(n, df = 3)
  data.frame(a = a, x1 = x1, x2 = x2)
}

# The spread a density family reads its residuals at, computed here without the
# package's estimators. The root mean square is uncentered, because the
# residuals of a fitted model already average to zero. The Laplace scale is the
# mean absolute residual. The t scale is the root of the estimating equation
# `mean((df + 1) r^2 / (df s^2 + r^2)) = 1`, bracketed between a scale far
# below the root mean square, where the left side tends to `df + 1`, and one
# far above it, where it tends to zero.
dose_trim_sigma <- function(
  residuals,
  family = c("normal", "t", "laplace", "kernel", "function"),
  sigma_method = "rms",
  df = NULL
) {
  family <- match.arg(family)
  r <- residuals[!is.na(residuals)]
  rms <- sqrt(mean(r^2))

  if (!identical(sigma_method, "mle")) {
    return(rms)
  }

  switch(
    family,
    laplace = mean(abs(r)),
    t = stats::uniroot(
      function(s) mean((df + 1) * r^2 / (df * s^2 + r^2)) - 1,
      interval = rms * c(1e-3, 10),
      tol = 1e-14,
      maxiter = 1000
    )$root,
    stop("No maximum likelihood scale for this family.")
  )
}

# The standardized density `g` of each family, written out.
dose_trim_g <- function(z, family, df = NULL, fn = NULL) {
  switch(
    family,
    normal = stats::dnorm(z),
    t = stats::dt(z, df = df),
    laplace = exp(-abs(z)) / 2,
    "function" = fn(z)
  )
}

# A Gaussian kernel estimate of the standardized residuals that are present, fit
# over their range at the `stats::density()` defaults and interpolated to each.
dose_trim_kernel <- function(z) {
  present <- !is.na(z)
  fit <- stats::density(
    z[present],
    bw = "nrd0",
    adjust = 1,
    kernel = "gaussian",
    n = 512,
    from = min(z[present]),
    to = max(z[present])
  )
  out <- rep(NA_real_, length(z))
  out[present] <- stats::approxfun(fit$x, fit$y)(z[present])
  out
}

# What a density trim of the dose `a` around the conditional means `mu` keeps,
# worked out by hand. `"density"` keeps the units whose conditional density
# `g(z) / sigma` is at least its `lower` quantile; `"resid"` keeps those whose
# absolute standardized residual is at most `upper`, and records the density
# floor `g(upper) / sigma` that bound is equivalent to. A unit with a missing
# mean or dose takes no part in either.
dose_trim_oracle <- function(
  a,
  mu,
  method = c("density", "resid"),
  family = "normal",
  sigma_method = "rms",
  df = NULL,
  fn = NULL,
  sigma = NULL,
  lower = 0.01,
  upper = NULL
) {
  method <- match.arg(method)
  a <- as.numeric(a)
  mu <- as.numeric(mu)
  r <- a - mu
  if (is.null(sigma)) {
    sigma <- dose_trim_sigma(r, family, sigma_method, df)
  }
  z <- r / sigma

  g <- if (identical(family, "kernel")) {
    dose_trim_kernel(z)
  } else {
    dose_trim_g(z, family, df = df, fn = fn)
  }
  f <- g / sigma
  observed <- !is.na(f)

  if (method == "density") {
    threshold <- unname(stats::quantile(f, probs = lower, na.rm = TRUE))
    keep_idx <- which(f >= threshold)
  } else {
    threshold <- dose_trim_g(upper, family, df = df, fn = fn) / sigma
    keep_idx <- which(abs(z) <= upper)
  }

  list(
    mu = mu,
    sigma = sigma,
    z = z,
    f = f,
    threshold = threshold,
    keep_idx = keep_idx,
    trimmed_idx = setdiff(which(observed), keep_idx)
  )
}

# Holds a density trim against its oracle: the retained and trimmed positions,
# the retained values, which are the conditional means, and the realized floor.
expect_dose_trim <- function(trimmed, oracle, n) {
  meta <- ps_trim_meta(trimmed)

  testthat::expect_s3_class(trimmed, "ps_trim")
  testthat::expect_identical(as.integer(meta$keep_idx), oracle$keep_idx)
  testthat::expect_identical(
    as.integer(meta$trimmed_idx),
    as.integer(oracle$trimmed_idx)
  )
  testthat::expect_identical(meta$n_obs, n)
  testthat::expect_equal(meta$sigma, oracle$sigma, tolerance = 1e-8)
  testthat::expect_equal(meta$threshold, oracle$threshold, tolerance = 1e-8)

  values <- as.numeric(trimmed)
  testthat::expect_true(all(is.na(values[oracle$trimmed_idx])))
  testthat::expect_equal(
    values[oracle$keep_idx],
    oracle$mu[oracle$keep_idx],
    tolerance = 1e-12
  )

  # A fixture that trims nothing, or everything, would pass the comparisons
  # above without exercising the rule.
  testthat::expect_gt(length(oracle$trimmed_idx), 0)
  testthat::expect_gt(length(oracle$keep_idx), 0)
}
