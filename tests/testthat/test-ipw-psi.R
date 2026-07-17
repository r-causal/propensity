# Tests for the internal psi-building layer that backs the M-estimation path of
# ipw(): the weight registry (ipw_weight_fn), the theta layout (ipw_theta_layout),
# and the stacked estimating-function builder (build_ipw_psi). All three are
# unexported and reached through propensity:::. Specs are constructed directly as
# plain lists following the M-estimation design contract.

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
  data.frame(x1, x2, a, y)
}

sim_continuous <- function(seed = 77, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 0.5 + 0.7 * x1 - 0.3 * x2 + rnorm(n)
  y <- 1 + 0.6 * A + 0.4 * x1 + rnorm(n)
  data.frame(x1, x2, A, y)
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
  stab_score = NULL
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

  out_mod <- glm(
    y ~ a + x1,
    data = dat,
    family = quasibinomial(),
    weights = as.double(wts)
  )
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
      y = dat$y,
      family = "binomial",
      link = "logit",
      coefs = coef(out_mod),
      X_counterfactual = x_cf
    ),
    contrasts = c("rd", "log(rr)", "log(or)"),
    focal_level = focal_level,
    reference_level = lev[1]
  )
}

continuous_spec <- function(dat, stabilize = TRUE, stab_score = NULL) {
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  fitted_ps <- as.double(fitted(ps_mod))
  A <- dat$A
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      fitted_ps,
      A,
      exposure_type = "continuous",
      stabilize = stabilize,
      stabilization_score = stab_score
    )
  )
  msm <- lm(y ~ A, data = dat, weights = as.double(wts))

  list(
    exposure_type = "continuous",
    estimand = "ate",
    n = nrow(dat),
    exposure = A,
    ps = list(
      X = model.matrix(ps_mod),
      link = NULL,
      coefs = coef(ps_mod),
      k = NULL
    ),
    stab = list(stabilized = isTRUE(stabilize), score = stab_score),
    outcome = list(
      X = model.matrix(msm),
      y = dat$y,
      family = "gaussian",
      link = "identity",
      coefs = coef(msm),
      X_counterfactual = NULL
    ),
    contrasts = NULL,
    focal_level = NULL,
    reference_level = NULL
  )
}

# ---- accessors for the functions under test ---------------------------------

# The three targets are unexported, so the tests reach them through `:::`. Each
# reference is wrapped once here with a targeted jarl suppression: the
# internal_function rule guards against depending on another package's internals,
# which does not apply to a package exercising its own.
new_weight_fn <- function(exposure_type, estimand) {
  # jarl-ignore internal_function: unexported target exercised by the package's own tests
  propensity:::ipw_weight_fn(exposure_type, estimand)
}

new_layout <- function(spec) {
  # jarl-ignore internal_function: unexported target exercised by the package's own tests
  propensity:::ipw_theta_layout(spec)
}

new_psi <- function(spec, layout) {
  # jarl-ignore internal_function: unexported target exercised by the package's own tests
  propensity:::build_ipw_psi(spec, layout)
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
  fn <- new_weight_fn(exposure_type, estimand)
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
  block_order <- c("ps", "stab", "out", "mu", "contrast")
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
  layout <- new_layout(spec)
  psi <- new_psi(spec, layout)
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
  for (estimand in c("ate", "att", "atu", "atm", "ato", "entropy")) {
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
    wt_ate(fitted_ps, A, exposure_type = "continuous")
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
    new_weight_fn("continuous", "att"),
    class = "propensity_ipw_estimand_error"
  )
  expect_error(
    new_weight_fn("continuous", "ato"),
    class = "propensity_ipw_estimand_error"
  )
})

# ---- theta layout -----------------------------------------------------------

test_that("ipw_theta_layout partitions theta for binary specs", {
  # stabilized ate: stab block present
  spec_stab <- binary_spec(sim_binary(), "ate", stabilize = TRUE)
  layout_stab <- new_layout(spec_stab)
  expect_layout_partition(
    layout_stab,
    list(
      ps = ncol(spec_stab$ps$X),
      stab = 1,
      out = ncol(spec_stab$outcome$X),
      mu = 2,
      contrast = 3
    )
  )

  # unstabilized ate: no stab block
  spec_unstab <- binary_spec(sim_binary(), "ate")
  layout_unstab <- new_layout(spec_unstab)
  expect_layout_partition(
    layout_unstab,
    list(
      ps = ncol(spec_unstab$ps$X),
      stab = 0,
      out = ncol(spec_unstab$outcome$X),
      mu = 2,
      contrast = 3
    )
  )

  # gaussian outcome: single linear contrast
  spec_lin <- binary_spec(sim_binary(), "att", outcome_family = "gaussian")
  layout_lin <- new_layout(spec_lin)
  expect_layout_partition(
    layout_lin,
    list(
      ps = ncol(spec_lin$ps$X),
      stab = 0,
      out = ncol(spec_lin$outcome$X),
      mu = 2,
      contrast = 1
    )
  )
})

test_that("ipw_theta_layout partitions theta for categorical specs", {
  skip_if_not_installed("nnet")

  spec <- categorical_spec(sim_categorical(), "ate")
  k <- spec$ps$k
  layout <- new_layout(spec)
  expect_layout_partition(
    layout,
    list(
      ps = ncol(spec$ps$X) * (k - 1),
      stab = 0,
      out = ncol(spec$outcome$X),
      mu = k,
      contrast = (k - 1) * length(spec$contrasts)
    )
  )

  spec_stab <- categorical_spec(sim_categorical(), "ate", stabilize = TRUE)
  layout_stab <- new_layout(spec_stab)
  expect_layout_partition(
    layout_stab,
    list(
      ps = ncol(spec_stab$ps$X) * (k - 1),
      stab = k - 1,
      out = ncol(spec_stab$outcome$X),
      mu = k,
      contrast = (k - 1) * length(spec_stab$contrasts)
    )
  )
})

test_that("ipw_theta_layout partitions theta for a continuous stabilized spec", {
  spec <- continuous_spec(sim_continuous(), stabilize = TRUE)
  layout <- new_layout(spec)
  expect_layout_partition(
    layout,
    list(
      ps = ncol(spec$ps$X) + 1,
      stab = 2,
      out = ncol(spec$outcome$X),
      mu = 0,
      contrast = 0
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

test_that("build_ipw_psi is root-seeded at init for continuous stabilized ate MSM", {
  expect_root_seeded(continuous_spec(sim_continuous(), stabilize = TRUE))
})
