# ipw() refuses weights whose values were bounded by wt_trunc().
#
# A weight truncation is a hard bound on the weights, not a smooth function of
# the propensity score model's parameters, so no sandwich `ipw()` can build
# accounts for it. The refusal has its own class, fires on the `wt_truncated`
# flag alone, and fires on every exposure route before the estimand is read,
# since the truncation label would otherwise surface as an unrecognized
# estimand. Weights that are also trimmed, truncated on the score scale, or
# calibrated meet those refusals first.

# ---- fixtures ---------------------------------------------------------------

wt_trunc_ipw_binary_data <- function(seed = 2024, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.5 * x2))
  y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.4 * x1 - 0.3 * x2))
  data.frame(x1, x2, z, y)
}

wt_trunc_ipw_continuous_data <- function(seed = 2024, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  A <- 0.5 + 0.8 * x1 - 0.4 * x2 + rnorm(n)
  yc <- 1 + 0.6 * A + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  data.frame(x1, x2, A, yc)
}

wt_trunc_ipw_categorical_data <- function(seed = 2024, n = 500) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  eta_b <- -0.2 + 0.5 * x1 + 0.3 * x2
  eta_c <- 0.1 - 0.4 * x1 + 0.6 * x2
  denom <- 1 + exp(eta_b) + exp(eta_c)
  u <- runif(n)
  lab <- ifelse(
    u < 1 / denom,
    "a",
    ifelse(u < (1 + exp(eta_b)) / denom, "b", "c")
  )
  a <- factor(lab, levels = c("a", "b", "c"))
  y <- rbinom(n, 1, plogis(-0.3 + 0.4 * (a == "b") + 0.5 * x1))
  data.frame(x1, x2, a, y)
}

wt_trunc_ipw_joint_data <- function(seed = 6401, n = 600) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * x2))
  e <- rbinom(n, 1, plogis(-0.2 + 0.5 * x1 + 0.3 * x2 - 0.8 * a))
  y <- rbinom(n, 1, plogis(-0.5 + 0.7 * a + 0.5 * e + 0.6 * x1))
  data.frame(x1, x2, a, e, y)
}

# Fit the weighted outcome model, recording the formula itself in the call so a
# printed call reads the same on every R version. The formula is rebound to this
# frame so the model frame finds `wts` here, whatever the caller named them.
wt_trunc_ipw_outcome <- function(fmla, dat, wts, family = "binomial") {
  environment(fmla) <- environment()
  mod <- if (family == "binomial") {
    glm(fmla, data = dat, family = quasibinomial(), weights = wts)
  } else {
    lm(fmla, data = dat, weights = wts)
  }
  mod$call$formula <- fmla
  mod
}

wt_trunc_ipw_binary <- function(modify = c("none", "trim", "trunc", "calib")) {
  modify <- rlang::arg_match(modify)
  withr::local_options(propensity.quiet = TRUE)
  dat <- wt_trunc_ipw_binary_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  ps <- predict(ps_mod, type = "response")

  score <- switch(
    modify,
    none = ps,
    trim = ps_refit(
      ps_trim(ps, method = "ps", lower = 0.3, upper = 0.7),
      model = ps_mod
    ),
    trunc = ps_trunc(ps, method = "ps", lower = 0.3, upper = 0.7),
    calib = ps_calibrate(ps, .exposure = dat$z, .focal_level = 1)
  )
  base <- wt_ate(
    score,
    .exposure = dat$z,
    exposure_type = "binary",
    .focal_level = 1
  )
  wts <- wt_trunc(base, method = "wt", upper = 2)

  list(
    dat = dat,
    ps_mod = ps_mod,
    base = base,
    wts = wts,
    outcome_mod = wt_trunc_ipw_outcome(y ~ z, dat, wts)
  )
}

wt_trunc_ipw_continuous <- function(stabilize = TRUE) {
  withr::local_options(propensity.quiet = TRUE)
  dat <- wt_trunc_ipw_continuous_data()
  ps_mod <- lm(A ~ x1 + x2, data = dat)
  base <- wt_ate(
    as.double(fitted(ps_mod)),
    dat$A,
    exposure_type = "continuous",
    stabilize = stabilize
  )
  wts <- wt_trunc(base, method = "count", upper = 5)

  list(
    dat = dat,
    ps_mod = ps_mod,
    wts = wts,
    outcome_mod = wt_trunc_ipw_outcome(yc ~ A, dat, wts, family = "gaussian")
  )
}

wt_trunc_ipw_categorical <- function() {
  withr::local_options(propensity.quiet = TRUE)
  dat <- wt_trunc_ipw_categorical_data()
  ps_mod <- nnet::multinom(
    a ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- ps_mod$lev
  base <- wt_ate(ps, dat$a, exposure_type = "categorical")
  wts <- wt_trunc(base, method = "count", upper = 5)

  list(
    dat = dat,
    ps_mod = ps_mod,
    wts = wts,
    outcome_mod = wt_trunc_ipw_outcome(y ~ a + x1, dat, wts)
  )
}

wt_trunc_ipw_joint <- function() {
  withr::local_options(propensity.quiet = TRUE)
  dat <- wt_trunc_ipw_joint_data()
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- glm(e ~ a * x1 + x2, data = dat, family = binomial())
  base <- wt_joint(
    wt_ate(ps_a),
    wt_ate(ps_e),
    exposure_type = c("binary", "binary")
  )
  wts <- wt_trunc(base, method = "count", upper = 5)

  list(
    dat = dat,
    models = joint_wt_models(a = ps_a, e = ps_e),
    wts = wts,
    outcome_mod = wt_trunc_ipw_outcome(y ~ a * e + x1, dat, wts)
  )
}

# ---- the fixtures bound weights as the tests assume --------------------------

test_that("every fixture carries weights that wt_trunc() moved", {
  fixtures <- list(
    binary = wt_trunc_ipw_binary(),
    continuous_stabilized = wt_trunc_ipw_continuous(stabilize = TRUE),
    continuous_unstabilized = wt_trunc_ipw_continuous(stabilize = FALSE),
    categorical = wt_trunc_ipw_categorical(),
    joint = wt_trunc_ipw_joint()
  )

  for (fx in fixtures) {
    expect_true(is_wt_truncated(fx$wts))
    expect_gt(sum(is_unit_wt_truncated(fx$wts)), 0)
    expect_false(is_ps_trimmed(fx$wts))
    expect_false(is_ps_truncated(fx$wts))
    expect_false(is_ps_calibrated(fx$wts))
    expect_true(is_wt_truncated(extract_weights(fx$outcome_mod)))
  }
})

# ---- one refusal on every route ---------------------------------------------

test_that("ipw() refuses weight-truncated binary weights on both SE paths", {
  fx <- wt_trunc_ipw_binary()

  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, .data = fx$dat),
    class = "propensity_ipw_wt_truncated_error"
  )
  expect_error(
    ipw(
      fx$ps_mod,
      fx$outcome_mod,
      .data = fx$dat,
      se_method = "linearization"
    ),
    class = "propensity_ipw_wt_truncated_error"
  )
})

test_that("ipw() refuses weight-truncated stabilized continuous weights", {
  fx <- wt_trunc_ipw_continuous(stabilize = TRUE)
  expect_true(is_stabilized(fx$wts))

  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, .data = fx$dat),
    class = "propensity_ipw_wt_truncated_error"
  )
})

test_that("ipw() refuses weight-truncated unstabilized continuous weights", {
  fx <- wt_trunc_ipw_continuous(stabilize = FALSE)
  expect_false(is_stabilized(fx$wts))

  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, .data = fx$dat),
    class = "propensity_ipw_wt_truncated_error"
  )
})

test_that("ipw() refuses weight-truncated categorical weights", {
  fx <- wt_trunc_ipw_categorical()

  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, .data = fx$dat),
    class = "propensity_ipw_wt_truncated_error"
  )
})

test_that("ipw() refuses weight-truncated joint weights", {
  fx <- wt_trunc_ipw_joint()

  expect_error(
    ipw(fx$models, fx$outcome_mod),
    class = "propensity_ipw_wt_truncated_error"
  )
})

# ---- the refusal comes before the estimand is read ---------------------------

test_that("the weight-truncation refusal fires before estimand parsing", {
  # A deliberately wrong estimand argument must not win over the guard.
  binary <- wt_trunc_ipw_binary()
  expect_error(
    ipw(
      binary$ps_mod,
      binary$outcome_mod,
      .data = binary$dat,
      estimand = "att"
    ),
    class = "propensity_ipw_wt_truncated_error"
  )

  continuous <- wt_trunc_ipw_continuous()
  expect_error(
    ipw(
      continuous$ps_mod,
      continuous$outcome_mod,
      .data = continuous$dat,
      estimand = "att"
    ),
    class = "propensity_ipw_wt_truncated_error"
  )

  categorical <- wt_trunc_ipw_categorical()
  expect_error(
    ipw(
      categorical$ps_mod,
      categorical$outcome_mod,
      .data = categorical$dat,
      estimand = "ato"
    ),
    class = "propensity_ipw_wt_truncated_error"
  )
})

# ---- the flag alone is enough ------------------------------------------------

test_that("a weight-truncated flag whose record was dropped is still refused", {
  fx <- wt_trunc_ipw_binary()
  n <- length(fx$wts)
  recombined <- c(fx$wts[seq_len(n / 2)], fx$wts[seq(n / 2 + 1, n)])

  # Recombining drops the positional record but keeps the flag.
  expect_true(is_wt_truncated(recombined))
  expect_error(is_unit_wt_truncated(recombined))

  outcome_mod <- wt_trunc_ipw_outcome(y ~ z, fx$dat, recombined)
  expect_error(
    ipw(fx$ps_mod, outcome_mod, .data = fx$dat),
    class = "propensity_ipw_wt_truncated_error"
  )
})

test_that("the flag is refused even under an unlabeled estimand", {
  # Weights that match the propensity score model exactly and carry the plain
  # "ate" label: only the flag says they were truncated, and it must suffice.
  fx <- wt_trunc_ipw_binary()
  flagged <- psw(
    as.double(fx$base),
    estimand = "ate",
    wt_truncated = TRUE
  )
  expect_true(is_wt_truncated(flagged))

  outcome_mod <- wt_trunc_ipw_outcome(y ~ z, fx$dat, flagged)
  expect_error(
    ipw(fx$ps_mod, outcome_mod, .data = fx$dat),
    class = "propensity_ipw_wt_truncated_error"
  )
})

# ---- order against the other modified-weight refusals ------------------------

test_that("trimmed weights that were also weight-truncated meet the trim refusal", {
  fx <- wt_trunc_ipw_binary("trim")
  expect_true(is_ps_trimmed(fx$wts))
  expect_true(is_wt_truncated(fx$wts))

  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, .data = fx$dat),
    class = "propensity_ipw_trimmed_error"
  )
})

test_that("score-truncated weights that were also weight-truncated meet the score refusal", {
  fx <- wt_trunc_ipw_binary("trunc")
  expect_true(is_ps_truncated(fx$wts))
  expect_true(is_wt_truncated(fx$wts))

  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, .data = fx$dat),
    class = "propensity_ipw_truncated_error"
  )
})

test_that("calibrated weights that were also weight-truncated meet the calibration refusal", {
  fx <- wt_trunc_ipw_binary("calib")
  expect_true(is_ps_calibrated(fx$wts))
  expect_true(is_wt_truncated(fx$wts))

  expect_error(
    ipw(fx$ps_mod, fx$outcome_mod, .data = fx$dat),
    class = "propensity_ipw_calibrated_error"
  )
})

# ---- the message -------------------------------------------------------------

test_that("the weight-truncation refusal names both honest routes", {
  fx <- wt_trunc_ipw_binary()

  cnd <- rlang::catch_cnd(ipw(fx$ps_mod, fx$outcome_mod, .data = fx$dat))
  expect_s3_class(cnd, "propensity_ipw_wt_truncated_error")
  msg <- cli::ansi_strip(conditionMessage(cnd))
  expect_match(msg, "truncated propensity score weights", fixed = TRUE)
  expect_match(msg, "wt_trunc()", fixed = TRUE)
  expect_match(msg, "M-estimation", fixed = TRUE)
  expect_match(msg, "fixed-weight", fixed = TRUE)

  # Record the snapshot only once the dedicated refusal is the one raised, so
  # the message of another refusal is never written down as this one's.
  skip_if_not(inherits(cnd, "propensity_ipw_wt_truncated_error"))
  expect_propensity_error(ipw(fx$ps_mod, fx$outcome_mod, .data = fx$dat))
})

test_that("the weight-truncation refusal reads the same on every route", {
  continuous <- wt_trunc_ipw_continuous()
  categorical <- wt_trunc_ipw_categorical()
  joint <- wt_trunc_ipw_joint()

  cnds <- list(
    rlang::catch_cnd(ipw(
      continuous$ps_mod,
      continuous$outcome_mod,
      .data = continuous$dat
    )),
    rlang::catch_cnd(ipw(
      categorical$ps_mod,
      categorical$outcome_mod,
      .data = categorical$dat
    )),
    rlang::catch_cnd(ipw(joint$models, joint$outcome_mod))
  )
  for (cnd in cnds) {
    expect_s3_class(cnd, "propensity_ipw_wt_truncated_error")
  }

  skip_if_not(all(vapply(
    cnds,
    inherits,
    logical(1),
    what = "propensity_ipw_wt_truncated_error"
  )))
  expect_propensity_error(ipw(
    continuous$ps_mod,
    continuous$outcome_mod,
    .data = continuous$dat
  ))
  expect_propensity_error(ipw(
    categorical$ps_mod,
    categorical$outcome_mod,
    .data = categorical$dat
  ))
  expect_propensity_error(ipw(joint$models, joint$outcome_mod))
})
