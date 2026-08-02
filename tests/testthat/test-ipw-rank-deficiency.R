# Designs the fits could not separate, and solutions the solver could not pin
# down
#
# `ipw()` multiplies each fitted coefficient vector against its own design
# positionally, so a design whose columns are linearly dependent has no
# coefficient for the columns the fit dropped. `glm()` and `lm()` record that as
# an `NA` coefficient rather than as an error, and the `NA` travels: into the
# propensity scores the M-estimator rebuilds, into the seeded starting values,
# and into the stacked estimating functions. Both models are checked so the
# report names the model and the columns rather than arriving as a missing value
# in a condition, or in the M-estimator's own vocabulary.
#
# Separately, a fit whose estimating equations have no unique root can still
# return numbers. The solver says so; these tests pin that the report the user
# reads is in `ipw()`'s terms and names the effects whose standard errors are
# not meaningful.

rank_data <- function(seed = 2024, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  z <- rbinom(n, 1, plogis(0.3 * x1 - 0.6 * x2))
  y <- rbinom(n, 1, plogis(-0.4 + 1.1 * z + 0.5 * x1 - 0.3 * x2))
  yc <- 1.5 + 0.8 * z + 0.6 * x1 - 0.4 * x2 + rnorm(n)
  dat <- data.frame(x1, x2, z, y, yc)
  # an exact copy and a copy perturbed far below the fits' pivot tolerance; both
  # leave the design rank deficient and the copy without a coefficient
  dat$x1_copy <- dat$x1
  dat$x1_near <- dat$x1 + 1e-12 * rnorm(n)
  dat
}

rank_weights <- function(ps_mod) {
  withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_mod))
}

# Default IRLS control throughout: several of these outcome models are
# deliberately rank deficient, and a tightened tolerance leaves the fit unable
# to reach it and warning about that at fitting time, which is not the condition
# under test.
rank_outcome <- function(formula, data, wts) {
  glm(formula, data = data, family = quasibinomial(), weights = wts)
}

rank_msg <- function(cnd) {
  gsub("[[:space:]]+", " ", conditionMessage(cnd))
}

# ---- a propensity score design the fit could not separate -------------------
#
# The dropped coefficient reached the rebuilt propensity scores as an `NA`,
# which reached the count of saturated scores, which reached an `if`: the call
# stopped with base R's "missing value where TRUE/FALSE needed" and named
# neither model.

test_that("ipw() rejects a propensity score model with an unseparable design", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x1_near + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z + x1, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out, se_method = "mestimation"),
    class = "propensity_ipw_rank_error"
  )
  msg <- rank_msg(err)
  expect_match(msg, "wt_mod", fixed = TRUE)
  expect_match(msg, "x1_near", fixed = TRUE)
  expect_false(grepl("missing value where TRUE/FALSE", msg, fixed = TRUE))
})

test_that("ipw() rejects an exactly duplicated propensity score column", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x1_copy + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z + x1, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out, se_method = "mestimation"),
    class = "propensity_ipw_rank_error"
  )
  expect_match(rank_msg(err), "x1_copy", fixed = TRUE)
})

test_that("the linearization path rejects the same propensity score model", {
  skip_if_not_installed("deli")
  # The guard sits in the shared design extraction, so both standard error
  # methods reach it. This route failed inside LAPACK as an exactly singular
  # system.
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x1_copy + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out, se_method = "linearization"),
    class = "propensity_ipw_rank_error"
  )
  msg <- rank_msg(err)
  expect_match(msg, "wt_mod", fixed = TRUE)
  expect_false(grepl("dgesv", msg, fixed = TRUE))
})

test_that("ipw() rejects an unseparable continuous propensity score design", {
  skip_if_not_installed("deli")
  # An lm propensity model drops a dependent column the same way, and the report
  # is the same one. This reached the weight-consistency preflight, which
  # reported a disagreement about the estimand and the focal level.
  dat <- rank_data()
  withr::local_seed(11)
  dat$a <- 0.5 + 0.8 * dat$x1 + rnorm(nrow(dat))
  ps_mod <- lm(a ~ x1 + x1_copy, data = dat)
  wts <- wt_ate(
    fitted(ps_mod),
    dat$a,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  out <- lm(yc ~ a, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, out, se_method = "mestimation"),
    class = "propensity_ipw_rank_error"
  )
  msg <- rank_msg(err)
  expect_match(msg, "x1_copy", fixed = TRUE)
  expect_false(grepl("weights recomputed", msg, fixed = TRUE))
})

# ---- an outcome design the fit could not separate ---------------------------
#
# These passed every propensity-side guard and stopped inside the M-estimator,
# which reported that `stacked_equations` returned non-finite values at `init`.
# Neither name appears in anything the user wrote.

test_that("ipw() rejects an outcome model with an exactly duplicated column", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z + x1 + x1_copy, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out, se_method = "mestimation"),
    class = "propensity_ipw_rank_error"
  )
  msg <- rank_msg(err)
  expect_match(msg, "outcome_mod", fixed = TRUE)
  expect_match(msg, "x1_copy", fixed = TRUE)
  expect_false(grepl("stacked_equations", msg, fixed = TRUE))
  expect_false(grepl("init", msg, fixed = TRUE))
})

test_that("ipw() rejects a near-collinear outcome design", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z + x1 + x1_near, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out, se_method = "mestimation"),
    class = "propensity_ipw_rank_error"
  )
  expect_match(rank_msg(err), "x1_near", fixed = TRUE)
})

test_that("ipw() rejects an unseparable gaussian outcome design", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- lm(yc ~ z + x1 + x1_copy, data = dat, weights = wts)

  err <- expect_error(
    ipw(ps_mod, out, se_method = "mestimation"),
    class = "propensity_ipw_rank_error"
  )
  expect_match(rank_msg(err), "outcome_mod", fixed = TRUE)
})

test_that("ipw() rejects an unseparable outcome design for a categorical exposure", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- rank_data()
  dat$g <- factor(
    ifelse(dat$x1 < -0.5, "a", ifelse(dat$x1 < 0.5, "b", "c")),
    levels = c("a", "b", "c")
  )
  ps_mod <- nnet::multinom(g ~ x2, data = dat, trace = FALSE)
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- levels(dat$g)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps, dat$g, exposure_type = "categorical")
  )
  out <- rank_outcome(y ~ g + x1 + x1_copy, dat, wts)

  err <- expect_error(
    ipw(ps_mod, out, se_method = "mestimation"),
    class = "propensity_ipw_rank_error"
  )
  expect_match(rank_msg(err), "x1_copy", fixed = TRUE)
})

test_that("the rank-deficient design errors read in the user's terms", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z + x1 + x1_copy, dat, wts)

  expect_snapshot(error = TRUE, ipw(ps_mod, out, se_method = "mestimation"))
})

test_that("the propensity rank error reads in the user's terms", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x1_copy + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z + x1, dat, wts)

  expect_snapshot(error = TRUE, ipw(ps_mod, out, se_method = "mestimation"))
})

# ---- fits the guards must leave alone ---------------------------------------

test_that("a correlated but separable design is still accepted", {
  skip_if_not_installed("deli")
  # Correlation is not rank deficiency. A copy perturbed far enough for the
  # fits to keep a coefficient for it is a legitimate, if poorly conditioned,
  # model and must still run.
  dat <- rank_data()
  withr::local_seed(13)
  dat$x1_loose <- dat$x1 + 0.01 * rnorm(nrow(dat))
  ps_mod <- glm(z ~ x1 + x1_loose + x2, data = dat, family = binomial())
  expect_false(anyNA(coef(ps_mod)))
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z + x1 + x1_loose, dat, wts)
  expect_false(anyNA(coef(out)))

  res <- expect_no_warning(ipw(ps_mod, out, se_method = "mestimation"))
  expect_true(all(is.finite(res$estimates$std.err)))
})

test_that("a categorical propensity model keeps a duplicated column", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  # The rank guard reads the `NA` a fit records for a column it dropped, and
  # `nnet::multinom()` records none: it optimizes rather than pivots and returns
  # a finite coefficient for a dependent column. Its design and its coefficient
  # matrix still agree, which is what the guard is about, so the fit runs and
  # gives the estimates the deduplicated fit gives. What is left of the
  # dependence is reported by the solver.
  dat <- rank_data()
  dat$x2_copy <- dat$x2
  dat$g <- factor(
    ifelse(dat$x1 < -0.5, "a", ifelse(dat$x1 < 0.5, "b", "c")),
    levels = c("a", "b", "c")
  )
  ps_dup <- nnet::multinom(g ~ x2 + x2_copy, data = dat, trace = FALSE)
  ps_plain <- nnet::multinom(g ~ x2, data = dat, trace = FALSE)
  expect_false(anyNA(coef(ps_dup)))

  cat_outcome <- function(mod) {
    ps <- unname(predict(mod, type = "probs"))
    colnames(ps) <- levels(dat$g)
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_ate(ps, dat$g, exposure_type = "categorical")
    )
    glm(y ~ g + x1, data = dat, family = quasibinomial(), weights = wts)
  }
  out_dup <- cat_outcome(ps_dup)
  out_plain <- cat_outcome(ps_plain)

  seen <- character()
  res <- withCallingHandlers(
    ipw(ps_dup, out_dup, se_method = "mestimation"),
    warning = function(w) {
      seen <<- c(seen, class(w)[[1]])
      invokeRestart("muffleWarning")
    }
  )
  ref <- expect_no_warning(
    ipw(ps_plain, out_plain, se_method = "mestimation")
  )

  expect_setequal(seen, "propensity_ipw_solver_warning")
  expect_equal(res$estimates$estimate, ref$estimates$estimate, tolerance = 1e-6)
})

# ---- a solution the estimating equations do not pin down --------------------
#
# A binomial outcome with no events in one exposure arm has no finite fit within
# that arm. The M-estimator still returns numbers, and the standard errors it
# reports for them are many orders of magnitude below the estimates they
# accompany. The only signal was the solver's own warning about rank deficiency
# in the estimating equations, which names neither the effects nor `ipw()`.

zero_event_fit <- function() {
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  dat$y0 <- dat$y
  dat$y0[dat$z == 0] <- 0L
  out <- glm(
    y0 ~ z + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts
  )
  list(ps_mod = ps_mod, out = out)
}

test_that("ipw() warns in its own terms when the equations have no unique root", {
  skip_if_not_installed("deli")
  mods <- zero_event_fit()

  expect_warning(
    ipw(mods$ps_mod, mods$out, se_method = "mestimation"),
    class = "propensity_ipw_solver_warning"
  )
})

test_that("the unsolved-equations warning names the effects it applies to", {
  skip_if_not_installed("deli")
  mods <- zero_event_fit()

  cnd <- NULL
  withCallingHandlers(
    ipw(mods$ps_mod, mods$out, se_method = "mestimation"),
    propensity_ipw_solver_warning = function(w) {
      cnd <<- w
      invokeRestart("muffleWarning")
    }
  )

  msg <- rank_msg(cnd)
  expect_match(msg, "rd", fixed = TRUE)
  expect_match(msg, "log(rr)", fixed = TRUE)
  expect_match(msg, "log(or)", fixed = TRUE)
  expect_match(msg, "not meaningful", fixed = TRUE)
})

test_that("the unsolved-equations warning replaces the solver's own", {
  skip_if_not_installed("deli")
  mods <- zero_event_fit()

  seen <- character()
  withCallingHandlers(
    ipw(mods$ps_mod, mods$out, se_method = "mestimation"),
    warning = function(w) {
      seen <<- c(seen, class(w)[[1]])
      invokeRestart("muffleWarning")
    }
  )

  expect_setequal(seen, "propensity_ipw_solver_warning")
})

test_that("the unsolved-equations warning is raised once per fit", {
  skip_if_not_installed("deli")
  mods <- zero_event_fit()

  n_warnings <- 0L
  withCallingHandlers(
    ipw(mods$ps_mod, mods$out, se_method = "mestimation"),
    propensity_ipw_solver_warning = function(w) {
      n_warnings <<- n_warnings + 1L
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(n_warnings, 1L)
})

test_that("the unsolved-equations warning leaves the estimates in the output", {
  skip_if_not_installed("deli")
  mods <- zero_event_fit()

  res <- suppressWarnings(
    ipw(mods$ps_mod, mods$out, se_method = "mestimation")
  )

  expect_s3_class(res, "ipw")
  expect_identical(nrow(res$estimates), 3L)
  expect_true(all(is.finite(res$estimates$estimate)))
})

test_that("the unsolved-equations warning reads in the user's terms", {
  skip_if_not_installed("deli")
  mods <- zero_event_fit()

  expect_snapshot(ipw(mods$ps_mod, mods$out, se_method = "mestimation"))
})

# ---- a solve that produced no variance at all -------------------------------
#
# The solver can also return estimates the sandwich cannot be built around. When
# the derivatives of the stacked equations come back non-finite there is no
# bread to invert, and deli says so twice in its own vocabulary: a warning about
# a bread matrix holding NA values as it abandons the variance, and, when the
# standard errors are read, an error about a fit that has no variance to compute
# inference from. Neither names anything the user wrote, and the second discards
# the estimates on the way past. These pin the single report that stands in for
# both, and the structure it names when the fit carries one.

# An ill-conditioned design the fits still separate: the perturbation sits an
# order of magnitude above the tolerance the fits pivot at, so every coefficient
# is finite and the rank guard passes it, and the derivatives come back
# non-finite all the same. Nothing in this fit is degenerate in a way that can
# be named.
ill_conditioned_fit <- function() {
  dat <- rank_data()
  dat$x1_dup <- withr::with_seed(3, dat$x1 + 1e-10 * rnorm(nrow(dat)))
  ps_mod <- glm(z ~ x1 + x1_dup + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z + x1, dat, wts)
  list(ps_mod = ps_mod, out = out)
}

# A binomial outcome that is an event for every exposed observation.
all_event_fit <- function() {
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  dat$y1 <- dat$y
  dat$y1[dat$z == 1] <- 1L
  out <- rank_outcome(y1 ~ z + x1, dat, wts)
  list(ps_mod = ps_mod, out = out)
}

# The same shape for a categorical exposure: no events in one level.
no_event_level_fit <- function() {
  dat <- rank_data()
  dat$g <- factor(
    ifelse(dat$x1 < -0.5, "a", ifelse(dat$x1 < 0.5, "b", "c")),
    levels = c("a", "b", "c")
  )
  ps_mod <- nnet::multinom(g ~ x2, data = dat, trace = FALSE)
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- levels(dat$g)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps, dat$g, exposure_type = "categorical")
  )
  dat$y0 <- dat$y
  dat$y0[dat$g == "c"] <- 0L
  out <- rank_outcome(y0 ~ g + x1, dat, wts)
  list(ps_mod = ps_mod, out = out)
}

test_that("ipw() reports a fit whose variance could not be built", {
  skip_if_not_installed("deli")
  mods <- ill_conditioned_fit()
  expect_false(anyNA(coef(mods$ps_mod)))

  err <- expect_error(
    ipw(mods$ps_mod, mods$out, se_method = "mestimation"),
    class = "propensity_ipw_variance_error"
  )
  msg <- rank_msg(err)
  expect_match(msg, "no finite derivative", fixed = TRUE)
  expect_false(grepl("bread", msg, fixed = TRUE))
})

test_that("the no-variance report replaces the ones deli makes", {
  skip_if_not_installed("deli")
  mods <- ill_conditioned_fit()

  seen <- character()
  withCallingHandlers(
    tryCatch(
      ipw(mods$ps_mod, mods$out, se_method = "mestimation"),
      error = function(e) NULL
    ),
    warning = function(w) {
      seen <<- c(seen, class(w)[[1]])
      invokeRestart("muffleWarning")
    }
  )

  expect_length(seen, 0L)
})

test_that("the no-variance error names an exposure arm whose outcome never varies", {
  skip_if_not_installed("deli")
  mods <- all_event_fit()

  err <- suppressWarnings(expect_error(
    ipw(mods$ps_mod, mods$out, se_method = "mestimation"),
    class = "propensity_ipw_variance_error"
  ))
  msg <- rank_msg(err)
  expect_match(msg, "exposed group", fixed = TRUE)
  expect_match(msg, "nothing but events", fixed = TRUE)
})

test_that("the no-variance error names an exposure level whose outcome never varies", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- no_event_level_fit()

  # The level with no events also drives its marginal mean to the boundary,
  # which the contrast reporter reports in its own right.
  err <- suppressWarnings(expect_error(
    ipw(mods$ps_mod, mods$out, se_method = "mestimation"),
    class = "propensity_ipw_variance_error"
  ))
  msg <- rank_msg(err)
  expect_match(msg, "the c level", fixed = TRUE)
  expect_match(msg, "no events", fixed = TRUE)
})

test_that("the no-variance error reads in the user's terms", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  mods <- no_event_level_fit()

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, mods$out, se_method = "mestimation")
  )
})

test_that("the undiagnosable no-variance error reads in the user's terms", {
  skip_if_not_installed("deli")
  mods <- ill_conditioned_fit()

  expect_snapshot(
    error = TRUE,
    ipw(mods$ps_mod, mods$out, se_method = "mestimation")
  )
})

# ---- healthy fits stay quiet ------------------------------------------------

test_that("a healthy binary fit raises no solver warning", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z + x1, dat, wts)

  expect_no_warning(ipw(ps_mod, out, se_method = "mestimation"))
})

test_that("a healthy gaussian outcome fit raises no solver warning", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- lm(yc ~ z + x1, data = dat, weights = wts)

  expect_no_warning(ipw(ps_mod, out, se_method = "mestimation"))
})

test_that("a healthy binary fit raises no solver warning under linearization", {
  dat <- rank_data()
  ps_mod <- glm(z ~ x1 + x2, data = dat, family = binomial())
  wts <- rank_weights(ps_mod)
  out <- rank_outcome(y ~ z, dat, wts)

  expect_no_warning(ipw(ps_mod, out, se_method = "linearization"))
})

test_that("a healthy categorical fit raises no solver warning", {
  skip_if_not_installed("nnet")
  skip_if_not_installed("deli")
  dat <- rank_data()
  dat$g <- factor(
    ifelse(dat$x1 < -0.5, "a", ifelse(dat$x1 < 0.5, "b", "c")),
    levels = c("a", "b", "c")
  )
  ps_mod <- nnet::multinom(g ~ x2, data = dat, trace = FALSE)
  ps <- unname(predict(ps_mod, type = "probs"))
  colnames(ps) <- levels(dat$g)
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps, dat$g, exposure_type = "categorical")
  )
  out <- rank_outcome(y ~ g + x1, dat, wts)

  expect_no_warning(ipw(ps_mod, out, se_method = "mestimation"))
})

test_that("a healthy continuous fit raises no solver warning", {
  skip_if_not_installed("deli")
  dat <- rank_data()
  withr::local_seed(12)
  dat$a <- 0.5 + 0.8 * dat$x1 - 0.4 * dat$x2 + rnorm(nrow(dat))
  ps_mod <- lm(a ~ x1 + x2, data = dat)
  wts <- wt_ate(
    fitted(ps_mod),
    dat$a,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  out <- lm(yc ~ a, data = dat, weights = wts)

  expect_no_warning(ipw(ps_mod, out, se_method = "mestimation"))
})
