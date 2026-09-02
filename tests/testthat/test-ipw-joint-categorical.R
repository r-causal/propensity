# A joint intervention on a binary treatment and a categorical one.
#
# The two-model route already accepts a categorical component in its container:
# `joint_wt_models()` types a `nnet::multinom` fit as "categorical" and
# `wt_joint()` multiplies its weights like any other. What this file pins is
# `ipw()` over that pair, which is the part that has to carry a multinomial
# score in the stacked system rather than the binomial one a binary component
# takes.
#
# The surface is the one the declared-exposure route already reports for the
# same crossing, generalized in the second treatment's level count rather than
# redesigned: the counterfactual risk of each of the six cells under the effect
# label "mean", each treatment's simple effects within each level of the other,
# and the interaction of the first treatment against each non-reference level of
# the second. Twenty-four rows for a binary outcome, fifteen for a continuous
# one. As in the two-binary file the expectations are written out here against
# those conventions rather than imported, so a divergence between the routes
# shows up as a disagreement between two written-out contracts.
#
# The centrepiece is the same equivalence: a binomial model of the first
# treatment paired with a multinomial model of the second, and one multinomial
# model over the six cells, are two parameterizations of the same propensity
# score, so the two routes agree loosely in general and to the solvers when both
# parameterizations are saturated in the same covariate. Both gaps are measured
# on these fixtures and the tolerances are set from the measurements.
#
# The sections up to the stabilization slot are unstabilized, since a
# categorical component needs no stabilization to be weighted correctly. What
# the slot itself pins is the three numerators such a component may still carry:
# the marginal proportion of each level, a score the caller computed, and a
# fitted multinomial model, each with its own block of the stacked system.

# ---- data simulators --------------------------------------------------------

# The level a unit takes, drawn from a row of category probabilities. Written
# out rather than reached for, since the simulator needs the draw to be
# reproducible from one seed alongside everything else it draws.
draw_categorical <- function(probs, levels) {
  cumulative <- t(apply(probs, 1, cumsum))
  u <- stats::runif(nrow(probs))
  index <- 1L + rowSums(u > cumulative[, -ncol(cumulative), drop = FALSE])

  factor(levels[index], levels = levels)
}

# A binary treatment and a three-level one sharing confounders, the second
# depending on the first. As in the two-binary fixture the dependence is not
# additive, so the second model has something to model flexibly, and the cell
# counts are far enough from zero that no cell of the crossing is thin.
sim_joint_categorical <- function(seed = 7301, n = 1500) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * x2))

  eta_mid <- -0.1 + 0.5 * x1 + 0.3 * x2 - 0.7 * a + 0.5 * a * x1
  eta_hi <- -0.3 + 0.3 * x1 - 0.4 * x2 + 0.6 * a - 0.4 * a * x1
  denom <- 1 + exp(eta_mid) + exp(eta_hi)
  z <- draw_categorical(
    cbind(1, exp(eta_mid), exp(eta_hi)) / denom,
    c("lo", "mid", "hi")
  )

  z_mid <- as.integer(z == "mid")
  z_hi <- as.integer(z == "hi")
  y <- rbinom(
    n,
    1,
    plogis(
      -0.4 +
        0.7 * a +
        0.5 * z_mid +
        0.3 * z_hi +
        0.6 * x1 -
        0.3 * x2 +
        0.8 * a * z_mid -
        0.5 * a * z_hi
    )
  )
  yc <- 1 +
    0.6 * a +
    0.4 * z_mid +
    0.2 * z_hi +
    0.5 * x1 -
    0.2 * x2 +
    0.7 * a * z_mid -
    0.4 * a * z_hi +
    rnorm(n)

  data.frame(x1, x2, a, z, y, yc)
}

# The same two treatments assigned the other way round, so the factorization the
# container asks for puts the categorical treatment first: the three-level
# treatment is assigned from the covariates and the binary one conditions on it.
sim_joint_categorical_first <- function(seed = 7303, n = 1500) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)

  eta_mid <- -0.1 + 0.5 * x1 + 0.3 * x2
  eta_hi <- -0.3 + 0.3 * x1 - 0.4 * x2
  denom <- 1 + exp(eta_mid) + exp(eta_hi)
  z <- draw_categorical(
    cbind(1, exp(eta_mid), exp(eta_hi)) / denom,
    c("lo", "mid", "hi")
  )

  z_mid <- as.integer(z == "mid")
  z_hi <- as.integer(z == "hi")
  a <- rbinom(
    n,
    1,
    plogis(
      0.2 + 0.3 * x1 - 0.4 * x2 - 0.6 * z_mid + 0.5 * z_hi + 0.4 * z_hi * x1
    )
  )
  y <- rbinom(
    n,
    1,
    plogis(
      -0.4 +
        0.7 * a +
        0.5 * z_mid +
        0.3 * z_hi +
        0.6 * x1 -
        0.3 * x2 +
        0.8 * a * z_mid -
        0.5 * a * z_hi
    )
  )

  data.frame(x1, x2, a, z, y)
}

# One binary covariate and nothing else, so that a binomial model of the first
# treatment paired with a multinomial model of the second, both saturated in it,
# and a multinomial over the six cells saturated in it describe the same joint
# distribution of the two treatments exactly rather than approximately. Both
# parameterizations carry ten free coefficients over the same model, which is
# what makes the tight pin below a statement about the routes rather than about
# two models that happen to be close.
sim_joint_categorical_saturated <- function(seed = 7302, n = 2400) {
  withr::local_seed(seed)
  x2 <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, plogis(0.4 - 0.8 * x2))

  eta_mid <- -0.2 + 0.6 * x2 - 0.8 * a + 0.7 * a * x2
  eta_hi <- -0.4 - 0.5 * x2 + 0.5 * a + 0.6 * a * x2
  denom <- 1 + exp(eta_mid) + exp(eta_hi)
  z <- draw_categorical(
    cbind(1, exp(eta_mid), exp(eta_hi)) / denom,
    c("lo", "mid", "hi")
  )

  z_mid <- as.integer(z == "mid")
  z_hi <- as.integer(z == "hi")
  y <- rbinom(
    n,
    1,
    plogis(
      -0.4 +
        0.7 * a +
        0.5 * z_mid +
        0.3 * z_hi -
        0.3 * x2 +
        0.8 * a * z_mid -
        0.5 * a * z_hi
    )
  )

  data.frame(x2, a, z, y)
}

# ---- model fitting ----------------------------------------------------------

# Every multinomial fit here is solved past its default tolerance. The joint
# comparisons below are between two solvers, and a multinomial stopped at the
# default `reltol` leaves a residual larger than the differences being measured.
# Each is written where its data lives rather than through a fitting helper,
# since a multinomial keeps no model frame and the route rebuilds its design by
# re-evaluating the fitting call in the formula's own environment.

# The categorical component's weights, built from the fitted probabilities under
# the fit's own level names so the weight reads the column belonging to the
# level each unit took rather than the column in that position.
categorical_component_weights <- function(fit, exposure) {
  probabilities <- unname(stats::predict(fit, type = "probs"))
  colnames(probabilities) <- fit$lev

  wt_ate(probabilities, exposure, exposure_type = "categorical")
}

# The two-model route: a binomial model of the binary treatment, a multinomial
# model of the categorical one conditioning on it, and the product of their ate
# weights. The formulas are built here rather than passed in so the weights
# resolve against this frame when the outcome model looks for them.
fit_joint_categorical_route <- function(
  dat,
  a_rhs = c("x1", "x2"),
  z_rhs = "a * x1 + x2",
  outcome_rhs = "a * z + x1",
  outcome_family = "binomial"
) {
  ps_a <- glm(
    stats::reformulate(a_rhs, response = "a"),
    data = dat,
    family = binomial()
  )
  ps_z <- nnet::multinom(
    stats::reformulate(z_rhs, response = "z"),
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      categorical_component_weights(ps_z, dat$z),
      exposure_type = c("binary", "categorical")
    )
  )

  outcome_var <- if (outcome_family == "binomial") "y" else "yc"
  fmla <- stats::reformulate(outcome_rhs, response = outcome_var)
  outcome_mod <- if (outcome_family == "binomial") {
    glm(
      fmla,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  } else {
    lm(fmla, data = dat, weights = wts)
  }

  list(
    models = joint_wt_models(a = ps_a, z = ps_z),
    ps_a = ps_a,
    ps_z = ps_z,
    outcome_mod = outcome_mod,
    wts = wts
  )
}

# The declared-exposure route over the same data: the crossing as one factor,
# one multinomial propensity score model over the six cells, categorical ate
# weights, and an outcome model over the cells. `y ~ joint + x1` and
# `y ~ a * z + x1` span the same seven columns, so given the same weights the
# two outcome fits coincide and the only thing separating the routes is the
# propensity score parameterization.
fit_joint_categorical_exposure_route <- function(
  dat,
  ps_rhs = c("x1", "x2"),
  outcome_rhs = "joint + x1",
  outcome_family = "binomial"
) {
  dat$joint <- causalgenerics::joint_exposure(a = dat$a, z = dat$z)
  ps_mod <- nnet::multinom(
    stats::reformulate(ps_rhs, response = "joint"),
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    categorical_component_weights(ps_mod, dat$joint)
  )

  outcome_var <- if (outcome_family == "binomial") "y" else "yc"
  fmla <- stats::reformulate(outcome_rhs, response = outcome_var)
  outcome_mod <- if (outcome_family == "binomial") {
    glm(
      fmla,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  } else {
    lm(fmla, data = dat, weights = wts)
  }

  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts, dat = dat)
}

# ---- the reported surface ---------------------------------------------------

# The cells, in the order the crossing varies them: the first treatment fastest,
# the reference cell first.
joint_categorical_cells <- c(
  "a = 0, z = lo",
  "a = 1, z = lo",
  "a = 0, z = mid",
  "a = 1, z = mid",
  "a = 0, z = hi",
  "a = 1, z = hi"
)

# The simple effects and the interaction, in the order the surface reports them:
# the first treatment varied within each level of the second, then the second
# within each level of the first, the held level outer and the compared level
# inner; then the interaction, once under the first treatment's framing, against
# each non-reference level of the second.
#
# `hi` and `lo` index `joint_categorical_cells` for a simple effect and index
# this list for an interaction, which differences two simple effects rather than
# four cell means.
joint_categorical_entries <- list(
  list(contrast = "a: 1 vs 0", group = "z = lo", hi = 2L, lo = 1L),
  list(contrast = "a: 1 vs 0", group = "z = mid", hi = 4L, lo = 3L),
  list(contrast = "a: 1 vs 0", group = "z = hi", hi = 6L, lo = 5L),
  list(contrast = "z: mid vs lo", group = "a = 0", hi = 3L, lo = 1L),
  list(contrast = "z: hi vs lo", group = "a = 0", hi = 5L, lo = 1L),
  list(contrast = "z: mid vs lo", group = "a = 1", hi = 4L, lo = 2L),
  list(contrast = "z: hi vs lo", group = "a = 1", hi = 6L, lo = 2L),
  list(
    contrast = "a: 1 vs 0",
    group = "z = mid vs z = lo",
    hi = 2L,
    lo = 1L
  ),
  list(
    contrast = "a: 1 vs 0",
    group = "z = hi vs z = lo",
    hi = 3L,
    lo = 1L
  )
)

joint_categorical_simple <- 7L

joint_categorical_expected_rows <- function(forms = c("rd", "log(rr)")) {
  entries <- joint_categorical_entries

  data.frame(
    effect = c(
      rep("mean", length(joint_categorical_cells)),
      rep(forms, times = length(entries))
    ),
    contrast = c(
      joint_categorical_cells,
      rep(
        vapply(entries, function(x) x$contrast, character(1)),
        each = length(forms)
      )
    ),
    group = c(
      rep("overall", length(joint_categorical_cells)),
      rep(
        vapply(entries, function(x) x$group, character(1)),
        each = length(forms)
      )
    ),
    stringsAsFactors = FALSE
  )
}

joint_categorical_labels <- function(estimates) {
  paste(estimates$effect, estimates$contrast, estimates$group)
}

# The six counterfactual risks by g-computation on a weighted outcome model that
# reads the two treatments separately, standardized over the whole sample. The
# ate tilt is the constant one, so the averaging population is the whole sample
# and nothing is weighted here.
joint_categorical_cell_means <- function(outcome_mod, data) {
  settings <- list(
    c("0", "lo"),
    c("1", "lo"),
    c("0", "mid"),
    c("1", "mid"),
    c("0", "hi"),
    c("1", "hi")
  )

  out <- vapply(
    settings,
    function(cell) {
      d <- data
      d$a <- as.numeric(cell[[1]])
      d$z <- factor(cell[[2]], levels = levels(data$z))
      mean(stats::predict(outcome_mod, newdata = d, type = "response"))
    },
    numeric(1)
  )

  stats::setNames(out, joint_categorical_cells)
}

# The same six risks off the declared route's outcome model, which reads the
# crossing as one factor rather than the two treatments separately. The cell is
# set on that factor directly: a crossing holding one cell is not a crossing at
# all and cannot be rebuilt as a declared exposure, so what a counterfactual
# setting has to reach here is the column the outcome model reads.
joint_categorical_declared_cell_means <- function(outcome_mod, data) {
  out <- vapply(
    joint_categorical_cells,
    function(cell) {
      d <- data
      d$joint <- factor(cell, levels = levels(data$joint))
      mean(stats::predict(outcome_mod, newdata = d, type = "response"))
    },
    numeric(1)
  )

  stats::setNames(out, joint_categorical_cells)
}

joint_categorical_transform <- function(form, mu_hi, mu_lo) {
  switch(
    form,
    rd = mu_hi - mu_lo,
    diff = mu_hi - mu_lo,
    "log(rr)" = log(mu_hi) - log(mu_lo)
  )
}

# Every reported estimate, built from the six cell means in the row order above.
# The interaction rows difference two simple effects on their own scale rather
# than differencing four means, which is what the surface reports them as.
joint_categorical_expected_estimates <- function(
  mu,
  forms = c("rd", "log(rr)")
) {
  entries <- joint_categorical_entries
  simple <- lapply(entries[seq_len(joint_categorical_simple)], function(entry) {
    vapply(
      forms,
      function(f) {
        joint_categorical_transform(f, mu[[entry$hi]], mu[[entry$lo]])
      },
      numeric(1)
    )
  })
  interaction <- lapply(
    entries[-seq_len(joint_categorical_simple)],
    function(entry) simple[[entry$hi]] - simple[[entry$lo]]
  )

  unname(c(mu, unlist(simple), unlist(interaction)))
}

# The names the multinomial denominator block takes, which is the layout of
# `as.vector(t(coef(fit)))`: one run of terms per non-reference level, in the
# fit's own level order, each written `<level>:<term>`.
joint_categorical_ps_names <- function(fit) {
  coefs <- stats::coef(fit)
  terms <- colnames(coefs)

  as.vector(vapply(
    rownames(coefs),
    function(level) paste0(level, ":", terms),
    character(length(terms))
  ))
}

# ---- baselines: the halves that stand today ---------------------------------

test_that("the fixture fits both routes and the container accepts the pair", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()

  # Every cell of the crossing is populated, which the declared-exposure route
  # requires and which the two-model route needs a mean for.
  expect_true(all(table(dat$a, dat$z) > 0))

  two <- fit_joint_categorical_route(dat)
  expect_true(is_joint_wt_models(two$models))
  expect_identical(two$models$names, c("a", "z"))
  expect_identical(
    two$models$exposure_type,
    c(a = "binary", z = "categorical")
  )

  # The weights are the product the container's models imply, and they say so.
  expect_true(is_joint_wt(two$wts))
  expect_identical(estimand(two$wts), "ate")
  expect_equal(
    as.double(two$wts),
    as.double(wt_ate(two$ps_a)) *
      as.double(withr::with_options(
        list(propensity.quiet = TRUE),
        categorical_component_weights(two$ps_z, dat$z)
      )),
    tolerance = 1e-12
  )

  # The declared-exposure route reports this crossing today, so its half of the
  # equivalence comparison is a fact about the package rather than an
  # aspiration, and the twenty-four rows are the surface it already reports.
  joint <- fit_joint_categorical_exposure_route(dat)
  res <- ipw(joint$ps_mod, joint$outcome_mod)
  expected <- joint_categorical_expected_rows()
  expect_identical(nrow(res$estimates), 24L)
  expect_identical(res$estimates$effect, expected$effect)
  expect_identical(res$estimates$contrast, expected$contrast)
  expect_identical(res$estimates$group, expected$group)

  # The written-out surface is the g-computation it claims to be, which is what
  # makes it usable as the comparator the two-model route is held to.
  expect_equal(
    res$estimates$estimate,
    joint_categorical_expected_estimates(
      joint_categorical_declared_cell_means(joint$outcome_mod, joint$dat)
    ),
    tolerance = 1e-8
  )
})

test_that("the two parameterizations disagree about the weights they imply", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  two <- fit_joint_categorical_route(dat)
  joint <- fit_joint_categorical_exposure_route(dat)

  # What makes the loose equivalence below a statement about estimators rather
  # than about arithmetic: a binomial model paired with a multinomial one and a
  # multinomial over the six cells are different models of the same conditional
  # distribution, and on this fixture the weights they imply differ by up to 48
  # percent. Two routes agreeing to a few percent in the estimates after that is
  # the claim being made.
  gap <- max(
    abs(as.double(two$wts) - as.double(joint$wts)) / as.double(joint$wts)
  )
  expect_gt(gap, 0.1)
})

test_that("saturated parameterizations imply the same weights", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_saturated()
  two <- fit_joint_categorical_route(
    dat,
    a_rhs = "x2",
    z_rhs = "a * x2",
    outcome_rhs = "a * z + x2"
  )
  joint <- fit_joint_categorical_exposure_route(
    dat,
    ps_rhs = "x2",
    outcome_rhs = "joint + x2"
  )

  # `a ~ x2` with `z ~ a * x2` and a multinomial over the six cells on `x2` are
  # two parameterizations of the same saturated model, so they fit the same cell
  # probabilities and the weights coincide. Measured max relative gap on this
  # fixture is 3.6e-7, which is where the two multinomial solves stop.
  gap <- max(
    abs(as.double(two$wts) - as.double(joint$wts)) / as.double(joint$wts)
  )
  expect_lt(gap, 1e-5)
})

# ---- the type gate ----------------------------------------------------------

test_that("the joint two-model route admits a categorical component in either slot", {
  # The stacked system carries a multinomial score for a categorical treatment,
  # and a multinomial score reads no order into the pair: unlike a dose, whose
  # density the first factor of the factorization does not carry, a categorical
  # treatment is the same block whichever position it sits in.
  expect_silent(
    check_ipw_joint_models_types(c(a = "binary", z = "categorical"))
  )
  expect_silent(
    check_ipw_joint_models_types(c(z = "categorical", a = "binary"))
  )

  # The boundary the admission is measured against: a dose in the first slot is
  # still refused, and the refusal still names the component that is wrong and
  # the position it sits in.
  cnd <- tryCatch(
    check_ipw_joint_models_types(c(d = "continuous", a = "binary")),
    error = function(e) e
  )
  expect_s3_class(cnd, "propensity_ipw_exposure_error")
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`d`", fixed = TRUE)
  expect_match(message, "first", fixed = TRUE)
})

test_that("ipw() no longer refuses a categorical component by type", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  two <- fit_joint_categorical_route(dat)

  # The refusal this route used to raise for the pair, which is the one thing
  # every expectation below depends on not happening.
  expect_no_error(ipw(two$models, two$outcome_mod))
})

test_that("ipw() reports the joint surface with the categorical treatment first", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_first()

  # The treatments assigned the other way round, so the factorization the
  # container records is f(z | L) f(a | z, L) and the categorical component is
  # the one the second model conditions on.
  ps_z <- nnet::multinom(
    z ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_a <- glm(a ~ z * x1 + x2, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      categorical_component_weights(ps_z, dat$z),
      wt_ate(ps_a),
      exposure_type = c("categorical", "binary")
    )
  )
  outcome_mod <- glm(
    y ~ a * z + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  models <- joint_wt_models(z = ps_z, a = ps_a)
  expect_identical(models$exposure_type, c(z = "categorical", a = "binary"))

  res <- ipw(models, outcome_mod)

  # The crossing has the same six cells whichever treatment is varied fastest,
  # so the surface is the same size; what moves is which treatment the
  # interaction rows are framed under, which this ordering puts on `z`. What
  # this test pins is the position rather than those labels: a multinomial
  # score is the same block in either slot, so the route reports the crossing
  # here rather than refusing it.
  expect_s3_class(res, "ipw")
  expect_identical(nrow(res$estimates), 24L)
  expect_true(all(is.finite(res$estimates$std.err)))
  expect_true(all(res$estimates$std.err > 0))
  expect_false(any(res$estimates$effect == "log(or)"))
})

# ---- the reported surface ---------------------------------------------------

test_that("ipw() over a binary and a categorical treatment reports the joint surface", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  two <- fit_joint_categorical_route(dat)
  res <- ipw(two$models, two$outcome_mod)
  est <- res$estimates

  expect_s3_class(res, "ipw")
  expect_identical(res$estimand, "ate")
  expect_identical(res$se_method, "mestimation")

  expect_identical(
    names(est),
    c(
      "effect",
      "contrast",
      "group",
      "estimate",
      "std.err",
      "z",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )

  expected <- joint_categorical_expected_rows()
  expect_identical(nrow(est), 24L)
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)

  # No odds ratio anywhere: it is noncollapsible, so neither a simple effect
  # reported beside one nor a difference of two of them says what it appears to.
  # A three-level second treatment adds rows to the surface and changes nothing
  # about that argument.
  expect_false(any(est$effect == "log(or)"))

  # Twenty-four rows over three identity columns is enough room for two rows to
  # collide, and every reader of the fit keys its rows by the label.
  labels <- joint_categorical_labels(est)
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(names(coef(res)), labels)
  expect_identical(dimnames(vcov(res)), list(labels, labels))
})

test_that("ipw() over a binary and a categorical treatment hand-computes its cell risks and contrasts", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  two <- fit_joint_categorical_route(dat)
  res <- ipw(two$models, two$outcome_mod)

  # Every reported quantity is g-computation on the weighted outcome model,
  # standardized over the whole sample, with the two treatments set together.
  mu <- joint_categorical_cell_means(two$outcome_mod, dat)
  expect_equal(
    res$estimates$estimate,
    joint_categorical_expected_estimates(mu),
    tolerance = 1e-8
  )

  # The two interactions written out over the cells, so the identity the labels
  # claim is pinned against arithmetic rather than against the helper above.
  # Each is a double difference against the reference level of the second
  # treatment; every other difference-in-differences on this crossing is a
  # difference of these two.
  interaction <- res$estimates$estimate[
    res$estimates$effect == "rd" &
      res$estimates$group %in% c("z = mid vs z = lo", "z = hi vs z = lo")
  ]
  expect_equal(
    interaction,
    unname(c(
      mu[[4]] - mu[[3]] - (mu[[2]] - mu[[1]]),
      mu[[6]] - mu[[5]] - (mu[[2]] - mu[[1]])
    )),
    tolerance = 1e-8
  )
})

test_that("ipw() over a binary and a categorical treatment reports diffs for a continuous outcome", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  two <- fit_joint_categorical_route(dat, outcome_family = "gaussian")
  res <- ipw(two$models, two$outcome_mod)
  est <- res$estimates

  expected <- joint_categorical_expected_rows(forms = "diff")
  expect_identical(nrow(est), 15L)
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)

  mu <- joint_categorical_cell_means(two$outcome_mod, dat)
  expect_equal(
    est$estimate,
    joint_categorical_expected_estimates(mu, forms = "diff"),
    tolerance = 1e-8
  )
})

# ---- the denominator block --------------------------------------------------

test_that("the stacked system carries the multinomial denominator block", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  two <- fit_joint_categorical_route(dat)
  res <- ipw(two$models, two$outcome_mod)
  theta <- coef(res$fit)

  # A multinomial fit contributes one parameter per term per non-reference
  # level, so the block is (k - 1) times wider than the design rather than as
  # wide as it. Three coefficients for `a ~ x1 + x2`, ten for `z ~ a * x1 + x2`
  # over two non-reference levels, seven for the outcome model, six cell means,
  # and eighteen contrast rows.
  expect_identical(as.integer(res$fit@n_params), 44L)
  expect_identical(nobs(res), 1500L)
  expect_identical(df.residual(res), 1500L - 44L)

  # The block is named the way the single-treatment categorical route names its
  # own: `<level>:<term>`, level-major and term-minor, which is the order
  # `as.vector(t(coef()))` puts the coefficients in and the order
  # `nnet::multinom()` names its own covariance in.
  ps_names <- joint_categorical_ps_names(two$ps_z)
  expect_identical(
    ps_names,
    c(
      "mid:(Intercept)",
      "mid:a",
      "mid:x1",
      "mid:x2",
      "mid:a:x1",
      "hi:(Intercept)",
      "hi:a",
      "hi:x1",
      "hi:x2",
      "hi:a:x1"
    )
  )
  expect_identical(ps_names, rownames(vcov(two$ps_z)))

  # The two denominator blocks come first and in the components' own order, the
  # binary one under its coefficient names and the categorical one under the
  # composed ones.
  n_a <- length(coef(two$ps_a))
  expect_identical(names(theta)[seq_len(n_a)], names(coef(two$ps_a)))
  expect_identical(
    names(theta)[n_a + seq_along(ps_names)],
    ps_names
  )

  # The block is solved at the multinomial the caller fit rather than carried at
  # whatever it was seeded with, which is what makes the reported standard
  # errors account for having estimated it.
  expect_equal(
    unname(theta[n_a + seq_along(ps_names)]),
    as.vector(t(coef(two$ps_z))),
    tolerance = 1e-6
  )
  expect_equal(
    unname(theta[seq_len(n_a)]),
    unname(coef(two$ps_a)),
    tolerance = 1e-6
  )

  # The block's width tracks the multinomial rather than carrying a fixed
  # allowance for it: an additive second model has two coefficients fewer, one
  # per non-reference level, and the system is two parameters smaller.
  additive <- fit_joint_categorical_route(dat, z_rhs = "a + x1 + x2")
  res_additive <- ipw(additive$models, additive$outcome_mod)
  expect_identical(as.integer(res_additive$fit@n_params), 42L)
})

test_that("a categorical joint fit reports a usable standard error for every row", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  two <- fit_joint_categorical_route(dat)
  res <- ipw(two$models, two$outcome_mod)
  est <- res$estimates

  expect_true(all(is.finite(est$std.err)))
  expect_true(all(est$std.err > 0))
  expect_true(all(est$ci.lower < est$ci.upper))
  expect_equal(sqrt(diag(vcov(res))), est$std.err, ignore_attr = TRUE)
  expect_equal(vcov(res), t(vcov(res)), tolerance = 1e-12)
})

# ---- the weights the system starts from -------------------------------------

test_that("the categorical joint seed rebuilds the weights the outcome model was fit under", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  two <- fit_joint_categorical_route(dat)

  # The preflight rebuilds the product from the spec's own blocks at the seed,
  # so it is the check that the multinomial block evaluates to the probabilities
  # the categorical component's weights were built from. A block that read the
  # softmax the wrong way round, or that read the column in the level's position
  # rather than the column named for it, would reach a different product here.
  spec <- ipw_spec_joint_models(two$models, two$outcome_mod)
  layout <- ipw_theta_layout(spec)
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(two$wts),
    tolerance = 1e-8
  )

  # Every equation in the stacked system sits at its root at the seed, which is
  # what makes the solve report the models the caller fit.
  mat <- build_ipw_psi(spec, layout)(layout$init)
  expect_false(anyNA(mat))
  expect_true(all(abs(rowSums(mat)) / spec$n < 1e-8))
})

# ---- equivalence with the declared-exposure route ---------------------------

test_that("the two routes agree as estimators of the same joint effects", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()

  two <- fit_joint_categorical_route(dat)
  joint <- fit_joint_categorical_exposure_route(dat)
  res_two <- ipw(two$models, two$outcome_mod)
  res_joint <- ipw(joint$ps_mod, joint$outcome_mod)

  # The same rows in the same order, so the comparison is row for row and label
  # for label.
  expect_identical(res_two$estimates$effect, res_joint$estimates$effect)
  expect_identical(res_two$estimates$contrast, res_joint$estimates$contrast)
  expect_identical(res_two$estimates$group, res_joint$estimates$group)
  expect_identical(
    joint_categorical_labels(res_two$estimates),
    joint_categorical_labels(res_joint$estimates)
  )

  # They target the same estimand through different models of the same
  # propensity score, so what separates them is the difference between those
  # models rather than solver noise. Measured on this fixture the reported rows
  # differ by at most 1.6e-2 in absolute terms, on cell risks of 0.37 to 0.75
  # and contrasts up to 0.38, while the weights behind them differ by up to 48
  # percent. The bound is set from that measurement with room, and what it pins
  # is that the two routes answer the same question, not that they compute the
  # same number.
  expect_lt(
    max(abs(res_two$estimates$estimate - res_joint$estimates$estimate)),
    0.05
  )

  # The cell risks are the primary quantities and they are all probabilities, so
  # they carry a tighter bound of their own; measured 8e-3.
  means <- res_two$estimates$effect == "mean"
  expect_lt(
    max(abs(
      res_two$estimates$estimate[means] - res_joint$estimates$estimate[means]
    )),
    0.03
  )
})

test_that("the two routes agree for a continuous outcome", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()

  two <- fit_joint_categorical_route(dat, outcome_family = "gaussian")
  joint <- fit_joint_categorical_exposure_route(
    dat,
    outcome_family = "gaussian"
  )
  res_two <- ipw(two$models, two$outcome_mod)
  res_joint <- ipw(joint$ps_mod, joint$outcome_mod)

  expect_identical(nrow(res_two$estimates), 15L)
  expect_identical(res_two$estimates$effect, res_joint$estimates$effect)
  expect_identical(res_two$estimates$contrast, res_joint$estimates$contrast)
  expect_identical(res_two$estimates$group, res_joint$estimates$group)

  # Measured 2.0e-2 on this fixture, on means of 0.93 to 2.50 and contrasts up
  # to 1.16.
  expect_lt(
    max(abs(res_two$estimates$estimate - res_joint$estimates$estimate)),
    0.06
  )
})

test_that("saturated parameterizations make the two routes agree tightly", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_saturated()

  two <- fit_joint_categorical_route(
    dat,
    a_rhs = "x2",
    z_rhs = "a * x2",
    outcome_rhs = "a * z + x2"
  )
  joint <- fit_joint_categorical_exposure_route(
    dat,
    ps_rhs = "x2",
    outcome_rhs = "joint + x2"
  )
  res_two <- ipw(two$models, two$outcome_mod)
  res_joint <- ipw(joint$ps_mod, joint$outcome_mod)

  # With both parameterizations saturated in the one covariate they describe the
  # same joint distribution exactly, so nothing separates the routes but the two
  # solves. This is the strongest available statement that the two-model route
  # computes the same estimand rather than something that resembles it: the
  # loose comparison above would be satisfied by an estimator that was merely
  # close, and this one would not.
  #
  # Measured on this fixture the rows agree to 3.0e-8 absolute against weights
  # agreeing to 3.6e-7, which is where the two multinomial solves stop. The
  # bound is set well above the measurement to leave room for them.
  expect_equal(
    res_two$estimates$estimate,
    res_joint$estimates$estimate,
    tolerance = 1e-5
  )

  # The two saturated parameterizations carry the same number of free
  # propensity score parameters over the same model, so one is a smooth
  # reparameterization of the other and the sandwich for a functional of the
  # cell probabilities is invariant to which is written down. The standard
  # errors therefore agree as well as the estimates do, which is what makes this
  # a statement about the stacked system rather than about the point estimates
  # it happens to report.
  expect_equal(
    res_two$estimates$std.err,
    res_joint$estimates$std.err,
    tolerance = 1e-4
  )
})

# ---- a categorical treatment in both slots ----------------------------------

# Two three-level treatments sharing confounders, the second depending on the
# first. Nothing about the stacked system distinguishes the pair from the one
# above: each component carries a multinomial score, and the crossing is nine
# cells rather than six. What the fixture is for is the surface, which grows in
# both level counts at once, and the level counts are kept apart from the
# outcome model's own so a row read off the wrong axis would not line up.
sim_joint_categorical_pair <- function(seed = 7311, n = 2400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)

  eta_mid <- -0.1 + 0.5 * x1 + 0.3 * x2
  eta_hi <- -0.3 + 0.3 * x1 - 0.4 * x2
  denom <- 1 + exp(eta_mid) + exp(eta_hi)
  z1 <- draw_categorical(
    cbind(1, exp(eta_mid), exp(eta_hi)) / denom,
    c("lo", "mid", "hi")
  )
  z1_mid <- as.integer(z1 == "mid")
  z1_hi <- as.integer(z1 == "hi")

  eta_q <- -0.2 +
    0.4 * x1 +
    0.2 * x2 -
    0.6 * z1_mid +
    0.5 * z1_hi +
    0.3 * z1_mid * x1
  eta_r <- -0.3 -
    0.3 * x1 +
    0.3 * x2 +
    0.5 * z1_mid -
    0.4 * z1_hi -
    0.3 * z1_hi * x1
  denom2 <- 1 + exp(eta_q) + exp(eta_r)
  z2 <- draw_categorical(
    cbind(1, exp(eta_q), exp(eta_r)) / denom2,
    c("p", "q", "r")
  )
  z2_q <- as.integer(z2 == "q")
  z2_r <- as.integer(z2 == "r")

  y <- rbinom(
    n,
    1,
    plogis(
      -0.4 +
        0.6 * z1_mid +
        0.3 * z1_hi +
        0.5 * z2_q +
        0.2 * z2_r +
        0.6 * x1 -
        0.3 * x2 +
        0.7 * z1_mid * z2_q -
        0.4 * z1_hi * z2_r
    )
  )

  data.frame(x1, x2, z1, z2, y)
}

# The two-model route over the pair: a multinomial model of each treatment, the
# second conditioning on the first, and the product of their ate weights.
fit_joint_categorical_pair_route <- function(dat) {
  ps_z1 <- nnet::multinom(
    z1 ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_z2 <- nnet::multinom(
    z2 ~ z1 * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      categorical_component_weights(ps_z1, dat$z1),
      categorical_component_weights(ps_z2, dat$z2),
      exposure_type = c("categorical", "categorical")
    )
  )

  # `y ~ z1 * z2 + x1` spans the nine cells and the covariate, which is what
  # `y ~ joint + x1` spans on the other route.
  outcome_mod <- glm(
    y ~ z1 * z2 + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(
    models = joint_wt_models(z1 = ps_z1, z2 = ps_z2),
    ps_z1 = ps_z1,
    ps_z2 = ps_z2,
    outcome_mod = outcome_mod,
    wts = wts
  )
}

# The declared-exposure route over the same crossing: one multinomial model over
# the nine cells and an outcome model that reads them as one factor.
fit_joint_categorical_pair_exposure_route <- function(dat) {
  dat$joint <- causalgenerics::joint_exposure(z1 = dat$z1, z2 = dat$z2)
  ps_mod <- nnet::multinom(
    joint ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    categorical_component_weights(ps_mod, dat$joint)
  )
  outcome_mod <- glm(
    y ~ joint + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts, dat = dat)
}

# The nine cells, the first treatment varying fastest and the reference cell
# first, as for any other crossing.
joint_categorical_pair_cells <- c(
  "z1 = lo, z2 = p",
  "z1 = mid, z2 = p",
  "z1 = hi, z2 = p",
  "z1 = lo, z2 = q",
  "z1 = mid, z2 = q",
  "z1 = hi, z2 = q",
  "z1 = lo, z2 = r",
  "z1 = mid, z2 = r",
  "z1 = hi, z2 = r"
)

# Twelve simple effects and four interactions. The simple effects vary the first
# treatment within each level of the second and then the second within each
# level of the first, held level outer and compared level inner; the four
# interactions cross each non-reference level of one treatment with each of the
# other, the second treatment's level outer, and are reported under the first
# treatment's framing.
joint_categorical_pair_entries <- list(
  list(contrast = "z1: mid vs lo", group = "z2 = p", hi = 2L, lo = 1L),
  list(contrast = "z1: hi vs lo", group = "z2 = p", hi = 3L, lo = 1L),
  list(contrast = "z1: mid vs lo", group = "z2 = q", hi = 5L, lo = 4L),
  list(contrast = "z1: hi vs lo", group = "z2 = q", hi = 6L, lo = 4L),
  list(contrast = "z1: mid vs lo", group = "z2 = r", hi = 8L, lo = 7L),
  list(contrast = "z1: hi vs lo", group = "z2 = r", hi = 9L, lo = 7L),
  list(contrast = "z2: q vs p", group = "z1 = lo", hi = 4L, lo = 1L),
  list(contrast = "z2: r vs p", group = "z1 = lo", hi = 7L, lo = 1L),
  list(contrast = "z2: q vs p", group = "z1 = mid", hi = 5L, lo = 2L),
  list(contrast = "z2: r vs p", group = "z1 = mid", hi = 8L, lo = 2L),
  list(contrast = "z2: q vs p", group = "z1 = hi", hi = 6L, lo = 3L),
  list(contrast = "z2: r vs p", group = "z1 = hi", hi = 9L, lo = 3L),
  list(
    contrast = "z1: mid vs lo",
    group = "z2 = q vs z2 = p",
    hi = 3L,
    lo = 1L
  ),
  list(
    contrast = "z1: hi vs lo",
    group = "z2 = q vs z2 = p",
    hi = 4L,
    lo = 2L
  ),
  list(
    contrast = "z1: mid vs lo",
    group = "z2 = r vs z2 = p",
    hi = 5L,
    lo = 1L
  ),
  list(
    contrast = "z1: hi vs lo",
    group = "z2 = r vs z2 = p",
    hi = 6L,
    lo = 2L
  )
)

joint_categorical_pair_simple <- 12L

joint_categorical_pair_expected_rows <- function(forms = c("rd", "log(rr)")) {
  entries <- joint_categorical_pair_entries

  data.frame(
    effect = c(
      rep("mean", length(joint_categorical_pair_cells)),
      rep(forms, times = length(entries))
    ),
    contrast = c(
      joint_categorical_pair_cells,
      rep(
        vapply(entries, function(x) x$contrast, character(1)),
        each = length(forms)
      )
    ),
    group = c(
      rep("overall", length(joint_categorical_pair_cells)),
      rep(
        vapply(entries, function(x) x$group, character(1)),
        each = length(forms)
      )
    ),
    stringsAsFactors = FALSE
  )
}

# The nine counterfactual risks by g-computation on the weighted outcome model,
# standardized over the whole sample, with both treatments set together.
joint_categorical_pair_cell_means <- function(outcome_mod, data) {
  settings <- expand.grid(
    z1 = levels(data$z1),
    z2 = levels(data$z2),
    stringsAsFactors = FALSE
  )

  out <- vapply(
    seq_len(nrow(settings)),
    function(i) {
      d <- data
      d$z1 <- factor(settings$z1[[i]], levels = levels(data$z1))
      d$z2 <- factor(settings$z2[[i]], levels = levels(data$z2))
      mean(stats::predict(outcome_mod, newdata = d, type = "response"))
    },
    numeric(1)
  )

  stats::setNames(out, joint_categorical_pair_cells)
}

joint_categorical_pair_expected_estimates <- function(
  mu,
  forms = c("rd", "log(rr)")
) {
  entries <- joint_categorical_pair_entries
  simple <- lapply(
    entries[seq_len(joint_categorical_pair_simple)],
    function(entry) {
      vapply(
        forms,
        function(f) {
          joint_categorical_transform(f, mu[[entry$hi]], mu[[entry$lo]])
        },
        numeric(1)
      )
    }
  )
  interaction <- lapply(
    entries[-seq_len(joint_categorical_pair_simple)],
    function(entry) simple[[entry$hi]] - simple[[entry$lo]]
  )

  unname(c(mu, unlist(simple), unlist(interaction)))
}

test_that("ipw() over two categorical treatments reports the whole crossing", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_pair()
  expect_true(all(table(dat$z1, dat$z2) > 0))

  two <- fit_joint_categorical_pair_route(dat)
  expect_identical(
    two$models$exposure_type,
    c(z1 = "categorical", z2 = "categorical")
  )

  res <- ipw(two$models, two$outcome_mod)
  est <- res$estimates

  # Nine cell means, twelve simple effects on two scales, four interactions on
  # two scales. The surface grows in both level counts at once and the identity
  # columns are what index it, so the rows are written out rather than counted.
  expected <- joint_categorical_pair_expected_rows()
  expect_identical(nrow(est), 41L)
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)
  expect_false(any(est$effect == "log(or)"))

  labels <- joint_categorical_labels(est)
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(names(coef(res)), labels)

  # Every reported quantity is g-computation on the weighted outcome model, the
  # interactions differencing two simple effects on their own scale.
  mu <- joint_categorical_pair_cell_means(two$outcome_mod, dat)
  expect_equal(
    est$estimate,
    joint_categorical_pair_expected_estimates(mu),
    tolerance = 1e-8
  )

  # Six coefficients for `z1 ~ x1 + x2` over two non-reference levels, fourteen
  # for `z2 ~ z1 * x1 + x2` over two more, ten for the outcome model, nine cell
  # means, and thirty-two contrast rows.
  expect_identical(as.integer(res$fit@n_params), 71L)
  expect_true(all(is.finite(est$std.err)))
  expect_true(all(est$std.err > 0))
})

test_that("the two routes agree on a three by three crossing", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_pair()

  two <- fit_joint_categorical_pair_route(dat)
  joint <- fit_joint_categorical_pair_exposure_route(dat)
  res_two <- ipw(two$models, two$outcome_mod)
  res_joint <- ipw(joint$ps_mod, joint$outcome_mod)

  # The whole surface, row for row and label for label: two multinomial models
  # of the two treatments and one over the nine cells report the same crossing
  # under the same names.
  expect_identical(nrow(res_joint$estimates), 41L)
  expect_identical(res_two$estimates$effect, res_joint$estimates$effect)
  expect_identical(res_two$estimates$contrast, res_joint$estimates$contrast)
  expect_identical(res_two$estimates$group, res_joint$estimates$group)
  expect_identical(
    joint_categorical_labels(res_two$estimates),
    joint_categorical_labels(res_joint$estimates)
  )

  # Two parameterizations of the same propensity score, so what separates them
  # is the difference between the models rather than solver noise. Measured on
  # this fixture the reported rows differ by at most 2.3e-2 in absolute terms,
  # on cell risks of 0.34 to 0.84 and contrasts up to 0.51, while the weights
  # behind them differ by up to 32 percent. The bounds are set from those
  # measurements with room.
  expect_lt(
    max(abs(res_two$estimates$estimate - res_joint$estimates$estimate)),
    0.06
  )

  # The cell risks are probabilities and carry a tighter bound; measured 8.3e-3.
  means <- res_two$estimates$effect == "mean"
  expect_lt(
    max(abs(
      res_two$estimates$estimate[means] - res_joint$estimates$estimate[means]
    )),
    0.03
  )
})

# ---- the categorical treatment in the first slot ----------------------------

test_that("the categorical-first surface is the declared route's label for label", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_first()

  ps_z <- nnet::multinom(
    z ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_a <- glm(a ~ z * x1 + x2, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      categorical_component_weights(ps_z, dat$z),
      wt_ate(ps_a),
      exposure_type = c("categorical", "binary")
    )
  )
  outcome_mod <- glm(
    y ~ a * z + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  res_two <- ipw(joint_wt_models(z = ps_z, a = ps_a), outcome_mod)

  dat$joint <- causalgenerics::joint_exposure(z = dat$z, a = dat$a)
  ps_joint <- nnet::multinom(
    joint ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  wts_joint <- withr::with_options(
    list(propensity.quiet = TRUE),
    categorical_component_weights(ps_joint, dat$joint)
  )
  outcome_joint <- glm(
    y ~ joint + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts_joint,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  res_joint <- ipw(ps_joint, outcome_joint)

  # Which treatment is declared first decides which one the cells vary fastest
  # and which one the interaction rows are framed under. Both routes read the
  # crossing the same way round, so the surface is the same labels in the same
  # order rather than the same shape: the cells vary `z`, the simple effects run
  # `z` within each level of `a` before `a` within each level of `z`, and the
  # interactions are `z`'s contrasts against `a = 1 vs a = 0`.
  expected <- data.frame(
    effect = c(
      rep("mean", 6),
      rep(c("rd", "log(rr)"), times = 9)
    ),
    contrast = c(
      "z = lo, a = 0",
      "z = mid, a = 0",
      "z = hi, a = 0",
      "z = lo, a = 1",
      "z = mid, a = 1",
      "z = hi, a = 1",
      rep(
        c(
          "z: mid vs lo",
          "z: hi vs lo",
          "z: mid vs lo",
          "z: hi vs lo",
          "a: 1 vs 0",
          "a: 1 vs 0",
          "a: 1 vs 0",
          "z: mid vs lo",
          "z: hi vs lo"
        ),
        each = 2
      )
    ),
    group = c(
      rep("overall", 6),
      rep(
        c(
          "a = 0",
          "a = 0",
          "a = 1",
          "a = 1",
          "z = lo",
          "z = mid",
          "z = hi",
          "a = 1 vs a = 0",
          "a = 1 vs a = 0"
        ),
        each = 2
      )
    ),
    stringsAsFactors = FALSE
  )

  expect_identical(nrow(res_two$estimates), 24L)
  expect_identical(res_two$estimates$effect, expected$effect)
  expect_identical(res_two$estimates$contrast, expected$contrast)
  expect_identical(res_two$estimates$group, expected$group)
  expect_identical(
    joint_categorical_labels(res_two$estimates),
    joint_categorical_labels(res_joint$estimates)
  )

  # Measured on this fixture the rows differ by at most 1.1e-2, on cell risks of
  # 0.38 to 0.81 and contrasts up to 0.46, against weights differing by up to 18
  # percent; the cell risks alone by 2.9e-3.
  expect_lt(
    max(abs(res_two$estimates$estimate - res_joint$estimates$estimate)),
    0.05
  )
  means <- res_two$estimates$effect == "mean"
  expect_lt(
    max(abs(
      res_two$estimates$estimate[means] - res_joint$estimates$estimate[means]
    )),
    0.03
  )
})

# ---- a two-level categorical component --------------------------------------

# A binary treatment whose second treatment is a two-level factor, so the pair
# can be fit either way: `nnet::multinom()` types the second component
# categorical and `stats::glm()` types it binary, and the two fits are of the
# same model.
sim_joint_categorical_two_level <- function(seed = 7305, n = 1200) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * x2))
  e <- rbinom(
    n,
    1,
    plogis(-0.2 + 0.5 * x1 + 0.3 * x2 - 0.8 * a + 0.6 * a * x1)
  )
  y <- rbinom(
    n,
    1,
    plogis(-0.5 + 0.7 * a + 0.5 * e + 0.6 * x1 - 0.3 * x2 + 0.9 * a * e)
  )

  data.frame(x1, x2, a, e = factor(e, levels = c(0, 1)), y)
}

test_that("the flattened coefficients of a two-level multinomial keep the level's name", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_two_level()
  fit <- nnet::multinom(
    e ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  # At two levels `coef()` on a multinomial fit gives a bare vector rather than
  # the one-row matrix the rest of its level counts give, so the level the
  # coefficients belong to is nowhere in the object the flattening reads. The
  # name has to come from the fit's own levels instead, which is the branch
  # every route taking a categorical component from a fit shares.
  expect_false(is.matrix(stats::coef(fit)))
  expect_identical(
    names(ipw_multinom_coefs(fit)),
    c("1:(Intercept)", "1:a", "1:x1", "1:x2", "1:a:x1")
  )
  expect_equal(
    unname(ipw_multinom_coefs(fit)),
    unname(stats::coef(fit)),
    tolerance = 1e-12
  )
})

test_that("a two-level multinomial component is carried or refused in the package's own terms", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_two_level()

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e_multinom <- nnet::multinom(
    e ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_e_glm <- glm(
    e ~ a * x1 + x2,
    data = dat,
    family = binomial(),
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  # The premise: the two fits are the same model solved by two optimizers, and
  # they land on the same coefficients to 6.7e-8 and the same weights to 1.3e-6.
  # Whatever the route does with the multinomial fit, it is not being handed a
  # different propensity score.
  expect_equal(
    unname(stats::coef(ps_e_multinom)),
    unname(stats::coef(ps_e_glm)),
    tolerance = 1e-6
  )
  weights_multinom <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_e_multinom)
  )
  weights_glm <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_e_glm)
  )
  expect_equal(
    as.double(weights_multinom),
    as.double(weights_glm),
    tolerance = 1e-5
  )

  # The single-treatment route meets the same fit with a refusal that names what
  # is wrong and what to fit instead, so the pair below is measured against a
  # route that already has an answer for a two-level multinomial.
  outcome_alone <- glm(
    y ~ e + x1,
    data = dat,
    family = quasibinomial(),
    weights = weights_multinom,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  expect_error(
    ipw(ps_e_multinom, outcome_alone),
    class = "propensity_model_family_error"
  )

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      weights_multinom,
      exposure_type = c("binary", "categorical")
    )
  )
  outcome_mod <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  models <- joint_wt_models(a = ps_a, e = ps_e_multinom)
  expect_identical(models$exposure_type, c(a = "binary", e = "categorical"))

  # Two resolutions are open here and either would be legitimate: the pair takes
  # the multinomial block and reports the fourteen-row crossing the equivalent
  # `glm` pair reports, or the route refuses it toward `glm` the way the
  # single-treatment route above does. What is not legitimate is the third
  # thing: a multinomial score is the block the route hands a categorical
  # component, and a two-level component has no rows for it, so the pair reaches
  # the estimating equation and fails inside it. The error the user meets then
  # is about a `y` with too few categories and about fitting an outcome, neither
  # of which is a treatment model they wrote.
  cnd <- tryCatch(ipw(models, outcome_mod), error = identity)
  if (inherits(cnd, "error")) {
    expect_s3_class(cnd, "propensity_error")
  } else {
    expect_identical(nrow(cnd$estimates), 14L)
  }
})

# The pair of fits the two tests below both build: the same second treatment
# model written twice, once for each optimizer that solves it, beside one first
# treatment model they share. `wt_joint()` is left to read each component's own
# recorded type, which is what a caller who names none gets, and which is
# `"binary"` on both sides: `wt_ate()` records the type of the weight it built,
# and the weight it builds from a two-level multinomial is the binary one.
joint_categorical_two_level_pair <- function() {
  dat <- sim_joint_categorical_two_level()

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e_multinom <- nnet::multinom(
    e ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_e_glm <- glm(
    e ~ a * x1 + x2,
    data = dat,
    family = binomial(),
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  outcome_mod <- function(w) {
    glm(
      y ~ a * e + x1,
      data = dat,
      family = quasibinomial(),
      weights = w,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  }

  withr::with_options(list(propensity.quiet = TRUE), {
    w_a <- wt_ate(ps_a)

    list(
      multinom = ipw(
        joint_wt_models(a = ps_a, e = ps_e_multinom),
        outcome_mod(wt_joint(w_a, wt_ate(ps_e_multinom)))
      ),
      glm = ipw(
        joint_wt_models(a = ps_a, e = ps_e_glm),
        outcome_mod(wt_joint(w_a, wt_ate(ps_e_glm)))
      )
    )
  })
}

test_that("a two-level multinomial component crosses into the same surface a binary one does", {
  skip_if_not_installed("nnet")
  fits <- joint_categorical_two_level_pair()

  # Four cells and ten contrasts over them, which is the surface a pair of
  # binary treatments reports. A treatment observed at two levels has that
  # surface whichever function fit it, so the crossing is the same crossing:
  # same cells, same reference, same contrasts, and the labels that name them.
  expect_identical(nrow(fits$multinom$estimates), 14L)
  expect_identical(
    fits$multinom$estimates[c("effect", "contrast", "group")],
    fits$glm$estimates[c("effect", "contrast", "group")]
  )

  # The block the system wrote for the component, read where the layout is
  # visible. A multinomial block names its parameters for the level they belong
  # to, so a component stacked as the multinomial it is classed as would carry
  # `1:(Intercept)` where the binary block carries `(Intercept)`.
  expect_identical(
    names(stats::coef(fits$multinom$fit)),
    names(stats::coef(fits$glm$fit))
  )
})

test_that("a two-level multinomial component agrees with its glm twin", {
  skip_if_not_installed("nnet")
  fits <- joint_categorical_two_level_pair()

  # The two fits are one model solved by two optimizers, and the stacked system
  # re-solves each component's score from the coefficients it was seeded at, so
  # what the pair reports does not inherit the seed's own precision: the
  # estimates agree to 1.5e-13 relative, well below the 6.7e-8 the coefficients
  # themselves differ by. The standard errors are read from a numerically
  # differentiated bread evaluated at those roots, which is what leaves them
  # agreeing to 2.4e-7 rather than to the estimates' own precision.
  expect_equal(
    fits$multinom$estimates$estimate,
    fits$glm$estimates$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    fits$multinom$estimates$std.err,
    fits$glm$estimates$std.err,
    tolerance = 1e-5
  )
})

# ---- refusals ---------------------------------------------------------------

test_that("a dose is refused in the first slot beside a categorical treatment", {
  # The position is what the dose is refused in rather than the company it
  # keeps: the first factor of the factorization carries no density, whichever
  # discrete treatment the second one would be.
  cnd_first <- tryCatch(
    check_ipw_joint_models_types(c(d = "continuous", z = "categorical")),
    error = identity
  )
  expect_s3_class(cnd_first, "propensity_ipw_exposure_error")
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(cnd_first)),
    "first",
    fixed = TRUE
  )

  # Two categorical treatments are the pair this refusal is not about.
  expect_silent(
    check_ipw_joint_models_types(c(z1 = "categorical", z2 = "categorical"))
  )
})

# ---- separation -------------------------------------------------------------

test_that("the joint preflight counts a categorical component's separation at the assigned level", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  two <- fit_joint_categorical_route(dat)

  spec <- ipw_spec_joint_models(two$models, two$outcome_mod)
  layout <- ipw_theta_layout(spec)

  # A fitted multinomial cannot reach the state this guard is for. Separation
  # drives the probability of the level a unit took toward one, not toward zero:
  # the covariate pattern that predicts a level without error predicts it for
  # the units that took it, and the columns that underflow belong to levels
  # those units did not take. Measured on a quasi-separated fit whose separating
  # coefficient reaches 19, the smallest assigned-level probability is 0.26 and
  # no entry of the softmax is an exact zero or one. What can reach the corner
  # is a theta the solver visits on its way, so the guard is evaluated where the
  # solver would evaluate it: at an init whose categorical block is scaled up
  # until the assigned level underflows.
  n_a <- length(coef(two$ps_a))
  categorical_block <- layout$idx$ps[n_a + seq_len(10)]

  # Scaled by two hundred, 1069 entries of the softmax are an exact zero or one
  # and not one of them is at a level its unit took. The binary rule would
  # refuse this fit; the assigned-level rule is what makes it pass.
  near <- layout
  near$init[categorical_block] <- layout$init[categorical_block] * 200
  blocks <- ipw_joint_models_blocks(
    spec,
    near$init[layout$idx$ps],
    near$init[layout$idx$stab]
  )
  saturated <- blocks[[2]]$ps
  expect_gt(sum(saturated == 0 | saturated == 1), 1000)
  expect_identical(sum(rowSums(spec$exposure[[2]] * saturated) == 0), 0L)
  expect_no_error(ipw_weights_at_init(spec, near))

  # Scaled by a thousand, 159 units have an exact zero at the level they took,
  # and those are the units whose weight has no denominator.
  over <- layout
  over$init[categorical_block] <- layout$init[categorical_block] * 1000
  blocks_over <- ipw_joint_models_blocks(
    spec,
    over$init[layout$idx$ps],
    over$init[layout$idx$stab]
  )
  assigned_zero <- sum(
    rowSums(spec$exposure[[2]] * blocks_over[[2]]$ps) == 0
  )
  expect_identical(assigned_zero, 159L)

  cnd <- tryCatch(ipw_weights_at_init(spec, over), error = identity)
  expect_s3_class(cnd, "propensity_ipw_separation_error")

  # The count the refusal reports is the assigned-level count rather than every
  # saturated entry of the softmax, which is the whole of the difference between
  # the rule this branch applies and the one a binary component takes.
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(
    message,
    paste("for", assigned_zero, "observations"),
    fixed = TRUE
  )
  expect_false(grepl(
    as.character(sum(
      blocks_over[[2]]$ps == 0 | blocks_over[[2]]$ps == 1
    )),
    message,
    fixed = TRUE
  ))

  # The binary component is left at its own seed throughout, so nothing it
  # carries is what fires.
  expect_identical(
    sum(blocks_over[[1]]$ps == 0 | blocks_over[[1]]$ps == 1),
    0L
  )
})

# ---- the stabilization slot -------------------------------------------------
#
# A categorical component needs no stabilization, and the sections above are
# written without it. A caller may still ask for it, and the three numerators a
# categorical component can carry are the three the single-treatment categorical
# route carries: the marginal probability of each level, a score the caller
# computed, and a fitted multinomial model of the treatment. What each is worth
# in the stacked system is what these tests pin.
#
# The arrangement is the binary component's, read for a treatment with more than
# two levels. A marginal numerator is k - 1 free proportions where a binary one
# is the single proportion `pi`; a score is a known multiplier and widens the
# system by nothing at all; and a fitted model contributes its own coefficients,
# which for a multinomial is one run of terms per non-reference level rather
# than one run in all. The names those parameters take compose the two
# conventions that already exist: the joint route's `stab_<component>_`, which
# says whose numerator a parameter belongs to, and the categorical block's
# `<level>:<term>`, which says which level's equation it sits in.
#
# The numerator here reads `x2`, which the propensity score model also reads and
# which the marginal structural model is given so the fixtures are silent about
# coverage. The one test that is about coverage takes it away again.

# The names a categorical component's numerator model fills its slice with. The
# joint convention's prefix composed with the categorical block's own layout:
# level-major and term-minor, which is the order `as.vector(t(coef()))` puts the
# coefficients in and the order the denominator's block is named in.
joint_categorical_stab_names <- function(component, model) {
  coefs <- stats::coef(model)
  terms <- colnames(coefs)

  as.vector(vapply(
    rownames(coefs),
    function(level) paste0("stab_", component, "_", level, ":", terms),
    character(length(terms))
  ))
}

# The two-model route with the categorical component stabilized. Only the
# numerator moves between the three fixtures: the binary component is left
# unstabilized throughout and the two propensity score models are the ones the
# sections above fit, so a difference between any two of these fits is a
# difference between numerators rather than between routes.
#
# The numerator model is fit inline against the frame this function holds, for
# the reason every multinomial in this file is: the design behind it is
# recovered by re-evaluating the fitting call rather than read off a model frame
# it does not keep.
fit_joint_categorical_stab_route <- function(
  dat,
  numerator = c("marginal", "score", "model"),
  outcome_rhs = "a * z + x1 + x2"
) {
  numerator <- match.arg(numerator)

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_z <- nnet::multinom(
    z ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  # A proper subset of the propensity score model's design, which is what keeps
  # the block the numerator contributes from being the denominator's block under
  # another name.
  num_z <- nnet::multinom(
    z ~ x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  probabilities <- unname(stats::predict(ps_z, type = "probs"))
  colnames(probabilities) <- ps_z$lev
  taken <- cbind(
    seq_len(nrow(dat)),
    match(as.character(dat$z), ps_z$lev)
  )

  fitted_numerator <- stats::fitted(num_z)
  score <- fitted_numerator[cbind(
    seq_len(nrow(dat)),
    match(as.character(dat$z), colnames(fitted_numerator))
  )]

  w_a <- withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_a))
  w_z <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      probabilities,
      dat$z,
      exposure_type = "categorical",
      stabilize = if (identical(numerator, "model")) num_z else TRUE,
      stabilization_score = if (identical(numerator, "score")) score
    )
  )
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(w_a, w_z, exposure_type = c("binary", "categorical"))
  )
  outcome_mod <- glm(
    stats::reformulate(outcome_rhs, response = "y"),
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(
    models = joint_wt_models(a = ps_a, z = ps_z),
    ps_a = ps_a,
    ps_z = ps_z,
    num_z = num_z,
    w_a = w_a,
    # The probability the component's own model gives the level each unit took,
    # which is the denominator of its factor of the product.
    denominator = probabilities[taken],
    score = score,
    outcome_mod = outcome_mod,
    wts = wts
  )
}

# The marginal probability of the level each unit took, which is what the
# default stabilizer's numerator evaluates to.
joint_categorical_marginal_numerator <- function(exposure) {
  proportions <- as.numeric(prop.table(table(exposure)))

  proportions[match(as.character(exposure), levels(exposure))]
}

test_that("ipw() estimates a categorical component stabilized on the marginal proportions", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  fx <- fit_joint_categorical_stab_route(dat, numerator = "marginal")

  # The product written out, so the rebuild below is held to the arithmetic
  # rather than to the preflight's own tolerance.
  expect_equal(
    as.double(fx$wts),
    as.double(fx$w_a) *
      (joint_categorical_marginal_numerator(dat$z) / fx$denominator),
    tolerance = 1e-12
  )

  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  layout <- ipw_theta_layout(spec)

  # A binary component's marginal numerator is one proportion and a categorical
  # one's is k - 1 of them, the reference level's completing the set, so the
  # slice the component takes is as wide as its treatment has non-reference
  # levels. The unstabilized binary component takes none.
  expect_identical(spec$stab$widths, c(0L, 2L))

  # The block is seeded at the level proportions of the fit's own rows, which is
  # where the weights were built, and it is named for the component and the
  # level rather than for the role alone: two components mean a name saying only
  # that a parameter is a proportion would not say whose.
  indicator <- ipw_joint_models_fit_exposure(spec, 2L)
  expect_identical(dim(indicator), c(nrow(dat), 3L))

  seed <- ipw_init_joint_models_stab(spec)
  expect_identical(names(seed), c("stab_z_mid", "stab_z_hi"))
  expect_equal(
    unname(seed),
    unname(colMeans(indicator)[-1]),
    tolerance = 1e-12
  )

  # The seed rebuilds the weights the outcome model was fit under, and every
  # equation sits at its root there, which is what makes the solve report the
  # numerator the weights were built from rather than one the solver reached.
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(fx$wts),
    tolerance = 1e-10
  )
  mat <- build_ipw_psi(spec, layout)(layout$init)
  expect_false(anyNA(mat))
  expect_true(all(abs(rowSums(mat)) / spec$n < 1e-8))

  res <- ipw(fx$models, fx$outcome_mod)
  theta <- coef(res$fit)

  expect_s3_class(res, "ipw")
  expect_identical(
    names(theta)[grepl("^stab_", names(theta))],
    c("stab_z_mid", "stab_z_hi")
  )

  # The solved proportions are the proportions of the rows the system analyzes,
  # which is what the equations written for them say: one row per non-reference
  # level, the indicator of that level less the parameter.
  expect_equal(
    unname(theta[c("stab_z_mid", "stab_z_hi")]),
    c(mean(dat$z == "mid"), mean(dat$z == "hi")),
    tolerance = 1e-8
  )

  # Two parameters more than the unstabilized pair carries over the same
  # models: three coefficients for `a ~ x1 + x2`, ten for `z ~ a * x1 + x2`, two
  # marginal proportions, eight for the outcome model, six cell means, and
  # eighteen contrast rows.
  expect_identical(as.integer(res$fit@n_params), 47L)

  # Every reported quantity is g-computation on the weighted outcome model, as
  # it is without a numerator: stabilization changes the weights and so the fit,
  # and changes nothing about what the surface reports off it.
  expect_identical(nrow(res$estimates), 24L)
  expect_equal(
    res$estimates$estimate,
    joint_categorical_expected_estimates(
      joint_categorical_cell_means(fx$outcome_mod, dat)
    ),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res$estimates$std.err)))
  expect_true(all(res$estimates$std.err > 0))
})

test_that("ipw() rebuilds a categorical component stabilized on a fixed score", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  fx <- fit_joint_categorical_stab_route(dat, numerator = "score")

  # The product written out: the score multiplies the component's factor as it
  # stands, so the weights are the unstabilized ones times a vector the caller
  # computed.
  expect_equal(
    as.double(fx$wts),
    as.double(fx$w_a) * (fx$score / fx$denominator),
    tolerance = 1e-12
  )

  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  layout <- ipw_theta_layout(spec)

  # A score is the numerator itself rather than a description of one, so the
  # system multiplies by it and estimates nothing, exactly as a binary
  # component's score is treated: no block, no parameters, no width.
  expect_identical(spec$stab$widths, c(0L, 0L))
  expect_identical(ipw_init_joint_models_stab(spec), numeric(0))

  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(fx$wts),
    tolerance = 1e-12
  )
  mat <- build_ipw_psi(spec, layout)(layout$init)
  expect_false(anyNA(mat))
  expect_true(all(abs(rowSums(mat)) / spec$n < 1e-8))

  res <- ipw(fx$models, fx$outcome_mod)

  expect_s3_class(res, "ipw")
  expect_false(any(grepl("^stab_", names(coef(res$fit)))))
  expect_identical(as.integer(res$fit@n_params), 45L)

  expect_identical(nrow(res$estimates), 24L)
  expect_equal(
    res$estimates$estimate,
    joint_categorical_expected_estimates(
      joint_categorical_cell_means(fx$outcome_mod, dat)
    ),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res$estimates$std.err)))
  expect_true(all(res$estimates$std.err > 0))
})

test_that("ipw() stacks a categorical component's numerator model", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  fx <- fit_joint_categorical_stab_route(dat, numerator = "model")

  # A model and the probabilities it evaluates to build one and the same set of
  # weights, which is what makes everything below a statement about the block
  # rather than about the fixture.
  expect_equal(
    as.double(fx$wts),
    as.double(fx$w_a) * (fx$score / fx$denominator),
    tolerance = 1e-12
  )

  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod)
  layout <- ipw_theta_layout(spec)

  # A multinomial numerator contributes one parameter per term per
  # non-reference level, so the slice is p times k - 1 rather than p: two terms
  # for `z ~ x2` over two non-reference levels.
  expect_identical(spec$stab$widths, c(0L, 4L))

  stab_names <- joint_categorical_stab_names("z", fx$num_z)
  expect_identical(
    stab_names,
    c(
      "stab_z_mid:(Intercept)",
      "stab_z_mid:x2",
      "stab_z_hi:(Intercept)",
      "stab_z_hi:x2"
    )
  )

  seed <- ipw_init_joint_models_stab(spec)
  expect_identical(names(seed), stab_names)
  expect_equal(
    unname(seed),
    as.vector(t(coef(fx$num_z))),
    tolerance = 1e-12
  )

  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(fx$wts),
    tolerance = 1e-10
  )
  mat <- build_ipw_psi(spec, layout)(layout$init)
  expect_false(anyNA(mat))
  # The numerator block is seeded at the multinom fit's own coefficients, and
  # nnet's optimizer leaves a mean score near 2e-8 there whatever reltol asks
  # for, so the fit's convergence floor bounds this pin rather than 1e-8.
  expect_true(all(abs(rowSums(mat)) / spec$n < 1e-7))

  res <- ipw(fx$models, fx$outcome_mod)
  theta <- coef(res$fit)

  expect_s3_class(res, "ipw")
  expect_identical(names(theta)[grepl("^stab_", names(theta))], stab_names)

  # The block is solved at the model the caller fit rather than carried at
  # whatever it was seeded with, which is what makes the reported standard
  # errors account for having estimated it. The layout is the denominator's:
  # one run of terms per non-reference level, in the fit's own level order.
  expect_equal(
    unname(theta[stab_names]),
    as.vector(t(coef(fx$num_z))),
    tolerance = 1e-6
  )

  expect_identical(nrow(res$estimates), 24L)
  expect_equal(
    res$estimates$estimate,
    joint_categorical_expected_estimates(
      joint_categorical_cell_means(fx$outcome_mod, dat)
    ),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(res$estimates$std.err)))
  expect_true(all(res$estimates$std.err > 0))
})

test_that("a categorical component's numerator model reports the score's estimates and other standard errors", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  model <- fit_joint_categorical_stab_route(dat, numerator = "model")
  score <- fit_joint_categorical_stab_route(dat, numerator = "score")

  # The two fixtures weight the units identically, which is what makes the
  # standard errors comparable at all.
  expect_equal(
    as.double(model$wts),
    as.double(score$wts),
    tolerance = 1e-12
  )

  res_model <- ipw(model$models, model$outcome_mod)
  res_score <- ipw(score$models, score$outcome_mod)

  expect_equal(
    res_model$estimates$estimate,
    res_score$estimates$estimate,
    tolerance = 1e-8
  )

  # A supplied score is a constant of the stacked system and a supplied model is
  # not, so the difference between the two is the numerator having been
  # estimated rather than known. The marginal structural model reads `z`, `x2`,
  # and `x1` without the interaction of the first two, so it is not saturated in
  # the numerator's cells and the reading reaches the reported rows.
  expect_true(all(is.finite(res_model$estimates$std.err)))
  expect_false(isTRUE(all.equal(
    res_model$estimates$std.err,
    res_score$estimates$std.err,
    tolerance = 1e-6
  )))
})

test_that("a categorical component's numerator model is asked about coverage", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()

  # A numerator conditioning on a variable the marginal structural model does
  # not read leaves that variable predicting the treatment in the weighted
  # sample, so the fit reports the effect within its levels rather than the
  # marginal one. The report is type-agnostic and a categorical component's
  # numerator model is a numerator model, so it is asked the same question.
  omitted <- fit_joint_categorical_stab_route(
    dat,
    numerator = "model",
    outcome_rhs = "a * z + x1"
  )
  cnd <- expect_warning(
    res <- ipw(omitted$models, omitted$outcome_mod),
    class = "propensity_ipw_stabilizer_coverage_warning"
  )
  expect_match(gsub("[[:space:]]+", " ", conditionMessage(cnd)), "x2")
  expect_s3_class(res, "ipw")
  expect_true(all(is.finite(res$estimates$estimate)))

  # The same numerator beside a marginal structural model that reads `x2` has
  # nothing to be missing, which is the fixture every other test here uses.
  covered <- fit_joint_categorical_stab_route(dat, numerator = "model")
  expect_no_warning(
    ipw(covered$models, covered$outcome_mod),
    class = "propensity_ipw_stabilizer_coverage_warning"
  )
})

test_that("the joint record keeps a categorical component's numerator model", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  fx <- fit_joint_categorical_stab_route(dat, numerator = "model")
  meta <- joint_wt_meta(fx$wts)

  # The estimator reads the model off the product rather than being handed it,
  # so what the product records is what it stacks. A component whose numerator
  # is the marginal proportion records no model, and the slot is per component
  # rather than one slot the pair shares.
  expect_identical(meta$stabilized, c(FALSE, TRUE))
  models <- joint_wt_numerator_models(meta)
  expect_length(models, 2L)
  expect_null(models[[1]])
  expect_identical(models[[2]], fx$num_z)

  marginal <- fit_joint_categorical_stab_route(dat, numerator = "marginal")
  expect_identical(
    joint_wt_numerator_models(joint_wt_meta(marginal$wts)),
    list(NULL, NULL)
  )
})

test_that("a stabilized categorical component is no longer refused by kind", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()

  # The three numerators a categorical component can carry are the three a
  # binary one carries, and the stacked system now writes a block for each: the
  # k - 1 free proportions, the score it multiplies by and estimates nothing
  # for, and the multinomial coefficients of a model the caller fit. None of
  # them is refused, so nothing about a stabilized categorical component is left
  # unsupported.
  for (numerator in c("marginal", "score", "model")) {
    fx <- fit_joint_categorical_stab_route(dat, numerator = numerator)
    res <- tryCatch(ipw(fx$models, fx$outcome_mod), error = identity)

    expect_false(
      inherits(res, "propensity_ipw_joint_models_stabilize_error"),
      label = numerator
    )
    expect_s3_class(res, "ipw")
  }
})

# ---- what a stacked categorical numerator model has to be -------------------

test_that("a categorical component's numerator fit to the levels in another order is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_z <- nnet::multinom(
    z ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  # Everything the block does is positional: the coefficients are flattened in
  # the fit's level order, the indicator matrix the score reads is laid out in
  # the component's own, and the softmax the numerator is rebuilt from takes the
  # first level as the reference. A numerator fit to the same levels in another
  # order is a different parameterization rather than a permutation of this one,
  # and reordering its columns cannot repair it. The weight side gathers by name
  # and takes such a fit, so the refusal is the estimator's to make.
  relevelled <- dat
  relevelled$z <- factor(as.character(dat$z), levels = c("hi", "lo", "mid"))
  num_z <- nnet::multinom(
    z ~ x2,
    data = relevelled,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  expect_identical(num_z$lev, c("hi", "lo", "mid"))
  expect_identical(ps_z$lev, c("lo", "mid", "hi"))

  probabilities <- unname(stats::predict(ps_z, type = "probs"))
  colnames(probabilities) <- ps_z$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      wt_ate(
        probabilities,
        dat$z,
        exposure_type = "categorical",
        stabilize = num_z
      ),
      exposure_type = c("binary", "categorical")
    )
  )
  outcome_mod <- glm(
    y ~ a * z + x1 + x2,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  # The refusal is the single-treatment categorical route's, read for a
  # component: the same class and the same cause, since what is wrong with the
  # fit does not depend on how many treatments were intervened on.
  cnd <- tryCatch(
    ipw(joint_wt_models(a = ps_a, z = ps_z), outcome_mod),
    error = identity
  )
  expect_s3_class(cnd, "propensity_ipw_numerator_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "level order", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})

test_that("a categorical component's numerator model of another response is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()

  # A second three-level factor over the same levels the treatment carries.
  # Nothing models it and it is not a treatment, so it is what a numerator of
  # the wrong response is fit to.
  dat$g <- factor(
    c("lo", "mid", "hi")[1 + (rank(dat$x1) %% 3)],
    levels = c("lo", "mid", "hi")
  )

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_z <- nnet::multinom(
    z ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  num_g <- nnet::multinom(
    g ~ x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  # The block's equations are the score of a model of the component's own
  # treatment, so a model of something else sits away from the root of the rows
  # seeded for it and the solve would move it, reporting a numerator nobody fit.
  # The levels agree, so the weight side gathers a probability for every unit
  # and has nothing to object to; what is wrong is which variable those
  # probabilities describe.
  probabilities <- unname(stats::predict(ps_z, type = "probs"))
  colnames(probabilities) <- ps_z$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      wt_ate(
        probabilities,
        dat$z,
        exposure_type = "categorical",
        stabilize = num_g
      ),
      exposure_type = c("binary", "categorical")
    )
  )
  outcome_mod <- glm(
    y ~ a * z + x1 + x2,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  cnd <- tryCatch(
    ipw(joint_wt_models(a = ps_a, z = ps_z), outcome_mod),
    error = identity
  )
  expect_s3_class(cnd, "propensity_ipw_numerator_error")

  # Both responses are named, so the message says which fit was handed over and
  # which treatment the component it was handed for models.
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "\"g\"", fixed = TRUE)
  expect_match(message, "\"z\"", fixed = TRUE)
})

# The refusals above, read for what they say about which numerator was refused.
# Both components' numerators arrive in the same argument, one `stabilize` per
# component, and `wt_mod` here is the container the two treatment models arrived
# in rather than a model anyone could refit, so a message naming those two
# arguments alone leaves a caller with two stabilized components nothing to act
# on. The route's other refusals name the component by the name its treatment
# model was given, and these name it the same way.
#
# The scenario is the one the refusals above are raised from: the same fits, the
# same weights, and a numerator of the caller's handed to the categorical
# component.
joint_categorical_numerator_refusal <- function(dat, ps_a, ps_z, numerator) {
  probabilities <- unname(stats::predict(ps_z, type = "probs"))
  colnames(probabilities) <- ps_z$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      wt_ate(
        probabilities,
        dat$z,
        exposure_type = "categorical",
        stabilize = numerator
      ),
      exposure_type = c("binary", "categorical")
    )
  )
  outcome_mod <- glm(
    y ~ a * z + x1 + x2,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  tryCatch(
    ipw(joint_wt_models(a = ps_a, z = ps_z), outcome_mod),
    error = identity
  )
}

test_that("the level order refusal names the component the numerator was built for", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_z <- nnet::multinom(
    z ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  relevelled <- dat
  relevelled$z <- factor(as.character(dat$z), levels = c("hi", "lo", "mid"))
  num_z <- nnet::multinom(
    z ~ x2,
    data = relevelled,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  cnd <- joint_categorical_numerator_refusal(dat, ps_a, ps_z, num_z)
  expect_s3_class(cnd, "propensity_ipw_numerator_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`stabilize` for `z`", fixed = TRUE)
  # The order the numerator has to declare is the component's own, so the
  # component names it here too, and `wt_mod` names nothing a reader could act
  # on.
  expect_match(message, "the level order `z` declares them in", fixed = TRUE)
  expect_no_match(message, "wt_mod", fixed = TRUE)
})

test_that("the wrong-response refusal names the component the numerator was built for", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()
  dat$g <- factor(
    c("lo", "mid", "hi")[1 + (rank(dat$x1) %% 3)],
    levels = c("lo", "mid", "hi")
  )

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_z <- nnet::multinom(
    z ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  num_g <- nnet::multinom(
    g ~ x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  cnd <- joint_categorical_numerator_refusal(dat, ps_a, ps_z, num_g)
  expect_s3_class(cnd, "propensity_ipw_numerator_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`stabilize` for `z`", fixed = TRUE)
  expect_match(message, "`z` models \"z\"", fixed = TRUE)
  expect_no_match(message, "wt_mod", fixed = TRUE)
})

# ---- a stabilized two-level multinomial component ---------------------------
#
# A `nnet::multinom()` fit to two levels is a logistic regression solved by a
# different optimizer, and the route stacks it as one. What that resolution owes
# the stabilization slot is the binary numerator rather than the multinomial
# one: the marginal proportion where a caller asks for the default, and a
# binomial coefficient block where a caller fits a model. Both are measured
# against the pair the same treatment model fit with `stats::glm()` reports.

# The same pair of fits stabilized, once for each optimizer that solves the
# second treatment's model. `wt_joint()` is left to read each component's own
# recorded type, which is `"binary"` on both sides, and the numerator is the
# binary one either way.
joint_categorical_two_level_stabilized_pair <- function(
  numerator = c("marginal", "model")
) {
  numerator <- match.arg(numerator)
  dat <- sim_joint_categorical_two_level()

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e_multinom <- nnet::multinom(
    e ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_e_glm <- glm(
    e ~ a * x1 + x2,
    data = dat,
    family = binomial(),
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  # A covariate the marginal structural model reads, which is what a numerator
  # may condition on without changing what the fit reports.
  stabilize <- if (identical(numerator, "marginal")) {
    TRUE
  } else {
    glm(
      e ~ x1,
      data = dat,
      family = binomial(),
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  }

  fit_one <- function(ps_e) {
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_joint(wt_ate(ps_a), wt_ate(ps_e, stabilize = stabilize))
    )
    outcome_mod <- glm(
      y ~ a * e + x1,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )

    ipw(joint_wt_models(a = ps_a, e = ps_e), outcome_mod)
  }

  list(multinom = fit_one(ps_e_multinom), glm = fit_one(ps_e_glm))
}

test_that("a stabilized two-level multinomial component agrees with its glm twin", {
  skip_if_not_installed("nnet")

  for (numerator in c("marginal", "model")) {
    fits <- joint_categorical_two_level_stabilized_pair(numerator)

    # The same fourteen-row crossing under the same labels, the same block for
    # the numerator, and the same numbers: the two fits are one model solved by
    # two optimizers, and the stacked system re-solves each component's score
    # from the coefficients it was seeded at, so what the pair reports does not
    # inherit the seed's own precision. The tolerances are the ones the
    # unstabilized pair above is held to, the standard errors being read from a
    # numerically differentiated bread rather than solved for.
    expect_identical(nrow(fits$multinom$estimates), 14L, label = numerator)
    expect_identical(
      fits$multinom$estimates[c("effect", "contrast", "group")],
      fits$glm$estimates[c("effect", "contrast", "group")],
      label = numerator
    )
    expect_identical(
      names(stats::coef(fits$multinom$fit)),
      names(stats::coef(fits$glm$fit)),
      label = numerator
    )
    expect_equal(
      fits$multinom$estimates$estimate,
      fits$glm$estimates$estimate,
      tolerance = 1e-8,
      label = numerator
    )
    expect_equal(
      fits$multinom$estimates$std.err,
      fits$glm$estimates$std.err,
      tolerance = 1e-5,
      label = numerator
    )
  }
})

# ---- a two-level multinomial component beside a dose ------------------------

# The two-level fixture with a dose added, so the factor can be the first
# component of a factorization whose second is continuous. The dose depends on
# it, which is what the container requires of a second model, and the outcome is
# continuous, which is what the dose surface's coefficient rows are read off.
sim_joint_two_level_dose <- function(seed = 7307, n = 1200) {
  dat <- sim_joint_categorical_two_level(seed = seed, n = n)
  withr::local_seed(seed + 1L)

  e_num <- as.numeric(dat$e == "1")
  dat$d <- 0.4 + 0.5 * dat$x1 - 0.3 * dat$x2 - 0.8 * e_num + rnorm(n)
  dat$yd <- 1 + 0.7 * e_num + 0.5 * dat$d + 0.6 * e_num * dat$d + rnorm(n)

  dat
}

test_that("a two-level multinomial component sits in the first slot beside a dose", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_two_level_dose()

  ps_e_multinom <- nnet::multinom(
    e ~ x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_e_glm <- glm(
    e ~ x1 + x2,
    data = dat,
    family = binomial(),
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  ps_d <- lm(d ~ e + x1 + x2, data = dat)

  fit_one <- function(ps_e) {
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_joint(
        wt_ate(ps_e),
        wt_ate(
          as.double(stats::fitted(ps_d)),
          dat$d,
          exposure_type = "continuous",
          stabilize = TRUE
        ),
        exposure_type = c("binary", "continuous")
      )
    )
    outcome_mod <- lm(yd ~ e * d, data = dat, weights = wts)

    ipw(joint_wt_models(e = ps_e, d = ps_d), outcome_mod)
  }

  res_multinom <- fit_one(ps_e_multinom)
  res_glm <- fit_one(ps_e_glm)

  # A dose is reported beside a first treatment with two levels, and a
  # multinomial fit to two levels is such a treatment's model: the pair takes
  # the dose's own coefficient surface rather than the cells a discrete pair
  # reports, and takes it whichever optimizer fit the first component.
  expect_identical(nrow(res_multinom$estimates), 3L)
  expect_identical(
    res_multinom$estimates[c("effect", "contrast", "group")],
    res_glm$estimates[c("effect", "contrast", "group")]
  )
  expect_identical(res_multinom$estimates$effect, c("diff", "slope", "diff"))

  expect_equal(
    res_multinom$estimates$estimate,
    res_glm$estimates$estimate,
    tolerance = 1e-8
  )
  expect_equal(
    res_multinom$estimates$std.err,
    res_glm$estimates$std.err,
    tolerance = 1e-5
  )
})

# ---- the categorical joint route over the rows `.data` keeps ----------------
#
# `.data` is where the designs this route rebuilds come from, and the rows that
# frame leaves need not be the rows the treatment models were fit over: an
# outcome model reading a covariate with a hole in it analyzes fewer rows than
# the models beside it did, and the stacked system solves over those. What a
# categorical component owes that restriction is what the single-treatment
# routes owe it. Its denominator design, its counterfactual designs, and its
# numerator model's design are rebuilt from the frame rather than read off fits
# that kept other rows; its stabilization block is seeded at the moments the
# weights were built at, which are the fit's, and solved at the moments of the
# rows analyzed, which are not; and the columns those rebuilds read are asked
# for, typed, and levelled before any of them runs.
#
# Two readings of the same request run through every test here: supplying the
# frame the fits were given and supplying the rows `ipw()` restricts it to are
# the same request, so the two forms report the same thing.

# The three-level fixture with two more columns. `v` is read by nothing but a
# numerator model, so it is the column the guards below take away or supply as
# something else. `w` is read by nothing but the outcome model and is withheld
# at a hundred and twenty of the rows that took the `hi` level, so the rows the
# outcome model analyzes are fewer than the rows the treatment models were fit
# over, and the level proportions over the two row sets differ by percents
# rather than by rounding. Withholding it at random would move them by parts in
# a thousand, and a pin written against either reading would then pass against
# the other.
sim_joint_categorical_gap <- function(seed = 7311, n = 1500) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  v <- rnorm(n)
  a <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * x2))

  eta_mid <- -0.1 + 0.5 * x1 + 0.3 * x2 - 0.7 * a + 0.5 * a * x1 + 0.4 * v
  eta_hi <- -0.3 + 0.3 * x1 - 0.4 * x2 + 0.6 * a - 0.4 * a * x1 - 0.3 * v
  denom <- 1 + exp(eta_mid) + exp(eta_hi)
  z <- draw_categorical(
    cbind(1, exp(eta_mid), exp(eta_hi)) / denom,
    c("lo", "mid", "hi")
  )

  z_mid <- as.integer(z == "mid")
  z_hi <- as.integer(z == "hi")
  y <- rbinom(
    n,
    1,
    plogis(
      -0.4 +
        0.7 * a +
        0.5 * z_mid +
        0.3 * z_hi +
        0.6 * x1 -
        0.3 * x2 +
        0.8 * a * z_mid -
        0.5 * a * z_hi
    )
  )

  dat <- data.frame(x1, x2, v, a, z, y)
  dat$w <- rev(dat$x1)
  dat$w[which(dat$z == "hi")[seq_len(120)]] <- NA

  dat
}

# The two-model route over the gapped frame, the categorical component
# stabilized and the outcome model reading the gapped covariate so that it
# analyzes fewer rows than either treatment model was fit over.
#
# The numerator model reads `v`, which no other model here does, and it is not
# saturated in what it reads. A numerator saturated in a binary covariate leaves
# the equations of the levels the gap does not touch at the same root over both
# row sets, since their proportions within each cell of that covariate are what
# solves them and the gap does not move those; half the pins below would then
# pass whichever rows the block had been solved over.
fit_joint_categorical_gap_route <- function(
  dat,
  numerator = c("marginal", "model"),
  outcome_rhs = "a * z + x1 + x2 + w"
) {
  numerator <- match.arg(numerator)

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_z <- nnet::multinom(
    z ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  num_z <- nnet::multinom(
    z ~ v + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  probabilities <- unname(stats::predict(ps_z, type = "probs"))
  colnames(probabilities) <- ps_z$lev

  w_a <- withr::with_options(list(propensity.quiet = TRUE), wt_ate(ps_a))
  w_z <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      probabilities,
      dat$z,
      exposure_type = "categorical",
      stabilize = if (identical(numerator, "model")) num_z else TRUE
    )
  )
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(w_a, w_z, exposure_type = c("binary", "categorical"))
  )
  outcome_mod <- glm(
    stats::reformulate(outcome_rhs, response = "y"),
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(
    models = joint_wt_models(a = ps_a, z = ps_z),
    ps_a = ps_a,
    ps_z = ps_z,
    num_z = num_z,
    outcome_mod = outcome_mod,
    wts = wts
  )
}

# The proportion of each non-reference level over a set of rows, which is what
# both readings of a categorical component's marginal numerator are and what
# tells the seed from the solution.
joint_categorical_level_shares <- function(exposure) {
  vapply(
    levels(exposure)[-1],
    function(level) mean(exposure == level),
    numeric(1)
  )
}

test_that("a categorical component's marginal numerator is seeded at the fit's proportions and solved at the analyzed rows'", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_gap()
  kept <- !is.na(dat$w)
  fx <- fit_joint_categorical_gap_route(dat, numerator = "marginal")

  # The two readings the block could take, and the gap between them. Everything
  # below is a statement about which of the two a seed is and which a solution
  # is, so the fixture has to separate them by more than any tolerance here
  # tolerates.
  fitted_shares <- joint_categorical_level_shares(dat$z)
  analyzed_shares <- joint_categorical_level_shares(dat$z[kept])
  expect_gt(min(abs(analyzed_shares - fitted_shares)), 0.02)

  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod, .data = dat)
  layout <- ipw_theta_layout(spec)

  expect_identical(spec$n, sum(kept))
  expect_identical(spec$stab$widths, c(0L, 2L))

  # The exposure moment the seed reads is the fit's own indicator matrix, which
  # is as long as the frame the fit was given rather than as the rows the system
  # goes on to solve over.
  indicator <- ipw_joint_models_fit_exposure(spec, 2L)
  expect_identical(dim(indicator), c(nrow(dat), 3L))

  seed <- ipw_init_joint_models_stab(spec)
  expect_identical(names(seed), c("stab_z_mid", "stab_z_hi"))
  expect_equal(unname(seed), unname(fitted_shares), tolerance = 1e-12)

  # Those are the moments the weights were built at, so the seed rebuilds the
  # weights the outcome model was fit under over exactly the rows it kept. A
  # seed already at the analyzed rows' proportions would rebuild weights the
  # caller was never given, and this is where that would be reported.
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(fx$wts)[kept],
    tolerance = 1e-10
  )

  # The rows written for the proportions are not at their root there, and what
  # they are off by is exactly the gap between the two readings: the seed is
  # where the weights were built and the solve is what moves it to the rows
  # analyzed.
  mat <- build_ipw_psi(spec, layout)(layout$init)
  expect_false(anyNA(mat))
  stab_rows <- match(c("stab_z_mid", "stab_z_hi"), names(layout$init))
  expect_equal(
    unname(rowMeans(mat)[stab_rows]),
    unname(analyzed_shares - fitted_shares),
    tolerance = 1e-10
  )

  res <- ipw(fx$models, fx$outcome_mod, .data = dat)
  theta <- coef(res$fit)

  expect_s3_class(res, "ipw")
  expect_identical(
    names(theta)[grepl("^stab_", names(theta))],
    c("stab_z_mid", "stab_z_hi")
  )

  # One row per non-reference level, the indicator of that level less the
  # parameter, written over the rows the spec analyzes and so solved by their
  # proportions rather than by the fit's.
  expect_equal(
    unname(theta[c("stab_z_mid", "stab_z_hi")]),
    unname(analyzed_shares),
    tolerance = 1e-8
  )

  # Supplying the frame the fits were given and supplying the rows `ipw()`
  # restricts it to are the same request.
  res_kept <- ipw(fx$models, fx$outcome_mod, .data = dat[kept, ])
  expect_equal(
    res$estimates$estimate,
    res_kept$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    res$estimates$std.err,
    res_kept$estimates$std.err,
    tolerance = 1e-10
  )
  expect_identical(nrow(res$estimates), 24L)
  expect_true(all(is.finite(res$estimates$std.err)))
  expect_true(all(res$estimates$std.err > 0))
})

test_that("a categorical component's numerator model is seeded at its own coefficients and solved at the analyzed rows' refit", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_gap()
  kept <- !is.na(dat$w)
  fx <- fit_joint_categorical_gap_route(dat, numerator = "model")

  # The block's rows are the multinomial score of the numerator's own design,
  # whose root over a set of rows is the refit on them. The numerator arrived
  # fit to the whole frame and the system reads it over the rows the outcome
  # model kept, so the two readings are the fit and the refit, and they differ
  # in every equation rather than in the one the gap sits in.
  refit <- nnet::multinom(
    z ~ v + x2,
    data = dat[kept, ],
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  fitted_coefs <- as.vector(t(coef(fx$num_z)))
  refit_coefs <- as.vector(t(coef(refit)))
  expect_false(isTRUE(all.equal(fitted_coefs, refit_coefs, tolerance = 1e-3)))

  spec <- ipw_spec_joint_models(fx$models, fx$outcome_mod, .data = dat)
  layout <- ipw_theta_layout(spec)

  expect_identical(spec$n, sum(kept))
  # Three terms over two non-reference levels, the design rebuilt from the
  # frame rather than recovered by re-evaluating the fitting call.
  expect_identical(spec$stab$widths, c(0L, 6L))
  expect_identical(
    dim(spec$stab$components[[2]]$model$X),
    c(sum(kept), 3L)
  )

  stab_names <- joint_categorical_stab_names("z", fx$num_z)
  seed <- ipw_init_joint_models_stab(spec)
  expect_identical(names(seed), stab_names)
  expect_equal(unname(seed), fitted_coefs, tolerance = 1e-12)

  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(fx$wts)[kept],
    tolerance = 1e-10
  )
  expect_false(anyNA(build_ipw_psi(spec, layout)(layout$init)))

  res <- muffle_coverage_warning(ipw(fx$models, fx$outcome_mod, .data = dat))
  theta <- coef(res$fit)

  expect_s3_class(res, "ipw")
  expect_identical(names(theta)[grepl("^stab_", names(theta))], stab_names)
  expect_equal(unname(theta[stab_names]), refit_coefs, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(
    unname(theta[stab_names]),
    fitted_coefs,
    tolerance = 1e-3
  )))

  res_kept <- muffle_coverage_warning(ipw(
    fx$models,
    fx$outcome_mod,
    .data = dat[kept, ]
  ))
  expect_equal(
    res$estimates$estimate,
    res_kept$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    res$estimates$std.err,
    res_kept$estimates$std.err,
    tolerance = 1e-10
  )
  expect_identical(nrow(res$estimates), 24L)
  expect_true(all(is.finite(res$estimates$std.err)))
  expect_true(all(res$estimates$std.err > 0))
})

test_that("a categorical joint fit reports the same thing whether or not .data is supplied", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical()

  # The frame every model was fit to, with nothing to restrict: the designs
  # rebuilt from it are the designs read off the fits, so the two calls are two
  # ways of asking for one thing. Measured on this fixture the four routes agree
  # to the last bit, which is what a rebuild that reproduces the fits' own
  # designs gives; the pins leave room for a platform whose model matrices are
  # assembled in another order.
  unstabilized <- fit_joint_categorical_route(dat)
  res <- ipw(unstabilized$models, unstabilized$outcome_mod)
  res_data <- ipw(
    unstabilized$models,
    unstabilized$outcome_mod,
    .data = dat
  )

  expect_identical(
    res_data$estimates[c("effect", "contrast", "group")],
    res$estimates[c("effect", "contrast", "group")]
  )
  expect_equal(
    res_data$estimates$estimate,
    res$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    res_data$estimates$std.err,
    res$estimates$std.err,
    tolerance = 1e-8
  )

  # And with each numerator a categorical component can carry, since each adds
  # a block of its own to the system and the model numerator adds a design to
  # rebuild alongside the treatments'.
  for (numerator in c("marginal", "score", "model")) {
    fx <- fit_joint_categorical_stab_route(dat, numerator = numerator)
    fit <- ipw(fx$models, fx$outcome_mod)
    fit_data <- ipw(fx$models, fx$outcome_mod, .data = dat)

    expect_identical(
      names(coef(fit_data$fit)),
      names(coef(fit$fit)),
      label = numerator
    )
    expect_equal(
      fit_data$estimates$estimate,
      fit$estimates$estimate,
      tolerance = 1e-10,
      label = numerator
    )
    expect_equal(
      fit_data$estimates$std.err,
      fit$estimates$std.err,
      tolerance = 1e-8,
      label = numerator
    )
  }
})

# ---- what `.data` has to supply a categorical component ---------------------

test_that("a categorical component's numerator covariate absent from .data is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_gap()
  fx <- fit_joint_categorical_gap_route(dat, numerator = "model")

  # The columns every rebuild reads are asked for before any of them runs, and
  # a categorical component's numerator model joins that sweep like the other
  # models: `v` is read by nothing else here, so leaving it out is a column the
  # numerator's design alone would have gone looking for.
  supplied <- dat
  supplied$v <- NULL

  err <- expect_error(
    muffle_coverage_warning(ipw(
      fx$models,
      fx$outcome_mod,
      .data = supplied
    )),
    class = "propensity_columns_exist_error"
  )
  expect_match(
    gsub("[[:space:]]+", " ", conditionMessage(err)),
    "v",
    fixed = TRUE
  )
})

test_that("a categorical component's numerator covariate supplied as a factor is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_gap()
  fx <- fit_joint_categorical_gap_route(dat, numerator = "model")

  # The type sweep compares the class the fit recorded for a column with the
  # class `.data` supplies and refuses the pair before any design is rebuilt.
  # What it heads off is a factor of three levels taking two design columns
  # where the number it stands in for took one, which the multinomial
  # coefficients would then be multiplied against positionally.
  supplied <- dat
  supplied$v <- factor(
    c("a", "b", "c")[1 + (rank(dat$x1) %% 3)],
    levels = c("a", "b", "c")
  )

  err <- expect_error(
    muffle_coverage_warning(ipw(
      fx$models,
      fx$outcome_mod,
      .data = supplied
    )),
    class = "propensity_ipw_data_error"
  )

  message <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(message, "v", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)
  expect_match(message, "stabilize", fixed = TRUE)
})

# The gapped frame with the numerator's covariate written two more ways: as a
# factor of three labels, which is what a level guard compares, and as a matrix
# of the indicators those labels imply, which is what no guard before the width
# count can compare. Both are read by the numerator model alone, so what happens
# to either is a statement about the fit that reads it.
joint_categorical_numerator_covariates <- function(dat) {
  dat$vf <- factor(
    c("p", "q", "r")[1 + (rank(dat$v) %% 3)],
    levels = c("p", "q", "r")
  )
  dat$vm <- stats::model.matrix(~vf, data = dat)[, -1, drop = FALSE]

  dat
}

# The two-model route over that frame, the categorical component stabilized on a
# numerator model reading whichever spelling of the covariate a test hands it.
fit_joint_categorical_numerator_covariate <- function(dat, covariate) {
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_z <- nnet::multinom(
    z ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  num_z <- nnet::multinom(
    stats::reformulate(c(covariate, "x2"), response = "z"),
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )

  probabilities <- unname(stats::predict(ps_z, type = "probs"))
  colnames(probabilities) <- ps_z$lev

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      wt_ate(
        probabilities,
        dat$z,
        exposure_type = "categorical",
        stabilize = num_z
      ),
      exposure_type = c("binary", "categorical")
    )
  )
  outcome_mod <- glm(
    y ~ a * z + x1 + x2 + w,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(models = joint_wt_models(a = ps_a, z = ps_z), outcome_mod = outcome_mod)
}

test_that("a categorical component's numerator covariate `.data` levels beyond the fit is refused", {
  skip_if_not_installed("nnet")
  dat <- joint_categorical_numerator_covariates(sim_joint_categorical_gap())
  fx <- fit_joint_categorical_numerator_covariate(dat, "vf")

  # A factor column the numerator alone reads is re-leveled against the levels
  # that fit recorded before its design is rebuilt, so a fourth label is refused
  # where the re-leveling happens rather than reaching the width count as a
  # design one column wider than the coefficients it multiplies.
  supplied <- dat
  supplied$vf <- factor(
    c("p", "q", "r", "s")[1 + (rank(dat$x1) %% 4)],
    levels = c("p", "q", "r", "s")
  )

  err <- expect_error(
    muffle_coverage_warning(ipw(fx$models, fx$outcome_mod, .data = supplied)),
    class = "propensity_ipw_data_error"
  )
  message <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(message, "vf", fixed = TRUE)
  expect_match(message, "s", fixed = TRUE)
})

test_that("a categorical component's numerator design rebuilt to another width is refused", {
  skip_if_not_installed("nnet")
  dat <- joint_categorical_numerator_covariates(sim_joint_categorical_gap())
  fx <- fit_joint_categorical_numerator_covariate(dat, "vm")

  # A matrix column is the shape both earlier guards pass over: the type sweep
  # compares the codings a class implies and has no noun for a matrix, and a
  # matrix carries no levels to be re-leveled against. What is left to notice
  # that `.data` rebuilt a wider design than the numerator was fit to is the
  # count of its columns, and the refusal names the argument the fit arrived in
  # so that a caller with two stabilized components knows which one to look at.
  supplied <- dat
  supplied$vm <- cbind(dat$vm, s = as.double(dat$x1 > 0))

  err <- expect_error(
    muffle_coverage_warning(ipw(fx$models, fx$outcome_mod, .data = supplied)),
    class = "propensity_ipw_data_error"
  )
  message <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(message, "stabilize", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)
  expect_match(message, "5 columns", fixed = TRUE)
  expect_match(message, "4 coefficients per equation", fixed = TRUE)

  # The same frame with the column it was fit to rebuilds the design and the fit
  # runs, so the refusal is about the width rather than about the shape.
  res <- muffle_coverage_warning(ipw(fx$models, fx$outcome_mod, .data = dat))
  expect_s3_class(res, "ipw")
})

test_that("a categorical treatment `.data` levels another way is resolved rather than refused", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_gap()
  fx <- fit_joint_categorical_gap_route(dat, numerator = "marginal")

  # A binary component's treatment column is read positionally, so the order it
  # declares is a fault the route names. A categorical one is not: the column is
  # levelled against the fit's own `lev` before anything reads it, so the order
  # it arrives in says nothing and a column carrying the same labels in another
  # order is the same treatment written down differently. Both codings a factor
  # and a character column give are resolved that way, and all three report the
  # same thing.
  relevelled <- dat
  relevelled$z <- factor(as.character(dat$z), levels = c("hi", "lo", "mid"))
  as_character <- dat
  as_character$z <- as.character(dat$z)

  res <- ipw(fx$models, fx$outcome_mod, .data = dat)
  res_relevelled <- ipw(fx$models, fx$outcome_mod, .data = relevelled)
  res_character <- ipw(fx$models, fx$outcome_mod, .data = as_character)

  for (other in list(res_relevelled, res_character)) {
    expect_identical(
      other$estimates[c("effect", "contrast", "group")],
      res$estimates[c("effect", "contrast", "group")]
    )
    expect_equal(
      other$estimates$estimate,
      res$estimates$estimate,
      tolerance = 1e-10
    )
    expect_equal(
      other$estimates$std.err,
      res$estimates$std.err,
      tolerance = 1e-10
    )
  }
})

test_that("a categorical treatment `.data` holds a level the model never saw is refused", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_gap()
  fx <- fit_joint_categorical_gap_route(dat, numerator = "marginal")

  # The other half of that resolution: a label the fit was never given has no
  # column of coefficients to be read against, so it is refused where the
  # resolution happens rather than reaching the softmax as a level the
  # denominator cannot score.
  #
  # The label is written into a row the outcome model kept. `.data` is
  # restricted to those rows before any column is read, so a value on a row the
  # restriction drops is a value nothing goes on to score and there is nothing
  # there to refuse.
  supplied <- dat
  supplied$z <- as.character(dat$z)
  supplied$z[[which(!is.na(dat$w))[[1]]]] <- "elsewhere"

  err <- expect_error(
    ipw(fx$models, fx$outcome_mod, .data = supplied),
    class = "propensity_ipw_exposure_error"
  )

  message <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(message, "elsewhere", fixed = TRUE)
  expect_match(message, "z", fixed = TRUE)
})

test_that("a categorical component's numerator model whose data is gone asks for .data", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_categorical_gap()
  dat <- dat[!is.na(dat$w), ]

  # A multinomial keeps no model frame and a `glm` fit with `model = FALSE`
  # keeps none either, so both rebuild one by re-evaluating their fitting call,
  # which a fit made inside a function whose frame is gone cannot do. Both
  # components meet that failure through one report, so both are refused in the
  # same terms: the class is the one that names `.data` as the remedy rather
  # than the numerator class the level-order and wrong-response guards raise.
  # Those two are about the fit being the wrong fit, which no data supplies a
  # way out of; this one is about the data behind a fit that is right, which is
  # exactly what `.data` carries. Pinning the pair together is what keeps a
  # categorical component's answer from drifting from its binary sibling's.
  categorical_fmla <- z ~ v
  binary_fmla <- a ~ v
  fit_categorical_in_function <- function(fitting_data) {
    nnet::multinom(
      categorical_fmla,
      data = fitting_data,
      trace = FALSE,
      reltol = 1e-14,
      maxit = 2000
    )
  }
  fit_binary_in_function <- function(fitting_data) {
    glm(
      binary_fmla,
      data = fitting_data,
      family = binomial(),
      model = FALSE
    )
  }
  gone_z <- fit_categorical_in_function(dat)
  gone_a <- fit_binary_in_function(dat)

  expect_error(model.frame(gone_z))
  expect_error(model.frame(gone_a))

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_z <- nnet::multinom(
    z ~ a * x1 + x2,
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  probabilities <- unname(stats::predict(ps_z, type = "probs"))
  colnames(probabilities) <- ps_z$lev

  fit_one <- function(a_stabilize, z_stabilize) {
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_joint(
        wt_ate(ps_a, stabilize = a_stabilize),
        wt_ate(
          probabilities,
          dat$z,
          exposure_type = "categorical",
          stabilize = z_stabilize
        ),
        exposure_type = c("binary", "categorical")
      )
    )

    glm(
      y ~ a * z + x1 + x2,
      data = dat,
      family = quasibinomial(),
      weights = wts,
      control = glm.control(epsilon = 1e-14, maxit = 200)
    )
  }

  outcome_z <- fit_one(FALSE, gone_z)
  outcome_a <- fit_one(gone_a, FALSE)
  models <- joint_wt_models(a = ps_a, z = ps_z)

  cnd_z <- tryCatch(
    muffle_coverage_warning(ipw(models, outcome_z)),
    error = identity
  )
  cnd_a <- tryCatch(
    muffle_coverage_warning(ipw(models, outcome_a)),
    error = identity
  )

  expect_s3_class(cnd_z, "propensity_ipw_data_error")
  expect_identical(class(cnd_z), class(cnd_a))

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd_z))
  expect_match(message, "stabilize", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)

  # And the remedy it names is one: the frame rebuilds the numerator's design
  # from the fit's own terms, and the fit runs.
  res <- muffle_coverage_warning(ipw(models, outcome_z, .data = dat))
  expect_s3_class(res, "ipw")
  expect_identical(nrow(res$estimates), 24L)
  expect_true(all(is.finite(res$estimates$std.err)))
})
