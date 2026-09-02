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
# Everything here is unstabilized. A categorical component needs no
# stabilization, and the numerator blocks a stabilized one would carry are
# pinned separately.

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

test_that("a stabilized categorical component is refused toward an unstabilized one", {
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
  probabilities <- unname(stats::predict(ps_z, type = "probs"))
  colnames(probabilities) <- ps_z$lev
  stabilized <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      probabilities,
      dat$z,
      exposure_type = "categorical",
      stabilize = TRUE
    )
  )
  expect_true(is_stabilized(stabilized))

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      stabilized,
      exposure_type = c("binary", "categorical")
    )
  )
  outcome_mod <- glm(
    y ~ a * z + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  models <- joint_wt_models(a = ps_a, z = ps_z)

  # The denominator of a categorical component is carried and its numerator is
  # not, so a component whose weights record one is refused rather than weighted
  # by a numerator the stack cannot differentiate. The refusal names the
  # component and the argument that rebuilds it.
  cnd <- tryCatch(ipw(models, outcome_mod), error = identity)
  expect_s3_class(cnd, "propensity_ipw_joint_models_stabilize_error")
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`z`", fixed = TRUE)
  expect_match(message, "stabilize = FALSE", fixed = TRUE)

  # The refusal is raised before anything is estimated, so it cannot be reached
  # through a fixture that also reports on its stabilizer's coverage.
  expect_no_warning(tryCatch(ipw(models, outcome_mod), error = function(e) {
    NULL
  }))
})

test_that("a dose is refused beside a categorical first treatment", {
  # The rows a dose reports are the marginal structural model's own coefficients
  # written for the two levels of a binary first treatment, so a categorical one
  # beside it has nowhere to be reported. Admitting a categorical treatment into
  # the first slot is what makes this pair reachable and the refusal necessary,
  # and the refusal says which of the two positions the dose may take.
  cnd <- tryCatch(
    check_ipw_joint_models_types(c(z = "categorical", d = "continuous")),
    error = identity
  )
  expect_s3_class(cnd, "propensity_ipw_exposure_error")
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`d`", fixed = TRUE)
  expect_match(message, "second", fixed = TRUE)
  expect_match(message, "only beside a binary first one", fixed = TRUE)

  # The other direction, which the dose has always been refused in: the first
  # factor of the factorization carries no density.
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
