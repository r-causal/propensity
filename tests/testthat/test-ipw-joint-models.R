# The two-model route to a joint intervention.
#
# `joint_exposure()` crosses two treatments into one categorical exposure and
# weights it with one multinomial propensity score model. The route here weights
# the same intervention through the sequential factorization it really has,
# f(A | L) f(E | A, L): two treatment models, one per treatment, and the product
# of their weights. `joint_wt_models()` records the pair and `wt_joint()` builds
# the product; what this file pins is `ipw()` over that container.
#
# The reported surface is the one the declared-exposure route reports and not a
# second surface that resembles it: the four counterfactual cell risks under the
# effect label "mean", the simple effects keyed by the level the other treatment
# is held at, and the interaction on each collapsible scale. Fourteen rows for a
# binary outcome, nine for a continuous one, labelled identically. The
# expectations below are written out against those conventions rather than
# imported from the joint-exposure file, which keeps this file readable on its
# own and makes a divergence between the two routes show up as a disagreement
# between two written-out contracts.
#
# The centrepiece is the equivalence: the two routes estimate the same estimand
# through different parameterizations of the same propensity score, so they
# agree loosely in general and to the last bits when the two parameterizations
# are saturated in the same covariate. Both gaps are measured on these fixtures
# and the tolerances are set from the measurements.

# ---- data simulators --------------------------------------------------------

# Two binary treatments sharing confounders, the second depending on the first.
# The dependence of E on A is not additive: the sequential factorization is only
# as good as the second model, and the guardrail on `joint_wt_models()` asks for
# that dependence to be modeled flexibly rather than as a single term, so the
# simulator gives it something to model.
sim_joint_models <- function(seed = 6401, n = 900) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, plogis(0.3 * x1 - 0.4 * x2))
  e <- rbinom(n, 1, plogis(-0.2 + 0.5 * x1 + 0.3 * x2 - 0.8 * a + 0.6 * a * x1))
  y <- rbinom(
    n,
    1,
    plogis(-0.5 + 0.7 * a + 0.5 * e + 0.6 * x1 - 0.3 * x2 + 0.9 * a * e)
  )
  yc <- 1 +
    0.6 * a +
    0.4 * e +
    0.5 * x1 -
    0.2 * x2 +
    0.8 * a * e +
    rnorm(n)
  data.frame(x1, x2, a, e, y, yc)
}

# One binary covariate and nothing else, so that a sequential pair saturated in
# it and a multinomial saturated in it describe the same joint distribution of
# the two treatments exactly rather than approximately. This is the fixture the
# tight equivalence pin runs on.
sim_joint_models_saturated <- function(seed = 6402, n = 1200) {
  withr::local_seed(seed)
  x2 <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, plogis(0.4 - 0.8 * x2))
  e <- rbinom(n, 1, plogis(-0.2 + 0.6 * x2 - 0.9 * a + 0.7 * a * x2))
  y <- rbinom(
    n,
    1,
    plogis(-0.5 + 0.7 * a + 0.5 * e - 0.3 * x2 + 0.9 * a * e)
  )
  data.frame(x2, a, e, y)
}

# ---- model fitting ----------------------------------------------------------

# The two-model route. `e_rhs` carries the second model's right-hand side so a
# test can vary its dimension, and `outcome_rhs` the outcome model's; the
# formulas are built here rather than passed in so the weights resolve against
# this frame when `glm()` looks for them.
fit_joint_models_route <- function(
  dat,
  a_rhs = c("x1", "x2"),
  e_rhs = "a * x1 + x2",
  outcome_rhs = "a * e + x1",
  outcome_family = "binomial"
) {
  ps_a <- glm(
    stats::reformulate(a_rhs, response = "a"),
    data = dat,
    family = binomial()
  )
  ps_e <- glm(
    stats::reformulate(e_rhs, response = "e"),
    data = dat,
    family = binomial()
  )
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      wt_ate(ps_e),
      exposure_type = c("binary", "binary")
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
    models = joint_wt_models(a = ps_a, e = ps_e),
    ps_a = ps_a,
    ps_e = ps_e,
    outcome_mod = outcome_mod,
    wts = wts
  )
}

# The declared-exposure route over the same data: the crossing as one factor,
# one multinomial propensity score model, categorical ate weights, and an
# outcome model over the cells. `y ~ joint + x1` and `y ~ a * e + x1` span the
# same columns, so given the same weights the two outcome fits coincide and the
# only thing separating the routes is the propensity score parameterization.
fit_joint_exposure_route <- function(
  dat,
  ps_rhs = c("x1", "x2"),
  outcome_rhs = "joint + x1"
) {
  dat$joint <- causalgenerics::joint_exposure(a = dat$a, e = dat$e)
  ps_mod <- nnet::multinom(
    stats::reformulate(ps_rhs, response = "joint"),
    data = dat,
    trace = FALSE,
    reltol = 1e-14,
    maxit = 2000
  )
  ps_named <- unname(predict(ps_mod, type = "probs"))
  colnames(ps_named) <- ps_mod$lev
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(ps_named, dat$joint, exposure_type = "categorical")
  )
  outcome_mod <- glm(
    stats::reformulate(outcome_rhs, response = "y"),
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(ps_mod = ps_mod, outcome_mod = outcome_mod, wts = wts, dat = dat)
}

# ---- the reported surface ---------------------------------------------------

# The cells, in the order the crossing varies them: the first treatment fastest,
# the reference cell first.
joint_models_cells <- c(
  "a = 0, e = 0",
  "a = 1, e = 0",
  "a = 0, e = 1",
  "a = 1, e = 1"
)

# The identity columns of the surface both routes report, written out against
# the conventions rather than read off either of them.
joint_models_expected_rows <- function(forms = c("rd", "log(rr)")) {
  simple <- list(
    c(contrast = "a: 1 vs 0", group = "e = 0"),
    c(contrast = "a: 1 vs 0", group = "e = 1"),
    c(contrast = "e: 1 vs 0", group = "a = 0"),
    c(contrast = "e: 1 vs 0", group = "a = 1")
  )
  interaction <- list(
    c(contrast = "a: 1 vs 0", group = "e = 1 vs e = 0")
  )
  entries <- c(simple, interaction)

  data.frame(
    effect = c(
      rep("mean", length(joint_models_cells)),
      rep(forms, times = length(entries))
    ),
    contrast = c(
      joint_models_cells,
      rep(
        vapply(entries, function(x) x[["contrast"]], ""),
        each = length(forms)
      )
    ),
    group = c(
      rep("overall", length(joint_models_cells)),
      rep(vapply(entries, function(x) x[["group"]], ""), each = length(forms))
    ),
    stringsAsFactors = FALSE
  )
}

joint_models_labels <- function(estimates) {
  paste(estimates$effect, estimates$contrast, estimates$group)
}

# The four counterfactual cell risks by g-computation on a weighted outcome
# model that reads the two treatments separately, standardized over the whole
# sample. The ate tilt is the constant one, so the averaging population is the
# whole sample and nothing is weighted here.
joint_models_cell_means <- function(outcome_mod, data) {
  values <- list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  out <- vapply(
    values,
    function(cell) {
      d <- data
      d$a <- cell[[1]]
      d$e <- cell[[2]]
      mean(predict(outcome_mod, newdata = d, type = "response"))
    },
    numeric(1)
  )
  stats::setNames(out, joint_models_cells)
}

joint_models_transform <- function(form, mu_hi, mu_lo) {
  switch(
    form,
    rd = mu_hi - mu_lo,
    diff = mu_hi - mu_lo,
    "log(rr)" = log(mu_hi) - log(mu_lo)
  )
}

# Every reported estimate, built from the four cell means in the row order
# above. `mu[[1]]` is the reference cell, `mu[[2]]` is a = 1 e = 0, `mu[[3]]` is
# a = 0 e = 1, and `mu[[4]]` is a = 1 e = 1.
joint_models_expected_estimates <- function(mu, forms = c("rd", "log(rr)")) {
  pairs <- list(c(2L, 1L), c(4L, 3L), c(3L, 1L), c(4L, 2L))
  simple <- unlist(lapply(pairs, function(ix) {
    vapply(
      forms,
      function(f) joint_models_transform(f, mu[[ix[[1]]]], mu[[ix[[2]]]]),
      numeric(1)
    )
  }))
  interaction <- vapply(
    forms,
    function(f) {
      joint_models_transform(f, mu[[4]], mu[[3]]) -
        joint_models_transform(f, mu[[2]], mu[[1]])
    },
    numeric(1)
  )

  unname(c(mu, simple, interaction))
}

# ---- baselines: the halves that stand today ---------------------------------

test_that("the fixture fits both routes and the constructors accept them", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_models()

  # Every cell of the crossing is populated, which the declared-exposure route
  # requires and which the two-model route needs a mean for.
  expect_true(all(table(dat$a, dat$e) > 0))

  two <- fit_joint_models_route(dat)
  expect_true(is_joint_wt_models(two$models))
  expect_identical(two$models$names, c("a", "e"))
  expect_identical(
    two$models$exposure_type,
    c(a = "binary", e = "binary")
  )

  # The weights are the product the container's models imply, and they say so.
  expect_true(is_joint_wt(two$wts))
  expect_identical(estimand(two$wts), "ate")
  expect_equal(
    as.double(two$wts),
    as.double(wt_ate(two$ps_a)) * as.double(wt_ate(two$ps_e)),
    tolerance = 1e-12
  )

  # The declared-exposure route is merged, so its half of the equivalence
  # comparison is a fact about the package today rather than an aspiration.
  joint <- fit_joint_exposure_route(dat)
  res <- ipw(joint$ps_mod, joint$outcome_mod)
  expected <- joint_models_expected_rows()
  expect_identical(res$estimates$effect, expected$effect)
  expect_identical(res$estimates$contrast, expected$contrast)
  expect_identical(res$estimates$group, expected$group)
  expect_identical(nrow(res$estimates), 14L)
})

test_that("the two parameterizations disagree about the weights they imply", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)
  joint <- fit_joint_exposure_route(dat)

  # What makes the loose equivalence below a statement about estimators rather
  # than about arithmetic: a sequential pair of logistic models and a
  # multinomial logit over the cells are different models of the same
  # conditional distribution, and on this fixture the weights they imply differ
  # by up to 32 percent. Two routes agreeing to a few percent in the estimates
  # after that is the claim being made.
  gap <- max(
    abs(as.double(two$wts) - as.double(joint$wts)) / as.double(joint$wts)
  )
  expect_gt(gap, 0.1)
})

test_that("saturated parameterizations imply the same weights", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_models_saturated()
  two <- fit_joint_models_route(
    dat,
    a_rhs = "x2",
    e_rhs = "a * x2",
    outcome_rhs = "a * e + x2"
  )
  joint <- fit_joint_exposure_route(
    dat,
    ps_rhs = "x2",
    outcome_rhs = "joint + x2"
  )

  # `a ~ x2` with `e ~ a * x2` and a multinomial over the cells on `x2` are two
  # parameterizations of the same saturated model, so they fit the same cell
  # probabilities and the weights coincide. Measured max relative gap on this
  # fixture is 1.5e-7, which is where the multinomial's own convergence stops.
  gap <- max(
    abs(as.double(two$wts) - as.double(joint$wts)) / as.double(joint$wts)
  )
  expect_lt(gap, 1e-5)
})

# ---- the reported surface ---------------------------------------------------

test_that("ipw() over two treatment models reports the joint surface", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)
  res <- ipw(two$models, two$outcome_mod)
  est <- res$estimates

  expect_s3_class(res, "ipw")
  expect_identical(res$estimand, "ate")
  expect_identical(res$se_method, "mestimation")

  # The same three identity columns in the same order the declared-exposure
  # route uses, so a row means the same thing whichever route produced it.
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

  expected <- joint_models_expected_rows()
  expect_identical(nrow(est), 14L)
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)

  # No odds ratio anywhere: it is noncollapsible, so neither a simple effect
  # beside one nor a difference of two of them says what it appears to.
  expect_false(any(est$effect == "log(or)"))
})

test_that("ipw() over two treatment models hand-computes its cell risks and contrasts", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)
  res <- ipw(two$models, two$outcome_mod)

  # Every reported quantity is g-computation on the weighted outcome model,
  # standardized over the whole sample, with the two treatments set together.
  mu <- joint_models_cell_means(two$outcome_mod, dat)
  expect_equal(
    res$estimates$estimate,
    joint_models_expected_estimates(mu),
    tolerance = 1e-8
  )

  # The interaction written out over the cells, so the identity the labels claim
  # is pinned against arithmetic rather than against the helper above.
  interaction <- res$estimates$estimate[
    res$estimates$effect == "rd" & res$estimates$group == "e = 1 vs e = 0"
  ]
  expect_equal(
    interaction,
    unname(mu[[4]] - mu[[2]] - mu[[3]] + mu[[1]]),
    tolerance = 1e-8
  )
})

test_that("ipw() over two treatment models reports diffs for a continuous outcome", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat, outcome_family = "gaussian")
  res <- ipw(two$models, two$outcome_mod)
  est <- res$estimates

  expected <- joint_models_expected_rows(forms = "diff")
  expect_identical(nrow(est), 9L)
  expect_identical(est$effect, expected$effect)
  expect_identical(est$contrast, expected$contrast)
  expect_identical(est$group, expected$group)

  mu <- joint_models_cell_means(two$outcome_mod, dat)
  expect_equal(
    est$estimate,
    joint_models_expected_estimates(mu, forms = "diff"),
    tolerance = 1e-8
  )
})

# ---- equivalence with the declared-exposure route ---------------------------

test_that("the two routes agree as estimators of the same joint effects", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_models()

  two <- fit_joint_models_route(dat)
  joint <- fit_joint_exposure_route(dat)
  res_two <- ipw(two$models, two$outcome_mod)
  res_joint <- ipw(joint$ps_mod, joint$outcome_mod)

  # The same rows in the same order, so the comparison is row for row.
  expect_identical(res_two$estimates$effect, res_joint$estimates$effect)
  expect_identical(res_two$estimates$contrast, res_joint$estimates$contrast)
  expect_identical(res_two$estimates$group, res_joint$estimates$group)

  # They target the same estimand through different models of the same
  # propensity score, so what separates them is the difference between those
  # models rather than solver noise. Measured on this fixture the reported rows
  # differ by at most 1.5e-2 in absolute terms, on cell risks of 0.36 to 0.84
  # and contrasts of 0.09 to 0.63, while the weights behind them differ by up to
  # 32 percent. The bound is set from that measurement with room, and what it
  # pins is that the two routes answer the same question, not that they compute
  # the same number.
  expect_lt(
    max(abs(res_two$estimates$estimate - res_joint$estimates$estimate)),
    0.05
  )

  # The cell risks are the primary quantities and they are all probabilities, so
  # they carry a tighter bound of their own; measured 6e-3.
  means <- res_two$estimates$effect == "mean"
  expect_lt(
    max(abs(
      res_two$estimates$estimate[means] - res_joint$estimates$estimate[means]
    )),
    0.03
  )
})

test_that("saturated parameterizations make the two routes agree tightly", {
  skip_if_not_installed("nnet")
  dat <- sim_joint_models_saturated()

  two <- fit_joint_models_route(
    dat,
    a_rhs = "x2",
    e_rhs = "a * x2",
    outcome_rhs = "a * e + x2"
  )
  joint <- fit_joint_exposure_route(
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
  # Measured on this fixture the rows agree to 3.2e-14 absolute and 9.3e-14
  # relative, against weights agreeing to 1.5e-7. The bound is set well above
  # the measurement to leave room for the two solvers.
  expect_equal(
    res_two$estimates$estimate,
    res_joint$estimates$estimate,
    tolerance = 1e-6
  )
})

# ---- standard errors ---------------------------------------------------------

test_that("ipw() over two treatment models reports a usable standard error for every row", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)
  res <- ipw(two$models, two$outcome_mod)
  est <- res$estimates

  expect_identical(nrow(est), 14L)
  expect_true(all(is.finite(est$std.err)))
  expect_true(all(est$std.err > 0))
  expect_true(all(is.finite(est$ci.lower)))
  expect_true(all(is.finite(est$ci.upper)))
  expect_true(all(est$ci.lower < est$ci.upper))
  expect_equal(sqrt(diag(vcov(res))), est$std.err, ignore_attr = TRUE)
})

test_that("the stacked system carries both treatment models", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)
  res <- ipw(two$models, two$outcome_mod)

  # The claim this pins is that both propensity score models are estimated
  # alongside the outcome model rather than plugged in as fixed weights, which
  # is what makes the reported standard errors account for having estimated
  # them. It is pinned through the dimension of the stacked system rather than
  # by comparison against a fixed-weight fit: `ipw()` has no route that takes
  # product weights as given, and the outcome model's own standard errors live
  # on the coefficient scale rather than on the effect scale, so there is no
  # like-for-like comparator to make. The dimension is exact and it fails
  # loudly for the defect in question.
  #
  # Three coefficients for `a ~ x1 + x2`, five for `e ~ a * x1 + x2`, five for
  # the outcome model, four cell means, and ten contrast rows.
  expect_identical(as.integer(res$fit@n_params), 27L)
  expect_identical(nobs(res), 900L)
  expect_identical(df.residual(res), 900L - 27L)

  # Both models' coefficients are in there under names of their own.
  theta <- names(res$fit@theta)
  expect_true(all(c("x1", "x2") %in% theta))
  expect_true(any(grepl("a:x1", theta, fixed = TRUE)))

  # The second model's dimension is what moves when the second model changes, so
  # the count tracks the pair rather than carrying a fixed allowance for it. An
  # additive second model has one coefficient fewer and the system is one
  # parameter smaller.
  additive <- fit_joint_models_route(dat, e_rhs = "a + x1 + x2")
  res_additive <- ipw(additive$models, additive$outcome_mod)
  expect_identical(as.integer(res_additive$fit@n_params), 26L)
})

test_that("the covariance of a two-model joint fit couples its blocks", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)
  res <- ipw(two$models, two$outcome_mod)
  covariance <- vcov(res)

  # Every reported row is a function of the same four cell means, solved as one
  # stacked system on one sample, so the interaction covaries with each simple
  # effect it is built from and two simple effects sharing a cell covary with
  # each other. An assembly that computed any of these separately and stitched
  # it in would report exact zeros here while matching its own standard errors
  # everywhere. The floor is absolute and excludes the structural zero rather
  # than a small number.
  couples <- list(
    c("rd a: 1 vs 0 e = 1 vs e = 0", "rd a: 1 vs 0 e = 1"),
    c("rd a: 1 vs 0 e = 1 vs e = 0", "rd a: 1 vs 0 e = 0"),
    c("rd a: 1 vs 0 e = 0", "rd e: 1 vs 0 a = 0"),
    c("mean a = 1, e = 1 overall", "rd a: 1 vs 0 e = 1 vs e = 0")
  )
  for (pair in couples) {
    where <- paste(pair, collapse = " with ")
    present <- all(pair %in% rownames(covariance))
    expect_true(present, label = where)
    if (!present) {
      next
    }
    entry <- covariance[pair[[1]], pair[[2]]]
    expect_true(is.finite(entry), label = where)
    expect_gt(abs(entry), 1e-8, label = where)
  }

  expect_equal(covariance, t(covariance), tolerance = 1e-12)
})

test_that("a two-model joint fit labels every surface alike", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)
  res <- ipw(two$models, two$outcome_mod)

  labels <- joint_models_labels(res$estimates)
  expect_identical(anyDuplicated(labels), 0L)
  expect_identical(names(coef(res)), labels)
  expect_identical(dimnames(vcov(res)), list(labels, labels))
  expect_identical(rownames(confint(res)), labels)
  expect_identical(
    dimnames(attr(res$estimates, "ipw_vcov", exact = TRUE)),
    list(labels, labels)
  )

  frame <- as.data.frame(res)
  expect_identical(
    names(frame),
    c(
      "term",
      "contrast",
      "group",
      "estimate",
      "std.error",
      "statistic",
      "p.value"
    )
  )
  expect_identical(frame$group, res$estimates$group)

  tidied <- tidy(res)
  expect_s3_class(tidied, "tbl_df")
  expect_identical(names(tidied), names(frame))
})

test_that("the conditional reading of a two-model joint fit is the outcome model's", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)
  conditional <- causalgenerics::as_conditional(
    ipw(two$models, two$outcome_mod)
  )

  # The container changes what is reported, not what is fitted. The outcome
  # model reads the two treatments and a covariate, and the conditional reading
  # is its coefficients with the covariance the joint estimation implies.
  expect_identical(coef(conditional), coef(two$outcome_mod))
  expect_identical(
    rownames(vcov(conditional)),
    names(coef(two$outcome_mod))
  )
  # Read off the surface the conditional reading presents rather than off the
  # stored frame. `as_conditional()` flips which reading a result reports
  # without rebuilding anything, so `estimates` stays the marginal frame and
  # legitimately keeps the group column its rows are keyed by; so does
  # `as.data.frame()`, which tabulates that frame whichever reading is active.
  # `tidy()` is the surface that answers in the reading, and the coefficients it
  # reports belong to no subgroup.
  expect_false("group" %in% names(tidy(conditional)))
})

# ---- refusals ----------------------------------------------------------------

test_that("ipw() refuses an outcome model that does not read both treatments", {
  dat <- sim_joint_models()

  # A joint surface sets both treatments at once, so an outcome model reading
  # one of them has no counterfactual design for three of the four cells. Both
  # omissions are the same fault and are reported under the class the package
  # already uses for an outcome model that does not read the exposure.
  missing_e <- fit_joint_models_route(dat, outcome_rhs = "a + x1")
  expect_error(
    ipw(missing_e$models, missing_e$outcome_mod),
    class = "propensity_ipw_exposure_error"
  )

  missing_a <- fit_joint_models_route(dat, outcome_rhs = "e + x1")
  expect_error(
    ipw(missing_a$models, missing_a$outcome_mod),
    class = "propensity_ipw_exposure_error"
  )

  expect_propensity_error(ipw(missing_e$models, missing_e$outcome_mod))
})

test_that("ipw() refuses weights that are not a joint product", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)

  # The container says the analysis is a joint intervention on two treatments,
  # and the weights have to be the product that names. A single treatment's
  # weights are an ordinary psw of the right length and the right estimand, so
  # nothing else downstream would notice.
  #
  # What is pinned here is the minimum: weights carrying no `joint_wt_meta`
  # refuse. Whether the product is recomputed from the container's own models
  # and compared, the way the binary path's weight-consistency preflight does,
  # is a heavier guarantee and is left to the estimator rounds.
  single <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = withr::with_options(
      list(propensity.quiet = TRUE),
      wt_ate(two$ps_a)
    ),
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  expect_false(is_joint_wt(extract_weights(single)))
  expect_error(
    ipw(two$models, single),
    class = "propensity_ipw_joint_models_weights_error"
  )
  expect_propensity_error(ipw(two$models, single))
})

test_that("ipw() refuses a joint fit whose weights are not the ate", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)

  # `wt_joint()` refuses a non-ate component, so a non-ate product cannot be
  # built at all and the estimand refusal is reached only by handing over
  # weights that are not a product. That is the previous refusal, so what is
  # pinned here is the other half: the surface is defined for the ate, and the
  # result says so.
  res <- ipw(two$models, two$outcome_mod)
  expect_identical(res$estimand, "ate")
  expect_identical(estimand(two$wts), "ate")
  expect_error(
    wt_joint(
      withr::with_options(
        list(propensity.quiet = TRUE),
        wt_att(two$ps_a)
      ),
      withr::with_options(
        list(propensity.quiet = TRUE),
        wt_ate(two$ps_e)
      )
    ),
    class = "propensity_wt_joint_estimand_error"
  )
})

test_that("the weights mismatch on two binary treatments names no focal level", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)

  # A joint product scaled after it was built is still the product the
  # container asks for, so the route reaches the preflight that rebuilds the
  # weights from the two models and reports that they disagree. A joint
  # specification targets the joint ate and resolves no focal level, so naming
  # a focal level among the causes would send the reader after a setting this
  # route never read.
  scaled <- two$wts * 1.5
  expect_true(is_joint_wt(scaled))
  outcome_mod <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = scaled,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  err <- expect_error(
    ipw(two$models, outcome_mod),
    class = "propensity_ipw_weights_mismatch_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_no_match(msg, "focal level", fixed = TRUE)

  # On this route `wt_mod` is the container of the two treatment models rather
  # than one propensity score model, so the remedy names what it holds. A
  # reader told to refit the outcome model with weights from "this propensity
  # score model" is being sent to a model the call never had.
  expect_match(msg, "two treatment models", fixed = TRUE)
  expect_no_match(msg, "this propensity score model", fixed = TRUE)

  expect_propensity_error(ipw(two$models, outcome_mod))
})

test_that("ipw() refuses .by on the two-model route", {
  dat <- sim_joint_models()
  dat$grp <- factor(ifelse(dat$x1 > 0, "hi", "lo"), levels = c("lo", "hi"))
  # The outcome model carries the modifier interacted with both treatments, so
  # the only thing wrong with the request is the combination itself.
  two <- fit_joint_models_route(dat, outcome_rhs = "a * e * grp + x1")

  # Effect modification of a joint intervention is a three-way question whether
  # the crossing arrived declared or as two models, so the refusal is the one
  # the declared route raises.
  expect_no_warning(expect_error(
    ipw(two$models, two$outcome_mod, .by = grp),
    class = "propensity_ipw_joint_by_error"
  ))
  expect_propensity_error(ipw(two$models, two$outcome_mod, .by = grp))
})

test_that("ipw() refuses linearization standard errors on the two-model route", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)

  # The cell means and every contrast built from them are parameters of the
  # stacked system, and the linearization path solves no such system.
  expect_error(
    ipw(two$models, two$outcome_mod, se_method = "linearization"),
    class = "propensity_ipw_joint_models_method_error"
  )
  expect_propensity_error(
    ipw(two$models, two$outcome_mod, se_method = "linearization")
  )
})

test_that("ipw() refuses robust standard errors on the two-model route", {
  dat <- sim_joint_models()
  two <- fit_joint_models_route(dat)

  # The cell means of a joint treatment are parameters of the stacked system,
  # so the weights-fixed sandwich of a two-cell outcome model describes none of
  # them. The refusal names the method that was asked for.
  expect_error(
    ipw(two$models, two$outcome_mod, se_method = "robust"),
    class = "propensity_ipw_joint_models_method_error",
    regexp = "robust"
  )
  expect_propensity_error(
    ipw(two$models, two$outcome_mod, se_method = "robust")
  )
})

# ---- treatment models fit under case weights --------------------------------
#
# Every score this route stacks is unweighted, so a treatment model fit with
# prior case weights is not at the root of the equation written for it: the
# solve, seeded at its coefficients, would move away from them and report a
# treatment model nobody fit. The single-treatment routes refuse such a fit by
# name, and this route reads the same field of the same fitted objects, so the
# refusal is the same one and is made of each component in turn.

test_that("ipw() refuses a joint binary treatment model fit with case weights", {
  dat <- sim_joint_models()
  dat$cw <- rep(c(1, 2), length.out = nrow(dat))

  # Integer weights, so the binomial fit raises nothing of its own and the only
  # condition in the test is the refusal being asserted.
  weighted_route <- function(which) {
    a_weights <- if (identical(which, "a")) dat$cw else rep(1, nrow(dat))
    e_weights <- if (identical(which, "e")) dat$cw else rep(1, nrow(dat))
    ps_a <- glm(
      a ~ x1 + x2,
      data = dat,
      family = binomial(),
      weights = a_weights
    )
    ps_e <- glm(
      e ~ a * x1 + x2,
      data = dat,
      family = binomial(),
      weights = e_weights
    )
    wts <- withr::with_options(
      list(propensity.quiet = TRUE),
      wt_joint(
        wt_ate(ps_a),
        wt_ate(ps_e),
        exposure_type = c("binary", "binary")
      )
    )
    outcome_mod <- lm(yc ~ a * e + x1, data = dat, weights = wts)

    list(
      models = joint_wt_models(a = ps_a, e = ps_e),
      outcome_mod = outcome_mod
    )
  }

  # The message of whatever the call raised, or the empty string when it raised
  # nothing, so that one component failing to refuse does not stop the other
  # from being asserted.
  refusal_message <- function(expr) {
    cnd <- tryCatch(expr, error = function(e) e)
    if (!inherits(cnd, "condition")) {
      return("")
    }
    gsub("[[:space:]]+", " ", conditionMessage(cnd))
  }

  for (which in c("a", "e")) {
    fx <- weighted_route(which)
    expect_error(
      ipw(fx$models, fx$outcome_mod),
      class = "propensity_ipw_ps_weights_error",
      label = which
    )

    # The refusal names the component that carries the weights rather than the
    # container the two models arrived in, which is not a model to refit.
    msg <- refusal_message(ipw(fx$models, fx$outcome_mod))
    expect_match(msg, paste0("`", which, "`"), fixed = TRUE)
    expect_no_match(msg, "wt_mod", fixed = TRUE)
  }
})

# ---- stabilized discrete components -----------------------------------------
#
# Both components of a discrete pair may be stabilized, and each numerator is
# estimated in a block of its own, whether it is the marginal proportion the
# default stabilizer is or a model the caller fit. What the route used to do
# instead was rebuild every discrete component unstabilized, which reached a
# different product and was reported as a weights mismatch the caller had not
# caused.

test_that("ipw() stacks a stabilized discrete component's own numerator", {
  dat <- sim_joint_models()

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- glm(e ~ a * x1 + x2, data = dat, family = binomial())
  # One component takes the default numerator, the marginal proportion, and the
  # other a fitted model of its treatment on a covariate the outcome model
  # reads, which is what a numerator may condition on.
  num_e <- glm(e ~ x1, data = dat, family = binomial())

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a, stabilize = TRUE),
      wt_ate(ps_e, stabilize = num_e),
      exposure_type = c("binary", "binary")
    )
  )
  outcome_mod <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  models <- joint_wt_models(a = ps_a, e = ps_e)

  # The product written out, so the rebuild is held to the arithmetic rather
  # than to the preflight's tolerance.
  p1 <- mean(dat$a)
  p_e <- as.numeric(fitted(num_e))
  expected <- (((dat$a * p1) / fitted(ps_a)) +
    (((1 - dat$a) * (1 - p1)) / (1 - fitted(ps_a)))) *
    (((dat$e / fitted(ps_e)) + ((1 - dat$e) / (1 - fitted(ps_e)))) *
      (dat$e * p_e + (1 - dat$e) * (1 - p_e)))
  expect_equal(as.numeric(wts), unname(expected), tolerance = 1e-12)

  spec <- ipw_spec_joint_models(models, outcome_mod)
  layout <- ipw_theta_layout(spec)
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(wts),
    tolerance = 1e-12
  )

  # Every equation in the stacked system sits at its root at the seed, which is
  # what makes the solve report the models the caller fit.
  mat <- build_ipw_psi(spec, layout)(layout$init)
  expect_false(anyNA(mat))
  expect_true(all(abs(rowSums(mat)) / spec$n < 1e-8))

  res <- ipw(models, outcome_mod)
  expect_s3_class(res, "ipw")
  expect_equal(
    res$estimates$estimate,
    joint_models_expected_estimates(joint_models_cell_means(outcome_mod, dat)),
    tolerance = 1e-6
  )
  expect_true(all(is.finite(res$estimates$std.err)))

  # One proportion for the component that took the default numerator, and one
  # parameter per coefficient for the component that took a model, each named
  # for the component it belongs to.
  theta <- coef(res$fit)
  stab <- theta[grepl("^stab_", names(theta))]
  expect_length(stab, 1L + length(coef(num_e)))
  expect_equal(unname(theta[["stab_a_pi"]]), p1, tolerance = 1e-8)
  expect_equal(
    unname(stab[grepl("^stab_e_", names(stab))]),
    unname(coef(num_e)),
    tolerance = 1e-6
  )
})

# ---- a component stabilized on a fixed score --------------------------------

test_that("ipw() rebuilds a discrete component stabilized on a fixed score", {
  dat <- sim_joint_models()

  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- glm(e ~ a * x1 + x2, data = dat, family = binomial())

  # A numerator the caller computed rather than one the estimator can fit. It
  # is away from the marginal proportion, which is the record a stabilized
  # discrete component used to leave whatever it was built with, so a system
  # that stood the proportion in would reach a different product and be
  # reported as a mismatch the caller did not cause.
  score <- 0.3 + 0.4 * plogis(dat$x2)

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a, stabilize = TRUE, stabilization_score = score),
      wt_ate(ps_e),
      exposure_type = c("binary", "binary")
    )
  )
  outcome_mod <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )
  models <- joint_wt_models(a = ps_a, e = ps_e)

  # The product written out, so the rebuild is held to the arithmetic rather
  # than to the preflight's tolerance.
  expected <- (((dat$a / fitted(ps_a)) +
    ((1 - dat$a) / (1 - fitted(ps_a)))) *
    score) *
    ((dat$e / fitted(ps_e)) + ((1 - dat$e) / (1 - fitted(ps_e))))
  expect_equal(as.numeric(wts), unname(expected), tolerance = 1e-12)

  spec <- ipw_spec_joint_models(models, outcome_mod)
  layout <- ipw_theta_layout(spec)
  expect_equal(
    as.double(ipw_weights_at_init(spec, layout)),
    as.double(wts),
    tolerance = 1e-12
  )

  mat <- build_ipw_psi(spec, layout)(layout$init)
  expect_false(anyNA(mat))
  expect_true(all(abs(rowSums(mat)) / spec$n < 1e-8))

  res <- ipw(models, outcome_mod)
  expect_s3_class(res, "ipw")
  expect_equal(
    res$estimates$estimate,
    joint_models_expected_estimates(joint_models_cell_means(outcome_mod, dat)),
    tolerance = 1e-6
  )
  expect_true(all(is.finite(res$estimates$std.err)))

  # A score is a known multiplier, so it widens the system by nothing at all
  # where the marginal proportion the default stabilizer estimates would have
  # taken one parameter.
  expect_false(any(grepl("^stab_", names(coef(res$fit)))))
})

# ---- a component's numerator model whose data is gone -----------------------
#
# A numerator model's coefficients are estimated in the stacked system, so the
# design they were fit over is what the route needs rather than the
# probabilities they evaluate to. A `glm` usually keeps its model frame, so that
# design is usually there to be read. One fit with `model = FALSE` keeps none
# and rebuilds it by re-evaluating the fitting call, which a fit made inside a
# wrapper whose frame is gone cannot do. What that owes the caller is the
# request every other route makes of a design it cannot recover: name the model
# that cannot be rebuilt and ask for `.data`, which rebuilds it.

# The pair of fits this section weights, with the second component's numerator
# left to the caller so the same product can be built from a fit that kept its
# frame and from one that did not.
joint_models_numerator_gone_fits <- function(dat, numerator) {
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- glm(e ~ a * x1 + x2, data = dat, family = binomial())

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a, stabilize = TRUE),
      wt_ate(ps_e, stabilize = numerator),
      exposure_type = c("binary", "binary")
    )
  )
  outcome_mod <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(models = joint_wt_models(a = ps_a, e = ps_e), outcome_mod = outcome_mod)
}

# The formula is written in the calling frame and the fitting call names a
# variable that lives only inside the wrapper, so nothing can rebuild the
# numerator's design once the wrapper has returned. Its fitted probabilities are
# still readable, and the weights were built from them.
joint_models_numerator_frame_gone <- function(dat) {
  fmla <- e ~ x1
  fit_in_function <- function(fitting_data) {
    glm(fmla, data = fitting_data, family = binomial(), model = FALSE)
  }

  fit_in_function(dat)
}

test_that("a joint component's numerator model whose data is gone asks for .data", {
  dat <- sim_joint_models()
  gone <- joint_models_numerator_frame_gone(dat)

  expect_error(model.matrix(gone))

  fits <- joint_models_numerator_gone_fits(dat, gone)

  cnd <- tryCatch(
    ipw(fits$models, fits$outcome_mod),
    error = function(e) e
  )
  expect_s3_class(cnd, "propensity_ipw_data_error")

  # Both treatment models here are rebuildable, so a message naming one of them
  # would send the caller to the wrong fit.
  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "stabilize", fixed = TRUE)
  expect_match(message, "Supply `.data`", fixed = TRUE)
})

test_that("a joint component's numerator model whose data is gone is rebuilt from .data", {
  dat <- sim_joint_models()
  gone <- joint_models_numerator_frame_gone(dat)

  # The other half of the contract: the `.data` the refusal above asks for
  # rebuilds the numerator's design from the fit's own terms and contrasts, and
  # what comes back is what the same numerator reports when its frame was never
  # lost.
  kept <- glm(e ~ x1, data = dat, family = binomial())
  expect_equal(unname(coef(gone)), unname(coef(kept)), tolerance = 1e-10)

  gone_fits <- joint_models_numerator_gone_fits(dat, gone)
  kept_fits <- joint_models_numerator_gone_fits(dat, kept)

  res_data <- ipw(gone_fits$models, gone_fits$outcome_mod, .data = dat)
  res_recovered <- ipw(kept_fits$models, kept_fits$outcome_mod)

  expect_s3_class(res_data, "ipw")
  expect_equal(
    res_data$estimates$estimate,
    res_recovered$estimates$estimate,
    tolerance = 1e-6
  )
  expect_equal(
    res_data$estimates$std.err,
    res_recovered$estimates$std.err,
    tolerance = 1e-6
  )
})

# ---- a binary component's numerator fit on an unsupported link ---------------
#
# The block written for a binary numerator is the binomial score of a fit on one
# of the three links deli writes that score for, so a fit on another link is
# refused rather than stacked. Both components' numerators arrive in the same
# argument, so the refusal names the component whose weights this one was built
# for, the way the route's other refusals name a component.

test_that("a binary component's numerator on an unsupported link names the component", {
  dat <- sim_joint_models()

  # A cauchit fit is a real fit of the treatment whose fitted probabilities the
  # weight layer takes, so the numerator reaches the estimator recorded on the
  # product and the refusal is the estimator's to make.
  num_e <- glm(e ~ x1, data = dat, family = binomial("cauchit"))
  fits <- joint_models_numerator_gone_fits(dat, num_e)

  cnd <- tryCatch(
    ipw(fits$models, fits$outcome_mod),
    error = identity
  )
  expect_s3_class(cnd, "propensity_ipw_link_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "cauchit", fixed = TRUE)
  expect_match(message, "`stabilize` for `e`", fixed = TRUE)
})

# ---- the refusals a numerator shares with the other stacked models ----------
#
# A numerator model meets the guards every model this route stacks meets: the
# refusal of case weights, the requirement of a coefficient per design column,
# the recovery of a design from a fit that keeps no model frame, the width a
# `.data` rebuild has to come back at, and the refusal of a penalized fit. Each
# of those names the model by the argument it arrived in, which on this route is
# one `stabilize` per component, so each names the component beside it.

test_that("a numerator fit with case weights names the component", {
  dat <- sim_joint_models()
  dat$case_wt <- rep(c(1, 2), length.out = nrow(dat))
  weighted <- glm(e ~ x1, data = dat, family = binomial(), weights = case_wt)

  # `wt_ate()` refuses such a model where it arrives, so weights recording one
  # were assembled by hand or by a version that took it. The record is written
  # here rather than built, which is what those weights are, and the estimator
  # reads it the same way either way.
  ps_a <- glm(a ~ x1 + x2, data = dat, family = binomial())
  ps_e <- glm(e ~ a * x1 + x2, data = dat, family = binomial())
  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a, stabilize = TRUE),
      wt_ate(ps_e, stabilize = glm(e ~ x1, data = dat, family = binomial())),
      exposure_type = c("binary", "binary")
    )
  )
  meta <- joint_wt_meta(wts)
  meta$numerator_model[[2]] <- weighted
  attr(wts, "joint_wt_meta") <- meta

  outcome_mod <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  cnd <- tryCatch(
    ipw(joint_wt_models(a = ps_a, e = ps_e), outcome_mod),
    error = identity
  )
  expect_s3_class(cnd, "propensity_ipw_ps_weights_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`stabilize` for `e`", fixed = TRUE)
  expect_no_match(message, "wt_mod", fixed = TRUE)
})

test_that("a numerator missing a coefficient names the component", {
  dat <- sim_joint_models()

  # A duplicated column leaves the fit without a coefficient for it, which is
  # the state the block cannot multiply its design by.
  duplicated <- dat
  duplicated$x1_again <- duplicated$x1
  num_e <- glm(e ~ x1 + x1_again, data = duplicated, family = binomial())
  fits <- joint_models_numerator_gone_fits(duplicated, num_e)

  cnd <- tryCatch(
    muffle_coverage_warning(ipw(fits$models, fits$outcome_mod)),
    error = identity
  )
  expect_s3_class(cnd, "propensity_ipw_rank_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`stabilize` for `e`", fixed = TRUE)
  expect_match(message, "x1_again", fixed = TRUE)
})

test_that("a numerator whose data is gone names the component", {
  dat <- sim_joint_models()
  gone <- joint_models_numerator_frame_gone(dat)
  fits <- joint_models_numerator_gone_fits(dat, gone)

  cnd <- tryCatch(ipw(fits$models, fits$outcome_mod), error = identity)
  expect_s3_class(cnd, "propensity_ipw_data_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`stabilize` for `e`", fixed = TRUE)
})

test_that("a numerator design rebuilt to another width names the component", {
  dat <- sim_joint_models()

  # A matrix column is the shape the guards before the width count pass over,
  # so a frame holding one of another width reaches the count itself.
  dat$vm <- stats::model.matrix(~ factor(x2), data = dat)[, -1, drop = FALSE]
  num_e <- glm(e ~ vm, data = dat, family = binomial())
  fits <- joint_models_numerator_gone_fits(dat, num_e)

  supplied <- dat
  supplied$vm <- cbind(dat$vm, s = as.double(dat$x1 > 0))

  cnd <- tryCatch(
    muffle_coverage_warning(
      ipw(fits$models, fits$outcome_mod, .data = supplied)
    ),
    error = identity
  )
  expect_s3_class(cnd, "propensity_ipw_data_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`stabilize` for `e`", fixed = TRUE)
  expect_match(message, ".data", fixed = TRUE)
})

test_that("a numerator fit as an additive model names the component", {
  skip_if_not_installed("mgcv")
  dat <- sim_joint_models()

  # An additive fit's coefficients are the root of a penalized score, which is
  # not the score the block stacks, so the fit is refused rather than stacked.
  num_e <- mgcv::gam(e ~ s(x1), data = dat, family = binomial())
  fits <- joint_models_numerator_gone_fits(dat, num_e)

  cnd <- tryCatch(ipw(fits$models, fits$outcome_mod), error = identity)
  expect_s3_class(cnd, "propensity_ipw_se_method_unavailable_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "`stabilize` for `e`", fixed = TRUE)
})

# ---- a treatment model whose data is gone -----------------------------------
#
# A treatment model's coefficients are estimated in the stacked system, so the
# design they were fit over is what the route needs. A `glm` usually keeps its
# model frame, and one fit with `model = FALSE` rebuilds the design by
# re-evaluating its fitting call, which a fit made inside a wrapper whose frame
# is gone cannot do. Off the fit alone there is nothing left to read, and the
# route says so. With a `.data` the caller supplied there is: every design and
# every treatment column on this route is rebuilt from that frame rather than
# read off the fits, so the frame the caller still holds stands in for the one
# the fit lost, and the same product is reported either way.

# The fit nothing can rebuild. The formula is written in the wrapper's frame and
# the fitting call names a variable that lives only inside it, so re-evaluating
# the call once the wrapper has returned finds nothing. Its coefficients and
# fitted values are still readable.
joint_models_treatment_frame_gone <- function(dat) {
  fmla <- a ~ x1 + x2
  fit_in_function <- function(fitting_data) {
    glm(fmla, data = fitting_data, family = binomial(), model = FALSE)
  }

  fit_in_function(dat)
}

# The pair this section weights, with the first component's treatment model left
# to the caller so the same product can be built from a fit that kept its frame
# and from one that did not. `wt_ate()` reads a fit's exposure out of its model
# frame, so the weights a frame-gone fit implies are built through `weight_mod`,
# an intact fit holding the same coefficients, rather than through the fit
# handed to `joint_wt_models()`.
joint_models_treatment_gone_fits <- function(
  dat,
  treatment_mod,
  weight_mod = treatment_mod
) {
  ps_e <- glm(e ~ a * x1 + x2, data = dat, family = binomial())

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(weight_mod),
      wt_ate(ps_e),
      exposure_type = c("binary", "binary")
    )
  )
  outcome_mod <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(
    models = joint_wt_models(a = treatment_mod, e = ps_e),
    outcome_mod = outcome_mod
  )
}

test_that("a joint treatment model whose data is gone is rebuilt from .data", {
  dat <- sim_joint_models()
  gone <- joint_models_treatment_frame_gone(dat)
  kept <- glm(a ~ x1 + x2, data = dat, family = binomial())

  expect_error(model.matrix(gone))
  expect_equal(unname(coef(gone)), unname(coef(kept)), tolerance = 1e-10)

  gone_fits <- joint_models_treatment_gone_fits(dat, gone, weight_mod = kept)
  kept_fits <- joint_models_treatment_gone_fits(dat, kept)

  # The two calls weight the same product from the same coefficients and differ
  # only in where the first component's design comes from, so the frame the
  # caller supplied has to reach the design the lost frame held.
  res_data <- ipw(gone_fits$models, gone_fits$outcome_mod, .data = dat)
  res_recovered <- ipw(kept_fits$models, kept_fits$outcome_mod)

  expect_s3_class(res_data, "ipw")
  expect_equal(
    res_data$estimates$estimate,
    res_recovered$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    res_data$estimates$std.err,
    res_recovered$estimates$std.err,
    tolerance = 1e-10
  )
})

test_that("a joint treatment model whose data is gone asks to be refit", {
  dat <- sim_joint_models()
  gone <- joint_models_treatment_frame_gone(dat)
  kept <- glm(a ~ x1 + x2, data = dat, family = binomial())

  gone_fits <- joint_models_treatment_gone_fits(dat, gone, weight_mod = kept)

  # Without `.data` there is no frame to rebuild the design from, and the only
  # remedy left is the fit itself, so the refusal asks for that rather than for
  # an argument that is not there.
  cnd <- tryCatch(
    ipw(gone_fits$models, gone_fits$outcome_mod),
    error = function(e) e
  )
  expect_s3_class(cnd, "propensity_ipw_data_error")

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "data behind a treatment model", fixed = TRUE)
  expect_match(
    message,
    "Refit the treatment models where the data they were fit to is still available.",
    fixed = TRUE
  )
})

# ---- a treatment column .data carries under the wrong order -----------------
#
# With a `.data` the caller supplied, every treatment on this route is read out
# of that frame by the name its component is registered under. A column that
# carries its own level order is compared against the order its model was fit
# with, so a re-levelled factor supplied the other way round is refused by name.
# A column that carries no order of its own is not compared at all: it is
# levelled alphabetically wherever it is read, and when the fit's order is not
# alphabetical the indicator the route rebuilds is the fit's indicator inverted.
#
# Nothing about that is silent, but every refusal it reaches is about something
# else. The weights recomputed under the inverted indicator no longer match the
# ones that fit the outcome model, so the weights-at-init preflight refuses and
# blames the weights; where the counterfactual designs collapse first the
# refusal blames the outcome model's contrast coding. Neither sends the reader
# to the column that is actually miscoded, which is what the two sections below
# ask for, on the class the sibling guards already refuse a miscoded `.data`
# exposure with.

# `sim_joint_models()` with the first treatment re-expressed as a factor whose
# levels run "t" then "c". The fit therefore holds "t" as its reference and "c"
# as its exposed level, which is the reverse of what levelling the same values
# alphabetically would say.
sim_joint_models_relevelled <- function(...) {
  dat <- sim_joint_models(...)
  dat$af <- factor(ifelse(dat$a == 1L, "t", "c"), levels = c("t", "c"))

  dat
}

# The pair this section weights. Both treatment models, the weights and the
# outcome model read `af` under the order the frame declares, so the fit is
# internally consistent and only a `.data` that declares a different order can
# move it.
joint_models_relevelled_fits <- function(dat) {
  ps_af <- glm(af ~ x1 + x2, data = dat, family = binomial())
  ps_e <- glm(e ~ af * x1 + x2, data = dat, family = binomial())

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_af),
      wt_ate(ps_e),
      exposure_type = c("binary", "binary")
    )
  )
  outcome_mod <- glm(
    y ~ af * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(
    models = joint_wt_models(af = ps_af, e = ps_e),
    outcome_mod = outcome_mod
  )
}

test_that("a re-levelled joint treatment supplied as that factor is accepted", {
  dat <- sim_joint_models_relevelled()
  fits <- joint_models_relevelled_fits(dat)

  # The boundary the two refusals below sit against: a `.data` whose treatment
  # column declares exactly the order the fit holds is the frame the fits kept,
  # and the route reports the surface it reports without one.
  res_data <- ipw(fits$models, fits$outcome_mod, .data = dat)
  res_frames <- ipw(fits$models, fits$outcome_mod)

  expect_s3_class(res_data, "ipw")
  expect_equal(
    res_data$estimates$estimate,
    res_frames$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    res_data$estimates$std.err,
    res_frames$estimates$std.err,
    tolerance = 1e-10
  )
})

test_that("a joint treatment supplied as a character column is refused by name", {
  dat <- sim_joint_models_relevelled()
  fits <- joint_models_relevelled_fits(dat)

  # The same values, stripped of the order that told the route which level the
  # fit treats as exposed.
  as_character <- dat
  as_character$af <- as.character(dat$af)

  cnd <- tryCatch(
    ipw(fits$models, fits$outcome_mod, .data = as_character),
    error = function(e) e
  )

  expect_s3_class(cnd, "propensity_ipw_data_error")

  # The two refusals this reaches instead, each of which names something the
  # caller has no reason to change.
  expect_false(inherits(cnd, "propensity_ipw_weights_mismatch_error"))
  expect_false(inherits(cnd, "propensity_ipw_msm_error"))

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "\"af\"", fixed = TRUE)
  expect_match(message, "levels", fixed = TRUE)
})

# The same fault reached through a model that transforms its response in place.
# `joint_wt_models()` names a component by `all.vars()` of its response, so a
# fit written as `factor(a, levels = c("t", "c")) ~ .` registers as the
# component named "a" and passes the name check the pair is built under. What
# the `.data` rebuild then reads is the raw column `a`, which carries the two
# strings and no order, so the route levels it alphabetically and inverts the
# indicator the fit was made under.
sim_joint_models_transformed <- function(...) {
  dat <- sim_joint_models(...)
  dat$a <- ifelse(dat$a == 1L, "t", "c")

  dat
}

joint_models_transformed_fits <- function(dat) {
  ps_a <- glm(
    factor(a, levels = c("t", "c")) ~ x1 + x2,
    data = dat,
    family = binomial()
  )
  ps_e <- glm(e ~ a * x1 + x2, data = dat, family = binomial())

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      wt_ate(ps_e),
      exposure_type = c("binary", "binary")
    )
  )
  outcome_mod <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(
    models = joint_wt_models(a = ps_a, e = ps_e),
    outcome_mod = outcome_mod
  )
}

test_that("a joint treatment model that levels its own response is refused by name", {
  dat <- sim_joint_models_transformed()
  fits <- joint_models_transformed_fits(dat)

  # `.data` here is the frame every model was fit to, so nothing about the rows
  # or the values differs. What differs is the order the raw column implies
  # against the order the response transform declared.
  cnd <- tryCatch(
    ipw(fits$models, fits$outcome_mod, .data = dat),
    error = function(e) e
  )

  expect_s3_class(cnd, "propensity_ipw_data_error")

  expect_false(inherits(cnd, "propensity_ipw_weights_mismatch_error"))
  expect_false(inherits(cnd, "propensity_ipw_msm_error"))

  message <- gsub("[[:space:]]+", " ", conditionMessage(cnd))
  expect_match(message, "\"a\"", fixed = TRUE)
  expect_match(message, "levels", fixed = TRUE)
})

# The other side of the same boundary: a column that carries no order of its own
# is not the fault, so one whose values imply exactly the order its model was
# fit under is read the way it always was. Only the disagreement is refused.
sim_joint_models_alphabetical <- function(...) {
  dat <- sim_joint_models(...)
  dat$af <- factor(ifelse(dat$a == 1L, "t", "c"))

  dat
}

test_that("a joint treatment fit on the order its values imply is accepted", {
  # The pair reads `af` whatever order the column declares, so the fits this
  # section weights are the ones the re-levelled frame is weighted by.
  dat <- sim_joint_models_alphabetical()
  fits <- joint_models_relevelled_fits(dat)

  as_character <- dat
  as_character$af <- as.character(dat$af)

  res_data <- ipw(fits$models, fits$outcome_mod, .data = as_character)
  res_frames <- ipw(fits$models, fits$outcome_mod)

  expect_s3_class(res_data, "ipw")
  expect_equal(
    res_data$estimates$estimate,
    res_frames$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    res_data$estimates$std.err,
    res_frames$estimates$std.err,
    tolerance = 1e-10
  )
})

# A rename is not a re-ordering. `factor(a, labels = c("ctl", "trt"))` over a
# 0/1 column assigns its labels in sorted-value order, so a fit written that way
# codes the treatment exactly as the raw column implies and differs from it only
# in what the two groups are called. The label sets share nothing, which is
# precisely why the order comparison has nothing to say: the fitted labels are
# evidence of an order only when the column carries the same two labels. Where
# they do not, the weight-consistency preflight is the check that matters,
# because it compares the weights the rebuilt column produces rather than the
# names it produces them under.
joint_models_relabelled_fits <- function(dat, ps_a) {
  ps_e <- glm(e ~ a * x1 + x2, data = dat, family = binomial())

  wts <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_joint(
      wt_ate(ps_a),
      wt_ate(ps_e),
      exposure_type = c("binary", "binary")
    )
  )
  outcome_mod <- glm(
    y ~ a * e + x1,
    data = dat,
    family = quasibinomial(),
    weights = wts,
    control = glm.control(epsilon = 1e-14, maxit = 200)
  )

  list(
    models = joint_wt_models(a = ps_a, e = ps_e),
    outcome_mod = outcome_mod
  )
}

test_that("a joint treatment model that relabels its response agrees with the fitted frames", {
  dat <- sim_joint_models()
  ps_a <- glm(
    factor(a, labels = c("ctl", "trt")) ~ x1 + x2,
    data = dat,
    family = binomial()
  )
  fits <- joint_models_relabelled_fits(dat, ps_a)

  # `.data` here is the frame every model was fit to. The labels the response
  # transform assigned are the ones the column's sorted values carry, so the
  # route reads the coding the fit was made under and the two routes are the
  # same computation.
  res_data <- ipw(fits$models, fits$outcome_mod, .data = dat)
  res_frames <- ipw(fits$models, fits$outcome_mod)

  expect_s3_class(res_data, "ipw")
  expect_equal(
    res_data$estimates$estimate,
    res_frames$estimates$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    res_data$estimates$std.err,
    res_frames$estimates$std.err,
    tolerance = 1e-10
  )
})

test_that("a joint treatment model that relabels its response the other way round is caught by the weights", {
  dat <- sim_joint_models()
  # The same two labels assigned to the values the other way round, which does
  # invert the indicator the fit was made under.
  ps_a <- glm(
    factor(a, levels = c(1, 0), labels = c("trt", "ctl")) ~ x1 + x2,
    data = dat,
    family = binomial()
  )
  fits <- joint_models_relabelled_fits(dat, ps_a)

  cnd <- tryCatch(
    ipw(fits$models, fits$outcome_mod, .data = dat),
    error = function(e) e
  )

  # Nothing about the labels distinguishes this from the agreeing relabel, so
  # the refusal comes from the weights the inverted coding rebuilds rather than
  # from the level comparison.
  expect_s3_class(cnd, "propensity_ipw_weights_mismatch_error")
})
