# Product weights for two treatments intervened on at once.
#
# A joint intervention on A and E needs a weight for the joint exposure, which
# factorizes sequentially as f(A | L) f(E | A, L). `wt_joint()` builds it from
# the two component weight vectors, and `joint_wt_models()` records the two
# fitted treatment models the estimator later reads.
#
# Both constructors exist mostly to refuse. Multiplying two weight vectors is
# one line; what is worth a function is the set of things that make the product
# not a joint weight, and every refusal below is one of them:
#
#   * the product of two MARGINAL weights, f(A | L) f(E | L), is not a joint
#     weight at all. The second model must condition on the first treatment, and
#     that dependence should be modeled flexibly. This is the load-bearing
#     guardrail;
#   * a continuous component must be stabilized, or the density ratio has no
#     usable variance and the product inherits that;
#   * both components must target the ate, since the product targets the joint
#     ate and a tilted component targets a population the other does not.
#
# This file pins the two constructors and their validation. What the estimator
# does with them is a later concern.

# ---- data simulator ---------------------------------------------------------

# Two binary treatments, the second depending on the first; a continuous dose
# also depending on the first; a second dose depending on the first; and a
# three-level discrete treatment, for the multinomial arm of the factorization
# check. Every treatment depends on the covariate as well, so a model that omits
# the first treatment is still a plausible model rather than an empty one.
sim_joint_wt <- function(seed = 5301, n = 400) {
  withr::local_seed(seed)
  x1 <- rnorm(n)
  a <- rbinom(n, 1, plogis(0.3 * x1))
  e <- rbinom(n, 1, plogis(-0.2 + 0.5 * x1 - 0.8 * a))
  d <- 0.5 + 0.6 * x1 - 0.7 * a + rnorm(n)
  d2 <- -0.3 + 0.4 * x1 + 0.5 * d + rnorm(n)
  g <- factor(sample(c("p", "q", "r"), n, TRUE), levels = c("p", "q", "r"))
  data.frame(x1, a, e, d, d2, g)
}

# ---- model and weight fixtures ----------------------------------------------

# The treatment models the container records. Named for the treatment each one
# models and for whether it conditions on the first treatment: a `_marginal`
# model reads the covariate alone, which is the shape the guardrail refuses.
joint_wt_models_fixture <- function(dat) {
  list(
    a = glm(a ~ x1, data = dat, family = binomial()),
    e = glm(e ~ a * x1, data = dat, family = binomial()),
    e_marginal = glm(e ~ x1, data = dat, family = binomial()),
    # Each conditioning on the other, which specifies no factorization.
    a_on_e = glm(a ~ e * x1, data = dat, family = binomial()),
    d = lm(d ~ a * x1, data = dat),
    d_marginal = lm(d ~ x1, data = dat),
    d2 = lm(d2 ~ d + x1, data = dat),
    # A model whose response is not the treatment it would be named for.
    mislabelled = glm(e ~ a * x1, data = dat, family = binomial())
  )
}

joint_wt_multinom_fixture <- function(dat) {
  list(
    g = nnet::multinom(g ~ a + x1, data = dat, trace = FALSE),
    g_marginal = nnet::multinom(g ~ x1, data = dat, trace = FALSE)
  )
}

# The component weight vectors. Continuous components are stabilized, which the
# guardrail requires; `d_unstabilized` is the one that is not, for the refusal.
joint_wt_weights <- function(dat, mods) {
  quiet <- function(expr) {
    withr::with_options(list(propensity.quiet = TRUE), expr)
  }
  list(
    a = quiet(wt_ate(mods$a)),
    a_stabilized = quiet(wt_ate(mods$a, stabilize = TRUE)),
    a_att = quiet(wt_att(mods$a)),
    e = quiet(wt_ate(mods$e)),
    d = quiet(wt_ate(
      as.double(fitted(mods$d)),
      dat$d,
      exposure_type = "continuous",
      stabilize = TRUE
    )),
    d_unstabilized = quiet(wt_ate(
      as.double(fitted(mods$d)),
      dat$d,
      exposure_type = "continuous",
      stabilize = FALSE
    )),
    d2 = quiet(wt_ate(
      as.double(fitted(mods$d2)),
      dat$d2,
      exposure_type = "continuous",
      stabilize = TRUE
    ))
  )
}

# Everything a test needs, built once per test so nothing is shared across them.
joint_wt_fixture <- function(seed = 5301, n = 400) {
  dat <- sim_joint_wt(seed = seed, n = n)
  mods <- joint_wt_models_fixture(dat)
  list(dat = dat, mods = mods, w = joint_wt_weights(dat, mods))
}

# ---- the fixtures the guardrail tests rest on -------------------------------

test_that("the component weights and treatment models are what the guardrails read", {
  skip_if_not_installed("nnet")
  fx <- joint_wt_fixture()
  multinom_mods <- joint_wt_multinom_fixture(fx$dat)

  # Nothing below can be trusted if the fixture does not hold what the refusals
  # are written against, and none of it depends on the two constructors, so it
  # is checked here where a fixture fault shows up as a fixture fault.
  expect_true(all(vapply(fx$w, is_psw, logical(1))))
  expect_true(all(vapply(fx$w, length, integer(1)) == nrow(fx$dat)))
  expect_identical(
    vapply(fx$w, estimand, character(1))[c("a", "e", "d", "d2")],
    c(a = "ate", e = "ate", d = "ate", d2 = "ate")
  )
  expect_identical(estimand(fx$w$a_att), "att")

  # The stabilization the continuous guardrail turns on.
  expect_true(is_stabilized(fx$w$d))
  expect_true(is_stabilized(fx$w$d2))
  expect_false(is_stabilized(fx$w$d_unstabilized))
  expect_false(is_stabilized(fx$w$a))
  expect_true(is_stabilized(fx$w$a_stabilized))

  # Each model's response is the treatment it is named for, except the one built
  # to fail that check.
  responses <- function(mod) all.vars(formula(mod)[[2]])
  expect_identical(responses(fx$mods$a), "a")
  expect_identical(responses(fx$mods$e), "e")
  expect_identical(responses(fx$mods$d), "d")
  expect_identical(responses(multinom_mods$g), "g")
  expect_identical(responses(fx$mods$mislabelled), "e")

  # Which models read the first treatment and which do not, which is the whole
  # of the factorization check.
  reads <- function(mod, name) name %in% all.vars(formula(mod)[[3]])
  expect_true(reads(fx$mods$e, "a"))
  expect_true(reads(fx$mods$d, "a"))
  expect_true(reads(multinom_mods$g, "a"))
  expect_false(reads(fx$mods$e_marginal, "a"))
  expect_false(reads(fx$mods$d_marginal, "a"))
  expect_false(reads(multinom_mods$g_marginal, "a"))
  expect_false(reads(fx$mods$a, "e"))
  expect_true(reads(fx$mods$a_on_e, "e"))
})

# ---- wt_joint(): the product ------------------------------------------------

test_that("wt_joint() multiplies two binary ate weights", {
  fx <- joint_wt_fixture()
  joint <- wt_joint(fx$w$a, fx$w$e)

  # The values are the elementwise product and nothing else: the weight for the
  # cell a unit actually took is the product of the two treatment weights.
  expect_true(is_psw(joint))
  expect_equal(
    as.double(joint),
    as.double(fx$w$a) * as.double(fx$w$e),
    tolerance = 1e-12
  )
  expect_length(joint, nrow(fx$dat))

  # The product of two ate weights targets the joint ate.
  expect_identical(estimand(joint), "ate")

  # The provenance the estimator reads: that these are product weights, and
  # what each component was.
  expect_true(is_joint_wt(joint))
  expect_identical(
    joint_wt_meta(joint),
    list(
      exposure_type = c("binary", "binary"),
      stabilized = c(FALSE, FALSE)
    )
  )

  # A component weight is not a joint weight, which is what the flag is for.
  expect_false(is_joint_wt(fx$w$a))
  expect_null(joint_wt_meta(fx$w$a))
})

test_that("wt_joint() multiplies a binary weight by a stabilized continuous weight", {
  fx <- joint_wt_fixture()
  joint <- wt_joint(
    fx$w$a,
    fx$w$d,
    exposure_type = c("binary", "continuous")
  )

  expect_equal(
    as.double(joint),
    as.double(fx$w$a) * as.double(fx$w$d),
    tolerance = 1e-12
  )
  expect_identical(estimand(joint), "ate")
  expect_identical(
    joint_wt_meta(joint),
    list(
      exposure_type = c("binary", "continuous"),
      stabilized = c(FALSE, TRUE)
    )
  )
})

test_that("wt_joint() marks the product stabilized when any component is", {
  fx <- joint_wt_fixture()

  # A stabilizing numerator on either factor is a stabilizing numerator on the
  # product, so the flag is the disjunction rather than the conjunction the
  # combining rules use. The record keeps the per-component truth, so the coarse
  # flag hides nothing: a reader who needs to know which factor carried it can
  # still find out.
  neither <- wt_joint(fx$w$a, fx$w$e)
  expect_false(is_stabilized(neither))
  expect_identical(joint_wt_meta(neither)$stabilized, c(FALSE, FALSE))

  one <- wt_joint(
    fx$w$a,
    fx$w$d,
    exposure_type = c("binary", "continuous")
  )
  expect_true(is_stabilized(one))
  expect_identical(joint_wt_meta(one)$stabilized, c(FALSE, TRUE))

  both <- wt_joint(
    fx$w$a_stabilized,
    fx$w$d,
    exposure_type = c("binary", "continuous")
  )
  expect_true(is_stabilized(both))
  expect_identical(joint_wt_meta(both)$stabilized, c(TRUE, TRUE))
})

test_that("wt_joint() accepts two stabilized continuous components", {
  fx <- joint_wt_fixture()

  # A dose crossed with a dose is a joint intervention like any other, and the
  # stabilization requirement is met by both, so nothing here is refused.
  joint <- wt_joint(
    fx$w$d,
    fx$w$d2,
    exposure_type = c("continuous", "continuous")
  )
  expect_true(is_joint_wt(joint))
  expect_identical(estimand(joint), "ate")
  expect_true(is_stabilized(joint))
  expect_identical(
    joint_wt_meta(joint)$exposure_type,
    c("continuous", "continuous")
  )
  expect_equal(
    as.double(joint),
    as.double(fx$w$d) * as.double(fx$w$d2),
    tolerance = 1e-12
  )
})

test_that("the product keeps its joint record through subsetting", {
  fx <- joint_wt_fixture()
  joint <- wt_joint(
    fx$w$a,
    fx$w$d,
    exposure_type = c("binary", "continuous")
  )

  # The record names the two components rather than the observations, so like
  # the categorical attributes it holds at any length and survives every
  # operation that keeps the vector a psw.
  subset <- joint[1:10]
  expect_true(is_psw(subset))
  expect_length(subset, 10L)
  expect_true(is_joint_wt(subset))
  expect_identical(joint_wt_meta(subset), joint_wt_meta(joint))
  expect_identical(estimand(subset), "ate")
  expect_identical(is_stabilized(subset), is_stabilized(joint))
  expect_equal(as.double(subset), as.double(joint)[1:10], tolerance = 1e-12)
})

# ---- wt_joint(): the refusals -----------------------------------------------

test_that("wt_joint() refuses a component that does not target the ate", {
  fx <- joint_wt_fixture()

  # The product of an att weight and an ate weight targets neither the joint ate
  # nor any single population, so there is no estimand to record on the result.
  expect_error(
    wt_joint(fx$w$a_att, fx$w$e),
    class = "propensity_wt_joint_estimand_error"
  )
  expect_error(
    wt_joint(fx$w$a, fx$w$a_att),
    class = "propensity_wt_joint_estimand_error"
  )
  expect_propensity_error(wt_joint(fx$w$a_att, fx$w$e))
})

test_that("wt_joint() refuses an unstabilized continuous component", {
  fx <- joint_wt_fixture()

  # An unstabilized density ratio 1 / f(D | L) has a heavy right tail on its
  # own, and multiplying it by a second weight inherits that tail. A binary
  # component needs no stabilization, so the requirement is on the continuous
  # ones alone rather than on every component.
  expect_error(
    wt_joint(
      fx$w$a,
      fx$w$d_unstabilized,
      exposure_type = c("binary", "continuous")
    ),
    class = "propensity_wt_joint_stabilize_error"
  )
  expect_error(
    wt_joint(
      fx$w$d_unstabilized,
      fx$w$d,
      exposure_type = c("continuous", "continuous")
    ),
    class = "propensity_wt_joint_stabilize_error"
  )

  # The same unstabilized vector as a binary component is accepted, which is
  # what makes this a rule about the exposure type rather than about the flag.
  expect_true(is_joint_wt(wt_joint(fx$w$a, fx$w$e)))

  expect_propensity_error(
    wt_joint(
      fx$w$a,
      fx$w$d_unstabilized,
      exposure_type = c("binary", "continuous")
    )
  )
})

test_that("wt_joint() refuses components of different lengths", {
  fx <- joint_wt_fixture()

  expect_error(
    wt_joint(fx$w$a, fx$w$e[1:10]),
    class = "propensity_wt_joint_length_error"
  )
  expect_propensity_error(wt_joint(fx$w$a, fx$w$e[1:10]))
})

test_that("wt_joint() refuses a component that is not a propensity score weight", {
  fx <- joint_wt_fixture()

  # A bare numeric carries no estimand and no stabilization status, so none of
  # the checks above can be applied to it and the result could record nothing
  # about where it came from.
  expect_error(
    wt_joint(as.double(fx$w$a), fx$w$e),
    class = "propensity_wt_joint_class_error"
  )
  expect_error(
    wt_joint(fx$w$a, as.double(fx$w$e)),
    class = "propensity_wt_joint_class_error"
  )
  expect_propensity_error(wt_joint(as.double(fx$w$a), fx$w$e))
})

test_that("wt_joint() refuses an exposure_type that does not name both components", {
  fx <- joint_wt_fixture()

  # The exposure types are supplied rather than read off the weights: a psw
  # records its estimand and its stabilization but not the kind of exposure it
  # weights, and the stabilization rule above needs to tell a continuous
  # component from a binary one.
  expect_error(
    wt_joint(fx$w$a, fx$w$e, exposure_type = "binary"),
    class = "propensity_wt_joint_exposure_type_error"
  )
  expect_error(
    wt_joint(fx$w$a, fx$w$e, exposure_type = c("binary", "ordinal")),
    class = "propensity_wt_joint_exposure_type_error"
  )
  expect_propensity_error(
    wt_joint(fx$w$a, fx$w$e, exposure_type = c("binary", "ordinal"))
  )
})

# ---- joint_wt_models(): the container ---------------------------------------

test_that("joint_wt_models() records the two treatment models in order", {
  fx <- joint_wt_fixture()
  models <- joint_wt_models(a = fx$mods$a, e = fx$mods$e)

  expect_s3_class(models, "joint_wt_models")
  expect_true(is_joint_wt_models(models))

  # The order is the order the arguments were written in, and it is the
  # factorization's order: the first treatment is the one the second model
  # conditions on.
  expect_identical(models$names, c("a", "e"))
  expect_identical(names(models$models), c("a", "e"))
  expect_identical(models$models$a, fx$mods$a)
  expect_identical(models$models$e, fx$mods$e)
})

test_that("joint_wt_models() records each model's exposure type", {
  skip_if_not_installed("nnet")
  fx <- joint_wt_fixture()
  multinom_mods <- joint_wt_multinom_fixture(fx$dat)

  # The type is read off the model, in the vocabulary the rest of the package
  # uses for exposures, so the estimator does not have to work it out again.
  binary <- joint_wt_models(a = fx$mods$a, e = fx$mods$e)
  expect_identical(binary$exposure_type, c(a = "binary", e = "binary"))

  continuous <- joint_wt_models(a = fx$mods$a, d = fx$mods$d)
  expect_identical(
    continuous$exposure_type,
    c(a = "binary", d = "continuous")
  )

  categorical <- joint_wt_models(a = fx$mods$a, g = multinom_mods$g)
  expect_identical(
    categorical$exposure_type,
    c(a = "binary", g = "categorical")
  )
})

test_that("joint_wt_models() requires a discrete second model to condition on the first treatment", {
  fx <- joint_wt_fixture()

  # The load-bearing guardrail. f(A | L) f(E | L) is a product of two marginal
  # weights and is not the joint weight f(A | L) f(E | A, L) for any data in
  # which E depends on A. Nothing downstream can detect the difference: the
  # product is a perfectly ordinary vector of positive numbers.
  expect_error(
    joint_wt_models(a = fx$mods$a, e = fx$mods$e_marginal),
    class = "propensity_wt_joint_factorization_error"
  )

  # The message has to say what is wrong and what to do. Three things it must
  # carry: that the product of marginal weights is not a joint weight, that the
  # second model must condition on the first treatment, and that the dependence
  # should be modeled flexibly rather than as a single additive term.
  expect_error(
    joint_wt_models(a = fx$mods$a, e = fx$mods$e_marginal),
    regexp = "marginal"
  )
  expect_error(
    joint_wt_models(a = fx$mods$a, e = fx$mods$e_marginal),
    regexp = "condition"
  )
  expect_error(
    joint_wt_models(a = fx$mods$a, e = fx$mods$e_marginal),
    regexp = "flexib"
  )

  # An additive term satisfies the check, since the check is on the
  # factorization rather than on how well the dependence is modeled; the
  # flexibility advice is advice.
  additive <- glm(e ~ a + x1, data = fx$dat, family = binomial())
  expect_true(is_joint_wt_models(
    joint_wt_models(a = fx$mods$a, e = additive)
  ))

  expect_propensity_error(
    joint_wt_models(a = fx$mods$a, e = fx$mods$e_marginal)
  )
})

test_that("joint_wt_models() requires a multinomial second model to condition on the first treatment", {
  skip_if_not_installed("nnet")
  fx <- joint_wt_fixture()
  multinom_mods <- joint_wt_multinom_fixture(fx$dat)

  expect_true(is_joint_wt_models(
    joint_wt_models(a = fx$mods$a, g = multinom_mods$g)
  ))
  expect_error(
    joint_wt_models(a = fx$mods$a, g = multinom_mods$g_marginal),
    class = "propensity_wt_joint_factorization_error"
  )
})

test_that("joint_wt_models() requires a continuous second model to condition on the first treatment", {
  fx <- joint_wt_fixture()

  # The guardrail is about the factorization, not about the response type.
  # f(D | A, L) is what the second factor has to be whether D is a dose or a
  # coin flip, so a dose model that reads the covariate alone is refused on the
  # same terms.
  expect_true(is_joint_wt_models(
    joint_wt_models(a = fx$mods$a, d = fx$mods$d)
  ))
  expect_error(
    joint_wt_models(a = fx$mods$a, d = fx$mods$d_marginal),
    class = "propensity_wt_joint_factorization_error"
  )
})

test_that("joint_wt_models() refuses models that condition on each other", {
  fx <- joint_wt_fixture()

  # A sequential factorization has a first factor and a second one. Two models
  # that each read the other's treatment are not the two factors of any
  # factorization: there is no order in which the first one is marginal in the
  # second treatment. Refused rather than resolved by argument order, since the
  # order the caller wrote is exactly what such a pair contradicts.
  expect_error(
    joint_wt_models(a = fx$mods$a_on_e, e = fx$mods$e),
    class = "propensity_wt_joint_circular_error"
  )
  expect_propensity_error(joint_wt_models(a = fx$mods$a_on_e, e = fx$mods$e))
})

test_that("joint_wt_models() refuses a pair in which neither model conditions on the other", {
  fx <- joint_wt_fixture()

  # Written the other way round, the pair that would be accepted in one order is
  # refused in the other: the first model is the one that has to be marginal in
  # the second treatment, so naming them backwards leaves the second model
  # reading no treatment at all.
  expect_error(
    joint_wt_models(e = fx$mods$e_marginal, a = fx$mods$a),
    class = "propensity_wt_joint_factorization_error"
  )
})

test_that("joint_wt_models() refuses arguments that do not name exactly two treatments", {
  fx <- joint_wt_fixture()

  # The names are the treatment names, so an unnamed argument names no
  # treatment and the container could record none.
  expect_error(
    joint_wt_models(fx$mods$a, e = fx$mods$e),
    class = "propensity_wt_joint_models_error"
  )
  expect_error(
    joint_wt_models(a = fx$mods$a),
    class = "propensity_wt_joint_models_error"
  )
  expect_error(
    joint_wt_models(a = fx$mods$a, e = fx$mods$e, d = fx$mods$d),
    class = "propensity_wt_joint_models_error"
  )
  expect_propensity_error(joint_wt_models(fx$mods$a, e = fx$mods$e))
})

test_that("joint_wt_models() refuses two models named for the same treatment", {
  fx <- joint_wt_fixture()

  # Two factors of one treatment's density is not a crossing, and the container
  # keyed by name could not tell the two apart.
  #
  # The call is assembled rather than written out because the linter refuses a
  # literal duplicated argument name, which is the thing under test. `...`
  # accepts one, so what reaches the constructor is the call a user could write.
  duplicated_names <- list(a = fx$mods$a, a = fx$mods$a)
  expect_error(
    do.call("joint_wt_models", duplicated_names),
    class = "propensity_wt_joint_models_error"
  )
})

test_that("joint_wt_models() refuses a model whose response is not its treatment", {
  fx <- joint_wt_fixture()

  # Each model is the treatment model for the treatment it is named for, so its
  # response has to be that treatment. A mismatch means the container records a
  # crossing of two treatments one of which nothing models, and the
  # factorization check above would then look for the wrong variable.
  expect_error(
    joint_wt_models(a = fx$mods$mislabelled, e = fx$mods$e),
    class = "propensity_wt_joint_response_error"
  )
  expect_propensity_error(
    joint_wt_models(a = fx$mods$mislabelled, e = fx$mods$e)
  )
})
