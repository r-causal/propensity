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
    )),
    # The same dose weighted by a heavier-tailed density and at a spread the
    # caller fixed, which are two of the things the density record carries and
    # which the product has to carry per component.
    d_t = quiet(wt_ate(
      as.double(fitted(mods$d)),
      dat$d,
      exposure_type = "continuous",
      stabilize = TRUE,
      .density = dens_t(4)
    )),
    d_fixed_sigma = quiet(wt_ate(
      as.double(fitted(mods$d)),
      dat$d,
      exposure_type = "continuous",
      stabilize = TRUE,
      .sigma = 1.25
    )),
    # The same dose again, differing from `d` in the numerator of the ratio
    # alone, which is the finest difference the density record carries.
    d_integrated = quiet(wt_ate(
      as.double(fitted(mods$d)),
      dat$d,
      exposure_type = "continuous",
      stabilize = TRUE,
      numerator = "integrated"
    ))
  )
}

# A weight vector carrying an estimand and a stabilization status and nothing
# else, which is what `psw()` builds and what a caller who assembled a weight by
# hand would hold. It records no exposure type, so it is the component
# `wt_joint()` cannot read one off.
recordless_psw <- function(w, stabilized = FALSE) {
  psw(as.double(w), estimand = "ate", stabilized = stabilized)
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
  expect_true(all(lengths(fx$w) == nrow(fx$dat)))
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
      stabilized = c(FALSE, FALSE),
      density = list(NULL, NULL),
      numerator_model = list(NULL, NULL),
      stabilization_score = list(NULL, NULL)
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
      stabilized = c(FALSE, TRUE),
      density = list(NULL, density_meta(fx$w$d)),
      numerator_model = list(NULL, NULL),
      stabilization_score = list(NULL, NULL)
    )
  )
})

test_that("wt_joint() records each component's numerator model", {
  fx <- joint_wt_fixture()

  # A component stabilized on a fitted model records that model, and a component
  # stabilized on anything else records none. The slot is what an estimator
  # reads to estimate the numerator again rather than to take the value it
  # evaluated to as a constant.
  num_a <- glm(a ~ x1, data = fx$dat, family = binomial())
  w_a <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(fx$mods$a, stabilize = num_a)
  )

  joint <- wt_joint(w_a, fx$w$d, exposure_type = c("binary", "continuous"))
  models <- joint_wt_meta(joint)$numerator_model

  expect_length(models, 2L)
  expect_identical(coef(models[[1]]), coef(num_a))
  expect_null(models[[2]])

  # A dose keeps its model in its own density record, where the single
  # continuous route has always kept it, so the slot holds none for it. One
  # model in one place is what makes two copies unable to disagree.
  num_d <- lm(d ~ x1, data = fx$dat)
  w_d <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(fitted(fx$mods$d)),
      fx$dat$d,
      exposure_type = "continuous",
      stabilize = num_d
    )
  )

  dose_joint <- wt_joint(fx$w$a, w_d, exposure_type = c("binary", "continuous"))
  dose_meta <- joint_wt_meta(dose_joint)

  expect_null(dose_meta$numerator_model[[1]])
  expect_null(dose_meta$numerator_model[[2]])
  expect_identical(
    coef(dose_meta$density[[2]]$numerator_model),
    coef(num_d)
  )

  # The reader is what an estimator asks, and it answers for every component
  # alike, falling through to the density record for the ones that have one.
  read <- joint_wt_numerator_models(dose_meta)
  expect_null(read[[1]])
  expect_identical(coef(read[[2]]), coef(num_d))
})

test_that("a record holding a dose's model in the slot still reads", {
  fx <- joint_wt_fixture()

  # A product built before the density record was the canonical home holds the
  # dose's model in the slot and in the record both. The reader takes the slot
  # where it holds one, so such a record answers as it always did.
  num_d <- lm(d ~ x1, data = fx$dat)
  w_d <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(
      as.double(fitted(fx$mods$d)),
      fx$dat$d,
      exposure_type = "continuous",
      stabilize = num_d
    )
  )

  joint <- wt_joint(fx$w$a, w_d, exposure_type = c("binary", "continuous"))
  meta <- joint_wt_meta(joint)
  meta$numerator_model[[2]] <- num_d

  expect_identical(coef(joint_wt_numerator_models(meta)[[2]]), coef(num_d))

  # And one written before the slot existed at all reads the dose's model out
  # of the density record, which is where it has always been.
  meta$numerator_model <- NULL
  expect_identical(coef(joint_wt_numerator_models(meta)[[2]]), coef(num_d))
  expect_null(joint_wt_numerator_models(meta)[[1]])
})

test_that("wt_joint() records each component's stabilization score", {
  fx <- joint_wt_fixture()

  # A numerator the caller computed is a vector rather than a model, and the
  # product records it per component the way it records the model. Without it a
  # discrete component stabilized on a score leaves the record the default
  # stabilizer leaves, and a dose's score is named by the density record and
  # then nowhere to be found.
  score_a <- rep(0.4, nrow(fx$dat))
  score_d <- stats::dnorm(fx$dat$d, mean(fx$dat$d), stats::sd(fx$dat$d))

  quiet <- function(expr) {
    withr::with_options(list(propensity.quiet = TRUE), expr)
  }

  w_a <- quiet(wt_ate(
    fx$mods$a,
    stabilize = TRUE,
    stabilization_score = score_a
  ))
  w_d <- quiet(wt_ate(
    as.double(fitted(fx$mods$d)),
    fx$dat$d,
    exposure_type = "continuous",
    stabilize = TRUE,
    stabilization_score = score_d
  ))

  scores <- joint_wt_meta(
    wt_joint(w_a, w_d, exposure_type = c("binary", "continuous"))
  )$stabilization_score

  expect_length(scores, 2L)
  expect_identical(scores[[1]], score_a)
  expect_identical(scores[[2]], score_d)

  # A component stabilized on anything else records no score, which is what
  # says the numerator is one the estimator can rebuild for itself.
  expect_identical(
    joint_wt_meta(wt_joint(
      fx$w$a,
      fx$w$d,
      c("binary", "continuous")
    ))$stabilization_score,
    list(NULL, NULL)
  )
})

test_that("a length-changing operation drops a per-observation joint score", {
  fx <- joint_wt_fixture()

  # The joint record names the components rather than the units, with the one
  # exception a score is: it holds a value per observation, and a subset it no
  # longer describes would be rebuilt from a numerator belonging to units the
  # result does not hold. The score on the weights themselves is dropped for
  # that reason, and the one inside the record is the same vector.
  score_a <- seq(0.3, 0.5, length.out = nrow(fx$dat))
  w_a <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(fx$mods$a, stabilize = TRUE, stabilization_score = score_a)
  )
  joint <- wt_joint(w_a, fx$w$e)

  expect_identical(joint_wt_meta(joint)$stabilization_score[[1]], score_a)

  subset <- expect_warning(
    joint[1:10],
    class = "propensity_stabilization_score_warning"
  )
  expect_null(joint_wt_meta(subset)$stabilization_score[[1]])

  # A single value scales every weight at any length, so it is carried the way
  # the score on the weights themselves is.
  one <- withr::with_options(
    list(propensity.quiet = TRUE),
    wt_ate(fx$mods$a, stabilize = TRUE, stabilization_score = 0.4)
  )
  short <- expect_silent(wt_joint(one, fx$w$e)[1:10])
  expect_identical(joint_wt_meta(short)$stabilization_score[[1]], 0.4)
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

# ---- wt_joint(): combining product weights ----------------------------------

test_that("the product keeps its joint record through combination", {
  fx <- joint_wt_fixture()
  joint <- wt_joint(fx$w$a, fx$w$e)

  # Appending one set of product weights to another leaves every observation a
  # product of the same two components, so the record describes the result as
  # well as it describes either input and travels with it.
  combined <- expect_silent(c(joint, joint))
  expect_true(is_psw(combined))
  expect_length(combined, 2L * length(joint))
  expect_true(is_joint_wt(combined))
  expect_identical(joint_wt_meta(combined), joint_wt_meta(joint))
  expect_identical(estimand(combined), "ate")
  expect_identical(is_stabilized(combined), is_stabilized(joint))
  expect_equal(
    as.double(combined),
    rep(as.double(joint), 2),
    tolerance = 1e-12
  )

  # `c()` settles its type through `vec_ptype2()`, so the vctrs entry point and
  # the base one answer alike.
  through_vctrs <- expect_silent(vctrs::vec_c(joint, joint))
  expect_true(is_joint_wt(through_vctrs))
  expect_identical(joint_wt_meta(through_vctrs), joint_wt_meta(joint))

  # A column of weights stacked with the rest of its data frame settles the
  # same way, which is how the record reaches a result assembled by row.
  df <- tibble::tibble(id = seq_along(joint), wt = joint)
  stacked <- expect_silent(vctrs::vec_rbind(df, df))
  expect_s3_class(stacked$wt, "psw")
  expect_true(is_joint_wt(stacked$wt))
  expect_identical(joint_wt_meta(stacked$wt), joint_wt_meta(joint))
})

test_that("two products built the same way agree on the record", {
  first <- wt_joint(joint_wt_fixture()$w$a, joint_wt_fixture()$w$d)
  second <- wt_joint(joint_wt_fixture()$w$a, joint_wt_fixture()$w$d)

  # The density slot holds the specification the ratio was read in, whose
  # function is written afresh on every call, so two products built the same way
  # hold two records that say the same thing and are not the same object. A
  # merge comparing them by identity would report a disagreement between weights
  # nothing distinguishes.
  combined <- expect_silent(c(first, second))
  expect_true(is_joint_wt(combined))
  expect_identical(
    joint_wt_meta(combined)$exposure_type,
    c("binary", "continuous")
  )
  expect_identical(joint_wt_meta(combined)$stabilized, c(FALSE, TRUE))
  expect_null(joint_wt_meta(combined)$density[[1]])
  expect_identical(joint_wt_meta(combined)$density[[2]]$numerator, "marginal")
  expect_identical(joint_wt_meta(combined)$density[[2]]$sigma, "pooled")
})

test_that("combining products read in different densities drops the record", {
  fx <- joint_wt_fixture()
  normal <- wt_joint(fx$w$a, fx$w$d)
  heavier <- wt_joint(fx$w$a, fx$w$d_t)

  # The two products weight the same pair of exposures and stabilize the same
  # component, and differ in the density the dose's ratio was read in. That is
  # the ratio an estimator would rebuild the weights from, so neither record
  # describes the combination.
  out <- NULL
  cnd <- expect_warning(
    out <- c(normal, heavier),
    class = "propensity_metadata_conflict_warning"
  )

  expect_s3_class(out, "psw")
  expect_length(out, 2L * length(normal))
  expect_identical(estimand(out), "ate")
  expect_false(is_joint_wt(out))
  expect_null(joint_wt_meta(out))
  expect_true(grepl("joint_wt_meta", conditionMessage(cnd), fixed = TRUE))
})

test_that("combining products over different exposures drops the record", {
  fx <- joint_wt_fixture()
  binary_dose <- wt_joint(fx$w$a, fx$w$d)
  two_doses <- wt_joint(fx$w$d, fx$w$d2)

  # A product of a binary weight and a dose weight is not the same kind of
  # object as a product of two dose weights, and the record is what says which
  # one a set of weights is.
  out <- NULL
  cnd <- expect_warning(
    out <- c(binary_dose, two_doses),
    class = "propensity_metadata_conflict_warning"
  )

  expect_s3_class(out, "psw")
  expect_false(is_joint_wt(out))
  expect_null(joint_wt_meta(out))
  expect_true(grepl("joint_wt_meta", conditionMessage(cnd), fixed = TRUE))
})

test_that("combining products that stabilize different components drops the record", {
  fx <- joint_wt_fixture()
  plain <- wt_joint(fx$w$a, fx$w$d)
  both <- wt_joint(fx$w$a_stabilized, fx$w$d)

  # The two products weight the same pair of exposures with the same density,
  # and differ in whether the binary component's weight is stabilized. The
  # product is marked stabilized when any component is, so both are, and the
  # per-component record is the only thing that says which of the two the
  # stabilization belongs to.
  expect_true(is_stabilized(plain))
  expect_true(is_stabilized(both))
  expect_identical(joint_wt_meta(plain)$stabilized, c(FALSE, TRUE))
  expect_identical(joint_wt_meta(both)$stabilized, c(TRUE, TRUE))

  out <- NULL
  cnd <- expect_warning(
    out <- c(plain, both),
    class = "propensity_metadata_conflict_warning"
  )

  expect_s3_class(out, "psw")
  expect_length(out, 2L * length(plain))
  expect_identical(estimand(out), "ate")
  expect_false(is_joint_wt(out))
  expect_null(joint_wt_meta(out))
  expect_true(grepl("joint_wt_meta", conditionMessage(cnd), fixed = TRUE))
})

test_that("combining products read with different numerators drops the record", {
  fx <- joint_wt_fixture()
  marginal <- wt_joint(fx$w$a, fx$w$d)
  integrated <- wt_joint(fx$w$a, fx$w$d_integrated)

  # The finest difference the density record carries: the same dose, the same
  # density family, the same spread, and a numerator read one way rather than
  # the other. It is still the ratio an estimator would rebuild the weights
  # from, so the record does not survive the disagreement.
  expect_identical(
    joint_wt_meta(marginal)$density[[2]]$numerator,
    "marginal"
  )
  expect_identical(
    joint_wt_meta(integrated)$density[[2]]$numerator,
    "integrated"
  )
  expect_identical(
    joint_wt_meta(marginal)$stabilized,
    joint_wt_meta(integrated)$stabilized
  )

  out <- NULL
  cnd <- expect_warning(
    out <- c(marginal, integrated),
    class = "propensity_metadata_conflict_warning"
  )

  expect_s3_class(out, "psw")
  expect_length(out, 2L * length(marginal))
  expect_false(is_joint_wt(out))
  expect_null(joint_wt_meta(out))
  expect_true(grepl("joint_wt_meta", conditionMessage(cnd), fixed = TRUE))
})

test_that("the conflict warning names the record it dropped", {
  fx <- joint_wt_fixture()
  normal <- wt_joint(fx$w$a, fx$w$d)
  heavier <- wt_joint(fx$w$a, fx$w$d_t)

  # The warning is the whole account a caller gets of a record the result no
  # longer carries, so it names that record by the name the accessor reads it
  # under.
  out <- expect_propensity_warning(c(normal, heavier))
  expect_false(is_joint_wt(out))
})

test_that("combining a product with a plain psw carries the record", {
  fx <- joint_wt_fixture()
  joint <- wt_joint(fx$w$a, fx$w$e)
  plain <- recordless_psw(fx$w$a)

  # Weights recording nothing about a product have nothing to disagree with,
  # which is the rule every carried attribute follows: the one operand that
  # records the thing speaks for the result. A disagreement takes two records
  # that say different things.
  combined <- expect_silent(c(joint, plain))
  expect_true(is_joint_wt(combined))
  expect_identical(joint_wt_meta(combined), joint_wt_meta(joint))

  # The rule does not depend on the order the two were written in.
  reversed <- expect_silent(c(plain, joint))
  expect_true(is_joint_wt(reversed))
  expect_identical(joint_wt_meta(reversed), joint_wt_meta(joint))
})

# ---- wt_joint(): what the components record ---------------------------------

test_that("wt_joint() reads each component's recorded exposure type", {
  fx <- joint_wt_fixture()

  # A psw built by a weight function records the exposure type it was built
  # for, so the caller no longer has to repeat it. What the product records is
  # what the components record, in the order they were given.
  joint <- wt_joint(fx$w$a, fx$w$d)
  expect_identical(
    joint_wt_meta(joint)$exposure_type,
    c("binary", "continuous")
  )
  expect_identical(joint_wt_meta(joint)$stabilized, c(FALSE, TRUE))
  expect_true(is_joint_wt(joint))
  expect_equal(
    as.double(joint),
    as.double(fx$w$a) * as.double(fx$w$d),
    tolerance = 1e-12
  )

  # The record is per component, so the joint record is the pair and the
  # product itself records no single type: it weights two exposures.
  expect_null(exposure_type(joint))

  two_binary <- wt_joint(fx$w$a, fx$w$e)
  expect_identical(
    joint_wt_meta(two_binary)$exposure_type,
    c("binary", "binary")
  )

  two_doses <- wt_joint(fx$w$d, fx$w$d2)
  expect_identical(
    joint_wt_meta(two_doses)$exposure_type,
    c("continuous", "continuous")
  )

  # The stabilization requirement is decided from the types that were read, so
  # a dose component that is not stabilized is refused without the caller
  # having said it was a dose.
  expect_error(
    wt_joint(fx$w$a, fx$w$d_unstabilized),
    class = "propensity_wt_joint_stabilize_error"
  )
})

test_that("wt_joint() refuses a component that records no exposure type", {
  fx <- joint_wt_fixture()
  hand_built <- recordless_psw(fx$w$e)

  # A weight assembled by hand records an estimand and a stabilization status
  # and nothing about the exposure, so there is no type to read and the
  # stabilization requirement could not be applied to it.
  expect_error(
    wt_joint(fx$w$a, hand_built),
    class = "propensity_wt_joint_exposure_type_error"
  )

  # The refusal names the component it could not read, since the remedy is to
  # say what that one weights.
  expect_error(
    wt_joint(hand_built, fx$w$e),
    class = "propensity_wt_joint_exposure_type_error",
    regexp = "w_a"
  )
  expect_error(
    wt_joint(fx$w$a, hand_built),
    class = "propensity_wt_joint_exposure_type_error",
    regexp = "w_e"
  )

  expect_propensity_error(wt_joint(fx$w$a, hand_built))

  # Saying it is what makes the same pair acceptable.
  named <- wt_joint(fx$w$a, hand_built, exposure_type = c("binary", "binary"))
  expect_true(is_joint_wt(named))
  expect_identical(
    joint_wt_meta(named)$exposure_type,
    c("binary", "binary")
  )
})

test_that("an explicit exposure_type wins over what the components record", {
  fx <- joint_wt_fixture()

  # Read off the record, the unstabilized dose is a dose and is refused. Named
  # as a binary component it is accepted, and the record the product carries is
  # the one the caller named: the argument decides, and the reading is what
  # happens when there is no argument.
  expect_error(
    wt_joint(fx$w$a, fx$w$d_unstabilized),
    class = "propensity_wt_joint_stabilize_error"
  )
  declared <- wt_joint(
    fx$w$a,
    fx$w$d_unstabilized,
    exposure_type = c("binary", "binary")
  )
  expect_identical(
    joint_wt_meta(declared)$exposure_type,
    c("binary", "binary")
  )

  # An explicit value is validated exactly as before.
  expect_error(
    wt_joint(fx$w$a, fx$w$e, exposure_type = c("binary", "ordinal")),
    class = "propensity_wt_joint_exposure_type_error"
  )
  expect_error(
    wt_joint(fx$w$a, fx$w$e, exposure_type = "binary"),
    class = "propensity_wt_joint_exposure_type_error"
  )
})

test_that("the product records each component's density", {
  fx <- joint_wt_fixture()

  # A dose's weights are a ratio of densities, and which ratio they are is what
  # the estimator has to rebuild. The product keeps each component's record in
  # the order the components were given, so the dose's family, numerator, and
  # spread survive the multiplication.
  joint <- wt_joint(fx$w$a, fx$w$d)
  density <- joint_wt_meta(joint)$density

  expect_length(density, 2L)
  expect_null(density[[1]])
  expect_s3_class(density[[2]], "propensity_density_meta")
  expect_true(density_specs_agree(density[[2]]$density, dens_normal()))
  expect_identical(density[[2]]$numerator, "marginal")
  expect_identical(density[[2]]$sigma, "pooled")
  expect_null(density[[2]]$sigma_value)

  # A binary component weights no density, so its slot is empty rather than
  # absent: the record is one element per component either way.
  expect_identical(
    joint_wt_meta(wt_joint(fx$w$a, fx$w$e))$density,
    list(NULL, NULL)
  )

  # The two choices that change the ratio, each read back off the product.
  heavy <- wt_joint(fx$w$a, fx$w$d_t)
  expect_true(density_specs_agree(
    joint_wt_meta(heavy)$density[[2]]$density,
    dens_t(4)
  ))

  fixed <- wt_joint(fx$w$a, fx$w$d_fixed_sigma)
  expect_identical(joint_wt_meta(fixed)$density[[2]]$sigma, "supplied")
  expect_equal(
    joint_wt_meta(fixed)$density[[2]]$sigma_value,
    1.25,
    tolerance = 1e-12
  )

  # A dose in the first position records in the first slot, so the pair is
  # positional rather than a single record the product happens to carry.
  doses <- wt_joint(fx$w$d, fx$w$d2)
  expect_identical(
    joint_wt_meta(doses)$density,
    list(density_meta(fx$w$d), density_meta(fx$w$d2))
  )
})

test_that("multiplying two psw vectors does not build a joint weight", {
  fx <- joint_wt_fixture()

  # `wt_joint()` is how a joint weight is built. Arithmetic on two psw vectors
  # gives their product and none of the record, which is the whole reason the
  # constructor exists: the checks it applies cannot be applied afterwards.
  expect_warning(
    product <- fx$w$a * fx$w$d,
    class = "propensity_metadata_conflict_warning"
  )

  expect_true(is_psw(product))
  expect_false(is_joint_wt(product))
  expect_null(joint_wt_meta(product))

  # The two components record different exposure types and only one of them
  # records a density, so neither describes the product and the combining rules
  # drop both.
  expect_null(exposure_type(product))
  expect_null(density_meta(product))

  # The values are the same numbers `wt_joint()` would give, which is what
  # makes the missing record the only difference and the reason nothing
  # downstream could notice it.
  expect_equal(
    as.double(product),
    as.double(wt_joint(fx$w$a, fx$w$d)),
    tolerance = 1e-12
  )
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

  # A value supplied here is used in place of what the components record, so it
  # is validated on its own terms: the stabilization rule above needs to tell a
  # continuous component from a binary one, and a vector naming one type, or a
  # type this package has no weights for, tells it neither.
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

test_that("the unsupported-type refusal is causalgenerics' refusal wrapped", {
  fx <- joint_wt_fixture()

  # Exposure-type vocabulary belongs to causalgenerics, which owns the list of
  # types and refuses one it does not know. This function keeps its own class,
  # because the same class also answers a vector of the wrong length and a
  # component that records no type at all, neither of which causalgenerics has
  # anything to say about. The refusal it raised is carried as the parent, so
  # the reader is told what was wrong with the type as well as which argument
  # was wrong.
  err <- expect_error(
    wt_joint(fx$w$a, fx$w$e, exposure_type = c("binary", "ordinal")),
    class = "propensity_wt_joint_exposure_type_error"
  )
  expect_s3_class(err, "propensity_error")
  expect_s3_class(err$parent, "condition")

  # The parent is raised outside this package: causalgenerics' own condition
  # where it recognizes the type, and the `arg_match()` refusal it answers a
  # type it has never heard of with. Either way it names the type, which is the
  # half of the report this package has nothing to say about.
  expect_false(inherits(err$parent, "propensity_error"))
  expect_true(
    inherits(err$parent, "causalgenerics_error") ||
      inherits(err$parent, "rlang_error")
  )
  expect_match(conditionMessage(err$parent), "ordinal", fixed = TRUE)

  # The length of the vector is this function's own requirement, and the refusal
  # for it stays its own: there is no single type for causalgenerics to read.
  short <- expect_error(
    wt_joint(fx$w$a, fx$w$e, exposure_type = "binary"),
    class = "propensity_wt_joint_exposure_type_error"
  )
  expect_null(short$parent)

  # A pair of types the package supports still passes.
  expect_s3_class(
    wt_joint(fx$w$a, fx$w$e, exposure_type = c("binary", "binary")),
    "psw"
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

test_that("is_joint_wt_models() is FALSE for anything else", {
  fx <- joint_wt_fixture()

  # The predicate answers for whatever it is handed, so a caller can branch on
  # it without checking the class first.
  expect_false(is_joint_wt_models(fx$mods$a))
  expect_false(is_joint_wt_models(list(a = fx$mods$a, e = fx$mods$e)))
  expect_false(is_joint_wt_models(fx$w$a))
  expect_false(is_joint_wt_models(NULL))
  expect_false(is_joint_wt_models("a"))
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

test_that("joint_wt_models() reads a dose from a robust or an additive model", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("mgcv")
  fx <- joint_wt_fixture()

  # The classes a dose density can be read from are named rather than reached
  # by inheritance: an `rlm` is an `lm` and a `gam` is a `glm`, so both would
  # otherwise be typed by the branch written for the class they inherit from,
  # and a class the package cannot read would be typed by it too.
  robust <- MASS::rlm(d ~ a * x1, data = fx$dat, acc = 1e-10)
  additive <- mgcv::gam(d ~ a + s(x1), data = fx$dat)

  expect_identical(joint_wt_model_type(robust), "continuous")
  expect_identical(joint_wt_model_type(additive), "continuous")

  expect_identical(
    joint_wt_models(a = fx$mods$a, d = robust)$exposure_type,
    c(a = "binary", d = "continuous")
  )
  expect_identical(
    joint_wt_models(a = fx$mods$a, d = additive)$exposure_type,
    c(a = "binary", d = "continuous")
  )
})

test_that("the unsupported-model refusal names the dose models it reads", {
  fx <- joint_wt_fixture()
  withr::local_seed(5302)
  dat <- fx$dat
  dat$k <- rpois(nrow(dat), 2)

  # A count model is not a treatment model this package reads a density from,
  # and the refusal is where a caller learns which classes it does read. The
  # robust and additive fits belong on that list, since both are accepted.
  counted <- glm(k ~ a + x1, data = dat, family = poisson())
  err <- expect_error(
    joint_wt_models(a = fx$mods$a, k = counted),
    class = "propensity_wt_joint_models_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "rlm", fixed = TRUE)
  expect_match(msg, "gam", fixed = TRUE)

  # What is refused is the family rather than the class: a `glm` is on the list
  # of classes the container reads, and every model it rejects here is a class
  # that is on that list fit to a family that is not. Naming only the class
  # leaves the reader looking at a class the same message says is supported.
  expect_match(msg, "poisson()", fixed = TRUE)
})

test_that("the unsupported-model refusal names the family of an additive fit", {
  skip_if_not_installed("mgcv")
  fx <- joint_wt_fixture()
  withr::local_seed(5302)
  dat <- fx$dat
  dat$k <- rpois(nrow(dat), 2)

  # An additive model is read as a dose model when it is gaussian and refused
  # when it is not, and it is the family that decides, so it is the family the
  # refusal names.
  counted <- mgcv::gam(k ~ a + s(x1), data = dat, family = poisson())
  err <- expect_error(
    joint_wt_models(a = fx$mods$a, k = counted),
    class = "propensity_wt_joint_models_error"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(msg, "poisson()", fixed = TRUE)

  expect_propensity_error(joint_wt_models(a = fx$mods$a, k = counted))
})

test_that("the unsupported-model refusal survives an argument that is no model", {
  fx <- joint_wt_fixture()

  # An argument that is not a fitted model at all has no family to read, and
  # reading one anyway raises a base error about subscripts in place of the
  # refusal that would have named what arrived.
  err <- expect_error(
    joint_wt_models(a = fx$mods$a, k = "not a model"),
    class = "propensity_wt_joint_models_error"
  )
  expect_match(conditionMessage(err), "character", fixed = TRUE)

  err_fn <- expect_error(
    joint_wt_models(a = fx$mods$a, k = function(x) x),
    class = "propensity_wt_joint_models_error"
  )
  expect_match(conditionMessage(err_fn), "function", fixed = TRUE)

  # A data frame is a list, so the guard cannot be a test of that alone: `[[`
  # asks it for a column of that name and has none to return.
  err_df <- expect_error(
    joint_wt_models(a = fx$mods$a, k = fx$dat),
    class = "propensity_wt_joint_models_error"
  )
  expect_match(conditionMessage(err_df), "data.frame", fixed = TRUE)
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

  # What the check reads is the variables the right-hand side mentions, not the
  # term labels it produces. The two agree wherever the treatment is written as
  # a bare main effect, which is every spelling above, so these two are what
  # separate them.
  #
  # A transformed treatment produces the term labels `factor(a)`, `x1`, and
  # `factor(a):x1`, none of which is `a`, and an interaction-only dependence
  # produces `x1` and `x1:a`, neither of which is `a` either. Both models
  # condition on the first treatment, and a check written against term labels
  # would refuse them.
  transformed <- glm(e ~ factor(a) * x1, data = fx$dat, family = binomial())
  expect_false("a" %in% attr(stats::terms(transformed), "term.labels"))
  expect_true(is_joint_wt_models(
    joint_wt_models(a = fx$mods$a, e = transformed)
  ))

  interaction_only <- glm(e ~ x1 + a:x1, data = fx$dat, family = binomial())
  expect_false("a" %in% attr(stats::terms(interaction_only), "term.labels"))
  expect_true(is_joint_wt_models(
    joint_wt_models(a = fx$mods$a, e = interaction_only)
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

test_that("the factorization refusal offers the swapped order when the pair is one reversed", {
  fx <- joint_wt_fixture()

  # `a_on_e` reads "e" and `e_marginal` does not read "a", so this pair is a
  # factorization written backwards rather than a pair of marginal models. The
  # refusal stands, since the order the caller wrote is not a factorization, but
  # adding "a" to the second model is not the only way out: the same two fits
  # supplied the other way round are f(e | L) f(a | e, L) and need no refitting.
  # Which of the two is right is a modelling choice about the order the
  # treatments are assigned in, so the message offers it rather than deciding.
  expect_error(
    joint_wt_models(a = fx$mods$a_on_e, e = fx$mods$e_marginal),
    class = "propensity_wt_joint_factorization_error"
  )
  expect_error(
    joint_wt_models(a = fx$mods$a_on_e, e = fx$mods$e_marginal),
    regexp = "other order"
  )

  # The swapped pair really is accepted, which is what makes the advice advice.
  expect_true(is_joint_wt_models(
    joint_wt_models(e = fx$mods$e_marginal, a = fx$mods$a_on_e)
  ))

  # The other quadrant keeps the message it has: neither model reads the other's
  # treatment, so there is no order in which this pair is a factorization and
  # nothing to offer beyond refitting.
  plain <- expect_error(
    joint_wt_models(a = fx$mods$a, e = fx$mods$e_marginal),
    class = "propensity_wt_joint_factorization_error"
  )
  expect_false(grepl("other order", conditionMessage(plain), fixed = TRUE))

  expect_propensity_error(
    joint_wt_models(a = fx$mods$a_on_e, e = fx$mods$e_marginal)
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

  # A call naming none of its arguments carries no name vector at all rather
  # than one of empty strings, so the count of unnamed arguments has to come
  # from the models. Counted off the names alone it collapses to zero by
  # zero-length recycling, and the refusal reports that none of them is the
  # thing it is refusing them for.
  expect_error(
    joint_wt_models(fx$mods$a, fx$mods$e),
    class = "propensity_wt_joint_models_error"
  )
  expect_error(
    joint_wt_models(fx$mods$a, fx$mods$e),
    class = "propensity_error"
  )
  expect_error(
    joint_wt_models(fx$mods$a, fx$mods$e),
    regexp = "2 of them are unnamed"
  )
})

test_that("joint_wt_models() refuses two models named for the same treatment", {
  fx <- joint_wt_fixture()

  # Two factors of one treatment's density is not a crossing, and the container
  # keyed by name could not tell the two apart.
  #
  # The names are attached after the list is built rather than written into it,
  # because the linter refuses a literal duplicated argument name, which is the
  # thing under test. `...` accepts one, so what reaches the constructor is the
  # call a user could write.
  duplicated_names <- stats::setNames(
    list(fx$mods$a, fx$mods$a),
    c("a", "a")
  )
  expect_identical(names(duplicated_names), c("a", "a"))
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
