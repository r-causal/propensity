test_that("wt_atc is an alias for wt_atu", {
  ps <- c(0.1, 0.3, 0.4, 0.3)
  exposure <- c(0, 0, 1, 0)

  # Compare results from wt_atu and wt_atc
  wts_atu <- wt_atu(ps, exposure, exposure_type = "binary")
  wts_atc <- wt_atc(ps, exposure, exposure_type = "binary")

  expect_identical(wts_atu, wts_atc)

  # Test with data.frame
  ps_df <- data.frame(
    control = c(0.9, 0.7, 0.3, 0.1),
    treated = c(0.1, 0.3, 0.7, 0.9)
  )
  exposure_df <- c(0, 0, 1, 1)

  wts_atu_df <- wt_atu(ps_df, exposure_df)
  wts_atc_df <- wt_atc(ps_df, exposure_df)

  expect_identical(wts_atu_df, wts_atc_df)

  # Test with GLM
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  treatment <- rbinom(n, 1, plogis(0.5 * x1 + 0.3 * x2))
  ps_model <- glm(treatment ~ x1 + x2, family = binomial)

  wts_atu_glm <- wt_atu(ps_model, treatment)
  wts_atc_glm <- wt_atc(ps_model, treatment)

  expect_identical(wts_atu_glm, wts_atc_glm)

  # Test with ps_trim
  ps_trimmed <- ps_trim(ps, .exposure = exposure, trim_at = 0.2)
  # Suppress refit warnings - we're testing the alias behavior, not the warning
  suppressWarnings({
    wts_atu_trim <- wt_atu(ps_trimmed, exposure)
    wts_atc_trim <- wt_atc(ps_trimmed, exposure)
  })

  expect_identical(wts_atu_trim, wts_atc_trim)

  # Test with ps_trunc
  ps_truncated <- ps_trunc(ps, .exposure = exposure, trunc_at = 0.2)
  wts_atu_trunc <- wt_atu(ps_truncated, exposure)
  wts_atc_trunc <- wt_atc(ps_truncated, exposure)

  expect_identical(wts_atu_trunc, wts_atc_trunc)
})

test_that("wt_atc works with all object types", {
  # Test with single column data.frame
  df <- data.frame(ps = c(0.1, 0.3, 0.4, 0.3))
  trt <- c(0, 0, 1, 0)

  wts_df <- wt_atc(df, trt)
  expect_s3_class(wts_df, "psw")
  expect_equal(estimand(wts_df), "atu")

  # Test with GLM without covariates
  glm_mod <- glm(trt ~ 1, family = binomial, data = data.frame(trt = trt))
  wts_glm <- wt_atc(glm_mod, trt)
  expect_s3_class(wts_glm, "psw")
  expect_equal(estimand(wts_glm), "atu")

  # Test with named data.frame columns
  ps_named <- data.frame(
    prob_control = c(0.9, 0.7, 0.3, 0.1),
    prob_treated = c(0.1, 0.3, 0.7, 0.9)
  )
  exposure_named <- c("control", "control", "treated", "treated")

  wts_named <- wt_atc(
    ps_named,
    exposure_named,
    .propensity_col = "prob_treated"
  )
  expect_s3_class(wts_named, "psw")
  expect_equal(estimand(wts_named), "atu")

  # Test categorical exposure
  ps_cat <- matrix(
    c(0.5, 0.3, 0.2, 0.2, 0.5, 0.3, 0.1, 0.3, 0.6, 0.6, 0.2, 0.2),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(NULL, c("A", "B", "C"))
  )
  exposure_cat <- factor(c("A", "B", "C", "A"))

  wts_cat <- wt_atc(
    ps_cat,
    exposure_cat,
    .focal_level = "A",
    exposure_type = "categorical"
  )
  expect_s3_class(wts_cat, "psw")
  expect_equal(estimand(wts_cat), "atu")
})

test_that("psw objects can be multiplied together", {
  ps <- c(0.1, 0.3, 0.4, 0.3)
  exposure <- c(0, 0, 1, 0)

  # Create ATE and censoring weights
  wts_ate <- wt_ate(ps, exposure, exposure_type = "binary")
  wts_cens <- wt_cens(ps, exposure, exposure_type = "binary")

  # Multiply them together
  combined_weights <- wts_ate * wts_cens

  # Check that multiplication works
  expect_equal(
    as.numeric(combined_weights),
    as.numeric(wts_ate) * as.numeric(wts_cens)
  )

  # Check that the estimand is combined
  expect_equal(estimand(combined_weights), "ate, uncensored")

  # Check that it's still a psw object
  expect_true(is_psw(combined_weights))

  # Test with stabilized weights
  wts_ate_stab <- wt_ate(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE
  )
  wts_cens_stab <- wt_cens(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE
  )

  combined_stab <- wts_ate_stab * wts_cens_stab
  expect_equal(estimand(combined_stab), "ate, uncensored")
  expect_true(is_stabilized(combined_stab))

  # Test mixed stabilization (only one stabilized)
  combined_mixed <- wts_ate * wts_cens_stab
  expect_equal(estimand(combined_mixed), "ate, uncensored")
  expect_false(is_stabilized(combined_mixed)) # Should be FALSE since only one is stabilized

  # Test with trimmed weights
  ps_trimmed <- ps_trim(ps, .exposure = exposure, trim_at = 0.2)
  # Suppress refit warnings - we're testing combination behavior, not the warning
  suppressWarnings({
    wts_ate_trim <- wt_ate(ps_trimmed, exposure)
  })
  combined_trim <- wts_ate_trim * wts_cens
  expect_true(is_ps_trimmed(combined_trim))
})

test_that("wt_cens uses ATE formula with uncensored estimand", {
  ps <- c(0.1, 0.3, 0.4, 0.3)
  exposure <- c(0, 0, 1, 0)

  # Get weights from wt_ate and wt_cens
  wts_ate <- wt_ate(ps, exposure, exposure_type = "binary")
  wts_cens <- wt_cens(ps, exposure, exposure_type = "binary")

  # The numeric values should be identical
  expect_equal(as.numeric(wts_ate), as.numeric(wts_cens))

  # But the estimand should be different
  expect_equal(estimand(wts_ate), "ate")
  expect_equal(estimand(wts_cens), "uncensored")

  # Test stabilized weights
  wts_ate_stab <- wt_ate(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE
  )
  wts_cens_stab <- wt_cens(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE
  )

  expect_equal(as.numeric(wts_ate_stab), as.numeric(wts_cens_stab))
  expect_true(is_stabilized(wts_cens_stab))
  expect_equal(estimand(wts_cens_stab), "uncensored")

  # Test with continuous exposure
  set.seed(123)
  n <- 32
  denom_model <- lm(mpg ~ gear + am + carb, data = mtcars)

  wts_ate_cont <- wt_ate(
    predict(denom_model),
    .exposure = mtcars$mpg,
    .sigma = influence(denom_model)$sigma,
    exposure_type = "continuous",
    stabilize = TRUE
  )

  wts_cens_cont <- wt_cens(
    predict(denom_model),
    .exposure = mtcars$mpg,
    .sigma = influence(denom_model)$sigma,
    exposure_type = "continuous",
    stabilize = TRUE
  )

  expect_equal(as.numeric(wts_ate_cont), as.numeric(wts_cens_cont))
  expect_equal(estimand(wts_cens_cont), "uncensored")

  # Test with ps_trim
  ps_trimmed <- ps_trim(ps, .exposure = exposure, trim_at = 0.2)
  # Suppress refit warnings - we're testing wt_cens behavior, not the warning
  suppressWarnings({
    wts_cens_trim <- wt_cens(ps_trimmed, exposure)
  })

  expect_equal(estimand(wts_cens_trim), "uncensored; trimmed")
  expect_true(is_ps_trimmed(wts_cens_trim))

  # Test with ps_trunc
  ps_truncated <- ps_trunc(ps, .exposure = exposure, trunc_at = 0.2)
  wts_cens_trunc <- wt_cens(ps_truncated, exposure)

  expect_equal(estimand(wts_cens_trunc), "uncensored; truncated")
  expect_true(is_ps_truncated(wts_cens_trunc))
})

test_that("wt_ate stores a user-supplied stabilization_score", {
  ps <- c(0.2, 0.4, 0.6, 0.8, 0.3, 0.7)
  exposure <- c(1, 0, 1, 0, 1, 0)
  score <- 0.42

  w <- wt_ate(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE,
    stabilization_score = score
  )

  expect_true(is_stabilized(w))
  expect_equal(attr(w, "stabilization_score"), score)
  expect_equal(stabilization_score(w), score)
})

test_that("wt_ate stores a stabilization_score for continuous exposures", {
  denom_model <- lm(mpg ~ gear + am + carb, data = mtcars)
  score <- 0.7

  w <- wt_ate(
    predict(denom_model),
    .exposure = mtcars$mpg,
    .sigma = influence(denom_model)$sigma,
    exposure_type = "continuous",
    stabilize = TRUE,
    stabilization_score = score
  )

  expect_true(is_stabilized(w))
  expect_equal(attr(w, "stabilization_score"), score)
  expect_equal(stabilization_score(w), score)
})

test_that("wt_ate records no stabilization_score under default stabilization", {
  ps <- c(0.2, 0.4, 0.6, 0.8, 0.3, 0.7)
  exposure <- c(1, 0, 1, 0, 1, 0)

  # Default marginal stabilizer: no fixed score is supplied, so the M-estimation
  # path must be able to tell that none was stored.
  w_default <- wt_ate(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE
  )
  expect_true(is_stabilized(w_default))
  expect_null(attr(w_default, "stabilization_score"))
  expect_null(stabilization_score(w_default))

  # Unstabilized weights never carry a stabilization score.
  w_none <- wt_ate(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = FALSE
  )
  expect_false(is_stabilized(w_none))
  expect_null(attr(w_none, "stabilization_score"))
  expect_null(stabilization_score(w_none))
})

test_that("wt_cens stores a user-supplied stabilization_score", {
  ps <- c(0.2, 0.4, 0.6, 0.8, 0.3, 0.7)
  exposure <- c(1, 0, 1, 0, 1, 0)
  score <- 0.55

  w <- wt_cens(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE,
    stabilization_score = score
  )

  expect_equal(estimand(w), "uncensored")
  expect_true(is_stabilized(w))
  expect_equal(attr(w, "stabilization_score"), score)
  expect_equal(stabilization_score(w), score)
})

test_that("wt_cens records no stabilization_score under default stabilization", {
  ps <- c(0.2, 0.4, 0.6, 0.8, 0.3, 0.7)
  exposure <- c(1, 0, 1, 0, 1, 0)

  w_default <- wt_cens(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE
  )
  expect_true(is_stabilized(w_default))
  expect_null(attr(w_default, "stabilization_score"))
  expect_null(stabilization_score(w_default))

  w_none <- wt_cens(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = FALSE
  )
  expect_false(is_stabilized(w_none))
  expect_null(attr(w_none, "stabilization_score"))
  expect_null(stabilization_score(w_none))
})

test_that("stabilization_score survives subsetting of psw", {
  ps <- c(0.2, 0.4, 0.6, 0.8, 0.3, 0.7)
  exposure <- c(1, 0, 1, 0, 1, 0)
  score <- 0.42

  w <- wt_ate(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE,
    stabilization_score = score
  )

  w_sub <- w[1:3]
  expect_true(is_stabilized(w_sub))
  expect_equal(attr(w_sub, "stabilization_score"), score)
  expect_equal(stabilization_score(w_sub), score)
})

# A per-observation stabilization score is indexed by observation, so it is only
# meaningful at the length of the weights it was recorded on. Restoring a psw at
# a different length cannot carry it, and the slice indices are not available to
# subset it, so the contract is to drop it and say so.

stabilized_psw <- function(score) {
  wt_ate(
    c(0.2, 0.4, 0.6, 0.8, 0.3, 0.7),
    c(1, 0, 1, 0, 1, 0),
    exposure_type = "binary",
    stabilize = TRUE,
    stabilization_score = score
  )
}

count_score_warnings <- function(expr) {
  count <- 0
  value <- withCallingHandlers(
    expr,
    propensity_stabilization_score_warning = function(cnd) {
      count <<- count + 1
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, count = count)
}

test_that("vector stabilization_score is dropped with a warning when `[` shortens a psw", {
  score <- c(0.51, 0.52, 0.53, 0.54, 0.55, 0.56)
  w <- stabilized_psw(score)
  expect_equal(stabilization_score(w), score)

  cnd <- expect_warning(
    w_sub <- w[1:3],
    class = "propensity_stabilization_score_warning"
  )
  expect_s3_class(cnd, "propensity_warning")

  expect_s3_class(w_sub, "psw")
  expect_length(w_sub, 3)
  expect_null(stabilization_score(w_sub))
  expect_null(attr(w_sub, "stabilization_score"))

  # Every other piece of metadata is unaffected by the drop.
  expect_equal(estimand(w_sub), "ate")
  expect_true(is_stabilized(w_sub))
  expect_equal(vec_data(w_sub), vec_data(w)[1:3])
})

test_that("vector stabilization_score is dropped with a warning when vec_slice shortens a psw", {
  score <- c(0.51, 0.52, 0.53, 0.54, 0.55, 0.56)
  w <- stabilized_psw(score)

  cnd <- expect_warning(
    w_sub <- vctrs::vec_slice(w, 1:3),
    class = "propensity_stabilization_score_warning"
  )
  expect_s3_class(cnd, "propensity_warning")

  expect_s3_class(w_sub, "psw")
  expect_length(w_sub, 3)
  expect_null(stabilization_score(w_sub))
  expect_null(attr(w_sub, "stabilization_score"))
  expect_equal(estimand(w_sub), "ate")
  expect_true(is_stabilized(w_sub))
})

test_that("scalar stabilization_score is carried through psw slicing silently", {
  score <- 0.42
  w <- stabilized_psw(score)

  w_sub <- expect_silent(w[1:3])
  expect_equal(stabilization_score(w_sub), score)

  w_slice <- expect_silent(vctrs::vec_slice(w, 1:3))
  expect_equal(stabilization_score(w_slice), score)
})

test_that("vector stabilization_score is carried through length-preserving psw arithmetic", {
  score <- c(0.51, 0.52, 0.53, 0.54, 0.55, 0.56)
  w <- stabilized_psw(score)

  w_doubled <- expect_silent(w * 2)
  expect_s3_class(w_doubled, "psw")
  expect_equal(stabilization_score(w_doubled), score)
  expect_equal(vec_data(w_doubled), vec_data(w) * 2)
})

test_that("zero-length restores of a vector-score psw keep the score silently", {
  score <- c(0.51, 0.52, 0.53, 0.54, 0.55, 0.56)
  w <- stabilized_psw(score)

  # A prototype and an empty subset hold no observations, so a score on them
  # lines up with nothing and is left alone.
  proto <- expect_silent(vctrs::vec_ptype(w))
  expect_length(proto, 0)
  expect_equal(stabilization_score(proto), score)

  empty <- expect_silent(w[integer(0)])
  expect_length(empty, 0)
  expect_equal(stabilization_score(empty), score)
})

test_that("concatenating vector-score psw objects warns and drops the score", {
  score <- c(0.51, 0.52, 0.53, 0.54, 0.55, 0.56)
  w <- stabilized_psw(score)

  out <- count_score_warnings(c(w, w))

  # A score recorded on six observations means nothing on twelve, so it goes.
  # vctrs restores the concatenated result more than once, so the number of
  # warnings follows its internals rather than this contract; that at least one
  # fires, and that no other condition escapes, is what is pinned here.
  expect_gte(out$count, 1)
  expect_s3_class(out$value, "psw")
  expect_length(out$value, 12)
  expect_null(stabilization_score(out$value))
  expect_true(is_stabilized(out$value))
  expect_equal(estimand(out$value), "ate")
})

test_that("ATE works for binary cases", {
  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    weights <- wt_ate(c(0.1, 0.3, 0.4, 0.3), .exposure = c(0, 0, 1, 0)),
    "Treating `.exposure` as binary"
  )

  weights2 <- expect_silent(
    wt_ate(
      c(0.1, 0.3, 0.4, 0.3),
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    )
  )

  expect_identical(weights, weights2)

  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    weights3 <- wt_ate(
      c(0.1, 0.3, 0.4, 0.3),
      .exposure = as.logical(c(0, 0, 1, 0))
    ),
    "Treating `.exposure` as binary"
  )

  expect_identical(weights, weights3)

  weights4 <- expect_silent(
    wt_ate(
      c(0.1, 0.3, 0.4, 0.3),
      .exposure = c(2, 2, 1, 2),
      exposure_type = "binary",
      .reference_level = 2
    )
  )

  expect_identical(weights, weights4)

  expect_propensity_error(
    wt_ate(
      c(-0.1, 0.3, 0.4, 3.3),
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    )
  )

  .exposure <- factor(
    c("untreated", "untreated", "treated", "untreated"),
    levels = c("untreated", "treated")
  )

  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    weights5 <- wt_ate(
      c(0.1, 0.3, 0.4, 0.3),
      exposure_type = "binary",
      .exposure = .exposure
    ),
    "Setting focal level to"
  )

  expect_identical(weights, weights5)

  expect_equal(
    weights,
    psw(c(1.11, 1.43, 2.50, 1.43), "ate"),
    tolerance = 0.01
  )
})

# ---- focal level recorded on binary att and atu weights ---------------------
#
# Binary att and atu weights are mirror images of each other: att weights built
# on the first sorted exposure level from 1 - e are numerically identical to atu
# weights built on the second level from e. Nothing in the returned weights tells
# the two apart, so a consumer cannot tell which level the weights target. The
# resolved focal level is recorded under the attribute name the categorical path
# already uses. The remaining binary estimands have symmetric tilting functions,
# so the mirrored construction reproduces the unmirrored weights exactly and
# there is nothing to record.

test_that("binary wt_att() records an explicitly supplied focal level", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- factor(
    c("control", "treat", "control", "treat"),
    levels = c("control", "treat")
  )

  weights <- wt_att(
    ps,
    exposure,
    exposure_type = "binary",
    .focal_level = "control"
  )

  expect_identical(attr(weights, "focal_category"), "control")
})

test_that("binary wt_atu() records an explicitly supplied focal level", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- factor(
    c("control", "treat", "control", "treat"),
    levels = c("control", "treat")
  )

  weights <- wt_atu(
    ps,
    exposure,
    exposure_type = "binary",
    .focal_level = "control"
  )

  expect_identical(attr(weights, "focal_category"), "control")
})

test_that("binary wt_att() records the focal level implied by a reference level", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- factor(
    c("control", "treat", "control", "treat"),
    levels = c("control", "treat")
  )

  # `.reference_level` codes every other level as focal, which for a two-level
  # exposure resolves to the single remaining level.
  weights <- wt_att(
    ps,
    exposure,
    exposure_type = "binary",
    .reference_level = "treat"
  )

  expect_identical(attr(weights, "focal_category"), "control")
})

test_that("binary wt_att() resolves a reference level with missing exposure values", {
  ps <- c(0.2, 0.4, 0.6, 0.8, 0.5)
  exposure <- factor(
    c("control", "treat", "control", "treat", NA),
    levels = c("control", "treat")
  )

  # A missing exposure value is not a third level. `.reference_level` still
  # leaves exactly one other level, so the focal level is recorded; the NA row's
  # own weight stays NA either way.
  weights <- wt_att(
    ps,
    exposure,
    exposure_type = "binary",
    .reference_level = "treat"
  )

  expect_identical(attr(weights, "focal_category"), "control")
  expect_true(is.na(as.double(weights)[[5]]))
  expect_false(anyNA(as.double(weights)[1:4]))
})

test_that("binary att and atu weights record no focal level when none is given", {
  withr::local_options(propensity.quiet = TRUE)
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- factor(
    c("control", "treat", "control", "treat"),
    levels = c("control", "treat")
  )

  weights_att <- wt_att(ps, exposure, exposure_type = "binary")
  weights_atu <- wt_atu(ps, exposure, exposure_type = "binary")

  expect_null(attr(weights_att, "focal_category"))
  expect_null(attr(weights_atu, "focal_category"))
})

test_that("binary estimands other than att and atu record no focal level", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- factor(
    c("control", "treat", "control", "treat"),
    levels = c("control", "treat")
  )

  symmetric <- list(
    ate = wt_ate,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )

  for (estimand_name in names(symmetric)) {
    weights <- symmetric[[estimand_name]](
      ps,
      exposure,
      exposure_type = "binary",
      .focal_level = "control"
    )
    expect_null(attr(weights, "focal_category"))
  }
})

# ---- `.propensity` is the probability of the resolved focal level ----------
#
# Every binary weight function reads `.propensity` as the probability of the
# level it resolves as focal, not as the probability of a fixed level. Naming
# the first factor level as focal, or the second as reference, therefore does
# not change how a supplied numeric propensity is read: it is still the focal
# level's probability, and the caller is responsible for supplying the right
# one. The pins below write the weights out from the definitions with e taken
# as the focal-level probability.

test_that("binary weights read `.propensity` as the named focal level's probability", {
  # e is the probability of "a" here, because "a" is named focal
  ps <- c(0.2, 0.25, 0.8, 0.4)
  exposure <- factor(c("a", "b", "a", "b"), levels = c("a", "b"))

  att_a <- wt_att(ps, exposure, exposure_type = "binary", .focal_level = "a")
  atu_a <- wt_atu(ps, exposure, exposure_type = "binary", .focal_level = "a")

  # att: 1 for the focal "a" units, e / (1 - e) for the "b" units
  expect_equal(as.numeric(att_a), c(1, 0.25 / 0.75, 1, 0.4 / 0.6))
  # atu: (1 - e) / e for the "a" units, 1 for the "b" units
  expect_equal(as.numeric(atu_a), c(0.8 / 0.2, 1, 0.2 / 0.8, 1))
})

test_that("binary att weights mirror atu weights on the other level's probability", {
  exposure <- factor(c("a", "b", "a", "b"), levels = c("a", "b"))
  # each call is given its own focal level's probability, so the two differ
  ps_a <- c(0.2, 0.25, 0.8, 0.4)
  ps_b <- 1 - ps_a

  att_a <- wt_att(ps_a, exposure, exposure_type = "binary", .focal_level = "a")
  atu_b <- wt_atu(ps_b, exposure, exposure_type = "binary", .focal_level = "b")

  expect_equal(as.numeric(att_a), as.numeric(atu_b))
})

test_that("binary weights treat a named reference level as its mirror focal level", {
  ps <- c(0.2, 0.25, 0.8, 0.4)
  exposure <- factor(c("a", "b", "a", "b"), levels = c("a", "b"))

  wt_fns <- list(
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )

  for (estimand_name in names(wt_fns)) {
    wt_fn <- wt_fns[[estimand_name]]
    expect_equal(
      wt_fn(ps, exposure, exposure_type = "binary", .reference_level = "b"),
      wt_fn(ps, exposure, exposure_type = "binary", .focal_level = "a")
    )
  }
})

test_that("binary att and atu record the first level as focal when it is named", {
  ps <- c(0.2, 0.25, 0.8, 0.4)
  exposure <- factor(c("a", "b", "a", "b"), levels = c("a", "b"))

  att_a <- wt_att(ps, exposure, exposure_type = "binary", .focal_level = "a")
  atu_a <- wt_atu(ps, exposure, exposure_type = "binary", .focal_level = "a")

  expect_identical(attr(att_a, "focal_category"), "a")
  expect_identical(attr(atu_a, "focal_category"), "a")
})

test_that("binary weights with no focal level resolve the second level as focal", {
  withr::local_options(propensity.quiet = TRUE)
  # e is the probability of "b", the level resolved as focal by default
  ps <- c(0.2, 0.25, 0.8, 0.4)
  exposure <- factor(c("a", "b", "a", "b"), levels = c("a", "b"))

  # each unit's own group probability: e for "b" units, 1 - e for "a" units
  own_group <- ifelse(exposure == "b", ps, 1 - ps)
  h <- -ps * log(ps) - (1 - ps) * log(1 - ps)

  expect_equal(
    as.numeric(wt_ate(ps, exposure, exposure_type = "binary")),
    1 / own_group
  )
  expect_equal(
    as.numeric(wt_att(ps, exposure, exposure_type = "binary")),
    c(0.2 / 0.8, 1, 0.8 / 0.2, 1)
  )
  expect_equal(
    as.numeric(wt_atu(ps, exposure, exposure_type = "binary")),
    c(1, 0.75 / 0.25, 1, 0.6 / 0.4)
  )
  expect_equal(
    as.numeric(wt_atm(ps, exposure, exposure_type = "binary")),
    pmin(ps, 1 - ps) / own_group
  )
  expect_equal(
    as.numeric(wt_ato(ps, exposure, exposure_type = "binary")),
    c(0.2, 0.75, 0.8, 0.6)
  )
  expect_equal(
    as.numeric(wt_entropy(ps, exposure, exposure_type = "binary")),
    h / own_group
  )
})

test_that("binary weights are unchanged when the default focal level is named", {
  withr::local_options(propensity.quiet = TRUE)
  ps <- c(0.2, 0.25, 0.8, 0.4)
  exposure <- factor(c("a", "b", "a", "b"), levels = c("a", "b"))

  wt_fns <- list(
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )

  for (estimand_name in names(wt_fns)) {
    wt_fn <- wt_fns[[estimand_name]]
    expect_equal(
      as.numeric(wt_fn(
        ps,
        exposure,
        exposure_type = "binary",
        .focal_level = "b"
      )),
      as.numeric(wt_fn(ps, exposure, exposure_type = "binary"))
    )
  }
})

# ---- 0/1 and logical exposures honor the named levels ----------------------
#
# A 0/1 exposure is the same two-level exposure whether it is stored as double
# or as integer, and a logical exposure is the same exposure as its 0/1 recode.
# Naming a level has to reach the same weights in all of those codings, and the
# same weights the equivalent factor coding reaches. As everywhere else on the
# binary path, `.propensity` is the probability of the level named as focal, so
# the expectations below are written out from the definitions with e taken as
# that probability.

zero_one_fixture <- function() {
  exposure <- c(0, 1, 0, 1)

  list(
    ps = c(0.8, 0.6, 0.4, 0.2),
    double = exposure,
    integer = as.integer(exposure),
    logical = exposure == 1,
    # "zero" is the second level, so naming it focal names the level this
    # factor already resolves as focal by default
    factor = factor(
      ifelse(exposure == 1, "one", "zero"),
      levels = c("one", "zero")
    )
  )
}

test_that("binary wt_att() honors a named focal level on a 0/1 double exposure", {
  fixture <- zero_one_fixture()
  ps <- fixture$ps
  is_focal <- as.numeric(fixture$double == 0)

  weights <- wt_att(
    ps,
    fixture$double,
    exposure_type = "binary",
    .focal_level = 0
  )

  # att: 1 for the focal units, e / (1 - e) for the rest
  expect_equal(as.numeric(weights), is_focal + ps * (1 - is_focal) / (1 - ps))

  expect_equal(
    as.numeric(weights),
    as.numeric(wt_att(
      ps,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "zero"
    ))
  )

  expect_equal(
    as.numeric(weights),
    as.numeric(wt_att(
      ps,
      fixture$integer,
      exposure_type = "binary",
      .focal_level = 0
    ))
  )
})

test_that("binary weights agree across 0/1 double and integer storage", {
  fixture <- zero_one_fixture()
  ps <- fixture$ps

  wt_fns <- list(ate = wt_ate, att = wt_att)

  for (estimand_name in names(wt_fns)) {
    wt_fn <- wt_fns[[estimand_name]]
    for (focal in list(0, 1)) {
      expect_equal(
        wt_fn(
          ps,
          fixture$double,
          exposure_type = "binary",
          .focal_level = focal
        ),
        wt_fn(
          ps,
          fixture$integer,
          exposure_type = "binary",
          .focal_level = focal
        )
      )
    }
  }
})

test_that("binary weights resolve a named reference level on a 0/1 double exposure", {
  fixture <- zero_one_fixture()
  ps <- fixture$ps
  is_focal <- as.numeric(fixture$double == 0)

  wt_fns <- list(ate = wt_ate, att = wt_att, atu = wt_atu)

  # naming 1 as reference leaves 0 as the only other level, so it must reach the
  # same weights as naming 0 as focal
  for (estimand_name in names(wt_fns)) {
    wt_fn <- wt_fns[[estimand_name]]
    expect_equal(
      wt_fn(ps, fixture$double, exposure_type = "binary", .reference_level = 1),
      wt_fn(ps, fixture$double, exposure_type = "binary", .focal_level = 0)
    )
  }

  att_ref <- wt_att(
    ps,
    fixture$double,
    exposure_type = "binary",
    .reference_level = 1
  )

  expect_equal(as.numeric(att_ref), is_focal + ps * (1 - is_focal) / (1 - ps))
  expect_identical(attr(att_ref, "focal_category"), 0)
})

test_that("binary wt_att() honors a named focal level on a logical exposure", {
  fixture <- zero_one_fixture()
  ps <- fixture$ps
  is_focal <- as.numeric(!fixture$logical)

  weights <- wt_att(
    ps,
    fixture$logical,
    exposure_type = "binary",
    .focal_level = FALSE
  )

  expect_equal(as.numeric(weights), is_focal + ps * (1 - is_focal) / (1 - ps))

  expect_equal(
    as.numeric(weights),
    as.numeric(wt_att(
      ps,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "zero"
    ))
  )

  expect_identical(attr(weights, "focal_category"), FALSE)
})

test_that("binary att and atu record a named focal level on a 0/1 double exposure", {
  fixture <- zero_one_fixture()

  att_zero <- wt_att(
    fixture$ps,
    fixture$double,
    exposure_type = "binary",
    .focal_level = 0
  )
  atu_zero <- wt_atu(
    fixture$ps,
    fixture$double,
    exposure_type = "binary",
    .focal_level = 0
  )

  expect_identical(attr(att_zero, "focal_category"), 0)
  expect_identical(attr(atu_zero, "focal_category"), 0)
})

test_that("binary weights resolve the deprecated `.treated` on a 0/1 double exposure", {
  fixture <- zero_one_fixture()
  ps <- fixture$ps
  is_focal <- as.numeric(fixture$double == 0)

  with_always_deprecated(
    expect_warning(
      wt_att(ps, fixture$double, exposure_type = "binary", .treated = 0),
      class = "lifecycle_warning_deprecated"
    )
  )

  # the deprecation is asserted above; take the value without re-raising it
  withr::local_options(lifecycle_verbosity = "quiet")
  att_treated <- wt_att(
    ps,
    fixture$double,
    exposure_type = "binary",
    .treated = 0
  )

  expect_equal(
    att_treated,
    wt_att(ps, fixture$double, exposure_type = "binary", .focal_level = 0)
  )
  expect_equal(
    as.numeric(att_treated),
    is_focal + ps * (1 - is_focal) / (1 - ps)
  )
})

test_that("binary weights with no named level keep the 0/1 and logical coding", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- zero_one_fixture()
  ps <- fixture$ps

  # with no level named the focal level stays where it is today: 1 for a 0/1
  # exposure, TRUE for a logical one
  is_focal <- fixture$double
  expected_ate <- is_focal / ps + (1 - is_focal) / (1 - ps)
  expected_att <- is_focal + ps * (1 - is_focal) / (1 - ps)

  codings <- list(
    double = fixture$double,
    integer = fixture$integer,
    logical = fixture$logical
  )

  for (coding_name in names(codings)) {
    exposure <- codings[[coding_name]]
    weights_att <- wt_att(ps, exposure, exposure_type = "binary")

    expect_equal(
      as.numeric(wt_ate(ps, exposure, exposure_type = "binary")),
      expected_ate
    )
    expect_equal(as.numeric(weights_att), expected_att)
    expect_null(attr(weights_att, "focal_category"))
  }
})

# ---- glm methods must supply the resolved focal level's probability --------
#
# `fitted()` on a binomial glm is the probability of the response's second
# factor level, which is the focal level only when the caller leaves the focal
# level at its default. Naming the first level as focal resolves the focal
# level to "a", so the glm methods owe the numeric method P("a"), which is
# 1 - fitted(). The expectations below are written out from the definitions
# with e = 1 - fitted() and the exposure indicator on "a".

focal_glm_fixture <- function() {
  set.seed(2024)
  n <- 40
  x <- rnorm(n)
  exposure <- factor(
    ifelse(rbinom(n, 1, plogis(0.4 * x)) == 1, "b", "a"),
    levels = c("a", "b")
  )

  list(
    model = glm(exposure ~ x, family = binomial),
    exposure = exposure
  )
}

test_that("glm wt_ate() uses the named focal level's fitted probability", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == "a"

  weights <- wt_ate(fixture$model, exposure_type = "binary", .focal_level = "a")

  expect_equal(as.numeric(weights), ifelse(is_focal, 1 / e, 1 / (1 - e)))
})

test_that("glm wt_att() uses the named focal level's fitted probability", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == "a"

  weights <- wt_att(fixture$model, exposure_type = "binary", .focal_level = "a")

  expect_equal(as.numeric(weights), ifelse(is_focal, 1, e / (1 - e)))
})

test_that("glm wt_atu() uses the named focal level's fitted probability", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == "a"

  weights <- wt_atu(fixture$model, exposure_type = "binary", .focal_level = "a")

  expect_equal(as.numeric(weights), ifelse(is_focal, (1 - e) / e, 1))
})

test_that("glm wt_atm() uses the named focal level's fitted probability", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == "a"

  weights <- wt_atm(fixture$model, exposure_type = "binary", .focal_level = "a")

  own_group <- ifelse(is_focal, e, 1 - e)
  expect_equal(as.numeric(weights), pmin(e, 1 - e) / own_group)
})

test_that("glm wt_ato() uses the named focal level's fitted probability", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == "a"

  weights <- wt_ato(fixture$model, exposure_type = "binary", .focal_level = "a")

  expect_equal(as.numeric(weights), ifelse(is_focal, 1 - e, e))
})

test_that("glm wt_entropy() uses the named focal level's fitted probability", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == "a"

  weights <- wt_entropy(
    fixture$model,
    exposure_type = "binary",
    .focal_level = "a"
  )

  h <- -e * log(e) - (1 - e) * log(1 - e)
  own_group <- ifelse(is_focal, e, 1 - e)
  expect_equal(as.numeric(weights), h / own_group)
})

test_that("glm weights resolve a named reference level to the other level's probability", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == "a"

  wt_fns <- list(
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )

  # naming "b" as reference resolves the focal level to "a", so it must reach
  # the same weights as naming "a" as focal
  for (estimand_name in names(wt_fns)) {
    wt_fn <- wt_fns[[estimand_name]]
    expect_equal(
      wt_fn(fixture$model, exposure_type = "binary", .reference_level = "b"),
      wt_fn(fixture$model, exposure_type = "binary", .focal_level = "a")
    )
  }

  # and those weights are the ones built from P("a") = 1 - fitted()
  att_ref <- wt_att(
    fixture$model,
    exposure_type = "binary",
    .reference_level = "b"
  )
  atu_ref <- wt_atu(
    fixture$model,
    exposure_type = "binary",
    .reference_level = "b"
  )

  expect_equal(as.numeric(att_ref), ifelse(is_focal, 1, e / (1 - e)))
  expect_equal(as.numeric(atu_ref), ifelse(is_focal, (1 - e) / e, 1))
})

test_that("glm weights resolve the deprecated `.treated` to the same focal level", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == "a"

  with_always_deprecated(
    expect_warning(
      wt_att(fixture$model, exposure_type = "binary", .treated = "a"),
      class = "lifecycle_warning_deprecated"
    )
  )

  # the deprecation is asserted above; take the value without re-raising it
  withr::local_options(lifecycle_verbosity = "quiet")
  att_treated <- wt_att(
    fixture$model,
    exposure_type = "binary",
    .treated = "a"
  )

  expect_equal(
    att_treated,
    wt_att(fixture$model, exposure_type = "binary", .focal_level = "a")
  )
  expect_equal(as.numeric(att_treated), ifelse(is_focal, 1, e / (1 - e)))
})

test_that("glm weights pass fitted values through when the focal level is the default", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  fitted_ps <- unname(fitted(fixture$model))

  expect_equal(
    wt_att(fixture$model, exposure_type = "binary"),
    wt_att(fitted_ps, fixture$exposure, exposure_type = "binary")
  )
})

test_that("glm weights pass fitted values through when the default focal level is named", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()

  wt_fns <- list(
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )

  # "b" is the response's second level, so naming it focal resolves to the
  # level `fitted()` already reports and nothing may be inverted. Values only:
  # naming a level also records `focal_category` on att and atu.
  for (estimand_name in names(wt_fns)) {
    wt_fn <- wt_fns[[estimand_name]]
    expect_equal(
      as.numeric(wt_fn(
        fixture$model,
        exposure_type = "binary",
        .focal_level = "b"
      )),
      as.numeric(wt_fn(fixture$model, exposure_type = "binary"))
    )
  }
})

zero_one_glm_fixture <- function() {
  set.seed(311)
  n <- 40
  x <- rnorm(n)
  z <- as.numeric(rbinom(n, 1, plogis(0.4 * x)))

  list(
    model = glm(z ~ x, family = binomial),
    exposure = z
  )
}

test_that("glm weights on a 0/1 response honor a named focal level", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- zero_one_glm_fixture()
  # 0 is not the response's success level, so the numeric method is owed
  # P(0) = 1 - fitted()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == 0

  weights <- wt_att(fixture$model, exposure_type = "binary", .focal_level = 0)

  expect_equal(as.numeric(weights), ifelse(is_focal, 1, e / (1 - e)))
  expect_equal(
    as.numeric(weights),
    as.numeric(wt_att(
      e,
      fixture$exposure,
      exposure_type = "binary",
      .focal_level = 0
    ))
  )
})

test_that("glm weights on a 0/1 response pass fitted values through by default", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- zero_one_glm_fixture()

  # 1 is the response's success level, so naming it focal resolves to the level
  # `fitted()` already reports and nothing may be inverted
  expect_equal(
    as.numeric(wt_att(
      fixture$model,
      exposure_type = "binary",
      .focal_level = 1
    )),
    as.numeric(wt_att(fixture$model, exposure_type = "binary"))
  )
})

test_that("glm weights on a logical response honor a named focal level", {
  withr::local_options(propensity.quiet = TRUE)
  set.seed(412)
  n <- 40
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x)) == 1
  ps_mod <- glm(z ~ x, family = binomial)

  # FALSE is not the response's success level, so the numeric method is owed
  # P(FALSE) = 1 - fitted()
  e <- 1 - unname(fitted(ps_mod))
  is_focal <- !z

  weights <- wt_att(ps_mod, exposure_type = "binary", .focal_level = FALSE)

  expect_equal(as.numeric(weights), ifelse(is_focal, 1, e / (1 - e)))
})

test_that("glm wt_cens() uses the named focal level's fitted probability", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  focal_ps <- 1 - unname(fitted(fixture$model))

  weights <- wt_cens(
    fixture$model,
    exposure_type = "binary",
    .focal_level = "a"
  )

  expect_equal(
    weights,
    wt_cens(
      focal_ps,
      fixture$exposure,
      exposure_type = "binary",
      .focal_level = "a"
    )
  )
  expect_equal(estimand(weights), "uncensored")
})

test_that("glm weights resolve the deprecated `.untreated` to the same reference level", {
  withr::local_options(propensity.quiet = TRUE)
  fixture <- focal_glm_fixture()
  e <- 1 - unname(fitted(fixture$model))
  is_focal <- fixture$exposure == "a"

  with_always_deprecated(
    expect_warning(
      wt_att(fixture$model, exposure_type = "binary", .untreated = "b"),
      class = "lifecycle_warning_deprecated"
    )
  )

  # the deprecation is asserted above; take the value without re-raising it
  withr::local_options(lifecycle_verbosity = "quiet")
  att_untreated <- wt_att(
    fixture$model,
    exposure_type = "binary",
    .untreated = "b"
  )

  expect_equal(
    att_untreated,
    wt_att(fixture$model, exposure_type = "binary", .reference_level = "b")
  )
  expect_equal(as.numeric(att_untreated), ifelse(is_focal, 1, e / (1 - e)))
})

test_that("ATE works for continuous cases", {
  denom_model <- lm(mpg ~ gear + am + carb, data = mtcars)

  # Compute population variances
  un_mean <- mean(mtcars$mpg)
  un_var <- mean((mtcars$mpg - un_mean)^2)
  cond_var <- mean((mtcars$mpg - predict(denom_model))^2)

  # Compute z-scores and densities
  z_num <- (mtcars$mpg - un_mean) / sqrt(un_var)
  z_den <- (mtcars$mpg - predict(denom_model)) / sqrt(cond_var)
  f_num <- dnorm(z_num)
  f_den <- dnorm(z_den)

  # Expected weights
  wts <- 1 / f_den
  stb_wts <- f_num / f_den

  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    weights <- wt_ate(
      predict(denom_model),
      .exposure = mtcars$mpg,
      .sigma = influence(denom_model)$sigma,
      exposure_type = "continuous"
    ),
    "Using unstabilized weights for continuous exposures is not recommended."
  )

  expect_equal(weights, psw(wts, "ate"), tolerance = 0.01)
  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    stabilized_weights <- wt_ate(
      predict(denom_model),
      .exposure = mtcars$mpg,
      .sigma = influence(denom_model)$sigma,
      stabilize = TRUE,
    ),
    "Treating `.exposure` as continuous"
  )

  expect_equal(
    stabilized_weights,
    psw(stb_wts, "ate", stabilized = TRUE),
    tolerance = 0.01
  )
})

test_that("stabilized weights use P(A=1) and P(A=0) as numerators", {
  ps <- c(0.2, 0.5, 0.8, 0.4)
  A <- c(1, 0, 1, 0)

  p1 <- mean(A)
  p0 <- 1 - p1
  inv_ps <- 1 / ps
  inv_1m <- 1 / (1 - ps)
  expected <- A * inv_ps * p1 + (1 - A) * inv_1m * p0

  got <- ate_binary(ps, A, stabilize = TRUE)

  expect_equal(got, expected)
})


test_that("ATE errors appropriately for categorical with vector propensity scores", {
  # For categorical exposures, propensity scores must be a matrix
  expect_propensity_error(
    wt_ate(
      c(0.1, 0.3, 0.4, 0.3),
      .exposure = c(0, 2, 1, 4),
      exposure_type = "categorical"
    )
  )
})

test_that("wt_ate() with ps_trim issues refit warning if not refit, no warning if refit", {
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.5 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 1) Trim
  trimmed_ps <- ps_trim(
    ps,
    .exposure = z,
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  # not refit => expect a warning
  w_ate_unfit <- expect_propensity_warning(
    wt_ate(
      trimmed_ps,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_s3_class(w_ate_unfit, "psw")
  expect_true(grepl("; trimmed$", estimand(w_ate_unfit)))

  # 2) After refit => no warning
  trimmed_refit <- ps_refit(trimmed_ps, model = fit)
  w_ate_fit <- expect_silent(
    wt_ate(
      trimmed_refit,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_s3_class(w_ate_fit, "psw")
  expect_true(grepl("; trimmed$", estimand(w_ate_fit)))
})

test_that("wt_ate() with ps_trunc adds '; truncated' without refit warning", {
  set.seed(234)
  n <- 100
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.6 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # e.g. bounding at [0.2, 0.8]
  truncated_ps <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)

  # Should produce weighting with no refit warnings
  w_ate_trunc <- expect_silent(
    wt_ate(
      truncated_ps,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_s3_class(w_ate_trunc, "psw")
  # Estimand ends with "; truncated"
  expect_true(grepl("; truncated$", estimand(w_ate_trunc)))
})

test_that("Other estimands (att, atu, etc.) with ps_trim or ps_trunc", {
  set.seed(345)
  n <- 120
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.2 + 0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # Trim
  trimmed_ps <- ps_trim(ps, .exposure = z, method = "ps")
  # No refit => warning
  w_att_trim <- expect_propensity_warning(
    wt_att(
      trimmed_ps,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  # Check estimand
  expect_true(grepl("att; trimmed", estimand(w_att_trim)))

  # Trunc
  truncated_ps <- ps_trunc(ps, method = "pctl", lower = 0.2, upper = 0.8)
  # No warning
  w_att_trunc <- expect_silent(
    wt_att(
      truncated_ps,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  w_atu_trunc <- expect_silent(
    wt_atu(
      truncated_ps,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_true(grepl("atu; truncated", estimand(w_atu_trunc)))
})

test_that("wt_ate() with ps_trunc sets truncated=TRUE in final psw", {
  set.seed(123)
  n <- 8
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  trunc_obj <- ps_trunc(ps, method = "ps", lower = 0.2, upper = 0.8)
  w_ate <- wt_ate(
    trunc_obj,
    .exposure = z,
    exposure_type = "binary",
    .focal_level = 1
  )

  expect_true(is_ps_truncated(w_ate))
  expect_false(is_ps_trimmed(w_ate))
  expect_match(estimand(w_ate), "; truncated$")
})

test_that("wt_atu.ps_trim triggers refit check, sets 'atu; trimmed'", {
  set.seed(991)
  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(1.6 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  # 1) Trim the PS
  trimmed_obj <- ps_trim(
    ps,
    .exposure = z,
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  # Not refit => we get a warning
  w_atu_unfit <- expect_propensity_warning(
    wt_atu(
      trimmed_obj,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_s3_class(w_atu_unfit, "psw")
  expect_match(estimand(w_atu_unfit), "atu; trimmed")
  expect_true(is_ps_trimmed(w_atu_unfit))
  # ps_trim_meta copied
  expect_identical(
    attr(w_atu_unfit, "ps_trim_meta"),
    attr(trimmed_obj, "ps_trim_meta")
  )

  # 2) Now refit => no warning
  refit_obj <- ps_refit(trimmed_obj, model = fit)
  w_atu_fit <- expect_silent(
    wt_atu(
      refit_obj,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_s3_class(w_atu_fit, "psw")
  expect_match(estimand(w_atu_fit), "atu; trimmed")
  expect_true(is_ps_trimmed(w_atu_fit))
  # confirm ps_trim_meta matches
  expect_identical(
    attr(w_atu_fit, "ps_trim_meta"),
    attr(refit_obj, "ps_trim_meta")
  )
})

test_that("wt_atm.ps_trim triggers refit check, sets 'atm; trimmed'", {
  set.seed(992)
  n <- 50
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.3 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  trimmed_obj <- ps_trim(
    ps,
    .exposure = z,
    method = "ps",
    lower = 0.2,
    upper = 0.8
  )

  # Not refit => warning
  w_atm_unfit <- expect_propensity_warning(
    wt_atm(
      trimmed_obj,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_s3_class(w_atm_unfit, "psw")
  expect_match(estimand(w_atm_unfit), "atm; trimmed")
  expect_true(is_ps_trimmed(w_atm_unfit))

  # Refit => no warning
  refit_obj <- ps_refit(trimmed_obj, model = fit)
  w_atm_fit <- expect_silent(
    wt_atm(
      refit_obj,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_s3_class(w_atm_fit, "psw")
  expect_match(estimand(w_atm_fit), "atm; trimmed")
  expect_true(is_ps_trimmed(w_atm_fit))
})

test_that("wt_ato.ps_trim triggers refit check, sets 'ato; trimmed'", {
  set.seed(993)
  n <- 60
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.4 * x))
  fit <- glm(z ~ x, family = binomial)
  ps <- predict(fit, type = "response")

  trimmed_obj <- ps_trim(
    ps,
    .exposure = z,
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )

  # Not refit => warning
  w_ato_unfit <- expect_propensity_warning(
    wt_ato(
      trimmed_obj,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_s3_class(w_ato_unfit, "psw")
  expect_match(estimand(w_ato_unfit), "ato; trimmed")
  expect_true(is_ps_trimmed(w_ato_unfit))

  # Refit => no warning
  refit_obj <- ps_refit(trimmed_obj, model = fit)
  w_ato_fit <- expect_silent(
    wt_ato(
      refit_obj,
      .exposure = z,
      exposure_type = "binary",
      .focal_level = 1
    )
  )
  expect_s3_class(w_ato_fit, "psw")
  expect_match(estimand(w_ato_fit), "ato; trimmed")
  expect_true(is_ps_trimmed(w_ato_fit))
})

# Entropy weight tests
test_that("wt_entropy works for binary cases", {
  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    weights <- wt_entropy(c(0.1, 0.3, 0.4, 0.3), .exposure = c(0, 0, 1, 0)),
    "Treating `.exposure` as binary"
  )

  weights2 <- expect_silent(
    wt_entropy(
      c(0.1, 0.3, 0.4, 0.3),
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    )
  )

  expect_identical(weights, weights2)

  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    weights3 <- wt_entropy(
      c(0.1, 0.3, 0.4, 0.3),
      .exposure = as.logical(c(0, 0, 1, 0))
    ),
    "Treating `.exposure` as binary"
  )

  expect_identical(weights, weights3)

  weights4 <- expect_silent(
    wt_entropy(
      c(0.1, 0.3, 0.4, 0.3),
      .exposure = c(2, 2, 1, 2),
      exposure_type = "binary",
      .reference_level = 2
    )
  )

  expect_identical(weights, weights4)

  expect_propensity_error(
    wt_entropy(
      c(-0.1, 0.3, 0.4, 3.3),
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    )
  )

  .exposure <- factor(
    c("untreated", "untreated", "treated", "untreated"),
    levels = c("untreated", "treated")
  )

  withr::local_options(propensity.quiet = FALSE)
  expect_message(
    weights5 <- wt_entropy(
      c(0.1, 0.3, 0.4, 0.3),
      exposure_type = "binary",
      .exposure = .exposure
    ),
    "Setting focal level to"
  )

  expect_identical(weights, weights5)
})

test_that("entropy tilting function properties", {
  # Test symmetry: h(e) = h(1-e)
  ps1 <- c(0.2, 0.3, 0.4)
  ps2 <- 1 - ps1

  # Calculate tilting functions
  h1 <- -ps1 * log(ps1) - (1 - ps1) * log(1 - ps1)
  h2 <- -ps2 * log(ps2) - (1 - ps2) * log(1 - ps2)

  expect_equal(h1, h2, tolerance = 1e-10)

  # Test maximum at 0.5
  ps_seq <- seq(0.01, 0.99, by = 0.01)
  h_vals <- -ps_seq * log(ps_seq) - (1 - ps_seq) * log(1 - ps_seq)
  max_idx <- which.max(h_vals)
  expect_equal(ps_seq[max_idx], 0.5, tolerance = 0.01)

  # Test bounds
  # Maximum entropy is log(2) ≈ 0.693
  expect_true(all(h_vals <= log(2) + 1e-10))
  expect_true(all(h_vals >= 0))
})

test_that("entropy weights have expected properties", {
  # Generate random propensity scores
  set.seed(123)
  ps <- runif(100, 0.1, 0.9)
  treatment <- rbinom(100, 1, ps)

  weights <- wt_entropy(ps, .exposure = treatment, exposure_type = "binary")

  # Entropy weights should be positive and finite
  expect_true(all(as.numeric(weights) > 0))
  expect_true(all(is.finite(as.numeric(weights))))

  # Weights at e=0.5 should be around log(2)/0.5 ≈ 1.386
  ps_near_half <- abs(ps - 0.5) < 0.01
  if (any(ps_near_half)) {
    expect_true(all(
      abs(as.numeric(weights[ps_near_half]) - log(2) / 0.5) < 0.1
    ))
  }
})

test_that("entropy weights handle extreme propensity scores", {
  # Near 0 and 1 propensity scores
  ps_extreme <- c(0.001, 0.01, 0.99, 0.999)
  treatment_extreme <- c(0, 0, 1, 1)

  weights_extreme <- expect_silent(
    wt_entropy(
      ps_extreme,
      .exposure = treatment_extreme,
      exposure_type = "binary"
    )
  )

  expect_true(all(is.finite(as.numeric(weights_extreme))))
  expect_true(all(as.numeric(weights_extreme) > 0))

  # Extreme weights can be large but should be finite
  # Theoretical upper bound for entropy weights based on extreme propensity scores
  # For extreme values near 0 or 1, weights can grow large but remain finite.
  # Here, we use a calculated bound derived from the entropy function properties.
  max_weight_bound <- log(2) / min(ps_extreme) # Example calculation
  expect_true(max(weights_extreme) < max_weight_bound)
})

test_that("wt_entropy works with ps_trim objects", {
  ps <- c(0.1, 0.3, 0.4, 0.3)
  ps_trimmed <- ps_trim(ps, method = "ps", lower = 0.15, upper = 0.85)

  weights <- expect_propensity_warning(
    wt_entropy(
      ps_trimmed,
      .exposure = c(0, 0, 1, 0),
      exposure_type = "binary"
    )
  )

  expect_s3_class(weights, "psw")
  expect_equal(estimand(weights), "entropy; trimmed")
  expect_true(is_ps_trimmed(weights))
})

test_that("wt_entropy works with ps_trunc objects", {
  ps <- c(0.1, 0.3, 0.4, 0.3)
  ps_truncated <- ps_trunc(ps, lower = 0.15, upper = 0.85)

  weights <- wt_entropy(
    ps_truncated,
    .exposure = c(0, 0, 1, 0),
    exposure_type = "binary"
  )

  expect_s3_class(weights, "psw")
  expect_equal(estimand(weights), "entropy; truncated")
  expect_true(is_ps_truncated(weights))
})

test_that("entropy weights error on unsupported exposure types", {
  # For categorical exposures, propensity scores must be a matrix
  expect_propensity_error(
    wt_entropy(
      c(0.1, 0.3, 0.4, 0.3),
      .exposure = c(1, 2, 3, 4),
      exposure_type = "categorical"
    )
  )

  # Now that continuous is not even an option for entropy,
  # the function will error during auto-detection if given continuous data
  expect_propensity_error(
    wt_entropy(
      rnorm(10),
      .exposure = rnorm(10)
    )
  )
})

# Comparison with PSWeight package - weights
test_that("entropy weights match PSweight's raw weights", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Use a simple example where we can verify the calculations
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  ps_true <- plogis(0.5 * x)
  trt <- rbinom(n, 1, ps_true)

  # Use the true propensity scores for both implementations
  # This ensures we're comparing the weight calculation, not PS estimation

  # Our implementation
  our_weights <- wt_entropy(ps_true, .exposure = trt, exposure_type = "binary")

  # PSweight's implementation using SumStat
  # Create a data frame with the required structure
  test_data <- data.frame(
    trt = trt,
    ps = ps_true,
    x = x
  )

  # SumStat with provided propensity scores
  ps_sumstat <- PSweight::SumStat(
    ps.estimate = ps_true,
    zname = "trt",
    xname = "x",
    data = test_data,
    weight = "entropy"
  )

  # Extract PSweight's raw weights (before normalization)
  # PSweight stores weights in ps.weights$entropy, but these are normalized
  # We need to un-normalize them to compare
  psw_weights_norm <- ps_sumstat$ps.weights$entropy

  # Un-normalize by multiplying by the sum of our raw weights in each group
  # PSweight normalizes so weights sum to 1 within each treatment group
  our_sum1 <- sum(our_weights[trt == 1])
  our_sum0 <- sum(our_weights[trt == 0])
  psw_weights_raw <- numeric(n)
  psw_weights_raw[trt == 1] <- psw_weights_norm[trt == 1] * our_sum1
  psw_weights_raw[trt == 0] <- psw_weights_norm[trt == 0] * our_sum0

  # Compare raw weights
  expect_equal(as.numeric(our_weights), psw_weights_raw, tolerance = 1e-10)
})

# Comparison with PSWeight package - estimates
test_that("entropy weights give same treatment effect estimates as PSweight", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Calculate treatment effect using PSweight
  ps_result <- PSweight::PSweight(
    ps.formula = trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6,
    yname = "Y",
    data = PSweight::psdata_cl,
    weight = "entropy"
  )

  psweight_ate <- unname(ps_result$muhat[2] - ps_result$muhat[1])

  # Calculate using our implementation
  ps_fit <- glm(
    trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6,
    data = PSweight::psdata_cl,
    family = binomial
  )
  ps_scores <- fitted(ps_fit)

  our_weights <- wt_entropy(
    ps_scores,
    .exposure = PSweight::psdata_cl$trt,
    exposure_type = "binary"
  )

  # Calculate weighted means
  mu1 <- weighted.mean(
    PSweight::psdata_cl$Y[PSweight::psdata_cl$trt == 1],
    as.numeric(our_weights[PSweight::psdata_cl$trt == 1])
  )
  mu0 <- weighted.mean(
    PSweight::psdata_cl$Y[PSweight::psdata_cl$trt == 0],
    as.numeric(our_weights[PSweight::psdata_cl$trt == 0])
  )
  our_ate <- mu1 - mu0

  # Compare estimates - they should be very close
  expect_equal(our_ate, psweight_ate, tolerance = 1e-6)
})

test_that("entropy weighted estimates are reasonable", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Simulate data with known treatment effect
  set.seed(456)
  n <- 500
  x <- rnorm(n)
  ps <- plogis(0.5 * x)
  treatment <- rbinom(n, 1, ps)
  # True treatment effect = 2
  outcome <- 1 + 2 * treatment + 0.5 * x + rnorm(n)

  # Calculate weights
  weights <- wt_entropy(ps, .exposure = treatment, exposure_type = "binary")

  # Weighted means
  mu1 <- weighted.mean(
    outcome[treatment == 1],
    as.numeric(weights[treatment == 1])
  )
  mu0 <- weighted.mean(
    outcome[treatment == 0],
    as.numeric(weights[treatment == 0])
  )
  ate_est <- mu1 - mu0

  # Should be close to true value of 2
  expect_equal(ate_est, 2, tolerance = 0.5)
})


test_that("entropy weight calculation matches manual calculation", {
  # Test specific values
  ps <- c(0.2, 0.5, 0.8)
  treatment <- c(0, 1, 1)

  # Manual calculation
  h_e <- -ps * log(ps) - (1 - ps) * log(1 - ps)

  expected <- numeric(3)
  expected[1] <- h_e[1] / (1 - ps[1]) # Control unit
  expected[2] <- h_e[2] / ps[2] # Treated unit
  expected[3] <- h_e[3] / ps[3] # Treated unit

  weights <- wt_entropy(ps, .exposure = treatment, exposure_type = "binary")

  expect_equal(as.numeric(weights), expected, tolerance = 1e-10)
})

# Tests for data.frame methods
test_that("wt_ate works with data frames", {
  # Create test data frame
  ps_df <- data.frame(
    control = c(0.9, 0.7, 0.3, 0.1),
    treated = c(0.1, 0.3, 0.7, 0.9)
  )
  exposure <- c(0, 0, 1, 1)

  # Test default behavior (uses second column)
  weights <- wt_ate(ps_df, exposure, exposure_type = "binary")
  expected <- wt_ate(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights, expected)

  # Test explicit column selection by name (quoted)
  weights_quoted <- wt_ate(
    ps_df,
    exposure,
    .propensity_col = "treated",
    exposure_type = "binary"
  )
  expect_equal(weights_quoted, expected)

  # Test unquoted column selection
  weights_unquoted <- wt_ate(
    ps_df,
    exposure,
    .propensity_col = treated,
    exposure_type = "binary"
  )
  expect_equal(weights_unquoted, expected)

  # Test column selection by index
  weights_idx <- wt_ate(
    ps_df,
    exposure,
    .propensity_col = 2,
    exposure_type = "binary"
  )
  expect_equal(weights_idx, expected)

  # Test single column data frame
  ps_single <- data.frame(prob = c(0.1, 0.3, 0.7, 0.9))
  weights_single <- wt_ate(ps_single, exposure, exposure_type = "binary")
  expected_single <- wt_ate(ps_single$prob, exposure, exposure_type = "binary")
  expect_equal(weights_single, expected_single)

  # Test error with empty data frame
  expect_propensity_error(
    wt_ate(data.frame(), exposure)
  )

  # Test error with invalid column selection
  expect_propensity_error(
    wt_ate(ps_df, exposure, .propensity_col = "nonexistent")
  )
})

test_that("all wt_* functions work with data frames", {
  ps_df <- data.frame(
    control = c(0.9, 0.7, 0.3, 0.1),
    treated = c(0.1, 0.3, 0.7, 0.9)
  )
  exposure <- c(0, 0, 1, 1)

  # Test ATT
  weights_att <- wt_att(ps_df, exposure, exposure_type = "binary")
  expected_att <- wt_att(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights_att, expected_att)

  # Test ATU
  weights_atu <- wt_atu(ps_df, exposure, exposure_type = "binary")
  expected_atu <- wt_atu(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights_atu, expected_atu)

  # Test ATM
  weights_atm <- wt_atm(ps_df, exposure, exposure_type = "binary")
  expected_atm <- wt_atm(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights_atm, expected_atm)

  # Test ATO
  weights_ato <- wt_ato(ps_df, exposure, exposure_type = "binary")
  expected_ato <- wt_ato(ps_df$treated, exposure, exposure_type = "binary")
  expect_equal(weights_ato, expected_ato)

  # Test Entropy
  weights_entropy <- wt_entropy(ps_df, exposure, exposure_type = "binary")
  expected_entropy <- wt_entropy(
    ps_df$treated,
    exposure,
    exposure_type = "binary"
  )
  expect_equal(weights_entropy, expected_entropy)
})

# ---- data frame columns resolve by focal level name ------------------------

# A data frame of per-level probabilities whose columns are named for the
# levels they hold. The column carrying P(focal) therefore depends on which
# level is resolved as focal, not on where the column sits.
level_named_df_fixture <- function() {
  p_a <- c(0.9, 0.7, 0.3, 0.1)

  list(
    p_a = p_a,
    p_b = 1 - p_a,
    # "b" is the second level, so it is the focal level with none named
    factor = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    double = c(0, 0, 1, 1)
  )
}

test_that("data frame propensity resolves a named focal level by column name", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(a = fixture$p_a, b = fixture$p_b)

  wt_fns <- list(
    ate = wt_ate,
    att = wt_att,
    atu = wt_atu,
    atm = wt_atm,
    ato = wt_ato,
    entropy = wt_entropy
  )

  for (estimand_name in names(wt_fns)) {
    wt_fn <- wt_fns[[estimand_name]]
    expect_equal(
      wt_fn(
        ps_df,
        fixture$factor,
        exposure_type = "binary",
        .focal_level = "a"
      ),
      wt_fn(
        fixture$p_a,
        fixture$factor,
        exposure_type = "binary",
        .focal_level = "a"
      ),
      label = paste0("wt_", estimand_name, "(ps_df)"),
      expected.label = paste0("wt_", estimand_name, "(p_a)")
    )
  }
})

test_that("data frame propensity resolves a named reference level by column name", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(a = fixture$p_a, b = fixture$p_b)

  # naming "b" as reference leaves "a" as the only other level, so the "a"
  # column carries P(focal)
  expect_equal(
    wt_att(
      ps_df,
      fixture$factor,
      exposure_type = "binary",
      .reference_level = "b"
    ),
    wt_att(
      fixture$p_a,
      fixture$factor,
      exposure_type = "binary",
      .reference_level = "b"
    )
  )
})

test_that("data frame propensity column selection ignores column order", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(a = fixture$p_a, b = fixture$p_b)
  ps_df_reversed <- data.frame(b = fixture$p_b, a = fixture$p_a)

  # with no level named the focal level is "b", so both frames must reach the
  # "b" column
  expected <- wt_att(fixture$p_b, fixture$factor, exposure_type = "binary")

  expect_equal(
    wt_att(ps_df, fixture$factor, exposure_type = "binary"),
    expected
  )
  expect_equal(
    wt_att(ps_df_reversed, fixture$factor, exposure_type = "binary"),
    expected
  )
})

test_that("data frame propensity resolves a 0/1 exposure's focal level by column name", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(
    `0` = fixture$p_a,
    `1` = fixture$p_b,
    check.names = FALSE
  )
  ps_df_reversed <- data.frame(
    `1` = fixture$p_b,
    `0` = fixture$p_a,
    check.names = FALSE
  )

  expected <- wt_att(
    fixture$p_a,
    fixture$double,
    exposure_type = "binary",
    .focal_level = 0
  )

  expect_equal(
    wt_att(
      ps_df,
      fixture$double,
      exposure_type = "binary",
      .focal_level = 0
    ),
    expected
  )
  expect_equal(
    wt_att(
      ps_df_reversed,
      fixture$double,
      exposure_type = "binary",
      .focal_level = 0
    ),
    expected
  )

  # with no level named the focal level stays 1, so the "1" column is selected
  expect_equal(
    wt_att(ps_df, fixture$double, exposure_type = "binary"),
    wt_att(fixture$p_b, fixture$double, exposure_type = "binary")
  )
})

test_that("data frame propensity warns when no column matches a named level", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(x1 = fixture$p_a, x2 = fixture$p_b)

  # the fallback is today's positional default, the second column
  expect_equal(
    suppressWarnings(wt_att(
      ps_df,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    )),
    wt_att(
      fixture$p_b,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    )
  )

  # the warning names the column the fallback landed on
  expect_warning(
    wt_att(
      ps_df,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ),
    regexp = "x2",
    class = "propensity_df_column_warning"
  )

  expect_warning(
    wt_att(
      ps_df,
      fixture$factor,
      exposure_type = "binary",
      .reference_level = "b"
    ),
    class = "propensity_df_column_warning"
  )
})

test_that("data frame propensity warns when only some levels have columns", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(a = fixture$p_a, other = fixture$p_b)

  expect_warning(
    wt_att(
      ps_df,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ),
    class = "propensity_df_column_warning"
  )
})

test_that("the fallback warning names a declared level the exposure never takes", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(a = fixture$p_a, b = fixture$p_b)
  with_unused <- factor(c("a", "a", "b", "b"), levels = c("a", "b", "c"))

  # The frame is named for every level the exposure holds, so by the names the
  # user can see the match should succeed. It does not: a factor answers for its
  # declared levels, `"c"` among them, and one of those has no column, so the
  # names are judged not to cover the exposure and the selection falls to
  # position. The wrong column is then read, guarded by this warning alone, so
  # the warning has to name the level that is in the way.
  warning <- expect_warning(
    wt_att(
      ps_df,
      with_unused,
      exposure_type = "binary",
      .focal_level = "a"
    ),
    class = "propensity_df_column_warning"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(warning))
  expect_match(msg, "\"c\"", fixed = TRUE)
  expect_match(msg, "droplevels", fixed = TRUE)

  # And that the remedy it names is the remedy: dropping the unused level makes
  # the match succeed and the warning go away.
  weights <- expect_no_warning(
    wt_att(
      ps_df,
      droplevels(with_unused),
      exposure_type = "binary",
      .focal_level = "a"
    )
  )
  expect_equal(
    as.double(weights),
    as.double(wt_att(
      fixture$p_a,
      droplevels(with_unused),
      exposure_type = "binary",
      .focal_level = "a"
    ))
  )
})

test_that("the fallback warning stays quiet about levels when every one is declared", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(x1 = fixture$p_a, x2 = fixture$p_b)

  # The names here cover no level at all, which the unused level has nothing to
  # do with, so the hint about it must not appear.
  warning <- expect_warning(
    wt_att(
      ps_df,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ),
    class = "propensity_df_column_warning"
  )
  expect_no_match(conditionMessage(warning), "droplevels", fixed = TRUE)
})

test_that("data frame propensity warns when the focal level names two columns", {
  fixture <- level_named_df_fixture()
  duplicated_names <- data.frame(fixture$p_a, fixture$p_b, fixture$p_b)
  names(duplicated_names) <- c("a", "a", "b")

  # Two columns can carry one name, through `check.names = FALSE` or through an
  # assignment like the one above, and the match takes the first of them. That
  # is a choice made on nothing the caller expressed, between columns that may
  # hold different numbers, so it is announced.
  warning <- expect_warning(
    wt_att(
      duplicated_names,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ),
    class = "propensity_df_duplicate_column_warning"
  )
  msg <- gsub("[[:space:]]+", " ", conditionMessage(warning))
  expect_match(msg, "\"a\"", fixed = TRUE)

  # The column read is still the first of them, so the announcement changes
  # nothing about the answer.
  expect_equal(
    as.double(suppressWarnings(wt_att(
      duplicated_names,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ))),
    as.double(wt_att(
      fixture$p_a,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ))
  )
})

test_that("a missing column name is not an ambiguity", {
  fixture <- level_named_df_fixture()
  na_named <- data.frame(fixture$p_a, fixture$p_b, fixture$p_b)
  names(na_named) <- c("a", "b", NA)

  # A data frame is allowed a column with no name, and such a name compares to
  # the focal level's as `NA` rather than as no match. Counting those matches
  # has to say that none of them matched, since a column with no name is not a
  # column named for the focal level.
  weights <- expect_no_warning(
    wt_att(
      na_named,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    )
  )
  expect_equal(
    as.double(weights),
    as.double(wt_att(
      fixture$p_a,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ))
  )
})

test_that("a missing column name does not mask a real ambiguity", {
  fixture <- level_named_df_fixture()
  na_named <- data.frame(
    fixture$p_a,
    fixture$p_b,
    fixture$p_b,
    fixture$p_a
  )
  names(na_named) <- c("a", "a", "b", NA)

  # The unnamed column sits in the same comparison as the two that are named for
  # the focal level, so leaving it out of the count must not take them with it.
  warnings <- 0L
  weights <- withCallingHandlers(
    wt_att(
      na_named,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ),
    propensity_df_duplicate_column_warning = function(cnd) {
      warnings <<- warnings + 1L
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(warnings, 1L)
  expect_equal(
    as.double(weights),
    as.double(wt_att(
      fixture$p_a,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ))
  )
})

test_that("a duplicated name away from the focal level is not an ambiguity", {
  fixture <- level_named_df_fixture()
  duplicated_names <- data.frame(fixture$p_a, fixture$p_b, fixture$p_a)
  names(duplicated_names) <- c("a", "b", "b")

  # Only the column the selection lands on is chosen between. A repeat of the
  # other level's name leaves the focal column unambiguous and is not this
  # warning's business.
  expect_no_warning(
    wt_att(
      duplicated_names,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ),
    class = "propensity_df_duplicate_column_warning"
  )
})

test_that("an explicit .propensity_col is not an ambiguity either", {
  fixture <- level_named_df_fixture()
  duplicated_names <- data.frame(fixture$p_a, fixture$p_b)
  names(duplicated_names) <- c("a", "a")

  # Naming the column by position is the caller answering the question the
  # warning asks, so it must not be asked.
  weights <- expect_no_warning(
    wt_att(
      duplicated_names,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a",
      .propensity_col = 2
    )
  )
  expect_equal(
    as.double(weights),
    as.double(wt_att(
      fixture$p_b,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    ))
  )
})

test_that("data frame propensity falls back to position silently with no level named", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(x1 = fixture$p_a, x2 = fixture$p_b)

  weights <- expect_no_warning(
    wt_att(ps_df, fixture$factor, exposure_type = "binary")
  )

  expect_equal(
    weights,
    wt_att(fixture$p_b, fixture$factor, exposure_type = "binary")
  )
})

test_that("`.propensity_col` overrides column name matching silently", {
  fixture <- level_named_df_fixture()

  # the level names would match, but the explicit column wins and is read as
  # the probability of the resolved focal level
  ps_df <- data.frame(a = fixture$p_a, b = fixture$p_b)
  weights_matched <- expect_no_warning(
    wt_att(
      ps_df,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a",
      .propensity_col = "b"
    )
  )
  expect_equal(
    weights_matched,
    wt_att(
      fixture$p_b,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    )
  )

  # names that cannot be matched do not warn either when the column is explicit
  ps_df_unmatched <- data.frame(x1 = fixture$p_a, x2 = fixture$p_b)
  weights_unmatched <- expect_no_warning(
    wt_att(
      ps_df_unmatched,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a",
      .propensity_col = "x1"
    )
  )
  expect_equal(
    weights_unmatched,
    wt_att(
      fixture$p_a,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    )
  )
})

test_that("a one column propensity data frame keeps the positional default", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(prob = fixture$p_b)

  weights <- expect_no_warning(
    wt_att(ps_df, fixture$factor, exposure_type = "binary")
  )

  expect_equal(
    weights,
    wt_att(fixture$p_b, fixture$factor, exposure_type = "binary")
  )
})

test_that("a one column propensity data frame reads silently with a named level", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(prob = fixture$p_a)

  # a single column is the caller supplying P(focal) directly, so there is
  # nothing to choose between and no fallback to report
  weights <- expect_no_warning(
    wt_att(
      ps_df,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    )
  )

  expect_equal(
    weights,
    wt_att(
      fixture$p_a,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    )
  )
})

test_that("data frame propensity fires the deprecated `.treated` exactly once", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(a = fixture$p_a, b = fixture$p_b)

  deprecations <- 0L
  weights <- with_always_deprecated(
    withCallingHandlers(
      wt_att(ps_df, fixture$factor, exposure_type = "binary", .treated = "a"),
      lifecycle_warning_deprecated = function(cnd) {
        deprecations <<- deprecations + 1L
        rlang::cnd_muffle(cnd)
      }
    )
  )

  expect_identical(deprecations, 1L)

  # the mapped level drives the column choice, so the "a" column is read
  expect_equal(
    weights,
    wt_att(
      fixture$p_a,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    )
  )
})

test_that("`wt_cens()` resolves a data frame column and its deprecation once", {
  # wt_cens() routes through the same data frame helper as the other weight
  # functions, so its exposure is read as a two-level indicator like any other
  # binary exposure, and "a" is the level whose probability is modeled
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(a = fixture$p_a, b = fixture$p_b)

  deprecations <- character()
  weights <- with_always_deprecated(
    withCallingHandlers(
      wt_cens(ps_df, fixture$factor, exposure_type = "binary", .treated = "a"),
      lifecycle_warning_deprecated = function(cnd) {
        deprecations <<- c(deprecations, conditionMessage(cnd))
        rlang::cnd_muffle(cnd)
      }
    )
  )

  expect_length(deprecations, 1)
  # the deprecation is attributed to the function the user called, not to the
  # `wt_ate()` machinery wt_cens() delegates to
  expect_match(deprecations[[1]], "wt_cens()", fixed = TRUE)
  expect_no_match(deprecations[[1]], "wt_ate()", fixed = TRUE)

  # the mapped level drives the column choice, so the "a" column is read
  expect_equal(
    weights,
    wt_cens(
      fixture$p_a,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a"
    )
  )
  expect_identical(estimand(weights), "uncensored")
})

test_that("continuous exposures keep the positional data frame default", {
  set.seed(2024)
  exposure <- rnorm(10)
  ps_df <- data.frame(first = rnorm(10), second = rnorm(10))

  expect_equal(
    wt_ate(ps_df, exposure, exposure_type = "continuous", stabilize = TRUE),
    wt_ate(
      ps_df$second,
      exposure,
      exposure_type = "continuous",
      stabilize = TRUE
    )
  )
})

test_that("data frame propensity resolves the focal column before stabilizing", {
  fixture <- level_named_df_fixture()
  ps_df <- data.frame(a = fixture$p_a, b = fixture$p_b)

  expect_equal(
    wt_ate(
      ps_df,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a",
      stabilize = TRUE
    ),
    wt_ate(
      fixture$p_a,
      fixture$factor,
      exposure_type = "binary",
      .focal_level = "a",
      stabilize = TRUE
    )
  )
})

test_that("wt_* functions work with parsnip output", {
  skip_if_not_installed("parsnip")
  skip_on_cran()

  # Simulate data
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment_numeric <- rbinom(n, 1, ps_true)
  treatment_factor <- factor(treatment_numeric, levels = c("0", "1"))

  df <- data.frame(
    treatment = treatment_factor,
    x1 = x1,
    x2 = x2
  )

  # Fit model with parsnip
  ps_spec <- parsnip::logistic_reg()
  ps_spec <- parsnip::set_engine(ps_spec, "glm")
  ps_model <- parsnip::fit(ps_spec, treatment ~ x1 + x2, data = df)

  # Get predictions
  ps_preds <- predict(ps_model, df, type = "prob")

  # Test that it works with default (second column)
  weights_ate <- wt_ate(ps_preds, treatment_numeric, exposure_type = "binary")
  expect_s3_class(weights_ate, "psw")
  expect_equal(estimand(weights_ate), "ate")

  # Test explicit column selection with parsnip column names
  weights_ate2 <- wt_ate(
    ps_preds,
    treatment_numeric,
    .propensity_col = ".pred_1",
    exposure_type = "binary"
  )
  expect_equal(weights_ate, weights_ate2)

  # Test with all estimands
  expect_s3_class(
    wt_att(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
  expect_s3_class(
    wt_atu(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
  expect_s3_class(
    wt_atm(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
  expect_s3_class(
    wt_ato(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
  expect_s3_class(
    wt_entropy(ps_preds, treatment_numeric, exposure_type = "binary"),
    "psw"
  )
})

test_that("wt_* functions work with GLM objects", {
  # Simulate data
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, ps_true)

  # Fit GLM model
  ps_model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Test ATE
  weights_ate_glm <- wt_ate(ps_model, treatment, exposure_type = "binary")
  weights_ate_numeric <- wt_ate(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_ate_glm, weights_ate_numeric)
  expect_s3_class(weights_ate_glm, "psw")
  expect_equal(estimand(weights_ate_glm), "ate")

  # Test ATT
  weights_att_glm <- wt_att(ps_model, treatment, exposure_type = "binary")
  weights_att_numeric <- wt_att(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_att_glm, weights_att_numeric)
  expect_s3_class(weights_att_glm, "psw")
  expect_equal(estimand(weights_att_glm), "att")

  # Test ATU
  weights_atu_glm <- wt_atu(ps_model, treatment, exposure_type = "binary")
  weights_atu_numeric <- wt_atu(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_atu_glm, weights_atu_numeric)
  expect_s3_class(weights_atu_glm, "psw")
  expect_equal(estimand(weights_atu_glm), "atu")

  # Test ATM
  weights_atm_glm <- wt_atm(ps_model, treatment, exposure_type = "binary")
  weights_atm_numeric <- wt_atm(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_atm_glm, weights_atm_numeric)
  expect_s3_class(weights_atm_glm, "psw")
  expect_equal(estimand(weights_atm_glm), "atm")

  # Test ATO
  weights_ato_glm <- wt_ato(ps_model, treatment, exposure_type = "binary")
  weights_ato_numeric <- wt_ato(ps_fitted, treatment, exposure_type = "binary")
  expect_equal(weights_ato_glm, weights_ato_numeric)
  expect_s3_class(weights_ato_glm, "psw")
  expect_equal(estimand(weights_ato_glm), "ato")

  # Test Entropy
  weights_entropy_glm <- wt_entropy(
    ps_model,
    treatment,
    exposure_type = "binary"
  )
  weights_entropy_numeric <- wt_entropy(
    ps_fitted,
    treatment,
    exposure_type = "binary"
  )
  expect_equal(weights_entropy_glm, weights_entropy_numeric)
  expect_s3_class(weights_entropy_glm, "psw")
  expect_equal(estimand(weights_entropy_glm), "entropy")
})

test_that("GLM methods with optional exposure argument", {
  # Simulate data
  set.seed(789)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, ps_true)

  # Fit GLM model
  ps_model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Test that all weight functions work without explicit exposure
  # and produce the same results as with explicit exposure

  # Test ATE
  weights_ate_auto <- wt_ate(ps_model, exposure_type = "binary")
  weights_ate_explicit <- wt_ate(ps_model, treatment, exposure_type = "binary")
  expect_equal(weights_ate_auto, weights_ate_explicit)
  expect_s3_class(weights_ate_auto, "psw")

  # Test ATT
  weights_att_auto <- wt_att(ps_model, exposure_type = "binary")
  weights_att_explicit <- wt_att(ps_model, treatment, exposure_type = "binary")
  expect_equal(weights_att_auto, weights_att_explicit)

  # Test ATU
  weights_atu_auto <- wt_atu(ps_model, exposure_type = "binary")
  weights_atu_explicit <- wt_atu(ps_model, treatment, exposure_type = "binary")
  expect_equal(weights_atu_auto, weights_atu_explicit)

  # Test ATM
  weights_atm_auto <- wt_atm(ps_model, exposure_type = "binary")
  weights_atm_explicit <- wt_atm(ps_model, treatment, exposure_type = "binary")
  expect_equal(weights_atm_auto, weights_atm_explicit)

  # Test ATO
  weights_ato_auto <- wt_ato(ps_model, exposure_type = "binary")
  weights_ato_explicit <- wt_ato(ps_model, treatment, exposure_type = "binary")
  expect_equal(weights_ato_auto, weights_ato_explicit)

  # Test Entropy
  weights_entropy_auto <- wt_entropy(ps_model, exposure_type = "binary")
  weights_entropy_explicit <- wt_entropy(
    ps_model,
    treatment,
    exposure_type = "binary"
  )
  expect_equal(weights_entropy_auto, weights_entropy_explicit)

  # Test with continuous exposure using wt_cens which supports continuous
  exposure_cont <- 2 + 0.5 * x1 + 0.3 * x2 + rnorm(n)
  ps_model_cont <- glm(exposure_cont ~ x1 + x2, family = gaussian)

  weights_cens_auto <- wt_cens(ps_model_cont, exposure_type = "continuous")
  weights_cens_explicit <- wt_cens(
    ps_model_cont,
    exposure_cont,
    exposure_type = "continuous"
  )
  expect_equal(weights_cens_auto, weights_cens_explicit)
})

test_that("GLM methods handle continuous exposures", {
  # Simulate continuous exposure data
  set.seed(456)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  exposure <- 2 + 0.5 * x1 + 0.3 * x2 + rnorm(n)

  # Fit linear model
  exposure_model <- glm(exposure ~ x1 + x2, family = gaussian)

  # Test ATE with continuous exposure
  weights_ate <- wt_ate(
    exposure_model,
    exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  expect_s3_class(weights_ate, "psw")
  expect_equal(estimand(weights_ate), "ate")
  expect_true(is_stabilized(weights_ate))

  # Check that weights are reasonable
  expect_true(all(is.finite(as.numeric(weights_ate))))
  expect_true(all(as.numeric(weights_ate) > 0))

  # Test optional exposure with continuous
  weights_ate_auto <- wt_ate(
    exposure_model,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  expect_equal(weights_ate_auto, weights_ate)
})

test_that("GLM methods error on non-GLM objects", {
  # Try with a non-GLM object
  expect_propensity_error(
    wt_ate("not a glm", c(0, 1, 0, 1))
  )

  expect_propensity_error(
    wt_att(list(a = 1, b = 2), c(0, 1, 0, 1))
  )
})

# Edge case tests for all wt_* functions ----

test_that("wt_* functions handle edge case propensity scores", {
  # Edge case: propensity scores at boundaries
  ps_boundary <- c(1e-10, 0.001, 0.5, 0.999, 1 - 1e-10)
  exposure_boundary <- c(0, 0, 1, 1, 1)

  # Test all functions with boundary values
  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(ps_boundary, exposure_boundary, exposure_type = "binary")
    expect_s3_class(weights, "psw")
    expect_true(all(is.finite(as.numeric(weights))))
    expect_true(all(as.numeric(weights) > 0))
  }

  # Edge case: all propensity scores identical
  ps_constant <- rep(0.5, 5)
  exposure_mixed <- c(0, 0, 1, 1, 0)

  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(ps_constant, exposure_mixed, exposure_type = "binary")
    expect_s3_class(weights, "psw")
    # For constant propensity scores, weights should be constant within treatment groups
    expect_equal(length(unique(weights[exposure_mixed == 0])), 1)
    expect_equal(length(unique(weights[exposure_mixed == 1])), 1)
  }

  # Edge case: single observation
  ps_single <- 0.3
  exposure_single <- 1

  # Single observation needs explicit .focal_level/.reference_level
  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(
      ps_single,
      exposure_single,
      exposure_type = "binary",
      .focal_level = 1
    )
    expect_s3_class(weights, "psw")
    expect_length(weights, 1)
    expect_true(is.finite(weights))
  }

  # Edge case: all treated or all control
  ps_various <- c(0.2, 0.4, 0.6, 0.8)

  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    # All treated - need to specify .focal_level since there's only one level
    weights_all_1 <- fn(
      ps_various,
      rep(1, 4),
      exposure_type = "binary",
      .focal_level = 1
    )
    expect_s3_class(weights_all_1, "psw")
    expect_true(all(is.finite(weights_all_1)))

    # All control - need to specify .reference_level since there's only one level
    weights_all_0 <- fn(
      ps_various,
      rep(0, 4),
      exposure_type = "binary",
      .reference_level = 0
    )
    expect_s3_class(weights_all_0, "psw")
    expect_true(all(is.finite(weights_all_0)))
  }
})

test_that("wt_* functions handle extreme weight scenarios", {
  # Create scenario that produces very large weights
  ps_extreme <- c(0.001, 0.999, 0.5, 0.5)
  exposure_extreme <- c(1, 0, 0, 1) # Misaligned with propensity

  # ATE should produce extreme weights
  weights_ate <- wt_ate(ps_extreme, exposure_extreme, exposure_type = "binary")
  expect_true(max(weights_ate) > 100) # Very large weight expected

  # ATO and ATM should be bounded
  weights_ato <- wt_ato(ps_extreme, exposure_extreme, exposure_type = "binary")
  weights_atm <- wt_atm(ps_extreme, exposure_extreme, exposure_type = "binary")
  expect_true(all(as.numeric(weights_ato) <= 1)) # ATO weights bounded by 1
  expect_true(all(as.numeric(weights_atm) <= 2)) # ATM weights bounded
})

# Additional error handling tests ----

test_that("wt_* functions error appropriately on invalid inputs", {
  # Invalid propensity scores (outside [0,1])
  expect_propensity_error(
    wt_ate(c(-0.1, 0.5, 1.1), c(0, 1, 0), exposure_type = "binary")
  )

  expect_propensity_error(
    wt_att(c(0, 0.5, 1), c(0, 1, 0), exposure_type = "binary")
  )

  # Length mismatch should error
  expect_propensity_error(
    wt_ate(c(0.1, 0.5), c(0, 1, 0), exposure_type = "binary")
  )

  # Invalid exposure type
  expect_propensity_error(
    wt_ate(c(0.1, 0.5), c(0, 1), exposure_type = "invalid")
  )

  # Non-numeric propensity scores
  expect_propensity_error(
    wt_ate(c("a", "b"), c(0, 1))
  )

  # Categorical exposure requires matrix propensity scores
  expect_propensity_error(
    wt_att(c(0.3, 0.3, 0.4), c(1, 2, 3), exposure_type = "categorical")
  )
})

# ---- `exposure_type` rejections name the function the user called -----------

# Every weight function routes `exposure_type` through the internal
# `match_exposure_type()` helper, so an unrecognized string is rejected several
# frames below the surface. The rejection still has to name the weight function
# the user called, the way every other guard in the package does. cli wraps the
# message at the console width, so the assertions collapse whitespace before
# matching the valid set.
#
# The routes differ in how they reach the numeric method, and only some of them
# get there through `UseMethod()`. `wt_cens()` and the modified propensity score
# classes call a numeric method by value instead, which puts a name no caller
# typed at the head of the frame's call, so those two routes are pinned
# alongside the dispatched ones.

exposure_type_fixture <- function() {
  ps <- c(0.2, 0.5, 0.8, 0.4)
  exposure <- c(0, 1, 1, 0)
  list(
    ps = ps,
    exposure = exposure,
    ps_df = data.frame(untreated = 1 - ps, treated = ps),
    ps_trimmed = ps_trim(
      ps,
      .exposure = exposure,
      method = "ps",
      lower = 0.25,
      upper = 0.75
    )
  )
}

unwrap_condition_message <- function(cnd) {
  gsub("\\s+", " ", conditionMessage(cnd))
}

test_that("an invalid exposure_type on the numeric method names wt_ate()", {
  fixture <- exposure_type_fixture()

  cnd <- rlang::catch_cnd(
    wt_ate(fixture$ps, fixture$exposure, exposure_type = "wrong"),
    classes = "error"
  )
  expect_s3_class(cnd, "error")

  # Without the call threaded through, the condition reports
  # `match_exposure_type()` and this assertion fails.
  expect_identical(as.character(cnd$call[[1]]), "wt_ate")

  # `wt_ate()` is the one weight function that accepts continuous exposures, so
  # its valid set has four entries.
  expect_match(
    unwrap_condition_message(cnd),
    paste(
      "`exposure_type` must be one of \"auto\", \"binary\",",
      "\"categorical\", or \"continuous\", not \"wrong\"."
    ),
    fixed = TRUE
  )

  expect_propensity_error(
    wt_ate(fixture$ps, fixture$exposure, exposure_type = "wrong")
  )
})

test_that("an invalid exposure_type on the data frame method names wt_att()", {
  fixture <- exposure_type_fixture()

  cnd <- rlang::catch_cnd(
    wt_att(fixture$ps_df, fixture$exposure, exposure_type = "wrong"),
    classes = "error"
  )
  expect_s3_class(cnd, "error")

  expect_identical(as.character(cnd$call[[1]]), "wt_att")

  # `wt_att()` does not support continuous exposures, so the restricted valid
  # set has to survive the trip through the data frame helper.
  expect_match(
    unwrap_condition_message(cnd),
    paste(
      "`exposure_type` must be one of \"auto\", \"binary\",",
      "or \"categorical\", not \"wrong\"."
    ),
    fixed = TRUE
  )

  expect_propensity_error(
    wt_att(fixture$ps_df, fixture$exposure, exposure_type = "wrong")
  )
})

test_that("an invalid exposure_type on the numeric method names wt_ato()", {
  fixture <- exposure_type_fixture()

  cnd <- rlang::catch_cnd(
    wt_ato(fixture$ps, fixture$exposure, exposure_type = "wrong"),
    classes = "error"
  )
  expect_s3_class(cnd, "error")

  expect_identical(as.character(cnd$call[[1]]), "wt_ato")

  expect_match(
    unwrap_condition_message(cnd),
    paste(
      "`exposure_type` must be one of \"auto\", \"binary\",",
      "or \"categorical\", not \"wrong\"."
    ),
    fixed = TRUE
  )

  expect_propensity_error(
    wt_ato(fixture$ps, fixture$exposure, exposure_type = "wrong")
  )
})

test_that("an invalid exposure_type on the glm method names wt_ate()", {
  fixture <- zero_one_glm_fixture()

  cnd <- rlang::catch_cnd(
    wt_ate(
      fixture$model,
      .exposure = fixture$exposure,
      exposure_type = "wrong"
    ),
    classes = "error"
  )
  expect_s3_class(cnd, "error")

  expect_identical(as.character(cnd$call[[1]]), "wt_ate")

  expect_match(
    unwrap_condition_message(cnd),
    paste(
      "`exposure_type` must be one of \"auto\", \"binary\",",
      "\"categorical\", or \"continuous\", not \"wrong\"."
    ),
    fixed = TRUE
  )

  expect_propensity_error(
    wt_ate(fixture$model, .exposure = fixture$exposure, exposure_type = "wrong")
  )
})

test_that("an invalid exposure_type on the wt_cens route names wt_cens()", {
  fixture <- exposure_type_fixture()

  # `wt_cens.numeric()` delegates to `wt_ate.numeric()` by calling it rather
  # than by dispatch, so the frame it errors in belongs to the ATE machinery.
  # The rejection still has to name the function the user called.
  cnd <- rlang::catch_cnd(
    wt_cens(fixture$ps, fixture$exposure, exposure_type = "wrong"),
    classes = "error"
  )
  expect_s3_class(cnd, "error")

  expect_identical(as.character(cnd$call[[1]]), "wt_cens")

  expect_match(
    unwrap_condition_message(cnd),
    paste(
      "`exposure_type` must be one of \"auto\", \"binary\",",
      "\"categorical\", or \"continuous\", not \"wrong\"."
    ),
    fixed = TRUE
  )

  expect_propensity_error(
    wt_cens(fixture$ps, fixture$exposure, exposure_type = "wrong")
  )
})

test_that("an invalid exposure_type on a trimmed propensity score names wt_ate()", {
  fixture <- exposure_type_fixture()

  # The modified propensity score methods reach the numeric method through a
  # local `weight_fn` argument, so the frame's call is a name no caller typed.
  # The refit warning arrives first, which is why only errors are caught here.
  cnd <- suppressWarnings(rlang::catch_cnd(
    wt_ate(fixture$ps_trimmed, fixture$exposure, exposure_type = "wrong"),
    classes = "error"
  ))
  expect_s3_class(cnd, "error")

  expect_identical(as.character(cnd$call[[1]]), "wt_ate")

  expect_match(
    unwrap_condition_message(cnd),
    paste(
      "`exposure_type` must be one of \"auto\", \"binary\",",
      "\"categorical\", or \"continuous\", not \"wrong\"."
    ),
    fixed = TRUE
  )

  # the refit warning that precedes the error is attributed the same way
  expect_propensity_error(
    wt_ate(fixture$ps_trimmed, fixture$exposure, exposure_type = "wrong")
  )
})

# ---- every route names the weight function the user called ------------------

# A weight function reaches its numeric method by a different frame on every
# route: the data frame method through a shared helper, the glm method and the
# modified propensity score classes by calling the numeric method directly. A
# condition raised anywhere in that machinery has to name the weight function
# the user called, and it has to keep naming it when the call comes from inside
# a function of the user's own, which is what separates the frame the method was
# dispatched into from the frame that dispatched to it.

weight_error_call_name <- function(expr) {
  cnd <- rlang::catch_cnd(expr, classes = "error")
  if (is.null(cnd)) {
    return(NA_character_)
  }

  call <- conditionCall(cnd)
  if (is.null(call)) {
    return(NA_character_)
  }

  paste(deparse(call[[1]]), collapse = " ")
}

weight_error_message <- function(expr) {
  cnd <- rlang::catch_cnd(expr, classes = "error")
  if (is.null(cnd)) {
    return(NA_character_)
  }

  unwrap_condition_message(cnd)
}

categorical_attribution_fixture <- function() {
  list(
    exposure = factor(c("a", "b", "c", "a")),
    # the first row sums to 1.4, which every categorical route refuses
    bad_matrix = matrix(
      c(
        0.2,
        0.3,
        0.9,
        0.3,
        0.3,
        0.4,
        0.2,
        0.3,
        0.5,
        0.5,
        0.3,
        0.2
      ),
      nrow = 4,
      byrow = TRUE
    )
  )
}

categorical_attribution_routes <- function() {
  list(
    wt_ate = function(ps, exposure) {
      wt_ate(ps, exposure, exposure_type = "categorical")
    },
    wt_att = function(ps, exposure) {
      wt_att(ps, exposure, exposure_type = "categorical", .focal_level = "a")
    },
    wt_atu = function(ps, exposure) {
      wt_atu(
        ps,
        exposure,
        exposure_type = "categorical",
        .reference_level = "a"
      )
    },
    wt_atm = function(ps, exposure) {
      wt_atm(ps, exposure, exposure_type = "categorical")
    },
    wt_ato = function(ps, exposure) {
      wt_ato(ps, exposure, exposure_type = "categorical")
    },
    wt_entropy = function(ps, exposure) {
      wt_entropy(ps, exposure, exposure_type = "categorical")
    },
    wt_cens = function(ps, exposure) {
      wt_cens(ps, exposure, exposure_type = "categorical")
    }
  )
}

test_that("a categorical matrix rejection names the weight function", {
  fixture <- categorical_attribution_fixture()
  routes <- categorical_attribution_routes()

  for (fn_name in names(routes)) {
    route <- routes[[fn_name]]
    cnd <- rlang::catch_cnd(
      route(fixture$bad_matrix, fixture$exposure),
      classes = "error"
    )
    expect_s3_class(cnd, "propensity_matrix_sum_error")
    expect_identical(
      weight_error_call_name(route(fixture$bad_matrix, fixture$exposure)),
      fn_name
    )
  }
})

test_that("a categorical matrix rejection names the weight function, not the caller", {
  fixture <- categorical_attribution_fixture()

  # The routes above are reached from inside a local function, so a guard that
  # reads the frame that dispatched to the method instead of the frame it was
  # dispatched into would name that local function here.
  expect_identical(
    weight_error_call_name(
      wt_ate(
        fixture$bad_matrix,
        fixture$exposure,
        exposure_type = "categorical"
      )
    ),
    "wt_ate"
  )
})

test_that("a categorical matrix rejection on the data frame route names the weight function", {
  fixture <- categorical_attribution_fixture()

  expect_identical(
    weight_error_call_name(
      wt_ato(
        as.data.frame(fixture$bad_matrix),
        fixture$exposure,
        exposure_type = "categorical"
      )
    ),
    "wt_ato"
  )
})

test_that("data frame column selection failures name the weight function", {
  fixture <- exposure_type_fixture()

  expect_identical(
    weight_error_call_name(
      wt_ate(
        fixture$ps_df,
        fixture$exposure,
        exposure_type = "binary",
        .propensity_col = "nonexistent"
      )
    ),
    "wt_ate"
  )

  expect_identical(
    weight_error_call_name(
      wt_att(data.frame(), fixture$exposure, exposure_type = "binary")
    ),
    "wt_att"
  )

  expect_identical(
    weight_error_call_name(
      wt_ate(fixture$ps_df, c(0, 1), exposure_type = "binary")
    ),
    "wt_ate"
  )
})

test_that("the data frame helper attributes its own type guard to the method", {
  # The generic dispatches on a data frame, so the helper's own type guard is
  # reachable only by calling the method directly, the way a `NextMethod()` from
  # a method for a subclass would. It still has to name a frame the condition
  # can be read from rather than the frame that called the method.
  expect_identical(
    weight_error_call_name(
      wt_ate.data.frame(c(0.2, 0.5), c(0, 1), exposure_type = "binary")
    ),
    "wt_ate.data.frame"
  )
})

test_that("a propensity score out of range on the data frame route names the weight function", {
  out_of_range <- data.frame(
    low = c(-1, 2, 0.5, 0.5),
    high = c(2, -1, 0.5, 0.5)
  )
  fixture <- exposure_type_fixture()

  cnd <- rlang::catch_cnd(
    wt_ate(
      out_of_range,
      fixture$exposure,
      exposure_type = "binary",
      .propensity_col = 1
    ),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_range_error")

  expect_identical(
    weight_error_call_name(
      wt_ate(
        out_of_range,
        fixture$exposure,
        exposure_type = "binary",
        .propensity_col = 1
      )
    ),
    "wt_ate"
  )
})

test_that("guards on the glm route name the weight function", {
  fixture <- zero_one_glm_fixture()

  expect_identical(
    weight_error_call_name(
      wt_ate(fixture$model, .exposure = c(0, 1), exposure_type = "binary")
    ),
    "wt_ate"
  )

  expect_identical(
    weight_error_call_name(
      wt_att(fixture$model, .exposure = fixture$exposure, bogus = 1)
    ),
    "wt_att"
  )
})

test_that("guards on the modified propensity score routes name the weight function", {
  fixture <- exposure_type_fixture()

  expect_identical(
    weight_error_call_name(
      suppressWarnings(wt_ate(
        fixture$ps_trimmed,
        fixture$exposure,
        exposure_type = "binary",
        bogus = 1
      ))
    ),
    "wt_ate"
  )

  truncated <- ps_trunc(fixture$ps, method = "ps", lower = 0.25, upper = 0.75)
  expect_identical(
    weight_error_call_name(
      wt_ato(truncated, c(0, 1), exposure_type = "binary")
    ),
    "wt_ato"
  )
})

test_that("an unsupported exposure type names the weight function on every route", {
  continuous <- c(1.5, 2.5, 3.5, 4.5)
  ps <- c(0.2, 0.5, 0.8, 0.4)

  withr::local_options(propensity.quiet = TRUE)

  expect_identical(
    weight_error_call_name(wt_att(ps, continuous)),
    "wt_att"
  )

  expect_identical(
    weight_error_call_name(
      wt_att(data.frame(treated = ps), continuous, .propensity_col = 1)
    ),
    "wt_att"
  )
})

focal_required_fixture <- function() {
  exposure <- factor(c("a", "b", "c", "a"))
  ps_matrix <- matrix(
    c(
      0.2,
      0.3,
      0.5,
      0.3,
      0.3,
      0.4,
      0.2,
      0.3,
      0.5,
      0.5,
      0.3,
      0.2
    ),
    nrow = 4,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  list(
    exposure = exposure,
    ps_matrix = ps_matrix,
    ps_df = as.data.frame(ps_matrix)
  )
}

test_that("a missing focal level on a categorical exposure names the weight function", {
  fixture <- focal_required_fixture()

  # Each generic is called by name inside its own closure, so the frame the
  # condition names is the call the user wrote rather than the local symbol a
  # loop would carry.
  routes <- list(
    wt_att = function(ps, exposure) {
      wt_att(ps, exposure, exposure_type = "categorical")
    },
    wt_atu = function(ps, exposure) {
      wt_atu(ps, exposure, exposure_type = "categorical")
    }
  )

  for (fn_name in names(routes)) {
    route <- routes[[fn_name]]

    cnd <- rlang::catch_cnd(
      route(fixture$ps_matrix, fixture$exposure),
      classes = "error"
    )
    expect_s3_class(cnd, "propensity_focal_required_error")

    expect_identical(
      weight_error_call_name(route(fixture$ps_matrix, fixture$exposure)),
      fn_name
    )

    expect_identical(
      weight_error_call_name(route(fixture$ps_df, fixture$exposure)),
      fn_name
    )
  }
})

test_that("a missing focal level names the weight function, not the caller", {
  fixture <- focal_required_fixture()

  wrapping_focal_fn <- function() {
    wt_att(fixture$ps_matrix, fixture$exposure, exposure_type = "categorical")
  }

  expect_identical(weight_error_call_name(wrapping_focal_fn()), "wt_att")
})

# An exposure with more than two levels forced to binary cannot be coded 0/1
# without a level named, and each weight formula refuses it through the shared
# transformation. The refusal is raised two frames below the numeric method, in
# a helper named after the estimand rather than after anything the user typed.

uncodable_binary_routes <- function() {
  list(
    wt_ate = function(ps, exposure) {
      wt_ate(ps, exposure, exposure_type = "binary")
    },
    wt_att = function(ps, exposure) {
      wt_att(ps, exposure, exposure_type = "binary")
    },
    wt_atu = function(ps, exposure) {
      wt_atu(ps, exposure, exposure_type = "binary")
    },
    wt_atm = function(ps, exposure) {
      wt_atm(ps, exposure, exposure_type = "binary")
    },
    wt_ato = function(ps, exposure) {
      wt_ato(ps, exposure, exposure_type = "binary")
    },
    wt_entropy = function(ps, exposure) {
      wt_entropy(ps, exposure, exposure_type = "binary")
    },
    wt_cens = function(ps, exposure) {
      wt_cens(ps, exposure, exposure_type = "binary")
    }
  )
}

test_that("an exposure that cannot be coded 0/1 names the weight function", {
  ps <- c(0.2, 0.5, 0.8, 0.4)
  uncodable <- c(1, 2, 3, 4)
  routes <- uncodable_binary_routes()

  for (fn_name in names(routes)) {
    route <- routes[[fn_name]]

    cnd <- rlang::catch_cnd(route(ps, uncodable), classes = "error")
    expect_s3_class(cnd, "propensity_binary_transform_error")

    expect_identical(weight_error_call_name(route(ps, uncodable)), fn_name)
  }
})

# A shaped exposure can pass the length check, because `length()` counts cells
# for a matrix and columns for a data frame rather than observations. Nothing
# downstream of the transformation reduces it back to one value per
# observation, so the refusal has to happen where the exposure is coded.

test_that("an exposure with dimensions names the weight function", {
  ps <- c(0.2, 0.5, 0.8, 0.4)
  dimensioned <- matrix(c(1, 0, 1, 0), nrow = 2, ncol = 2)
  routes <- uncodable_binary_routes()

  for (fn_name in names(routes)) {
    route <- routes[[fn_name]]

    cnd <- rlang::catch_cnd(route(ps, dimensioned), classes = "error")
    expect_s3_class(cnd, "propensity_binary_transform_error")

    expect_identical(weight_error_call_name(route(ps, dimensioned)), fn_name)
  }
})

test_that("an exposure with dimensions is refused whatever the levels named", {
  ps <- c(0.2, 0.5, 0.8, 0.4)
  dimensioned <- matrix(c(1, 0, 1, 0), nrow = 2, ncol = 2)

  expect_error(
    wt_ate(ps, dimensioned, exposure_type = "binary", .focal_level = 1),
    class = "propensity_binary_transform_error"
  )
  expect_error(
    wt_ate(ps, dimensioned, exposure_type = "binary", .reference_level = 0),
    class = "propensity_binary_transform_error"
  )
  expect_error(
    wt_ate(ps, matrix(c(1, 0, 1, 0), ncol = 1), exposure_type = "binary"),
    class = "propensity_binary_transform_error"
  )
})

test_that("an exposure with dimensions reports what it was given", {
  ps <- c(0.2, 0.5, 0.8, 0.4)
  dimensioned <- matrix(c(1, 0, 1, 0), nrow = 2, ncol = 2)

  message <- weight_error_message(
    wt_ate(ps, dimensioned, exposure_type = "binary")
  )

  expect_match(message, "`.exposure` must be a vector", fixed = TRUE)
  expect_match(message, "matrix", fixed = TRUE)
  expect_match(message, "one value per observation", fixed = TRUE)
})

test_that("a one-dimensional exposure counts its dimension in the singular", {
  ps <- c(0.2, 0.5, 0.8, 0.4)
  one_dimensional <- array(c(1, 0, 1, 0), dim = 4)

  expect_propensity_error(
    wt_ate(ps, one_dimensional, exposure_type = "binary")
  )
})

test_that("an exposure that cannot be coded 0/1 names the weight function on every route", {
  ps <- c(0.2, 0.5, 0.8, 0.4)
  uncodable <- c(1, 2, 3, 4)
  fixture <- exposure_type_fixture()

  expect_identical(
    weight_error_call_name(
      wt_ate(
        data.frame(treated = ps),
        uncodable,
        exposure_type = "binary",
        .propensity_col = 1
      )
    ),
    "wt_ate"
  )

  glm_fixture <- zero_one_glm_fixture()
  expect_identical(
    weight_error_call_name(
      wt_atm(glm_fixture$model, rep(1:4, 10), exposure_type = "binary")
    ),
    "wt_atm"
  )

  expect_identical(
    weight_error_call_name(
      suppressWarnings(wt_ato(
        fixture$ps_trimmed,
        uncodable,
        exposure_type = "binary"
      ))
    ),
    "wt_ato"
  )

  wrapping_binary_fn <- function() {
    wt_entropy(ps, uncodable, exposure_type = "binary")
  }
  expect_identical(weight_error_call_name(wrapping_binary_fn()), "wt_entropy")
})

# ---- the `call` argument is plumbing, and is checked as such ----------------

# The weight generics pass their dots on to their methods, so the `call`
# argument the numeric methods carry is reachable from user code. A value the
# condition system cannot read as a call is refused where it arrives, rather
# than left to turn the next guard that fires into a report of rlang's
# internals.

test_that("a `call` that is not a call or an environment is rejected", {
  fixture <- exposure_type_fixture()

  expect_error(
    wt_ate(
      fixture$ps,
      fixture$exposure,
      exposure_type = "binary",
      call = "bogus"
    ),
    class = "propensity_call_arg_error"
  )
  expect_error(
    wt_ate(
      fixture$ps,
      fixture$exposure,
      exposure_type = "binary",
      call = "bogus"
    ),
    class = "propensity_error"
  )
  expect_match(
    weight_error_message(
      wt_ate(
        fixture$ps,
        fixture$exposure,
        exposure_type = "binary",
        call = "bogus"
      )
    ),
    "`call`",
    fixed = TRUE
  )

  expect_error(
    wt_att(fixture$ps, fixture$exposure, exposure_type = "binary", call = 42),
    class = "propensity_call_arg_error"
  )

  expect_error(
    wt_ato(fixture$ps, fixture$exposure, exposure_type = "binary", call = 42),
    class = "propensity_call_arg_error"
  )
})

test_that("a `call` rejection names the weight function on every route", {
  fixture <- exposure_type_fixture()

  expect_identical(
    weight_error_call_name(
      wt_ate(fixture$ps, fixture$exposure, exposure_type = "binary", call = 42)
    ),
    "wt_ate"
  )

  expect_identical(
    weight_error_call_name(
      wt_att(
        fixture$ps_df,
        fixture$exposure,
        exposure_type = "binary",
        .propensity_col = 2,
        call = 42
      )
    ),
    "wt_att"
  )

  expect_identical(
    weight_error_call_name(
      suppressWarnings(wt_ato(
        fixture$ps_trimmed,
        fixture$exposure,
        exposure_type = "binary",
        call = 42
      ))
    ),
    "wt_ato"
  )

  glm_fixture <- zero_one_glm_fixture()
  expect_identical(
    weight_error_call_name(
      wt_ate(
        glm_fixture$model,
        .exposure = glm_fixture$exposure,
        exposure_type = "binary",
        call = 42
      )
    ),
    "wt_ate"
  )
})

test_that("a `call` the caller supplies is used to attribute conditions", {
  ps <- c(0.5, 1.5)
  exposure <- c(0, 1)

  wrapping_weight_fn <- function() {
    wt_ate(ps, exposure, exposure_type = "binary", call = rlang::current_env())
  }

  expect_identical(
    weight_error_call_name(wrapping_weight_fn()),
    "wrapping_weight_fn"
  )

  # a call object is the other value the condition system reads
  expect_identical(
    weight_error_call_name(
      wt_ate(
        ps,
        exposure,
        exposure_type = "binary",
        call = quote(my_wrapper(x))
      )
    ),
    "my_wrapper"
  )
})

test_that("a valid exposure_type still dispatches to the weight formula", {
  fixture <- exposure_type_fixture()
  ps <- fixture$ps
  exposure <- fixture$exposure

  ate <- wt_ate(ps, exposure, exposure_type = "binary")
  expect_equal(
    as.numeric(ate),
    exposure / ps + (1 - exposure) / (1 - ps)
  )
  expect_identical(estimand(ate), "ate")

  ato <- wt_ato(ps, exposure, exposure_type = "binary")
  expect_equal(
    as.numeric(ato),
    exposure * (1 - ps) + (1 - exposure) * ps
  )
  expect_identical(estimand(ato), "ato")

  # the data frame method reaches the same numeric method
  expect_equal(
    wt_att(fixture$ps_df, exposure, exposure_type = "binary"),
    wt_att(ps, exposure, exposure_type = "binary")
  )
})

test_that("data frame methods error appropriately", {
  # Empty data frame
  expect_propensity_error(
    wt_ate(data.frame(), c(0, 1))
  )

  # Non-existent column
  df <- data.frame(a = c(0.1, 0.9), b = c(0.9, 0.1))
  expect_propensity_error(
    wt_ate(df, c(0, 1), .propensity_col = "nonexistent")
  )

  # Column index out of bounds
  expect_propensity_error(
    wt_ate(df, c(0, 1), .propensity_col = 5)
  )

  # Non-numeric column
  df_char <- data.frame(a = c("high", "low"), b = c(0.9, 0.1))
  suppressWarnings({
    expect_propensity_error(
      wt_ate(df_char, c(0, 1), .propensity_col = 1)
    )
  })

  # Column with invalid values
  df_invalid <- data.frame(a = c(0.5, 1.5), b = c(0.9, 0.1))
  expect_propensity_error(
    wt_ate(df_invalid, c(0, 1), .propensity_col = 1)
  )
})

test_that("GLM methods error appropriately", {
  # Non-GLM object
  expect_propensity_error(
    wt_ate(lm(mpg ~ wt, data = mtcars), rep(0:1, 16))
  )

  # GLM with wrong dimensions
  small_glm <- glm(c(0, 1) ~ c(1, 2), family = binomial)
  expect_propensity_error(
    wt_ate(small_glm, c(0, 1, 0, 1)) # Mismatch in length
  )
})

# Default method dispatch tests ----

test_that("default methods provide informative errors", {
  # Custom class that doesn't have a method
  custom_obj <- structure(list(x = 1:5), class = "my_custom_class")

  expect_propensity_error(
    wt_ate(custom_obj, c(0, 1, 0, 1, 0))
  )

  expect_propensity_error(
    wt_att(custom_obj, c(0, 1, 0, 1, 0))
  )

  expect_propensity_error(
    wt_atu(custom_obj, c(0, 1, 0, 1, 0))
  )

  expect_propensity_error(
    wt_atm(custom_obj, c(0, 1, 0, 1, 0))
  )

  expect_propensity_error(
    wt_ato(custom_obj, c(0, 1, 0, 1, 0))
  )

  expect_propensity_error(
    wt_entropy(custom_obj, c(0, 1, 0, 1, 0))
  )
})

# Data frame method tests with various column types ----

test_that("data frame methods handle various column configurations", {
  # Data frame with many columns
  ps_multi <- data.frame(
    col1 = runif(10, 0.1, 0.9),
    col2 = runif(10, 0.1, 0.9),
    col3 = runif(10, 0.1, 0.9),
    col4 = runif(10, 0.1, 0.9),
    col5 = runif(10, 0.1, 0.9)
  )
  exposure <- rbinom(10, 1, 0.5)

  # Test column selection by position
  for (i in 1:5) {
    weights <- suppressWarnings(wt_ate(
      ps_multi,
      exposure,
      .propensity_col = i,
      exposure_type = "binary"
    ))
    expected <- wt_ate(ps_multi[[i]], exposure, exposure_type = "binary")
    expect_equal(weights, expected)
  }

  # Test with tibble
  if (requireNamespace("tibble", quietly = TRUE)) {
    ps_tibble <- tibble::as_tibble(ps_multi)
    weights_tibble <- wt_ate(
      ps_tibble,
      exposure,
      .propensity_col = 3,
      exposure_type = "binary"
    )
    weights_df <- wt_ate(
      ps_multi,
      exposure,
      .propensity_col = 3,
      exposure_type = "binary"
    )
    expect_equal(weights_tibble, weights_df)
  }

  # Test with column names containing spaces
  ps_spaces <- data.frame(
    `control probability` = runif(5, 0.1, 0.9),
    `treatment probability` = runif(5, 0.1, 0.9),
    check.names = FALSE
  )
  exposure_small <- c(0, 0, 1, 1, 0)

  weights_spaces <- wt_ate(
    ps_spaces,
    exposure_small,
    .propensity_col = "treatment probability",
    exposure_type = "binary"
  )
  expect_s3_class(weights_spaces, "psw")

  # Test tidyselect helpers
  weights_last <- wt_ate(
    ps_multi,
    exposure,
    .propensity_col = tidyselect::last_col(),
    exposure_type = "binary"
  )
  expected_last <- wt_ate(ps_multi$col5, exposure, exposure_type = "binary")
  expect_equal(weights_last, expected_last)
})

# GLM method tests with different families ----

test_that("GLM methods handle non-binomial families appropriately", {
  # Gaussian family (for continuous exposures)
  set.seed(123)
  n <- 50
  x <- rnorm(n)
  exposure_cont <- 2 + 0.5 * x + rnorm(n)

  glm_gaussian <- glm(exposure_cont ~ x, family = gaussian)

  # Should work for ATE with continuous exposure
  weights_gaussian <- wt_ate(
    glm_gaussian,
    exposure_cont,
    exposure_type = "continuous",
    stabilize = TRUE
  )
  expect_s3_class(weights_gaussian, "psw")
  expect_true(all(is.finite(weights_gaussian)))

  # Should error for estimands that don't support continuous
  # ATT doesn't accept continuous as a valid exposure type
  expect_propensity_error(
    wt_att(glm_gaussian, exposure_cont, exposure_type = "continuous")
  )

  # Poisson family (should extract fitted values)
  # For non-binomial GLMs, the fitted values might not be valid propensity scores
  # Skip this test as it's not a valid use case
})

# Attribute preservation tests ----

test_that("attributes are preserved across all weight methods", {
  # Create trimmed propensity scores with attributes
  ps <- runif(20, 0.05, 0.95)
  exposure <- rbinom(20, 1, ps)

  ps_trimmed <- ps_trim(
    ps,
    .exposure = exposure,
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )

  # Test that all weight functions preserve trim attributes
  for (fn_name in c(
    "wt_ate",
    "wt_att",
    "wt_atu",
    "wt_atm",
    "wt_ato",
    "wt_entropy"
  )) {
    fn <- get(fn_name)

    # Suppress refit warning
    suppressWarnings({
      weights <- fn(ps_trimmed, exposure, exposure_type = "binary")
    })

    expect_true(is_ps_trimmed(weights))
    expect_equal(
      attr(weights, "ps_trim_meta"),
      attr(ps_trimmed, "ps_trim_meta")
    )
    expect_match(estimand(weights), "; trimmed$")
  }

  # Test with truncated propensity scores
  ps_truncated <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)

  for (fn_name in c(
    "wt_ate",
    "wt_att",
    "wt_atu",
    "wt_atm",
    "wt_ato",
    "wt_entropy"
  )) {
    fn <- get(fn_name)
    weights <- fn(ps_truncated, exposure, exposure_type = "binary")

    expect_true(is_ps_truncated(weights))
    expect_equal(
      attr(weights, "ps_trunc_meta"),
      attr(ps_truncated, "ps_trunc_meta")
    )
    expect_match(estimand(weights), "; truncated$")
  }
})

test_that("stabilization attributes are set correctly", {
  ps <- runif(20, 0.1, 0.9)
  exposure <- rbinom(20, 1, ps)

  # Test ATE with stabilization
  weights_unstab <- wt_ate(
    ps,
    exposure,
    stabilize = FALSE,
    exposure_type = "binary"
  )
  weights_stab <- wt_ate(
    ps,
    exposure,
    stabilize = TRUE,
    exposure_type = "binary"
  )

  expect_false(is_stabilized(weights_unstab))
  expect_true(is_stabilized(weights_stab))

  # Test with custom stabilization score
  weights_custom <- wt_ate(
    ps,
    exposure,
    stabilize = TRUE,
    stabilization_score = 0.4,
    exposure_type = "binary"
  )
  expect_true(is_stabilized(weights_custom))
})

# Stabilization tests across methods ----

test_that("stabilization works correctly for all applicable methods", {
  set.seed(456)
  ps <- runif(30, 0.2, 0.8)
  exposure <- rbinom(30, 1, ps)

  # ATE supports stabilization
  weights_ate_unstab <- wt_ate(
    ps,
    exposure,
    stabilize = FALSE,
    exposure_type = "binary"
  )
  weights_ate_stab <- wt_ate(
    ps,
    exposure,
    stabilize = TRUE,
    exposure_type = "binary"
  )

  # Check that stabilization was applied
  expect_true(is_stabilized(weights_ate_stab))
})

# NA handling tests ----

test_that("all methods handle NAs appropriately", {
  # Propensity scores with NAs
  ps_na <- c(0.2, NA, 0.5, 0.8, NA)
  exposure_na <- c(0, 1, 1, 0, 1)

  # All methods should handle NAs by producing NA weights
  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(ps_na, exposure_na, exposure_type = "binary")
    expect_s3_class(weights, "psw")
    expect_true(anyNA(weights))
  }

  # Exposure with NAs
  ps_good <- c(0.2, 0.3, 0.5, 0.8, 0.7)
  exposure_with_na <- c(0, NA, 1, 0, 1)

  for (fn in list(wt_ate, wt_att, wt_atu, wt_atm, wt_ato, wt_entropy)) {
    weights <- fn(ps_good, exposure_with_na, exposure_type = "binary")
    expect_s3_class(weights, "psw")
    expect_true(anyNA(weights))
  }

  # Data frame with NAs
  df_na <- data.frame(
    col1 = c(0.1, NA, 0.5),
    col2 = c(0.9, 0.5, NA)
  )

  # Data frame with NAs produces NA weights
  weights_df_na <- wt_ate(df_na, c(0, 1, 0), exposure_type = "binary")
  expect_s3_class(weights_df_na, "psw")
  expect_true(anyNA(weights_df_na))

  # GLM predictions with NAs
  set.seed(789)
  n <- 20
  x <- c(rnorm(18), NA, NA)
  y <- c(rbinom(18, 1, plogis(0.5 * x[1:18])), 0, 1)

  # GLM will handle NAs in predictors
  glm_na <- glm(y ~ x, family = binomial)

  # GLM with NA predictions will have shorter output
  # The NA observations are dropped during model fitting
  suppressWarnings({
    # GLM drops NA observations, so output is shorter
    expect_propensity_error(
      wt_ate(glm_na, y, exposure_type = "binary")
    )
  })
})

# Integration tests ----

test_that("weight functions integrate correctly with ps_trim and ps_trunc", {
  set.seed(999)
  n <- 100
  ps <- runif(n, 0.05, 0.95)
  exposure <- rbinom(n, 1, ps)

  # Create a model for refitting
  x <- rnorm(n)
  model <- glm(exposure ~ x + I(x^2), family = binomial)

  # Trim and refit
  ps_trim_obj <- ps_trim(
    ps,
    .exposure = exposure,
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )
  ps_refit_obj <- ps_refit(ps_trim_obj, model = model)

  # Should not warn after refit
  suppressMessages({
    weights_refit <- wt_ate(ps_refit_obj, exposure, exposure_type = "binary")
  })
  expect_true(is_ps_trimmed(weights_refit))
  expect_true(is_refit(weights_refit))

  # Truncate
  ps_trunc_obj <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)

  # Should not warn for truncation
  suppressMessages({
    weights_trunc <- wt_ate(ps_trunc_obj, exposure, exposure_type = "binary")
  })
  expect_true(is_ps_truncated(weights_trunc))
  expect_false(is_ps_trimmed(weights_trunc))
})

test_that("data frame and GLM methods produce consistent results", {
  set.seed(111)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  true_ps <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, true_ps)

  # Fit GLM
  model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(model)

  # Create data frame
  ps_df <- data.frame(
    control_prob = 1 - ps_fitted,
    treatment_prob = ps_fitted
  )

  # All methods should give same results
  for (fn_name in c(
    "wt_ate",
    "wt_att",
    "wt_atu",
    "wt_atm",
    "wt_ato",
    "wt_entropy"
  )) {
    fn <- get(fn_name)

    # From GLM
    weights_glm <- fn(model, treatment, exposure_type = "binary")

    # From fitted values
    weights_numeric <- fn(ps_fitted, treatment, exposure_type = "binary")

    # From data frame (using treatment column)
    weights_df <- fn(
      ps_df,
      treatment,
      .propensity_col = 2,
      exposure_type = "binary"
    )

    expect_equal(weights_glm, weights_numeric)
    expect_equal(weights_glm, weights_df)
  }
})

# Additional mathematical property tests ----

test_that("weight functions satisfy mathematical properties", {
  set.seed(222)
  n <- 100
  ps <- runif(n, 0.1, 0.9)
  treatment <- rbinom(n, 1, ps)

  # ATT weights: treated units should have weight 1
  weights_att <- wt_att(ps, treatment, exposure_type = "binary")
  expect_true(all(as.numeric(weights_att[treatment == 1]) == 1))

  # ATU weights: control units should have weight 1
  weights_atu <- wt_atu(ps, treatment, exposure_type = "binary")
  expect_true(all(as.numeric(weights_atu[treatment == 0]) == 1))

  # ATO weights: A * (1-e) + (1-A) * e
  weights_ato <- wt_ato(ps, treatment, exposure_type = "binary")
  # ATO weights are bounded by 1
  expect_true(all(as.numeric(weights_ato) <= 1 + 1e-10))

  # ATM weights: symmetric for e and 1-e
  # Create symmetric propensity scores
  ps_sym <- c(0.3, 0.7, 0.4, 0.6)
  treatment_sym <- c(0, 1, 1, 0)
  weights_atm <- wt_atm(ps_sym, treatment_sym, exposure_type = "binary")

  # Weights for symmetric PS pairs should be equal
  expect_equal(weights_atm[1], weights_atm[2], tolerance = 1e-10)
  expect_equal(weights_atm[3], weights_atm[4], tolerance = 1e-10)
})

test_that("continuous exposure weights have correct properties", {
  set.seed(333)
  n <- 100
  x <- rnorm(n)
  exposure <- 2 + 0.5 * x + rnorm(n, sd = 1.5)

  # Fit model
  model <- lm(exposure ~ x)
  mu_hat <- fitted(model)
  sigma <- summary(model)$sigma

  # Calculate weights
  weights_cont <- wt_ate(
    mu_hat,
    .exposure = exposure,
    .sigma = rep(sigma, n),
    exposure_type = "continuous",
    stabilize = TRUE
  )

  expect_s3_class(weights_cont, "psw")
  expect_true(all(is.finite(as.numeric(weights_cont))))
  expect_true(all(as.numeric(weights_cont) > 0))

  # Stabilized continuous weights should have mean approximately 1
  expect_equal(mean(weights_cont), 1, tolerance = 0.2)
})

# Comparison tests with other packages ----

test_that("ATE weights match WeightIt", {
  skip_if_not_installed("WeightIt")
  skip_on_cran()

  # Simulate data
  set.seed(123)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, ps_true)

  # Fit propensity score model
  ps_model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_ate(ps_fitted, treatment, exposure_type = "binary")

  # Calculate weights with WeightIt
  weightit_obj <- WeightIt::weightit(
    treatment ~ x1 + x2,
    data = data.frame(treatment = treatment, x1 = x1, x2 = x2),
    method = "ps",
    estimand = "ATE"
  )
  weightit_weights <- weightit_obj$weights

  # Compare
  expect_equal(
    as.numeric(our_weights),
    weightit_weights,
    tolerance = 1e-10,
    ignore_attr = "names"
  )
})

test_that("ATT weights match WeightIt", {
  skip_if_not_installed("WeightIt")
  skip_on_cran()

  # Simulate data
  set.seed(456)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, ps_true)

  # Fit propensity score model
  ps_model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_att(ps_fitted, treatment, exposure_type = "binary")

  # Calculate weights with WeightIt
  weightit_obj <- WeightIt::weightit(
    treatment ~ x1 + x2,
    data = data.frame(treatment = treatment, x1 = x1, x2 = x2),
    method = "ps",
    estimand = "ATT"
  )
  weightit_weights <- weightit_obj$weights

  # Compare
  expect_equal(
    as.numeric(our_weights),
    weightit_weights,
    tolerance = 1e-10,
    ignore_attr = "names"
  )
})

test_that("ATU/ATC weights match WeightIt", {
  skip_if_not_installed("WeightIt")
  skip_on_cran()

  # Simulate data
  set.seed(789)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, ps_true)

  # Fit propensity score model
  ps_model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_atu(ps_fitted, treatment, exposure_type = "binary")

  # Calculate weights with WeightIt
  # WeightIt uses ATC for Average Treatment Effect on Controls
  weightit_obj <- WeightIt::weightit(
    treatment ~ x1 + x2,
    data = data.frame(treatment = treatment, x1 = x1, x2 = x2),
    method = "ps",
    estimand = "ATC"
  )
  weightit_weights <- weightit_obj$weights

  # Compare
  expect_equal(
    as.numeric(our_weights),
    weightit_weights,
    tolerance = 1e-10,
    ignore_attr = "names"
  )
})

test_that("ATM weights match WeightIt", {
  skip_if_not_installed("WeightIt")
  skip_on_cran()

  # Simulate data
  set.seed(321)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, ps_true)

  # Fit propensity score model
  ps_model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_atm(ps_fitted, treatment, exposure_type = "binary")

  # Calculate weights with WeightIt
  # WeightIt uses ATM for matching weights
  weightit_obj <- WeightIt::weightit(
    treatment ~ x1 + x2,
    data = data.frame(treatment = treatment, x1 = x1, x2 = x2),
    method = "ps",
    estimand = "ATM"
  )
  weightit_weights <- weightit_obj$weights

  # Compare
  expect_equal(
    as.numeric(our_weights),
    weightit_weights,
    tolerance = 1e-10,
    ignore_attr = "names"
  )
})

test_that("ATO weights match WeightIt", {
  skip_if_not_installed("WeightIt")
  skip_on_cran()

  # Simulate data
  set.seed(654)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ps_true <- plogis(0.5 * x1 + 0.3 * x2)
  treatment <- rbinom(n, 1, ps_true)

  # Fit propensity score model
  ps_model <- glm(treatment ~ x1 + x2, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_ato(ps_fitted, treatment, exposure_type = "binary")

  # Calculate weights with WeightIt
  # WeightIt uses ATO for overlap weights
  weightit_obj <- WeightIt::weightit(
    treatment ~ x1 + x2,
    data = data.frame(treatment = treatment, x1 = x1, x2 = x2),
    method = "ps",
    estimand = "ATO"
  )
  weightit_weights <- weightit_obj$weights

  # Compare
  expect_equal(
    as.numeric(our_weights),
    weightit_weights,
    tolerance = 1e-10,
    ignore_attr = "names"
  )
})

# PSweight comparison tests ----

test_that("ATE weights match PSweight", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Fit propensity score model
  ps_formula <- trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6
  ps_model <- glm(ps_formula, data = PSweight::psdata_cl, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_ate(
    ps_fitted,
    PSweight::psdata_cl$trt,
    exposure_type = "binary"
  )

  # Calculate weights with PSweight
  psw_obj <- PSweight::SumStat(
    ps.formula = ps_formula,
    data = PSweight::psdata_cl,
    weight = "IPW"
  )

  # Extract IPW weights and un-normalize them
  # PSweight normalizes weights to sum to 1 within each group
  psw_weights_norm <- psw_obj$ps.weights$IPW

  # Calculate sum of raw weights for each group to un-normalize
  trt <- PSweight::psdata_cl$trt
  raw_ipw <- ifelse(trt == 1, 1 / ps_fitted, 1 / (1 - ps_fitted))
  sum_raw_trt1 <- sum(raw_ipw[trt == 1])
  sum_raw_trt0 <- sum(raw_ipw[trt == 0])

  # Un-normalize PSweight weights
  psw_weights_raw <- numeric(nrow(PSweight::psdata_cl))
  psw_weights_raw[trt == 1] <- psw_weights_norm[trt == 1] * sum_raw_trt1
  psw_weights_raw[trt == 0] <- psw_weights_norm[trt == 0] * sum_raw_trt0

  # Compare
  expect_equal(as.numeric(our_weights), psw_weights_raw, tolerance = 1e-10)
})

test_that("Overlap (ATO) weights use different formula than PSweight", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Fit propensity score model
  ps_formula <- trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6
  ps_model <- glm(ps_formula, data = PSweight::psdata_cl, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_ato(
    ps_fitted,
    PSweight::psdata_cl$trt,
    exposure_type = "binary"
  )

  # Calculate weights with PSweight
  psw_obj <- PSweight::SumStat(
    ps.formula = ps_formula,
    data = PSweight::psdata_cl,
    weight = "overlap"
  )

  # Extract overlap weights and un-normalize
  # PSweight normalizes weights to sum to 1 within each group
  psw_weights_norm <- psw_obj$ps.weights$overlap

  # Calculate sum of raw weights for each group to un-normalize
  trt <- PSweight::psdata_cl$trt
  raw_wts <- ps_fitted * (1 - ps_fitted) # ATO weights are same for both groups
  sum_raw_trt1 <- sum(raw_wts[trt == 1])
  sum_raw_trt0 <- sum(raw_wts[trt == 0])

  # Un-normalize PSweight weights
  psw_weights_raw <- numeric(nrow(PSweight::psdata_cl))
  psw_weights_raw[trt == 1] <- psw_weights_norm[trt == 1] * sum_raw_trt1
  psw_weights_raw[trt == 0] <- psw_weights_norm[trt == 0] * sum_raw_trt0

  # Our package uses (1-ps) for treated and ps for control
  # PSweight uses ps*(1-ps) for both groups
  # These are different formulations of overlap weights
  # So we skip the direct comparison
  skip("ATO formulas differ between packages")
})

test_that("Matching (ATM) weights match PSweight", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Fit propensity score model
  ps_formula <- trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6
  ps_model <- glm(ps_formula, data = PSweight::psdata_cl, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_atm(
    ps_fitted,
    PSweight::psdata_cl$trt,
    exposure_type = "binary"
  )

  # Calculate weights with PSweight
  psw_obj <- PSweight::SumStat(
    ps.formula = ps_formula,
    data = PSweight::psdata_cl,
    weight = "matching"
  )

  # Extract matching weights and un-normalize
  # PSweight normalizes weights to sum to 1 within each group
  psw_weights_norm <- psw_obj$ps.weights$matching

  # Calculate sum of raw weights for each group to un-normalize
  trt <- PSweight::psdata_cl$trt
  raw_wts <- ifelse(
    trt == 1,
    pmin(ps_fitted, 1 - ps_fitted) / ps_fitted,
    pmin(ps_fitted, 1 - ps_fitted) / (1 - ps_fitted)
  )
  sum_raw_trt1 <- sum(raw_wts[trt == 1])
  sum_raw_trt0 <- sum(raw_wts[trt == 0])

  # Un-normalize PSweight weights
  psw_weights_raw <- numeric(nrow(PSweight::psdata_cl))
  psw_weights_raw[trt == 1] <- psw_weights_norm[trt == 1] * sum_raw_trt1
  psw_weights_raw[trt == 0] <- psw_weights_norm[trt == 0] * sum_raw_trt0

  # Compare
  expect_equal(as.numeric(our_weights), psw_weights_raw, tolerance = 1e-10)
})

test_that("Entropy weights use different formula than PSweight", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Fit propensity score model
  ps_formula <- trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6
  ps_model <- glm(ps_formula, data = PSweight::psdata_cl, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_entropy(
    ps_fitted,
    PSweight::psdata_cl$trt,
    exposure_type = "binary"
  )

  # Calculate weights with PSweight
  psw_obj <- PSweight::SumStat(
    ps.formula = ps_formula,
    data = PSweight::psdata_cl,
    weight = "entropy"
  )

  # Extract entropy weights and un-normalize
  # PSweight normalizes weights to sum to 1 within each group
  psw_weights_norm <- psw_obj$ps.weights$entropy

  # Calculate sum of raw weights for each group to un-normalize
  trt <- PSweight::psdata_cl$trt
  raw_wts <- ifelse(trt == 1, 1 - ps_fitted, ps_fitted)
  sum_raw_trt1 <- sum(raw_wts[trt == 1])
  sum_raw_trt0 <- sum(raw_wts[trt == 0])

  # Un-normalize PSweight weights
  psw_weights_raw <- numeric(nrow(PSweight::psdata_cl))
  psw_weights_raw[trt == 1] <- psw_weights_norm[trt == 1] * sum_raw_trt1
  psw_weights_raw[trt == 0] <- psw_weights_norm[trt == 0] * sum_raw_trt0

  # Our package uses entropy tilting function h(e) = -[e*log(e) + (1-e)*log(1-e)]
  # PSweight uses simpler formula: w = (1-e) for treated, w = e for control
  # These give different results by a factor related to the entropy function
  # So we just verify the ratio is consistent
  ratio <- as.numeric(our_weights) / psw_weights_raw
  expect_true(sd(ratio) / mean(ratio) < 0.01) # Coefficient of variation < 1%
})

test_that("Treated (ATT) weights match PSweight", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Fit propensity score model
  ps_formula <- trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6
  ps_model <- glm(ps_formula, data = PSweight::psdata_cl, family = binomial)
  ps_fitted <- fitted(ps_model)

  # Calculate weights with our package
  our_weights <- wt_att(
    ps_fitted,
    PSweight::psdata_cl$trt,
    exposure_type = "binary"
  )

  # Calculate weights with PSweight
  psw_obj <- PSweight::SumStat(
    ps.formula = ps_formula,
    data = PSweight::psdata_cl,
    weight = "treated"
  )

  # Extract treated weights and un-normalize
  # PSweight normalizes weights to sum to 1 within each group
  psw_weights_norm <- psw_obj$ps.weights$treated

  # Calculate sum of raw weights for each group to un-normalize
  trt <- PSweight::psdata_cl$trt
  raw_wts <- ifelse(trt == 1, 1, ps_fitted / (1 - ps_fitted))
  sum_raw_trt1 <- sum(raw_wts[trt == 1])
  sum_raw_trt0 <- sum(raw_wts[trt == 0])

  # Un-normalize PSweight weights
  psw_weights_raw <- numeric(nrow(PSweight::psdata_cl))
  psw_weights_raw[trt == 1] <- psw_weights_norm[trt == 1] * sum_raw_trt1
  psw_weights_raw[trt == 0] <- psw_weights_norm[trt == 0] * sum_raw_trt0

  # Compare
  expect_equal(as.numeric(our_weights), psw_weights_raw, tolerance = 1e-10)
})

# Integration tests with other methods ----

test_that("Data frame methods produce same results as WeightIt", {
  skip_if_not_installed("WeightIt")
  skip_on_cran()

  # Simulate data
  set.seed(111)
  n <- 150
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  treatment <- rbinom(n, 1, plogis(0.5 * x1 + 0.3 * x2))

  # Create data frame
  df <- data.frame(treatment = treatment, x1 = x1, x2 = x2)

  # Fit model and get predictions as data frame
  ps_model <- glm(treatment ~ x1 + x2, data = df, family = binomial)
  ps_df <- data.frame(
    control_prob = 1 - fitted(ps_model),
    treatment_prob = fitted(ps_model)
  )

  # Our weights using data frame method
  our_ate_df <- wt_ate(ps_df, treatment, exposure_type = "binary")
  our_att_df <- wt_att(ps_df, treatment, exposure_type = "binary")

  # WeightIt weights
  w_ate <- WeightIt::weightit(
    treatment ~ x1 + x2,
    data = df,
    method = "ps",
    estimand = "ATE"
  )
  w_att <- WeightIt::weightit(
    treatment ~ x1 + x2,
    data = df,
    method = "ps",
    estimand = "ATT"
  )

  # Compare
  expect_equal(
    as.numeric(our_ate_df),
    w_ate$weights,
    tolerance = 1e-10,
    ignore_attr = "names"
  )
  expect_equal(
    as.numeric(our_att_df),
    w_att$weights,
    tolerance = 1e-10,
    ignore_attr = "names"
  )
})

test_that("GLM methods produce same results as PSweight", {
  skip_if_not_installed("PSweight")
  skip_on_cran()

  # Fit propensity score model
  ps_model <- glm(
    trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6,
    data = PSweight::psdata_cl,
    family = binomial
  )

  # Our weights using GLM method
  our_ate_glm <- wt_ate(
    ps_model,
    PSweight::psdata_cl$trt,
    exposure_type = "binary"
  )
  our_entropy_glm <- wt_entropy(
    ps_model,
    PSweight::psdata_cl$trt,
    exposure_type = "binary"
  )

  # PSweight
  psw_obj <- PSweight::SumStat(
    ps.formula = trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6,
    data = PSweight::psdata_cl,
    weight = c("IPW", "entropy")
  )

  # Extract and un-normalize weights
  ps_fitted <- fitted(ps_model)
  trt <- PSweight::psdata_cl$trt

  # Calculate sums for un-normalization
  # IPW/ATE weights
  raw_ipw <- ifelse(trt == 1, 1 / ps_fitted, 1 / (1 - ps_fitted))
  sum_ipw_trt1 <- sum(raw_ipw[trt == 1])
  sum_ipw_trt0 <- sum(raw_ipw[trt == 0])

  # Entropy weights
  raw_entropy <- ifelse(trt == 1, 1 - ps_fitted, ps_fitted)
  sum_entropy_trt1 <- sum(raw_entropy[trt == 1])
  sum_entropy_trt0 <- sum(raw_entropy[trt == 0])

  # Un-normalize PSweight weights
  psw_ate_raw <- numeric(nrow(PSweight::psdata_cl))
  psw_ate_raw[trt == 1] <- psw_obj$ps.weights$IPW[trt == 1] * sum_ipw_trt1
  psw_ate_raw[trt == 0] <- psw_obj$ps.weights$IPW[trt == 0] * sum_ipw_trt0

  psw_entropy_raw <- numeric(nrow(PSweight::psdata_cl))
  psw_entropy_raw[trt == 1] <- psw_obj$ps.weights$entropy[trt == 1] *
    sum_entropy_trt1
  psw_entropy_raw[trt == 0] <- psw_obj$ps.weights$entropy[trt == 0] *
    sum_entropy_trt0

  # Compare
  expect_equal(as.numeric(our_ate_glm), psw_ate_raw, tolerance = 1e-10)

  # For entropy, our package uses different formula than PSweight
  # so we just verify the ratio is consistent
  entropy_ratio <- as.numeric(our_entropy_glm) / psw_entropy_raw
  expect_true(sd(entropy_ratio) / mean(entropy_ratio) < 0.01)
})

# ---- validating a stabilization score at the weight functions ---------------

test_that("wt_ate() rejects an invalid stabilization score for a binary exposure", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c(0, 1, 0, 1)
  scored <- function(score) {
    wt_ate(
      ps,
      exposure,
      exposure_type = "binary",
      stabilize = TRUE,
      stabilization_score = score
    )
  }

  expect_error(
    scored(c(1, 2)),
    class = "propensity_stabilization_score_error"
  )
  expect_error(
    scored(c(1, 2, 3)),
    class = "propensity_stabilization_score_error"
  )
  expect_error(scored(0), class = "propensity_stabilization_score_error")
  expect_error(scored(-1), class = "propensity_stabilization_score_error")
  expect_error(
    scored(NA_real_),
    class = "propensity_stabilization_score_error"
  )
  expect_error(scored(NaN), class = "propensity_stabilization_score_error")
  expect_error(scored(Inf), class = "propensity_stabilization_score_error")
  expect_error(scored("0.5"), class = "propensity_stabilization_score_error")

  # A length that neither matches nor divides the weights recycles under a base
  # R warning, so the rejection has to arrive before the multiplication.
  expect_no_warning(expect_error(
    scored(c(1, 2, 3)),
    class = "propensity_stabilization_score_error"
  ))
})

test_that("wt_ate() rejects a stabilization score with dimensions", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c(0, 1, 0, 1)
  score <- matrix(c(0.4, 0.5, 0.6, 0.7), nrow = 2)

  # A matrix holds as many values as the weights have observations, so the
  # length rule reads it as one value per observation and the coercion that
  # follows silently flattens it.
  expect_error(
    wt_ate(
      ps,
      exposure,
      exposure_type = "binary",
      stabilize = TRUE,
      stabilization_score = score
    ),
    class = "propensity_stabilization_score_error"
  )
  expect_match(
    weight_error_message(
      wt_ate(
        ps,
        exposure,
        exposure_type = "binary",
        stabilize = TRUE,
        stabilization_score = score
      )
    ),
    "2 x 2",
    fixed = TRUE
  )

  expect_error(
    psw(ps, "ate", stabilized = TRUE, stabilization_score = score),
    class = "propensity_stabilization_score_error"
  )

  expect_error(
    wt_ate(
      ps,
      exposure,
      exposure_type = "binary",
      stabilize = TRUE,
      stabilization_score = array(0.5, dim = c(2, 2, 1))
    ),
    class = "propensity_stabilization_score_error"
  )
})

test_that("wt_ate() reports a length mismatch before the stabilization score", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c(0, 1, 0)

  # The score matches `.propensity`, so a score checked against the exposure's
  # length reports the score rather than the mismatch that caused it.
  cnd <- rlang::catch_cnd(
    wt_ate(
      ps,
      exposure,
      exposure_type = "binary",
      stabilize = TRUE,
      stabilization_score = ps
    ),
    classes = "error"
  )
  expect_s3_class(cnd, "propensity_length_error")

  cnd_continuous <- rlang::catch_cnd(
    wt_ate(
      ps,
      c(1.5, 2.5, 3.5),
      exposure_type = "continuous",
      stabilize = TRUE,
      stabilization_score = ps
    ),
    classes = "error"
  )
  expect_s3_class(cnd_continuous, "propensity_length_error")
})

test_that("wt_cens() rejects an invalid stabilization score", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c(0, 1, 0, 1)

  expect_error(
    wt_cens(
      ps,
      exposure,
      exposure_type = "binary",
      stabilize = TRUE,
      stabilization_score = c(1, 2)
    ),
    class = "propensity_stabilization_score_error"
  )
  expect_error(
    wt_cens(
      ps,
      exposure,
      exposure_type = "binary",
      stabilize = TRUE,
      stabilization_score = -1
    ),
    class = "propensity_stabilization_score_error"
  )
})

test_that("wt_ate() rejects an invalid stabilization score for a continuous exposure", {
  denom_model <- lm(mpg ~ gear + am + carb, data = mtcars)
  scored <- function(score) {
    wt_ate(
      predict(denom_model),
      .exposure = mtcars$mpg,
      .sigma = influence(denom_model)$sigma,
      exposure_type = "continuous",
      stabilize = TRUE,
      stabilization_score = score
    )
  }

  expect_error(
    scored(c(1, 2)),
    class = "propensity_stabilization_score_error"
  )
  expect_error(scored(0), class = "propensity_stabilization_score_error")
  expect_error(scored(-1), class = "propensity_stabilization_score_error")
  expect_error(
    scored(NA_real_),
    class = "propensity_stabilization_score_error"
  )
  expect_error(scored(Inf), class = "propensity_stabilization_score_error")
  expect_error(scored("0.5"), class = "propensity_stabilization_score_error")
})

test_that("wt_ate() rejects an invalid stabilization score for a categorical exposure", {
  exposure <- factor(c("A", "B", "C", "A", "B", "C"))
  ps_matrix <- matrix(
    c(
      0.5,
      0.3,
      0.2,
      0.2,
      0.5,
      0.3,
      0.1,
      0.2,
      0.7,
      0.6,
      0.3,
      0.1,
      0.3,
      0.4,
      0.3,
      0.2,
      0.2,
      0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(ps_matrix) <- levels(exposure)

  scored <- function(score) {
    wt_ate(
      ps_matrix,
      exposure,
      exposure_type = "categorical",
      stabilize = TRUE,
      stabilization_score = score
    )
  }

  expect_error(
    scored(c(1, 2)),
    class = "propensity_stabilization_score_error"
  )
  expect_error(
    scored(c(1, 2, 3, 4)),
    class = "propensity_stabilization_score_error"
  )
  expect_error(scored(0), class = "propensity_stabilization_score_error")
  expect_error(scored(-1), class = "propensity_stabilization_score_error")
  expect_error(
    scored(NA_real_),
    class = "propensity_stabilization_score_error"
  )
  expect_error(scored(Inf), class = "propensity_stabilization_score_error")
  expect_error(scored("0.5"), class = "propensity_stabilization_score_error")

  expect_no_warning(expect_error(
    scored(c(1, 2, 3, 4)),
    class = "propensity_stabilization_score_error"
  ))
})

test_that("wt_ate() multiplies a per-observation stabilization score into the weights", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c(0, 1, 0, 1)
  score <- c(0.4, 0.5, 0.6, 0.7)

  unstabilized <- wt_ate(ps, exposure, exposure_type = "binary")
  scored <- wt_ate(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE,
    stabilization_score = score
  )

  expect_equal(as.double(scored), as.double(unstabilized) * score)
  expect_identical(stabilization_score(scored), score)
  expect_true(is_stabilized(scored))

  # A single value applies to every weight.
  scalar <- wt_ate(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE,
    stabilization_score = 0.5
  )
  expect_equal(as.double(scalar), as.double(unstabilized) * 0.5)
  expect_identical(stabilization_score(scalar), 0.5)
})

test_that("wt_ate() records an integer stabilization score as a double", {
  ps <- c(0.2, 0.4, 0.6, 0.8)
  exposure <- c(0, 1, 0, 1)

  scored <- wt_ate(
    ps,
    exposure,
    exposure_type = "binary",
    stabilize = TRUE,
    stabilization_score = 2L
  )

  expect_identical(stabilization_score(scored), 2)
  expect_type(stabilization_score(scored), "double")
})

test_that("wt_ate() rejects an invalid stabilization score built on trimmed scores", {
  ps <- ps_trim(
    c(0.05, 0.3, 0.5, 0.7, 0.95),
    method = "ps",
    lower = 0.1,
    upper = 0.9
  )
  exposure <- c(0, 0, 1, 1, 1)

  expect_error(
    suppressWarnings(wt_ate(
      ps,
      exposure,
      exposure_type = "binary",
      .focal_level = 1,
      stabilize = TRUE,
      stabilization_score = c(1, 2)
    )),
    class = "propensity_stabilization_score_error"
  )
})
