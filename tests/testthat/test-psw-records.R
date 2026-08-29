# What a `psw` records about the exposure it was built for: `exposure_type` on
# every set of weights, and `density_meta` on the continuous weights whose value
# is a ratio of densities. Both describe the exposure rather than the units, so
# neither follows a length rule, and both are read back through the exported
# accessors `exposure_type()` and `density_meta()`.

# A continuous exposure and the fitted conditional means of a model for it. The
# weights are a ratio of densities in these two, so neither has to be a
# probability.
records_exposure <- c(0.15, 0.32, 0.48, 0.61, 0.74, 0.88)
records_mu <- c(0.2, 0.3, 0.5, 0.6, 0.7, 0.85)

# Continuous ATE weights, stabilized on the marginal density unless the caller
# asks for something else. The arguments the record is written from are named
# rather than passed through dots so that a fixture cannot collide with them.
continuous_records_psw <- function(
  stabilize = TRUE,
  stabilization_score = NULL,
  .sigma = NULL
) {
  wt_ate(
    records_mu,
    records_exposure,
    .sigma = .sigma,
    exposure_type = "continuous",
    stabilize = stabilize,
    stabilization_score = stabilization_score
  )
}

records_binary_ps <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
records_binary_exposure <- c(0, 1, 0, 1, 0, 1)

records_categorical_exposure <- factor(c("A", "B", "C", "A", "B", "C"))

records_categorical_ps <- local({
  ps <- matrix(
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
  colnames(ps) <- levels(records_categorical_exposure)
  ps
})

# `wt_*()` on a trimmed propensity score warns that the model was never refit.
# These tests are about what the weights record, so that one warning is muffled
# by class rather than silencing the whole call.
muffle_refit_warning <- function(expr) {
  withCallingHandlers(
    expr,
    propensity_no_refit_warning = function(cnd) invokeRestart("muffleWarning")
  )
}

# Every warning raised while evaluating `expr`, by class, so a case that must
# report a disagreement once can assert that it reported it once and reported
# nothing else.
records_warning_classes <- function(expr) {
  classes <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(cnd) {
      classes <<- c(classes, class(cnd)[[1]])
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, classes = classes)
}

# A density record holds the specification that evaluates the density, and a
# function written the same way in two calls is not the same object, so a record
# built in one call is compared against one built in another by what it says
# rather than by identity.
density_meta_summary <- function(x) {
  meta <- density_meta(x)
  if (is.null(meta)) {
    return(NULL)
  }

  list(
    density = format(meta$density),
    numerator = meta$numerator,
    sigma = meta$sigma
  )
}

marginal_normal_record <- list(
  density = "normal",
  numerator = "marginal",
  sigma = "pooled"
)

# ---- the accessors ----------------------------------------------------------

test_that("the record accessors read what the weights carry", {
  w <- continuous_records_psw()

  expect_identical(exposure_type(w), attr(w, "exposure_type"))
  expect_identical(density_meta(w), attr(w, "density_meta"))

  expect_identical(exposure_type(w), "continuous")
  expect_type(density_meta(w), "list")
  expect_named(
    density_meta(w),
    c("density", "numerator", "sigma", "sigma_value")
  )
})

test_that("the record accessors report weights that carry nothing", {
  # Both records describe an exposure the weight functions were given. Weights
  # written by hand were given no exposure to describe, so there is nothing for
  # either record to say about them.
  w <- psw(c(1, 2, 3), estimand = "ate")

  expect_null(exposure_type(w))
  expect_null(density_meta(w))
})

# ---- what each weight function records --------------------------------------

test_that("every binary weight function records the exposure type", {
  weight_fns <- list(
    wt_ate = wt_ate,
    wt_att = wt_att,
    wt_atu = wt_atu,
    wt_atm = wt_atm,
    wt_ato = wt_ato,
    wt_entropy = wt_entropy,
    wt_cens = wt_cens
  )

  for (name in names(weight_fns)) {
    w <- weight_fns[[name]](
      records_binary_ps,
      records_binary_exposure,
      exposure_type = "binary",
      .focal_level = 1
    )

    expect_identical(exposure_type(w), "binary", label = name)

    # A binary exposure is weighted by propensity scores rather than by a
    # density, so there is no density for the second record to describe.
    expect_null(density_meta(w), label = name)
  }
})

test_that("every categorical weight function records the exposure type", {
  unfocused <- list(
    wt_ate = wt_ate,
    wt_atm = wt_atm,
    wt_ato = wt_ato,
    wt_entropy = wt_entropy
  )

  for (name in names(unfocused)) {
    w <- unfocused[[name]](
      records_categorical_ps,
      records_categorical_exposure,
      exposure_type = "categorical"
    )

    expect_identical(exposure_type(w), "categorical", label = name)
    expect_null(density_meta(w), label = name)
  }

  focused <- list(wt_att = wt_att, wt_atu = wt_atu)

  for (name in names(focused)) {
    w <- focused[[name]](
      records_categorical_ps,
      records_categorical_exposure,
      .focal_level = "A",
      exposure_type = "categorical"
    )

    expect_identical(exposure_type(w), "categorical", label = name)
    expect_null(density_meta(w), label = name)
  }
})

test_that("continuous ate weights record the density they were built from", {
  w <- continuous_records_psw()

  expect_identical(exposure_type(w), "continuous")
  expect_s3_class(density_meta(w)$density, "propensity_density")
  expect_identical(density_meta_summary(w), marginal_normal_record)
})

test_that("continuous weights record which numerator stabilized them", {
  # Stabilizing on the marginal density and stabilizing on a score the user
  # supplied are different weights, and weights that were not stabilized at all
  # are a third. The record names which of the three these are, so that a reader
  # of the weights, and the sandwich that rebuilds them, can tell them apart.
  marginal <- continuous_records_psw()
  expect_identical(density_meta(marginal)$numerator, "marginal")

  scored <- continuous_records_psw(stabilization_score = 0.5)
  expect_identical(density_meta(scored)$numerator, "score")

  unstabilized <- continuous_records_psw(stabilize = FALSE)
  expect_identical(density_meta(unstabilized)$numerator, "none")

  # The exposure is continuous whichever numerator the weights took.
  expect_identical(exposure_type(unstabilized), "continuous")
  expect_identical(density_meta(unstabilized)$sigma, "pooled")
})

test_that("continuous weights record where the residual spread came from", {
  pooled <- continuous_records_psw()
  expect_identical(density_meta(pooled)$sigma, "pooled")

  supplied <- continuous_records_psw(.sigma = 0.12)
  expect_identical(density_meta(supplied)$sigma, "supplied")

  # A spread supplied for each observation is supplied just as a single one is.
  per_observation <- continuous_records_psw(
    .sigma = c(0.1, 0.12, 0.14, 0.11, 0.13, 0.15)
  )
  expect_identical(density_meta(per_observation)$sigma, "supplied")
})

test_that("continuous weights record a single supplied spread as its value", {
  # A spread that is one number describes the conditional density of every unit,
  # so it is a constant the ratio can be rebuilt from and the record keeps it.
  # A pooled spread is a function of the data rather than a number the caller
  # chose, and a spread supplied per observation is a different density for each
  # unit; neither is one constant, so neither leaves a value behind.
  expect_null(density_meta(continuous_records_psw())$sigma_value)

  supplied <- continuous_records_psw(.sigma = 0.12)
  expect_identical(density_meta(supplied)$sigma_value, 0.12)

  per_observation <- continuous_records_psw(
    .sigma = c(0.1, 0.12, 0.14, 0.11, 0.13, 0.15)
  )
  expect_null(density_meta(per_observation)$sigma_value)
})

test_that("combining weights spread by different numbers drops the density", {
  # Two records that agree on everything but the number the spread was fixed at
  # describe two different ratios, so the combination is described by neither.
  first <- continuous_records_psw(.sigma = 0.12)
  second <- continuous_records_psw(.sigma = 0.2)

  out <- NULL
  expect_warning(
    out <- c(first, second),
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(density_meta(out))

  # The same number twice is the same ratio twice, and the record survives.
  same <- expect_silent(c(first, continuous_records_psw(.sigma = 0.12)))
  expect_identical(density_meta(same)$sigma_value, 0.12)
})

test_that("continuous censoring weights record both", {
  w <- wt_cens(
    records_mu,
    records_exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  )

  expect_identical(estimand(w), "uncensored")
  expect_identical(exposure_type(w), "continuous")
  expect_identical(density_meta_summary(w), marginal_normal_record)
})

# ---- printing the density record --------------------------------------------

test_that("a density record prints the density, numerator, and spread", {
  # The record is three choices, and the print reads them back in the terms the
  # arguments that made them are written in, one to a line and aligned so that
  # the three values read down the same column.
  expect_identical(
    capture.output(print(density_meta(continuous_records_psw()))),
    c(
      "density:   normal",
      "numerator: marginal",
      "sigma:     pooled"
    )
  )

  scored <- continuous_records_psw(stabilization_score = 0.5, .sigma = 0.12)
  expect_identical(
    capture.output(print(density_meta(scored))),
    c(
      "density:   normal",
      "numerator: score",
      "sigma:     supplied"
    )
  )
})

# ---- the records hold at any length -----------------------------------------

test_that("the records survive operations that change the length", {
  w <- continuous_records_psw()
  meta <- density_meta(w)

  shortened <- list(
    `[` = expect_silent(w[1:3]),
    head = expect_silent(utils::head(w, 3)),
    vec_slice = expect_silent(vctrs::vec_slice(w, 1:3)),
    empty = expect_silent(w[integer(0)])
  )

  for (label in names(shortened)) {
    out <- shortened[[label]]
    expect_identical(exposure_type(out), "continuous", label = label)
    expect_identical(density_meta(out), meta, label = label)
  }
})

test_that("the records survive a model frame", {
  # An outcome model shortens its weights column for the rows it drops, in C,
  # and re-attaches the original variable's attributes to the result. Both
  # records describe the exposure rather than the units, so both mean the same
  # thing on the rows the model kept as they did on the rows it was given.
  w <- continuous_records_psw()
  data <- data.frame(
    y = c(1, 2, 3, 4, 5, NA),
    x = c(1, 2, 3, 4, 5, 6),
    wt = w
  )

  fit <- stats::lm(y ~ x, data = data, weights = wt)
  weights <- stats::model.weights(stats::model.frame(fit))

  expect_length(weights, 5)
  expect_identical(exposure_type(weights), "continuous")
  expect_identical(density_meta(weights), density_meta(w))
})

test_that("the records survive arithmetic", {
  w <- continuous_records_psw()
  meta <- density_meta(w)

  # Arithmetic on two sets of weights merges what each of them records, which is
  # the workflow a weight against confounding multiplied by a weight against
  # censoring takes. These two were built by separate calls, so their density
  # records hold separate copies of the function that evaluates the density, and
  # a merge that compared them by identity would report a disagreement between
  # weights built the same way.
  twin <- continuous_records_psw()

  results <- list(
    `w * 2` = expect_silent(w * 2),
    `w + 1` = expect_silent(w + 1),
    `-w` = expect_silent(-w),
    `w / sum(w)` = expect_silent(w / sum(w)),
    `w * twin` = expect_silent(w * twin),
    `w + twin` = expect_silent(w + twin)
  )

  for (label in names(results)) {
    out <- results[[label]]
    expect_s3_class(out, "psw")
    expect_identical(exposure_type(out), "continuous", label = label)
    expect_identical(
      density_meta_summary(out),
      density_meta_summary(w),
      label = label
    )
  }

  # An operand's own record travels through unchanged where there is only one of
  # them to travel.
  expect_identical(density_meta(expect_silent(w * 2)), meta)
})

# ---- combining --------------------------------------------------------------

test_that("combining weights that record the same thing keeps both records", {
  w <- continuous_records_psw()

  out <- expect_silent(c(w, w))
  expect_length(out, 12)
  expect_identical(exposure_type(out), "continuous")
  expect_identical(density_meta_summary(out), marginal_normal_record)

  # Two objects built by separate calls record the same exposure and the same
  # density, so they agree as well as one object does with itself. The density
  # specification carries the function that evaluates it, and two functions
  # written the same way in two calls are different objects, so a merge that
  # compared the records by identity would report a disagreement between weights
  # that were built the same way.
  separate <- expect_silent(c(w, continuous_records_psw()))
  expect_identical(exposure_type(separate), "continuous")
  expect_identical(density_meta_summary(separate), marginal_normal_record)

  combined <- expect_silent(vctrs::vec_c(w, continuous_records_psw()))
  expect_identical(density_meta_summary(combined), marginal_normal_record)

  # Weights for a binary exposure have only the one record to agree on.
  binary <- wt_ate(
    records_binary_ps,
    records_binary_exposure,
    exposure_type = "binary",
    .focal_level = 1
  )
  both_binary <- expect_silent(c(binary, binary))
  expect_identical(exposure_type(both_binary), "binary")
  expect_null(density_meta(both_binary))
})

test_that("combining weights whose densities disagree drops the density", {
  # The two describe the same exposure and were built from the same family, and
  # disagree only on where the residual spread came from, so neither density
  # record describes the combination and only that record is dropped.
  pooled <- continuous_records_psw()
  supplied <- continuous_records_psw(.sigma = 0.12)

  out <- NULL
  cnd <- expect_warning(
    out <- c(pooled, supplied),
    class = "propensity_metadata_conflict_warning"
  )

  expect_s3_class(out, "psw")
  expect_length(out, 12)
  expect_null(density_meta(out))
  expect_identical(exposure_type(out), "continuous")
  expect_true(
    grepl("density_meta", conditionMessage(cnd), fixed = TRUE)
  )
})

test_that("combining a binary and a continuous psw drops both records", {
  binary <- wt_ate(
    records_binary_ps,
    records_binary_exposure,
    exposure_type = "binary",
    .focal_level = 1
  )
  continuous <- continuous_records_psw(stabilize = FALSE)

  out <- NULL
  cnd <- expect_warning(
    out <- c(binary, continuous),
    class = "propensity_metadata_conflict_warning"
  )

  expect_s3_class(out, "psw")
  expect_length(out, 12)
  expect_null(exposure_type(out))

  # The density record describes a continuous exposure. The combination is not
  # described as one, so the record has nothing left to be about, whether or not
  # the other input recorded a density of its own to disagree with it.
  expect_null(density_meta(out))
  expect_true(
    grepl("exposure_type", conditionMessage(cnd), fixed = TRUE)
  )
})

test_that("a dropped record stays dropped however the inputs are ordered", {
  # vctrs settles the common type two inputs at a time, so a third input meets a
  # prototype the earlier pair has already taken the record off. Carrying it
  # back from an input that has nothing left to disagree with would make the
  # result depend on the order the inputs were written in.
  pooled <- continuous_records_psw()
  supplied <- continuous_records_psw(.sigma = 0.12)

  forward <- records_warning_classes(c(pooled, supplied, supplied))
  backward <- records_warning_classes(c(supplied, supplied, pooled))

  # One operation answers for the disagreement once, whichever order the inputs
  # were written in: the pair that meets the record already dropped reports
  # nothing further.
  expect_identical(forward$classes, "propensity_metadata_conflict_warning")
  expect_identical(backward$classes, "propensity_metadata_conflict_warning")

  expect_length(forward$value, 18)
  expect_length(backward$value, 18)
  expect_null(density_meta(forward$value))
  expect_null(density_meta(backward$value))
  expect_identical(exposure_type(forward$value), "continuous")
  expect_identical(exposure_type(backward$value), "continuous")
})

# ---- the records are not part of the type -----------------------------------

test_that("weights that differ only in these records cast to each other", {
  # A cast asks whether one set of weights can be read as another's type. The
  # type is what the weights are: their estimand, their stabilization, and which
  # modification of the propensity scores they were built from. What exposure
  # they describe and which density they came from are recorded alongside that
  # rather than as part of it, so a disagreement about either is not a reason to
  # refuse the conversion.
  pooled <- continuous_records_psw()
  supplied <- continuous_records_psw(.sigma = 0.12)
  binary <- wt_ate(
    records_binary_ps,
    records_binary_exposure,
    exposure_type = "binary",
    .focal_level = 1
  )

  out <- expect_silent(vctrs::vec_cast(pooled, to = supplied))
  expect_s3_class(out, "psw")
  expect_identical(vctrs::vec_data(out), vctrs::vec_data(pooled))

  reversed <- expect_silent(vctrs::vec_cast(supplied, to = pooled))
  expect_s3_class(reversed, "psw")
  expect_identical(vctrs::vec_data(reversed), vctrs::vec_data(supplied))

  # Weights for two different exposures are each other's type as well, as long
  # as they agree about what the weights themselves are. These two are both
  # unstabilized ATE weights built from unmodified propensity scores.
  unstabilized <- continuous_records_psw(stabilize = FALSE)
  crossed <- expect_silent(vctrs::vec_cast(binary, to = unstabilized))
  expect_s3_class(crossed, "psw")
  expect_identical(vctrs::vec_data(crossed), vctrs::vec_data(binary))
})

# ---- weights built from a modified propensity score -------------------------

test_that("weights from a trimmed propensity score carry the records", {
  trimmed <- ps_trim(records_mu, method = "ps", lower = 0.25, upper = 0.8)

  w <- muffle_refit_warning(wt_ate(
    trimmed,
    records_exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  ))

  expect_identical(estimand(w), "ate; trimmed")
  expect_identical(exposure_type(w), "continuous")
  expect_identical(density_meta_summary(w), marginal_normal_record)

  binary <- muffle_refit_warning(wt_ate(
    ps_trim(records_binary_ps, method = "ps", lower = 0.25, upper = 0.65),
    records_binary_exposure,
    exposure_type = "binary",
    .focal_level = 1
  ))

  expect_identical(exposure_type(binary), "binary")
  expect_null(density_meta(binary))
})

test_that("weights from a truncated propensity score carry the records", {
  truncated <- ps_trunc(records_mu, method = "ps", lower = 0.25, upper = 0.8)

  w <- wt_ate(
    truncated,
    records_exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  )

  expect_identical(estimand(w), "ate; truncated")
  expect_identical(exposure_type(w), "continuous")
  expect_identical(density_meta_summary(w), marginal_normal_record)

  binary <- wt_ate(
    ps_trunc(records_binary_ps, method = "ps", lower = 0.25, upper = 0.65),
    records_binary_exposure,
    exposure_type = "binary",
    .focal_level = 1
  )

  expect_identical(exposure_type(binary), "binary")
  expect_null(density_meta(binary))
})

test_that("weights from a calibrated propensity score carry the records", {
  # Calibration maps propensity scores onto the observed rate of a two-level
  # exposure, so the calibrated route is a binary one and has only the exposure
  # type to record.
  exposure <- c(0, 1, 0, 0, 1, 0, 1, 1, 0, 1)
  calibrated <- ps_calibrate(
    c(0.14, 0.22, 0.31, 0.4, 0.48, 0.55, 0.62, 0.7, 0.78, 0.86),
    exposure
  )

  w <- wt_ate(
    calibrated,
    exposure,
    exposure_type = "binary",
    .focal_level = 1
  )

  expect_identical(estimand(w), "ate; calibrated")
  expect_identical(exposure_type(w), "binary")
  expect_null(density_meta(w))
})
