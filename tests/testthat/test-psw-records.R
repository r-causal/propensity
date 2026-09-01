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
#
# The fitted means track the exposure closely, so the weights sit past the
# finite-variance boundary and say so. These tests are written about what the
# weights record rather than about how well behaved they are, so the report is
# muffled here.
continuous_records_psw <- function(
  stabilize = TRUE,
  stabilization_score = NULL,
  .sigma = NULL,
  .density = "normal"
) {
  muffle_variance_warning(wt_ate(
    records_mu,
    records_exposure,
    .sigma = .sigma,
    exposure_type = "continuous",
    stabilize = stabilize,
    stabilization_score = stabilization_score,
    .density = .density
  ))
}

# A numerator model for the same exposure: the exposure on a baseline covariate
# the conditional means above do not read. `shift` moves the covariate, which
# moves the coefficients the model is fit at without changing the terms it
# reads.
records_numerator_model <- function(shift = 0) {
  exposure <- records_exposure
  covariate <- c(0.1, 0.5, 0.2, 0.9, 0.4, 0.7) + shift

  stats::lm(exposure ~ covariate)
}

# A numerator model whose formula is too long for `deparse()` to write on one
# line. The covariate names are what make it wrap, so the record has several
# deparsed lines to read back as one.
records_long_numerator_model <- function() {
  frame <- data.frame(
    exposure = records_exposure,
    baseline_covariate_measured_at_enrollment = c(
      0.1,
      0.5,
      0.2,
      0.9,
      0.4,
      0.7
    ),
    baseline_covariate_measured_at_randomization = c(
      0.3,
      0.2,
      0.8,
      0.1,
      0.6,
      0.5
    ),
    baseline_covariate_measured_at_the_first_visit = c(
      0.4,
      0.7,
      0.1,
      0.3,
      0.9,
      0.2
    )
  )

  stats::lm(
    exposure ~
      baseline_covariate_measured_at_enrollment +
      baseline_covariate_measured_at_randomization +
      baseline_covariate_measured_at_the_first_visit,
    data = frame
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
    c("density", "numerator", "sigma", "sigma_value", "numerator_model")
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

test_that("combining weights stabilized on different models drops the density", {
  # A numerator model is compared by what it says rather than by identity, the
  # way the density specification is: a fit carries the frame it was made in
  # and the call that made it, so the same numerator fit in two calls would
  # otherwise read as a disagreement.
  first <- continuous_records_psw(stabilize = records_numerator_model())
  again <- continuous_records_psw(stabilize = records_numerator_model())

  same <- expect_silent(c(first, again))
  expect_identical(
    stats::coef(density_meta(same)$numerator_model),
    stats::coef(records_numerator_model())
  )

  # Two models that read the same terms and were fit to different values of
  # them estimate two different numerators, so the combination is stabilized on
  # neither.
  other <- continuous_records_psw(
    stabilize = records_numerator_model(shift = 0.4)
  )

  out <- NULL
  expect_warning(
    out <- c(first, other),
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(density_meta(out))

  # So do a model and the marginal density it was supplied instead of.
  marginal <- NULL
  expect_warning(
    marginal <- c(first, continuous_records_psw()),
    class = "propensity_metadata_conflict_warning"
  )
  expect_null(density_meta(marginal))
})

test_that("continuous censoring weights record both", {
  w <- muffle_variance_warning(wt_cens(
    records_mu,
    records_exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  ))

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

  # A numerator of "model" says that a model the caller fit estimated the
  # numerator, but not which one, so the record reads the model's formula back
  # under the three choices, in the terms of the argument it arrived in.
  modeled <- continuous_records_psw(stabilize = records_numerator_model())
  expect_identical(
    capture.output(print(density_meta(modeled))),
    c(
      "density:   normal",
      "numerator: model",
      "sigma:     pooled",
      "stabilize: exposure ~ covariate"
    )
  )
})

test_that("a wrapped numerator formula reads back as one line of formula spacing", {
  # `deparse()` breaks a long formula across lines and indents the
  # continuation, so collapsing those lines leaves runs of spaces inside the
  # formula. The record reads it back spaced the way a formula is written.
  modeled <- continuous_records_psw(stabilize = records_long_numerator_model())
  printed <- capture.output(print(density_meta(modeled)))

  line <- printed[startsWith(printed, "stabilize:")]
  expect_length(line, 1L)

  read_back <- sub("^stabilize: +", "", line)
  expect_false(grepl("  ", read_back, fixed = TRUE))
  expect_identical(
    read_back,
    paste(
      "exposure ~ baseline_covariate_measured_at_enrollment +",
      "baseline_covariate_measured_at_randomization +",
      "baseline_covariate_measured_at_the_first_visit"
    )
  )
})

# ---- the density record under the printed weights ---------------------------

# What printing a set of weights writes, at a width no line of it wraps at.
psw_printed <- function(x) {
  withr::local_options(width = 200)

  utils::capture.output(print(x))
}

test_that("printing continuous weights writes the density record under them", {
  w <- continuous_records_psw()
  printed <- psw_printed(w)

  # The footer is the record's own `format()` method rather than a second
  # rendering of the same three lines, so the printed weights and the printed
  # record cannot come to disagree.
  expect_identical(utils::tail(printed, 3L), format(density_meta(w)))

  expect_match(printed, "^density:\\s+normal$", all = FALSE)
  expect_match(printed, "^numerator:\\s+marginal$", all = FALSE)
  expect_match(printed, "^sigma:\\s+pooled$", all = FALSE)

  # The header and the weights themselves are what they always were: the record
  # is written under them and nothing else moves.
  expect_length(printed, 5L)
  expect_identical(printed[[1]], "<psw{estimand = ate; stabilized}[6]>")

  # Weights stabilized on a model carry a fourth line naming it, and it reaches
  # the printed weights the same way: through the record's own `format()`.
  modeled <- continuous_records_psw(stabilize = records_numerator_model())
  modeled_printed <- psw_printed(modeled)

  expect_identical(
    utils::tail(modeled_printed, 4L),
    format(density_meta(modeled))
  )
  expect_match(modeled_printed, "^numerator:\\s+model$", all = FALSE)
  expect_match(
    modeled_printed,
    "^stabilize:\\s+exposure ~ covariate$",
    all = FALSE
  )
  expect_length(modeled_printed, 6L)
})

test_that("the footer tells apart weights whose values print alike", {
  # Two sets of weights built from different families carry the same estimand
  # and, to a reader of the numbers, much the same weights. What each is a ratio
  # of is the difference between them, and it is the footer that says so.
  normal <- psw_printed(continuous_records_psw())
  heavy <- psw_printed(continuous_records_psw(.density = dens_t(4)))

  expect_identical(normal[[1]], heavy[[1]])
  expect_match(heavy, "^density:\\s+t\\(df = 4\\)$", all = FALSE)
  expect_false(any(grepl("^density:\\s+t", normal)))
})

test_that("binary and categorical weights print with no footer", {
  weights <- list(
    binary = wt_ate(records_binary_ps, records_binary_exposure),
    categorical = wt_ate(
      records_categorical_ps,
      records_categorical_exposure
    )
  )

  for (name in names(weights)) {
    printed <- psw_printed(weights[[name]])

    # Weights that are not a ratio of densities record no ratio and print none:
    # the header and the values, exactly as they always have.
    expect_identical(length(printed), 2L, label = name)
    expect_false(
      any(grepl("^(density|numerator|sigma):", printed)),
      label = name
    )
  }
})

test_that("the footer follows the record through a shorter slice", {
  w <- continuous_records_psw()
  footer <- format(density_meta(w))

  # The record describes the exposure the weights were built for rather than the
  # units, so it survives a slice, and what is printed under the weights follows
  # it.
  expect_identical(utils::tail(psw_printed(w[1:3]), 3L), footer)
  expect_identical(utils::tail(psw_printed(utils::head(w, 2L)), 3L), footer)
})

test_that("weights that lost the record print no footer", {
  first <- continuous_records_psw(.sigma = 0.12)
  second <- continuous_records_psw(.sigma = 0.2)

  out <- NULL
  expect_warning(
    out <- c(first, second),
    class = "propensity_metadata_conflict_warning"
  )

  expect_null(density_meta(out))
  expect_false(any(grepl("^(density|numerator|sigma):", psw_printed(out))))
})

test_that("the printed weights read as they are pinned", {
  # The whole of what printing writes, for weights that carry a record and for
  # weights that do not: the footer is user-visible output, and the header and
  # the values under it are what they were before there was one.
  expect_snapshot({
    print(continuous_records_psw())

    print(continuous_records_psw(.density = dens_t(4)))

    print(continuous_records_psw(stabilize = records_numerator_model()))

    print(wt_ate(records_binary_ps, records_binary_exposure))
  })
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

test_that("arithmetic that leaves the result unstabilized drops the numerator", {
  # A numerator is what the weights were divided into, and a product only half
  # of which was divided by it was not. The merge already drops the
  # stabilization status and the score there, and the numerator fields of the
  # density record say the same thing a third time: left standing, the record
  # reads back as a ratio the product is not, and prints a `stabilize:` line
  # naming a model that never multiplied it.
  #
  # The other operand carries no record of its own, so the modeled record
  # reaches the result by the single-operand route rather than by agreeing with
  # anything. Two records that disagree about the numerator are a conflict, and
  # are dropped whole with a report of their own.
  modeled <- continuous_records_psw(stabilize = records_numerator_model())
  plain <- psw(rep(1.5, length(records_exposure)), estimand = "ate")

  results <- list(
    `modeled * plain` = expect_silent(modeled * plain),
    `plain * modeled` = expect_silent(plain * modeled),
    `modeled + plain` = expect_silent(modeled + plain)
  )

  for (label in names(results)) {
    out <- results[[label]]
    expect_false(is_stabilized(out), label = label)
    expect_null(numerator_model(out), label = label)
    expect_null(density_meta(out)$numerator_model, label = label)
    expect_identical(density_meta(out)$numerator, "none", label = label)

    # The denominator is the conditional density the weights divide by, which
    # the arithmetic leaves alone, so the family and the spread stay.
    expect_identical(
      format(density_meta(out)$density),
      format(density_meta(modeled)$density),
      label = label
    )
    expect_identical(density_meta(out)$sigma, "pooled", label = label)
    expect_false(
      any(grepl("stabilize:", format(density_meta(out)), fixed = TRUE)),
      label = label
    )
  }
})

test_that("arithmetic that stays stabilized keeps the numerator model", {
  # The drop above is the unstabilized result's, not every product's: two sets
  # of weights stabilized on the same model were both divided by it, and so was
  # their product.
  modeled <- continuous_records_psw(stabilize = records_numerator_model())
  twin <- continuous_records_psw(stabilize = records_numerator_model())

  out <- expect_silent(modeled * twin)
  expect_true(is_stabilized(out))
  expect_identical(density_meta(out)$numerator, "model")
  expect_identical(
    stats::coef(numerator_model(out)),
    stats::coef(records_numerator_model())
  )
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

  w <- muffle_variance_warning(muffle_refit_warning(wt_ate(
    trimmed,
    records_exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  )))

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

  w <- muffle_variance_warning(wt_ate(
    truncated,
    records_exposure,
    exposure_type = "continuous",
    stabilize = TRUE
  ))

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

test_that("weights from a calibrated score that was trimmed carry the calibration", {
  # Trimming and truncation read the scores a calibrated object holds and
  # record the calibration in their own metadata, so the weights built from one
  # were built from calibrated scores and say so, exactly as the weights built
  # from the calibrated score itself do.
  exposure <- c(0, 1, 0, 0, 1, 0, 1, 1, 0, 1)
  calibrated <- ps_calibrate(
    c(0.14, 0.22, 0.31, 0.4, 0.48, 0.55, 0.62, 0.7, 0.78, 0.86),
    exposure
  )

  trimmed <- muffle_refit_warning(wt_ate(
    ps_trim(calibrated, method = "ps", lower = 0.05, upper = 0.95),
    exposure,
    exposure_type = "binary",
    .focal_level = 1
  ))

  expect_identical(estimand(trimmed), "ate; trimmed; calibrated")
  expect_true(is_ps_calibrated(trimmed))
  expect_true(is_ps_trimmed(trimmed))

  # The record is the score's rather than the estimand's, so it reaches the
  # weights whichever formula built them.
  truncated <- wt_att(
    ps_trunc(calibrated, method = "ps", lower = 0.05, upper = 0.95),
    exposure,
    exposure_type = "binary",
    .focal_level = 1
  )

  expect_identical(estimand(truncated), "att; truncated; calibrated")
  expect_true(is_ps_calibrated(truncated))
  expect_true(is_ps_truncated(truncated))

  # A score that was never calibrated records nothing that says it was.
  plain <- wt_att(
    ps_trunc(records_binary_ps, method = "ps", lower = 0.05, upper = 0.95),
    records_binary_exposure,
    exposure_type = "binary",
    .focal_level = 1
  )

  expect_identical(estimand(plain), "att; truncated")
  expect_false(is_ps_calibrated(plain))
})
