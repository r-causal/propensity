# The fitted models a categorical exposure's weights can be read from.
#
# The weights are read off a probability for every level, so what the model has
# to supply is one fitted column per level of the exposure, named after the
# level it belongs to. `nnet::multinom()` is the fit that reports that, and it
# is the class named here: methods are resolved in class order, so a `multinom`
# is read as a `multinom` rather than as the `nnet` it inherits from, and the
# set of models a categorical exposure's weights can be built from is the set
# of methods.
#
# A `multinom` of two levels is not a categorical fit. It reports a single
# probability, the shape a binomial `glm` reports, and a two-level exposure
# resolves as binary, so such a fit is read on the binary path instead. The
# methods below hold both halves of that: what the categorical path reads, and
# what the binary path accepts and reads.
extract_categorical_ps <- function(
  model,
  exposure_levels = NULL,
  call = rlang::caller_env()
) {
  UseMethod("extract_categorical_ps")
}

#' @export
extract_categorical_ps.default <- function(
  model,
  exposure_levels = NULL,
  call = rlang::caller_env()
) {
  abort_categorical_model(model, call = call)
}

# `fitted()` is the n by K matrix of level probabilities, carrying the levels as
# its column names in the order the response was coded in. The order the
# exposure was supplied in is the order the weights are read in, and
# `check_ps_matrix()` matches the columns to it by name, so nothing is reordered
# here.
#
# The response check that every model method runs first has already refused a
# fit with no levels, so the columns are named whenever this is reached.
#' @export
extract_categorical_ps.multinom <- function(
  model,
  exposure_levels = NULL,
  call = rlang::caller_env()
) {
  check_multinom_levels(model, exposure_levels, call = call)

  stats::fitted(model)
}

# The levels a `multinom` was fit to are the columns it reports, so a fit made
# to some other set of levels than those of the exposure being weighted has no
# column to read as the probability of each level. Comparing the two sets here
# reports that as what it is, a model of something other than the exposure,
# rather than leaving `check_ps_matrix()` to report the width of a matrix the
# caller never built. Order is not part of the comparison, because the matrix
# check matches the columns to the levels by name, as does the gather a
# numerator model's fitted probabilities are read with.
#
# An exposure of fewer than three levels is not a categorical exposure at all,
# whatever the model was fit to, and `transform_exposure_categorical()` reports
# that; the comparison leaves such an exposure to it.
#
# `arg` and `remedy` name the argument the fit arrived in and what to do about
# it, because the same question is asked of two models: the propensity score
# model the weights are read off, and the numerator model they are stabilized
# on. Both are read by level name, so both are held to the same set.
check_multinom_levels <- function(
  model,
  exposure_levels,
  arg = ".propensity",
  remedy = "Fit the propensity score model to the exposure being weighted.",
  call = rlang::caller_env()
) {
  if (length(exposure_levels) < 3L) {
    return(invisible(NULL))
  }

  if (setequal(as.character(exposure_levels), model$lev)) {
    return(invisible(NULL))
  }

  n_levels <- length(model$lev)
  n_exposure <- length(exposure_levels)

  # The mismatch an unused factor level makes: `nnet::multinom()` drops a level
  # no unit took, so a fit made to an exposure with an empty level reports
  # every level the exposure has units at and no other. The fit is of the
  # exposure being weighted and only the declaration disagrees, which is fixed
  # without refitting anything, so it is worth naming. A fit reporting a level
  # the exposure does not declare is not that, whatever else it reports: it is
  # a fit of something else, and dropping levels cannot supply the column it is
  # missing.
  exposure_levels_chr <- as.character(exposure_levels)
  unfitted <- setdiff(exposure_levels_chr, model$lev)
  unused_level_remedy <- if (
    length(unfitted) > 0 && all(model$lev %in% exposure_levels_chr)
  ) {
    c(
      i = "{.fun nnet::multinom} drops a level no unit took, so an exposure
           declaring {.val {unfitted}} with no unit there is fit one column
           short; call {.fun droplevels} on {.arg .exposure} if that is the
           disagreement."
    )
  }

  abort(
    c(
      "Weights for a categorical exposure need a probability for every level
       of the exposure being weighted.",
      x = "{.arg {arg}} was fit to {n_levels} level{?s}
           ({.val {model$lev}}), and the exposure has {n_exposure}
           ({.val {exposure_levels}}).",
      unused_level_remedy,
      i = remedy
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

# Whether a fitted model reports a probability for every level of the exposure
# rather than a single probability for each unit. That is what decides which
# route a fit takes through the propensity score modifiers, whose vector and
# matrix methods read the two shapes. The weight functions settle the same
# question from the exposure they were given, which a modifier such as
# `ps_tilt()` does not have, so the answer is read from the fit itself.
model_fits_levels <- function(model) {
  UseMethod("model_fits_levels")
}

#' @export
model_fits_levels.default <- function(model) {
  FALSE
}

# A two-level fit reports one column, the shape a binomial `glm` reports, and a
# two-level exposure is binary however it was fit, so such a fit is read on the
# vector route. `check_multinom_response()` has already refused a fit with no
# levels wherever this is reached.
#' @export
model_fits_levels.multinom <- function(model) {
  length(model$lev) > 2L
}

# The levels a fitted model reports a probability for, read off the fit itself.
# A modification that reads no exposure still has to know which levels the
# columns it is given belong to, and the fit records them, so nothing has to be
# recovered from the data the model was fit to. This is the pair of
# `model_fits_levels()`: a class that answers `TRUE` there reports a column for
# every level and so has levels to name here.
model_levels <- function(model) {
  UseMethod("model_levels")
}

# `lev` is the response's levels in the order the columns of `fitted()` are laid
# out in, which is the order the columns are named in.
#' @export
model_levels.multinom <- function(model) {
  model$lev
}

# A `multinom` fits a probability for every level and so has a single
# probability to give only when it was fit to two of them. More levels than
# that leave nothing to read as the probability of the exposure, which is what a
# binary exposure needs.
#
# `remedy` goes unread here. What is wrong with a `multinom` of more than two
# levels is the number of levels rather than the family, so the remedy is the
# one written below whichever argument the model arrived in. A `multinom`
# supplied as a numerator model never reaches this method at all:
# `check_numerator_model()` refuses one for a binary exposure before the family
# is read, since a two-level fit would pass here and the binary numerator is
# then read off a response and a family a `multinom` does not carry.
#' @export
check_binary_model_family.multinom <- function(
  model,
  arg = ".propensity",
  remedy = NULL,
  call = rlang::caller_env()
) {
  n_levels <- length(model$lev)

  if (n_levels == 2L) {
    return(invisible(NULL))
  }

  abort(
    c(
      "A binary propensity score needs the probability of one of the
       exposure's two levels.",
      x = "{.arg {arg}} was fit to {n_levels} levels
           ({.val {model$lev}}), so it fits no single probability to read
           against a binary exposure.",
      i = "Fit the propensity score model to the exposure being weighted, or
           weight the exposure the model was fit to with
           {.code exposure_type = \"categorical\"}."
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

# A two-level fit reports one column, the probability of the second level, which
# is the level the default coding treats as focal. Flattening it to a bare
# vector gives the binary formulas the shape they read and is what
# `prepare_model_weight_args()` inverts when the caller names the other level.
#' @export
extract_binary_ps.multinom <- function(model, call = rlang::caller_env()) {
  as.numeric(stats::fitted(model))
}

# A `multinom` keeps no model frame, so `model.frame()` rebuilds one by
# re-evaluating the fitting call in the environment the formula came from. That
# environment no longer holds the arguments of a fit made inside a function, and
# the rebuild fails there, which would leave the exposure unreadable for a model
# fit by an ordinary wrapper.
#
# The fit records the response itself: `residuals` is the indicator matrix of
# the observed level less the fitted probabilities, laid out in the order of
# `lev`, so adding the fitted values back recovers the indicators exactly. A
# two-level fit indicates the second level in its single column; more levels get
# a column each, and the level indicated is the one whose column holds the one.
#' @export
extract_model_exposure.multinom <- function(model) {
  indicators <- model$residuals + model$fitted.values

  indicated <- if (ncol(indicators) == 1L) {
    round(indicators[, 1L]) + 1L
  } else {
    max.col(indicators)
  }

  pad_dropped_rows(model, factor(model$lev[indicated], levels = model$lev))
}

# What a `multinom` was fit to. A matrix response is read as counts rather than
# as a factor, so the fit records no levels, and its fitted values name no level
# that could be matched to an exposure. The model methods run this before
# anything reads the model, including before the exposure is taken off it, so
# that a fit the route cannot read is reported as such rather than as whatever
# the exposure detector makes of a matrix of counts.
#
# `arg` names the argument the fit arrived as, since the estimator takes it as
# `wt_mod` and the weight functions as `.propensity`, and the refusal is the
# same one in both places. `what` names which side of the ratio the fit was
# supplied for, a numerator model being held to the same requirement as the
# model the scores themselves are read off.
check_multinom_response <- function(
  model,
  arg = ".propensity",
  what = "A propensity score model",
  call = rlang::caller_env()
) {
  if (length(model$lev) > 0) {
    return(invisible(NULL))
  }

  abort(
    c(
      "{what} must be fit to the levels of the exposure.",
      x = "{.arg {arg}} was fit to a matrix response, which
           {.fun nnet::multinom} reads as counts rather than as levels.",
      i = "Refit the model with the exposure factor on the left-hand side, as
           in {.code nnet::multinom(exposure ~ x)}."
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

# Whether a fitted model reports a probability for every level of a categorical
# exposure. This is the categorical counterpart of `check_binary_model_family()`
# and is asked of the numerator model a categorical exposure's weights are
# stabilized on; the propensity score model of one is read by
# `extract_categorical_ps()`, whose methods name the same set of classes. Which
# classes answer is which classes have a method, so a model fitting one value
# per unit is refused by the default whatever family it was fit with.
#
# `arg` names the argument the model arrived in and `remedy` says what to do
# about it, both written the way the binary check takes them: the two sides of
# the ratio arrive in different arguments and take different remedies.
check_categorical_model_family <- function(
  model,
  arg = ".propensity",
  remedy = "Pass a fitted {.fun nnet::multinom} of the exposure, or a matrix or
            data frame with one column per level.",
  call = rlang::caller_env()
) {
  UseMethod("check_categorical_model_family")
}

#' @export
check_categorical_model_family.default <- function(
  model,
  arg = ".propensity",
  remedy = "Pass a fitted {.fun nnet::multinom} of the exposure, or a matrix or
            data frame with one column per level.",
  call = rlang::caller_env()
) {
  abort(
    c(
      "Weights for a categorical exposure need a probability for every level.",
      x = "{.arg {arg}} is {.cls {class(model)[[1]]}}, which fits one value for
           each unit rather than one for each level.",
      i = remedy
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

# A `multinom` reports a column per level, which is the shape both sides of the
# ratio are read in. A matrix response is the one fit that does not: it is read
# as counts rather than as a factor, so the fit records no levels and its
# columns name none for the levels of the exposure to be matched against.
#
# The number of levels is no part of the question here. Two of them make a
# binary exposure however the fit was made, and the exposure's own levels decide
# which route the weights take, so a fit reporting too few columns for the
# exposure is reported against those levels rather than against its family.
#' @export
check_categorical_model_family.multinom <- function(
  model,
  arg = ".propensity",
  remedy = "Pass a fitted {.fun nnet::multinom} of the exposure, or a matrix or
            data frame with one column per level.",
  call = rlang::caller_env()
) {
  check_multinom_response(
    model,
    arg = arg,
    what = categorical_model_role(arg),
    call = call
  )
}

# Which of the two models a refusal is about, read from the argument the fit
# arrived in: the propensity score model of a categorical exposure is supplied
# as `.propensity` and its numerator model as `stabilize`, so the argument name
# is what tells the two apart wherever the same refusal serves both.
categorical_model_role <- function(arg) {
  if (identical(arg, "stabilize")) {
    "A numerator model"
  } else {
    "A propensity score model"
  }
}

# A categorical exposure needs a probability for every level, which a model
# reporting one fitted value per unit cannot give. The classes that fit one
# value per unit reach this through the default method.
abort_categorical_model <- function(model, call = rlang::caller_env()) {
  abort(
    c(
      "Weights for a categorical exposure need a probability for every level.",
      x = "{.arg .propensity} is {.cls {class(model)[[1]]}}, which fits one
           value for each unit rather than one for each level.",
      i = "Pass a fitted {.fun nnet::multinom}, or a matrix or data frame with
           one column per level, such as that model's fitted values."
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}
