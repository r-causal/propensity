# The fitted models a continuous exposure's weights can be read from.
#
# The weights are a ratio of densities, so what the model has to supply is the
# center of the conditional density: one fitted conditional mean for each unit,
# on the scale of the exposure, around which a single spread describes the
# residuals. Every class here reports that mean through `fitted()`, whatever
# link it was fit under, and every one of them is spread by the pooled residual
# root mean square unless the caller supplies `.sigma`.
#
# Methods are resolved in class order, so a `gam` is read as a `gam` rather than
# as the `glm` it inherits from, and an `rlm` as an `rlm` rather than as the
# `lm` it inherits from. Each supported class is named here rather than reached
# by inheritance, so that the set of models the weights can be built from is
# the set of methods.
extract_continuous_ps <- function(model, call = rlang::caller_env()) {
  UseMethod("extract_continuous_ps")
}

#' @export
extract_continuous_ps.default <- function(model, call = rlang::caller_env()) {
  abort_no_method(model, call = call)
}

#' @export
extract_continuous_ps.lm <- function(model, call = rlang::caller_env()) {
  list(mu = stats::fitted(model))
}

# `MASS::rlm()` fits the same conditional mean by a weighted least squares
# iteration, so its fitted values are read the way an `lm`'s are. Its spread is
# not: `rlm` reports a robust scale estimate in `fit$s`, which resists the
# extreme residuals rather than pooling all of them, while the density these
# weights are a ratio of is spread by the pooled residual root mean square, as
# it is for every other class. That is deliberate, and the robust scale is
# available by passing `.sigma = fit$s`.
#' @export
extract_continuous_ps.rlm <- function(model, call = rlang::caller_env()) {
  extract_continuous_ps.lm(model, call = call)
}

#' @export
extract_continuous_ps.glm <- function(model, call = rlang::caller_env()) {
  check_continuous_model_family(model, call = call)

  # `fitted()` is on the scale of the response whatever the link is, so a log or
  # inverse link never has to be undone here.
  list(mu = stats::fitted(model))
}

# An `mgcv::gam()` is a `glm` whose conditional mean is a sum of smooth terms,
# and reports that mean the same way, so it is held to the same family check.
#' @export
extract_continuous_ps.gam <- function(model, call = rlang::caller_env()) {
  extract_continuous_ps.glm(model, call = call)
}

# The families whose fitted values are a conditional mean with one spread for
# every observation. A model fit by least squares carries no family and is one
# of them. `quasi()` is one when its variance is constant, which is the
# gaussian variance under another name; any other variance function describes a
# spread that changes with the fitted value, which the conditional density has
# no way to take.
check_continuous_model_family <- function(model, call = rlang::caller_env()) {
  family <- model[["family"]]

  no_family <- is.null(family)
  gaussian_variance <- identical(family$family, "gaussian") ||
    (identical(family$family, "quasi") && identical(family$varfun, "constant"))

  if (no_family || gaussian_variance) {
    return(invisible(NULL))
  }

  abort(
    c(
      "Weights for a continuous exposure need a model of its conditional mean
       with a single spread.",
      x = "{.arg .propensity} was fit with {.code {model_family_label(family)}},
           whose spread changes with its fitted values.",
      i = "Fit the propensity score model with {.fun gaussian}, {.fun lm},
           {.fun mgcv::gam}, or {.fun MASS::rlm}, or pass fitted conditional
           means to {.arg .propensity} directly."
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

# Whether a fitted model reports the probability of a binary exposure. Most
# classes answer that from their family, which is what the default method
# reads; a class that carries no family and still fits a probability, such as
# `nnet::multinom()`, answers it from what it was fit to and has a method of its
# own.
check_binary_model_family <- function(model, call = rlang::caller_env()) {
  UseMethod("check_binary_model_family")
}

# The families whose fitted values are the probability of the exposure. A model
# of a conditional mean is not one of them, however close to a probability its
# fitted values happen to fall: a linear probability model is not held to the
# unit interval, so this reads the model rather than the values it fitted, and
# reads it before the range of a propensity score is checked.
#' @export
check_binary_model_family.default <- function(
  model,
  call = rlang::caller_env()
) {
  family <- model[["family"]]

  if (!is.null(family) && family$family %in% c("binomial", "quasibinomial")) {
    return(invisible(NULL))
  }

  fitted_by <- if (is.null(family)) {
    "{.arg .propensity} is {.cls {class(model)[[1]]}}, whose fitted values are
     conditional means rather than probabilities."
  } else {
    "{.arg .propensity} was fit with {.code {model_family_label(family)}}, whose
     fitted values are conditional means rather than probabilities."
  }

  abort(
    c(
      "Weights for a binary exposure need a model of the probability of that
       exposure.",
      x = fitted_by,
      i = "Fit the propensity score model with {.fun binomial}, or pass fitted
           probabilities to {.arg .propensity} directly."
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

# How a family reads in a message. A `quasi()` family is named by its variance
# function as well, that being what decides whether its fitted values have one
# spread or many.
model_family_label <- function(family) {
  if (identical(family$family, "quasi")) {
    return(paste0("quasi(variance = \"", family$varfun, "\")"))
  }

  paste0(family$family, "()")
}
