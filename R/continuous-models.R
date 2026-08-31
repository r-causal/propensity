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
# `lm` it inherits from: both are named so that they are reached before the
# classes they inherit from. Other subclasses of `lm` are reached by
# inheritance, which is what a fit that reports one conditional mean per unit
# through `fitted()` needs, and all the `lm` method asks of a model.
extract_continuous_ps <- function(model, call = rlang::caller_env()) {
  UseMethod("extract_continuous_ps")
}

#' @export
extract_continuous_ps.default <- function(model, call = rlang::caller_env()) {
  abort_no_method(model, call = call)
}

#' @export
extract_continuous_ps.lm <- function(model, call = rlang::caller_env()) {
  check_continuous_model_response(model, call = call)

  list(mu = stats::fitted(model))
}

# A model of several responses at once, such as `lm(cbind(a, b) ~ x)`, fits a
# conditional mean for each of them and reports the whole set as a matrix. One
# of those columns may well be the mean of the exposure, but nothing on this
# route says which, so the fit is refused here rather than downstream, where a
# matrix with a row for each unit is only ever described by a length that the
# caller never wrote.
check_continuous_model_response <- function(
  model,
  call = rlang::caller_env()
) {
  fitted_values <- stats::fitted(model)

  if (!inherits(model, "mlm") && !is.matrix(fitted_values)) {
    return(invisible(NULL))
  }

  dims <- dim(as.matrix(fitted_values))

  abort(
    c(
      "Weights for a continuous exposure need a model of one conditional mean
       for each unit.",
      x = "{.arg .propensity} is {.cls {class(model)[[1]]}}, a fit of
           {dims[[2]]} response{?s}, whose fitted values are {dims[[1]]} by
           {dims[[2]]}.",
      i = "Fit the exposure on its own, or pass the conditional means of this
           exposure to {.arg .propensity} as a numeric vector."
    ),
    error_class = "propensity_ps_shape_error",
    call = call
  )
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
#
# `arg` is the argument the model arrived in, and `remedy` the advice that goes
# with it: a propensity score model's fitted means can be passed directly
# instead, and a numerator model a caller stabilizes on is a model or nothing.
check_continuous_model_family <- function(
  model,
  arg = ".propensity",
  remedy = "Fit the propensity score model with {.fun gaussian}, {.fun lm},
            {.fun mgcv::gam}, or {.fun MASS::rlm}, or pass fitted conditional
            means to {.arg .propensity} directly.",
  call = rlang::caller_env()
) {
  family <- model[["family"]]

  no_family <- is.null(family)

  # A family object is a list, and every test below reads a component of it.
  # An element that is not a list is not a family, and `$` on an atomic vector
  # raises an error of base R's rather than one of this package's, so what such
  # an element is gets read before anything is asked of it.
  family_object <- is.list(family)
  gaussian_variance <- family_object &&
    (identical(family$family, "gaussian") ||
      (identical(family$family, "quasi") &&
        identical(family$varfun, "constant")))

  if (no_family || gaussian_variance) {
    return(invisible(NULL))
  }

  fit_with <- if (!family_object) {
    "{.arg {arg}} holds {.obj_type_friendly {family}} where a family
     object goes, so it names no spread at all."
  } else {
    "{.arg {arg}} was fit with {.code {model_family_label(family)}},
     whose spread changes with its fitted values."
  }

  abort(
    c(
      "Weights for a continuous exposure need a model of its conditional mean
       with a single spread.",
      x = fit_with,
      i = remedy
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

  # A family object is a list. An element that is not one is not a family, and
  # `$` on an atomic vector raises an error of base R's rather than one of this
  # package's, so what the element is gets read before anything is asked of it.
  #
  # `isTRUE()` because a family object that carries no name answers the test
  # with nothing at all, which `&&` can read as neither true nor false. A
  # family like that is not one of the two, and is refused below by what it is
  # missing.
  family_object <- is.list(family)
  binomial_family <- family_object &&
    isTRUE(family$family %in% c("binomial", "quasibinomial"))

  if (binomial_family) {
    return(invisible(NULL))
  }

  fitted_by <- if (is.null(family)) {
    "{.arg .propensity} is {.cls {class(model)[[1]]}}, whose fitted values are
     conditional means rather than probabilities."
  } else if (!family_object) {
    "{.arg .propensity} holds {.obj_type_friendly {family}} where a family
     object goes, so it names nothing that fits a probability."
  } else if (!family_is_named(family)) {
    "{.arg .propensity} was fit with an unnamed family, whose fitted values are
     conditional means rather than probabilities."
  } else {
    "{.arg .propensity} was fit with {.code {model_family_label(family)}}, whose
     fitted values are conditional means rather than probabilities."
  }

  abort(
    c(
      "A binary propensity score needs a model of the probability of the
       exposure.",
      x = fitted_by,
      i = "Fit the propensity score model with {.fun binomial}, or pass fitted
           probabilities to {.arg .propensity} directly."
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

# Whether a family object names itself. Every family built by `binomial()` and
# its neighbors holds its name in `$family`, and a message about one is written
# from that string; an object that carries a link and no name has nothing to be
# written from, and is described by what it is rather than by what it is called.
# A name of no characters is no name either: written into a message it leaves a
# pair of backticks with nothing between them.
family_is_named <- function(family) {
  rlang::is_string(family$family) && nzchar(family$family)
}

# How a family reads in a message. A `quasi()` family is named by its variance
# function as well, that being what decides whether its fitted values have one
# spread or many. An extended family names itself with the parameters it was fit
# at, such as mgcv's `"Scaled t(Inf,1.012)"`, so a name is written as a call only
# when it is a bare one: appending parentheses to a name that already reads as a
# call writes one that could not be run.
model_family_label <- function(family) {
  # A family with no name has none to write. A refusal that has a fuller
  # sentence for that case tests `family_is_named()` itself and never arrives
  # here; this keeps the label readable for any that does not.
  if (!family_is_named(family)) {
    return("an unnamed family")
  }

  if (identical(family$family, "quasi")) {
    return(paste0("quasi(variance = \"", family$varfun, "\")"))
  }

  if (grepl("^[A-Za-z.][A-Za-z0-9._]*$", family$family)) {
    return(paste0(family$family, "()"))
  }

  family$family
}
