# Product weights for a joint intervention on two treatments, and the container
# recording the two treatment models the weights were built from.
#
# A joint intervention on A and E needs a weight for the joint exposure, which
# factorizes sequentially as f(A | L) f(E | A, L). Multiplying the two component
# weight vectors is one line; what these constructors are for is the set of
# things that make the product not a joint weight. Every refusal below is one of
# them, and the factorization guardrail is the load-bearing one: the product of
# two marginal weights is an ordinary vector of positive numbers that nothing
# downstream can tell from the real thing.

# The exposure types a component may declare. A joint weight is a product over
# two treatments, and each of them is weighted the way its own type is.
joint_wt_exposure_types <- c("binary", "categorical", "continuous")

# The model classes a treatment model may be, and the exposure type each one
# implies, in the vocabulary the rest of the package uses for exposures.
#
# Each class is named rather than reached by inheritance, because the two the
# package reads a dose from are subclasses of the two it reads other things
# from: an `rlm` is an `lm` and a `gam` is a `glm`. Written as inheritance the
# branch for the parent would type them, and it would type an unfamiliar
# subclass of either parent along with them, which is how a model this package
# cannot read a density from would pass silently.
joint_wt_model_type <- function(model) {
  classes <- class(model)

  if (inherits(model, "multinom")) {
    return("categorical")
  }

  if (inherits(model, "gam")) {
    # An additive model of a dose is a gaussian one. Any other family fits
    # something this package reads no treatment density from.
    if (identical(model$family$family, "gaussian")) {
      return("continuous")
    }
    return(NULL)
  }

  if (inherits(model, "rlm")) {
    return("continuous")
  }

  if (inherits(model, "glm")) {
    family <- model$family$family
    if (family %in% c("binomial", "quasibinomial")) {
      return("binary")
    }
    if (family %in% c("gaussian", "quasi")) {
      return("continuous")
    }
    return(NULL)
  }

  if (identical(classes, "lm")) {
    return("continuous")
  }

  NULL
}

# The variables a model's right-hand side reads. Membership rather than term
# labels: a treatment written as `factor(a) * x1` produces the labels
# `factor(a)`, `x1`, and `factor(a):x1`, and one written as `x1 + a:x1` produces
# `x1` and `x1:a`. Neither produces the label `a`, and both condition on `a`, so
# a check written against term labels would refuse two models that factorize
# correctly.
joint_wt_reads <- function(model) {
  all.vars(stats::formula(model)[[3L]])
}

# The variables a model's response reads, which is how a model is matched
# against the treatment it is recorded under.
joint_wt_response <- function(model) {
  all.vars(stats::formula(model)[[2L]])
}

# ---- wt_joint() -------------------------------------------------------------

#' Product weights for a joint intervention on two treatments
#'
#' @description
#' `wt_joint()` builds the weight for a joint intervention on two treatments by
#' multiplying the two component weight vectors. A joint exposure needs a weight
#' for the pair, and that weight factorizes sequentially: the density of the
#' first treatment given the covariates, times the density of the second given
#' the first treatment and the covariates. `wt_joint()` is the product, plus the
#' checks that make the product a joint weight rather than an arbitrary one.
#'
#' `is_joint_wt()` reports whether a set of weights was built this way, and
#' `joint_wt_meta()` reads back what each component was.
#'
#' @details
#' # The factorization
#'
#' The weight for a joint intervention on `A` and `E` is
#' \eqn{1 / [f(A | L) f(E | A, L)]}. The second factor conditions on the first
#' treatment. The product of two *marginal* weights,
#' \eqn{1 / [f(A | L) f(E | L)]}, is a different quantity, and it is not the
#' joint weight for any data in which `E` depends on `A`.
#'
#' Nothing downstream can tell the two apart. Both are ordinary vectors of
#' positive numbers, both fit an outcome model without complaint, and both
#' return estimates with standard errors. That is why [joint_wt_models()]
#' refuses a second model that does not condition on the first treatment, and
#' why the dependence should be modeled flexibly rather than as a single
#' additive term: an additive term satisfies the factorization but may model it
#' badly, and a badly modeled dependence biases the joint weight in a way that
#' is equally invisible.
#'
#' # Stabilization
#'
#' A continuous component must be stabilized. The unstabilized density ratio
#' \eqn{1 / f(A | L)} has a heavy right tail on its own, and multiplying it by a
#' second weight inherits that tail, leaving the product with no usable
#' variance. Build a continuous component with `stabilize = TRUE`. A binary or
#' categorical component needs no stabilization and is accepted either way.
#'
#' Any numerator that stabilizes the ratio satisfies the requirement, and any
#' density the ratio is read in is accepted: a component built with
#' `numerator = "integrated"`, with a heavier-tailed `.density`, or with a
#' `stabilization_score` of its own is stabilized as much as the default
#' marginal one is. What each component was built with is recorded per component
#' rather than merged, so the product carries both answers.
#'
#' # What the product records
#'
#' The result is a [psw()] with estimand `"ate"`, since the product of two `ate`
#' weights targets the joint `ate`. It is marked stabilized when *either*
#' component is, because a stabilizing numerator on either factor is a
#' stabilizing numerator on the product; `joint_wt_meta()` keeps the
#' per-component truth, so the coarse flag hides nothing.
#'
#' The product records nothing about trimming, truncation, or calibration of the
#' propensity scores the components were built from. Those records name the
#' observations of one component, and the product is not that component.
#'
#' @param w_a,w_e The two component weight vectors, as [psw()] objects. Both
#'   must target the `ate`, both must be the same length, and a component
#'   weighting a continuous exposure must be stabilized. The order is the
#'   factorization's order: `w_a` weights the first treatment and `w_e` weights
#'   the second, whose model conditions on the first.
#' @param exposure_type A character vector of length two naming each component's
#'   exposure type, one of `"binary"`, `"categorical"`, or `"continuous"`, in
#'   the order the components were given. Defaults to the types the components
#'   record, which is what a weight function writes on the weights it builds. A
#'   value given here is used instead of what they record, and a component
#'   recording no type, such as one assembled by hand, is refused unless both
#'   types are named.
#' @param x An object to test, or the weights to read the record from.
#'
#' @return
#' `wt_joint()` returns a [psw()] vector of the elementwise product, carrying a
#' `joint_wt_meta` record.
#'
#' `is_joint_wt()` returns a single logical.
#'
#' `joint_wt_meta()` returns the record as a list of three elements, each one
#' per component and in the order the components were given, or `NULL` for
#' weights that are not a product:
#' \describe{
#'   \item{`exposure_type`}{The two components' exposure types, in order.}
#'   \item{`stabilized`}{Whether each component was stabilized, in order.}
#'   \item{`density`}{Each component's density record, as [density_meta()]
#'     returns it, in order. A component weighting a discrete exposure is
#'     `NULL`, since its weights are no ratio of densities.}
#' }
#'
#' The record names the components rather than the observations, so it survives
#' subsetting and every other operation that keeps the vector a [psw()].
#'
#' @seealso [joint_wt_models()] for the container recording the two treatment
#'   models, and [wt_ate()] for building the components. These are used by
#'   [ipw()] methods for joint treatment models.
#'
#' @examples
#' set.seed(4)
#' n <- 200
#' x1 <- rnorm(n)
#' a <- rbinom(n, 1, plogis(0.3 * x1))
#' e <- rbinom(n, 1, plogis(-0.2 + 0.5 * x1 - 0.8 * a))
#' dat <- data.frame(x1, a, e)
#'
#' # The second model conditions on the first treatment, and does so flexibly
#' mod_a <- glm(a ~ x1, data = dat, family = binomial())
#' mod_e <- glm(e ~ a * x1, data = dat, family = binomial())
#'
#' w <- wt_joint(wt_ate(mod_a), wt_ate(mod_e))
#' head(w)
#' estimand(w)
#' is_joint_wt(w)
#' joint_wt_meta(w)
#'
#' # A continuous component must be stabilized
#' d <- 0.5 + 0.6 * x1 - 0.7 * a + rnorm(n)
#' mod_d <- lm(d ~ a * x1, data = dat)
#' w_d <- wt_ate(
#'   fitted(mod_d),
#'   d,
#'   exposure_type = "continuous",
#'   stabilize = TRUE
#' )
#' # Each component records the exposure type it weights, so the product knows
#' # which of them is the dose without being told
#' joint_wt_meta(wt_joint(wt_ate(mod_a), w_d))
#'
#' @name wt_joint
NULL

#' @rdname wt_joint
#' @export
wt_joint <- function(w_a, w_e, exposure_type = NULL) {
  # The class check comes first because everything after it reads a record off
  # the components, and a bare numeric carries none of them.
  check_wt_joint_class(w_a, w_e)
  exposure_type <- wt_joint_exposure_type(w_a, w_e, exposure_type)
  check_wt_joint_exposure_type(exposure_type)
  check_wt_joint_length(w_a, w_e)
  check_wt_joint_estimand(w_a, w_e)

  stabilized <- c(is_stabilized(w_a), is_stabilized(w_e))
  check_wt_joint_stabilize(exposure_type, stabilized)

  out <- new_psw(
    as.double(w_a) * as.double(w_e),
    estimand = "ate",
    # A stabilizing numerator on either factor is a stabilizing numerator on the
    # product, so the flag is the disjunction rather than the conjunction the
    # combining rules use. The record below keeps the per-component truth.
    stabilized = any(stabilized)
  )
  attr(out, "joint_wt_meta") <- list(
    exposure_type = unname(exposure_type),
    stabilized = unname(stabilized),
    # A dose's weights are a ratio of densities, and which ratio they are is
    # what an estimator has to rebuild them from. The slot is one per component,
    # empty for a component weighting no density, so the record is read
    # positionally either way.
    density = list(density_meta(w_a), density_meta(w_e))
  )

  out
}

# The exposure types the product is built under: the ones the caller named, or
# the ones the components recorded when the caller named none. A weight function
# records the type it was given, so repeating it is the exception rather than
# the rule, and a caller who does name it is taken at their word: the argument
# decides and the reading is what happens in its absence.
#
# A component recording no type is refused rather than guessed at. The
# stabilization requirement is a rule about doses, so a component whose type is
# unknown is one the rule could not be applied to.
wt_joint_exposure_type <- function(
  w_a,
  w_e,
  declared,
  call = rlang::caller_env()
) {
  if (!is.null(declared)) {
    return(declared)
  }

  recorded <- list(exposure_type(w_a), exposure_type(w_e))
  unrecorded <- vapply(recorded, is.null, logical(1))

  if (!any(unrecorded)) {
    return(unlist(recorded))
  }

  bad <- c("w_a", "w_e")[unrecorded]

  abort(
    c(
      "{.fun wt_joint} requires each component to record the exposure type it \\
      weights.",
      x = "{.arg {bad}} record{?s/} none.",
      i = "A weight built by hand, or by a version of the package that did \\
      not record it, carries an estimand and a stabilization status and \\
      nothing about the exposure, so the requirement that a continuous \\
      component be stabilized could not be applied to it.",
      i = "Rebuild {cli::qty(bad)}{?it/them} with a weight function such as \\
      {.fun wt_ate}, or name both types in {.arg exposure_type}, for example \\
      {.code exposure_type = c(\"binary\", \"continuous\")}."
    ),
    error_class = "propensity_wt_joint_exposure_type_error",
    call = call
  )
}

#' @rdname wt_joint
#' @export
is_joint_wt <- function(x) {
  !is.null(joint_wt_meta(x))
}

#' @rdname wt_joint
#' @export
joint_wt_meta <- function(x) {
  attr(x, "joint_wt_meta")
}

check_wt_joint_exposure_type <- function(
  exposure_type,
  call = rlang::caller_env()
) {
  problem <- if (!is.character(exposure_type) || length(exposure_type) != 2) {
    "{.arg exposure_type} names {length(exposure_type)} exposure type{?s}."
  } else if (!all(exposure_type %in% joint_wt_exposure_types)) {
    unsupported <- setdiff(exposure_type, joint_wt_exposure_types)
    "{.val {unsupported}} {?is/are} not {?a/} supported exposure type{?s}."
  } else {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.arg exposure_type} must name the exposure type of each component.",
      x = problem,
      i = "Supported types: {.val {joint_wt_exposure_types}}.",
      i = "Supply one per component, in the order the components were given, \\
      for example {.code exposure_type = c(\"binary\", \"continuous\")}, or \\
      leave {.arg exposure_type} unset to read each component's own record."
    ),
    error_class = "propensity_wt_joint_exposure_type_error",
    call = call
  )
}

check_wt_joint_class <- function(w_a, w_e, call = rlang::caller_env()) {
  bad <- c("w_a", "w_e")[!c(is_psw(w_a), is_psw(w_e))]

  if (length(bad) == 0) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun wt_joint} requires both components to be {.cls psw} vectors.",
      x = "{.arg {bad}} {?is/are} not.",
      i = "A bare numeric records neither an estimand nor a stabilization \\
      status, so the checks the product depends on cannot be applied to it and \\
      the result could record nothing about where it came from.",
      i = "Build each component with a weight function such as {.fun wt_ate}."
    ),
    error_class = "propensity_wt_joint_class_error",
    call = call
  )
}

check_wt_joint_length <- function(w_a, w_e, call = rlang::caller_env()) {
  if (identical(length(w_a), length(w_e))) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun wt_joint} requires its two components to be the same length.",
      x = "{.arg w_a} is length {length(w_a)}; {.arg w_e} is length \\
      {length(w_e)}.",
      i = "The product is elementwise, so each unit needs one weight from each \\
      component.",
      i = "Build both components from the same observations."
    ),
    error_class = "propensity_wt_joint_length_error",
    call = call
  )
}

check_wt_joint_estimand <- function(w_a, w_e, call = rlang::caller_env()) {
  estimands <- c(joint_wt_estimand_label(w_a), joint_wt_estimand_label(w_e))

  if (all(estimands == "ate")) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun wt_joint} requires both components to target the \\
      {.val {'ate'}} estimand.",
      x = "{.arg w_a} targets {.val {estimands[[1]]}} and {.arg w_e} targets \\
      {.val {estimands[[2]]}}.",
      i = "The product targets the joint {.val {'ate'}}, and a tilted \\
      component targets a population the other one does not, so the product \\
      targets no single population.",
      i = "Rebuild both components with {.fun wt_ate}."
    ),
    error_class = "propensity_wt_joint_estimand_error",
    call = call
  )
}

# The estimand a component records, or a word standing for the absence of one,
# so a message naming what each component targets says something either way.
joint_wt_estimand_label <- function(w) {
  label <- estimand(w)

  if (is.null(label)) "none" else label
}

check_wt_joint_stabilize <- function(
  exposure_type,
  stabilized,
  call = rlang::caller_env()
) {
  continuous <- exposure_type == "continuous"
  bad <- c("w_a", "w_e")[continuous & !stabilized]

  if (length(bad) == 0) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun wt_joint} requires a continuous component to be stabilized.",
      x = "{.arg {bad}} weight{?s/} a continuous exposure and {?is/are} not \\
      stabilized.",
      i = "The unstabilized density ratio has a heavy right tail on its own, \\
      and multiplying it by a second weight inherits that tail, leaving the \\
      product with no usable variance.",
      i = "Rebuild {cli::qty(bad)}{?it/them} with {.code stabilize = TRUE}. A \\
      binary or categorical component needs no stabilization."
    ),
    error_class = "propensity_wt_joint_stabilize_error",
    call = call
  )
}

# ---- joint_wt_models() ------------------------------------------------------

#' Record the two treatment models of a joint exposure
#'
#' @description
#' `joint_wt_models()` records the two fitted treatment models a joint weight is
#' built from, in the order the joint density factorizes: the first treatment's
#' model, then the second treatment's model, which conditions on the first.
#' `is_joint_wt_models()` reports whether an object is such a record.
#'
#' The container exists to refuse. It checks that the pair describes a
#' sequential factorization at all, which is the one thing about a joint weight
#' that cannot be detected from the weight vector afterwards.
#'
#' @details
#' # What is checked
#'
#' Each argument is named for the treatment its model fits, and the name has to
#' be that model's response: the factorization check looks for the first
#' treatment by name in the second model's formula, so a model recorded under
#' the wrong name would be checked for the wrong variable.
#'
#' The second model must condition on the first treatment. The product of two
#' marginal models, \eqn{f(A | L) f(E | L)}, is not the joint weight
#' \eqn{f(A | L) f(E | A, L)} for any data in which the second treatment depends
#' on the first, and nothing downstream can tell the two apart. What is checked
#' is whether the second model's right-hand side *mentions* the first treatment,
#' so `e ~ factor(a) * x1` and `e ~ x1 + a:x1` both pass, as does the additive
#' `e ~ a + x1`.
#'
#' Passing is not the same as modeling the dependence well. An additive term
#' satisfies the factorization and may still model it badly, and a badly modeled
#' dependence biases the joint weight invisibly. Model the dependence on the
#' first treatment flexibly, with an interaction or a spline, unless you have a
#' reason to believe it is additive.
#'
#' Neither may the two models condition on each other. A sequential
#' factorization has a first factor that is marginal in the second treatment, so
#' a pair that each read the other's treatment are not the two factors of any
#' factorization, in either order.
#'
#' @param ... Exactly two fitted treatment models, each named for the treatment
#'   it fits. The order is the factorization's order. Supported models are a
#'   binomial [stats::glm()] for a binary treatment, a [nnet::multinom()] for a
#'   categorical one, and an [stats::lm()], a gaussian [stats::glm()], a
#'   gaussian `mgcv::gam()`, or a `MASS::rlm()` for a continuous one. Each class
#'   is recognized by name rather than by inheritance, so a subclass of one of
#'   them is not read as its parent.
#' @param x An object to test.
#'
#' @return
#' `joint_wt_models()` returns an S3 object of class `joint_wt_models`, a list
#' with three fields:
#' \describe{
#'   \item{`names`}{The two treatment names, in the factorization's order.}
#'   \item{`models`}{The two fitted models, named for their treatments.}
#'   \item{`exposure_type`}{Each model's exposure type, named for its treatment:
#'     `"binary"`, `"categorical"`, or `"continuous"`.}
#' }
#'
#' `is_joint_wt_models()` returns a single logical.
#'
#' @seealso [wt_joint()] for the product weight itself. These are used by
#'   [ipw()] methods for joint treatment models.
#'
#' @examples
#' set.seed(4)
#' n <- 200
#' x1 <- rnorm(n)
#' a <- rbinom(n, 1, plogis(0.3 * x1))
#' e <- rbinom(n, 1, plogis(-0.2 + 0.5 * x1 - 0.8 * a))
#' dat <- data.frame(x1, a, e)
#'
#' models <- joint_wt_models(
#'   a = glm(a ~ x1, data = dat, family = binomial()),
#'   e = glm(e ~ a * x1, data = dat, family = binomial())
#' )
#' models$names
#' models$exposure_type
#' is_joint_wt_models(models)
#'
#' @name joint_wt_models
NULL

#' @rdname joint_wt_models
#' @export
joint_wt_models <- function(...) {
  models <- list(...)
  names <- names(models)

  check_joint_wt_models_named(models, names)

  exposure_type <- check_joint_wt_models_supported(models, names)
  check_joint_wt_models_response(models, names)
  check_joint_wt_models_factorization(models, names)

  structure(
    list(
      names = names,
      models = models,
      exposure_type = stats::setNames(exposure_type, names)
    ),
    class = "joint_wt_models"
  )
}

#' @rdname joint_wt_models
#' @export
is_joint_wt_models <- function(x) {
  inherits(x, "joint_wt_models")
}

# Exactly two models, each under a name of its own. The names are the treatment
# names the factorization is written in, so an unnamed argument names no
# treatment and a repeated name names one treatment twice, which is not a
# crossing and which the container keyed by name could not tell apart.
check_joint_wt_models_named <- function(
  models,
  names,
  call = rlang::caller_env()
) {
  n <- length(models)

  # A call that names none of its arguments has no name vector at all rather
  # than a vector of empty strings, so the count has to come from the models
  # there. Read off the names alone it collapses to zero by zero-length
  # recycling, and the message reports that none of them is unnamed.
  n_unnamed <- if (is.null(names)) {
    n
  } else {
    sum(is.na(names) | names == "")
  }

  problem <- if (n != 2) {
    "{n} model{?s} {?was/were} supplied."
  } else if (n_unnamed > 0) {
    "{n_unnamed} of them {?is/are} unnamed."
  } else if (anyDuplicated(names)) {
    duplicated_name <- names[[anyDuplicated(names)]]
    "{.val {duplicated_name}} names both of them."
  } else {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun joint_wt_models} must be given exactly two models, each named \\
      for the treatment it fits.",
      x = problem,
      i = "The names are the treatment names, and their order is the \\
      factorization's order: the first treatment is the one the second model \\
      conditions on."
    ),
    error_class = "propensity_wt_joint_models_error",
    call = call
  )
}

# Each model's exposure type, which is also the check that the model is one this
# package can read a treatment density from at all.
check_joint_wt_models_supported <- function(
  models,
  names,
  call = rlang::caller_env()
) {
  types <- lapply(models, joint_wt_model_type)
  unsupported <- vapply(types, is.null, logical(1))

  if (any(unsupported)) {
    bad <- names[unsupported]
    bad_class <- class(models[[which(unsupported)[[1]]]])[[1]]
    abort(
      c(
        "{.fun joint_wt_models} must be given models it can read a treatment \\
        density from.",
        x = "The model named {.arg {bad[[1]]}} is {.cls {bad_class}}.",
        i = "Supported models: a binomial {.fun glm} for a binary treatment, \\
        a {.fun nnet::multinom} for a categorical one, and an {.fun lm}, a \\
        gaussian {.fun glm}, a gaussian {.fun mgcv::gam}, or a \\
        {.fun MASS::rlm} for a continuous one."
      ),
      error_class = "propensity_wt_joint_models_error",
      call = call
    )
  }

  unlist(types)
}

check_joint_wt_models_response <- function(
  models,
  names,
  call = rlang::caller_env()
) {
  mismatched <- !vapply(
    seq_along(models),
    function(i) identical(joint_wt_response(models[[i]]), names[[i]]),
    logical(1)
  )

  if (!any(mismatched)) {
    return(invisible(TRUE))
  }

  first <- which(mismatched)[[1]]
  name <- names[[first]]
  response <- deparse1(stats::formula(models[[first]])[[2L]])

  abort(
    c(
      "{.fun joint_wt_models} requires each model's response to be the \\
      treatment it is named for.",
      x = "The model named {.arg {name}} has response {.code {response}}.",
      i = "The factorization check looks for the first treatment by name in \\
      the second model's formula, so a model recorded under the wrong name \\
      would be checked for the wrong variable.",
      i = "Rename the argument to the treatment the model fits, or supply the \\
      model that fits {.val {name}}."
    ),
    error_class = "propensity_wt_joint_response_error",
    call = call
  )
}

# The load-bearing guardrail, and the pair of refusals it splits into. The two
# are mutually exclusive: a circular pair is one in which the second model does
# read the first treatment, which is the condition the factorization check
# passes on. They carry their own classes because the remedies differ, one
# adding a term and the other removing one.
check_joint_wt_models_factorization <- function(
  models,
  names,
  call = rlang::caller_env()
) {
  first <- names[[1]]
  second <- names[[2]]
  first_reads_second <- second %in% joint_wt_reads(models[[1]])
  second_reads_first <- first %in% joint_wt_reads(models[[2]])

  if (first_reads_second && second_reads_first) {
    abort(
      c(
        "{.fun joint_wt_models} requires a sequential factorization, and each \\
        of these models conditions on the other's treatment.",
        x = "{.arg {first}} reads {.val {second}}, and {.arg {second}} reads \\
        {.val {first}}.",
        i = "A sequential factorization has a first factor that is marginal \\
        in the second treatment, so such a pair are not the two factors of \\
        any factorization, in either order.",
        i = "Refit the first treatment's model without {.val {second}}."
      ),
      error_class = "propensity_wt_joint_circular_error",
      call = call
    )
  }

  if (second_reads_first) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.fun joint_wt_models} requires the second model to condition on the \\
      first treatment.",
      x = "{.arg {second}} does not read {.val {first}} on its right-hand \\
      side.",
      i = "A joint weight factorizes as f({first} | L) f({second} | {first}, \\
      L). The product of two marginal models, f({first} | L) f({second} | L), \\
      is a different quantity, and it is not the joint weight wherever \\
      {.val {second}} depends on {.val {first}}.",
      i = "Nothing downstream can tell the two apart: the product is an \\
      ordinary vector of positive numbers either way.",
      i = "Add {.val {first}} to the formula of {.arg {second}}, and model \\
      that dependence flexibly rather than as a single additive term."
    ),
    error_class = "propensity_wt_joint_factorization_error",
    call = call
  )
}
