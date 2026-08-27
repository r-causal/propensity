# A continuous-exposure fit whose outcome model reads the exposure through more
# than one design column, which is the boundary the reported rows are named by:
# one column is a slope, several are coefficients of a curve.
ipw_is_dose_basis <- function(spec) {
  length(spec$coefficients$col) > 1L
}

# The claim the constructor and every accessor make about such a fit, written
# once. A curve has a different slope at every dose, so no coefficient of it is
# the effect of the exposure and there is no row of the surface that answers the
# question the marginal reading is asked.
ipw_dose_basis_lead <- function() {
  "{.fun ipw} reports only the conditional reading when the exposure enters \\
  {.arg outcome_mod} through several columns."
}

# Why there is one reading, and where the estimand this package does not compute
# belongs. Marginalizing a dose response over the observed doses is a separate
# quantity, and assembling it from the coefficients by hand is the part a user
# should not be left to do.
ipw_dose_basis_reason <- function() {
  c(
    "The coefficient surface is the outcome model's own, and no single row of \\
    it is a causal effect.",
    "Marginalizing over the dose is left to the {.pkg marginaleffects} \\
    package: call {.fun avg_slopes} or {.fun avg_comparisons} on the \\
    conditional result. See \\
    {.url https://marginaleffects.com/chapters/interactions.html}."
  )
}

# The reading such a fit records is a default rather than a fact the caller
# asked for, so it is announced once at construction. A caller who named the
# reading has already been told, and the announcement stops.
ipw_dose_basis_announce <- function() {
  reason <- ipw_dose_basis_reason()

  alert_info(ipw_dose_basis_lead())
  alert_info(reason[[1]])
  alert_info(reason[[2]])
  alert_info("Set {.code effects = \"conditional\"} to silence this message.")

  invisible(NULL)
}

# The refusal every door into the marginal reading of a dose-basis fit makes.
# `NULL` is a caller who named no reading, and the conditional one is the
# reading the result records, so only the marginal reading is refused.
check_ipw_dose_basis_effects <- function(effects, call = rlang::caller_env()) {
  if (!identical(effects, "marginal")) {
    return(invisible(effects))
  }

  reason <- ipw_dose_basis_reason()

  abort(
    c(
      ipw_dose_basis_lead(),
      x = "{.code effects = \"marginal\"} asks for a reading this model has \\
      none of.",
      i = reason[[1]],
      i = reason[[2]],
      i = "Use {.code effects = \"conditional\"} or omit {.arg effects}."
    ),
    error_class = "propensity_ipw_effects_error",
    call = call
  )
}

# The subclass carries the refusals and nothing else. Every generic without a
# method here resolves at `ipw`, which is what keeps a dose-basis result
# printable, poolable, and readable by everything written against the shared
# class. `as_conditional()` is among them: the reading such a result records is
# the one it presents, so the upstream method already answers.

#' @exportS3Method causalgenerics::as_marginal ipw_dose_basis
as_marginal.ipw_dose_basis <- function(x, ...) {
  check_ipw_dose_basis_effects("marginal")
}

#' @exportS3Method stats::coef ipw_dose_basis
coef.ipw_dose_basis <- function(object, ..., effects = NULL) {
  check_ipw_dose_basis_effects(effects)
  NextMethod()
}

#' @exportS3Method stats::vcov ipw_dose_basis
vcov.ipw_dose_basis <- function(object, ..., effects = NULL) {
  check_ipw_dose_basis_effects(effects)
  NextMethod()
}

#' @exportS3Method stats::confint ipw_dose_basis
confint.ipw_dose_basis <- function(
  object,
  parm,
  level = 0.95,
  ...,
  effects = NULL
) {
  check_ipw_dose_basis_effects(effects)
  NextMethod()
}

#' @exportS3Method generics::tidy ipw_dose_basis
tidy.ipw_dose_basis <- function(x, ..., effects = NULL) {
  check_ipw_dose_basis_effects(effects)
  NextMethod()
}

# The coercion surface upstream takes no reading at all, so the argument is
# added here rather than passed on: a caller who names the marginal reading is
# refused, and a call that names none behaves exactly as the upstream method
# does.
#' @export
as.data.frame.ipw_dose_basis <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  effects = NULL
) {
  check_ipw_dose_basis_effects(effects)
  NextMethod()
}
