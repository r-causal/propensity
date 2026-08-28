# A continuous-exposure fit whose outcome model reads the exposure through more
# than one design column, which is the boundary the reported rows are named by:
# one column is a slope, several are coefficients of a curve.
ipw_is_dose_basis <- function(spec) {
  length(spec$coefficients$col) > 1L
}

# The claim `ipw()` makes about such a fit, written once. A curve has a
# different slope at every dose, so no coefficient of it is the effect of the
# exposure and there is no row of the surface that answers the question the
# marginal reading is asked.
ipw_dose_basis_lead <- function() {
  "{.fun ipw} reports only the conditional reading because the exposure \\
  enters {.arg outcome_mod} through more than one term, such as a spline or \\
  polynomial."
}

# Why there is one reading, and where the estimand this package does not compute
# belongs. Marginalizing a dose response over the observed doses is a separate
# quantity, and assembling it from the coefficients by hand is the part a user
# should not be left to do.
ipw_dose_basis_reason <- function() {
  c(
    "With a nonlinear dose-response, no single coefficient is the effect of \\
    the exposure, so there is no marginal effect to report.",
    "Use the {.pkg marginaleffects} package to marginalize over the dose: \\
    {.fun avg_slopes} for slopes, {.fun avg_comparisons} for contrasts, and \\
    {.fun avg_predictions} for causal dose-response functions. See \\
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

# The refusal `ipw()` makes when it is asked to build such a fit in a reading
# the fit does not have. `NULL` is a caller who named no reading, and the
# conditional one is the reading the result records, so only the marginal
# reading is refused. Every other door into that reading is refused by the
# result class, which reads the set of readings the fit declares.
check_ipw_dose_basis_effects <- function(effects, call = rlang::caller_env()) {
  if (!identical(effects, "marginal")) {
    return(invisible(effects))
  }

  reason <- ipw_dose_basis_reason()

  abort(
    c(
      ipw_dose_basis_lead(),
      x = "{.code effects = \"marginal\"} is not available for this model.",
      i = reason[[1]],
      i = reason[[2]],
      i = "Use {.code effects = \"conditional\"} or omit {.arg effects}."
    ),
    error_class = "propensity_ipw_effects_error",
    call = call
  )
}
