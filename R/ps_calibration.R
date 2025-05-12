#' Calibrate propensity scores
#'
#' This function calibrates propensity scores to improve their accuracy
#' using Platt scaling (logistic regression).
#' It preserves the attributes of causal weight objects when applicable.
#'
#' @param ps Numeric vector of propensity scores between 0 and 1
#' @param treat A binary vector of treatment assignments
#' @param .treated Value that represents the treated units in the `treat` vector
#' @param .untreated Value that represents the untreated units in the `treat` vector
#' @param estimand Character indicating the estimand type.
#'
#' @return A calibrated propensity score object (`psw`)
#'
#' @examples
#' # Generate example data
#' ps <- runif(100)
#' treat <- rbinom(100, 1, ps)
#'
#' # Calibrate using Platt scaling
#' calibrated <- ps_calibrate(ps, treat)
#' @export
ps_calibrate <- function(ps, treat,
                         .treated = 1, .untreated = 0,
                         estimand = NULL) {
  # Check that ps is numeric and in valid range
  if (!is.numeric(ps)) {
    abort("`ps` must be a numeric vector.")
  }

  if (is_causal_wt(ps) && is_ps_calibrated(ps)) {
    abort("`ps` already calibrated.")
  }

  if (any(ps < 0 | ps > 1, na.rm = TRUE)) {
    abort("`ps` values must be between 0 and 1.")
  }

  # Transform treatment to binary if needed
  treat <- transform_exposure_binary(
    treat,
    .treated = .treated,
    .untreated = .untreated
  )

  if (length(ps) != length(treat)) {
    abort("Propensity score vector `ps` must be the same length as `treat`.")
  }

  # Extract attributes from causal weight objects if applicable
  if (is_causal_wt(ps)) {
    if (is.null(estimand)) estimand <- estimand(ps)
    stabilized <- is_stabilized(ps)
    trimmed <- is_ps_trimmed(ps)
    truncated <- is_ps_truncated(ps)
  } else {
    if (is.null(estimand)) estimand <- "unknown"
    stabilized <- FALSE
    trimmed <- FALSE
    truncated <- FALSE
  }

  # Fit logistic regression: treat ~ ps (using ps as sole predictor)
  calib_model <- stats::glm(treat ~ ps, family = stats::binomial())

  # Calibrated probabilities are the fitted values from this model
  calib_ps <- stats::fitted(calib_model)

  # Return calibrated scores as a psw object with preserved attributes
  psw(
    x = calib_ps,
    estimand = estimand,
    stabilized = stabilized,
    trimmed = trimmed,
    truncated = truncated,
    calibrated = TRUE
  )
}

#' @rdname psw
#' @export
is_ps_calibrated <- function(wt) {
  isTRUE(attr(wt, "calibrated"))
}

