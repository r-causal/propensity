#' Calibrate propensity scores
#'
#' This function calibrates propensity scores to improve their accuracy using
#' either Platt scaling (logistic regression) or isotonic regression. It
#' preserves the attributes of causal weight objects when applicable.
#'
#' @param ps Numeric vector of propensity scores between 0 and 1
#' @param treat A binary vector of treatment assignments
#' @param method Calibration method, either "platt" (default) for logistic
#'   calibration or "isoreg" for isotonic regression calibration
#' @param .treated Value that represents the treated units in the `treat` vector
#' @param .untreated Value that represents the untreated units in the `treat`
#'   vector
#' @param estimand Character indicating the estimand type.
#'
#' @return A calibrated propensity score object (`psw`)
#'
#' @examples
#' # Generate example data
#' ps <- runif(100)
#' treat <- rbinom(100, 1, ps)
#'
#' # Platt scaling (default)
#' calibrated_platt <- ps_calibrate(ps, treat)
#'
#' # Isotonic regression
#' calibrated_iso <- ps_calibrate(ps, treat, method = "isoreg")
#' @importFrom stats glm fitted isoreg binomial
#' @export
ps_calibrate <- function(
  ps,
  treat,
  method = c("platt", "isoreg"),
  .treated = 1,
  .untreated = 0,
  estimand = NULL
) {
  method <- match.arg(method)
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

  # Handle NA values
  na_idx <- is.na(ps) | is.na(treat)

  # Perform calibration based on method
  if (method == "platt") {
    # Fit logistic regression: treat ~ ps (using ps as sole predictor)
    calib_model <- stats::glm(treat ~ ps, family = stats::binomial())

    if (!calib_model$converged) {
      warn(
        "Calibration model did not converge",
        warning_class = "propensity_convergence_warning"
      )
    }

    # Calibrated probabilities are the fitted values from this model
    # This will be shorter if there are NAs, so we need to expand back
    fitted_vals <- stats::fitted(calib_model)
    calib_ps <- numeric(length(ps))
    calib_ps[!na_idx] <- fitted_vals
    calib_ps[na_idx] <- NA
  } else if (method == "isoreg") {
    # Isotonic regression calibration
    if (any(na_idx)) {
      # Work with non-NA values only
      ps_valid <- ps[!na_idx]
      treat_valid <- treat[!na_idx]

      # Order by propensity scores for isotonic regression
      ord <- order(ps_valid)
      ps_ordered <- ps_valid[ord]
      treat_ordered <- treat_valid[ord]

      # Fit isotonic regression
      iso_fit <- stats::isoreg(ps_ordered, treat_ordered)

      # Get calibrated values and map back to original order
      calib_ps_ordered <- iso_fit$yf
      calib_ps_valid <- numeric(length(ps_valid))
      calib_ps_valid[ord] <- calib_ps_ordered

      # Map back to full vector with NAs
      calib_ps <- numeric(length(ps))
      calib_ps[!na_idx] <- calib_ps_valid
      calib_ps[na_idx] <- NA
    } else {
      # No NAs, proceed normally
      ord <- order(ps)
      ps_ordered <- ps[ord]
      treat_ordered <- treat[ord]

      # Fit isotonic regression
      iso_fit <- stats::isoreg(ps_ordered, treat_ordered)

      # Get calibrated values and map back to original order
      calib_ps_ordered <- iso_fit$yf
      calib_ps <- numeric(length(ps))
      calib_ps[ord] <- calib_ps_ordered
    }

    # Ensure calibrated values are in [0, 1]
    calib_ps[!na_idx] <- pmax(0, pmin(1, calib_ps[!na_idx]))
  }

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
