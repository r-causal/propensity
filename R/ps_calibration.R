#' Calibrate propensity scores
#'
#' This function calibrates propensity scores to improve their accuracy using
#' either Platt scaling (logistic regression) or isotonic regression. It
#' preserves the attributes of causal weight objects when applicable.
#'
#' @param ps Numeric vector of propensity scores between 0 and 1
#' @param treat A binary vector of treatment assignments
#' @param method Calibration method:
#'   \describe{
#'     \item{`"logistic"`}{(Default) Logistic calibration (also known as Platt scaling).
#'       Assumes a sigmoid relationship between observed and true probabilities.
#'       Best when: propensity scores follow a logistic pattern but are
#'       systematically biased. Provides smooth, parametric calibration.
#'       Faster and more stable with small samples.}
#'     \item{`"isoreg"`}{Isotonic regression calibration. Uses a non-parametric
#'       monotonic transformation. Best when: the relationship between observed
#'       and true probabilities is non-linear or when you want to preserve
#'       the rank order without assuming a specific functional form.
#'       More flexible but requires larger samples for stable estimates.}
#'   }
#' @param smooth Logical. For `method = "logistic"`, whether to use a smoothed
#'   logistic spline model (`smooth = TRUE`, default) or simple logistic
#'   regression (`smooth = FALSE`). When `TRUE`, uses `mgcv::gam()` with
#'   spline smoothing. When `FALSE`, uses `stats::glm()`. Ignored for
#'   `method = "isoreg"`.
#' @param .treated The value representing the treatment group. If not provided,
#'   `ps_calibrate()` will attempt to automatically determine the treatment coding.
#' @param .untreated The value representing the control group. If not provided,
#'   `ps_calibrate()` will attempt to automatically determine the control coding.
#' @param estimand Character indicating the estimand type.
#'
#' @return A calibrated propensity score object (`psw`)
#'
#' @examples
#' # Generate example data
#' ps <- runif(100)
#' treat <- rbinom(100, 1, ps)
#'
#' # Logistic calibration with smoothing (default)
#' calibrated_smooth <- ps_calibrate(ps, treat)
#'
#' # Logistic calibration without smoothing (simple logistic regression)
#' calibrated_simple <- ps_calibrate(ps, treat, smooth = FALSE)
#'
#' # Isotonic regression
#' calibrated_iso <- ps_calibrate(ps, treat, method = "isoreg")
#' @importFrom stats glm fitted isoreg binomial
#' @export
ps_calibrate <- function(
  ps,
  treat,
  method = c("logistic", "isoreg"),
  smooth = TRUE,
  .treated = NULL,
  .untreated = NULL,
  estimand = NULL
) {
  method <- match.arg(method)
  # Check that ps is numeric and in valid range
  if (!is.numeric(ps)) {
    abort("`ps` must be a numeric vector.")
  }

  if (is_causal_wt(ps) && is_ps_calibrated(ps)) {
    abort(
      "`ps` is already calibrated. Cannot calibrate already calibrated propensity scores."
    )
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
  calib_model <- NULL
  if (method == "logistic") {
    if (smooth) {
      # Check if mgcv is available for smooth calibration
      rlang::check_installed(
        "mgcv",
        reason = "for smooth calibration using GAM (Generalized Additive Model)."
      )

      # Create data frame for GAM fitting (only non-NA values)
      calib_data <- data.frame(
        treat = treat[!na_idx],
        ps = ps[!na_idx]
      )

      # Check if we have enough unique values for smoothing (like probably does)
      n_unique <- length(unique(calib_data$ps))
      if (n_unique < 10) {
        # Fall back to simple logistic regression if too few unique values
        smooth <- FALSE
        calib_model <- stats::glm(
          treat ~ ps,
          data = calib_data,
          family = stats::binomial()
        )
      } else {
        # Fit GAM with smooth spline: treat ~ s(ps)
        calib_model <- mgcv::gam(
          treat ~ s(ps),
          data = calib_data,
          family = "binomial"
        )
      }
    }

    if (!smooth) {
      # For simple logistic regression, fit on original data (not data frame)
      # This handles the case where smooth was originally FALSE or was set to FALSE due to fallback
      if (is.null(calib_model)) {
        calib_model <- stats::glm(treat ~ ps, family = stats::binomial())
      }
    }

    if (isFALSE(calib_model$converged)) {
      warn(
        "Calibration model did not converge",
        warning_class = "propensity_convergence_warning"
      )
    }

    # Calibrated probabilities are the fitted values from this model
    # This will be shorter if there are NAs, so we need to expand back
    if (smooth) {
      # For GAM models, predict on the full ps vector (including NAs)
      # Create prediction data frame
      pred_data <- data.frame(ps = ps)
      fitted_vals <- as.numeric(predict(
        calib_model,
        newdata = pred_data,
        type = "response"
      ))
      calib_ps <- fitted_vals
    } else {
      # For GLM models, use fitted values and expand back
      fitted_vals <- stats::fitted(calib_model)
      calib_ps <- numeric(length(ps))
      calib_ps[!na_idx] <- fitted_vals
      calib_ps[na_idx] <- NA
    }
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
