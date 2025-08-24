#' Calibrate propensity scores
#'
#' This function calibrates propensity scores to improve their accuracy using
#' either Platt scaling (logistic regression) or isotonic regression. It
#' preserves the attributes of causal weight objects when applicable.
#'
#' @param ps Numeric vector of propensity scores between 0 and 1
#' @param .exposure A binary vector of treatment assignments
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
#' @param .focal_level The value representing the focal group (typically treatment).
#'   If not provided, `ps_calibrate()` will attempt to automatically determine the coding.
#' @param .reference_level The value representing the reference group (typically control).
#'   If not provided, `ps_calibrate()` will attempt to automatically determine the coding.
#' @param .treated `r lifecycle::badge("deprecated")` Use `.focal_level` instead.
#' @param .untreated `r lifecycle::badge("deprecated")` Use `.reference_level` instead.
#'
#' @return A calibrated propensity score object (`ps_calib`)
#'
#' @examples
#' # Generate example data
#' ps <- runif(100)
#' exposure <- rbinom(100, 1, ps)
#'
#' # Logistic calibration with smoothing (default)
#' calibrated_smooth <- ps_calibrate(ps, exposure)
#'
#' # Logistic calibration without smoothing (simple logistic regression)
#' calibrated_simple <- ps_calibrate(ps, exposure, smooth = FALSE)
#'
#' # Isotonic regression
#' calibrated_iso <- ps_calibrate(ps, exposure, method = "isoreg")
#' @importFrom stats glm fitted isoreg binomial
#' @export
ps_calibrate <- function(
  ps,
  .exposure,
  method = c("logistic", "isoreg"),
  smooth = TRUE,
  .focal_level = NULL,
  .reference_level = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  method <- rlang::arg_match(method)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "ps_calibrate"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  # Check that ps is numeric and in valid range
  if (!is.numeric(ps)) {
    abort(
      "`ps` must be a numeric vector.",
      error_class = "propensity_type_error"
    )
  }

  if (is_ps_calibrated(ps)) {
    abort(
      "`ps` is already calibrated. Cannot calibrate already calibrated propensity scores.",
      error_class = "propensity_already_calibrated_error"
    )
  }

  # Extract numeric values for comparison if ps is a causal weight object
  ps_numeric <- if (is_causal_wt(ps)) {
    as.numeric(ps)
  } else if (inherits(ps, "ps_calib")) {
    vec_data(ps)
  } else {
    ps
  }

  if (any(ps_numeric < 0 | ps_numeric > 1, na.rm = TRUE)) {
    abort(
      "`ps` values must be between 0 and 1.",
      error_class = "propensity_range_error"
    )
  }

  # Transform treatment to binary if needed
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  if (length(ps) != length(.exposure)) {
    abort(
      "Propensity score vector `ps` must be the same length as `.exposure`.",
      error_class = "propensity_length_error"
    )
  }

  # ps_calib objects should not have weight-specific attributes like estimand
  # Only store calibration-specific metadata

  # Handle NA values
  na_idx <- is.na(ps) | is.na(.exposure)

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
        treat = .exposure[!na_idx],
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
        calib_model <- stats::glm(.exposure ~ ps, family = stats::binomial())
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
    # Smoothing is not applicable to isotonic regression
    smooth <- FALSE
    # Isotonic regression calibration
    if (any(na_idx)) {
      # Work with non-NA values only
      ps_valid <- ps[!na_idx]
      treat_valid <- .exposure[!na_idx]

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
      treat_ordered <- .exposure[ord]

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

  # Return calibrated scores as a ps_calib object
  new_ps_calib(
    x = calib_ps,
    ps_calib_meta = list(
      method = method,
      smooth = smooth,
      calib_model = calib_model
    )
  )
}

#' Constructor for ps_calib objects
#' @noRd
new_ps_calib <- function(x, ps_calib_meta = list()) {
  if (is.matrix(x)) {
    # For matrices, we don't use vctrs
    structure(
      x,
      ps_calib_meta = ps_calib_meta,
      class = c("ps_calib_matrix", "ps_calib", "matrix")
    )
  } else {
    vec_assert(x, ptype = double())
    new_vctr(
      x,
      ps_calib_meta = ps_calib_meta,
      class = "ps_calib",
      inherit_base_type = TRUE
    )
  }
}

#' Extract metadata from ps_calib object
#' @noRd
ps_calib_meta <- function(x) {
  attr(x, "ps_calib_meta")
}

#' Check if object is calibrated
#'
#' @description
#' `is_ps_calibrated()` is an S3 generic that returns `TRUE` if its argument represents
#' calibrated propensity scores.
#'
#' @param x An R object.
#' @return A logical scalar (`TRUE` or `FALSE`).
#' @export
is_ps_calibrated <- function(x) {
  UseMethod("is_ps_calibrated")
}

#' @export
is_ps_calibrated.default <- function(x) {
  FALSE
}

#' @export
is_ps_calibrated.psw <- function(x) {
  isTRUE(attr(x, "calibrated"))
}

#' @export
is_ps_calibrated.ps_calib <- function(x) {
  TRUE
}

#' @export
is_ps_calibrated.ps_calib_matrix <- function(x) {
  TRUE
}

# vctrs machinery for ps_calib

#' @export
vec_ptype_abbr.ps_calib <- function(x, ...) "ps_calib"

#' @export
vec_ptype_full.ps_calib <- function(x, ...) {
  meta <- ps_calib_meta(x)
  paste("ps_calib", meta$method)
}

#' @export
#' @method vec_arith ps_calib
vec_arith.ps_calib <- function(op, x, y, ...) {
  UseMethod("vec_arith.ps_calib", y)
}

#' @export
#' @method vec_arith.ps_calib default
vec_arith.ps_calib.default <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_calib ps_calib
vec_arith.ps_calib.ps_calib <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_calib MISSING
vec_arith.ps_calib.MISSING <- function(op, x, y, ...) {
  switch(
    op,
    `-` = -1 * vec_data(x), # Returns numeric
    `+` = vec_data(x), # Returns numeric
    stop_incompatible_op(op, x, y)
  )
}

#' @export
#' @method vec_arith.ps_calib numeric
vec_arith.ps_calib.numeric <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_calib
vec_arith.numeric.ps_calib <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
vec_ptype2.ps_calib.ps_calib <- function(x, y, ...) {
  x_meta <- ps_calib_meta(x)
  y_meta <- ps_calib_meta(y)

  # Check if calibration parameters match
  if (
    !identical(x_meta$method, y_meta$method) ||
      !identical(x_meta$smooth, y_meta$smooth)
  ) {
    warn_incompatible_metadata(
      x,
      y,
      "different calibration parameters"
    )
    return(double())
  }

  # If parameters match, return ps_calib prototype
  new_ps_calib(double(), ps_calib_meta = x_meta)
}

#' @export
vec_ptype2.ps_calib.double <- function(x, y, ...) {
  warn_class_downgrade("ps_calib")
  double()
}

#' @export
vec_ptype2.double.ps_calib <- function(x, y, ...) {
  warn_class_downgrade("ps_calib")
  double()
}

#' @export
vec_cast.ps_calib.ps_calib <- function(x, to, ...) {
  # Check if metadata matches
  x_meta <- ps_calib_meta(x)
  to_meta <- ps_calib_meta(to)

  if (
    !identical(x_meta$method, to_meta$method) ||
      !identical(x_meta$smooth, to_meta$smooth)
  ) {
    vctrs::stop_incompatible_cast(x, to, x_arg = "", to_arg = "")
  }

  # Return x as-is if metadata matches
  x
}

#' @export
vec_cast.double.ps_calib <- function(x, to, ...) {
  # degrade to numeric
  vec_data(x)
}

#' @export
as.double.ps_calib <- function(x, ...) {
  vec_data(x)
}

#' @export
vec_cast.ps_calib.double <- function(x, to, ...) {
  # create a default ps_calib
  new_ps_calib(
    x,
    ps_calib_meta = list(
      method = "unknown",
      smooth = FALSE
    )
  )
}

#' @export
vec_ptype2.psw.ps_calib <- function(x, y, ...) {
  warn_class_downgrade(c("psw", "ps_calib"))
  double()
}

#' @export
vec_ptype2.ps_calib.psw <- function(x, y, ...) {
  warn_class_downgrade(c("ps_calib", "psw"))
  double()
}

#' @export
vec_ptype2.ps_calib.ps_trim <- function(x, y, ...) {
  warn_class_downgrade(c("ps_calib", "ps_trim"))
  double()
}

#' @export
vec_ptype2.ps_trim.ps_calib <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trim", "ps_calib"))
  double()
}

#' @export
vec_ptype2.ps_calib.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade(c("ps_calib", "ps_trunc"))
  double()
}

#' @export
vec_ptype2.ps_trunc.ps_calib <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trunc", "ps_calib"))
  double()
}

#' @export
vec_math.ps_calib <- function(.fn, .x, ...) {
  vec_math_base(.fn, vec_data(.x), ...)
}

#' @export
Summary.ps_calib <- function(..., na.rm = FALSE) {
  .fn <- .Generic
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call(.fn, c(numeric_args, list(na.rm = na.rm)))
}

#' @export
min.ps_calib <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("min", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
max.ps_calib <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("max", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
range.ps_calib <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("range", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
median.ps_calib <- function(x, na.rm = FALSE, ...) {
  median(vec_data(x), na.rm = na.rm, ...)
}

#' @export
`[.ps_calib` <- function(x, i, ...) {
  # If i is missing, just call NextMethod
  if (missing(i)) {
    return(NextMethod())
  }

  # Get the subset of data using NextMethod to handle vctrs subsetting
  result <- NextMethod()

  # Preserve metadata
  attr(result, "ps_calib_meta") <- ps_calib_meta(x)

  result
}

#' @export
print.ps_calib <- function(x, ..., n = NULL) {
  meta <- ps_calib_meta(x)
  n_vals <- length(x)

  # Create header
  cat(sprintf(
    "<ps_calib[%d]; method=%s%s>\n",
    n_vals,
    meta$method,
    if (isTRUE(meta$smooth)) " (smooth)" else ""
  ))

  # Determine how many values to show
  if (is.null(n)) {
    n_print <- getOption("propensity.print_max", default = 10)
  } else {
    n_print <- n
  }

  # Show values
  if (is.infinite(n_print) || n_print >= n_vals) {
    print(vec_data(x))
  } else {
    n_show <- min(n_print, n_vals)
    print(vec_data(x)[seq_len(n_show)])

    if (n_vals > n_show) {
      cat("# ... with", n_vals - n_show, "more values\n")
    }
  }

  invisible(x)
}

#' @export
print.ps_calib_matrix <- function(x, ..., n = NULL) {
  meta <- ps_calib_meta(x)
  n_rows <- nrow(x)
  k <- ncol(x)

  # Create header
  cat(sprintf(
    "<ps_calib_matrix[%d x %d]; method=%s%s>\n",
    n_rows,
    k,
    meta$method,
    if (isTRUE(meta$smooth)) " (smooth)" else ""
  ))

  # Determine how many rows to show
  if (is.null(n)) {
    n_print <- getOption("propensity.print_max", default = 10)
  } else {
    n_print <- n
  }

  # Show all rows if n is Inf or very large
  if (is.infinite(n_print) || n_print >= n_rows) {
    print(unclass(x))
  } else {
    # Show first n_print rows
    n_show <- min(n_print, n_rows)
    x_sub <- x[seq_len(n_show), , drop = FALSE]
    print(unclass(x_sub))

    if (n_rows > n_show) {
      cat("# ... with", n_rows - n_show, "more rows\n")
    }
  }

  invisible(x)
}
