#' Weighted pool-adjacent-violators algorithm (PAVA)
#'
#' Implements isotonic regression with optional observation weights. Unlike
#' [stats::isoreg()], this does not group tied x-values before fitting, so
#' tied inputs can receive different fitted values.
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values.
#' @param w Numeric vector of weights (default: equal weights).
#' @return Numeric vector of fitted values in the original order of `x`.
#' @noRd
pava_weighted <- function(x, y, w = rep(1, length(x))) {
  n <- length(x)
  if (n <= 1L) {
    return(y)
  }

  # Order by x
  ord <- order(x)
  y_ord <- y[ord]
  w_ord <- w[ord]

  # Initialize each observation as its own block
  # Store value (weighted mean) and weight (sum of weights) per block
  block_val <- y_ord
  block_wt <- w_ord
  # block_end[i] gives the last original index belonging to block i

  block_end <- seq_len(n)
  n_blocks <- n

  i <- 1L
  while (i < n_blocks) {
    if (block_val[i] > block_val[i + 1L]) {
      # Merge blocks i and i+1
      new_wt <- block_wt[i] + block_wt[i + 1L]
      new_val <- (block_val[i] *
        block_wt[i] +
        block_val[i + 1L] * block_wt[i + 1L]) /
        new_wt
      block_val[i] <- new_val
      block_wt[i] <- new_wt
      block_end[i] <- block_end[i + 1L]

      # Remove block i+1
      if (i + 1L < n_blocks) {
        idx_keep <- seq_len(n_blocks)[-c(i + 1L)]
        block_val <- block_val[idx_keep]
        block_wt <- block_wt[idx_keep]
        block_end <- block_end[idx_keep]
      } else {
        block_val <- block_val[seq_len(n_blocks - 1L)]
        block_wt <- block_wt[seq_len(n_blocks - 1L)]
        block_end <- block_end[seq_len(n_blocks - 1L)]
      }
      n_blocks <- n_blocks - 1L

      # Step back to recheck
      if (i > 1L) i <- i - 1L
    } else {
      i <- i + 1L
    }
  }

  # Reconstruct fitted values from blocks
  fitted <- numeric(n)
  start <- 1L
  for (j in seq_len(n_blocks)) {
    fitted[start:block_end[j]] <- block_val[j]
    start <- block_end[j] + 1L
  }

  # Return in original order
  fitted[order(ord)]
}

#' Calibrate propensity scores
#'
#' @description
#' `ps_calibrate()` adjusts estimated propensity scores so they better reflect
#' true treatment probabilities. This can improve the accuracy of inverse
#' probability weights derived from a misspecified propensity score model.
#'
#' @param ps A numeric vector of propensity scores between 0 and 1. Must not
#'   already be calibrated.
#' @param .exposure A binary vector of observed treatment assignments, the same
#'   length as `ps`.
#' @param method Calibration method. One of:
#'   \describe{
#'     \item{`"logistic"`}{(Default) Logistic calibration, also called Platt
#'       scaling. Fits a logistic regression of `.exposure` on `ps`, yielding
#'       a smooth, parametric correction. Works well with small samples and
#'       when the bias in `ps` is approximately monotone.}
#'     \item{`"isoreg"`}{Isotonic regression. Fits a non-parametric,
#'       monotone step function. More flexible than logistic calibration
#'       because it makes no distributional assumption, but needs larger
#'       samples for stable estimates.}
#'   }
#' @param smooth Logical. When `method = "logistic"`, controls the form of the
#'   calibration model. If `TRUE` (default), fits a GAM with a spline on
#'   `ps` via [mgcv::gam()]; if `FALSE`, fits a simple logistic regression
#'   via [stats::glm()]. Ignored when `method = "isoreg"`.
#' @param .focal_level The value of `.exposure` representing the focal group
#'   (typically the treated group). If `NULL` (default), coding is determined
#'   automatically.
#' @param .reference_level The value of `.exposure` representing the reference
#'   group (typically the control group). If `NULL` (default), coding is
#'   determined automatically.
#' @param .treated `r lifecycle::badge("deprecated")` Use `.focal_level` instead.
#' @param .untreated `r lifecycle::badge("deprecated")` Use `.reference_level` instead.
#'
#' @details
#' Calibration is useful when the propensity score model is correctly
#' specified in terms of variable selection but produces probabilities that
#' are systematically too high or too low. Unlike [ps_trim()] and
#' [ps_trunc()], which handle extreme scores by removing or bounding them,
#' calibration reshapes the entire distribution of scores.
#'
#' **Choosing a method:**
#' - Use `"logistic"` (the default) as a first choice. It is stable and
#'   fast, and the `smooth = TRUE` option adds flexibility via a spline.
#' - Use `"isoreg"` when you suspect a non-smooth or irregular relationship
#'   between estimated and true probabilities and have a sufficiently large
#'   sample.
#'
#' The calibrated scores are returned as a `ps_calib` object, which can be
#' passed directly to weight functions such as [wt_ate()].
#'
#' @return A `ps_calib` vector the same length as `ps`. The attribute
#'   `ps_calib_meta` stores calibration metadata (method and whether
#'   smoothing was applied). Use [is_ps_calibrated()] to test whether an
#'   object has been calibrated.
#'
#' @seealso [is_ps_calibrated()] to test for calibrated scores;
#'   [ps_trim()] and [ps_trunc()] for alternative approaches to extreme
#'   propensity scores; [wt_ate()] and other weight functions that accept
#'   `ps_calib` objects.
#'
#' @examples
#' # Simulate data
#' set.seed(42)
#' ps <- runif(200)
#' exposure <- rbinom(200, 1, ps)
#'
#' # Logistic calibration without smoothing (simple Platt scaling)
#' cal <- ps_calibrate(ps, exposure, smooth = FALSE)
#' cal
#'
#' # Use calibrated scores to calculate weights
#' wt_ate(cal, exposure)
#'
#' # Isotonic regression calibration
#' cal_iso <- ps_calibrate(ps, exposure, method = "isoreg")
#'
#' if (rlang::is_installed("mgcv")) {
#'   # Logistic calibration with spline smoothing (default)
#'   cal_smooth <- ps_calibrate(ps, exposure)
#' }
#'
#' @importFrom stats glm fitted binomial
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

    # Two-step isotonic regression calibration following van der Laan et al.
    # (2024, arXiv:2411.06342): fit separately for treated and control groups,
    # then combine. This matches WeightIt::calibrate(method = "isoreg").
    ps_valid <- ps[!na_idx]
    treat_valid <- .exposure[!na_idx]

    # Calibrate for controls: P(treat=0|ps) via isotonic regression on 1-ps
    p0 <- 1 - pava_weighted(1 - ps_valid, 1 - treat_valid)
    # Calibrate for treated: P(treat=1|ps) via isotonic regression on ps
    p1 <- pava_weighted(ps_valid, treat_valid)

    # Squish to prevent extrapolation beyond observed range
    is_ctrl <- treat_valid == 0
    is_trt <- treat_valid == 1
    if (any(is_ctrl)) {
      p0 <- pmax(p0, min(p0[is_ctrl]))
    }
    if (any(is_trt)) {
      p1 <- pmax(p1, min(p1[is_trt]))
    }

    # Combine: controls use p0, treated use p1
    calib_ps_valid <- p0
    calib_ps_valid[is_trt] <- p1[is_trt]

    # Map back to full vector with NAs
    calib_ps <- numeric(length(ps))
    calib_ps[!na_idx] <- calib_ps_valid
    if (any(na_idx)) {
      calib_ps[na_idx] <- NA
    }

    # Ensure calibrated values are in [0, 1]
    calib_ps[!na_idx] <- pmax(0, pmin(1, calib_ps[!na_idx]))
  }

  # Return calibrated scores as a ps_calib object
  new_ps_calib(
    x = calib_ps,
    ps_calib_meta = list(
      method = method,
      smooth = smooth
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

#' Check if propensity scores are calibrated
#'
#' `is_ps_calibrated()` tests whether `x` is a calibrated propensity score
#' object (class `ps_calib`) or a `psw` object derived from calibrated scores.
#'
#' @param x An object to test.
#' @return A single `TRUE` or `FALSE`.
#'
#' @seealso [ps_calibrate()] to calibrate propensity scores.
#'
#' @examples
#' ps <- runif(100)
#' exposure <- rbinom(100, 1, ps)
#'
#' is_ps_calibrated(ps)
#'
#' calibrated <- ps_calibrate(ps, exposure)
#' is_ps_calibrated(calibrated)
#'
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
