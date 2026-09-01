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
#' @param .propensity A numeric vector of propensity scores in `[0, 1]`. Unlike
#'   the rest of the package, calibration accepts scores of exactly 0 and
#'   exactly 1, since repairing scores a model pushed to the ends of the
#'   interval is part of what calibration is for. The logistic calibration curve
#'   maps them back inside the interval; isotonic calibration can return a score
#'   at an endpoint when its pooled block is pure, and such scores are rejected
#'   by the weight functions. Must not already be calibrated.
#'
#'   A data frame of predicted probabilities, one column per exposure level, is
#'   reduced to a single column the way [ps_trim()] and [ps_trunc()] reduce one:
#'   the second column of a pair, or the only column of a frame that has one.
#'   The column read is announced. A frame holding a `.pred_class` column, which
#'   a fitted tidymodels classification model returns when no prediction type is
#'   named, carries predicted levels rather than probabilities and is refused
#'   with an error of class `propensity_df_class_column_error`.
#' @param .exposure A binary vector of observed treatment assignments, the same
#'   length as `.propensity`.
#' @param method Calibration method. One of:
#'   \describe{
#'     \item{`"logistic"`}{(Default) Logistic calibration, also called Platt
#'       scaling. Fits a logistic regression of `.exposure` on `.propensity`,
#'       yielding a smooth, parametric correction. Works well with small
#'       samples and when the bias in `.propensity` is approximately monotone.}
#'     \item{`"isoreg"`}{Isotonic regression. Fits a non-parametric,
#'       monotone step function. More flexible than logistic calibration
#'       because it makes no distributional assumption, but needs larger
#'       samples for stable estimates.}
#'   }
#' @param smooth Logical. When `method = "logistic"`, controls the form of the
#'   calibration model. If `TRUE` (default), fits a GAM with a spline on
#'   `.propensity` via [mgcv::gam()]; if `FALSE`, fits a simple logistic
#'   regression via [stats::glm()]. Ignored when `method = "isoreg"`. A spline
#'   needs enough distinct scores to place its knots, so with fewer than 10
#'   distinct values of `.propensity` among the observations the fit reads,
#'   those with both a score and an exposure recorded, the fit falls back to
#'   logistic regression without a spline. The fallback is announced and
#'   recorded in the returned metadata.
#' @param .focal_level The value of `.exposure` representing the focal group
#'   (typically the treated group). Every binary coding honors it: 0/1 numeric,
#'   logical, two-level factor, and two-level character exposures are all coded
#'   with the named level as focal, and a level the exposure never takes is an
#'   error. With no level named, the exposure defaults to its higher level,
#'   which is `1` for a 0/1 exposure, `TRUE` for a logical one, and the second
#'   of the two levels a factor or character exposure takes. Levels a factor
#'   declares but never takes are not candidates. Naming any other level
#'   reverses the coding, so `.propensity` must then hold the probability of
#'   the named level.
#' @param .reference_level The value of `.exposure` representing the reference
#'   group (typically the control group). Naming it makes the exposure's other
#'   level focal, with the same consequence for `.propensity`, and a level the
#'   exposure never takes is an error. Automatically detected if not supplied.
#' @param .treated `r lifecycle::badge("deprecated")` Use `.focal_level` instead.
#' @param .untreated `r lifecycle::badge("deprecated")` Use `.reference_level` instead.
#' @param ... Additional arguments passed to methods.
#' @param ps `r lifecycle::badge("deprecated")` Use `.propensity` instead. A
#'   call that names `ps` must name the arguments after it as well, since a
#'   positional argument binds to `.propensity`.
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
#' **Fitted models:**
#' `.propensity` can be the fitted propensity score model instead of the scores
#' it reports, and the exposure is then read off the model as well unless
#' `.exposure` is supplied. The variable read is announced;
#' `options(propensity.quiet = TRUE)` silences the announcement. A binomial
#' [stats::glm()] and a two-level `nnet::multinom()` report the one score per
#' unit calibration reads. A `nnet::multinom()` of three or more levels reports
#' one column per level and no single score to calibrate, so it is refused with
#' an error of class `propensity_model_family_error`; calibrate the columns one
#' at a time against the exposure indicator each of them belongs to.
#'
#' **Missing values:**
#' A unit with a missing exposure tells the calibration model nothing, so it is
#' dropped from the fit. What happens to that unit afterwards depends on the
#' method. Logistic calibration, smoothed or not, fits a curve from propensity
#' score to calibrated score and reads every unit off it, so a unit with a
#' missing exposure but an observed score is still calibrated. Isotonic
#' calibration fits separately within each exposure group and reads each unit
#' from the fit for the group it belongs to; a unit with a missing exposure
#' belongs to neither group and is returned as `NA`. Under both methods a
#' missing propensity score yields a missing calibrated score.
#'
#' The calibrated scores are returned as a `ps_calib` object, which can be
#' passed directly to weight functions such as [wt_ate()].
#'
#' @return A `ps_calib` vector the same length as `.propensity`. The attribute
#'   `ps_calib_meta` stores calibration metadata (method and whether
#'   smoothing was applied). Use [is_ps_calibrated()] to test whether an
#'   object has been calibrated.
#'
#' @references
#' Platt, J. (1999). Probabilistic outputs for support vector machines and
#' comparisons to regularized likelihood methods. *Advances in Large Margin
#' Classifiers*, 61--74.
#'
#' Zadrozny, B., & Elkan, C. (2002). Transforming classifier scores into
#' accurate multiclass probability estimates. *Proceedings of the Eighth ACM
#' SIGKDD International Conference on Knowledge Discovery and Data Mining*,
#' 694--699. \doi{10.1145/775047.775151}
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
#' # Calibrate the scores a fitted model reports, reading the exposure off it
#' x <- rnorm(200)
#' fit <- glm(exposure ~ x, family = binomial())
#' ps_calibrate(fit, smooth = FALSE)
#'
#' if (rlang::is_installed("nnet")) {
#'   exposure_factor <- factor(exposure)
#'   multinomial_fit <- nnet::multinom(exposure_factor ~ x, trace = FALSE)
#'   ps_calibrate(multinomial_fit, smooth = FALSE)
#' }
#'
#' @importFrom stats glm binomial predict
#' @export
ps_calibrate <- function(
  .propensity,
  .exposure,
  method = c("logistic", "isoreg"),
  smooth = TRUE,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated()
) {
  .propensity <- handle_propensity_deprecation(
    rlang::maybe_missing(.propensity),
    ps,
    "ps_calibrate"
  )

  UseMethod("ps_calibrate", .propensity)
}

#' @export
ps_calibrate.default <- function(
  .propensity,
  .exposure,
  method = c("logistic", "isoreg"),
  smooth = TRUE,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  # Two frames arrive here because two condition systems read them.
  # `user_env` is the frame lifecycle reports a deprecation from, which decides
  # whether the reader is told to change their own call or to report an issue.
  # `call` is the frame every other condition is attributed to, which rlang
  # reads to name the function in the report. A route that reaches this method
  # by a call rather than by dispatch has to supply both.
  user_env = rlang::caller_env(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)
  check_ps_method(.propensity, call = call)
  method <- rlang::arg_match(method, error_call = call)

  if (rlang::is_missing(.exposure)) {
    abort(
      c(
        "{.arg .exposure} must be supplied.",
        i = "Calibration reads the propensity scores against the exposure
             they are the probability of."
      ),
      error_class = "propensity_missing_arg_error",
      call = call
    )
  }

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "ps_calibrate",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)
  # Check that the propensity scores are numeric and in valid range
  # A one-dimensional array holds one score per observation, which is the shape
  # calibration reads, so its dimension is dropped the way the weight functions
  # already drop it. Two or more dimensions hold something else.
  if (length(dim(.propensity)) == 1L) {
    dim(.propensity) <- NULL
  }

  if (!is.null(dim(.propensity))) {
    abort(
      c(
        "{.arg .propensity} must be a numeric vector.",
        x = "It is {.cls {class(.propensity)[[1]]}} with \\
        {length(dim(.propensity))} dimension{?s}.",
        i = "Calibration reads one propensity score per observation. Calibrate
             the columns of a matrix of scores one at a time."
      ),
      error_class = "propensity_type_error",
      call = call
    )
  }

  if (!is.numeric(.propensity)) {
    abort(
      "{.arg .propensity} must be a numeric vector.",
      error_class = "propensity_type_error",
      call = call
    )
  }

  if (is_ps_calibrated(.propensity)) {
    abort(
      "{.arg .propensity} is already calibrated. Cannot calibrate already calibrated propensity scores.",
      error_class = "propensity_already_calibrated_error",
      call = call
    )
  }

  # Calibration reads the values the scores hold, not the class holding them,
  # and trimmed or truncated scores reach the numeric routes through their
  # common type with a number, which announces a conversion the caller never
  # asked for. Reading them once here keeps the rest of the function numeric.
  ps_numeric <- as.numeric(.propensity)

  # The range accepted here is the closed interval, not the open one
  # `check_ps_range()` enforces everywhere else. Trimming and weighting divide
  # by a score, so a 0 or a 1 has no weight to give and is refused; calibration
  # repairs scores a model pushed to the ends of the interval, so refusing them
  # would turn away the very scores calibration exists to fix. What comes back
  # is inside the interval under logistic calibration, whose curve never
  # reaches either end. Isotonic calibration can return a score at an endpoint
  # when its pooled block is pure, and such scores are rejected by the weight
  # functions in turn.
  if (any(ps_numeric < 0 | ps_numeric > 1, na.rm = TRUE)) {
    abort(
      "{.arg .propensity} values must be between 0 and 1.",
      error_class = "propensity_range_error",
      call = call
    )
  }

  # Transform treatment to binary if needed
  .exposure <- transform_exposure_binary(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    call = call
  )

  if (length(.propensity) != length(.exposure)) {
    abort(
      "Propensity score vector {.arg .propensity} must be the same length as \\
      {.arg .exposure}.",
      error_class = "propensity_length_error",
      call = call
    )
  }

  # ps_calib objects should not have weight-specific attributes like estimand
  # Only store calibration-specific metadata

  # Handle NA values
  na_idx <- is.na(ps_numeric) | is.na(.exposure)

  # Perform calibration based on method
  if (method == "logistic") {
    # The units the fit reads: those with both a score and an exposure.
    calib_data <- data.frame(
      treat = .exposure[!na_idx],
      ps = ps_numeric[!na_idx]
    )

    if (smooth) {
      # Whether a spline can be placed is settled before mgcv is asked for.
      # Deciding the other way round refuses the call over a package the
      # calibration it goes on to perform takes no part in.
      n_unique <- length(unique(calib_data$ps))
      if (n_unique < 10) {
        # Fall back to simple logistic regression if too few unique values.
        # The caller asked for a spline and is getting a straight line, so the
        # substitution is announced rather than left to the metadata.
        smooth <- FALSE
        alert_info(
          "{n_unique} distinct propensity score{?s} {?is/are} too few to place the knots of a spline, so calibration falls back to logistic regression without one."
        )
      }
    }

    if (smooth) {
      rlang::check_installed(
        "mgcv",
        reason = "for smooth calibration using GAM (Generalized Additive Model)."
      )

      # Fit GAM with smooth spline: treat ~ s(ps)
      calib_model <- mgcv::gam(
        treat ~ s(ps),
        data = calib_data,
        family = "binomial"
      )
    } else {
      calib_model <- stats::glm(
        treat ~ ps,
        data = calib_data,
        family = stats::binomial()
      )
    }

    if (isFALSE(calib_model$converged)) {
      warn(
        "Calibration model did not converge",
        warning_class = "propensity_convergence_warning",
        call = call
      )
    }

    # Both logistic fits read the units with an observed exposure, and both
    # produce a curve from propensity score to calibrated score, so both are
    # read over every unit. A unit dropped from the fit for want of an exposure
    # still has a score the curve maps; only a missing score has no answer.
    calib_ps <- as.numeric(predict(
      calib_model,
      newdata = data.frame(ps = ps_numeric),
      type = "response"
    ))
  } else if (method == "isoreg") {
    # Smoothing is not applicable to isotonic regression
    smooth <- FALSE

    # Two-step isotonic regression calibration following van der Laan et al.
    # (2024, arXiv:2411.06342): fit separately for treated and control groups,
    # then combine. This matches WeightIt::calibrate(method = "isoreg").
    ps_valid <- ps_numeric[!na_idx]
    treat_valid <- .exposure[!na_idx]

    # Calibrate for controls: P(treat=0|ps) via isotonic regression on 1-ps
    p0 <- 1 - pava_weighted(1 - ps_valid, 1 - treat_valid)
    # Calibrate for treated: P(treat=1|ps) via isotonic regression on ps
    p1 <- pava_weighted(ps_valid, treat_valid)

    is_trt <- treat_valid == 1

    # Combine: controls use p0, treated use p1
    calib_ps_valid <- p0
    calib_ps_valid[is_trt] <- p1[is_trt]

    # Map back to full vector with NAs
    calib_ps <- numeric(length(.propensity))
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

# A frame of predicted probabilities holds one column per exposure level, and
# calibration reads one score per unit. The column is chosen the way trimming
# and truncation choose it, so the same frame reaches the same score whichever
# of the three it is handed to.
#' @export
ps_calibrate.data.frame <- function(
  .propensity,
  .exposure,
  method = c("logistic", "isoreg"),
  smooth = TRUE,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)
  check_predicted_class_column(.propensity, call = call)

  ps_vec <- binary_ps_column(.propensity, "ps_calibrate")

  # The default method is reached here by a call rather than by dispatch, so it
  # is handed a frame to report against. Left to its own, a refusal on this
  # route would name the method the caller never wrote.
  ps_calibrate.default(
    .propensity = ps_vec,
    .exposure = rlang::maybe_missing(.exposure),
    method = method,
    smooth = smooth,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .treated = .treated,
    .untreated = .untreated,
    user_env = rlang::caller_env(),
    call = call
  )
}

# The fitted propensity score models calibration reads, registered for the same
# classes the weight functions read: a `glm`, whose binomial families fit the
# probability of a binary exposure, and a `multinom`. A `lm` is not among them,
# its fitted values being conditional means rather than probabilities, and it
# reaches the default method, which reports that it has no scores to calibrate.
#' @export
ps_calibrate.glm <- function(
  .propensity,
  .exposure,
  method = c("logistic", "isoreg"),
  smooth = TRUE,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)

  calibrate_from_model(
    .propensity,
    .exposure = rlang::maybe_missing(.exposure, NULL),
    method = method,
    smooth = smooth,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    call = call,
    user_env = rlang::caller_env()
  )
}

# A `multinom` is `c("multinom", "nnet")` and inherits from neither `glm` nor
# `lm`, so it reaches a method of its own. What it was fit to is checked before
# anything reads it, so that a fit the route cannot read is reported as such
# rather than as whatever its fitted values are made of.
#' @export
ps_calibrate.multinom <- function(
  .propensity,
  .exposure,
  method = c("logistic", "isoreg"),
  smooth = TRUE,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)
  check_multinom_response(.propensity, call = call)

  calibrate_from_model(
    .propensity,
    .exposure = rlang::maybe_missing(.exposure, NULL),
    method = method,
    smooth = smooth,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    call = call,
    user_env = rlang::caller_env()
  )
}

# Calibration reads one propensity score for each unit and the exposure those
# scores are calibrated against, so a fit is always read for its probability of
# a binary exposure. A fit of three or more levels reports a column for each of
# them and no single probability to calibrate, which is refused here rather than
# on the binary reading, whose remedy is about the weights a categorical
# exposure takes and so has nothing to offer a caller calibrating scores.
#
# The default method is reached by a call rather than by dispatch, so it is
# handed the frame the route was entered on. Left to its own, a refusal here
# would name a method the caller never wrote.
calibrate_from_model <- function(
  model,
  .exposure,
  method,
  smooth,
  .focal_level,
  .reference_level,
  .treated,
  .untreated,
  call,
  user_env
) {
  if (model_fits_levels(model)) {
    abort_calibrate_levels(model, call = call)
  }

  args <- prepare_model_ps(
    model,
    .exposure,
    needs_exposure = TRUE,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "ps_calibrate",
    call = call,
    user_env = user_env
  )

  ps_calibrate.default(
    .propensity = args$propensity,
    .exposure = args$exposure,
    method = method,
    smooth = smooth,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    user_env = user_env,
    call = call
  )
}

# A fit of three or more levels holds one probability for each level and no
# single score to calibrate. The remedy is calibration's own: each column is the
# probability of one level, and calibrating it against the indicator for that
# level is the reading the model has to offer.
abort_calibrate_levels <- function(model, call = rlang::caller_env()) {
  fit_levels <- model_levels(model)
  n_levels <- length(fit_levels)

  abort(
    c(
      "Calibration reads one propensity score for each unit.",
      x = "{.arg .propensity} was fit to {n_levels} levels
           ({.val {fit_levels}}), so it fits one probability for each level rather
           than one for each unit.",
      i = "Calibrate the columns of {.code fitted(fit)} one at a time against
           the indicator for each level."
    ),
    error_class = "propensity_model_family_error",
    call = call
  )
}

#' Constructor for ps_calib objects
#' @noRd
new_ps_calib <- function(x, ps_calib_meta = list()) {
  vec_assert(x, ptype = double())
  new_vctr(
    x,
    ps_calib_meta = ps_calib_meta,
    class = "ps_calib",
    inherit_base_type = TRUE
  )
}

#' Extract metadata from ps_calib object
#' @noRd
ps_calib_meta <- function(x) {
  attr(x, "ps_calib_meta")
}

#' Check if propensity scores are calibrated
#'
#' `is_ps_calibrated()` tests whether `x` was calibrated, rather than whether it
#' is still a `ps_calib` object: it is `TRUE` for a calibrated propensity score,
#' for weights derived from one, and for scores that were calibrated before
#' [ps_trim()] or [ps_trunc()] bounded them, each of which records the
#' calibration in its own metadata.
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
#' calibrated <- ps_calibrate(ps, exposure, smooth = FALSE)
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

# Trimming and truncation read the scores a calibrated object holds and drop
# the class, because what comes back is a trimmed or bounded score rather than
# a calibrated one. The calibration is recorded in their own metadata, and this
# predicate has always answered from a record rather than from a class, as it
# does for the weights built from calibrated scores.
#' @export
is_ps_calibrated.ps_trim <- function(x) {
  isTRUE(ps_trim_meta(x)$calibrated)
}

#' @export
is_ps_calibrated.ps_trunc <- function(x) {
  isTRUE(ps_trunc_meta(x)$calibrated)
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

# A comparison asks about the values a vector holds, not about the type holding
# them, so it reads them directly instead of routing through the common type of
# the two sides, which for a calibrated vector and a number is the numeric
# downgrade. Recycling still goes through vctrs, so sizes with no common size
# have no answer rather than the one base R recycling would invent. The refusal
# is reported against the comparison the user wrote rather than against this
# helper, which they have no way to reach.
ps_calib_compare <- function(e1, e2, call = rlang::caller_env()) {
  args <- vctrs::vec_recycle_common(e1, e2, .call = call)

  if (inherits(args[[1]], "ps_calib")) {
    args[[1]] <- vctrs::vec_data(args[[1]])
  }
  if (inherits(args[[2]], "ps_calib")) {
    args[[2]] <- vctrs::vec_data(args[[2]])
  }

  args
}

#' @export
`<.ps_calib` <- function(e1, e2) {
  args <- ps_calib_compare(e1, e2)
  args[[1]] < args[[2]]
}

#' @export
`<=.ps_calib` <- function(e1, e2) {
  args <- ps_calib_compare(e1, e2)
  args[[1]] <= args[[2]]
}

#' @export
`==.ps_calib` <- function(e1, e2) {
  args <- ps_calib_compare(e1, e2)
  args[[1]] == args[[2]]
}

#' @export
`>.ps_calib` <- function(e1, e2) {
  args <- ps_calib_compare(e1, e2)
  args[[1]] > args[[2]]
}

#' @export
`>=.ps_calib` <- function(e1, e2) {
  args <- ps_calib_compare(e1, e2)
  args[[1]] >= args[[2]]
}

#' @export
`!=.ps_calib` <- function(e1, e2) {
  args <- ps_calib_compare(e1, e2)
  args[[1]] != args[[2]]
}

# How a calibration is described: the curve it was fit with and whether that
# curve was smoothed. A `ps_calib` carries no positional record, so this is the
# whole of its metadata, and it is read through one helper so that the combine
# and the cast compare the same thing.
calib_parameters <- function(meta) {
  fields <- c("method", "smooth")

  rlang::set_names(lapply(fields, function(field) meta[[field]]), fields)
}

#' @export
vec_ptype2.ps_calib.ps_calib <- function(x, y, ...) {
  x_meta <- ps_calib_meta(x)
  y_meta <- ps_calib_meta(y)

  # Check if calibration parameters match
  if (!identical(calib_parameters(x_meta), calib_parameters(y_meta))) {
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

# A cast returns the values it was handed in the type it was handed, and a
# `ps_calib`'s type is the whole description of the calibration. That is the
# comparison `vec_ptype2()` already makes when it refuses to find a common type,
# so the cast makes it too.
#
# `vec_ptype_full()` names the method but not whether the fit was smoothed, so
# two types that differ only in that render identically and the refusal would
# read as a type that cannot be converted to itself. What disagrees is named
# alongside them, the way the combine names it.
#' @export
vec_cast.ps_calib.ps_calib <- function(x, to, ...) {
  # Check if metadata matches
  x_meta <- ps_calib_meta(x)
  to_meta <- ps_calib_meta(to)

  if (!identical(calib_parameters(x_meta), calib_parameters(to_meta))) {
    vctrs::stop_incompatible_cast(
      x,
      to,
      x_arg = "",
      to_arg = "",
      details = "different calibration parameters"
    )
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

# A cast returns the values it was handed in the type it was handed, and how the
# calibration was performed is part of that type, so the description comes from
# the target rather than being invented here.
#' @export
vec_cast.ps_calib.double <- function(x, to, ...) {
  new_ps_calib(vec_cast(x, double()), ps_calib_meta = ps_calib_meta(to))
}

# Calibrated propensity scores lie strictly between 0 and 1, so meeting an
# integer in the integers would round every one of them away. The common type is
# the one that holds both sets of values.
#' @export
vec_ptype2.ps_calib.integer <- function(x, y, ...) {
  warn_class_downgrade("ps_calib")
  double()
}

#' @export
vec_ptype2.integer.ps_calib <- function(x, y, ...) {
  warn_class_downgrade("ps_calib")
  double()
}

#' @export
vec_cast.integer.ps_calib <- function(x, to, ...) {
  # A propensity score has no integer to be, so vctrs' own check reports the
  # loss rather than silently rounding it away.
  vec_cast(vec_data(x), integer(), x_arg = "ps_calib")
}

#' @export
vec_cast.ps_calib.integer <- function(x, to, ...) {
  new_ps_calib(as.double(x), ps_calib_meta = ps_calib_meta(to))
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
