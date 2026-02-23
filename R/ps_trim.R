#' Trim Propensity Scores
#'
#' @description
#' Trim observations with extreme propensity scores by replacing them with `NA`,
#' effectively removing those units from downstream analyses. The returned object
#' has the same length (or dimensions) as the input, with trimmed entries set to
#' `NA`. After trimming, refit the propensity score model on the retained
#' observations with [ps_refit()].
#'
#' @param ps A numeric vector of propensity scores in (0, 1) for binary
#'   exposures, or a matrix / data frame where each column gives the propensity
#'   score for one level of a categorical exposure.
#' @param method Trimming method. One of:
#'
#'   * **`"ps"`** (default): Fixed threshold. Observations with propensity
#'     scores outside `[lower, upper]` are trimmed. For categorical exposures,
#'     observations where *any* column falls below `lower` (the symmetric
#'     threshold delta) are trimmed.
#'   * **`"adaptive"`**: Data-driven threshold that minimizes the asymptotic
#'     variance of the IPW estimator (Crump et al., 2009). The `lower` and
#'     `upper` arguments are ignored.
#'   * **`"pctl"`**: Quantile-based. Observations outside the `[lower, upper]`
#'     quantiles of the propensity score distribution are trimmed. Defaults:
#'     `lower = 0.05`, `upper = 0.95`.
#'   * **`"pref"`**: Preference score trimming. Transforms propensity scores
#'     to the preference scale (Walker et al., 2013) and trims outside
#'     `[lower, upper]`. Requires `.exposure`. Binary exposures only. Defaults:
#'     `lower = 0.3`, `upper = 0.7`.
#'   * **`"cr"`**: Common range (clinical equipoise). Trims to the overlap
#'     region of the propensity score distributions across exposure groups.
#'     Requires `.exposure`. Binary exposures only. The `lower` and `upper`
#'     arguments are ignored.
#'   * **`"optimal"`**: Multi-category optimal trimming (Yang et al., 2016).
#'     Categorical exposures only. Requires `.exposure`.
#'
#'   For categorical exposures, only `"ps"` and `"optimal"` are supported.
#' @param lower,upper Numeric thresholds whose interpretation depends on
#'   `method`:
#'
#'   * `"ps"`: absolute propensity score bounds (defaults: 0.1, 0.9). For
#'     categorical exposures, only `lower` is used as the symmetric threshold.
#'   * `"pctl"`: quantile probabilities (defaults: 0.05, 0.95).
#'   * `"pref"`: preference score bounds (defaults: 0.3, 0.7).
#'   * `"adaptive"`, `"cr"`, `"optimal"`: ignored (thresholds are data-driven).
#' @param .exposure An exposure variable. Required for `"pref"`, `"cr"` (binary
#'   vector), and `"optimal"` (factor or character). Not required for other
#'   methods.
#' @inheritParams wt_ate
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' ## How trimming works
#'
#' Trimming identifies observations with extreme (near 0 or 1) propensity
#' scores and sets them to `NA`. These observations are excluded from
#' subsequent weight calculations and effect estimation. The goal is to
#' remove units that lack sufficient overlap between exposure groups, which
#' would otherwise receive extreme weights and destabilize estimates.
#'
#' ## Choosing a method
#'
#' * Use `"ps"` when you have a specific threshold in mind or want a simple
#'   default.
#' * Use `"adaptive"` for a principled, data-driven cutoff that targets
#'   variance reduction.
#' * Use `"pctl"` to trim a fixed percentage of extreme values from each tail.
#' * Use `"pref"` when you want to restrict to the region of clinical
#'   equipoise based on the preference score.
#' * Use `"cr"` to restrict to the common support region where both exposure
#'   groups have observed propensity scores.
#' * Use `"optimal"` for multi-category (3+) exposures; this is the only
#'   data-driven method available for categorical treatments.
#'
#' ## Typical workflow
#'
#' 1. Fit a propensity score model
#' 2. Apply `ps_trim()` to flag extreme values
#' 3. Call [ps_refit()] to re-estimate propensity scores on the retained sample
#' 4. Compute weights with [wt_ate()] or another weight function
#'
#' ## Object behavior
#'
#' Arithmetic operations on `ps_trim` objects return plain numeric vectors,
#' since transformed propensity scores (e.g., `1/ps`) are no longer propensity
#' scores. Trimmed values propagate as `NA` in calculations; use `na.rm = TRUE`
#' where appropriate.
#'
#' When combining `ps_trim` objects with [c()], trimming parameters must match.
#' Mismatched parameters trigger a warning and return a numeric vector.
#'
#' Use [ps_trim_meta()] to inspect the trimming metadata, including the method,
#' cutoffs, and which observations were retained or trimmed.
#'
#' @return A **`ps_trim`** object (a numeric vector with class `"ps_trim"`, or a
#'   matrix with class `"ps_trim_matrix"`). Trimmed observations are `NA`.
#'   Metadata is stored in the `"ps_trim_meta"` attribute and can be accessed
#'   with [ps_trim_meta()]. Key fields include:
#'
#'   * `method`: the trimming method used
#'   * `keep_idx`: integer indices of retained observations
#'   * `trimmed_idx`: integer indices of trimmed (NA) observations
#'   * Method-specific fields such as `cutoff` (adaptive), `q_lower`/`q_upper`
#'     (pctl), `cr_lower`/`cr_upper` (cr), `delta` (categorical ps),
#'     or `lambda` (optimal)
#'
#' @references
#' Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing
#' with limited overlap in estimation of average treatment effects.
#' *Biometrika*, 96(1), 187--199.
#'
#' Walker, A. M., Patrick, A. R., Lauer, M. S., et al. (2013). A tool for
#' assessing the feasibility of comparative effectiveness research.
#' *Comparative Effectiveness Research*, 3, 11--20.
#'
#' Yang, S., Imbens, G. W., Cui, Z., Faries, D. E., & Kadziola, Z. (2016).
#' Propensity score matching and subclassification in observational studies
#' with multi-level treatments. *Biometrics*, 72(4), 1055--1065.
#'
#' @seealso [ps_trunc()] for bounding (winsorizing) instead of discarding,
#'   [ps_refit()] to re-estimate propensity scores after trimming,
#'   [ps_calibrate()] for calibration-based adjustment,
#'   [ps_trim_meta()] to inspect trimming metadata,
#'   [is_ps_trimmed()] and [is_unit_trimmed()] for logical queries.
#'
#' @examples
#' set.seed(2)
#' n <- 300
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(1.3 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' # Fixed threshold trimming (default)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#' trimmed
#'
#' # How many observations were trimmed?
#' sum(is_unit_trimmed(trimmed))
#'
#' # Data-driven adaptive trimming
#' ps_trim(ps, method = "adaptive")
#'
#' # Quantile-based trimming at 5th and 95th percentiles
#' ps_trim(ps, method = "pctl")
#'
#' # Refit after trimming, then compute weights
#' trimmed <- ps_trim(ps, method = "adaptive")
#' refitted <- ps_refit(trimmed, fit)
#' wt_ate(refitted, .exposure = z)
#'
#' @export
ps_trim <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("ps_trim")
}

#' @export
ps_trim.default <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  method <- rlang::arg_match(method)
  check_ps_range(ps)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "ps_trim"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level

  if (method == "ps") {
    if (is.null(lower)) {
      lower <- 0.1
    }
    if (is.null(upper)) {
      upper <- 0.9
    }
    check_lower_upper(lower, upper)
  } else if (method == "adaptive") {
    if (!is.null(lower) || !is.null(upper)) {
      warn(
        "For {.code method = 'adaptive'}, {.code lower} and {.code upper} are ignored."
      )
    }
  } else if (method == "pctl") {
    if (is.null(lower)) {
      lower <- 0.05
    }
    if (is.null(upper)) upper <- 0.95
  } else if (method == "pref") {
    if (is.null(lower)) {
      lower <- 0.3
    }
    if (is.null(upper)) upper <- 0.7
  } else {
    if (!is.null(lower) || !is.null(upper)) {
      warn(
        "For {.code method = 'cr'}, {.code lower} and {.code upper} are ignored."
      )
    }
  }

  n <- length(ps)
  keep_idx <- integer(0)
  trimmed_idx <- integer(0)
  meta_list <- list(method = method, lower = lower, upper = upper)

  # Decide which indices are kept
  if (method == "ps") {
    keep_idx <- which(ps >= lower & ps <= upper)
  } else if (method == "adaptive") {
    sum_wt <- 1 / (ps * (1 - ps))
    k <- 2 * mean(sum_wt) - max(sum_wt)

    if (k >= 0) {
      cutoff <- 0
    } else {
      trim_fun <- function(x) {
        sum_wt_trim <- sum_wt[sum_wt <= x]
        x - 2 * mean(sum_wt_trim)
      }
      rng <- range(sum_wt)
      lambda <- uniroot(trim_fun, lower = rng[1], upper = rng[2])$root
      cutoff <- 0.5 - sqrt(0.25 - 1 / lambda)
    }
    meta_list$cutoff <- cutoff
    keep_idx <- which(pmin(ps, 1 - ps) > cutoff)
  } else if (method == "pctl") {
    q_lower <- quantile(ps, probs = lower)
    q_upper <- quantile(ps, probs = upper)
    meta_list$q_lower <- q_lower
    meta_list$q_upper <- q_upper
    keep_idx <- which(ps >= q_lower & ps <= q_upper)
  } else if (method == "pref") {
    if (is.null(.exposure)) {
      abort(
        "For {.code method = 'pref'}, must supply {.arg exposure}.",
        error_class = "propensity_missing_arg_error"
      )
    }
    .exposure <- transform_exposure_binary(
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
    prop_exposure <- mean(.exposure)
    pref_score <- plogis(qlogis(ps) - qlogis(prop_exposure))
    meta_list$P <- prop_exposure
    keep_idx <- which(pref_score >= lower & pref_score <= upper)
  } else {
    if (is.null(.exposure)) {
      abort(
        "For {.code method = 'cr'}, must supply {.arg exposure}.",
        error_class = "propensity_missing_arg_error"
      )
    }
    .exposure <- transform_exposure_binary(
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
    ps_treat <- ps[.exposure == 1]
    ps_untrt <- ps[.exposure == 0]
    cr_lower <- min(ps_treat)
    cr_upper <- max(ps_untrt)
    meta_list$cr_lower <- cr_lower
    meta_list$cr_upper <- cr_upper

    keep_idx <- which(ps >= cr_lower & ps <= cr_upper)
  }

  trimmed_idx <- setdiff(seq_len(n), keep_idx)

  # Replace trimmed entries with NA
  ps_na <- ps
  ps_na[trimmed_idx] <- NA_real_

  new_trimmed_ps(
    x = ps_na,
    ps_trim_meta = c(
      meta_list,
      list(
        keep_idx = keep_idx,
        trimmed_idx = trimmed_idx
      )
    )
  )
}

#' @export
ps_trim.matrix <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Only ps and optimal are valid for categorical
  method <- rlang::arg_match(method, values = c("ps", "optimal"))

  # Validate exposure for categorical
  if (is.null(.exposure)) {
    abort(
      "`.exposure` must be provided for categorical propensity score trimming.",
      error_class = "propensity_missing_arg_error"
    )
  }

  # Transform to factor and validate
  .exposure <- transform_exposure_categorical(.exposure)

  # Validate matrix
  ps <- check_ps_matrix(ps, .exposure, call = rlang::caller_env())

  n <- nrow(ps)
  k <- ncol(ps)

  # Initialize metadata
  meta_list <- list(method = method, is_matrix = TRUE)

  if (method == "ps") {
    # Symmetric trimming
    if (is.null(lower)) {
      lower <- 0.1
    }
    delta <- lower # Use lower as delta for consistency

    # Validate delta
    if (delta >= 1 / k) {
      warn(
        "Invalid trimming threshold (delta >= 1/k); returning original data",
        warning_class = "propensity_range_warning"
      )
      keep_idx <- seq_len(n)
    } else {
      # Apply symmetric trimming rule: keep if min(propensity scores) > delta
      keep_idx <- which(apply(ps, 1, function(x) min(x) > delta))

      # Check if all treatment groups are preserved
      if (length(unique(.exposure[keep_idx])) < k) {
        warn(
          "One or more groups removed after trimming; returning original data",
          warning_class = "propensity_no_data_warning"
        )
        keep_idx <- seq_len(n)
      }
    }

    meta_list$delta <- delta
  } else {
    # optimal
    # Multi-category optimal trimming (Yang et al., 2016)
    # Calculate sum of inverse propensity scores
    sum_inv_ps <- rowSums(1 / ps)

    # Define trimming function
    trim_fun <- function(x) {
      sum_trim <- sum_inv_ps[sum_inv_ps <= x]
      if (length(sum_trim) == 0) {
        return(x)
      }
      x - 2 * mean(sum_trim) / mean(sum_inv_ps <= x)
    }

    # Check if trimming is needed
    if (trim_fun(max(sum_inv_ps)) < 0) {
      # No valid solution, use maximum + 1
      lambda <- max(sum_inv_ps) + 1
      keep_idx <- seq_len(n) # Keep all
    } else {
      # Find optimal lambda
      result <- tryCatch(
        {
          uniroot(
            trim_fun,
            lower = min(sum_inv_ps),
            upper = max(sum_inv_ps)
          )$root
        },
        error = function(e) {
          warn(
            "Could not find optimal trimming threshold; using no trimming",
            warning_class = "propensity_convergence_warning"
          )
          NULL
        }
      )

      if (!is.null(result)) {
        lambda <- result
        keep_idx <- which(sum_inv_ps <= lambda)

        # Check if all treatment groups are preserved
        if (length(unique(.exposure[keep_idx])) < k) {
          warn(
            "One or more groups removed after trimming; returning original data",
            warning_class = "propensity_no_data_warning"
          )
          keep_idx <- seq_len(n)
          lambda <- NULL
        }
      } else {
        keep_idx <- seq_len(n)
        lambda <- NULL
      }
    }

    meta_list$lambda <- lambda
  }

  trimmed_idx <- setdiff(seq_len(n), keep_idx)

  # Replace trimmed entries with NA
  ps_na <- ps
  ps_na[trimmed_idx, ] <- NA_real_

  new_trimmed_ps(
    x = ps_na,
    ps_trim_meta = c(
      meta_list,
      list(
        keep_idx = keep_idx,
        trimmed_idx = trimmed_idx
      )
    )
  )
}

#' @export
ps_trim.data.frame <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # For categorical exposures, convert to matrix and call matrix method
  if (!is.null(.exposure)) {
    exposure_type <- detect_exposure_type(.exposure)
    if (exposure_type == "categorical") {
      ps_matrix <- as.matrix(ps)
      return(ps_trim.matrix(
        ps = ps_matrix,
        method = method,
        lower = lower,
        upper = upper,
        .exposure = .exposure,
        ...
      ))
    }
  }

  # For binary exposures, extract appropriate column and call default method
  # This is a simplified version - in practice you might want to handle
  # .propensity_col parameter like in weight functions
  if (ncol(ps) == 2) {
    # Use second column by default for binary
    ps_vec <- ps[[2]]
  } else {
    # Use first column
    ps_vec <- ps[[1]]
  }

  ps_trim.default(
    ps = ps_vec,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .treated = .treated,
    .untreated = .untreated
  )
}

#' @export
ps_trim.ps_trim <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  warn(
    "Propensity scores have already been trimmed. Returning original object.",
    warning_class = "propensity_already_modified_warning"
  )
  ps
}

new_trimmed_ps <- function(x, ps_trim_meta = list()) {
  if (is.matrix(x)) {
    # For matrices, we don't use vctrs
    structure(
      x,
      ps_trim_meta = ps_trim_meta,
      class = c("ps_trim_matrix", "ps_trim", "matrix")
    )
  } else {
    vec_assert(x, ptype = double())
    new_vctr(
      x,
      ps_trim_meta = ps_trim_meta,
      class = "ps_trim",
      inherit_base_type = TRUE
    )
  }
}

#' Test whether propensity scores have been trimmed
#'
#' @description `is_ps_trimmed()` returns `TRUE` if `x` is a `ps_trim` object
#'   or a `psw` object created from trimmed propensity scores, and `FALSE`
#'   otherwise. This tests whether the *object* carries trimming information, not
#'   which individual units were trimmed; see [is_unit_trimmed()] for that.
#'
#' @param x An object to test.
#' @return A logical scalar (`TRUE` or `FALSE`).
#'
#' @seealso [ps_trim()] for trimming propensity scores, [is_unit_trimmed()] to
#'   identify which units were trimmed, [ps_trim_meta()] to retrieve full
#'   trimming metadata.
#'
#' @examples
#' ps <- c(0.05, 0.3, 0.6, 0.95)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' is_ps_trimmed(trimmed)
#' is_ps_trimmed(ps)
#'
#' @export
is_ps_trimmed <- function(x) {
  UseMethod("is_ps_trimmed")
}

#' @export
is_ps_trimmed.default <- function(x) {
  FALSE
}

#' @export
is_ps_trimmed.ps_trim <- function(x) {
  TRUE
}

#' @export
is_ps_trimmed.ps_trim_matrix <- function(x) {
  TRUE
}

#' Identify which units were trimmed
#'
#' @description `is_unit_trimmed()` returns a logical vector indicating which
#'   observations were removed by trimming. This is a per-unit query, as opposed
#'   to [is_ps_trimmed()], which tests whether the object has been trimmed at
#'   all.
#'
#' @param x A `ps_trim` object created by [ps_trim()].
#' @return A logical vector the same length as `x`, where `TRUE` marks a
#'   trimmed unit.
#'
#' @seealso [ps_trim()] for trimming propensity scores, [is_ps_trimmed()] to
#'   test whether an object has been trimmed, [ps_trim_meta()] to retrieve full
#'   trimming metadata.
#'
#' @examples
#' ps <- c(0.05, 0.3, 0.6, 0.95)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' is_unit_trimmed(trimmed)
#'
#' # Use to subset data to retained observations
#' kept <- !is_unit_trimmed(trimmed)
#' ps[kept]
#'
#' @export
is_unit_trimmed <- function(x) {
  UseMethod("is_unit_trimmed")
}

#' @export
is_unit_trimmed.default <- function(x) {
  abort(
    "{.code is_unit_trimmed()} not supported for class {.val {class(x)}}",
    error_class = "propensity_method_error"
  )
}

#' @export
is_unit_trimmed.ps_trim <- function(x) {
  meta <- ps_trim_meta(x)
  out <- vector("logical", length = length(x))
  out[meta$trimmed_idx] <- TRUE

  out
}

#' @export
is_unit_trimmed.ps_trim_matrix <- function(x) {
  meta <- ps_trim_meta(x)
  out <- vector("logical", length = nrow(x))
  out[meta$trimmed_idx] <- TRUE

  out
}


#' @export
`[.ps_trim_matrix` <- function(x, i, j, ..., drop = TRUE) {
  # Get metadata
  meta <- ps_trim_meta(x)

  # Handle single index (matrix as vector) - bypass ps_trim method
  if (nargs() == 2) {
    return(unclass(x)[i])
  }

  # Handle missing i (all rows) - bypass ps_trim method
  if (missing(i)) {
    return(unclass(x)[, j, ..., drop = drop])
  }

  # For single element extraction, call base method directly to avoid ps_trim method
  # Check if this will result in a single element
  if (!missing(j) && length(i) == 1 && length(j) == 1) {
    return(as.numeric(unclass(x)[i, j, drop = TRUE]))
  }

  # Perform subsetting with base matrix method to avoid calling [.ps_trim
  result <- unclass(x)[i, j, ..., drop = drop]

  # If not a matrix anymore or dropped dimensions, return as-is (no metadata)
  if (!is.matrix(result)) {
    return(result)
  }

  # If result is a single element, return as numeric (no metadata)
  if (length(result) == 1) {
    return(as.numeric(result))
  }

  # Handle different index types for rows
  n <- nrow(x)

  if (is.logical(i)) {
    # Convert logical to positions
    i <- which(i)
  } else if (any(i < 0)) {
    # Handle negative indexing
    i <- setdiff(seq_len(n), -i)
  }

  # Map indices to new positions
  new_keep_idx <- integer(0)
  new_trimmed_idx <- integer(0)

  for (idx in seq_along(i)) {
    old_pos <- i[idx]
    if (old_pos %in% meta$keep_idx) {
      new_keep_idx <- c(new_keep_idx, idx)
    }
    if (old_pos %in% meta$trimmed_idx) {
      new_trimmed_idx <- c(new_trimmed_idx, idx)
    }
  }

  # Update metadata
  new_meta <- meta
  new_meta$keep_idx <- new_keep_idx
  new_meta$trimmed_idx <- new_trimmed_idx

  attr(result, "ps_trim_meta") <- new_meta
  class(result) <- c("ps_trim_matrix", "ps_trim", "matrix")
  result
}


# Print methods for ps_trim_matrix

#' @export
print.ps_trim_matrix <- function(x, ..., n = NULL) {
  meta <- ps_trim_meta(x)
  n_rows <- nrow(x)
  k <- ncol(x)
  n_trimmed <- length(meta$trimmed_idx)

  # Create header
  cat(sprintf(
    "<ps_trim_matrix[%d x %d]; trimmed %d of %d; method=%s>\n",
    n_rows,
    k,
    n_trimmed,
    n_rows,
    meta$method
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

# vctrs machinery for ps_trim

#' @export
vec_ptype_abbr.ps_trim <- function(x, ...) "ps_trim"

#' @export
vec_ptype_full.ps_trim <- function(x, ...) {
  paste(
    "ps_trim;",
    "trimmed",
    length(ps_trim_meta(x)$trimmed_idx),
    "of "
  )
}

#' @export
#' @method vec_arith ps_trim
vec_arith.ps_trim <- function(op, x, y, ...) {
  UseMethod("vec_arith.ps_trim", y)
}

#' @export
#' @method vec_arith.ps_trim default
vec_arith.ps_trim.default <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim ps_trim
vec_arith.ps_trim.ps_trim <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim MISSING
vec_arith.ps_trim.MISSING <- function(op, x, y, ...) {
  switch(
    op,
    `-` = -1 * vec_data(x), # Returns numeric
    `+` = vec_data(x), # Returns numeric
    stop_incompatible_op(op, x, y)
  )
}

#' @export
#' @method vec_arith.ps_trim numeric
vec_arith.ps_trim.numeric <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_trim
vec_arith.numeric.ps_trim <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim integer
vec_arith.ps_trim.integer <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim list
vec_arith.ps_trim.list <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
vec_ptype2.ps_trim.ps_trim <- function(x, y, ...) {
  x_meta <- ps_trim_meta(x)
  y_meta <- ps_trim_meta(y)

  # Check if trim parameters match
  if (
    !identical(x_meta$lower, y_meta$lower) ||
      !identical(x_meta$upper, y_meta$upper) ||
      !identical(x_meta$method, y_meta$method)
  ) {
    warn_incompatible_metadata(
      x,
      y,
      "different trimming parameters"
    )
    return(double())
  }

  # Check if refit status matches
  if (!identical(x_meta$refit, y_meta$refit)) {
    warn_incompatible_metadata(
      x,
      y,
      "different refit status"
    )
    return(double())
  }

  # If parameters match, return ps_trim prototype
  # The actual index combining will happen in vec_c
  # Handle missing metadata gracefully
  if (
    is.null(x_meta$method) || is.null(x_meta$lower) || is.null(x_meta$upper)
  ) {
    # Return basic ps_trim if metadata is incomplete
    new_trimmed_ps(double(), ps_trim_meta = x_meta)
  } else {
    ps_trim(
      double(),
      method = x_meta$method,
      lower = x_meta$lower,
      upper = x_meta$upper
    )
  }
}
#' @export
vec_ptype2.ps_trim.double <- function(x, y, ...) {
  warn_class_downgrade("ps_trim")
  double()
}
#' @export
vec_ptype2.double.ps_trim <- function(x, y, ...) {
  warn_class_downgrade("ps_trim")
  double()
}

#' @export
vec_cast.ps_trim.ps_trim <- function(x, to, ...) {
  # Check if metadata matches (excluding indices)
  x_meta <- ps_trim_meta(x)
  to_meta <- ps_trim_meta(to)

  if (
    !identical(x_meta$lower, to_meta$lower) ||
      !identical(x_meta$upper, to_meta$upper) ||
      !identical(x_meta$method, to_meta$method)
  ) {
    vctrs::stop_incompatible_cast(x, to, x_arg = "", to_arg = "")
  }

  # Return x as-is if metadata matches
  x
}

#' @export
vec_cast.double.ps_trim <- function(x, to, ...) {
  # degrade to numeric with NAs
  vec_data(x)
}

#' @export
vec_cast.ps_trim.double <- function(x, to, ...) {
  # create a default ps_trim with no trimming
  new_trimmed_ps(
    x,
    ps_trim_meta = list(
      method = "unknown",
      keep_idx = seq_along(x),
      trimmed_idx = integer(0)
    )
  )
}

#' @export
vec_ptype2.psw.ps_trim <- function(x, y, ...) {
  warn_class_downgrade(c("psw", "ps_trim"))
  double()
}

#' @export
vec_ptype2.ps_trim.psw <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trim", "psw"))
  double()
}

#' @export
vec_cast.character.ps_trim <- function(x, to, ...) as.character(vec_data(x))

#' @export
vec_ptype2.ps_trim.integer <- function(x, y, ...) {
  warn_class_downgrade("ps_trim", "integer")
  integer()
}
#' @export
vec_ptype2.integer.ps_trim <- function(x, y, ...) {
  warn_class_downgrade("ps_trim", "integer")
  integer()
}

#' @export
vec_ptype2.ps_trim.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trim", "ps_trunc"))
  double()
}
#' @export
vec_ptype2.ps_trunc.ps_trim <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trunc", "ps_trim"))
  double()
}

#' @export
vec_cast.integer.ps_trim <- function(x, to, ...) {
  as.integer(vec_data(x))
}

#' @export
vec_cast.ps_trim.integer <- function(x, to, ...) {
  xx <- as.double(x)
  new_trimmed_ps(
    xx,
    ps_trim_meta = list(
      method = "unknown",
      keep_idx = seq_along(xx),
      trimmed_idx = integer(0)
    )
  )
}

#' @export
vec_cast.ps_trim.ps_trunc <- function(x, to, ...) {
  # Convert ps_trunc to ps_trim (no trimming, just convert)
  new_trimmed_ps(
    vec_data(x),
    ps_trim_meta = list(
      method = "from_trunc",
      keep_idx = seq_along(x),
      trimmed_idx = integer(0)
    )
  )
}

#' @export
vec_cast.ps_trunc.ps_trim <- function(x, to, ...) {
  # Convert ps_trim to ps_trunc (ignore NAs)
  ps_trunc(vec_data(x), method = "ps", lower = 0, upper = 1)
}

#' @export
vec_math.ps_trim <- function(.fn, .x, ...) {
  vec_math_base(.fn, vec_data(.x), ...)
}

#' @export
Summary.ps_trim <- function(..., na.rm = FALSE) {
  .fn <- .Generic
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call(.fn, c(numeric_args, list(na.rm = na.rm)))
}

#' @export
min.ps_trim <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("min", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
max.ps_trim <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("max", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
range.ps_trim <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("range", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
median.ps_trim <- function(x, na.rm = FALSE, ...) {
  median(vec_data(x), na.rm = na.rm, ...)
}


#' @export
`[.ps_trim` <- function(x, i, ...) {
  # If i is missing, just call NextMethod
  if (missing(i)) {
    return(NextMethod())
  }

  # Get original metadata
  meta <- ps_trim_meta(x)

  # Convert i to positive integer indices if needed
  # For logical indices, only allow length 1 or length n
  if (is.logical(i) && length(i) != 1 && length(i) != length(x)) {
    abort(
      "Logical subscript `i` must be size 1 or {length(x)}, not {length(i)}.",
      error_class = "propensity_length_error"
    )
  }
  i <- vec_as_location(i, n = length(x))

  # Get the subset of data using NextMethod to handle vctrs subsetting
  result <- NextMethod()

  # Update indices: map old positions to new positions
  new_trimmed_idx <- match(meta$trimmed_idx, i)
  new_trimmed_idx <- new_trimmed_idx[!is.na(new_trimmed_idx)]

  new_keep_idx <- match(meta$keep_idx, i)
  new_keep_idx <- new_keep_idx[!is.na(new_keep_idx)]

  # Update metadata
  new_meta <- meta
  new_meta$trimmed_idx <- new_trimmed_idx
  new_meta$keep_idx <- new_keep_idx

  # Update the attributes on the result
  attr(result, "ps_trim_meta") <- new_meta

  result
}

#' @export
sort.ps_trim <- function(x, decreasing = FALSE, na.last = NA, ...) {
  # Get original metadata
  meta <- ps_trim_meta(x)

  # Get numeric data
  x_data <- vec_data(x)

  # Create a tracking vector to know which NAs are from trimming
  # TRUE = trimmed, FALSE = not trimmed (includes original NAs)
  is_trimmed <- logical(length(x))
  is_trimmed[meta$trimmed_idx] <- TRUE

  # Get the order
  ord <- order(x_data, na.last = na.last, decreasing = decreasing, ...)

  # Apply the ordering to both data and tracking vector
  sorted_data <- x_data[ord]
  sorted_is_trimmed <- is_trimmed[ord]

  # Find new positions of trimmed values
  new_trimmed_idx <- which(sorted_is_trimmed)
  new_keep_idx <- which(!sorted_is_trimmed & !is.na(sorted_data))

  # Create new metadata with updated indices
  new_meta <- meta
  new_meta$trimmed_idx <- new_trimmed_idx
  new_meta$keep_idx <- new_keep_idx

  # Create the sorted ps_trim object
  new_trimmed_ps(sorted_data, ps_trim_meta = new_meta)
}

#' @export
summary.ps_trim <- function(object, ...) {
  summary(vec_data(object), ...)
}

#' @importFrom stats na.omit
#' @export
na.omit.ps_trim <- function(object, ...) {
  # Get metadata
  meta <- ps_trim_meta(object)

  # Get non-NA values and their positions
  not_na <- !is.na(object)
  kept_positions <- which(not_na)

  # Get the clean data
  clean_data <- vec_data(object)[not_na]

  # Update metadata - only keep indices that are in kept_positions
  new_trimmed_idx <- match(meta$trimmed_idx, kept_positions)
  new_trimmed_idx <- new_trimmed_idx[!is.na(new_trimmed_idx)]

  new_keep_idx <- match(meta$keep_idx, kept_positions)
  new_keep_idx <- new_keep_idx[!is.na(new_keep_idx)]

  # Create new metadata
  new_meta <- meta
  new_meta$trimmed_idx <- new_trimmed_idx
  new_meta$keep_idx <- new_keep_idx

  # Create the result with na.action attribute
  result <- new_trimmed_ps(clean_data, ps_trim_meta = new_meta)

  # Add na.action attribute as base na.omit does
  attr(result, "na.action") <- which(!not_na)
  class(attr(result, "na.action")) <- "omit"

  result
}

#' @export
vec_restore.ps_trim <- function(x, to, ...) {
  # Get the prototype's metadata
  to_meta <- ps_trim_meta(to)

  # For combining multiple ps_trim objects, we need to reconstruct indices
  # based on which values are NA (these were trimmed)
  if (length(to_meta$trimmed_idx) == 0 && length(x) > 0) {
    # Identify which positions have NA values (these were trimmed)
    na_positions <- which(is.na(x))

    # Update metadata with the NA positions as trimmed indices
    new_meta <- to_meta
    new_meta$trimmed_idx <- na_positions
    new_meta$keep_idx <- setdiff(seq_along(x), na_positions)

    # Use the constructor to create a proper ps_trim object
    # vec_data in case x is already a vctr
    return(new_trimmed_ps(vec_data(x), ps_trim_meta = new_meta))
  }

  # Use the constructor with the prototype's metadata
  # vec_data in case x is already a vctr
  new_trimmed_ps(vec_data(x), ps_trim_meta = to_meta)
}

#' @export
quantile.ps_trim <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, ...) {
  quantile(vec_data(x), probs = probs, na.rm = na.rm, ...)
}

#' @export
anyDuplicated.ps_trim <- function(x, incomparables = FALSE, ...) {
  anyDuplicated(vec_data(x), incomparables = incomparables, ...)
}

#' @export
diff.ps_trim <- function(x, lag = 1L, differences = 1L, ...) {
  diff(vec_data(x), lag = lag, differences = differences, ...)
}

#' Refit a Propensity Score Model on Retained Observations
#'
#' @description
#' Re-estimates a propensity score model using only the observations retained
#' after trimming. This is the recommended intermediate step between
#' [ps_trim()] and weight calculation (e.g. [wt_ate()]):
#'
#' **`ps_trim()` -> `ps_refit()` -> `wt_*()`**
#'
#' Trimming changes the target population by removing observations with extreme
#' propensity scores. Refitting the model on the retained subset produces
#' propensity scores that better reflect this population, improving both model
#' fit and downstream weight estimation. Weight functions warn if a trimmed
#' propensity score has not been refit.
#'
#' @param trimmed_ps A `ps_trim` object returned by [ps_trim()].
#' @param model The original fitted model used to estimate the propensity
#'   scores (e.g. a [glm][stats::glm] or [multinom][nnet::multinom] object).
#'   The model is refit via [update()][stats::update] on the retained subset.
#' @param .data A data frame. If `NULL` (the default), the data are extracted
#'   from `model` via [model.frame()][stats::model.frame].
#' @param ... Additional arguments passed to [update()][stats::update].
#'
#' @return A `ps_trim` object with re-estimated propensity scores for retained
#'   observations and `NA` for trimmed observations. Use [is_refit()] to
#'   confirm refitting was applied.
#'
#' @seealso [ps_trim()] for the trimming step, [is_refit()] to check refit
#'   status, [wt_ate()] and other weight functions for the next step in the
#'   pipeline.
#'
#' @examples
#' set.seed(2)
#' n <- 200
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#'
#' # fit a propensity score model
#' ps_model <- glm(z ~ x, family = binomial)
#' ps <- predict(ps_model, type = "response")
#'
#' # trim -> refit -> weight pipeline
#' trimmed <- ps_trim(ps, lower = 0.1, upper = 0.9)
#' refit <- ps_refit(trimmed, ps_model)
#' wts <- wt_ate(refit, .exposure = z)
#'
#' is_refit(refit)
#'
#' @export
ps_refit <- function(trimmed_ps, model, .data = NULL, ...) {
  assert_class(trimmed_ps, "ps_trim")
  meta <- ps_trim_meta(trimmed_ps)

  if (length(meta$keep_idx) == 0) {
    abort(
      "No retained rows to refit on (all were trimmed).",
      error_class = "propensity_no_data_error"
    )
  }

  if (is.null(.data)) {
    .data <- model.frame(model)
  }

  # Get the number of observations
  n_obs <- if (is.matrix(trimmed_ps)) nrow(trimmed_ps) else length(trimmed_ps)

  if (nrow(.data) != n_obs) {
    abort(
      c(
        "{.arg .data} must have the same number of rows as observations in {.arg trimmed_ps}.",
        x = "{.arg .data} has {nrow(.data)} row{?s}.",
        x = "{.arg trimmed_ps} has {n_obs} observation{?s}."
      ),
      error_class = "propensity_length_error"
    )
  }

  # refit on untrimmed rows
  data_sub <- .data[meta$keep_idx, , drop = FALSE]
  refit_model <- stats::update(model, data = data_sub, ...)

  # predict new PS for all rows
  if (is.matrix(trimmed_ps)) {
    # For matrix propensity scores (categorical exposures)
    new_ps <- matrix(NA_real_, nrow = n_obs, ncol = ncol(trimmed_ps))
    colnames(new_ps) <- colnames(trimmed_ps)

    # Predict probabilities for retained observations
    if (inherits(refit_model, "multinom")) {
      # For multinomial models from nnet
      pred_probs <- stats::predict(
        refit_model,
        newdata = data_sub,
        type = "probs"
      )
      # Ensure it's a matrix
      if (!is.matrix(pred_probs)) {
        pred_probs <- matrix(pred_probs, nrow = 1)
      }
      new_ps[meta$keep_idx, ] <- pred_probs
    } else {
      # Generic prediction
      new_ps[meta$keep_idx, ] <- stats::predict(
        refit_model,
        newdata = data_sub,
        type = "response"
      )
    }
  } else {
    # For vector propensity scores (binary exposures)
    new_ps <- rep(NA_real_, n_obs)
    new_ps[meta$keep_idx] <- stats::predict(
      refit_model,
      newdata = data_sub,
      type = "response"
    )
  }

  meta$refit <- TRUE

  new_trimmed_ps(
    x = new_ps,
    ps_trim_meta = meta
  )
}

#' Check if propensity scores have been refit
#'
#' `is_refit()` tests whether `x` is a `ps_trim` object whose propensity
#' model has been refit on the retained (non-trimmed) observations via
#' [ps_refit()].
#'
#' @param x An object to test (typically a [ps_trim] vector).
#' @return A single `TRUE` or `FALSE`.
#'
#' @seealso [ps_refit()] to refit a propensity model after trimming,
#'   [ps_trim()] to trim propensity scores.
#'
#' @examples
#' set.seed(2)
#' n <- 30
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' trimmed <- ps_trim(ps, lower = 0.2, upper = 0.8)
#' is_refit(trimmed)
#'
#' refit <- ps_refit(trimmed, fit)
#' is_refit(refit)
#'
#' @export
is_refit <- function(x) {
  UseMethod("is_refit")
}

#' @export
is_refit.default <- function(x) {
  FALSE
}

#' @export
is_refit.ps_trim <- function(x) {
  meta <- ps_trim_meta(x)
  isTRUE(meta$refit)
}

#' Extract trimming metadata from a `ps_trim` object
#'
#' @description `ps_trim_meta()` returns the metadata list attached to a
#'   `ps_trim` object by [ps_trim()].
#'
#' @param x A `ps_trim` object.
#' @return A named list with elements:
#' \describe{
#'   \item{`method`}{Character string indicating the trimming method used.}
#'   \item{`keep_idx`}{Integer vector of retained observation indices.}
#'   \item{`trimmed_idx`}{Integer vector of trimmed observation indices.}
#'   \item{`lower`, `upper`}{Numeric cutoffs, when applicable.}
#'   \item{`refit`}{Logical, `TRUE` if the model was refit via [ps_refit()].}
#' }
#' Additional method-specific elements (e.g. `cutoff`, `delta`, `lambda`) may
#' also be present.
#'
#' @seealso [ps_trim()] for trimming propensity scores, [is_ps_trimmed()] and
#'   [is_unit_trimmed()] for predicate queries.
#'
#' @examples
#' ps <- c(0.05, 0.3, 0.6, 0.95)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' ps_trim_meta(trimmed)
#'
#' @export
ps_trim_meta <- function(x) {
  attr(x, "ps_trim_meta")
}
