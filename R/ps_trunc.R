#' Truncate (Winsorize) Propensity Scores
#'
#' `ps_trunc()` bounds extreme propensity scores to fixed limits, replacing
#' out-of-range values with the boundary value (a form of *winsorizing*). The
#' result is a vector or matrix of the same length and dimensions as `ps`, with
#' no observations removed. This contrasts with [ps_trim()], which sets extreme
#' values to `NA` (effectively removing those observations from analysis).
#'
#' @param ps A numeric vector of propensity scores between 0 and 1 (binary
#'   exposures), or a matrix/data.frame where each column contains propensity
#'   scores for one level of a categorical exposure.
#' @param .exposure An exposure vector. Required for method `"cr"` (binary
#'   exposure vector) and for categorical exposures (factor or character vector)
#'   with any method.
#' @param method One of `"ps"`, `"pctl"`, or `"cr"`:
#'   * `"ps"` (default): Truncate directly on propensity score values. Values
#'     outside `[lower, upper]` are set to the nearest bound. For categorical
#'     exposures, applies symmetric truncation using `lower` as the threshold
#'     (delta) and renormalizes rows to sum to 1.
#'   * `"pctl"`: Truncate at quantiles of the propensity score distribution.
#'     The `lower` and `upper` arguments specify quantile probabilities. For
#'     categorical exposures, quantiles are computed across all columns.
#'   * `"cr"`: Truncate to the common range of propensity scores across
#'     exposure groups (binary exposures only). Bounds are
#'     `[min(ps[focal]), max(ps[reference])]`. Requires `.exposure`.
#' @param lower,upper Bounds for truncation. Interpretation depends on `method`:
#'   * `method = "ps"`: Propensity score values (defaults: 0.1 and 0.9). For
#'     categorical exposures, `lower` is the truncation threshold delta
#'     (default: 0.01); `upper` is ignored.
#'   * `method = "pctl"`: Quantile probabilities (defaults: 0.05 and 0.95;
#'     categorical defaults: 0.01 and 0.99).
#'   * `method = "cr"`: Ignored; bounds are determined by the data.
#' @inheritParams wt_ate
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' Unlike [ps_trim()], truncation preserves all observations. No `NA` values
#' are introduced; out-of-range scores are replaced with the boundary value.
#'
#' For **binary exposures**, each propensity score \eqn{e_i} is bounded:
#'
#'  - If \eqn{e_i < l}, set \eqn{e_i = l} (the lower bound).
#'
#'  - If \eqn{e_i > u}, set \eqn{e_i = u} (the upper bound).
#'
#' For **categorical exposures**, values below the threshold are set to the
#' threshold and each row is renormalized to sum to 1.
#'
#' **Arithmetic behavior**: Arithmetic operations on `ps_trunc` objects return
#' plain numeric vectors. Once propensity scores are transformed (e.g., into
#' weights), the result is no longer a propensity score.
#'
#' **Combining behavior**: Combining `ps_trunc` objects with `c()` requires
#' matching truncation parameters. Mismatched parameters produce a warning and
#' return a plain numeric vector.
#'
#' @return A `ps_trunc` object (a numeric vector for binary exposures, or a
#'   matrix for categorical exposures). Use [ps_trunc_meta()] to inspect
#'   metadata including `method`, `lower_bound`, `upper_bound`, and
#'   `truncated_idx` (positions of modified values).
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
#' @seealso [ps_trim()] for removing (rather than bounding) extreme values,
#'   [ps_refit()] for refitting the propensity model after trimming,
#'   [is_ps_truncated()], [is_unit_truncated()], [ps_trunc_meta()]
#'
#' @examples
#' set.seed(2)
#' n <- 200
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(1.2 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' # Truncate to [0.1, 0.9]
#' ps_t <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)
#' ps_t
#'
#' # Truncate at the 1st and 99th percentiles
#' ps_trunc(ps, method = "pctl", lower = 0.01, upper = 0.99)
#'
#' # Use truncated scores to calculate weights
#' wt_ate(ps_t, .exposure = z)
#'
#' # Inspect which observations were truncated
#' is_unit_truncated(ps_t)
#'
#' @export
ps_trunc <- function(
  ps,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("ps_trunc")
}

#' @export
ps_trunc.default <- function(
  ps,
  method = c("ps", "pctl", "cr"),
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
  meta_list <- list(method = method)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "ps_trunc"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level

  check_ps_range(ps)

  if (method == "ps") {
    if (is.null(lower)) {
      lower <- 0.1
    }
    if (is.null(upper)) {
      upper <- 0.9
    }
    check_lower_upper(lower, upper)

    lb <- lower
    ub <- upper
  } else if (method == "pctl") {
    if (is.null(lower)) {
      lower <- 0.05
    }
    if (is.null(upper)) {
      upper <- 0.95
    }
    lb <- quantile(ps, probs = lower)
    ub <- quantile(ps, probs = upper)
    meta_list$lower_pctl <- lower
    meta_list$upper_pctl <- upper
  } else {
    .exposure <- transform_exposure_binary(
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
    ps_treat <- ps[.exposure == 1]
    ps_untrt <- ps[.exposure == 0]
    lb <- min(ps_treat)
    ub <- max(ps_untrt)
  }

  # winsorize
  pinned_low <- which(ps < lb)
  pinned_high <- which(ps > ub)
  truncated_idx <- sort(c(pinned_low, pinned_high))

  ps[pinned_low] <- lb
  ps[pinned_high] <- ub

  meta <- c(
    meta_list,
    list(
      lower_bound = lb,
      upper_bound = ub,
      truncated_idx = truncated_idx
    )
  )

  new_ps_trunc(ps, meta)
}

#' @export
ps_trunc.matrix <- function(
  ps,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # Only ps and pctl are valid for categorical
  method <- rlang::arg_match(method, values = c("ps", "pctl"))

  # Validate exposure for categorical
  if (is.null(.exposure)) {
    abort(
      "`.exposure` must be provided for categorical propensity score truncation.",
      error_class = "propensity_missing_arg_error"
    )
  }

  # Transform to factor and validate
  .exposure <- transform_exposure_categorical(.exposure)

  # Validate matrix
  ps <- check_ps_matrix(ps, .exposure, call = rlang::caller_env())

  n <- nrow(ps)
  k <- ncol(ps)

  if (method == "ps") {
    # Symmetric truncation
    if (is.null(lower)) {
      lower <- 0.01
    } # Default threshold
    delta <- lower # Use lower as delta for consistency

    # Validate delta
    if (delta >= 1 / k) {
      abort(
        "Invalid truncation threshold (delta >= 1/k).",
        error_class = "propensity_range_error"
      )
    }

    # Track which values were truncated
    truncated_idx <- which(apply(ps, 1, function(x) any(x < delta)))

    # Apply truncation and renormalize
    ps_trunc <- ps
    for (i in 1:n) {
      row_vals <- ps_trunc[i, ]
      # Clamp values below delta
      row_vals[row_vals < delta] <- delta
      # Renormalize to sum to 1
      ps_trunc[i, ] <- row_vals / sum(row_vals)
    }

    lower_bound <- delta
    upper_bound <- NA_real_
  } else {
    # pctl
    # Percentile-based truncation
    if (is.null(lower)) {
      lower <- 0.01
    } # Default percentile
    if (is.null(upper)) {
      upper <- 0.99
    } # Default percentile

    # Calculate thresholds based on the distribution of all propensity scores
    all_ps_vals <- as.vector(ps)
    lower_threshold <- quantile(all_ps_vals, probs = lower)
    upper_threshold <- quantile(all_ps_vals, probs = upper)

    # Track which rows had values truncated
    truncated_idx <- which(apply(
      ps,
      1,
      function(x) any(x < lower_threshold | x > upper_threshold)
    ))

    # Apply truncation and renormalize
    ps_trunc <- ps
    for (i in 1:n) {
      row_vals <- ps_trunc[i, ]
      # Clamp values outside thresholds
      row_vals[row_vals < lower_threshold] <- lower_threshold
      row_vals[row_vals > upper_threshold] <- upper_threshold
      # Renormalize to sum to 1
      ps_trunc[i, ] <- row_vals / sum(row_vals)
    }

    lower_bound <- lower_threshold
    upper_bound <- upper_threshold
  }

  meta <- list(
    method = method,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    truncated_idx = truncated_idx,
    is_matrix = TRUE
  )

  # Add percentile info if using pctl method
  if (method == "pctl") {
    meta$lower_pctl <- lower
    meta$upper_pctl <- upper
  }

  new_ps_trunc(ps_trunc, meta)
}

#' @export
ps_trunc.data.frame <- function(
  ps,
  method = c("ps", "pctl", "cr"),
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
      return(ps_trunc.matrix(
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
  if (ncol(ps) == 2) {
    # Use second column by default for binary
    ps_vec <- ps[[2]]
  } else {
    # Use first column
    ps_vec <- ps[[1]]
  }

  ps_trunc.default(
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
ps_trunc.ps_trunc <- function(
  ps,
  method = c("ps", "pctl", "cr"),
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
    "Propensity scores have already been truncated. Returning original object.",
    warning_class = "propensity_already_modified_warning"
  )
  ps
}

new_ps_trunc <- function(x, meta) {
  if (is.matrix(x)) {
    # For matrices, we don't use vctrs
    structure(
      x,
      ps_trunc_meta = meta,
      class = c("ps_trunc_matrix", "ps_trunc", "matrix")
    )
  } else {
    vec_assert(x, ptype = double())
    new_vctr(
      x,
      ps_trunc_meta = meta,
      class = "ps_trunc",
      inherit_base_type = TRUE
    )
  }
}

#' Extract truncation metadata from a `ps_trunc` object
#'
#' @description Returns the metadata list attached to a [`ps_trunc`][ps_trunc()]
#'   object. The list includes fields such as `method`, `lower_bound`,
#'   `upper_bound`, and `truncated_idx`.
#'
#' @param x A `ps_trunc` object created by [ps_trunc()].
#' @return A named list with truncation metadata, including:
#'   * `method` -- the truncation method used (`"ps"`, `"pctl"`, or `"cr"`)
#'   * `lower_bound`, `upper_bound` -- the applied bounds
#'   * `truncated_idx` -- integer positions of values that were winsorized
#'
#' @seealso [ps_trunc()], [is_ps_truncated()], [is_unit_truncated()]
#'
#' @examples
#' ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
#' ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
#' ps_trunc_meta(ps_t)
#'
#' @export
ps_trunc_meta <- function(x) {
  attr(x, "ps_trunc_meta")
}

#' Test whether propensity scores have been truncated
#'
#' @description `is_ps_truncated()` returns `TRUE` if `x` is a `ps_trunc`
#'   object or a `psw` object derived from truncated propensity scores.
#'   Use [is_unit_truncated()] to find out *which* observations were modified.
#'
#' @param x An object.
#' @return A single `TRUE` or `FALSE`.
#'
#' @seealso [ps_trunc()], [is_unit_truncated()], [ps_trunc_meta()]
#'
#' @examples
#' ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
#' is_ps_truncated(ps)
#'
#' ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
#' is_ps_truncated(ps_t)
#'
#' @export
is_ps_truncated <- function(x) {
  UseMethod("is_ps_truncated")
}

#' @export
is_ps_truncated.default <- function(x) {
  FALSE
}

#' @export
is_ps_truncated.ps_trunc <- function(x) {
  TRUE
}

#' @export
is_ps_truncated.ps_trunc_matrix <- function(x) {
  TRUE
}

#' Identify which units were truncated
#'
#' @description `is_unit_truncated()` returns a logical vector indicating which
#'   observations had their propensity scores modified by truncation. Use
#'   [is_ps_truncated()] to test whether an object has been truncated at all.
#'
#' @param x A `ps_trunc` object created by [ps_trunc()].
#' @return A logical vector the same length as `x` (or number of rows for
#'   matrix input). `TRUE` marks observations whose values were winsorized.
#'
#' @seealso [ps_trunc()], [is_ps_truncated()], [ps_trunc_meta()]
#'
#' @examples
#' ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
#' ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
#' is_unit_truncated(ps_t)
#'
#' @export
is_unit_truncated <- function(x) {
  UseMethod("is_unit_truncated")
}

#' @export
is_unit_truncated.default <- function(x) {
  abort(
    "{.code is_unit_truncated()} not supported for class {.val {class(x)}}",
    error_class = "propensity_method_error"
  )
}

#' @export
is_unit_truncated.ps_trunc <- function(x) {
  meta <- ps_trunc_meta(x)
  out <- vector("logical", length = length(x))
  out[meta$truncated_idx] <- TRUE

  out
}

#' @export
is_unit_truncated.ps_trunc_matrix <- function(x) {
  meta <- ps_trunc_meta(x)
  out <- vector("logical", length = nrow(x))
  out[meta$truncated_idx] <- TRUE

  out
}


#' @export
`[.ps_trunc_matrix` <- function(x, i, j, ..., drop = TRUE) {
  # Get metadata
  meta <- ps_trunc_meta(x)

  # Handle single index (matrix as vector) - bypass ps_trunc method
  if (nargs() == 2) {
    return(unclass(x)[i])
  }

  # Handle missing i (all rows) - bypass ps_trunc method
  if (missing(i)) {
    return(unclass(x)[, j, ..., drop = drop])
  }

  # For single element extraction, call base method directly to avoid ps_trunc method
  # Check if this will result in a single element
  if (!missing(j) && length(i) == 1 && length(j) == 1) {
    return(as.numeric(unclass(x)[i, j, drop = TRUE]))
  }

  # Perform subsetting with base matrix method to avoid calling [.ps_trunc
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
  new_truncated_idx <- integer(0)

  for (idx in seq_along(i)) {
    old_pos <- i[idx]
    if (old_pos %in% meta$truncated_idx) {
      new_truncated_idx <- c(new_truncated_idx, idx)
    }
  }

  # Update metadata
  new_meta <- meta
  new_meta$truncated_idx <- new_truncated_idx

  attr(result, "ps_trunc_meta") <- new_meta
  class(result) <- c("ps_trunc_matrix", "ps_trunc", "matrix")
  result
}


# Print methods for ps_trunc_matrix

#' @export
print.ps_trunc_matrix <- function(x, ..., n = NULL) {
  meta <- ps_trunc_meta(x)
  n_rows <- nrow(x)
  k <- ncol(x)
  n_truncated <- length(meta$truncated_idx)

  # Create header
  if (meta$method == "pctl") {
    cat(sprintf(
      "<ps_trunc_matrix[%d x %d]; truncated %d of %d; method=%s[%.2f,%.2f]>\n",
      n_rows,
      k,
      n_truncated,
      n_rows,
      meta$method,
      meta$lower_pctl,
      meta$upper_pctl
    ))
  } else {
    cat(sprintf(
      "<ps_trunc_matrix[%d x %d]; truncated %d of %d; method=%s[%.4f,%.4f]>\n",
      n_rows,
      k,
      n_truncated,
      n_rows,
      meta$method,
      meta$lower_bound,
      ifelse(is.na(meta$upper_bound), Inf, meta$upper_bound)
    ))
  }

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

# vctrs machinery for ps_trunc

#' @export
vec_ptype_abbr.ps_trunc <- function(x, ...) {
  "ps_trunc"
}

#' @export
vec_ptype_full.ps_trunc <- function(x, ...) {
  m <- ps_trunc_meta(x)
  paste0(
    "ps_trunc{[",
    m$lower_bound,
    ",",
    m$upper_bound,
    "], method=",
    m$method,
    "}"
  )
}

# Arithmetic disallowed
#' @export
#' @method vec_arith ps_trunc
vec_arith.ps_trunc <- function(op, x, y, ...) {
  UseMethod("vec_arith.ps_trunc", y)
}

#' @export
#' @method vec_arith.ps_trunc default
vec_arith.ps_trunc.default <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc ps_trunc
vec_arith.ps_trunc.ps_trunc <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc MISSING
vec_arith.ps_trunc.MISSING <- function(op, x, y, ...) {
  switch(
    op,
    `-` = -1 * vec_data(x), # Returns numeric
    `+` = vec_data(x), # Returns numeric
    stop_incompatible_op(op, x, y)
  )
}

#' @export
#' @method vec_arith.ps_trunc numeric
vec_arith.ps_trunc.numeric <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_trunc
vec_arith.numeric.ps_trunc <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc integer
vec_arith.ps_trunc.integer <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

# Combining / Casting
#' @export
vec_ptype2.ps_trunc.ps_trunc <- function(x, y, ...) {
  x_meta <- ps_trunc_meta(x)
  y_meta <- ps_trunc_meta(y)

  # Check if truncation parameters match
  if (
    !identical(x_meta$lower_bound, y_meta$lower_bound) ||
      !identical(x_meta$upper_bound, y_meta$upper_bound) ||
      !identical(x_meta$method, y_meta$method)
  ) {
    warn_incompatible_metadata(
      x,
      y,
      "different truncation parameters"
    )
    return(double())
  }

  # If parameters match, return ps_trunc prototype
  ps_trunc(
    double(),
    method = x_meta$method,
    lower = x_meta$lower_bound,
    upper = x_meta$upper_bound
  )
}

#' @export
vec_ptype2.ps_trunc.double <- function(x, y, ...) {
  warn_class_downgrade("ps_trunc")
  double()
}

#' @export
vec_ptype2.double.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade("ps_trunc")
  double()
}

#' @export
vec_cast.ps_trunc.ps_trunc <- function(x, to, ...) {
  # Check if metadata matches (excluding indices)
  x_meta <- ps_trunc_meta(x)
  to_meta <- ps_trunc_meta(to)

  if (
    !identical(x_meta$lower_bound, to_meta$lower_bound) ||
      !identical(x_meta$upper_bound, to_meta$upper_bound) ||
      !identical(x_meta$method, to_meta$method)
  ) {
    vctrs::stop_incompatible_cast(x, to, x_arg = "", to_arg = "")
  }

  # Return x as-is if metadata matches
  x
}

#' @export
vec_cast.double.ps_trunc <- function(x, to, ...) {
  vec_data(x)
}

#' @export
vec_cast.ps_trunc.double <- function(x, to, ...) {
  new_ps_trunc(
    x,
    meta = list(method = "unknown", lower_bound = NA, upper_bound = NA)
  )
}

#' @export
vec_ptype2.psw.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade(c("psw", "ps_trunc"))
  double()
}

#' @export
vec_ptype2.ps_trunc.psw <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trunc", "psw"))
  double()
}

#' @export
vec_cast.character.ps_trunc <- function(x, to, ...) as.character(vec_data(x))

#' @export
vec_ptype2.ps_trunc.integer <- function(x, y, ...) {
  warn_class_downgrade("ps_trunc", "integer")
  integer()
}

#' @export
vec_ptype2.integer.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade("ps_trunc", "integer")
  integer()
}

#' @export
vec_cast.integer.ps_trunc <- function(x, to, ...) as.integer(vec_data(x))

#' @export
vec_cast.ps_trunc.integer <- function(x, to, ...) {
  xx <- as.double(x)
  new_ps_trunc(
    xx,
    meta = list(method = "unknown", lower_bound = NA, upper_bound = NA)
  )
}

#' @export
vec_math.ps_trunc <- function(.fn, .x, ...) {
  vec_math_base(.fn, vec_data(.x), ...)
}

#' @export
Summary.ps_trunc <- function(..., na.rm = FALSE) {
  .fn <- .Generic
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call(.fn, c(numeric_args, list(na.rm = na.rm)))
}

#' @export
min.ps_trunc <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("min", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
max.ps_trunc <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("max", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
range.ps_trunc <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("range", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
median.ps_trunc <- function(x, na.rm = FALSE, ...) {
  median(vec_data(x), na.rm = na.rm, ...)
}

#' @export
quantile.ps_trunc <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, ...) {
  quantile(vec_data(x), probs = probs, na.rm = na.rm, ...)
}


# Note: We cannot implement custom vec_c for ps_trunc
# vctrs handles combination internally through vec_ptype2 and vec_cast
# When combining ps_trunc objects with same parameters, indices won't be preserved

#' @export
`[.ps_trunc` <- function(x, i, ...) {
  # If i is missing, just call NextMethod
  if (missing(i)) {
    return(NextMethod())
  }

  # Get original metadata
  meta <- ps_trunc_meta(x)

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
  new_truncated_idx <- match(meta$truncated_idx, i)
  new_truncated_idx <- new_truncated_idx[!is.na(new_truncated_idx)]

  # Update metadata
  new_meta <- meta
  new_meta$truncated_idx <- new_truncated_idx

  # Update the attributes on the result
  attr(result, "ps_trunc_meta") <- new_meta

  result
}

#' @export
summary.ps_trunc <- function(object, ...) {
  summary(vec_data(object), ...)
}

#' @export
sort.ps_trunc <- function(x, decreasing = FALSE, na.last = NA, ...) {
  # Get original metadata
  meta <- ps_trunc_meta(x)

  # Get numeric data
  x_data <- vec_data(x)

  # Create a tracking vector to know which values were truncated
  # TRUE = truncated, FALSE = not truncated
  was_truncated <- logical(length(x))
  was_truncated[meta$truncated_idx] <- TRUE

  # Get the order
  ord <- order(x_data, na.last = na.last, decreasing = decreasing, ...)

  # Apply the ordering to both data and tracking vector
  sorted_data <- x_data[ord]
  sorted_was_truncated <- was_truncated[ord]

  # Find new positions of truncated values
  new_truncated_idx <- which(sorted_was_truncated)

  # Create new metadata with updated indices
  new_meta <- meta
  new_meta$truncated_idx <- new_truncated_idx

  # Create the sorted ps_trunc object
  new_ps_trunc(sorted_data, new_meta)
}

#' @export
vec_restore.ps_trunc <- function(x, to, ...) {
  # Get the prototype's metadata
  to_meta <- ps_trunc_meta(to)

  # Extract numeric data for comparisons
  x_data <- vec_data(x)

  # For combining multiple ps_trunc objects, we need to reconstruct indices
  # For ps_trunc, values at boundaries are modified, not NA
  if (
    length(to_meta$truncated_idx) == 0 &&
      length(x_data) > 0 &&
      !is.null(to_meta$lower_bound) &&
      !is.null(to_meta$upper_bound)
  ) {
    # Identify which positions were truncated (at the bounds)
    truncated_positions <- which(
      x_data == to_meta$lower_bound | x_data == to_meta$upper_bound
    )

    # Update metadata with the truncated positions
    new_meta <- to_meta
    new_meta$truncated_idx <- truncated_positions

    # Use the constructor to create a proper ps_trunc object
    return(new_ps_trunc(x_data, new_meta))
  }

  # Use the constructor with the prototype's metadata
  new_ps_trunc(x_data, to_meta)
}

#' @export
anyDuplicated.ps_trunc <- function(x, incomparables = FALSE, ...) {
  anyDuplicated(vec_data(x), incomparables = incomparables, ...)
}

#' @export
diff.ps_trunc <- function(x, lag = 1L, differences = 1L, ...) {
  diff(vec_data(x), lag = lag, differences = differences, ...)
}
