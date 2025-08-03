#' Truncate (Winsorize) Propensity Scores
#'
#' **`ps_trunc()`** sets out‐of‐range propensity scores to fixed bounding values
#' (a form of *winsorizing*). This is an alternative to [ps_trim()], which removes
#' (sets `NA`) instead of bounding and is then refit with [ps_refit()]
#'
#' @param ps The propensity score, either a numeric vector between 0 and 1 for
#'   binary exposures, or a matrix/data.frame where each column represents
#'   propensity scores for each level of a categorical exposure.
#' @param .exposure For method "cr", a binary exposure vector. For categorical
#'   exposures, must be a factor or character vector.
#' @param method One of `"ps"`, `"pctl"`, or `"cr"`.
#'   * `"ps"`: directly cut on `[lower, upper]` of `ps`. For categorical, uses
#'     symmetric truncation with `lower` as the threshold.
#'   * `"pctl"`: use quantiles of `ps` as bounding values. For categorical,
#'     calculates quantiles across all propensity score values.
#'   * `"cr"`: the common range of `ps` given `.exposure`, bounding
#'     `[min(ps[treated]), max(ps[untreated])]` (binary only)
#' @param lower,upper Numeric or quantile bounds. If `NULL`, defaults vary by method.
#'   For categorical exposures with method `"ps"`, `lower` represents the
#'   truncation threshold (delta).
#' @inheritParams wt_ate
#' @param ... Additional arguments passed to methods
#'
#' @details
#' For binary exposures with each \eqn{ps[i]}:
#'  - If \eqn{ps[i] < lower\_bound}, we set \eqn{ps[i] = lower\_bound}.
#'  - If \eqn{ps[i] > upper\_bound}, we set \eqn{ps[i] = upper\_bound}.
#'
#' For categorical exposures:
#'  - Each value below the threshold is set to the threshold
#'  - Rows are renormalized to sum to 1
#'
#' This approach is often called *winsorizing*.
#'
#' @return A **`ps_trunc`** object (numeric vector or matrix). It has an attribute
#'   `ps_trunc_meta` storing fields like `method`, `lower_bound`, and
#'   `upper_bound`.
#' @seealso [ps_trim()] and [ps_refit()] for removing extreme values vs. bounding
#'
#' @examples
#' set.seed(2)
#' n <- 30
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' # truncate just the 99th percentile
#' ps_trunc(ps, method = "pctl", lower = 0, upper = .99)
#'
#' @export
ps_trunc <- function(
  ps,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .treated = NULL,
  .untreated = NULL,
  ...
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  method <- rlang::arg_match(method)
  meta_list <- list(method = method)

  check_ps_range(ps)

  if (method == "ps") {
    if (is.null(lower)) lower <- 0.1
    if (is.null(upper)) upper <- 0.9
    check_lower_upper(lower, upper)

    lb <- lower
    ub <- upper
  } else if (method == "pctl") {
    if (is.null(lower)) lower <- 0.05
    if (is.null(upper)) upper <- 0.95
    lb <- quantile(ps, probs = lower)
    ub <- quantile(ps, probs = upper)
    meta_list$lower_pctl <- lower
    meta_list$upper_pctl <- upper
  } else {
    .exposure <- transform_exposure_binary(
      .exposure,
      .treated = .treated,
      .untreated = .untreated
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
  .treated = NULL,
  .untreated = NULL,
  ...
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
  ps <- check_ps_matrix(ps, .exposure)

  n <- nrow(ps)
  k <- ncol(ps)

  if (method == "ps") {
    # Symmetric truncation
    if (is.null(lower)) lower <- 0.01 # Default threshold
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
    if (is.null(lower)) lower <- 0.01 # Default percentile
    if (is.null(upper)) upper <- 0.99 # Default percentile

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
  .treated = NULL,
  .untreated = NULL,
  ...
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
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
ps_trunc.ps_trunc <- function(
  ps,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .treated = NULL,
  .untreated = NULL,
  ...
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

#' @title Extract `ps_trunc` metadata
#' @description Returns the internal metadata list for a `ps_trunc` object.
#' @param x A **`ps_trunc`** object.
#' @return A named list of metadata.
#' @export
ps_trunc_meta <- function(x) {
  attr(x, "ps_trunc_meta")
}

#' Check if object is truncated
#'
#' @description `is_ps_truncated()` is an S3 generic that returns `TRUE` if its
#' argument represents a `ps_trunc` object or `psw` object created from
#' truncated propensity scores. `is_ps_truncated()` is a question about whether
#' the propensity scores *have* been truncated, as opposed to
#' [is_unit_truncated()], which is a question about which *units* have been
#' truncated.
#'
#' @param x An object.
#' @return A logical scalar (`TRUE` or `FALSE`).
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

#' Check if units have been truncated
#'
#' @description `is_ps_truncated()` is an S3 generic that returns a vector of
#'   `TRUE` or `FALSE`, representing if the element has been truncated.
#'   [is_unit_truncated()] is a question about which *units* have been
#'   truncated, as opposed to `is_ps_truncated()`, which is a question about
#'   whether the propensity scores *have* been truncated.
#'
#' @param x An object.
#' @return A logical vector.
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
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc ps_trunc
vec_arith.ps_trunc.ps_trunc <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc numeric
vec_arith.ps_trunc.numeric <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_trunc
vec_arith.numeric.ps_trunc <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc integer
vec_arith.ps_trunc.integer <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

# Combining / Casting
#' @export
vec_ptype2.ps_trunc.ps_trunc <- function(x, y, ...) {
  stop_incompatible_type(x, y, message = "Can't combine ps_trunc objects.")
}

#' @export
vec_ptype2.ps_trunc.double <- function(x, y, ...) double()

#' @export
vec_ptype2.double.ps_trunc <- function(x, y, ...) double()

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
vec_ptype2.psw.ps_trunc <- function(x, y, ...) character()

#' @export
vec_ptype2.ps_trunc.psw <- function(x, y, ...) character()

#' @export
vec_cast.character.ps_trunc <- function(x, to, ...) as.character(vec_data(x))

#' @export
vec_ptype2.ps_trunc.integer <- function(x, y, ...) integer()

#' @export
vec_ptype2.integer.ps_trunc <- function(x, y, ...) integer()

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
