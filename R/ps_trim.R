#' Trim Propensity Scores
#'
#' `ps_trim()` applies trimming methods to a propensity-score vector or matrix,
#' returning a new vector/matrix of the *same length/dimensions*, with trimmed
#' entries replaced by `NA.` You can inspect further metadata in `ps_trim_meta(x)`.
#' After running `ps_trim()`, you should refit the model with `ps_refit()`.
#'
#' @param ps The propensity score, either a numeric vector between 0 and 1 for
#'   binary exposures, or a matrix/data.frame where each column represents
#'   propensity scores for each level of a categorical exposure.
#' @param .exposure For methods like `"pref"` or `"cr"`, a vector for a binary
#'   exposure. For categorical exposures with method `"optimal"`, must be a
#'   factor or character vector.
#' @param method One of `c("ps", "adaptive", "pctl", "pref", "cr", "optimal")`.
#'   For categorical exposures, only `"ps"` and `"optimal"` are supported.
#' @param lower,upper Numeric cutoffs or quantiles. If `NULL`, defaults vary by
#'   method. For categorical exposures with method `"ps"`, `lower` represents the
#'   symmetric trimming threshold (delta).
#' @inheritParams wt_ate
#' @param ... Additional arguments passed to methods
#'
#' @details The returned object is a **`ps_trim`** vector/matrix of the same
#' length/dimensions as `ps`, but with trimmed entries replaced by `NA`. An
#' attribute `ps_trim_meta` contains:
#'
#' - `method`: Which trimming method was used
#' - `keep_idx`: Indices retained
#' - `trimmed_idx`: Indices replaced by `NA`
#' - Possibly other fields such as final cutoffs, etc.
#'
#' For categorical exposures:
#' - **Symmetric trimming** (`method = "ps"`): Removes observations where any
#'   propensity score falls below the threshold delta (specified via `lower`).
#' - **Optimal trimming** (`method = "optimal"`): Uses the Yang et al. (2016)
#'   approach for multi-category treatments.
#'
#' **Arithmetic behavior**: Arithmetic operations on `ps_trim` objects return
#' numeric vectors, not `ps_trim` objects. This is intentional - once you
#' transform propensity scores (e.g., `1/ps` for weights), the result is no
#' longer a propensity score.
#'
#' **NA handling**: Trimmed values are set to `NA`. Operations that don't handle
#' `NA` values will propagate them (e.g., `sum()` returns `NA` unless
#' `na.rm = TRUE`).
#'
#' **Metadata tracking**: The `trimmed_idx` and `keep_idx` are updated when
#' subsetting or reordering:
#' - Subsetting with `[` updates indices to new positions
#' - `sort()` reorders data and updates indices accordingly
#' - `unique()` may change lengths but preserves the class
#' - `na.omit()` removes trimmed values and updates indices
#'
#' **Combining behavior**: When combining `ps_trim` objects with `c()`, metadata
#' must match (same trimming parameters). Mismatched metadata triggers a warning
#' and returns a numeric vector.
#'
#' @return A `ps_trim` object (numeric vector or matrix). The attribute
#' `ps_trim_meta` stores metadata.
#'
#' @seealso [ps_trunc()] for bounding/winsorizing instead of discarding,
#'   [is_refit()], [is_ps_trimmed()]
#' @examples
#'
#' set.seed(2)
#' n <- 300
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(1.3 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' ps_trim(ps, method = "adaptive")
#'
#' # Coercion behavior with ps_trim objects
#' ps_trim1 <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#' ps_trim2 <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' # Compatible objects combine silently
#' c(ps_trim1[1:50], ps_trim2[51:100])  # Returns ps_trim object
#'
#' # Different trim parameters trigger warning
#' ps_trim3 <- ps_trim(ps, method = "ps", lower = 0.2, upper = 0.8)
#' c(ps_trim1[1:50], ps_trim3[51:100])  # Warning: returns numeric
#'
#' # Cross-class combinations warn and return numeric
#' psw_obj <- psw(ps[1:50], estimand = "ate")
#' c(ps_trim1[1:50], psw_obj)  # Warning: returns numeric
#'
#' @export
ps_trim <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .treated = NULL,
  .untreated = NULL,
  ...
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
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  method <- rlang::arg_match(method)
  check_ps_range(ps)

  if (method == "ps") {
    if (is.null(lower)) lower <- 0.1
    if (is.null(upper)) upper <- 0.9
    check_lower_upper(lower, upper)
  } else if (method == "adaptive") {
    if (!is.null(lower) || !is.null(upper)) {
      warn(
        "For {.code method = 'adaptive'}, {.code lower} and {.code upper} are ignored."
      )
    }
  } else if (method == "pctl") {
    if (is.null(lower)) lower <- 0.05
    if (is.null(upper)) upper <- 0.95
  } else if (method == "pref") {
    if (is.null(lower)) lower <- 0.3
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
      .treated = .treated,
      .untreated = .untreated
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
      .treated = .treated,
      .untreated = .untreated
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
  .treated = NULL,
  .untreated = NULL,
  ...
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
    if (is.null(lower)) lower <- 0.1
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
      if (length(sum_trim) == 0) return(x)
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
  .treated = NULL,
  .untreated = NULL,
  ...
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
    .treated = .treated,
    .untreated = .untreated,
    ...
  )
}

#' @export
ps_trim.ps_trim <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .treated = NULL,
  .untreated = NULL,
  ...
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

#' Check if object is trimmed
#'
#' @description `is_ps_trimmed()` is an S3 generic that returns `TRUE` if its
#'   argument represents a `ps_trim` object or `psw` object created from trimmed
#'   propensity scores. `is_ps_trimmed()` is a question about whether or not the
#'   propensity scores *have* been trimmed, as opposed to [is_unit_trimmed()],
#'   which is a question about which *units* have been trimmed.
#'
#' @param x An object.
#' @return A logical scalar (`TRUE` or `FALSE`).
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

#' Check if units have been trimmed
#'
#' @description `is_unit_trimmed()` is an that vector of `TRUE` or `FALSE`
#'   values, representing if the unit was trimmed. `is_unit_trimmed()` is a
#'   question about which *units* have been trimmed, as opposed to
#'   [is_ps_trimmed()], which is a question about whether or not the propensity
#'   scores *have* been trimmed.
#'
#' @param x An object.
#' @return A logical scalar (`TRUE` or `FALSE`).
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

#' Refit the Propensity Score Model on Retained Observations
#'
#' Takes a `ps_trim` object and the original model
#' used to calculate the propensity score, then:
#' 1. Retrieves data from the model (or from `.df` argument if provided)
#' 2. Subsets rows to the non‐trimmed indices
#' 3. Refits the model
#' 4. Predicts new propensity scores for all rows (trimmed rows -> `NA`)
#' 5. Returns a new `ps_trim` object with `refit = TRUE`.
#'
#' @param trimmed_ps A `ps_trim` object (same length as data, NAs for trimmed).
#' @param model The fitted model used to get the original PS (e.g. a glm).
#' @param .df Optional. A data frame. If `NULL`, we try to retrieve from `model`.
#' @param ... Additional arguments passed to `update()`.
#'
#' @return
#' A new `ps_trim` object with updated propensity scores and
#' `ps_trim_meta(x)$refit` set to `TRUE`.
#'
#' @seealso [ps_trim()], [is_refit()], [is_ps_trimmed()]
#'
#' @examples
#' set.seed(2)
#' n <- 30
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' # trim and refit
#' refit <- ps_trim(ps, lower = .2, upper = .8) |>
#'   ps_refit(fit)
#'
#' is_refit(refit)
#'
#' @export
ps_refit <- function(trimmed_ps, model, .df = NULL, ...) {
  assert_class(trimmed_ps, "ps_trim")
  meta <- ps_trim_meta(trimmed_ps)

  if (length(meta$keep_idx) == 0) {
    abort(
      "No retained rows to refit on (all were trimmed).",
      error_class = "propensity_no_data_error"
    )
  }

  if (is.null(.df)) {
    .df <- model.frame(model)
  }

  # Get the number of observations
  n_obs <- if (is.matrix(trimmed_ps)) nrow(trimmed_ps) else length(trimmed_ps)

  if (nrow(.df) != n_obs) {
    abort(
      c(
        "{.arg .df} must have the same number of rows as observations in {.arg trimmed_ps}.",
        x = "{.arg .df} has {nrow(.df)} row{?s}.",
        x = "{.arg trimmed_ps} has {n_obs} observation{?s}."
      ),
      error_class = "propensity_length_error"
    )
  }

  # refit on untrimmed rows
  data_sub <- .df[meta$keep_idx, , drop = FALSE]
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

#' Check if an object has been refit
#'
#' @description
#' **`is_refit()`** is an S3 generic that returns `TRUE` if its argument
#' represents a [ps_trim] object (or a weighting object) that has had the
#' propensity model refit on the retained subset.
#'
#' @param x An R object (e.g. a [ps_trim] or [psw]).
#' @return A logical scalar (`TRUE` or `FALSE`).
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

#' @title Extract `ps_trim` metadata
#' @description Returns the internal metadata list for a `ps_trim` object.
#' @param x A **`ps_trim`** object.
#' @return A named list of metadata.
#' @export
ps_trim_meta <- function(x) {
  attr(x, "ps_trim_meta")
}
