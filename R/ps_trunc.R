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
#' ## The truncation record
#'
#' A `ps_trunc` records which units were winsorized as positions among the
#' observations it was written for, along with how many observations that was.
#' Operations that hand this package the subscript re-index those positions onto
#' the result: subsetting with `[`, [sort()], [unique()], and [rep()] all return
#' a record written for what they return, and a subscript naming a position more
#' than once reports that unit at every place it now holds.
#'
#' Operations that change how many observations there are without supplying a
#' subscript cannot re-index the record, and it is dropped rather than worked
#' out from the values, since a score that arrived equal to a bound is
#' indistinguishable from one this function pinned there. [vctrs::vec_slice()],
#' which is how filtering, joining, and grouped summaries in dplyr reach a
#' column, is the usual route, and dropping the record there raises a warning
#' of class `propensity_trunc_record_warning`. Combining with [c()] drops it
#' without comment, because concatenation appends one set of observations to
#' another and the prototype it builds the result from holds no positions to
#' lose. The values, the class, and the method and its bounds are untouched
#' either way.
#'
#' A record can also outlive the observations it describes, because it travels
#' by routes vctrs does not see: growing a `ps_trunc` by subassignment carries
#' it across the length change. [is_unit_truncated()] therefore checks that the
#' record covers the object it is given and raises an error of class
#' `propensity_missing_meta_error` when it does not, rather than name truncated
#' units at stale positions.
#'
#' That check compares how many observations the record was written for against
#' how many the object holds, which a reordering does not change. An operation
#' that reorders the observations through vctrs, rather than through `[`,
#' therefore keeps a record written for the order they used to be in:
#' `vctrs::vec_slice(x, 5:1)` and `dplyr::arrange()` both return the values in
#' a new order under positions still naming the old one, and
#' [is_unit_truncated()] answers from those positions and names the wrong
#' units. Subsetting with `[` is handed the subscript and re-indexes, so
#' reorder with `[`, or put the propensity scores in the order you want before
#' truncating them.
#'
#' @return A `ps_trunc` object (a numeric vector for binary exposures, or a
#'   matrix for categorical exposures). Use [ps_trunc_meta()] to inspect
#'   metadata including `method`, `lower_bound`, `upper_bound`,
#'   `truncated_idx` (positions of modified values), and `n_obs` (the number of
#'   observations those positions describe).
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
      truncated_idx = truncated_idx,
      n_obs = length(ps)
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
    n_obs = n,
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

# The positional half of a truncation record. The rest of it, the method and its
# bounds, describes the truncation rather than the units and means the same
# thing at any length.
reindex_trunc_record <- function(meta, i) {
  meta$truncated_idx <- reindex_positions(meta$truncated_idx, i)
  meta$n_obs <- length(i)

  meta
}

drop_trunc_record <- function(meta) {
  meta[["truncated_idx"]] <- NULL
  meta[["n_obs"]] <- NULL

  meta
}

# A record that no longer describes the observations in front of it is dropped
# rather than guessed at. Nothing in the values says which units a shorter or a
# longer vector once had winsorized: a score that arrived equal to a bound is
# indistinguishable from one this package pinned there.
#
# A record over no observations names no unit, so replacing it costs the caller
# nothing and goes without comment. That is the record every prototype carries,
# which is what concatenation builds its result from.
discard_trunc_record <- function(meta, n) {
  recorded <- meta$n_obs

  if (!is.null(recorded) && recorded > 0) {
    warn(
      c(
        "Dropping the record of which units were truncated.",
        i = "The record describes {recorded} observation{?s} and this result
             has {n}, so its positions do not describe them.",
        i = "The values are unchanged and the result is still a
             {.cls ps_trunc}. Truncate the propensity scores you want to work
             with to get a record written for them."
      ),
      warning_class = "propensity_trunc_record_warning",
      # One of the routes here is vctrs' internal dispatch, whose call would be
      # reported and names nothing the caller wrote, so no call is attributed.
      call = NULL
    )
  }

  drop_trunc_record(meta)
}

# `i` holds the old positions the result is built from, in the order it holds
# them, among the `n_obs` observations it was taken from. A record that covers
# those observations re-indexes onto `i`; one that does not cannot be placed
# against them at all.
carry_trunc_record <- function(meta, n_obs, i) {
  if (record_covers(meta, n_obs)) {
    reindex_trunc_record(meta, i)
  } else {
    discard_trunc_record(meta, length(i))
  }
}

# A positional query reads its answer out of the record, so a record that does
# not cover the object in front of it has no answer to give: reporting every
# unit as untouched would be a wrong answer rather than a missing one.
check_trunc_record <- function(meta, n, fn, call = rlang::caller_env()) {
  if (record_covers(meta, n)) {
    return(invisible(meta))
  }

  recorded <- meta$n_obs
  problem <- if (is.null(recorded)) {
    "These propensity scores carry no record of which units were truncated."
  } else {
    "The record describes {recorded} observation{?s} and these propensity
     scores have {n}, so its positions do not describe them."
  }

  abort(
    c(
      "{.code {fn}} has no usable truncation record for these propensity
       scores.",
      x = problem,
      i = "An operation that changed the number of observations dropped it.
           Truncate the propensity scores you want to work with, or query the
           object the record was written for."
    ),
    error_class = "propensity_missing_meta_error",
    call = call
  )
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
#'   * `n_obs` -- the number of observations those positions describe
#'
#'   `truncated_idx` and `n_obs` are absent when an operation that changed the
#'   number of observations dropped the record; see [ps_trunc()].
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
#'   The answer comes from the truncation record, which is written for a fixed
#'   set of observations and can both be lost and outlive them:
#'   [vctrs::vec_slice()] and [c()] drop it, and subassignment that grows the
#'   vector carries it across a length change. `is_unit_truncated()` therefore
#'   checks that the record covers the object it is given, and raises an error
#'   of class `propensity_missing_meta_error` when it does not, rather than
#'   name truncated units at stale positions.
#'
#'   That check counts observations, which a reordering does not change, so it
#'   does not catch one. An operation that reorders through vctrs rather than
#'   through `[`, such as `vctrs::vec_slice(x, 5:1)` or `dplyr::arrange()`,
#'   keeps a record written for the old order, and `is_unit_truncated()`
#'   answers from it and names the wrong units. See [ps_trunc()] for the whole
#'   contract.
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
  # No observations, no answers. A record kept on an empty vector describes
  # observations it does not have, and indexing an empty logical by the
  # positions it names would grow one padded with `NA`.
  if (length(x) == 0) {
    return(logical(0))
  }

  meta <- ps_trunc_meta(x)
  check_trunc_record(meta, length(x), "is_unit_truncated()")

  out <- vector("logical", length = length(x))
  out[meta$truncated_idx] <- TRUE

  out
}

#' @export
is_unit_truncated.ps_trunc_matrix <- function(x) {
  if (nrow(x) == 0) {
    return(logical(0))
  }

  meta <- ps_trunc_meta(x)
  check_trunc_record(meta, nrow(x), "is_unit_truncated()")

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

  # The rows the result was built from, as positions in `x`, which is what
  # re-indexing the record takes. These are worked out from the subscript rather
  # than read off the result, so a form they do not reproduce, seen as a count
  # that does not match the rows base returned, degrades the record rather than
  # certifying it against rows it does not describe.
  rows <- subscript_row_positions(i, x)
  new_meta <- if (length(rows) == nrow(result)) {
    carry_trunc_record(meta, nrow(x), rows)
  } else {
    discard_trunc_record(meta, nrow(result))
  }

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
    meta = list(
      method = "unknown",
      lower_bound = NA,
      upper_bound = NA,
      truncated_idx = integer(0),
      n_obs = length(x)
    )
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
    meta = list(
      method = "unknown",
      lower_bound = NA,
      upper_bound = NA,
      truncated_idx = integer(0),
      n_obs = length(xx)
    )
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

  # The subscript is in hand here, which is what re-indexing a record takes and
  # what the restore behind `NextMethod()` lacks: that restore would see a
  # record written for a different number of observations and drop it. Building
  # the result here keeps ordinary subsetting describing what it returns.
  new_ps_trunc(vec_data(x)[i], carry_trunc_record(meta, length(x), i))
}

#' @export
unique.ps_trunc <- function(x, incomparables = FALSE, ...) {
  check_incomparables(incomparables, "ps_trunc")

  # `vec_unique_loc()` names the position each retained value came from, which
  # is the subscript re-indexing the record takes. Without this the restore
  # behind vctrs' own method sees only a shorter vector and drops the record.
  x[vec_unique_loc(x)]
}

#' @export
rep.ps_trunc <- function(x, ...) {
  # The positions `rep()` produces are the subscript the result is built from,
  # so the record follows them and a repeated unit is reported at each place it
  # now holds.
  x[rep(seq_along(x), ...)]
}

#' @export
summary.ps_trunc <- function(object, ...) {
  summary(vec_data(object), ...)
}

#' @export
sort.ps_trunc <- function(x, decreasing = FALSE, na.last = NA, ...) {
  meta <- ps_trunc_meta(x)
  x_data <- vec_data(x)

  # `order()` returns the old positions in their new order, so it is the
  # subscript the result is built from and the record follows it.
  ord <- order(x_data, na.last = na.last, decreasing = decreasing, ...)

  new_ps_trunc(x_data[ord], carry_trunc_record(meta, length(x), ord))
}

#' @export
vec_restore.ps_trunc <- function(x, to, ...) {
  # vec_data in case x is already a vctr
  data <- vec_data(x)
  meta <- ps_trunc_meta(to)

  # Nothing rebuilding a `ps_trunc` is handed the subscript behind a length
  # change, so a record written for a different number of observations cannot be
  # re-indexed onto the data arriving here. Zero-length data is exempt: a
  # prototype or an empty slice holds no observations, so no position in the
  # record contradicts it, and the record rides along to the restore that builds
  # the real result.
  if (length(data) > 0 && !record_covers(meta, length(data))) {
    meta <- discard_trunc_record(meta, length(data))
  }

  new_ps_trunc(data, meta)
}

#' @export
anyDuplicated.ps_trunc <- function(x, incomparables = FALSE, ...) {
  anyDuplicated(vec_data(x), incomparables = incomparables, ...)
}

#' @export
diff.ps_trunc <- function(x, lag = 1L, differences = 1L, ...) {
  diff(vec_data(x), lag = lag, differences = differences, ...)
}
