abort_unsupported <- function(exposure_type, what, call = rlang::caller_env()) {
  abort(
    "Exposure type {.val {exposure_type}} not currently supported for {.field {what}}",
    call = call,
    error_class = "propensity_wt_not_supported_error"
  )
}

abort_no_method <- function(.propensity, call = rlang::caller_env()) {
  abort(
    paste0(
      "No method for objects of class ",
      paste(class(.propensity), collapse = ", ")
    ),
    call = call,
    error_class = "propensity_method_error"
  )
}

match_exposure_type <- function(
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .exposure,
  valid_types = c("auto", "binary", "categorical", "continuous"),
  call = rlang::caller_env()
) {
  .exposure_type <- rlang::arg_match(
    exposure_type,
    valid_types,
    error_call = call
  )
  if (.exposure_type == "auto") {
    detect_exposure_type(.exposure)
  } else {
    .exposure_type
  }
}

detect_exposure_type <- function(.exposure) {
  exposure_type <- if (has_two_levels(.exposure)) {
    "binary"
  } else if (is.factor(.exposure) || is.character(.exposure)) {
    # Check number of unique values for factor/character
    if (length(unique(.exposure)) > 2) {
      "categorical"
    } else {
      "binary"
    }
  } else if (is_categorical(.exposure)) {
    "categorical"
  } else {
    "continuous"
  }

  alert_info("Treating {.arg .exposure} as {exposure_type}")

  exposure_type
}

handle_focal_deprecation <- function(
  .focal_level,
  .reference_level,
  .treated,
  .untreated,
  fn_name
) {
  # Handle deprecation warnings and parameter mapping
  if (!is.null(.treated)) {
    lifecycle::deprecate_warn(
      "0.1.0",
      paste0(fn_name, "(.treated)"),
      paste0(fn_name, "(.focal_level)")
    )
    if (is.null(.focal_level)) {
      .focal_level <- .treated
    }
  }

  if (!is.null(.untreated)) {
    lifecycle::deprecate_warn(
      "0.1.0",
      paste0(fn_name, "(.untreated)"),
      paste0(fn_name, "(.reference_level)")
    )
    if (is.null(.reference_level)) {
      .reference_level <- .untreated
    }
  }

  list(.focal_level = .focal_level, .reference_level = .reference_level)
}

transform_exposure_binary <- function(
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  if (!is.null(.focal_level)) {
    return(ifelse(.exposure == .focal_level, 1, 0))
  }

  if (!is.null(.reference_level)) {
    return(ifelse(.exposure != .reference_level, 1, 0))
  }

  # With no level named, an exposure that already carries the 0/1 coding, or a
  # logical one, is taken at face value and needs no announcement.
  if (is_binary(.exposure)) {
    return(.exposure)
  }

  if (is.logical(.exposure)) {
    return(as.numeric(.exposure))
  }

  if (has_two_levels(.exposure)) {
    levels <- if (is.factor(.exposure)) {
      levels(.exposure)
    } else {
      sort(unique(.exposure))
    }
    alert_info("Setting focal level to {.val {levels[[2]]}}")
    return(ifelse(.exposure == levels[[2]], 1, 0))
  } else {
    abort(
      c(
        "Don't know how to transform `.exposure` to 0/1 binary variable.",
        i = "Specify `.focal_level` and `.reference_level`."
      ),
      call = call,
      error_class = "propensity_binary_transform_error"
    )
  }
}

# The exposure level `transform_exposure_binary()` actually codes as 1. A named
# level always wins; only when neither is named does the coding fall back to the
# exposure's own ordering, which is 1 for a 0/1 exposure and `TRUE` for a
# logical one. Mirrors that function's branching in the same order.
effective_binary_focal_level <- function(
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  if (!is.null(.focal_level)) {
    return(.focal_level)
  }

  if (!is.null(.reference_level)) {
    # Every level other than the reference is coded as focal, so a single focal
    # level is named only when the exposure takes exactly one other value.
    # Missing values are dropped first: `!=` returns NA rather than FALSE for
    # them, so an NA would survive the subscript and count as a second remaining
    # value, leaving a two-level exposure with no focal level recorded.
    observed <- .exposure[!is.na(.exposure)]
    remaining <- unique(observed[observed != .reference_level])
    if (length(remaining) != 1) {
      return(NULL)
    }
    return(as_focal_label(remaining))
  }

  if (is_binary(.exposure)) {
    return(1)
  }

  if (is.logical(.exposure)) {
    return(TRUE)
  }

  if (!has_two_levels(.exposure)) {
    return(NULL)
  }

  if (is.factor(.exposure)) {
    levels(.exposure)[[2]]
  } else {
    sort(unique(.exposure))[[2]]
  }
}

# Subsetting a factor yields a factor; record its label, which is what the
# categorical path stores and what comparisons against exposure values expect.
as_focal_label <- function(x) {
  if (is.factor(x)) as.character(x) else x
}

# Whether the named levels leave the focal level where the default coding puts
# it. Both sides come from `effective_binary_focal_level()`, so this agrees with
# `transform_exposure_binary()` by construction. A level that cannot be resolved
# counts as the default, leaving the exposure to be handled downstream as it is
# today. Note that an NA-bearing factor or character exposure short-circuits
# this to TRUE regardless of the argument: `has_two_levels()` counts NA as a
# level, so the default focal level is NULL and the comparison below never runs.
resolved_focal_is_default <- function(
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  resolved <- effective_binary_focal_level(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )
  default <- effective_binary_focal_level(.exposure)

  if (is.null(resolved) || is.null(default)) {
    return(TRUE)
  }

  # `==` rather than `identical()`, matching how `transform_exposure_binary()`
  # compares the exposure against a named level: 1 and 1L are the same level.
  isTRUE(resolved == default)
}

# Binary att and atu weights are mirror images of one another: att weights built
# on one exposure level are numerically identical to atu weights built on the
# other, so nothing in the values records which level they target. Store the
# level under the attribute name the categorical path uses, but only when the
# caller named one; with no level named there is no intent to check against.
record_binary_focal_level <- function(
  psw_obj,
  .exposure,
  exposure_type,
  .focal_level = NULL,
  .reference_level = NULL
) {
  level_supplied <- !is.null(.focal_level) || !is.null(.reference_level)
  if (!identical(exposure_type, "binary") || !level_supplied) {
    return(psw_obj)
  }

  attr(psw_obj, "focal_category") <- effective_binary_focal_level(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )

  psw_obj
}

is_binary <- function(.exposure) {
  is.numeric(.exposure) &&
    identical(sort(unique(as.double(.exposure))), c(0, 1))
}

is_categorical <- function(.exposure) {
  # assumption: a variable where the proportion of unique values
  # to total number of observations is less than 20% is categorical
  n_non_na <- sum(!is.na(.exposure))
  if (n_non_na == 0) {
    return(FALSE)
  }

  ratio <- length(unique(.exposure)) / n_non_na
  # Handle NaN case explicitly
  if (is.nan(ratio)) {
    return(FALSE)
  }

  ratio < 0.2
}

has_two_levels <- function(.x) {
  length(unique(.x)) == 2
}

check_refit <- function(.propensity, call = rlang::caller_env()) {
  if (!is_refit(.propensity)) {
    warn(
      c(
        "It appears you trimmed your propensity score but did not refit the model.",
        i = "Use {.code ps_refit()} for more accurate re-estimation."
      ),
      warning_class = "propensity_no_refit_warning",
      call = call
    )
  }
}

check_ps_range <- function(ps, call = rlang::caller_env()) {
  if (is.matrix(ps) || is.data.frame(ps)) {
    # For matrices/data frames, check all values
    ps_vals <- as.numeric(as.matrix(ps))
    # Check only non-NA values
    non_na_vals <- ps_vals[!is.na(ps_vals)]
    if (
      length(non_na_vals) > 0 &&
        any(non_na_vals <= 0 | non_na_vals >= 1 | !is.finite(non_na_vals))
    ) {
      abort(
        c(
          "All propensity scores must be between 0 and 1.",
          i = "The range of values in {.arg ps} is \\
        {format(range(ps_vals, na.rm = TRUE), nsmall = 1, digits = 1)}"
        ),
        call = call,
        error_class = "propensity_range_error"
      )
    }
  } else {
    ps <- as.numeric(ps)
    # Check only non-NA values
    non_na_vals <- ps[!is.na(ps)]
    if (
      length(non_na_vals) > 0 &&
        any(non_na_vals <= 0 | non_na_vals >= 1 | !is.finite(non_na_vals))
    ) {
      abort(
        c(
          "The propensity score must be between 0 and 1.",
          i = "The range of {.arg ps} is \\
        {format(range(ps, na.rm = TRUE), nsmall = 1, digits = 1)}"
        ),
        call = call,
        error_class = "propensity_range_error"
      )
    }
  }

  invisible(TRUE)
}

check_lower_upper <- function(lower, upper, call = rlang::caller_env()) {
  if (lower >= upper) {
    abort(
      c(
        "{.arg lower} must be smaller than {.arg upper}",
        x = "{.arg lower} is {lower} and {.arg upper} is {upper}"
      ),
      call = call,
      error_class = "propensity_range_error"
    )
  }

  invisible(TRUE)
}

check_lengths_match <- function(
  .propensity,
  .exposure,
  call = rlang::caller_env()
) {
  # Handle matrix/data.frame inputs
  if (is.matrix(.propensity) || is.data.frame(.propensity)) {
    len_prop <- nrow(.propensity)
  } else {
    len_prop <- length(.propensity)
  }

  len_exp <- length(.exposure)

  if (len_prop != len_exp) {
    abort(
      c(
        "{.arg .propensity} and {.arg .exposure} must have the same length.",
        i = "{.arg .propensity} has {if (is.matrix(.propensity) || is.data.frame(.propensity)) 'rows' else 'length'} {len_prop}",
        i = "{.arg .exposure} has length {len_exp}"
      ),
      call = call,
      error_class = "propensity_length_error"
    )
  }

  invisible(TRUE)
}

transform_exposure_categorical <- function(
  .exposure,
  .focal_level = NULL,
  call = rlang::caller_env()
) {
  # Convert to factor if not already
  if (!is.factor(.exposure)) {
    .exposure <- as.factor(.exposure)
  }

  # Check if we have more than 2 levels
  n_levels <- nlevels(.exposure)
  if (n_levels <= 2) {
    abort(
      c(
        "Categorical exposure must have more than 2 levels.",
        i = "Found {n_levels} levels.",
        i = "Use binary exposure methods for 2-level exposures."
      ),
      call = call,
      error_class = "propensity_categorical_levels_error"
    )
  }

  # Validate focal category if provided
  if (!is.null(.focal_level)) {
    if (!.focal_level %in% levels(.exposure)) {
      abort(
        c(
          "Focal category must be one of the exposure levels.",
          i = "Focal category: {.val {(.focal_level)}}",
          i = "Available levels: {.val {levels(.exposure)}}"
        ),
        call = call,
        error_class = "propensity_focal_category_error"
      )
    }
  }

  .exposure
}

check_ps_matrix <- function(
  ps_matrix,
  .exposure,
  call = rlang::caller_env()
) {
  # Convert to matrix if data frame first
  if (is.data.frame(ps_matrix)) {
    ps_matrix <- as.matrix(ps_matrix)
  }

  # Check if it's a matrix
  if (!is.matrix(ps_matrix)) {
    abort(
      "For categorical exposures, `.propensity` must be a matrix or data frame.",
      call = call,
      error_class = "propensity_matrix_type_error"
    )
  }

  # Check dimensions
  n_obs <- length(.exposure)
  n_cats <- nlevels(.exposure)

  if (nrow(ps_matrix) != n_obs) {
    abort(
      c(
        "Number of rows in propensity score matrix must match number of observations.",
        i = "Matrix rows: {nrow(ps_matrix)}",
        i = "Observations: {n_obs}"
      ),
      call = call,
      error_class = "propensity_matrix_dims_error"
    )
  }

  if (ncol(ps_matrix) != n_cats) {
    abort(
      c(
        "Number of columns in propensity score matrix must match number of exposure categories.",
        i = "Matrix columns: {ncol(ps_matrix)}",
        i = "Categories: {n_cats}"
      ),
      call = call,
      error_class = "propensity_matrix_dims_error"
    )
  }

  # Check that rows sum to 1 (within tolerance)
  # Only check non-NA rows
  row_sums <- rowSums(ps_matrix, na.rm = FALSE)
  ROW_SUM_TOLERANCE <- 1e-6 # Tolerance for floating point comparison
  non_na_rows <- !is.na(row_sums)

  if (any(non_na_rows)) {
    # Check only the rows that don't have NA values
    if (any(abs(row_sums[non_na_rows] - 1) > ROW_SUM_TOLERANCE)) {
      bad_rows <- which(abs(row_sums - 1) > ROW_SUM_TOLERANCE & non_na_rows)
      abort(
        c(
          "Propensity score matrix rows must sum to 1.",
          i = "Problem rows: {bad_rows[1:min(5, length(bad_rows))]}{if (length(bad_rows) > 5) ' ...' else ''}"
        ),
        call = call,
        error_class = "propensity_matrix_sum_error"
      )
    }
  }

  # Check for valid probabilities
  if (any(ps_matrix < 0 | ps_matrix > 1, na.rm = TRUE)) {
    abort(
      "All propensity scores must be between 0 and 1.",
      call = call,
      error_class = "propensity_range_error"
    )
  }

  # Ensure columns are in the same order as factor levels
  # This is critical for correct weight calculation
  exp_levels <- levels(.exposure)

  # Check if columns have names
  if (!is.null(colnames(ps_matrix))) {
    # Try to match column names to factor levels
    # Handle both plain names (A, B, C) and parsnip-style names (.pred_A, .pred_B, .pred_C)
    col_names <- colnames(ps_matrix)

    # Remove common prefixes like ".pred_"
    clean_names <- gsub("^\\.pred_", "", col_names)

    # Check if clean names match factor levels
    if (setequal(clean_names, exp_levels)) {
      # Reorder columns to match factor levels
      if (!identical(clean_names, exp_levels)) {
        col_order <- match(exp_levels, clean_names)
        ps_matrix <- ps_matrix[, col_order, drop = FALSE]
        # Update column names to match
        colnames(ps_matrix) <- col_names[col_order]
      }
    } else {
      # Column names don't match factor levels
      abort(
        c(
          "Column names of propensity score matrix must match exposure levels.",
          i = "Expected levels: {.val {exp_levels}}",
          i = "Found columns: {.val {clean_names}}"
        ),
        call = call,
        error_class = "propensity_matrix_names_error"
      )
    }
  } else {
    # No column names - assume they're in factor level order
    # Issue a warning as this is risky
    warn(
      c(
        "Propensity score matrix has no column names.",
        i = "Assuming columns are in factor level order: {.val {exp_levels}}",
        i = "This may lead to incorrect results if columns are misaligned."
      ),
      warning_class = "propensity_matrix_no_names_warning",
      call = call
    )
  }

  ps_matrix
}

# Helper for ps_trim and ps_trunc methods
calculate_weight_from_modified_ps <- function(
  .propensity,
  .exposure,
  weight_fn,
  modification_type = c("trim", "trunc", "calib"),
  ...
) {
  modification_type <- rlang::arg_match(modification_type)

  # `weight_fn` is a numeric method invoked through a local symbol, so its own
  # frame deparses as `weight_fn()`. Bind the dispatching method's frame here
  # and hand it down, so conditions raised inside the weight machinery name the
  # weight function the user called.
  call <- rlang::caller_env()

  # Only check refit for trim
  if (modification_type == "trim") {
    check_refit(.propensity, call = call)
  }

  # Handle matrix or vector propensity scores
  if (
    inherits(
      .propensity,
      c("ps_trim_matrix", "ps_trunc_matrix", "ps_calib_matrix")
    )
  ) {
    # For matrix propensity scores, pass them directly
    # The weight function should handle the matrix appropriately
    base_wt <- weight_fn(
      .propensity,
      .exposure = .exposure,
      call = call,
      ...
    )
  } else {
    # Convert to numeric for vector propensity scores
    numeric_ps <- as.numeric(.propensity)

    # Call the weight function with the numeric propensity scores
    base_wt <- weight_fn(
      numeric_ps,
      .exposure = .exposure,
      call = call,
      ...
    )
  }

  # Update estimand
  if (modification_type == "trim") {
    old_est <- estimand(base_wt)
    estimand(base_wt) <- paste0(old_est, "; trimmed")
    attr(base_wt, "trimmed") <- TRUE
    attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")
  } else if (modification_type == "trunc") {
    estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")
    attr(base_wt, "truncated") <- TRUE
    attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)
  } else {
    # calib
    estimand(base_wt) <- paste0(estimand(base_wt), "; calibrated")
    attr(base_wt, "calibrated") <- TRUE
    attr(base_wt, "ps_calib_meta") <- ps_calib_meta(.propensity)
  }

  base_wt
}

# Helper to preserve categorical attributes on psw objects
preserve_categorical_attrs <- function(psw_obj, wts, exposure_type) {
  if (exposure_type == "categorical") {
    attr(psw_obj, "n_categories") <- attr(wts, "n_categories")
    attr(psw_obj, "category_names") <- attr(wts, "category_names")
    # focal_category might not always exist
    if (!is.null(attr(wts, "focal_category"))) {
      attr(psw_obj, "focal_category") <- attr(wts, "focal_category")
    }
  }
  psw_obj
}

# The exposure levels `transform_exposure_binary()` reads, derived the same way
# it derives them so that a column named for a level is named for the level the
# coding actually uses.
binary_exposure_levels <- function(.exposure) {
  if (is.factor(.exposure)) {
    levels(.exposure)
  } else {
    sort(unique(.exposure))
  }
}

# Position of the column holding the probability of the resolved focal level,
# or NULL when the column names cannot answer the question. The names must
# cover every exposure level before any of them is trusted, so that a lone
# coincidental match cannot redirect the selection. A reference level that
# leaves more than one candidate resolves to no focal level at all, which is
# likewise unanswerable.
match_focal_level_column <- function(
  .propensity,
  .exposure,
  .focal_level = NULL,
  .reference_level = NULL
) {
  exposure_levels <- as.character(binary_exposure_levels(.exposure))
  if (!all(exposure_levels %in% names(.propensity))) {
    return(NULL)
  }

  focal_level <- effective_binary_focal_level(
    .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level
  )
  if (is.null(focal_level)) {
    return(NULL)
  }

  col_pos <- match(as.character(focal_level), names(.propensity))
  if (is.na(col_pos)) NULL else col_pos
}

# Helper function to extract propensity scores from data frames
# This consolidates the logic used across multiple weight functions
extract_propensity_from_df <- function(
  .propensity,
  .propensity_col_quo = NULL,
  .exposure = NULL,
  exposure_type = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  call = rlang::caller_env()
) {
  if (!rlang::quo_is_null(.propensity_col_quo)) {
    col_pos <- tryCatch(
      tidyselect::eval_select(
        .propensity_col_quo,
        data = .propensity
      ),
      error = function(e) {
        abort(
          paste0("Column selection failed: ", e$message),
          call = call,
          error_class = "propensity_df_column_error"
        )
      }
    )

    if (length(col_pos) != 1) {
      abort(
        "`.propensity_col` must select exactly one column.",
        call = call,
        error_class = "propensity_df_column_error"
      )
    }

    return(.propensity[[col_pos]])
  }

  if (ncol(.propensity) == 0) {
    abort(
      "`.propensity` data frame must have at least one column.",
      call = call,
      error_class = "propensity_df_ncol_error"
    )
  }

  is_binary_exposure <- identical(exposure_type, "binary")

  if (is_binary_exposure) {
    col_pos <- match_focal_level_column(
      .propensity,
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
    if (!is.null(col_pos)) {
      return(.propensity[[col_pos]])
    }
  }

  # Default behavior: use second column if available, otherwise first
  col_pos <- if (ncol(.propensity) >= 2) 2L else 1L

  # A single column is the caller supplying the probability of the focal level
  # directly, so there is nothing to choose between and nothing to report.
  level_named <- !is.null(.focal_level) || !is.null(.reference_level)
  if (is_binary_exposure && level_named && ncol(.propensity) > 1) {
    warn(
      c(
        "Can't tell which column of {.arg .propensity} holds the probability of the focal level.",
        i = "Selected {.val {names(.propensity)[[col_pos]]}} by position.",
        i = "Name the columns after the levels of {.arg .exposure}, or set {.arg .propensity_col}, to select the column by name."
      ),
      warning_class = "propensity_df_column_warning",
      call = call
    )
  }

  .propensity[[col_pos]]
}

# Helper function to extract propensity scores from GLM objects
extract_propensity_from_glm <- function(
  .propensity,
  call = rlang::caller_env()
) {
  # Check if it's a valid GLM object
  if (!inherits(.propensity, "glm")) {
    abort(
      "`.propensity` must be a GLM object.",
      call = call,
      error_class = "propensity_glm_type_error"
    )
  }

  # Check if it's a binomial GLM for binary propensity scores
  if (
    !is.null(.propensity$family) &&
      .propensity$family$family == "binomial"
  ) {
    # Get predicted probabilities
    ps_vec <- stats::predict(.propensity, type = "response")
  } else {
    # For non-binomial GLMs, get linear predictor
    ps_vec <- stats::fitted(.propensity)
  }

  ps_vec
}

# Helper function to handle common data frame method pattern
# This encapsulates the logic used across all weight function data.frame methods
#
# The deprecated arguments are resolved here rather than downstream, for the
# same reason as in `prepare_glm_weight_args()`: the resolved focal level
# decides which column carries `.propensity`, and mapping it twice would emit
# the deprecation warning twice, so the numeric method receives the mapped
# levels and no deprecated arguments.
handle_data_frame_weight_calculation <- function(
  weight_fn_numeric,
  .propensity,
  .exposure,
  exposure_type,
  valid_exposure_types = c("auto", "binary", "categorical", "continuous"),
  .propensity_col_quo,
  .focal_level = NULL,
  .reference_level = NULL,
  .treated = NULL,
  .untreated = NULL,
  fn_name,
  ...
) {
  # Validate inputs
  if (!is.data.frame(.propensity)) {
    abort(
      "`.propensity` must be a data frame.",
      call = rlang::caller_env(2),
      error_class = "propensity_matrix_type_error"
    )
  }

  # Check exposure type
  exposure_type_check <- match_exposure_type(
    exposure_type,
    .exposure,
    valid_exposure_types,
    call = rlang::caller_env()
  )

  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    fn_name
  )

  if (exposure_type_check == "categorical") {
    # For categorical exposures, pass the whole data frame
    return(weight_fn_numeric(
      .propensity = .propensity,
      .exposure = .exposure,
      exposure_type = exposure_type,
      .focal_level = focal_params$.focal_level,
      .reference_level = focal_params$.reference_level,
      ...
    ))
  }

  # For non-categorical exposures, extract single column
  ps_vec <- extract_propensity_from_df(
    .propensity,
    .propensity_col_quo,
    .exposure = .exposure,
    exposure_type = exposure_type_check,
    .focal_level = focal_params$.focal_level,
    .reference_level = focal_params$.reference_level,
    call = rlang::caller_env(2)
  )

  # Call the numeric method
  weight_fn_numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .focal_level = focal_params$.focal_level,
    .reference_level = focal_params$.reference_level,
    ...
  )
}

# Helper functions for extracting information from GLM formulas
# (moved from ipw.R to be shared across the package)
fmla_extract_left_vctr <- function(mod) {
  .data <- mod |>
    model.frame()

  .data[[1]]
}

fmla_extract_left_chr <- function(mod) {
  as.character(formula(mod)[[2]])
}

# Variables a model's left-hand side reads. `fmla_extract_left_chr()` deparses
# the whole expression, so a transformed response such as `log(y)` gives back the
# function name alongside its argument; a rebuild from a supplied data frame needs
# the columns the response is computed from, which is what `all.vars()` gives and
# what the user can actually supply.
fmla_extract_left_vars <- function(mod) {
  all.vars(formula(mod)[[2]])
}

# Evaluate a model's left-hand side against a data frame. A response written as a
# transformation has to be computed rather than looked up by name, and the
# functions and constants it reads resolve from the formula's environment, the
# same way they did when the model was fit. `scale()` and friends return a
# one-column matrix; drop that dimension so the result matches what
# stats::model.response() hands back on the model-frame route.
fmla_eval_left <- function(mod, data) {
  fmla <- formula(mod)
  left <- eval(fmla[[2]], data, environment(fmla))

  if (is.matrix(left) && ncol(left) == 1L) {
    left <- drop(left)
  }

  left
}

# Helper function to handle optional exposure in GLM methods
extract_exposure_from_glm <- function(
  glm_obj,
  .exposure = NULL
) {
  if (is.null(.exposure)) {
    # Extract exposure from GLM
    .exposure <- fmla_extract_left_vctr(glm_obj)
    exposure_name <- fmla_extract_left_chr(glm_obj)
    alert_info("Using exposure variable {.val {exposure_name}} from GLM model")
  }
  .exposure
}

# Shared preparation for the GLM methods of the weight functions.
#
# Fitted values report the probability of the level the response's default
# coding treats as focal: the second factor level, or the second sorted value.
# The weight formulas read `.propensity` as the probability of the level
# actually resolved as focal, so naming the other level means the fitted values
# must be inverted before the numeric method sees them.
#
# The deprecated arguments are resolved here rather than downstream. The
# resolved focal level is needed to decide whether to invert, and mapping it
# twice would emit the deprecation warning twice, so the numeric method
# receives the mapped levels and no deprecated arguments.
prepare_glm_weight_args <- function(
  glm_obj,
  .exposure,
  exposure_type,
  valid_exposure_types,
  .focal_level,
  .reference_level,
  .treated,
  .untreated,
  fn_name,
  call = rlang::caller_env()
) {
  .exposure <- extract_exposure_from_glm(glm_obj, .exposure)
  exposure_type <- match_exposure_type(
    exposure_type,
    .exposure,
    valid_exposure_types,
    call = call
  )

  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    fn_name
  )

  ps_vec <- extract_propensity_from_glm(glm_obj, call = call)

  invert <- identical(exposure_type, "binary") &&
    !resolved_focal_is_default(
      .exposure,
      .focal_level = focal_params$.focal_level,
      .reference_level = focal_params$.reference_level
    )

  if (invert) {
    ps_vec <- 1 - ps_vec
  }

  list(
    propensity = ps_vec,
    exposure = .exposure,
    exposure_type = exposure_type,
    focal_level = focal_params$.focal_level,
    reference_level = focal_params$.reference_level
  )
}
