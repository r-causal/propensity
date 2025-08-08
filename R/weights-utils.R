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
  valid_types = c("auto", "binary", "categorical", "continuous")
) {
  .exposure_type <- rlang::arg_match(exposure_type, valid_types)
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

  alert_info("Treating `.exposure` as {exposure_type}")

  exposure_type
}

transform_exposure_binary <- function(
  .exposure,
  .treated = NULL,
  .untreated = NULL,
  call = rlang::caller_env()
) {
  if (is_binary(.exposure)) {
    return(.exposure)
  }

  if (is.logical(.exposure)) {
    return(as.numeric(.exposure))
  }

  if (!is.null(.treated)) {
    return(ifelse(.exposure == .treated, 1, 0))
  }

  if (!is.null(.untreated)) {
    return(ifelse(.exposure != .untreated, 1, 0))
  }

  if (is.null(.treated) && is.null(.untreated) && has_two_levels(.exposure)) {
    levels <- if (is.factor(.exposure)) levels(.exposure) else
      sort(unique(.exposure))
    alert_info("Setting treatment to {.var {levels[[2]]}}")
    return(ifelse(.exposure == levels[[2]], 1, 0))
  } else {
    abort(
      c(
        "Don't know how to transform `.exposure` to 0/1 binary variable.",
        i = "Specify `.treated` and `.untreated.`"
      ),
      call = call,
      error_class = "propensity_binary_transform_error"
    )
  }
}

is_binary <- function(.exposure) {
  identical(sort(unique(.exposure)), c(0, 1))
}

is_categorical <- function(.exposure) {
  # assumption: a variable where the proportion of unique values
  # to total number of observations is less than 20% is categorical
  n_non_na <- sum(!is.na(.exposure))
  if (n_non_na == 0) return(FALSE)

  ratio <- length(unique(.exposure)) / n_non_na
  # Handle NaN case explicitly
  if (is.nan(ratio)) return(FALSE)

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
  focal = NULL,
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
  if (!is.null(focal)) {
    if (!focal %in% levels(.exposure)) {
      abort(
        c(
          "Focal category must be one of the exposure levels.",
          i = "Focal category: {.val {focal}}",
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
  modification_type = c("trim", "trunc"),
  ...
) {
  modification_type <- rlang::arg_match(modification_type)

  # Only check refit for trim
  if (modification_type == "trim") {
    check_refit(.propensity, call = rlang::caller_env())
  }

  # Handle matrix or vector propensity scores
  if (inherits(.propensity, c("ps_trim_matrix", "ps_trunc_matrix"))) {
    # For matrix propensity scores, pass them directly
    # The weight function should handle the matrix appropriately
    base_wt <- weight_fn(
      .propensity,
      .exposure = .exposure,
      ...
    )
  } else {
    # Convert to numeric for vector propensity scores
    numeric_ps <- as.numeric(.propensity)

    # Call the weight function with the numeric propensity scores
    base_wt <- weight_fn(
      numeric_ps,
      .exposure = .exposure,
      ...
    )
  }

  # Update estimand
  if (modification_type == "trim") {
    old_est <- estimand(base_wt)
    estimand(base_wt) <- paste0(old_est, "; trimmed")
    attr(base_wt, "trimmed") <- TRUE
    attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")
  } else {
    estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")
    attr(base_wt, "truncated") <- TRUE
    attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)
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

# Helper function to extract propensity scores from data frames
# This consolidates the logic used across multiple weight functions
extract_propensity_from_df <- function(
  .propensity,
  .propensity_col_quo = NULL,
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

    ps_vec <- .propensity[[col_pos]]
  } else {
    # Default behavior: use second column if available, otherwise first
    if (ncol(.propensity) >= 2) {
      ps_vec <- .propensity[[2]]
    } else if (ncol(.propensity) == 1) {
      ps_vec <- .propensity[[1]]
    } else {
      abort(
        "`.propensity` data frame must have at least one column.",
        call = call,
        error_class = "propensity_df_ncol_error"
      )
    }
  }

  ps_vec
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
handle_data_frame_weight_calculation <- function(
  weight_fn_numeric,
  .propensity,
  .exposure,
  exposure_type,
  valid_exposure_types = c("auto", "binary", "categorical", "continuous"),
  .propensity_col_quo,
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
    valid_exposure_types
  )

  if (exposure_type_check == "categorical") {
    # For categorical exposures, pass the whole data frame
    return(weight_fn_numeric(
      .propensity = .propensity,
      .exposure = .exposure,
      exposure_type = exposure_type,
      ...
    ))
  }

  # For non-categorical exposures, extract single column
  ps_vec <- extract_propensity_from_df(
    .propensity,
    .propensity_col_quo,
    call = rlang::caller_env(2)
  )

  # Call the numeric method
  weight_fn_numeric(
    .propensity = ps_vec,
    .exposure = .exposure,
    exposure_type = exposure_type,
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

# Helper function to handle optional exposure in GLM methods
extract_exposure_from_glm <- function(
  glm_obj,
  .exposure = NULL,
  call = rlang::caller_env()
) {
  if (is.null(.exposure)) {
    # Extract exposure from GLM
    .exposure <- fmla_extract_left_vctr(glm_obj)
    exposure_name <- fmla_extract_left_chr(glm_obj)
    alert_info("Using exposure variable {.val {exposure_name}} from GLM model")
  }
  .exposure
}
