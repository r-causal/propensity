abort_unsupported <- function(exposure_type, what, call = rlang::caller_env()) {
  abort(
    "Exposure type {.val {exposure_type}} not currently supported for {.field {what}}",
    call = call,
    error_class = "propensity_wt_not_supported_error"
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
  length(unique(.exposure)) / sum(!is.na(.exposure)) < 0.2
}

has_two_levels <- function(.x) {
  length(unique(.x)) == 2
}

check_refit <- function(.propensity) {
  if (!is_refit(.propensity)) {
    warn(
      c(
        "It appears you trimmed your propensity score but did not refit the model.",
        i = "Use {.code ps_refit()} for more accurate re-estimation."
      ),
      warning_class = "propensity_no_refit_warning"
    )
  }
}

check_ps_range <- function(ps, call = rlang::caller_env()) {
  ps <- as.numeric(ps)
  if (any(ps <= 0 | ps >= 1, na.rm = TRUE)) {
    abort(
      c(
        "The propensity score must be between 0 and 1.",
        i = "The range of {.arg ps} is \\
      {format(range(ps), nsmall = 1, digits = 1)}"
      ),
      call = call,
      error_class = "propensity_range_error"
    )
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
  row_sums <- rowSums(ps_matrix)
  ROW_SUM_TOLERANCE <- 1e-6 # Tolerance for floating point comparison
  if (any(abs(row_sums - 1) > ROW_SUM_TOLERANCE)) {
    bad_rows <- which(abs(row_sums - 1) > ROW_SUM_TOLERANCE)
    abort(
      c(
        "Propensity score matrix rows must sum to 1.",
        i = "Problem rows: {bad_rows[1:min(5, length(bad_rows))]}{if (length(bad_rows) > 5) ' ...' else ''}"
      ),
      call = call,
      error_class = "propensity_matrix_sum_error"
    )
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
      warning_class = "propensity_matrix_no_names_warning"
    )
  }

  ps_matrix
}
