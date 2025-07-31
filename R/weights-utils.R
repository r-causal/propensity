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
