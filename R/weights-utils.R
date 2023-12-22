
abort_unsupported <- function(exposure_type, what) {
  rlang::abort("Exposure type {exposure_type} not currently supported for {what}")
}

match_exposure_type <- function(exposure_type = c("auto", "binary", "categorical", "continuous"), .exposure) {
  .exposure_type <- rlang::arg_match(exposure_type)
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

transform_exposure_binary <- function(.exposure, .treated = NULL, .untreated = NULL) {
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
    levels <- if (is.factor(.exposure)) levels(.exposure) else sort(unique(.exposure))
    alert_info("Setting treatment to {.var {levels[[2]]}}")
    return(ifelse(.exposure == levels[[2]], 1, 0))
  } else {
    rlang::abort(c(
      "Don't know how to transform `.exposure` to 0/1 binary variable.",
      i = "Specify `.treated` and `.untreated.`"
    ))
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
