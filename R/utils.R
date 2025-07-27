be_quiet <- function() {
  getOption("propensity.quiet", default = FALSE)
}

abort <- function(
  ...,
  error_class = NULL,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  cli::cli_abort(
    ...,
    class = c(error_class, "propensity_error"),
    call = call,
    .envir = .envir
  )
}

warn <- function(
  ...,
  warning_class = NULL,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  cli::cli_warn(
    ...,
    class = c(warning_class, "propensity_warning"),
    call = call,
    .envir = .envir
  )
}

alert_info <- function(.message, .envir = parent.frame()) {
  if (!be_quiet()) {
    cli::cli_alert_info(text = .message, .envir = .envir)
  }
}

assert_class <- function(
  x,
  classes,
  .length = NULL,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  classes <- as.character(classes)
  .stop <- FALSE
  .msg <- if (length(classes) == 1) {
    "{.arg {arg}} must be of class {.val {classes}}."
  } else {
    "{.arg {arg}} must be one of class {.val {classes}}."
  }
  .class_msg <- NULL
  .length_msg <- NULL

  if (!any(vapply(classes, function(cls) inherits(x, cls), logical(1)))) {
    .stop <- TRUE
    .class_msg <- "It has class {.val {class(x)}}."
  }

  if (!is.null(.length) && length(x) != .length) {
    .stop <- TRUE
    .msg <- if (length(classes) == 1) {
      "{.arg {arg}} must be of class {.val {classes}} and length {.val { .length}}."
    } else {
      "{.arg {arg}} must be one of class {.val {classes}} and length {.val { .length}}."
    }
    .length_msg <- "It has length {.val {length(x)}}."
  }

  if (.stop) {
    abort(
      c(
        .msg,
        x = .class_msg,
        x = .length_msg
      ),
      error_class = "propensity_class_error",
      call = call
    )
  }

  invisible(TRUE)
}

assert_columns_exist <- function(
  .df,
  names_vec,
  arg = rlang::caller_arg(.df),
  call = rlang::caller_env()
) {
  missing <- setdiff(names_vec, names(.df))
  if (length(missing) > 0) {
    abort(
      "The data frame {.arg {arg}} is missing the {.val {missing}} column{?s}.",
      error_class = "propensity_columns_exist_error",
      call = call
    )
  }

  invisible(TRUE)
}

#' Cut Numeric Vector into Quantile-Based Intervals
#'
#' `cut_quantile()` is a wrapper around [cut()] that uses quantiles as 
#' breakpoints. By default, it creates 10 intervals based on deciles 
#' (quantiles at 0, 0.1, 0.2, ..., 1.0).
#'
#' @param x A numeric vector to be cut into intervals.
#' @param probs A numeric vector of probabilities with values in [0,1] 
#'   giving the quantiles to use as breakpoints. Default creates deciles.
#' @param include.lowest Logical, indicating if an 'x[i]' equal to the 
#'   lowest (or highest, for right = FALSE) 'breaks' value should be 
#'   included. Default is `TRUE`.
#' @param ... Additional arguments passed to [cut()].
#'
#' @return A factor with levels corresponding to the quantile-based intervals.
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' 
#' # Default: cut into deciles
#' cut_quantile(x)
#' 
#' # Custom quantiles: cut into quartiles
#' cut_quantile(x, probs = seq(0, 1, 0.25))
#' 
#' # Cut into tertiles
#' cut_quantile(x, probs = c(0, 1/3, 2/3, 1))
#'
#' @export
cut_quantile <- function(x, probs = seq(0, 1, 0.1), include.lowest = TRUE, ...) {
  if (!is.numeric(x)) {
    abort("{.arg x} must be numeric.")
  }
  
  if (!is.numeric(probs) || any(probs < 0) || any(probs > 1)) {
    abort("{.arg probs} must be numeric values between 0 and 1.")
  }
  
  if (length(probs) < 2) {
    abort("{.arg probs} must have at least 2 values to create intervals.")
  }
  
  # Remove NA values for quantile calculation, but preserve them in result
  x_complete <- x[!is.na(x)]
  
  if (length(x_complete) == 0) {
    abort("All values in {.arg x} are NA.")
  }
  
  # Calculate quantile breaks
  breaks <- stats::quantile(x_complete, probs = probs, na.rm = TRUE)
  
  # Check if breaks are unique (handles case where all values are identical)
  unique_breaks <- unique(breaks)
  if (length(unique_breaks) < 2) {
    # If all values are the same, create a single interval
    warn("All values are identical; creating a single interval.")
    # Create a factor with all values in the same level
    result <- factor(rep(1, length(x)), levels = 1, 
                    labels = paste0("[", min(x_complete), ",", max(x_complete), "]"))
    result[is.na(x)] <- NA
    return(result)
  }
  
  # Use cut() with quantile breaks
  cut(x, breaks = breaks, include.lowest = include.lowest, ...)
}
