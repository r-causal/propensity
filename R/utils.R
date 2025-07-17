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
