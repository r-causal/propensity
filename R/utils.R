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

# `unique()` takes values to hold out of the comparison. These methods find
# duplicates with `vctrs::vec_unique_loc()`, which has no equivalent, so an
# `incomparables` argument would be compared like any other value and the
# caller would get an answer to a question they did not ask.
check_incomparables <- function(
  incomparables,
  class_name,
  call = rlang::caller_env()
) {
  if (isFALSE(incomparables)) {
    return(invisible(incomparables))
  }

  abort(
    c(
      "{.arg incomparables} is not supported for {.cls {class_name}} objects.",
      x = "Values named there would be compared along with the rest, so the
           result would not be the one asked for.",
      i = "Call {.fun unique} on {.code as.numeric(x)} to hold values out of
           the comparison. The result is a plain numeric vector and carries no
           record of which units were modified."
    ),
    error_class = "propensity_unsupported_arg_error",
    call = call
  )
}

# Coercion warning helpers
warn_incompatible_metadata <- function(
  x,
  y,
  reason,
  ...,
  call = rlang::caller_env()
) {
  x_class <- class(x)[1]
  y_class <- class(y)[1]

  warn(
    c(
      paste0("Converting ", x_class, " to numeric: ", reason),
      i = "Metadata cannot be preserved when combining incompatible objects",
      i = "Use identical objects or explicitly cast to numeric to avoid this warning"
    ),
    warning_class = "propensity_coercion_warning",
    call = call
  )
}

warn_class_downgrade <- function(
  from_classes,
  to_class = "numeric",
  ...,
  call = rlang::caller_env()
) {
  # Handle both single class and multiple classes
  if (length(from_classes) == 1) {
    from_text <- from_classes
  } else {
    from_text <- paste(from_classes, collapse = " and ")
  }

  warn(
    c(
      paste0("Converting ", from_text, " to ", to_class),
      i = "Class-specific attributes and metadata have been dropped",
      i = "Use explicit casting to numeric to avoid this warning"
    ),
    warning_class = "propensity_class_downgrade_warning",
    call = call
  )
}

# Modification record helpers

# `ps_trim()` and `ps_trunc()` each leave a record of which units the
# modification touched. A record names positions, and a position only means
# something against a known number of observations, so the record writes that
# number down rather than leaving it to be worked out: `truncated_idx` names
# only the units that changed, and a record whose positions have been dropped
# names none at all.
record_covers <- function(meta, n) {
  !is.null(meta$n_obs) && meta$n_obs == n
}

# `[` knows the subscript, which is what re-indexing a record takes. Every
# occurrence of a recorded position is mapped onto the position it now holds, so
# a subscript naming a position twice reports that unit twice. `NA` names no
# position, so an element taken by one falls in neither set.
reindex_positions <- function(positions, i) {
  which(i %in% positions)
}

# The rows `x[i, ]` is built from, as positions in `x`, which is the form a
# record is re-indexed against. A row subscript is not already in that form:
# character subscripts name rows rather than positions, a logical subscript is
# recycled over the rows and takes a whole row per element rather than per
# `TRUE`, and a zero selects nothing. `NA` stands for a row that names no
# observation in `x`, which is what base returns for it, and is left as `NA`
# rather than dropped so the positions stay aligned with the rows.
subscript_row_positions <- function(i, x) {
  n <- nrow(x)

  if (is.character(i)) {
    return(match(i, rownames(x)))
  }

  if (is.logical(i)) {
    i <- rep_len(i, n)
    taken <- i | is.na(i)
    return(replace(seq_len(n), is.na(i), NA_integer_)[taken])
  }

  i <- as.integer(i)

  # Negative and missing subscripts never mix, so a subscript holding an `NA`
  # names positions to take rather than positions to leave out.
  if (!anyNA(i) && any(i < 0)) {
    return(setdiff(seq_len(n), -i))
  }

  i[is.na(i) | i != 0]
}
