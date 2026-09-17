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

inform <- function(
  ...,
  message_class = NULL,
  call = rlang::caller_env(),
  .envir = parent.frame()
) {
  cli::cli_inform(
    ...,
    class = c(message_class, "propensity_message"),
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

# Predictions from a binary model arrive as one column per level, and the
# binary path works on a single column of scores. The second column of a pair is
# the probability of the second level in the layout those predictions come in,
# and the only column otherwise. The caller did not make that choice, so it is
# announced rather than left to be inferred from the result.
binary_ps_column <- function(ps, fn) {
  column <- if (ncol(ps) == 2) 2L else 1L
  name <- names(ps)[[column]]

  alert_info(
    "Using the {.val {name}} column as the propensity score for {.fun {fn}}."
  )

  ps[[column]]
}

# The names a frame of predictions offers for matching against exposure levels.
# Predicting probabilities from a fitted parsnip model returns one column per
# level, named `.pred_<level>`, so the level a column holds is the part of the
# name after that prefix. A frame whose columns are not all named that way is
# read under the names it has: stripping a prefix from some of them would match
# levels against a mixture of two conventions.
level_column_names <- function(nms) {
  if (length(nms) == 0 || anyNA(nms) || !all(grepl("^\\.pred_", nms))) {
    return(nms)
  }

  sub("^\\.pred_", "", nms)
}

# Predicting from a fitted classification model without naming a type returns a
# single `.pred_class` column, which holds the level predicted for each unit
# rather than a probability. It arrives in the same shape as a frame of scores
# and is no propensity score at all, so every route that reads a frame refuses
# it here. Left to run on, it is reported downstream as levels that are not
# numeric or not between zero and one, neither of which tells the caller that
# the prediction they need is a different one.
check_predicted_class_column <- function(
  .propensity,
  call = rlang::caller_env()
) {
  if (!is.data.frame(.propensity) || !".pred_class" %in% names(.propensity)) {
    return(invisible(FALSE))
  }

  abort(
    c(
      "{.arg .propensity} holds predicted classes rather than propensity
       scores.",
      x = "The {.val .pred_class} column holds the level predicted for each
           unit.",
      i = "Predict probabilities instead, with
           {.code predict(fit, new_data, type = \"prob\")}, and pass the
           columns that returns."
    ),
    error_class = "propensity_df_class_column_error",
    call = call
  )
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

# The scores a modified propensity score holds, one unit per element or, for a
# matrix, one unit per row. vctrs compares a bare matrix row by row, which is
# what `unique()` of a matrix of scores merges.
score_values <- function(x) {
  if (is.matrix(x)) {
    out <- unclass(x)[,, drop = FALSE]
    dimnames(out) <- NULL
    out
  } else {
    vctrs::vec_data(x)
  }
}

# The units at `loc`, taken through the class's own subsetting so the result
# keeps its class and a record re-indexed onto them. A matrix is subset by row;
# a single subscript would index it as one long vector.
subset_score_units <- function(x, loc) {
  if (is.matrix(x)) {
    x[loc, , drop = FALSE]
  } else {
    x[loc]
  }
}

# Whether every set of units `unique()` merges into one element shares a single
# status. `vec_group_id()` groups values the way `vec_unique_loc()` keeps them,
# missing values included, so each group is one element of the result.
merged_units_agree <- function(values, status) {
  groups <- vctrs::vec_group_id(values)
  pairs <- data.frame(group = as.integer(groups), status = status)

  vctrs::vec_unique_count(pairs) == attr(groups, "n")
}

# A modified score matrix as a data frame. A row of the matrix is a unit, and so
# is a row of the frame, so every column is a score vector of the matrix's class
# carrying the matrix's record unchanged: its positions name rows, which are the
# elements of each column. The frame's own rows and names are the ones base R
# gives the plain matrix.
modified_score_matrix_frame <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  record_attr,
  build_column,
  ...
) {
  meta <- attr(x, record_attr)
  values <- unclass(x)
  attr(values, record_attr) <- NULL
  out <- as.data.frame(values, row.names = row.names, optional = optional, ...)

  for (j in seq_along(out)) {
    out[[j]] <- build_column(unname(as.double(out[[j]])), meta)
  }

  out
}

# The score matrix a data frame of modified score columns stands for, or `NULL`
# when no column is modified. The frame describes one matrix only when every
# column is the same class and carries the same record, as `as.data.frame()`
# leaves them and as row subsetting keeps them; anything else mixes scores that
# were modified differently, or not at all, into one set of units.
frame_as_modified_score_matrix <- function(frame, call = rlang::caller_env()) {
  columns <- unclass(frame)
  kinds <- vapply(
    columns,
    function(column) {
      if (inherits(column, "ps_trim")) {
        "trim"
      } else if (inherits(column, "ps_trunc")) {
        "trunc"
      } else {
        "plain"
      }
    },
    character(1)
  )

  if (all(kinds == "plain")) {
    return(NULL)
  }

  kind <- kinds[[1]]
  record_attr <- if (identical(kind, "trim")) {
    "ps_trim_meta"
  } else {
    "ps_trunc_meta"
  }
  records <- lapply(columns, attr, record_attr)
  agree <- all(kinds == kind) &&
    all(vapply(records, identical, logical(1), records[[1]]))

  if (!agree) {
    abort(
      c(
        "The columns of {.arg .propensity} must be modified alike.",
        x = "Some columns are trimmed or truncated propensity scores and the
             others are not, or their records differ, so the columns do not
             describe one set of units.",
        i = "Build the data frame from a single trimmed or truncated score
             matrix with {.fn as.data.frame}, or pass the matrix itself."
      ),
      error_class = "propensity_matrix_type_error",
      call = call
    )
  }

  values <- do.call(cbind, lapply(columns, vctrs::vec_data))
  row_names <- if (.row_names_info(frame) > 0) rownames(frame)
  dimnames(values) <- list(row_names, names(frame))

  if (identical(kind, "trim")) {
    list(
      scores = new_trimmed_ps(values, ps_trim_meta = records[[1]]),
      modification_type = "trim"
    )
  } else {
    list(
      scores = new_ps_trunc(values, records[[1]]),
      modification_type = "trunc"
    )
  }
}

# The units a modified score holds, one per element or, for a matrix, one per
# row, and those of them with a missing score. A row with any missing score has
# no complete probability vector, so the whole unit is missing.
missing_unit_count <- function(x) {
  if (is.matrix(x)) nrow(x) else length(x)
}

missing_units <- function(x) {
  values <- score_values(x)
  if (is.matrix(values)) {
    which(rowSums(is.na(values)) > 0)
  } else {
    which(is.na(values))
  }
}

# `[` knows the subscript, which is what re-indexing a record takes. Every
# occurrence of a recorded position is mapped onto the position it now holds, so
# a subscript naming a position twice reports that unit twice. `NA` names no
# position, so an element taken by one falls in neither set, and neither does
# one taken from past the last recorded position.
#
# `i` holds positions as `vec_as_location()` returns them, positive or `NA`. A
# zero or negative value would drop or exclude elements of the mask rather than
# name a position.
reindex_positions <- function(positions, i) {
  mask <- logical(max(c(positions, 0L)))
  mask[positions] <- TRUE
  which(mask[i])
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

# One run of spaces where a string carries several. Text assembled from
# `deparse()` output carries the indentation of the continuation lines, which
# reads as a gap rather than as a space once the lines are joined.
squash_whitespace <- function(x) {
  gsub("[[:space:]]+", " ", x)
}
