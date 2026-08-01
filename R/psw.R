#' Propensity Score Weight Vectors
#'
#' @description
#' `psw` objects are numeric vectors that carry metadata about propensity score
#' weights, including the target estimand and whether the underlying propensity
#' scores were trimmed, truncated, or calibrated.
#'
#' Most users will encounter `psw` objects as return values from [wt_ate()] and
#' related weight functions. These constructor and helper functions are useful
#' for inspecting weight objects or for package developers extending propensity.
#'
#' @details
#'
#' ## Constructors
#'
#' * `psw()` is the **user-facing** constructor. It coerces `x` to double and
#'   validates inputs before creating the object.
#' * `new_psw()` is the **low-level** constructor intended for developers. It
#'   assumes `x` is already a double vector and performs minimal validation.
#' * `as_psw()` coerces an existing numeric vector to a `psw` object.
#'
#' ## Queries
#'
#' * `is_psw()` tests whether an object is a `psw` vector.
#' * `is_causal_wt()` tests whether an object inherits from the broader
#'   `causal_wts` class (which includes `psw` objects).
#' * `estimand()` and `estimand<-` get and set the estimand attribute.
#' * `is_stabilized()` returns `TRUE` if the weights are stabilized.
#' * `stabilization_score()` returns the user-supplied stabilization score, or
#'   `NULL` when none was recorded or when a per-observation score was dropped
#'   because an operation changed the length of the weights.
#'
#' ## Arithmetic and combining
#'
#' Arithmetic operations on `psw` objects preserve the class and attributes,
#' so operations like normalization (`weights / sum(weights)`) retain metadata.
#' Combining `psw` objects with [c()] preserves the class only when all
#' metadata matches; mismatched metadata produces a warning and falls back to a
#' plain numeric vector.
#'
#' Subsetting with `[` preserves class and attributes for vector subscripts.
#' Two kinds of attribute hold one value per observation and so cannot be
#' re-indexed for a subset: a `stabilization_score` with more than one value,
#' and the records left by a modified propensity score (`ps_trim_meta`,
#' `ps_trunc_meta`, and `ps_calib_meta`). Where an operation goes through
#' vctrs, these are carried when the result comes back at the length they were
#' recorded on and dropped when it does not. Any same-length operation keeps
#' them, a reordering or a subscript with duplicates included, so the positions
#' a record names can end up describing different observations than they did.
#'
#' Dropping the `stabilization_score` warns, because the score was supplied by
#' the user and the weights can be recomputed on the subset. Dropping a
#' modification record is silent, because these records also travel by routes
#' vctrs does not see, and a warning on the one route it controls would be
#' neither complete nor about anything the user wrote. Subassignment that grows
#' the vector carries a record across the length change, since `[<-` casts the
#' replacement and then leaves base R to preserve the attributes; and
#' `model.frame()` drops the `NA`-weighted rows from a weights column in C and
#' re-attaches the original variable's attributes to the shortened result, so
#' weights built on trimmed propensity scores come back out of every outcome
#' model fit on them still carrying a record written for rows that are no
#' longer there.
#'
#' Honesty therefore lives at query time. [is_unit_trimmed()] answers by
#' position, so it checks that the record covers the vector it is given and
#' raises an error of class `propensity_missing_meta_error` when it does not,
#' or when weights marked as trimmed carry no record at all, rather than name
#' trimmed units at stale positions. [is_refit()] reads a single flag rather
#' than a position, so it answers from any record present and refuses only when
#' the record is absent entirely.
#'
#' The result of any of these operations stays a `psw` and keeps every other
#' attribute, including its stabilized, trimmed, truncated, and calibrated
#' status, and the attributes describing a categorical exposure, which name the
#' exposure levels rather than the units and so mean the same thing at any
#' length.
#'
#' Matrix or array subscripts intentionally drop the `psw` class and return
#' a plain numeric vector via base R linear indexing; this is required so
#' `glm.fit()`-style internal indexing works on `psw`-weighted GLMs.
#' Summary functions ([sum()], [mean()], etc.) return plain numeric values.
#'
#' @import vctrs
#' @export
#'
#' @param x For `psw()` and `new_psw()`: a numeric vector of weights
#'   (default: `double()`). For `is_psw()`, `is_causal_wt()`, and `as_psw()`:
#'   an object to test or coerce.
#' @param estimand A character string identifying the target estimand (e.g.,
#'   `"ate"`, `"att"`, `"ato"`). Defaults to `NULL`.
#' @param stabilized Logical. Were the weights stabilized? Defaults to `FALSE`.
#' @param trimmed Logical. Were the weights derived from trimmed propensity
#'   scores? Defaults to `FALSE`.
#' @param truncated Logical. Were the weights derived from truncated propensity
#'   scores? Defaults to `FALSE`.
#' @param calibrated Logical. Were the weights derived from calibrated
#'   propensity scores? Defaults to `FALSE`.
#' @param stabilization_score Optional numeric stabilization score to record on
#'   the object, either a single value or one value per observation. Defaults to
#'   `NULL`, meaning no fixed score was supplied.
#' @param wt A `psw` or `causal_wts` object.
#' @param ... Additional attributes stored on the object (developer use only).
#'
#' @return
#' * `new_psw()`, `psw()`, `as_psw()`: A `psw` vector.
#' * `is_psw()`, `is_causal_wt()`, `is_stabilized()`: A single logical value.
#' * `estimand()`: A character string, or `NULL` if no estimand is set.
#' * `estimand<-`: The modified `psw` object (called for its side effect).
#' * `stabilization_score()`: A numeric value or vector, or `NULL` if none was
#'   recorded or a per-observation score was dropped.
#'
#' @seealso
#' [wt_ate()], [wt_att()], [wt_atu()], [wt_atm()], [wt_ato()] for
#' calculating propensity score weights (which return `psw` objects).
#'
#' [ps_trim()], [ps_trunc()], and [ps_calibrate()] for modifying propensity
#' scores before weight calculation.
#'
#' @examples
#' # Create psw objects directly
#' w <- psw(c(1.2, 0.8, 1.5), estimand = "ate")
#' w
#'
#' # Query metadata
#' is_psw(w)
#' estimand(w)
#' is_stabilized(w)
#'
#' # Coerce a plain numeric vector
#' as_psw(c(1.0, 2.0), estimand = "att")
#'
#' # Arithmetic preserves the psw class
#' w / sum(w)
#'
#' # Combining: compatible metadata is preserved
#' x <- psw(c(1.2, 0.8), estimand = "ate")
#' y <- psw(c(1.1, 0.9), estimand = "ate")
#' c(x, y)
#'
#' # Combining: incompatible metadata warns and returns numeric
#' x <- psw(c(1.2, 0.8), estimand = "ate")
#' y <- psw(c(1.1, 0.9), estimand = "att")
#' c(x, y)
#'
#' @import vctrs
#' @name psw
NULL

#' @rdname psw
#' @export
new_psw <- function(
  x = double(),
  estimand = NULL,
  stabilized = FALSE,
  trimmed = FALSE,
  truncated = FALSE,
  calibrated = FALSE,
  stabilization_score = NULL,
  ...
) {
  vec_assert(x, ptype = double())
  vec_assert(stabilized, ptype = logical(), size = 1)

  new_vctr(
    x,
    estimand = estimand,
    stabilized = stabilized,
    trimmed = trimmed,
    truncated = truncated,
    calibrated = calibrated,
    stabilization_score = stabilization_score,
    ...,
    class = c("psw", "causal_wts"),
    inherit_base_type = TRUE
  )
}


#' @rdname psw
#' @export
psw <- function(
  x = double(),
  estimand = NULL,
  stabilized = FALSE,
  trimmed = FALSE,
  truncated = FALSE,
  calibrated = FALSE,
  stabilization_score = NULL
) {
  x <- vec_cast(x, to = double())
  attributes(x) <- NULL
  new_psw(
    x,
    estimand = estimand,
    stabilized = stabilized,
    trimmed = trimmed,
    truncated = truncated,
    calibrated = calibrated,
    stabilization_score = stabilization_score
  )
}

#' @rdname psw
#' @export
is_psw <- function(x) {
  inherits(x, "psw")
}

#' @rdname psw
#' @export
is_stabilized <- function(wt) {
  isTRUE(attr(wt, "stabilized"))
}

#' @rdname psw
#' @export
stabilization_score <- function(wt) {
  attr(wt, "stabilization_score")
}

#' @rdname psw
#' @export
as_psw <- function(x, estimand = NULL) {
  x <- vec_cast(x, to = double())
  psw(x, estimand = estimand)
}

#' @export
is_ps_trimmed.psw <- function(x) {
  isTRUE(attr(x, "trimmed"))
}

# The number of observations a trimming record was written for. `keep_idx` and
# `trimmed_idx` partition those observations in every `ps_trim()` method, so
# their total is the length the record is entitled to speak for.
recorded_trim_length <- function(meta) {
  length(meta$keep_idx) + length(meta$trimmed_idx)
}

# `recorded` is the length the record covers, or NULL when there is no record.
abort_trim_meta <- function(fn, recorded, n, call) {
  problem <- if (is.null(recorded)) {
    "These weights are marked as built from trimmed propensity scores but
     carry no record of which units were trimmed."
  } else {
    "The record covers {recorded} observation{?s} and these weights have {n},
     so its positions do not describe them."
  }

  abort(
    c(
      "{.code {fn}} has no usable trimming record for these weights.",
      x = problem,
      i = "Call {.code {fn}} on the {.cls ps_trim} object the weights were
           built from, or rebuild the weights from it."
    ),
    error_class = "propensity_missing_meta_error",
    call = call
  )
}

# Both unit-level queries read their answer out of `ps_trim_meta`, so weights
# marked as trimmed with no record left have no answer to give: reporting every
# unit as retained, or no refit as having happened, would be a wrong answer
# rather than a missing one.
check_trim_meta <- function(x, fn, call = rlang::caller_env()) {
  if (!is_ps_trimmed(x) || !is.null(ps_trim_meta(x))) {
    return(invisible(x))
  }

  abort_trim_meta(fn, recorded = NULL, n = length(x), call = call)
}

# A positional query needs more than a record: it needs one written for the
# vector in front of it. A record can outlive the observations it describes,
# because it does not only travel by routes vctrs controls. `model.frame()`
# drops the `NA`-weighted rows from a weights column in C and re-attaches the
# original variable's attributes to the shortened result, so weights built on
# trimmed propensity scores come back out of every outcome model fit on them
# still carrying a record written for rows that are no longer there. Growing a
# psw by subassignment crosses a length change the same way. Answering from
# either would name trimmed units at positions holding retained ones.
check_trim_meta_covers <- function(x, fn, call = rlang::caller_env()) {
  check_trim_meta(x, fn, call = call)

  meta <- ps_trim_meta(x)
  if (is.null(meta) || recorded_trim_length(meta) == length(x)) {
    return(invisible(x))
  }

  abort_trim_meta(
    fn,
    recorded = recorded_trim_length(meta),
    n = length(x),
    call = call
  )
}

#' @export
is_unit_trimmed.psw <- function(x) {
  # No observations, no answers. A record kept on an empty vector describes
  # observations it does not have, and indexing an empty logical by the
  # positions it names would grow one padded with `NA`.
  if (length(x) == 0) {
    return(logical(0))
  }

  out <- vector("logical", length = length(x))
  if (!is_ps_trimmed(x)) {
    return(out)
  }

  check_trim_meta_covers(x, "is_unit_trimmed()")

  meta <- ps_trim_meta(x)
  out[meta$trimmed_idx] <- TRUE

  out
}

#' @export
is_ps_truncated.psw <- function(x) {
  isTRUE(attr(x, "truncated"))
}

#' @export
is_unit_truncated.psw <- function(x) {
  isTRUE(attr(x, "truncated"))
}


#' @export
is_refit.psw <- function(x) {
  # `refit` is one flag about the propensity model, not a fact about any
  # position, so a record that no longer covers these weights still answers it.
  check_trim_meta(x, "is_refit()")

  isTRUE(ps_trim_meta(x)$refit)
}

#' @export
vec_ptype_abbr.psw <- function(x, ...) {
  estimand <- estimand(x)
  if (is.null(estimand)) {
    "psw"
  } else {
    paste0("psw{", estimand, "}")
  }
}

#' @export
vec_ptype_full.psw <- function(x, ...) {
  estimand <- estimand(x)
  if (is_stabilized(x)) {
    stabilized <- "; stabilized"
  } else {
    stabilized <- NULL
  }

  if (is.null(estimand)) {
    paste0("psw{estimand = unknown", stabilized, "}")
  } else {
    paste0("psw{estimand = ", estimand, stabilized, "}")
  }
}

#' @export
#' @method vec_arith psw
vec_arith.psw <- function(op, x, y, ...) {
  UseMethod("vec_arith.psw", y)
}

#' @export
#' @method vec_arith.psw default
vec_arith.psw.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.psw psw
vec_arith.psw.psw <- function(op, x, y, ...) {
  estimand_x <- estimand(x)
  estimand_y <- estimand(y)
  if (!identical(estimand_x, estimand_y)) {
    estimand <- paste0(estimand_x, ", ", estimand_y)
  } else {
    estimand <- estimand_x
  }

  # Determine stabilized status: both must be stabilized for result to be stabilized
  stabilized <- is_stabilized(x) && is_stabilized(y)

  # For trimmed/truncated/calibrated: if either has the property, result has it
  trimmed <- is_ps_trimmed(x) || is_ps_trimmed(y)
  truncated <- is_ps_truncated(x) || is_ps_truncated(y)
  calibrated <- is_ps_calibrated(x) || is_ps_calibrated(y)

  rslts <- vec_arith_base(op, x, y)
  psw(
    rslts,
    estimand = estimand,
    stabilized = stabilized,
    trimmed = trimmed,
    truncated = truncated,
    calibrated = calibrated
  )
}

#' @export
#' @method vec_arith.psw MISSING
vec_arith.psw.MISSING <- function(op, x, y, ...) {
  switch(
    op,
    `-` = vec_restore(-1 * vec_data(x), x), # Returns psw (preserves class)
    `+` = x, # Returns psw unchanged
    stop_incompatible_op(op, x, y)
  )
}

#' @export
#' @method vec_arith.psw numeric
vec_arith.psw.numeric <- function(op, x, y, ...) {
  result <- vec_arith_base(op, x, y)
  vec_restore(result, x)
}

#' @export
#' @method vec_arith.numeric psw
vec_arith.numeric.psw <- function(op, x, y, ...) {
  result <- vec_arith_base(op, x, y)
  vec_restore(result, y)
}

#' @export
#' @method vec_arith.psw integer
vec_arith.psw.integer <- function(op, x, y, ...) {
  result <- vec_arith_base(op, x, y)
  vec_restore(result, x)
}

# A stabilization score of length > 1 is indexed by observation, so it is only
# meaningful at the length it was recorded on. Nothing rebuilding a psw is given
# the indices behind a length change, and the subscript method for this class is
# owned upstream, so per-observation metadata cannot be re-indexed here.
#
# Zero-length data is exempt. A prototype or empty subset carries no
# observations, so a score on it lines up with nothing and misleads no one, and
# keeping it lets the restore that builds the real result compare the score
# against the length that observations actually arrive at.
stabilization_score_aligns <- function(score, n) {
  n == 0 || length(score) <= 1 || length(score) == n
}

# The records a modified propensity score leaves on the weights built from it.
# Each holds indices into the observations it was recorded on, so like a
# per-observation stabilization score, each is only meaningful at that length.
psw_modification_meta <- c("ps_trim_meta", "ps_trunc_meta", "ps_calib_meta")

# Attributes describing a categorical exposure. These name the exposure levels
# rather than the units, so they carry no length of their own.
psw_categorical_attrs <- c("n_categories", "category_names", "focal_category")

# A modification record is judged against `to`, the object supplying the type:
# the record was written for `to`'s observations, so data arriving at `to`'s
# length is the only data those indices describe. Zero-length data is exempt for
# the reason a stabilization score is: it lines up with nothing and so
# contradicts nothing, and a prototype that keeps the record lets the restore
# building the real result carry it on to the observations.
#
# Any same-length operation keeps indices that may no longer point where they
# did, a reordering or a subscript with duplicates included. Nothing rebuilding
# a psw is handed the subscript, so neither can be told apart from an untouched
# vector here. `is_unit_trimmed()` checks the record covers the vector it is
# given, which catches a length change from any route, but a same-length
# rearrangement is beyond what either can see and is documented rather than
# guarded.
modification_meta_aligns <- function(x, to) {
  length(x) == 0 || length(x) == length(to)
}

# The psw type is every attribute the object carries, so anything that rebuilds
# a psw around another vector's data owes all of them to the object supplying
# the type. The attributes indexed by observation are the exception, and are
# dropped when the data arrives at a length they do not describe: a later read
# sees an absent record rather than one silently misaligned with the weights.
carry_psw_metadata <- function(x, to) {
  score <- stabilization_score(to)
  if (!stabilization_score_aligns(score, length(x))) {
    score <- NULL
  }

  out <- new_psw(
    x,
    estimand = estimand(to),
    stabilized = is_stabilized(to),
    trimmed = is_ps_trimmed(to),
    truncated = is_ps_truncated(to),
    calibrated = is_ps_calibrated(to),
    stabilization_score = score
  )

  if (modification_meta_aligns(x, to)) {
    for (meta in psw_modification_meta) {
      attr(out, meta) <- attr(to, meta)
    }
  }

  for (attribute in psw_categorical_attrs) {
    attr(out, attribute) <- attr(to, attribute)
  }

  out
}

#' @export
vec_restore.psw <- function(x, to, ...) {
  # Extract numeric data if needed
  if (inherits(x, "psw")) {
    x <- vec_data(x)
  }

  # Only a restore announces the drop, because only a restore has an audience
  # for it: `to` is the weights the operation was performed on, so a user who
  # recorded a per-observation score has something to act on. In a cast, `to`
  # supplies the type and its score describes `to`'s own observations rather
  # than the incoming data. Outside vctrs' subassignment intermediates, where
  # the restore that follows reattaches the target's score anyway, a cast's `to`
  # is a prototype, which cannot carry a score misaligned with itself.
  score <- stabilization_score(to)
  if (!stabilization_score_aligns(score, length(x))) {
    warn(
      c(
        "Dropping the per-observation {.arg stabilization_score}.",
        i = "A score of length {length(score)} cannot be carried through an
             operation that changes the length of the weights to
             {length(x)}.",
        i = "The result is still marked as stabilized. Recompute the weights
             on the subset if you need a score aligned with them."
      ),
      warning_class = "propensity_stabilization_score_warning",
      # Restore is reached through vctrs' internal dispatch, whose call would
      # be reported here and names nothing the caller wrote.
      call = NULL
    )
  }

  carry_psw_metadata(x, to)
}

#' @export
vec_ptype2.psw.psw <- function(x, y, ...) {
  # Check estimand compatibility
  if (!identical(estimand(x), estimand(y))) {
    warn_incompatible_metadata(
      x,
      y,
      paste0(
        "incompatible estimands '",
        estimand(x),
        "' and '",
        estimand(y),
        "'"
      )
    )
    return(double())
  }

  # Check stabilization status
  if (!identical(attr(x, "stabilized"), attr(y, "stabilized"))) {
    warn_incompatible_metadata(x, y, "different stabilization status")
    return(double())
  }

  # Check stabilization score
  if (!identical(stabilization_score(x), stabilization_score(y))) {
    warn_incompatible_metadata(x, y, "different stabilization scores")
    return(double())
  }

  # Check trimmed status
  if (!identical(is_ps_trimmed(x), is_ps_trimmed(y))) {
    warn_incompatible_metadata(x, y, "different trimming status")
    return(double())
  }

  # Check truncated status
  if (!identical(is_ps_truncated(x), is_ps_truncated(y))) {
    warn_incompatible_metadata(x, y, "different truncation status")
    return(double())
  }

  # Check calibrated status
  if (!identical(is_ps_calibrated(x), is_ps_calibrated(y))) {
    warn_incompatible_metadata(x, y, "different calibration status")
    return(double())
  }

  # If all metadata matches, return psw with all attributes preserved
  new_psw(
    estimand = estimand(x),
    stabilized = is_stabilized(x),
    trimmed = is_ps_trimmed(x),
    truncated = is_ps_truncated(x),
    calibrated = is_ps_calibrated(x),
    stabilization_score = stabilization_score(x)
  )
}

#' @export
vec_ptype2.psw.double <- function(x, y, ...) {
  warn_class_downgrade("psw")
  double()
}

#' @export
vec_ptype2.double.psw <- function(x, y, ...) {
  warn_class_downgrade("psw")
  double()
}

# A cast returns `x`'s data in `to`'s type, so the data is stripped back to bare
# doubles and the whole type, metadata included, comes from `to`.
cast_to_psw <- function(x, to) {
  x <- vec_cast(vec_data(x), to = double())
  attributes(x) <- NULL

  carry_psw_metadata(x, to)
}

#' @export
vec_cast.psw.psw <- function(x, to, ...) x

#' @export
vec_cast.psw.double <- function(x, to, ...) cast_to_psw(x, to)

#' @export
vec_cast.double.psw <- function(x, to, ...) vec_data(x)

#' @export
vec_ptype2.psw.character <- function(x, y, ...) {
  warn_class_downgrade("psw", "character")
  character()
}

#' @export
vec_ptype2.character.psw <- function(x, y, ...) {
  warn_class_downgrade("psw", "character")
  character()
}

#' @export
vec_cast.character.psw <- function(x, to, ...) as.character(vec_data(x))


#' @export
vec_ptype2.psw.integer <- function(x, y, ...) {
  warn_class_downgrade("psw", "integer")
  integer()
}

#' @export
vec_ptype2.integer.psw <- function(x, y, ...) {
  warn_class_downgrade("psw", "integer")
  integer()
}

#' @export
vec_cast.psw.integer <- function(x, to, ...) cast_to_psw(x, to)

#' @export
vec_cast.integer.psw <- function(x, to, ...) {
  vec_cast(vec_data(x), integer(), x_arg = "psw")
}

#' @export
vec_cast.psw.ps_trim <- function(x, to, ...) cast_to_psw(x, to)

#' @export
vec_cast.ps_trim.psw <- function(x, to, ...) {
  ps_trim(vec_data(x), method = "ps", lower = 0, upper = 1)
}

#' @export
vec_cast.psw.ps_trunc <- function(x, to, ...) cast_to_psw(x, to)

#' @export
vec_cast.ps_trunc.psw <- function(x, to, ...) {
  ps_trunc(vec_data(x), method = "ps", lower = 0, upper = 1)
}
