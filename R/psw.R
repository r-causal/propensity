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
#' * `is_stabilized()` returns `TRUE` if the weights are stabilized.
#' * `stabilization_score()` returns the user-supplied stabilization score, or
#'   `NULL` when none was recorded or when a per-observation score was dropped
#'   because an operation changed the length of the weights.
#'
#' ## The stabilization score
#'
#' A `stabilization_score` is the multiplier the weights were stabilized on. It
#' holds either a single value, which scales every weight, or one value per
#' observation, which scales each weight by its own. Both forms are checked
#' where the score is recorded: a score must be numeric, positive, and finite,
#' and hold a length the weights can use, and one that does not is refused with
#' an error of class `propensity_stabilization_score_error`. A score is recorded
#' as a plain double vector, with its storage type normalized and any names
#' dropped, so a score written `1L` and one written `1` are the same score and
#' combine without conflict.
#'
#' A zero-length `psw` is a prototype: it records what a result built from it
#' will carry rather than describing observations of its own, so a
#' per-observation score on one is recorded as given and checked for length when
#' the observations arrive.
#'
#' `psw` objects also inherit the broader `causal_wts` class. The accessors that
#' class carries, [causalgenerics::is_causal_wt()],
#' [causalgenerics::estimand()], and `estimand<-`, are re-exported by propensity
#' and documented at [causalgenerics::estimand()].
#'
#' ## Arithmetic and combining
#'
#' Arithmetic operations on `psw` objects preserve the class and attributes,
#' so operations like normalization (`weights / sum(weights)`) retain metadata.
#'
#' An operation between two `psw` objects merges what each of them records. Two
#' different estimands are pasted together, and an estimand only one operand
#' names stands for the result; the result is stabilized only when both operands
#' are, and it is marked as trimmed, truncated, or calibrated when either
#' operand is. The remaining attributes, the `stabilization_score`, the records
#' left by a modified propensity score, and the attributes describing a
#' categorical exposure, are carried by agreement: one only a single operand
#' records carries, and one both record with the same value carries. One they
#' record differently is dropped, since neither value describes the result, and
#' a single warning of class
#' `propensity_metadata_conflict_warning` names every attribute dropped that
#' way. The rule is applied per operation, so an attribute one operation drops
#' for a disagreement can be carried again by a later operation whose operands
#' agree.
#'
#' A `stabilization_score` is carried only when the result is stabilized, which
#' takes both operands. A result that is not stabilized drops the score without
#' comment.
#'
#' Combining `psw` objects with [c()] preserves the class only when all
#' metadata matches; mismatched metadata produces a warning and falls back to a
#' plain numeric vector. Concatenation appends one set of observations to
#' another, so the positions a modification record names would describe units
#' from the other input; those records are dropped from the result whether or
#' not the inputs agree on them. The categorical attributes name exposure levels
#' rather than positions and carry when the inputs agree.
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
#'   (default: `double()`). For `is_psw()` and `as_psw()`: an object to test or
#'   coerce.
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
#'   the object, either a single value or one value per observation. Every value
#'   must be positive and finite. Defaults to `NULL`, meaning no fixed score was
#'   supplied.
#' @param wt A `psw` or `causal_wts` object.
#' @param ... Additional attributes stored on the object (developer use only).
#'
#' @return
#' * `new_psw()`, `psw()`, `as_psw()`: A `psw` vector.
#' * `is_psw()`, `is_stabilized()`: A single logical value.
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
  stabilization_score <- check_stabilization_score(
    stabilization_score,
    length(x)
  )
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

# `recorded` is the number of observations the record covers, or NULL when
# there is nothing to read that number from.
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
#
# Coverage is the same question here as on a `ps_trim`, and it is asked of the
# same field: the record writes down how many observations its positions were
# written for. A record that does not write that number down, which nothing in
# this package produces, has positions with no sample to place them against and
# is refused like an absent one.
check_trim_meta_covers <- function(x, fn, call = rlang::caller_env()) {
  check_trim_meta(x, fn, call = call)

  meta <- ps_trim_meta(x)
  if (is.null(meta) || record_covers(meta, length(x))) {
    return(invisible(x))
  }

  abort_trim_meta(fn, recorded = meta$n_obs, n = length(x), call = call)
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
  estimand <- merge_psw_estimands(estimand(x), estimand(y))

  # Determine stabilized status: both must be stabilized for result to be stabilized
  stabilized <- is_stabilized(x) && is_stabilized(y)

  # For trimmed/truncated/calibrated: if either has the property, result has it
  trimmed <- is_ps_trimmed(x) || is_ps_trimmed(y)
  truncated <- is_ps_truncated(x) || is_ps_truncated(y)
  calibrated <- is_ps_calibrated(x) || is_ps_calibrated(y)

  rslts <- vec_cast(vec_arith_base(op, x, y), to = double())
  attributes(rslts) <- NULL

  merged <- merge_psw_attrs(x, y, length(rslts))
  attrs <- merged$attrs
  conflicts <- merged$conflicts

  # A score is the value the weights were stabilized on, so a result the merge
  # leaves unstabilized has nothing for it to describe and drops it with the
  # status. That is the same kind of reduction as a drop for length rather than
  # a disagreement between the operands, so it goes without comment and takes
  # any disagreement about the score with it.
  if (!stabilized) {
    attrs["stabilization_score"] <- list(NULL)
    conflicts <- setdiff(conflicts, "stabilization_score")
  }

  if (length(conflicts) > 0) {
    # Arithmetic reaches this method through vctrs' dispatch, whose call would be
    # reported here and names nothing the caller wrote.
    warn_conflicting_psw_attrs(conflicts, call = NULL)
  }

  build_psw(
    rslts,
    estimand = estimand,
    stabilized = stabilized,
    trimmed = trimmed,
    truncated = truncated,
    calibrated = calibrated,
    attrs = attrs
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

# The score a caller supplies, checked where it is recorded. It multiplies the
# weights, so a value that is not a positive finite number produces something
# that is not a set of weights, and a length that is neither one nor one per
# observation recycles into weights nobody asked for. Both are caught here
# rather than left to surface as a strange number much later.
#
# Zero-length data is exempt from the length rule for the reason
# `stabilization_score_aligns()` gives: a prototype carries metadata for
# observations that have not arrived, so there is no length to check against,
# and the restore that builds the real result checks the score there instead.
#
# Returns the score as a double. The metadata merge compares scores with
# `identical()`, which reads a difference in storage type as a difference in the
# score, so normalizing here is what lets a score written `1L` and one written
# `1` describe the same stabilization.
check_stabilization_score <- function(score, n, call = rlang::caller_env()) {
  if (is.null(score)) {
    return(NULL)
  }

  if (!is.numeric(score)) {
    abort(
      c(
        "{.arg stabilization_score} must be numeric.",
        x = "It has class {.cls {class(score)}}."
      ),
      error_class = "propensity_stabilization_score_error",
      call = call
    )
  }

  # A score holding no values passes the alignment rule, which reads any length
  # below two as the single value that scales every weight. There is no such
  # value here, so it is refused separately.
  if (length(score) == 0 || !stabilization_score_aligns(score, n)) {
    abort(
      c(
        "{.arg stabilization_score} must hold one value or one value per
         observation.",
        x = "It holds {length(score)} value{?s}.",
        x = "The weights have {n} observation{?s}.",
        i = "Supply a single value to scale every weight, or one value for
             each observation."
      ),
      error_class = "propensity_stabilization_score_error",
      call = call
    )
  }

  invalid <- !is.finite(score) | score <= 0
  if (any(invalid)) {
    first <- which(invalid)[[1]]
    detail <- if (length(score) == 1) {
      "It is {.val {score}}."
    } else {
      "{sum(invalid)} of its values {?is/are} not, the first at
       position {first}: {.val {score[[first]]}}."
    }

    abort(
      c(
        "{.arg stabilization_score} must be positive and finite.",
        x = detail,
        i = "The score multiplies the weights, so a value that is not positive
             and finite leaves them unusable."
      ),
      error_class = "propensity_stabilization_score_error",
      call = call
    )
  }

  as.double(score)
}

# The records a modified propensity score leaves on the weights built from it.
# Each holds indices into the observations it was recorded on, so like a
# per-observation stabilization score, each is only meaningful at that length.
psw_modification_meta <- c("ps_trim_meta", "ps_trunc_meta", "ps_calib_meta")

# Attributes describing a categorical exposure. These name the exposure levels
# rather than the units, so they carry no length of their own.
psw_categorical_attrs <- c("n_categories", "category_names", "focal_category")

# A modification record is judged against `to`, the object supplying it: the
# record was written for `to`'s observations, so data arriving at `to`'s length
# is the only data those indices describe. `n` is the number of observations the
# data arrives at. Zero-length data is exempt for the reason a stabilization
# score is: it lines up with nothing and so contradicts nothing, and a prototype
# that keeps the record lets the restore building the real result carry it on to
# the observations.
#
# Any same-length operation keeps indices that may no longer point where they
# did, a reordering or a subscript with duplicates included. Nothing rebuilding
# a psw is handed the subscript, so neither can be told apart from an untouched
# vector here. `is_unit_trimmed()` checks the record covers the vector it is
# given, which catches a length change from any route, but a same-length
# rearrangement is beyond what either can see and is documented rather than
# guarded.
modification_meta_aligns <- function(n, to) {
  n == 0 || n == length(to)
}

# Everything a psw carries beyond the six fields describing the weights as a
# whole. The score and the modification records are indexed by observation; the
# categorical attributes name exposure levels and so hold at any length.
psw_carried_attrs <- c(
  "stabilization_score",
  psw_modification_meta,
  psw_categorical_attrs
)

# The attributes a disagreement has already dropped from a prototype under
# construction. vctrs folds `vec_ptype2()` over the inputs two at a time, so a
# prototype meets the third input carrying no trace of the second: without this
# record, an attribute the first pair dropped is carried back by a later pair
# that happens to agree, the result depends on the order the inputs were written
# in, and the warning naming the attribute as dropped is untrue of the result
# the caller gets. An attribute named here is dropped by every remaining pair
# and reported by none of them, so one operation reports each disagreement once
# and answers for it once.
#
# The record belongs to a prototype in the middle of a fold and reaches no
# further: the only route that copies it is the one below, which stops at the
# first vector holding observations.
psw_conflicted_attr <- "psw_conflicted_attrs"

conflicted_psw_attrs <- function(x) {
  out <- attr(x, psw_conflicted_attr)
  if (is.null(out)) character() else out
}

# The carried attributes `to` can still speak for once data has arrived at `n`
# observations. One indexed by observation that does not describe `n` of them is
# reported as absent, so a later read sees nothing rather than a record silently
# misaligned with the weights.
aligned_psw_attrs <- function(to, n) {
  attrs <- lapply(psw_carried_attrs, function(attribute) attr(to, attribute))
  names(attrs) <- psw_carried_attrs

  if (!stabilization_score_aligns(attrs$stabilization_score, n)) {
    attrs["stabilization_score"] <- list(NULL)
  }

  if (!modification_meta_aligns(n, to)) {
    attrs[psw_modification_meta] <- list(NULL)
  }

  # The record of what a fold has already dropped belongs to the prototype being
  # folded and goes no further. Data arriving at observations is the result the
  # caller holds, and the fold that settled its type is over by then. A combined
  # result that holds no observations is the one thing vctrs cannot tell apart
  # from the prototype it was built from, so it keeps the record; it names
  # nothing about observations there are none of.
  if (n == 0) {
    attrs[psw_conflicted_attr] <- list(attr(to, psw_conflicted_attr))
  }

  attrs
}

# The psw type is every attribute the object carries, so anything that rebuilds
# a psw around another vector's data owes all of them: the six fields describing
# the weights as a whole, then the carried attributes, which the caller has
# already reduced to what the result is entitled to.
build_psw <- function(
  x,
  estimand,
  stabilized,
  trimmed,
  truncated,
  calibrated,
  attrs
) {
  out <- new_psw(
    x,
    estimand = estimand,
    stabilized = stabilized,
    trimmed = trimmed,
    truncated = truncated,
    calibrated = calibrated,
    stabilization_score = attrs$stabilization_score
  )

  for (attribute in setdiff(names(attrs), "stabilization_score")) {
    attr(out, attribute) <- attrs[[attribute]]
  }

  out
}

carry_psw_metadata <- function(x, to) {
  build_psw(
    x,
    estimand = estimand(to),
    stabilized = is_stabilized(to),
    trimmed = is_ps_trimmed(to),
    truncated = is_ps_truncated(to),
    calibrated = is_ps_calibrated(to),
    attrs = aligned_psw_attrs(to, length(x))
  )
}

# Two operands naming different estimands describe a result that targets both,
# which is what the pasted label says. An operand naming none has nothing to add
# to the other's label, so the label the one operand supplies stands for the
# result rather than being pasted against an empty half.
merge_psw_estimands <- function(x, y) {
  if (identical(x, y) || is.null(y)) {
    return(x)
  }

  if (is.null(x)) {
    return(y)
  }

  paste0(x, ", ", y)
}

# An operation on two psw objects has two sets of carried attributes to answer
# for, each already reduced to what it can still speak for at the result's
# length. One only a single operand records has nothing to disagree with, which
# is what the single-operand route carries through `w * 1`. One both record
# identically describes the result either way. One they record differently
# describes neither, so it is dropped and named in `conflicts` for the caller to
# report. One an earlier pair of the same operation already dropped stays
# dropped, and is named in `conflicted` rather than `conflicts` so it is not
# reported a second time.
merge_psw_attrs <- function(x, y, n, fields = psw_carried_attrs) {
  x_attrs <- aligned_psw_attrs(x, n)
  y_attrs <- aligned_psw_attrs(y, n)
  dropped <- union(conflicted_psw_attrs(x), conflicted_psw_attrs(y))

  attrs <- list()
  conflicts <- character()

  for (field in fields) {
    x_value <- x_attrs[[field]]
    y_value <- y_attrs[[field]]

    value <- if (field %in% dropped) {
      NULL
    } else if (is.null(x_value) || identical(x_value, y_value)) {
      y_value
    } else if (is.null(y_value)) {
      x_value
    } else {
      conflicts <- c(conflicts, field)
      NULL
    }

    attrs[field] <- list(value)
  }

  list(
    attrs = attrs,
    conflicts = conflicts,
    conflicted = union(dropped, conflicts)
  )
}

# Unlike a drop for length, which happens to records that also travel by routes
# vctrs does not see, a disagreement is about the two objects the caller named
# and is reported once for the operation.
warn_conflicting_psw_attrs <- function(fields, call = rlang::caller_env()) {
  warn(
    c(
      "Dropping the {.field {fields}} attribute{?s} from the result.",
      i = "The two sets of weights record {cli::qty(fields)}{?it/them}
           differently, so neither value describes the result.",
      i = "Every attribute the two agree on is carried through."
    ),
    warning_class = "propensity_metadata_conflict_warning",
    call = call
  )
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

  # The prototype is shared by inputs whose observations are appended one after
  # another, so the positions a modification record names would describe units
  # from the other input. Nothing rebuilding the combined vector is handed the
  # offsets, so the records are left off the prototype whether or not the inputs
  # agree on them. The categorical attributes name exposure levels rather than
  # positions and so mean the same thing at the combined length.
  merged <- merge_psw_attrs(x, y, 0, fields = psw_categorical_attrs)
  if (length(merged$conflicts) > 0) {
    warn_conflicting_psw_attrs(merged$conflicts)
  }

  out <- build_psw(
    double(),
    estimand = estimand(x),
    stabilized = is_stabilized(x),
    trimmed = is_ps_trimmed(x),
    truncated = is_ps_truncated(x),
    calibrated = is_ps_calibrated(x),
    attrs = c(
      list(stabilization_score = stabilization_score(x)),
      merged$attrs
    )
  )

  if (length(merged$conflicted) > 0) {
    attr(out, psw_conflicted_attr) <- merged$conflicted
  }

  out
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
