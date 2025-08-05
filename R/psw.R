#' Create and Manipulate `psw` Objects
#'
#' Functions to create and manipulate `psw` objects, which are specialized
#' vectors for propensity score weights with optional `estimand` attributes.
#' Most users should use [`wt_ate()`] and friends, but these functions can help
#' extend the functionality of `psw` objects.
#'
#' @import vctrs
#' @export
#'
#' @param x A numeric vector (default: `double()`).
#' @param estimand A character string representing the estimand (e.g., "ate",
#'   "att", "ato"). Default is `NULL`.
#' @param stabilized A logical `TRUE`
#' @param trimmed Logical, whether these weights came from a trimmed PS.
#' @param truncated Logical, whether these weights came from a truncated PS.
#' @param calibrated Logical, whether these weights came from a calibrated PS.
#' @param wt An object to check or convert.
#' @param value The value to add to the attribute.
#' @param ... Additional attributes to track in the weights.
#' @return
#' - `new_psw()`: A `psw` object.
#' - `psw()`: A `psw` object.
#' - `is_psw()`: `TRUE` if the object is a `psw`, otherwise `FALSE`.
#' - `as_psw()`: A `psw` object.
#' - `estimand()`: The `estimand` attribute of a `psw` object.
#' - `is_stabilized()`: The `stabilized` attribute of a `psw` object.
#'
#' @details
#' The `psw` class is a vctrs-based S3 class that represents propensity score
#' weights. It extends numeric vectors with additional metadata tracking the
#' estimand type, stabilization status, and source transformations.
#'
#' **Arithmetic behavior**: Unlike `ps_trim` and `ps_trunc` objects, arithmetic
#' operations on `psw` objects preserve the class and attributes. This allows
#' weight manipulations like normalization (`weights / sum(weights)`) while
#' maintaining metadata.
#'
#' **Combining behavior**: When combining `psw` objects with `c()`, the class
#' is preserved only if all metadata matches. Mismatched metadata triggers a
#' warning and returns a numeric vector.
#'
#' **Base R compatibility**: Most base R operations work seamlessly:
#' - Subsetting with `[` preserves class and attributes
#' - Summary functions (`sum()`, `mean()`, etc.) return numeric values
#' - Comparison operators return logical vectors
#' - Works in data frames and with tidyverse functions
#'
#' @examples
#' psw_weights <- new_psw(c(0.1, 0.2, 0.3), estimand = "ate")
#' is_psw(psw_weights)
#' estimand(psw_weights)
#'
#' psw_helper <- psw(c(0.5, 0.7), estimand = "att")
#' as_psw(c(0.1, 0.2), estimand = "ato")
#'
#' # Coercion behavior - compatible objects combine silently
#' x <- psw(c(0.5, 0.7), estimand = "ate")
#' y <- psw(c(0.3, 0.8), estimand = "ate")
#' c(x, y)  # Returns psw object
#'
#' # Incompatible metadata triggers warning and returns numeric
#' x <- psw(c(0.5, 0.7), estimand = "ate")
#' y <- psw(c(0.3, 0.8), estimand = "att")
#' c(x, y)  # Warning: returns numeric
#'
#' # Works with tidyr::pivot_longer for plotting
#' if (requireNamespace("tidyr", quietly = TRUE)) {
#'   df <- data.frame(
#'     id = 1:4,
#'     ate_wts = psw(c(0.5, 0.7, 0.3, 0.8), estimand = "ate"),
#'     att_wts = psw(c(0.4, 0.6, 0.2, 0.9), estimand = "att")
#'   )
#'   # This will warn but succeed, returning numeric in the pivoted column
#'   tidyr::pivot_longer(df, cols = c(ate_wts, att_wts))
#' }
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
  calibrated = FALSE
) {
  x <- vec_cast(x, to = double())
  attributes(x) <- NULL
  new_psw(
    x,
    estimand = estimand,
    stabilized = stabilized,
    trimmed = trimmed,
    truncated = truncated,
    calibrated = calibrated
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
is_causal_wt <- function(x) {
  inherits(x, "causal_wts")
}

#' @rdname psw
#' @export
as_psw <- function(x, estimand = NULL) {
  x <- vec_cast(x, to = double())
  psw(x, estimand = estimand)
}

#' @rdname psw
#' @export
estimand <- function(wt) {
  attr(wt, "estimand")
}

#' @rdname psw
#' @export
`estimand<-` <- function(wt, value) {
  assert_class(wt, "causal_wts")
  attr(wt, "estimand") <- value
  wt
}

#' @export
is_ps_trimmed.psw <- function(x) {
  isTRUE(attr(x, "trimmed"))
}

#' @export
is_unit_trimmed.psw <- function(x) {
  out <- vector("logical", length = length(x))
  if (!is_ps_trimmed(x)) {
    return(out)
  }

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
  meta <- ps_trim_meta(x)
  if (!is.null(meta)) {
    return(isTRUE(meta$refit))
  }
  FALSE
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

  rslts <- vec_arith_base(op, x, y)
  psw(rslts, estimand = estimand)
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

#' @export
vec_math.psw <- function(.fn, .x, ...) {
  # Some functions like cumsum/cumprod should preserve psw class
  if (.fn %in% c("cumsum", "cumprod", "cummin", "cummax")) {
    result <- vec_math_base(.fn, vec_data(.x), ...)
    return(vec_restore(result, .x))
  }
  # Other functions like log, sqrt return numeric
  vec_math_base(.fn, vec_data(.x), ...)
}

#' @export
Summary.psw <- function(..., na.rm = FALSE) {
  .fn <- .Generic
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call(.fn, c(numeric_args, list(na.rm = na.rm)))
}

#' @export
min.psw <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("min", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
max.psw <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("max", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
range.psw <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("range", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
median.psw <- function(x, na.rm = FALSE, ...) {
  median(vec_data(x), na.rm = na.rm, ...)
}

#' @export
vec_restore.psw <- function(x, to, ...) {
  # Extract numeric data if needed
  if (inherits(x, "psw")) {
    x <- vec_data(x)
  }

  # Preserve psw attributes
  new_psw(
    x,
    estimand = estimand(to),
    stabilized = is_stabilized(to),
    trimmed = is_ps_trimmed(to),
    truncated = is_ps_truncated(to),
    calibrated = isTRUE(attr(to, "calibrated"))
  )
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
    truncated = is_ps_truncated(x)
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

#' @export
vec_cast.psw.psw <- function(x, to, ...) x

#' @export
vec_cast.psw.double <- function(x, to, ...) psw(x, estimand = estimand(to))

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
vec_cast.psw.integer <- function(x, to, estimand = NULL, ...)
  psw(x, estimand = estimand(to))

#' @export
vec_cast.integer.psw <- function(x, to, ...) {
  vec_cast(vec_data(x), integer(), x_arg = "psw")
}

#' @export
vec_cast.psw.ps_trim <- function(x, to, ...) {
  psw(vec_data(x), estimand = estimand(to))
}

#' @export
vec_cast.ps_trim.psw <- function(x, to, ...) {
  ps_trim(vec_data(x), method = "ps", lower = 0, upper = 1)
}

#' @export
vec_cast.psw.ps_trunc <- function(x, to, ...) {
  psw(vec_data(x), estimand = estimand(to))
}

#' @export
vec_cast.ps_trunc.psw <- function(x, to, ...) {
  ps_trunc(vec_data(x), method = "ps", lower = 0, upper = 1)
}

#' @export
summary.psw <- function(object, ...) {
  summary(as.numeric(object), ...)
}

#' @export
quantile.psw <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, ...) {
  quantile(vec_data(x), probs = probs, na.rm = na.rm, ...)
}

#' @export
anyDuplicated.psw <- function(x, incomparables = FALSE, ...) {
  anyDuplicated(vec_data(x), incomparables = incomparables, ...)
}

#' @export
diff.psw <- function(x, lag = 1L, differences = 1L, ...) {
  diff(vec_data(x), lag = lag, differences = differences, ...)
}
