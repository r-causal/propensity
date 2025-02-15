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
#' @examples
#' psw_weights <- new_psw(c(0.1, 0.2, 0.3), estimand = "ate")
#' is_psw(psw_weights)
#' estimand(psw_weights)
#'
#' psw_helper <- psw(c(0.5, 0.7), estimand = "att")
#' as_psw(c(0.1, 0.2), estimand = "ato")
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
    class = c("psw", "causal_wts")
  )
}



#' @rdname psw
#' @export
psw <- function(
  x = double(),
  estimand = NULL,
  stabilized = FALSE,
  trimmed = FALSE,
  truncated = FALSE
) {
  x <- vec_cast(x, to = double())
  attributes(x) <- NULL
  new_psw(
    x,
    estimand   = estimand,
    stabilized = stabilized,
    trimmed    = trimmed,
    truncated  = truncated
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
  vec_cast(x, to = new_psw(estimand = estimand))
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
is_trimmed.psw <- function(x) {
  isTRUE(attr(x, "trimmed"))
}

#' @export
is_truncated.psw <- function(x) {
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

vec_ptype_full.psw <- function(x, ...) {
  estimand <- estimand(x)
  if (is.null(estimand)) {
    "psw{estimand = unknown}"
  } else {
    paste0("psw{estimand = ", estimand, "}")
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
#' @method vec_arith.psw numeric
vec_arith.psw.numeric <- function(op, x, y, ...) {
  vec_arith_base(op, x, y) |>
    psw(estimand = estimand(x))
}

#' @export
#' @method vec_arith.numeric psw
vec_arith.numeric.psw <- function(op, x, y, ...) {
  vec_arith_base(op, x, y) |>
    psw(estimand = estimand(y))
}

#' @export
#' @method vec_arith.psw integer
vec_arith.psw.integer <- function(op, x, y, ...) {
  estimand <- estimand(x)

  vec_arith_base(op, x, y) |>
    new_psw(estimand = estimand)
}

#' @export
vec_math.psw <- function(.fn, .x, ...) {
  vec_math_base(.fn, vec_data(.x), ...)
}

#' @export
vec_ptype2.psw.psw <- function(x, y, ...) {
  if (!identical(estimand(x), estimand(y))) {
    stop_incompatible_type(
      x,
      y,
      x_arg = "x",
      y_arg = "y",
      message = "Can't combine weights with different estimands"
    )
  }
  new_psw(estimand = estimand(x))
}

#' @export
vec_ptype2.psw.double <- function(x, y, ...) double()

#' @export
vec_ptype2.double.psw <- function(x, y, ...) double()

#' @export
vec_cast.psw.double <- function(x, to, ...) psw(x, estimand = estimand(to))

#' @export
vec_cast.double.psw <- function(x, to, ...) vec_data(x)

#' @export
vec_ptype2.psw.integer <- function(x, y, ...) integer()

#' @export
vec_ptype2.integer.psw <- function(x, y, ...) integer()

#' @export
vec_cast.psw.integer <- function(x, to, estimand = NULL, ...) psw(x, estimand = estimand(to))

#' @export
vec_cast.integer.psw <- function(x, to, ...) {
  vec_cast(vec_data(x), integer(), x_arg = "psw")
}
