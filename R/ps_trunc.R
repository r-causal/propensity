#' @title Winsorize (Truncate) Propensity Scores
#'
#' @description
#' **`ps_trunc()`** sets out‐of‐range propensity scores to fixed bounding values
#' (a form of *winsorizing*). This is an alternative to [ps_trim()], which removes
#' (sets `NA`) instead of bounding and is then refit with [ps_refit()]
#'
#' @param ps A numeric vector in (0, 1) that represents the probability of receiving `exposure`.
#' @param exposure For method "cr", a 0/1 exposure vector.
#' @param method One of `"ps"`, `"pctl"`, or `"cr"`.
#'   * `"ps"`: directly cut on `[lower, upper]`
#'   * `"pctl"`: use quantiles of `ps` as bounding values
#'   * `"cr"`: for binary exposure, bounding `[min(ps[treated]), max(ps[untreated])]`
#' @param lower,upper Numeric or quantile bounds. If `NULL`, defaults vary by method.
#'
#' @details
#' For each \eqn{ps[i]}:
#'  - If \eqn{ps[i] < lower\_bound}, we set \eqn{ps[i] = lower\_bound}.
#'  - If \eqn{ps[i] > upper\_bound}, we set \eqn{ps[i] = upper\_bound}.
#' This approach is often called *winsorizing*.
#'
#' @return A **`ps_trunc`** object (numeric vector). It has an attribute
#'   `ps_trunc_meta` storing fields like `method`, `lower_bound`, and
#'   `upper_bound`.
#' @seealso [ps_trim()] and [ps_refit()] for removing extreme values vs. bounding
#' @export
ps_trunc <- function(
  ps,
  exposure = NULL,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL
) {
  method <- match.arg(method)
  meta_list <- list(method = method)

  if (any(ps <= 0 | ps >= 1)) {
    stop("All propensity scores must be strictly between 0 and 1.")
  }

  if (method == "ps") {
    if (is.null(lower)) lower <- 0.1
    if (is.null(upper)) upper <- 0.9
    if (lower >= upper) {
      stop("For method='ps', need lower < upper. Got lower=", lower, ", upper=", upper)
    }
    lb <- lower
    ub <- upper
  } else if (method == "pctl") {
    if (is.null(lower)) lower <- 0.05
    if (is.null(upper)) upper <- 0.95
    lb <- quantile(ps, probs = lower)
    ub <- quantile(ps, probs = upper)
    meta_list$lower_pctl <- lower
    meta_list$upper_pctl <- upper
  } else { # method == "cr"
    if (is.null(exposure)) {
      stop("For method='cr', must supply a 0/1 'exposure'.")
    }
    if (!all(exposure %in% c(0, 1))) {
      stop("Exposure must be 0/1 for method='cr'.")
    }
    ps_treat <- ps[exposure == 1]
    ps_untrt <- ps[exposure == 0]
    lb <- min(ps_treat)
    ub <- max(ps_untrt)
  }

  # Winsorize
  ps2 <- pmin(pmax(ps, lb), ub)

  meta <- c(
    meta_list,
    list(
      lower_bound = lb,
      upper_bound = ub
    )
  )

  new_ps_trunc(ps2, meta)
}

new_ps_trunc <- function(x, meta) {
  vec_assert(x, double())
  new_vctr(
    x,
    ps_trunc_meta = meta,
    class = "ps_trunc"
  )
}

#' @title Extract `ps_trunc` metadata
#' @description Returns the internal metadata list for a `ps_trunc` object.
#' @param x A **`ps_trunc`** object.
#' @return A named list of metadata.
#' @export
ps_trunc_meta <- function(x) {
  attr(x, "ps_trunc_meta")
}

#' Check if object is truncated
#'
#' @description
#' `is_truncated()` is an S3 generic that returns `TRUE` if its argument represents a
#' ps_trunc object or psw object flagged as truncated.
#'
#' @param x An R object.
#' @return A logical scalar (`TRUE` or `FALSE`).
#' @export
is_truncated <- function(x) {
  UseMethod("is_truncated")
}


#' @export
is_truncated.default <- function(x) {
  FALSE
}

#' @export
is_truncated.ps_trunc <- function(x) {
  TRUE
}

# vctrs machinery for ps_trunc

#' @export
vec_ptype_abbr.ps_trunc <- function(x, ...) {
  "ps_trunc"
}

#' @export
vec_ptype_full.ps_trunc <- function(x, ...) {
  m <- ps_trunc_meta(x)
  paste0("ps_trunc{[", m$lower_bound, ",", m$upper_bound, "], method=", m$method, "}")
}

# Arithmetic disallowed
#' @export
#' @method vec_arith ps_trunc
vec_arith.ps_trunc <- function(op, x, y, ...) {
  UseMethod("vec_arith.ps_trunc", y)
}

#' @export
#' @method vec_arith.ps_trunc default
vec_arith.ps_trunc.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc ps_trunc
vec_arith.ps_trunc.ps_trunc <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc numeric
vec_arith.ps_trunc.numeric <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_trunc
vec_arith.numeric.ps_trunc <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc integer
vec_arith.ps_trunc.integer <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
vec_math.ps_trunc <- function(.fn, .x, ...) {
  stop_incompatible_op(.fn, .x, subtype = "math")
}

# Combining / Casting
#' @export
vec_ptype2.ps_trunc.ps_trunc <- function(x, y, ...) {
  stop_incompatible_type(x, y, message = "Can't combine ps_trunc objects.")
}

#' @export
vec_ptype2.ps_trunc.double <- function(x, y, ...) double()
#' @export
vec_ptype2.double.ps_trunc <- function(x, y, ...) double()

#' @export
vec_cast.double.ps_trunc <- function(x, to, ...) {
  vec_data(x)
}
#' @export
vec_cast.ps_trunc.double <- function(x, to, ...) {
  new_ps_trunc(
    x,
    meta = list(method = "unknown", lower_bound = NA, upper_bound = NA)
  )
}

#' @export
vec_ptype2.ps_trunc.integer <- function(x, y, ...) integer()
#' @export
vec_ptype2.integer.ps_trunc <- function(x, y, ...) integer()

#' @export
vec_cast.integer.ps_trunc <- function(x, to, ...) as.integer(vec_data(x))
#' @export
vec_cast.ps_trunc.integer <- function(x, to, ...) {
  xx <- as.double(x)
  new_ps_trunc(
    xx,
    meta = list(method = "unknown", lower_bound = NA, upper_bound = NA)
  )
}
