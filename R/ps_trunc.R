#' Truncate (Winsorize) Propensity Scores
#'
#' **`ps_trunc()`** sets out‐of‐range propensity scores to fixed bounding values
#' (a form of *winsorizing*). This is an alternative to [ps_trim()], which removes
#' (sets `NA`) instead of bounding and is then refit with [ps_refit()]
#'
#' @param ps The propensity score, a numeric vector between 0 and 1.
#' @param .exposure For method "cr", a binary exposure vector.
#' @param method One of `"ps"`, `"pctl"`, or `"cr"`.
#'   * `"ps"`: directly cut on `[lower, upper]` of `ps`.
#'   * `"pctl"`: use quantiles of `ps` as bounding values
#'   * `"cr"`: the common range of `ps` given `.exposure`, bounding `[min(ps[treated]), max(ps[untreated])]`
#' @param lower,upper Numeric or quantile bounds. If `NULL`, defaults vary by method.
#' @inheritParams wt_ate
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
#'
#' @examples
#' set.seed(2)
#' n <- 30
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' # truncate just the 99th percentile
#' ps_trunc(ps, method = "pctl", lower = 0, upper = .99)
#'
#' @export
ps_trunc <- function(
  ps,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .treated = NULL,
  .untreated = NULL
) {
  method <- rlang::arg_match(method)
  meta_list <- list(method = method)

  check_ps_range(ps)

  if (method == "ps") {
    if (is.null(lower)) lower <- 0.1
    if (is.null(upper)) upper <- 0.9
    check_lower_upper(lower, upper)

    lb <- lower
    ub <- upper
  } else if (method == "pctl") {
    if (is.null(lower)) lower <- 0.05
    if (is.null(upper)) upper <- 0.95
    lb <- quantile(ps, probs = lower)
    ub <- quantile(ps, probs = upper)
    meta_list$lower_pctl <- lower
    meta_list$upper_pctl <- upper
  } else {
    .exposure <- transform_exposure_binary(
      .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
    ps_treat <- ps[.exposure == 1]
    ps_untrt <- ps[.exposure == 0]
    lb <- min(ps_treat)
    ub <- max(ps_untrt)
  }

  # winsorize
  pinned_low <- which(ps < lb)
  pinned_high <- which(ps > ub)
  truncated_idx <- sort(c(pinned_low, pinned_high))

  ps[pinned_low] <- lb
  ps[pinned_high] <- ub

  meta <- c(
    meta_list,
    list(
      lower_bound = lb,
      upper_bound = ub,
      truncated_idx = truncated_idx
    )
  )

  new_ps_trunc(ps, meta)
}

new_ps_trunc <- function(x, meta) {
  vec_assert(x, double())
  new_vctr(
    x,
    ps_trunc_meta = meta,
    class = "ps_trunc",
    inherit_base_type = TRUE
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
#' @description `is_ps_truncated()` is an S3 generic that returns `TRUE` if its
#' argument represents a `ps_trunc` object or `psw` object created from
#' truncated propensity scores. `is_ps_truncated()` is a question about whether
#' the propensity scores *have* been truncated, as opposed to
#' [is_unit_truncated()], which is a question about which *units* have been
#' truncated.
#'
#' @param x An object.
#' @return A logical scalar (`TRUE` or `FALSE`).
#' @export
is_ps_truncated <- function(x) {
  UseMethod("is_ps_truncated")
}

#' @export
is_ps_truncated.default <- function(x) {
  FALSE
}

#' @export
is_ps_truncated.ps_trunc <- function(x) {
  TRUE
}

#' Check if units have been truncated
#'
#' @description `is_ps_truncated()` is an S3 generic that returns a vector of
#'   `TRUE` or `FALSE`, representing if the element has been truncated.
#'   [is_unit_truncated()] is a question about which *units* have been
#'   truncated, as opposed to `is_ps_truncated()`, which is a question about
#'   whether the propensity scores *have* been truncated.
#'
#' @param x An object.
#' @return A logical vector.
#' @export
is_unit_truncated <- function(x) {
  UseMethod("is_unit_truncated")
}

#' @export
is_unit_truncated.default <- function(x) {
  abort(
    "{.code is_unit_truncated()} not supported for class {.val {class(x)}}"
  )
}

#' @export
is_unit_truncated.ps_trunc <- function(x) {
  meta <- ps_trunc_meta(x)
  out <- vector("logical", length = length(x))
  out[meta$truncated_idx] <- TRUE

  out
}


#' Check if object is truncated
#'
#' @description
#' `is_ps_truncated()` is an S3 generic that returns `TRUE` if its argument represents a
#' ps_trunc object or psw object flagged as truncated.
#'
#' @param x An R object.
#' @return A logical scalar (`TRUE` or `FALSE`).
#' @export
is_ps_truncated <- function(x) {
  UseMethod("is_ps_truncated")
}

#' @export
is_ps_truncated.default <- function(x) {
  FALSE
}

#' @export
is_ps_truncated.ps_trunc <- function(x) {
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
  paste0(
    "ps_trunc{[",
    m$lower_bound,
    ",",
    m$upper_bound,
    "], method=",
    m$method,
    "}"
  )
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
vec_ptype2.psw.ps_trunc <- function(x, y, ...) character()

#' @export
vec_ptype2.ps_trunc.psw <- function(x, y, ...) character()

#' @export
vec_cast.character.ps_trunc <- function(x, to, ...) as.character(vec_data(x))

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

#' @export
vec_math.ps_trunc <- function(.fn, .x, ...) {
  vec_math_base(.fn, vec_data(.x), ...)
}
