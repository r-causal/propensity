#' Tabulate Truncated Weights Over a Grid of Bounds
#'
#' `wt_trunc_sensitivity()` truncates a set of weights with [wt_trunc()] at
#' each bound in a grid and tabulates what each truncation does to them: the
#' bound as given and as applied, the number of weights moved, and the range
#' and mean of the weights that result. The first row describes the weights
#' before any truncation, so every other row can be read against it. This is
#' the weight side of the sensitivity table Cole and Hernán (2008) recommend
#' for choosing a bound.
#'
#' @param .weights A [psw] vector, such as one returned by [wt_ate()], or a
#'   numeric vector of weights. Propensity scores, including those returned by
#'   [ps_trim()], [ps_trunc()], and [ps_calibrate()], are refused with an error
#'   of class `propensity_type_error`, and weights that have already been
#'   truncated with an error of class `propensity_already_modified_error`.
#' @param method The rule that reads each bound, one of `"pctl"` (the
#'   default), `"wt"`, or `"count"`, as in [wt_trunc()]. The `"adaptive"`
#'   method sets a single bound from the number of weights, so it has no grid
#'   to vary and is not accepted.
#' @param lower An optional grid of lower bounds, read according to `method`
#'   as `wt_trunc()` reads `lower`. A single value is used with every upper
#'   bound; otherwise `lower` must have the same length as `upper`, and the two
#'   are paired element by element.
#' @param upper The grid of upper bounds, read according to `method` as
#'   `wt_trunc()` reads `upper`. For `"pctl"` it defaults to
#'   `c(0.99, 0.975, 0.95, 0.90)`; for `"wt"` and `"count"` it is required.
#'
#' @details
#' Each row after the first is computed by calling [wt_trunc()] with that
#' row's bounds, so it describes exactly the weights that call returns. Every
#' bound in the grid must be one `wt_trunc()` accepts: a single invalid bound
#' refuses the whole grid, with the error `wt_trunc()` would raise for it,
#' rather than leaving a gap in the table. The rows follow the grid in the
#' order given.
#'
#' The summaries are taken over the weights that are present, so missing
#' weights, such as those of units set aside by [ps_trim()], take no part in
#' them and are never counted as truncated. A grid over weights of which none
#' are present is refused with an error of class `propensity_range_error`.
#' Weights can be zero, and a zero weight makes `range_ratio` infinite.
#'
#' Weights that have already been truncated are refused, because the first
#' row would then describe weights that were already bounded rather than the
#' weights before truncation. Pass the weights as they were before
#' `wt_trunc()` instead.
#'
#' ## What is left to the caller
#'
#' The table does not include effect estimates or effective sample sizes.
#' propensity fits no outcome models, so to add estimates to the table, map
#' over the grid, truncate the weights at each bound with `wt_trunc()`, and fit
#' the outcome model with each set of weights. For the effective sample size,
#' add the weights before and after truncation to a data frame as separate
#' columns and pass them to `halfmoon::check_ess()`.
#'
#' @return A [tibble][tibble::tibble] with one row for the untruncated weights
#'   followed by one row per bound, and the columns:
#'   * `lower`, `upper`: the bounds as given, `NA` in the first row and `lower`
#'     `NA` when no lower bound was given.
#'   * `lower_value`, `upper_value`: the bounds as applied, on the scale of the
#'     weights, recorded by `wt_trunc()`. Both are `NA` in the first row, and
#'     `lower_value` is `NA` when no lower bound was given.
#'   * `n_truncated`: the number of weights moved, 0 in the first row.
#'   * `min`, `max`, `mean`: the smallest, largest, and mean weight present.
#'   * `range_ratio`: `max / min`.
#'
#' @references
#' Cole, S. R., & Hernán, M. A. (2008). Constructing inverse probability
#' weights for marginal structural models. *American Journal of
#' Epidemiology*, 168(6), 656--664.
#'
#' @seealso [wt_trunc()] to truncate the weights at a chosen bound.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(2 * x))
#' fit <- glm(z ~ x, family = binomial)
#' w <- wt_ate(fit, .exposure = z)
#'
#' # The default percentile grid. The first row is the untruncated weights.
#' wt_trunc_sensitivity(w)
#'
#' # A grid of absolute bounds, and a lower tail bounded at every one of them
#' wt_trunc_sensitivity(w, method = "wt", upper = c(20, 10, 5))
#' wt_trunc_sensitivity(w, method = "wt", lower = 1.05, upper = c(20, 10, 5))
#'
#' # The largest one, three, and ten weights
#' wt_trunc_sensitivity(w, method = "count", upper = c(1, 3, 10))
#'
#' @export
wt_trunc_sensitivity <- function(
  .weights,
  method = c("pctl", "wt", "count"),
  lower = NULL,
  upper = NULL
) {
  call <- rlang::current_env()
  method <- rlang::arg_match(method, error_call = call)

  check_wt_trunc_sensitivity_weights(.weights, call = call)

  check_wt_trunc_sensitivity_grid(lower, "lower", call = call)
  check_wt_trunc_sensitivity_grid(upper, "upper", call = call)

  if (is.null(upper)) {
    upper <- wt_trunc_sensitivity_default_grid(method, call = call)
  }
  lower <- recycle_wt_trunc_sensitivity_lower(lower, upper, call = call)

  x <- vec_data(.weights)
  present <- x[!is.na(x)]
  if (length(present) == 0L) {
    abort(
      c(
        "{.arg .weights} must have at least one weight present.",
        x = "Of the {length(x)} weight{?s} given, none {?is/are}
             present, so there is no range or mean to report."
      ),
      error_class = "propensity_range_error",
      call = call
    )
  }

  rows <- lapply(seq_along(upper), function(i) {
    wt_trunc_sensitivity_row(
      .weights,
      i = i,
      method = method,
      lower = if (is.null(lower)) NULL else lower[[i]],
      upper = upper[[i]],
      call = call
    )
  })

  reference <- wt_trunc_sensitivity_summary(
    present,
    lower = NA_real_,
    upper = NA_real_,
    lower_value = NA_real_,
    upper_value = NA_real_,
    n_truncated = 0L
  )

  vec_rbind(reference, !!!rows)
}

check_wt_trunc_sensitivity_weights <- function(.weights, call) {
  if (inherits(.weights, c("ps_trim", "ps_trunc", "ps_calib"))) {
    abort_wt_trunc_scores(.weights, call = call)
  }

  if (!inherits(.weights, "psw") && !is.numeric(.weights)) {
    abort(
      c(
        "{.arg .weights} must be a numeric vector of weights.",
        x = "It is {.cls {class(.weights)[[1]]}}, which holds no weights to
             bound.",
        i = "Pass a {.cls psw} vector, such as one {.fun wt_ate} returns."
      ),
      error_class = "propensity_type_error",
      call = call
    )
  }

  if (!is.null(dim(.weights))) {
    abort_wt_trunc_dims(.weights, call = call)
  }

  if (inherits(.weights, "psw") && is_wt_truncated(.weights)) {
    abort(
      c(
        "{.arg .weights} must not already be truncated.",
        x = "The weights have already been bounded with {.fun wt_trunc}, so
             the first row of the grid would describe weights that were
             already bounded rather than the weights before truncation.",
        i = "Pass the weights as they were before {.fun wt_trunc}."
      ),
      error_class = "propensity_already_modified_error",
      call = call
    )
  }

  invisible(TRUE)
}

# A list would be split into bounds of any type by `[[`, so a grid must be an
# atomic vector before its elements are handed to `wt_trunc()`.
check_wt_trunc_sensitivity_grid <- function(value, arg, call) {
  if (is.null(value) || is.atomic(value)) {
    return(invisible(TRUE))
  }

  abort(
    c(
      "{.arg {arg}} must be a numeric vector of bounds.",
      x = "It is {.cls {class(value)[[1]]}}."
    ),
    error_class = "propensity_type_error",
    call = call
  )
}

wt_trunc_sensitivity_default_grid <- function(method, call) {
  if (method == "pctl") {
    return(c(0.99, 0.975, 0.95, 0.90))
  }

  abort(
    c(
      "For {.code method = '{method}'}, {.arg upper} is required.",
      i = "Supply the grid of upper bounds, or use {.code method = 'pctl'} to
           read a default grid of quantiles."
    ),
    error_class = "propensity_missing_arg_error",
    call = call
  )
}

# `lower` is recycled to the length of `upper` and never the other way round:
# `upper` is the grid, and `lower` a bound paired with it.
recycle_wt_trunc_sensitivity_lower <- function(lower, upper, call) {
  if (length(upper) == 0L) {
    abort(
      c(
        "{.arg upper} must have at least one bound.",
        x = "{.arg upper} has length 0, so there is no grid to tabulate."
      ),
      error_class = "propensity_length_error",
      call = call
    )
  }

  if (is.null(lower)) {
    return(NULL)
  }
  if (length(lower) == 1L) {
    return(rep(lower, length(upper)))
  }
  if (length(lower) != length(upper)) {
    abort(
      c(
        "{.arg lower} must have length 1 or the length of {.arg upper}.",
        x = "{.arg lower} has length {length(lower)} and {.arg upper} has
             length {length(upper)}."
      ),
      error_class = "propensity_length_error",
      call = call
    )
  }

  lower
}

# One truncated row. `wt_trunc()` validates the bounds and applies them, and a
# refusal it raises is reported against the grid the caller asked for.
wt_trunc_sensitivity_row <- function(
  .weights,
  i,
  method,
  lower,
  upper,
  call
) {
  truncated <- rlang::try_fetch(
    wt_trunc(.weights, method = method, lower = lower, upper = upper),
    propensity_error = function(cnd) {
      cnd$call <- call
      # The position rather than the argument: a recycled `lower` can be the
      # bound refused as well as `upper`.
      cnd$body <- c(
        cnd$body,
        i = cli::format_inline(
          "The refused bound is at position {i} of the grid."
        )
      )
      rlang::cnd_signal(cnd)
    }
  )

  meta <- attr(truncated, "psw_trunc_meta")
  values <- vec_data(truncated)

  wt_trunc_sensitivity_summary(
    values[!is.na(values)],
    lower = na_if_null(meta$lower),
    upper = meta$upper,
    lower_value = na_if_null(meta$lower_value),
    upper_value = meta$upper_value,
    n_truncated = length(meta$truncated_idx)
  )
}

wt_trunc_sensitivity_summary <- function(
  present,
  lower,
  upper,
  lower_value,
  upper_value,
  n_truncated
) {
  # The names a psw carries identify units, not rows of the table.
  present <- unname(present)
  tibble::tibble(
    lower = lower,
    upper = upper,
    lower_value = lower_value,
    upper_value = upper_value,
    n_truncated = n_truncated,
    min = min(present),
    max = max(present),
    range_ratio = max(present) / min(present),
    mean = mean(present)
  )
}

na_if_null <- function(x) {
  if (is.null(x)) NA_real_ else x
}
