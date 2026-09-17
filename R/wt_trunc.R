#' Truncate (Winsorize) Propensity Score Weights
#'
#' `wt_trunc()` bounds extreme weights at a limit read on the scale of the
#' weights themselves, replacing each weight beyond the limit with the limit.
#' No unit is removed. Where [ps_trunc()] bounds the propensity scores the
#' weights are built from, `wt_trunc()` bounds the weights after they are
#' built, so it applies to weights for any exposure type, including the
#' density-ratio weights of a continuous exposure, which have no propensity
#' score to bound.
#'
#' @param .weights A [psw] vector, such as one returned by [wt_ate()], or a
#'   numeric vector of weights. A numeric vector is returned as a [psw] vector
#'   with no estimand. Propensity scores, including those returned by
#'   [ps_trim()], [ps_trunc()], and [ps_calibrate()], are refused with an error
#'   of class `propensity_type_error`: build weights from them first.
#' @param method The rule that sets the bound:
#'   * `"adaptive"` (default): an upper bound of \eqn{\sqrt{n} \log(n) / 5}
#'     (Gruber et al., 2022), where \eqn{n} is the number of weights present.
#'   * `"wt"`: bounds given as weight values.
#'   * `"pctl"`: bounds read at sample quantiles of the weights.
#'   * `"count"`: bounds that move a given number of the most extreme weights.
#' @param lower,upper The bounds, read according to `method`:
#'   * `"adaptive"`: not used. Supplying either is ignored with a warning.
#'   * `"wt"`: `upper` is required and must be a positive finite number.
#'     `lower` is optional and must be a finite number of at least 0 and below
#'     `upper`.
#'   * `"pctl"`: `upper` is required and must lie in (0.5, 1). `lower` is
#'     optional and must lie in (0, 0.5).
#'   * `"count"`: `upper` is required and must be a whole number of at least 1;
#'     the `upper` largest weights are bounded at the next largest. `lower` is
#'     optional and bounds the `lower` smallest weights at the next smallest in
#'     the same way.
#' @param ... Not used. Any argument supplied here is an error, so that a
#'   misspelled bound is not silently ignored.
#'
#' @details
#' ## Methods
#'
#' `"adaptive"` is the default. Its bound depends on the sample size alone and
#' loosens as the sample grows: at \eqn{n = 1000} it is 43.7. Gruber et al.
#' (2022) derived it for the propensity score of a binary treatment, choosing a
#' bound that loosens with the sample so that its bias shrinks. Whether the
#' rule carries over to the density-ratio weights of a continuous exposure has
#' not been established in the literature. Here \eqn{n} counts the weights that
#' are present, so for weights built from a trimmed propensity score, whose
#' trimmed units are `NA`, it is the number of units kept by the trim. The
#' bound is below 1 for \eqn{n \le 6} (0.72 at \eqn{n = 5}), so for a sample
#' that small `"adaptive"` pulls every weight above that value down to a
#' constant below 1.
#'
#' `"wt"` winsorizes at the values given, and `"count"` at the order statistic
#' next to the weights it bounds. Under `"count"`, a weight tied with that
#' order statistic is not moved, so fewer weights than `upper` may change. The
#' two counts must leave at least two weights between them, and `upper` alone
#' must leave at least one weight below it.
#'
#' `"pctl"` reads its bounds with [stats::quantile()] at its default type, as
#' `ps_trunc(method = "pctl")` does. Because such a rule moves a fixed share of
#' the units however the weights are distributed, it is best used to check how
#' sensitive an analysis is to its bound rather than as the bound itself. An
#' `upper` below 0.5 is refused rather than read as the other tail: write
#' `upper = 0.9` to bound the upper tail at the 0.9 quantile, or
#' `lower = 0.1` to bound the lower tail at the 0.1 quantile.
#'
#' `"adaptive"` and `"pctl"` need at least two weights present, and every
#' method leaves a missing weight missing.
#'
#' ## One-sided by default
#'
#' Each method bounds only the upper tail unless `lower` is supplied, and
#' `"adaptive"`, which takes no `lower`, always bounds only the upper tail.
#' The largest weights are the ones that dominate an estimate; small weights
#' contribute little to it, so bounding them changes little. `lower` bounds the
#' lower tail when that is wanted.
#'
#' ## What is recorded
#'
#' The result keeps every attribute of `.weights`: the stabilization status and
#' any numerator model, the density record of a continuous exposure, the
#' attributes of a categorical exposure, and the records left by [ps_trim()],
#' [ps_trunc()], or [ps_calibrate()]. It gains a `wt_truncated` flag, read by
#' [is_wt_truncated()], and a record of the method, the bounds as given and as
#' applied, and the positions of the weights moved, read by
#' [is_unit_wt_truncated()]. The record is positional and is dropped by any
#' operation that changes the length of the weights; see [psw].
#'
#' The estimand gains the label `"; weights truncated"`, as in
#' `"ate; weights truncated"`. A truncated weight targets a different quantity
#' from the weight before truncation: any fixed bound introduces a bias (Ma and
#' Wang, 2020; Gruber et al., 2022), and an interval computed with the
#' truncated weights held fixed is an interval for that quantity. A weight
#' truncation does not describe the same operation as truncating the propensity
#' score, which is labeled `"; truncated"` and reported by [is_ps_truncated()].
#'
#' Weights that have already been truncated are returned unchanged with a
#' warning of class `propensity_already_modified_warning`. One record describes
#' one truncation, so to truncate at a different bound, call `wt_trunc()` on
#' the original weights.
#'
#' ## Weights from a trimmed, truncated, or calibrated score
#'
#' Weights built from a score modified by [ps_trim()], [ps_trunc()], or
#' [ps_calibrate()] can be truncated as well. The result carries both records,
#' and the labels stack, as in `"ate; trimmed; weights truncated"`. The package
#' does not prevent combining these tools.
#'
#' ## Effect estimation
#'
#' [ipw()] does not accept weights whose values were truncated, because its
#' standard errors do not account for the bound. An analysis with truncated
#' weights needs an interval built another way: one computed with the weights
#' held fixed, which conditions on the bound applied, or a bootstrap that
#' repeats the truncation in every resample, which a percentile bound in
#' particular needs, since that bound is itself estimated from the weights.
#'
#' To see how much a truncation changed the effective sample size, add the
#' weights before and after it to a data frame as separate columns and pass
#' both to `halfmoon::check_ess()`.
#'
#' @return A [psw] vector of the same length as `.weights`, with the bounded
#'   values, the estimand labeled as truncated, and a record of the truncation.
#'
#' @references
#' Gruber, S., Phillips, R. V., Lee, H., & van der Laan, M. J. (2022).
#' Data-adaptive selection of the propensity score truncation level for
#' inverse-probability-weighted and targeted maximum likelihood estimators of
#' marginal point treatment effects. *American Journal of Epidemiology*,
#' 191(9), 1640--1651.
#'
#' Ma, X., & Wang, J. (2020). Robust inference using inverse probability
#' weighting. *Journal of the American Statistical Association*, 115(532),
#' 1851--1860.
#'
#' @seealso [is_wt_truncated()] and [is_unit_wt_truncated()] to read the
#'   record, [ps_trunc()] to bound propensity scores instead, and [ps_trim()]
#'   to remove units.
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(2 * x))
#' fit <- glm(z ~ x, family = binomial)
#' w <- wt_ate(fit, .exposure = z)
#'
#' # The adaptive bound, on the upper tail
#' w_adaptive <- wt_trunc(w)
#' w_adaptive
#' is_wt_truncated(w_adaptive)
#' which(is_unit_wt_truncated(w_adaptive))
#'
#' # A bound given as a weight, on one tail or both
#' wt_trunc(w, method = "wt", upper = 5)
#' wt_trunc(w, method = "wt", lower = 1.05, upper = 5)
#'
#' # The 99th percentile, as a sensitivity check
#' wt_trunc(w, method = "pctl", upper = 0.99)
#'
#' # The three largest weights, bounded at the fourth largest
#' wt_trunc(w, method = "count", upper = 3)
#'
#' # Weights for a continuous exposure
#' a <- 1 + 0.5 * x + rnorm(n)
#' dose <- lm(a ~ x)
#' w_dose <- wt_ate(dose, .exposure = a)
#' wt_trunc(w_dose)
#'
#' # Weights from a trimmed propensity score carry both records
#' trimmed <- ps_refit(ps_trim(fit, lower = 0.05, upper = 0.95), fit)
#' wt_trunc(wt_ate(trimmed, .exposure = z), method = "wt", upper = 5)
#'
#' @export
wt_trunc <- function(
  .weights,
  method = c("adaptive", "wt", "pctl", "count"),
  lower = NULL,
  upper = NULL,
  ...
) {
  UseMethod("wt_trunc")
}

#' @export
wt_trunc.default <- function(
  .weights,
  method = c("adaptive", "wt", "pctl", "count"),
  lower = NULL,
  upper = NULL,
  ...
) {
  abort(
    c(
      "{.arg .weights} must be a numeric vector of weights.",
      x = "It is {.cls {class(.weights)[[1]]}}, which holds no weights to
           bound.",
      i = "Pass a {.cls psw} vector, such as one {.fun wt_ate} returns."
    ),
    error_class = "propensity_type_error",
    call = rlang::current_env()
  )
}

#' @export
wt_trunc.numeric <- function(
  .weights,
  method = c("adaptive", "wt", "pctl", "count"),
  lower = NULL,
  upper = NULL,
  ...
) {
  call <- rlang::current_env()
  rlang::check_dots_empty(call = call)
  method <- rlang::arg_match(method, error_call = call)

  # `psw()` drops every attribute of what it is given, names included, and the
  # names identify the units the weights belong to.
  weights <- psw(.weights)
  names(weights) <- names(.weights)

  truncate_weights(weights, method, lower, upper, call = call)
}

#' @export
wt_trunc.psw <- function(
  .weights,
  method = c("adaptive", "wt", "pctl", "count"),
  lower = NULL,
  upper = NULL,
  ...
) {
  call <- rlang::current_env()
  rlang::check_dots_empty(call = call)
  method <- rlang::arg_match(method, error_call = call)

  if (is_wt_truncated(.weights)) {
    warn(
      c(
        "These weights have already been truncated. Returning them unchanged.",
        i = "To truncate at a different bound, call {.fun wt_trunc} on the
             original weights."
      ),
      warning_class = "propensity_already_modified_warning",
      call = call
    )
    return(.weights)
  }

  truncate_weights(.weights, method, lower, upper, call = call)
}

#' @export
wt_trunc.ps_trim <- function(.weights, ...) {
  abort_wt_trunc_scores(.weights, call = rlang::current_env())
}

#' @export
wt_trunc.ps_trunc <- function(.weights, ...) {
  abort_wt_trunc_scores(.weights, call = rlang::current_env())
}

#' @export
wt_trunc.ps_calib <- function(.weights, ...) {
  abort_wt_trunc_scores(.weights, call = rlang::current_env())
}

# A modified propensity score is a vector of numbers between 0 and 1, and
# bounding it as though it were a weight would bound the wrong quantity without
# complaint.
abort_wt_trunc_scores <- function(.weights, call = rlang::caller_env()) {
  abort(
    c(
      "{.arg .weights} must be weights, not propensity scores.",
      x = "It is {.cls {class(.weights)[[1]]}}, a modified propensity score.",
      i = "Build weights from it with a weight function such as {.fun wt_ate}
           and truncate those, or bound the scores themselves with
           {.fun ps_trunc}."
    ),
    error_class = "propensity_type_error",
    call = call
  )
}

# Winsorize a psw at the bounds `method` reads from `lower` and `upper`, and
# write the result back over the original so that every attribute the
# truncation does not describe is kept.
truncate_weights <- function(.weights, method, lower, upper, call) {
  # A psw holds doubles, and `vec_data()` keeps the names the weights carry.
  x <- vec_data(.weights)
  present <- x[!is.na(x)]

  bounds <- switch(
    method,
    adaptive = wt_trunc_adaptive_bounds(present, lower, upper, call = call),
    wt = wt_trunc_wt_bounds(lower, upper, call = call),
    pctl = wt_trunc_pctl_bounds(present, lower, upper, call = call),
    count = wt_trunc_count_bounds(present, lower, upper, call = call)
  )

  # A missing weight compares as neither above nor below a bound, and `which()`
  # drops it, so it is left missing and is not among the weights moved.
  low <- integer()
  if (!is.null(bounds$lower_value)) {
    low <- which(x < bounds$lower_value)
    x[low] <- bounds$lower_value
  }
  high <- which(x > bounds$upper_value)
  x[high] <- bounds$upper_value

  meta <- new_psw_trunc_meta(
    method = method,
    lower = bounds$lower,
    upper = bounds$upper,
    lower_value = bounds$lower_value,
    upper_value = bounds$upper_value,
    truncated_idx = sort(c(low, high)),
    n_obs = length(x)
  )

  out <- vec_restore(x, .weights)
  estimand <- estimand(.weights)
  if (!is.null(estimand)) {
    attr(out, "estimand") <- paste0(estimand, "; weights truncated")
  }
  attr(out, "wt_truncated") <- TRUE
  attr(out, "psw_trunc_meta") <- meta

  # The restore rebuilds the attributes in its own order. Putting them back in
  # the order the weights arrived with leaves the attributes the truncation did
  # not touch exactly as they were, with the record following them.
  attrs <- attributes(out)
  arrived <- intersect(names(attributes(.weights)), names(attrs))
  attributes(out) <- attrs[union(arrived, names(attrs))]

  out
}

wt_trunc_adaptive_bounds <- function(present, lower, upper, call) {
  if (!is.null(lower) || !is.null(upper)) {
    warn(
      c(
        "For {.code method = 'adaptive'}, {.arg lower} and {.arg upper} are
         ignored.",
        i = "The adaptive bound is set by the number of weights. To give a
             bound yourself, use {.code method = 'wt'}."
      ),
      call = call
    )
  }

  check_wt_trunc_n_present(present, "adaptive", call = call)
  n <- length(present)

  list(
    lower = NULL,
    upper = NULL,
    lower_value = NULL,
    upper_value = sqrt(n) * log(n) / 5
  )
}

wt_trunc_wt_bounds <- function(lower, upper, call) {
  upper <- wt_trunc_bound_arg(upper, "upper", "wt", call = call)
  lower <- wt_trunc_bound_arg(lower, "lower", "wt", call = call)

  if (upper <= 0) {
    abort_wt_trunc_range(
      "For {.code method = 'wt'}, {.arg upper} must be a positive weight.",
      "{.arg upper} is {upper}.",
      call = call
    )
  }

  if (!is.null(lower)) {
    if (lower < 0) {
      abort_wt_trunc_range(
        "For {.code method = 'wt'}, {.arg lower} must be zero or more.",
        "{.arg lower} is {lower}, and a weight is never negative.",
        call = call
      )
    }
    if (lower >= upper) {
      abort_wt_trunc_range(
        "{.arg lower} must be smaller than {.arg upper}.",
        "{.arg lower} is {lower} and {.arg upper} is {upper}.",
        call = call
      )
    }
  }

  list(lower = lower, upper = upper, lower_value = lower, upper_value = upper)
}

wt_trunc_pctl_bounds <- function(present, lower, upper, call) {
  upper <- wt_trunc_bound_arg(upper, "upper", "pctl", call = call)
  lower <- wt_trunc_bound_arg(lower, "lower", "pctl", call = call)

  # Read as a lower tail, a probability below one half names a real request,
  # so both readings are offered rather than one chosen for the caller.
  if (upper > 0 && upper < 0.5) {
    abort(
      c(
        "For {.code method = 'pctl'}, {.arg upper} must be between 0.5 and 1.",
        x = "{.arg upper} is {upper}, which is in the lower tail.",
        i = "To bound the upper tail write {.code upper = {1 - upper}}; to
             bound the lower tail write {.code lower = {upper}}."
      ),
      error_class = "propensity_range_error",
      call = call
    )
  }
  if (!(upper > 0.5 && upper < 1)) {
    abort_wt_trunc_range(
      "For {.code method = 'pctl'}, {.arg upper} must be between 0.5 and 1.",
      "{.arg upper} is {upper}.",
      call = call
    )
  }
  if (!is.null(lower) && !(lower > 0 && lower < 0.5)) {
    abort_wt_trunc_range(
      "For {.code method = 'pctl'}, {.arg lower} must be between 0 and 0.5.",
      "{.arg lower} is {lower}.",
      call = call
    )
  }

  check_wt_trunc_n_present(present, "pctl", call = call)

  # `quantile()` names its result for the probability, which says nothing
  # about the bound and would reappear wherever the bound is compared.
  lower_value <- NULL
  if (!is.null(lower)) {
    lower_value <- unname(quantile(present, probs = lower))
  }

  list(
    lower = lower,
    upper = upper,
    lower_value = lower_value,
    upper_value = unname(quantile(present, probs = upper))
  )
}

wt_trunc_count_bounds <- function(present, lower, upper, call) {
  upper <- wt_trunc_bound_arg(upper, "upper", "count", call = call)
  lower <- wt_trunc_bound_arg(lower, "lower", "count", call = call)

  check_wt_trunc_count(upper, "upper", call = call)
  check_wt_trunc_count(lower, "lower", call = call)

  n <- length(present)
  if (is.null(lower) && upper > n - 1) {
    abort_wt_trunc_range(
      "For {.code method = 'count'}, {.arg upper} must leave a weight below
       the ones it bounds.",
      "{.arg upper} is {upper} and {n} weight{?s} {?is/are} present.",
      call = call
    )
  }
  if (!is.null(lower) && lower + upper > n - 2) {
    abort_wt_trunc_range(
      "For {.code method = 'count'}, {.arg lower} and {.arg upper} must leave
       at least two weights between them.",
      "{.arg lower} is {lower}, {.arg upper} is {upper}, and {n} weight{?s}
       {?is/are} present.",
      call = call
    )
  }

  sorted <- sort(present)
  lower_value <- NULL
  if (!is.null(lower)) {
    lower_value <- sorted[[lower + 1]]
  }

  list(
    lower = lower,
    upper = upper,
    lower_value = lower_value,
    upper_value = sorted[[n - upper]]
  )
}

check_wt_trunc_count <- function(value, arg, call) {
  if (is.null(value) || (value >= 1 && value == round(value))) {
    return(invisible(TRUE))
  }

  abort_wt_trunc_range(
    "For {.code method = 'count'}, {.arg {arg}} must be a whole number of at
     least 1.",
    "{.arg {arg}} is {value}.",
    call = call
  )
}

# A bound as the caller gave it, checked for the things every method needs and
# stored as a double: the record's bounds are compared with `identical()` when
# truncated weights are combined, so an integer and a double giving the same
# bound must be recorded the same way.
wt_trunc_bound_arg <- function(value, arg, method, call) {
  if (is.null(value)) {
    if (arg == "lower") {
      return(NULL)
    }
    abort(
      c(
        "For {.code method = '{method}'}, {.arg upper} is required.",
        i = "Supply the bound, or use {.code method = 'adaptive'} to have it
             set from the number of weights."
      ),
      error_class = "propensity_missing_arg_error",
      call = call
    )
  }

  if (anyNA(value)) {
    abort(
      c(
        "{.arg {arg}} must not be missing.",
        i = "A missing bound decides nothing about which weights to move.
             Supply a value, or leave {.arg lower} unset to bound only the
             upper tail."
      ),
      error_class = "propensity_missing_value_error",
      call = call
    )
  }

  if (!is.numeric(value)) {
    abort(
      c(
        "{.arg {arg}} must be a number.",
        x = "It is {.cls {class(value)[[1]]}}."
      ),
      error_class = "propensity_type_error",
      call = call
    )
  }

  if (length(value) == 0L) {
    abort_wt_trunc_range(
      "{.arg {arg}} must be a single finite number.",
      "{.arg {arg}} has length 0.",
      call = call
    )
  }
  if (length(value) != 1L || !is.finite(value)) {
    abort_wt_trunc_range(
      "{.arg {arg}} must be a single finite number.",
      "{.arg {arg}} has length {length(value)} and value{?s} {value}.",
      call = call
    )
  }

  as.double(value)
}

check_wt_trunc_n_present <- function(present, method, call) {
  n <- length(present)
  if (n >= 2) {
    return(invisible(TRUE))
  }

  abort_wt_trunc_range(
    "For {.code method = '{method}'}, at least two weights must be present.",
    "{n} weight{?s} {?is/are} present.",
    call = call
  )
}

abort_wt_trunc_range <- function(
  problem,
  detail,
  call,
  .envir = parent.frame()
) {
  abort(
    c(problem, x = detail),
    error_class = "propensity_range_error",
    call = call,
    .envir = .envir
  )
}

# One line describing the truncation for the printed footer, or `NULL` for
# weights that were not truncated. The counts come from the record's positions,
# so a record that no longer covers the weights is not read for them.
format_psw_trunc_line <- function(x) {
  if (!is_wt_truncated(x)) {
    return(NULL)
  }

  meta <- attr(x, "psw_trunc_meta")
  if (!record_covers(meta, length(x))) {
    return("truncation: weights truncated")
  }

  given <- c(meta$lower, meta$upper)
  label <- meta$method
  if (meta$method %in% c("pctl", "count")) {
    label <- paste(
      label,
      paste(vapply(given, format, character(1)), collapse = "/")
    )
  }

  realized <- c(
    if (!is.null(meta$lower_value)) {
      paste("lower", format(meta$lower_value, digits = 3))
    },
    paste("upper", format(meta$upper_value, digits = 3))
  )

  paste0(
    "truncation: ",
    label,
    " (",
    paste(realized, collapse = ", "),
    "), ",
    length(meta$truncated_idx),
    " of ",
    meta$n_obs,
    " weights truncated"
  )
}
