#' Tidy an inverse probability weighted result
#'
#' @description
#' `tidy()` returns the estimates of an [ipw()] result as a tibble using the
#' column names broom conventions use. There is one row per effect measure, and
#' one row per effect measure per comparison for a categorical exposure, in the
#' order the result stores them. Nothing is dropped.
#'
#' The values are the ones the result already holds: `tidy()` renames and
#' selects rather than re-estimating anything. The one exception is the
#' confidence interval, which is rebuilt from the estimate and its standard
#' error when the requested `conf.level` differs from the level the result was
#' fit at.
#'
#' @param x An `ipw` object, as returned by [ipw()].
#' @param conf.int Logical. Should the confidence interval bounds be returned in
#'   the `conf.low` and `conf.high` columns? Defaults to `FALSE`.
#' @param conf.level The confidence level of the interval, a single number
#'   between 0 and 1. Defaults to `0.95`, which is this method's own default
#'   rather than the level `x` was fit at. When the two differ, the bounds are
#'   recomputed as the normal approximation
#'   `estimate +/- qnorm(1 - (1 - conf.level) / 2) * std.error`, which is the
#'   interval [ipw()] itself reports. Ignored when `conf.int` is `FALSE`.
#' @param exponentiate Logical. Should the estimate and its bounds be
#'   exponentiated on the rows reported on the log scale? Defaults to `FALSE`.
#'   This behaves exactly as it does in
#'   [`as.data.frame()`][causalgenerics::new_ipw()]: the `log(rr)` and `log(or)`
#'   rows have their estimate and bounds exponentiated and their `term` relabeled
#'   to `rr` and `or`, while every other row is left alone. The standard error,
#'   the test statistic, and the p-value describe the log scale estimate and stay
#'   there.
#' @param ... These dots are for future extensions and must be empty.
#'
#' @return A [tibble][tibble::tibble] with one row per estimate and the columns:
#' \describe{
#'   \item{`term`}{The effect measure, such as `"rd"`, `"log(rr)"`, `"log(or)"`,
#'     `"diff"`, or `"slope"`.}
#'   \item{`comparison`}{The contrast the row reports, such as `"b vs a"`.
#'     Categorical exposures only.}
#'   \item{`estimate`}{The estimated effect.}
#'   \item{`std.error`}{The standard error of the estimate.}
#'   \item{`statistic`}{The z statistic, the estimate over its standard error.}
#'   \item{`p.value`}{The two-sided p-value of that statistic.}
#'   \item{`conf.low`, `conf.high`}{The interval bounds. Present only when
#'     `conf.int` is `TRUE`.}
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' x1 <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.5 * x1))
#' y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.3 * x1))
#' dat <- data.frame(x1, z, y)
#'
#' ps_mod <- glm(z ~ x1, data = dat, family = binomial())
#' wts <- wt_ate(ps_mod)
#' outcome_mod <- glm(y ~ z, data = dat, family = quasibinomial(), weights = wts)
#' result <- ipw(ps_mod, outcome_mod)
#'
#' tidy(result)
#'
#' # A 90% interval, with the ratio measures on their natural scale
#' tidy(result, conf.int = TRUE, conf.level = 0.9, exponentiate = TRUE)
#'
#' @seealso [ipw()] for the estimator, and
#'   [`as.data.frame()`][causalgenerics::new_ipw()] for the result's own
#'   columns.
#'
#' @exportS3Method generics::tidy ipw
tidy.ipw <- function(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  ...
) {
  rlang::check_dots_empty()
  check_conf_level(conf.level)

  if (conf.int) {
    x$estimates <- recompute_ipw_interval(x$estimates, conf.level)
  }

  # Exponentiation is the result's own contract, so the coercion method applies
  # it, to whichever interval this call settled on.
  estimates <- as.data.frame(x, exponentiate = exponentiate)

  out <- list(
    term = estimates$effect,
    estimate = estimates$estimate,
    std.error = estimates$std.err,
    statistic = estimates$z,
    p.value = estimates$p.value
  )

  # A categorical exposure reports each effect measure once per contrast, and
  # the column naming the contrast follows the term it qualifies.
  if (!is.null(estimates[["comparison"]])) {
    out <- append(out, list(comparison = estimates$comparison), after = 1L)
  }

  if (conf.int) {
    out$conf.low <- estimates$ci.lower
    out$conf.high <- estimates$ci.upper
  }

  tibble::as_tibble(out)
}

# The estimates table with its interval expressed at `conf_level`. The stored
# bounds are returned untouched when they already describe that level, so a
# request for the level the result was fit at reports the result's own numbers
# rather than a recomputation that merely agrees with them.
recompute_ipw_interval <- function(estimates, conf_level) {
  stored <- estimates[["conf.level"]]
  if (!is.null(stored) && isTRUE(all(stored == conf_level))) {
    return(estimates)
  }

  half_width <- stats::qnorm(1 - (1 - conf_level) / 2) * estimates$std.err
  estimates$ci.lower <- estimates$estimate - half_width
  estimates$ci.upper <- estimates$estimate + half_width
  estimates$conf.level <- conf_level

  estimates
}

check_conf_level <- function(
  conf.level,
  arg = rlang::caller_arg(conf.level),
  call = rlang::caller_env()
) {
  valid <- is.numeric(conf.level) &&
    length(conf.level) == 1L &&
    !is.na(conf.level) &&
    conf.level > 0 &&
    conf.level < 1

  if (!valid) {
    # An empty value has nothing for `{.val}` to print, so it is described by
    # its type instead.
    supplied <- if (length(conf.level) == 0) {
      "You supplied {.obj_type_friendly {conf.level}}."
    } else {
      "You supplied {.val {conf.level}}."
    }

    abort(
      c(
        "{.arg {arg}} must be a single number between 0 and 1.",
        x = supplied,
        i = "Use {.code {arg} = 0.95} for a 95% interval."
      ),
      error_class = "propensity_conf_level_error",
      call = call
    )
  }

  invisible(conf.level)
}
