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
#' @seealso [ipw()] for the estimator, [`glance()`][glance.ipw()] for a one-row
#'   summary of the fit, [`augment()`][augment.ipw()] for its per-observation
#'   columns, and [`as.data.frame()`][causalgenerics::new_ipw()] for the
#'   result's own columns.
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

#' Glance at an inverse probability weighted result
#'
#' @description
#' `glance()` describes an [ipw()] result rather than its estimates: one row
#' naming the estimand and counting the observations and the residual degrees of
#' freedom of the system the standard errors came from. A fit reporting several
#' effect measures, or several comparisons of a categorical exposure, still
#' returns exactly one row.
#'
#' Under M-estimation that system is the stacked estimating equations, which
#' hold the propensity score model, the outcome model, and the effect measures
#' at once. Its residual degrees of freedom are the observations it was solved
#' on less the parameters it solves for, so a multinomial propensity score model
#' leaves fewer of them than a binary one does on the same data. Linearization
#' stacks nothing and records no parameter count, so the observations are the
#' outcome model's and there is no count to subtract from them.
#'
#' The columns and their types are the same on every route [ipw()] takes, so the
#' rows of several results stack into one table.
#'
#' @param x An `ipw` object, as returned by [ipw()].
#' @param ... These dots are for future extensions and must be empty.
#'
#' @return A one-row [tibble][tibble::tibble] with the columns:
#' \describe{
#'   \item{`estimand`}{The estimand the weights target, such as `"ate"`.}
#'   \item{`nobs`}{The number of observations the standard errors were estimated
#'     from: the observations the stacked estimating equations were solved on
#'     under M-estimation, and the observations the outcome model was fit on
#'     under linearization.}
#'   \item{`df.residual`}{The residual degrees of freedom of the stacked
#'     estimating equations, `nobs` less the number of parameters the system
#'     solves for. `NA` under linearization, which records no parameter count.}
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
#' glance(result)
#'
#' @seealso [ipw()] for the estimator, [`tidy()`][tidy.ipw()] for its estimates,
#'   and [`augment()`][augment.ipw()] for its per-observation columns.
#'
#' @exportS3Method generics::glance ipw
glance.ipw <- function(x, ...) {
  rlang::check_dots_empty()

  # The counts describe whatever produced the standard errors. That is the
  # stacked estimating equations when the result holds a fit, and the outcome
  # model alone otherwise, which leaves no parameter count to spend the
  # observations against.
  if (is.null(x$fit)) {
    nobs <- as.integer(stats::nobs(x$outcome_mod))
    df_residual <- NA_integer_
  } else {
    nobs <- as.integer(stats::nobs(x$fit))
    df_residual <- as.integer(stats::df.residual(x$fit))
  }

  tibble::tibble(
    estimand = x$estimand,
    nobs = nobs,
    df.residual = df_residual
  )
}

#' Augment an inverse probability weighted result with per-observation columns
#'
#' @description
#' `augment()` works per observation rather than per estimate: the data an
#' [ipw()] result was produced from, carried through in full, with the propensity
#' score, the weights, the fitted values, and the residuals attached as
#' dot-prefixed columns. Nothing is dropped and no column of the source frame is
#' changed.
#'
#' The added columns are read from the two models the result holds, so they
#' describe the fit rather than the frame they are attached to. `data` says which
#' frame to carry them on, which is how broom uses the argument: it names the
#' data the fit was produced from, not new data to predict on. Supplying the
#' modeling data rather than leaving the default is how covariates the outcome
#' formula left out arrive beside the fit's own columns.
#'
#' @param x An `ipw` object, as returned by [ipw()].
#' @param data A data frame with one row for each observation the fit used, or
#'   `NULL`, the default, to use the outcome model's own model frame. That frame
#'   holds the response, the terms of the outcome formula, and the `(weights)`
#'   column [stats::model.frame()] records, and all of them are kept.
#' @param ... These dots are for future extensions and must be empty.
#'
#' @return A [tibble][tibble::tibble] with one row per observation, holding every
#' column of the source frame in its own order followed by:
#' \describe{
#'   \item{`.propensity`}{The propensity score, as the propensity score model
#'     predicts it: the probability of exposure for a binary exposure, and the
#'     conditional mean of the exposure for a continuous one. This is the
#'     propensity score the weight functions take, so feeding it back with the
#'     exposure returns the weights the outcome model was fit with.}
#'   \item{`.propensity_<level>`}{For a categorical exposure, which has a
#'     probability for every level rather than a single number, one column per
#'     level in place of `.propensity`, named for the level and ordered as the
#'     propensity score model orders its levels.}
#'   \item{`.weights`}{The weights the outcome model was fit with, as the
#'     [psw()] vector they were supplied as, so the estimand they record travels
#'     with them.}
#'   \item{`.fitted`}{The outcome model's fitted values, on the response scale.}
#'   \item{`.resid`}{The observed outcome less `.fitted`. Present only when the
#'     source frame holds the outcome the model was fit on, which the default
#'     frame always does. A factor or logical outcome is differenced on the 0/1
#'     scale its fitted values are on.}
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
#' augment(result)
#'
#' # The modeling data carries through the covariates the outcome formula left
#' # out
#' augment(result, data = dat)
#'
#' @seealso [ipw()] for the estimator, [`tidy()`][tidy.ipw()] for its estimates,
#'   and [`glance()`][glance.ipw()] for a one-row summary of the fit.
#'
#' @exportS3Method generics::augment ipw
augment.ipw <- function(x, data = NULL, ...) {
  rlang::check_dots_empty()

  outcome_frame <- stats::model.frame(x$outcome_mod)

  if (is.null(data)) {
    data <- outcome_frame
  } else {
    assert_class(data, "data.frame")
    check_augment_rows(data, nrow(outcome_frame))
  }

  fitted <- unname(stats::predict(x$outcome_mod, type = "response"))

  columns <- c(
    as.list(data),
    augment_propensity_columns(x$wt_mod),
    list(
      .weights = stats::model.weights(outcome_frame),
      .fitted = fitted
    )
  )

  # A residual is an observed outcome less a fitted value, so it belongs to a
  # frame that holds the outcome and to no other. The outcome goes through the
  # conversion the estimator itself uses, which is what puts the residual of a
  # factor or logical outcome on the 0/1 scale its fitted values are on.
  response <- names(outcome_frame)[[
    attr(stats::terms(x$outcome_mod), "response")
  ]]
  if (response %in% names(data)) {
    columns$.resid <- ipw_outcome_numeric(data[[response]]) - fitted
  }

  # The source frame is carried verbatim, `(weights)` included, so its names are
  # taken as they are rather than repaired into syntactic ones.
  tibble::as_tibble(columns, .name_repair = "minimal")
}

# The propensity score of each observation, as the columns it takes to hold it.
# A binary or continuous exposure has one number per observation and so one
# column; a categorical exposure has a probability for every level, and each
# level gets a column of its own.
augment_propensity_columns <- function(wt_mod) {
  if (!inherits(wt_mod, "multinom")) {
    return(list(
      .propensity = unname(stats::predict(wt_mod, type = "response"))
    ))
  }

  levels <- wt_mod$lev
  probs <- stats::predict(wt_mod, type = "probs")

  # nnet::multinom() predicts a two-level exposure as the probability of the
  # second level alone rather than as the matrix of level probabilities it
  # returns for more levels.
  if (!is.matrix(probs)) {
    probs <- cbind(1 - probs, probs)
    colnames(probs) <- levels
  }

  rlang::set_names(
    lapply(levels, function(level) unname(probs[, level])),
    paste0(".propensity_", levels)
  )
}

check_augment_rows <- function(
  data,
  n_obs,
  arg = rlang::caller_arg(data),
  call = rlang::caller_env()
) {
  if (nrow(data) == n_obs) {
    return(invisible(data))
  }

  abort(
    c(
      "{.arg {arg}} must have one row for each observation the fit used.",
      x = "{.arg {arg}} has {nrow(data)} row{?s}, but the fit used \\
      {n_obs} observation{?s}.",
      i = "{.arg {arg}} names the data the fit was produced from. Leave it as \\
      {.code NULL} to use the outcome model's own frame."
    ),
    error_class = "propensity_augment_data_error",
    call = call
  )
}
