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
#' A result reports its effects in one of two readings, and `tidy()` returns the
#' one the result records unless `effects` names the other for the call. The
#' marginal reading is the table of causal contrasts described above. The
#' conditional reading is the outcome model's coefficient surface: one row per
#' coefficient, with the standard errors of the block of the joint estimation
#' that carries the uncertainty of having estimated the weights from the same
#' data. Both readings return the same columns in the same order, so their rows
#' stack.
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
#'
#'   The conditional reading has no rows labeled as ratios to pick out, so the
#'   link of the outcome model settles the question there: a `logit` link puts
#'   every coefficient on the log odds scale and a `log` link puts every
#'   coefficient on the log risk scale, and both are scales an exponential
#'   undoes. The estimate and, when they were asked for, the interval bounds are
#'   exponentiated, no term is relabeled, and the standard error, the test
#'   statistic, and the p-value describe the link scale and stay there. Every
#'   other link errors rather than exponentiating coefficients that are not on
#'   such a scale.
#' @param ... These dots are for future extensions and must be empty.
#' @param effects The reading to report, either `"marginal"` or `"conditional"`.
#'   `NULL`, the default, reports the reading the result records; any other value
#'   overrides it for the one call and leaves the result as it is. The marginal
#'   reading reports the population-averaged causal contrasts; the conditional
#'   reading reports the outcome model's coefficient surface. [as_marginal()] and
#'   [as_conditional()] move a result between the two readings.
#'
#'   The covariance the conditional reading reports is the outcome block of the
#'   jointly estimated sandwich, which every route that stacks estimating
#'   equations attaches to the outcome model it stores:
#'   `se_method = "mestimation"` for a binary exposure, and the categorical and
#'   continuous routes, which run on M-estimation alone. A linearization fit
#'   stacks no such system and records no such block, so its conditional reading
#'   errors rather than reporting the covariance the outcome model computed for
#'   itself, which treats the estimated weights as fixed.
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
#' The conditional reading returns those columns in that order, holding one row
#' per coefficient of the outcome model: `term` is the coefficient's name,
#' `estimate` is the coefficient, `std.error` is the square root of the diagonal
#' of the corrected covariance, `statistic` is the estimate over that standard
#' error, and `p.value` is the two-sided normal p-value of the statistic. There
#' is no `comparison` column, because a coefficient is not a contrast of exposure
#' levels, and the bounds are built at the level `conf.level` asks for, the
#' stored ones belonging to the effects the marginal reading reports.
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
#' # The outcome model's coefficients, with the standard errors the joint
#' # estimation of the weights and the outcome implies
#' tidy(result, conf.int = TRUE, effects = "conditional")
#'
#' @seealso [ipw()] for the estimator, [`glance()`][glance.ipw()] for a one-row
#'   summary of the fit, [`augment()`][augment.ipw()] for its per-observation
#'   columns, [`as.data.frame()`][causalgenerics::new_ipw()] for the result's own
#'   columns, and [as_marginal()] and [as_conditional()] for the reading a result
#'   records.
#'
#' @exportS3Method generics::tidy ipw
tidy.ipw <- function(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  ...,
  effects = NULL
) {
  rlang::check_dots_empty()
  check_conf_level(conf.level)

  # The accessors own the `effects` argument, so the estimates of the reading
  # this call reports come from one before anything here branches on a reading.
  # That call resolves a `NULL` against the reading the result records and
  # refuses a value naming neither, on every route a result can have been built
  # by; validating the value here as well would leave two validators to keep in
  # step.
  estimate <- stats::coef(x, effects = effects)

  # The value has been through that check by the time this reads it, so this
  # settles which reading was asked for rather than checking it a second time.
  reading <- if (is.null(effects)) ipw_stored_effects(x) else effects
  if (reading == "conditional") {
    return(tidy_ipw_conditional(
      x,
      estimate = estimate,
      conf.int = conf.int,
      conf.level = conf.level,
      exponentiate = exponentiate
    ))
  }

  tidy_ipw_marginal(
    x,
    conf.int = conf.int,
    conf.level = conf.level,
    exponentiate = exponentiate
  )
}

# The reading a result records. A result built before the field existed carries
# none, and marginal is the reading every method produced then.
ipw_stored_effects <- function(x) {
  if (is.null(x$effects)) "marginal" else x$effects
}

# The marginal reading: the result's own table of causal contrasts, renamed and
# selected into the broom columns.
tidy_ipw_marginal <- function(x, conf.int, conf.level, exponentiate) {
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

# The conditional reading: the outcome model's coefficient surface, in the
# columns the marginal reading uses so that the two tables stack. The estimates
# arrive from the accessor that resolved the reading, and the rest of the row is
# read through the accessors beside it, which is what keeps this table and the
# printed reading of the same result the same numbers.
tidy_ipw_conditional <- function(
  x,
  estimate,
  conf.int,
  conf.level,
  exponentiate,
  call = rlang::caller_env()
) {
  # A row is its estimate and its inference together, so the corrected
  # covariance is read whether or not bounds were asked for. A result that
  # stacked no estimating equations records none, and the error the accessor
  # raises for that is the answer here rather than the covariance the outcome
  # model computed for itself, which treats the estimated weights as fixed.
  covariance <- stats::vcov(x, effects = "conditional")

  # A column of a tibble is addressed by its position, so the names the
  # accessors label their vectors with belong to `term` and to nothing else.
  std_error <- unname(sqrt(diag(covariance)))
  statistic <- unname(estimate) / std_error

  out <- list(
    term = names(estimate),
    estimate = unname(estimate),
    std.error = std_error,
    statistic = statistic,
    # The form causalgenerics prints this reading with, which keeps its
    # significant digits in the far tail where `2 * (1 - pnorm(abs(z)))` loses
    # them.
    p.value = 2 * stats::pnorm(-abs(statistic))
  )

  if (conf.int) {
    # The bounds the result stores describe the effects of the other reading, so
    # there are none to prefer here and the accessor builds them at whatever
    # level the call asked for.
    limits <- stats::confint(x, level = conf.level, effects = "conditional")
    out$conf.low <- unname(limits[, 1L])
    out$conf.high <- unname(limits[, 2L])
  }

  if (exponentiate) {
    check_exponentiate_link(x$outcome_mod, call = call)

    # Broom's convention, which every method taking the argument follows: the
    # estimate and, when they were asked for, the bounds move to the natural
    # scale, and the columns describing the link scale estimate stay there.
    out$estimate <- exp(out$estimate)
    if (conf.int) {
      out$conf.low <- exp(out$conf.low)
      out$conf.high <- exp(out$conf.high)
    }
  }

  tibble::as_tibble(out)
}

# The conditional reading has no rows labeled as ratios to pick out, so the link
# of the outcome model settles whether there is anything for an exponential to
# undo: a logit link puts every coefficient on the log odds scale and a log link
# puts every coefficient on the log risk scale, and a coefficient on any other
# scale exponentiates to a number describing nothing.
check_exponentiate_link <- function(outcome_mod, call = rlang::caller_env()) {
  # An lm is a gaussian identity model, which is how the outcome family check
  # reads one too.
  link <- if (inherits(outcome_mod, "glm")) {
    outcome_mod$family$link
  } else {
    "identity"
  }
  exponentiable <- c("logit", "log")

  if (!link %in% exponentiable) {
    abort(
      c(
        "{.arg exponentiate} needs coefficients on a scale an exponential \\
        undoes.",
        x = "The outcome model was fit with the {.val {link}} link, whose \\
        coefficients are not on such a scale.",
        i = "Only {.val {exponentiable}} links exponentiate in the conditional \\
        reading."
      ),
      error_class = "propensity_exponentiate_error",
      call = call
    )
  }

  invisible(link)
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
#'   \item{`nobs`}{The number of observations the outcome model was fit on,
#'     which is also what the stacked estimating equations are solved on under
#'     M-estimation. Reported by [stats::nobs()].}
#'   \item{`df.residual`}{The residual degrees of freedom of the stacked
#'     estimating equations, `nobs` less the number of parameters the system
#'     solves for. `NA` under linearization, which records no parameter count.
#'     Reported by [stats::df.residual()].}
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

  # Both counts are the result class' own accessors, so the number `glance()`
  # reports and the number `nobs()` or `df.residual()` reports are the same
  # number rather than two readings that agree by construction. The accessors
  # also settle what a result whose `fit` nothing can be read from reports: the
  # residual degrees of freedom are missing, and the rest of the summary stands.
  tibble::tibble(
    estimand = x$estimand,
    nobs = stats::nobs(x),
    df.residual = stats::df.residual(x)
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
#' The weights are reported once, as `.weights`, which is a deliberate departure
#' from broom. A broom `augment()` method carries the model frame's `(weights)`
#' column through and adds no weights column of its own; here the weights are the
#' central per-observation quantity of the analysis, and the [psw()] vector
#' reports them with the estimand they target attached, which the model frame's
#' plain numeric copy of the same values does not. The default frame therefore
#' leaves `(weights)` out. A frame supplied through `data` is the caller's own
#' and is carried as it arrives, so a weight column they hold stays where they
#' put it, whatever it is named.
#'
#' @param x An `ipw` object, as returned by [ipw()].
#' @param data A data frame with one row for each observation the fit used, or
#'   `NULL`, the default, to use the outcome model's own model frame. That frame
#'   holds the response and the terms of the outcome formula, and both are kept.
#'   The `(weights)` column [stats::model.frame()] records is the one thing left
#'   out of it, because those weights are reported as `.weights`.
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
  weights <- stats::model.weights(outcome_frame)

  if (is.null(data)) {
    # The weights are reported as `.weights`, so the frame this method builds
    # for itself leaves out the numeric copy of them `model.frame()` records.
    # The drop belongs to this branch alone: a frame the caller supplies is
    # theirs, and a column of theirs is not deleted for matching a name.
    data <- outcome_frame[names(outcome_frame) != "(weights)"]
  } else {
    assert_class(data, "data.frame")
    check_augment_rows(data, nrow(outcome_frame))
  }

  fitted <- unname(stats::predict(x$outcome_mod, type = "response"))

  columns <- c(
    as.list(data),
    augment_propensity_columns(x$wt_mod),
    list(
      .weights = weights,
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

  # A frame supplied through `data` is the caller's own and may hold
  # non-syntactic names, `(weights)` among them, so the names of the source
  # frame are taken as they are rather than repaired into syntactic ones.
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
