#' Tidy an inverse probability weighted result
#'
#' @description
#' `tidy()` returns the estimates of an [ipw()] result as a tibble using the
#' column names broom conventions use. For a binary or categorical exposure
#' there is one row per exposure level, under the term `"mean"`, and then one
#' row per effect measure per contrast of levels; for a continuous exposure
#' there is one row per reported coefficient. The rows arrive in the order the
#' result stores them and nothing is dropped.
#'
#' Those columns are the ones the result's own coercion surface reports, so the
#' marginal reading is [`as.data.frame()`][causalgenerics::new_ipw()] read as a
#' tibble rather than a second assembly of the same table. The values are the
#' ones the result already holds, nothing being re-estimated: the confidence
#' interval is the only thing rebuilt, and only when the requested `conf.level`
#' differs from the level the result was fit at. The one thing the frame carries
#' that the tibble does not is the covariance of the effects it attaches as an
#' attribute, a tidied table being its columns.
#'
#' `conf.int`, `conf.level`, and `exponentiate` are the arguments this method
#' shares with that surface, so the surface validates them, in both readings and
#' before either is assembled. A value it refuses is reported under the
#' causalgenerics condition naming the argument at fault, which is what keeps one
#' argument of one method from being well formed in one reading of a result and
#' refused in the other.
#'
#' A result reports its effects in one of two readings, and `tidy()` returns the
#' one the result records unless `effects` names the other for the call. The
#' marginal reading is the table of causal contrasts described above. The
#' conditional reading is the outcome model's coefficient surface: one row per
#' coefficient, with the standard errors of the block of the joint estimation
#' that carries the uncertainty of having estimated the weights from the same
#' data. The two readings return the same columns in the same order, with one
#' exception: the `contrast` column belongs to the marginal reading. A binary
#' and a categorical result both carry it there, naming the exposure level each
#' mean row belongs to and the pair of levels each effect measure compares, and
#' the conditional reading of either returns the same table one column narrower.
#' A coefficient names neither a level nor a pair of them.
#'
#' @param x An `ipw` object, as returned by [ipw()].
#' @param conf.int Logical. Should the confidence interval bounds be returned in
#'   the `conf.low` and `conf.high` columns? Defaults to `FALSE`.
#' @param conf.level The confidence level of the interval, a single number
#'   between 0 and 1. Defaults to `0.95`, which is this method's own default
#'   rather than the level `x` was fit at. When the two differ, the bounds are
#'   recomputed as the normal approximation
#'   `estimate +/- qnorm(1 - (1 - conf.level) / 2) * std.error`, which is the
#'   interval [ipw()] itself reports. Not used when `conf.int` is `FALSE`, though
#'   it must still be a valid level.
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
#'   `"fixed"` is also accepted, and asks for the reading the result records,
#'   exactly as `NULL` does. It is the name mixed models give the effects that
#'   are not random, and [mice::pool()] asks every tidier it calls for it. A
#'   result reports one reading at a time and neither of them is random, so the
#'   name resolves to whichever one the result records. The accessors do not
#'   accept it: it is spelled here, where multiple imputation reaches the
#'   result, rather than in the surface underneath.
#'
#'   The covariance the conditional reading reports is the outcome block of the
#'   jointly estimated sandwich, which every route that stacks estimating
#'   equations attaches to the outcome model it stores:
#'   `se_method = "mestimation"` for a binary exposure, and the categorical and
#'   continuous routes, which run on M-estimation alone. A linearization fit
#'   stacks no such system and records no such block, so its conditional reading
#'   errors rather than reporting the covariance the outcome model computed for
#'   itself, which treats the estimated weights as fixed.
#' @param parametric Accepted and ignored. [mice::pool()] passes
#'   `parametric = TRUE` to every tidier it calls, to ask the models it was
#'   written against for their parametric coefficient table rather than for a
#'   smooth term. An `ipw` result reports one table, and the inference it
#'   reports is already parametric, so there is nothing for the argument to
#'   select and nothing is read from it. It is here so that a result can be
#'   pooled across multiply imputed datasets.
#'
#' @return A [tibble][tibble::tibble] with one row per estimate and the columns:
#' \describe{
#'   \item{`term`}{The effect measure, such as `"mean"` for a counterfactual
#'     mean, or `"rd"`, `"log(rr)"`, `"log(or)"`, `"diff"`, `"slope"`, or
#'     `"coef"` for the measures built from them.}
#'   \item{`contrast`}{What the row is about: the exposure level a `"mean"` row
#'     belongs to, such as `"0"` or `"b"`; the pair of levels an effect measure
#'     compares, such as `"1 vs 0"` or `"b vs a"`; or the coefficient name on a
#'     surface whose rows are named after their coefficients. Present only where
#'     a row reports one.}
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
#' is no `contrast` column, whatever the exposure, because a coefficient names
#' neither an exposure level nor a pair of them, and the bounds are built at the
#' level `conf.level` asks for, the stored ones belonging to the effects the
#' marginal reading reports.
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
#'   records. For results fitted to multiply imputed data, see the Multiple
#'   imputation section of [ipw()], [pool_ipw()] for the pooling, and
#'   [`tidy()`][tidy.ipw_pooled()] for what it returns.
#'
#' @exportS3Method generics::tidy ipw
tidy.ipw <- function(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  ...,
  effects = NULL,
  parametric = NULL
) {
  rlang::check_dots_empty()

  # `"fixed"` is the reading the result records, which is what `NULL` asks for,
  # so it is resolved here and the accessors below never see it. Compared
  # exactly rather than matched loosely, so that a value close to it is still
  # refused by the accessors that own the argument.
  if (identical(effects, "fixed")) {
    effects <- NULL
  }

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

  # The conditional reading exponentiates by the link of the outcome model
  # rather than by the labels of the rows, and the refusal of a link no
  # exponential undoes is this package's, so it is raised before the surface
  # below reaches a refusal of its own. A value of `exponentiate` that is not a
  # flag is not read here; it is refused there, with the other two arguments.
  if (reading == "conditional" && isTRUE(exponentiate)) {
    check_exponentiate_link(x$outcome_mod)
  }

  # The result's own coercion surface, in the reading this call reports rather
  # than the one the result records: the three arguments the two readings share
  # are validated there, and one argument of one method cannot be well formed in
  # one reading of a result and refused in the other.
  #
  # The conditional reading answers from the accessors and so discards the
  # frame, which costs a table of a few rows: nothing is refit to build it.
  frame <- as.data.frame(
    x,
    conf.int = conf.int,
    conf.level = conf.level,
    exponentiate = exponentiate,
    effects = reading
  )

  if (reading == "conditional") {
    return(tidy_ipw_conditional(
      x,
      estimate = estimate,
      conf.int = conf.int,
      conf.level = conf.level,
      exponentiate = exponentiate
    ))
  }

  # The covariance of the effects the frame attaches is the one thing that does
  # not travel: it is an attribute rather than a column, and a tidied table is
  # its columns.
  out <- tibble::as_tibble(frame)
  attr(out, "ipw_vcov") <- NULL

  out
}

# The reading a result records. A result built before the field existed carries
# none, and marginal is the reading every method produced then.
ipw_stored_effects <- function(x) {
  if (is.null(x$effects)) "marginal" else x$effects
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
  exponentiate
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

#' Glance at an inverse probability weighted result
#'
#' @description
#' `glance()` describes an [ipw()] result rather than its estimates: one row
#' naming the estimand and counting the observations and the residual degrees of
#' freedom of the system the standard errors came from. A fit reporting several
#' effect measures, or several contrasts of a categorical exposure, still
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

#' Tidy a pooled inverse probability weighted result
#'
#' @description
#' `tidy()` returns the estimates of a [pool_ipw()] result as a tibble using the
#' column names broom conventions use. It holds the rows the results themselves
#' report: for a binary or categorical exposure one row per exposure level,
#' under the term `"mean"`, and then one row per effect measure per contrast of
#' levels; for a continuous exposure one row per reported coefficient. The rows
#' arrive in the order the pooled result stores them and nothing is dropped.
#'
#' Those columns are the ones the pooled result's own coercion surface reports,
#' so this method is that surface read as a tibble rather than a second assembly
#' of the same table. The covariance of the effects the frame attaches is the one
#' thing that does not travel: it is an attribute rather than a column, and a
#' tidied table is its columns.
#'
#' The pooling has already happened by the time this method sees the result.
#' Nothing here re-pools and no argument changes an estimate: they say how the
#' table already in the object is reported, which is which of the two readings
#' it is, whether an interval comes with it, the level that interval is reported
#' at, and the scale the ratio rows are put on.
#'
#' A pooled result reports its effects in one of two readings, and `tidy()`
#' returns the one the result records unless `effects` names the other for the
#' call. [pool_ipw()] pools both readings of the results it is given, so naming
#' one here reports a table the pooling already built rather than pooling again.
#' The marginal reading is the table of pooled causal contrasts described above.
#' The conditional reading is the outcome models' coefficient surface, pooled
#' the same way: one row per coefficient, from the block of the joint estimation
#' that carries the uncertainty of having estimated the weights from the same
#' data. The two readings return the same columns in the same order, with one
#' exception: the `contrast` column belongs to the marginal reading. A binary
#' and a categorical pool both carry it there, naming the exposure level each
#' mean row belongs to and the pair of levels each effect measure compares, and
#' the conditional reading of either returns the same table one column narrower.
#'
#' @param x An `ipw_pooled` object, as returned by [pool_ipw()].
#' @param conf.int Logical. Should the confidence interval bounds be returned in
#'   the `conf.low` and `conf.high` columns? Defaults to `FALSE`.
#' @param conf.level The confidence level of the interval. `NULL`, the default,
#'   reports the bounds at the level the pooled result stored, which is the level
#'   its own intervals were built at. A number between 0 and 1 rebuilds them at
#'   that level instead, from the estimate and its standard error referred to t
#'   on the pooled degrees of freedom. Not used when `conf.int` is `FALSE`,
#'   though it must still be a valid level.
#' @param exponentiate Logical. Should the estimate and its bounds be
#'   exponentiated on the rows reported on the log scale? Defaults to `FALSE`.
#'   This behaves exactly as it does on the pooled result's own coercion surface:
#'   in the marginal reading the `log(rr)` and `log(or)` rows have their estimate
#'   and bounds exponentiated and their `term` relabeled to `rr` and `or`, and in
#'   the conditional reading the outcome model's link settles the scale for every
#'   row at once and no term is relabeled. In both readings the standard error,
#'   the test statistic, and the p-value describe the log scale and stay there,
#'   which is where the inference was done.
#' @param ... These dots are for future extensions and must be empty.
#' @param effects The reading to report, either `"marginal"` or `"conditional"`.
#'   `NULL`, the default, reports the reading the pooled result records; naming
#'   a reading reports that one for the one call and leaves the result as it is.
#'   The marginal reading reports the pooled causal contrasts; the conditional
#'   reading reports the pooled coefficients of the outcome models, which is the
#'   reading `exponentiate` above puts on the natural scale every row at once.
#'   [as_marginal()] and [as_conditional()] move a pooled result between the two
#'   readings for good rather than for a call.
#'
#'   A reading the pooling could not compute is refused rather than reported
#'   under the other one's name. The conditional reading needs the covariance
#'   the joint estimation of the weights and the outcome implies, and a
#'   linearization fit stacks no such system and records no such block, so a set
#'   of those fits pools the marginal reading alone. A result pooled before both
#'   readings were kept holds only the one it was pooled in, and asking it for
#'   the other is refused the same way.
#'
#'   `"fixed"` is not among the values this method takes, though
#'   [`tidy()`][tidy.ipw()] on a single result accepts it. [mice::pool()] asks
#'   for that reading from the tidier of each imputation's own fit, which is
#'   where the name arrives and where the tolerance for it belongs; nothing
#'   tidies a pooled result on its behalf.
#'
#' @return A [tibble][tibble::tibble] with one row per pooled estimate and the
#'   columns:
#' \describe{
#'   \item{`term`}{The effect measure, such as `"mean"` for a counterfactual
#'     mean, or `"rd"`, `"log(rr)"`, `"log(or)"`, `"diff"`, `"slope"`, or
#'     `"coef"` for the measures built from them, or the coefficient name in the
#'     conditional reading.}
#'   \item{`contrast`}{What the row is about: the exposure level a `"mean"` row
#'     belongs to, such as `"0"` or `"b"`; the pair of levels an effect measure
#'     compares, such as `"1 vs 0"` or `"b vs a"`; or the coefficient name on a
#'     surface whose rows are named after their coefficients. Present only in
#'     the marginal reading, and there only where a row reports one: the
#'     conditional reading names each coefficient in `term` and returns no
#'     `contrast` column.}
#'   \item{`estimate`}{The pooled point estimate.}
#'   \item{`std.error`}{The pooled standard error.}
#'   \item{`statistic`}{The test statistic, the estimate over its standard
#'     error.}
#'   \item{`df`}{The pooled degrees of freedom the statistic is referred to.
#'     These carry the Barnard-Rubin small-sample adjustment, so the statistic is
#'     referred to t on them rather than to the normal.}
#'   \item{`p.value`}{The two-sided p-value of that statistic.}
#'   \item{`conf.low`, `conf.high`}{The interval bounds. Present only when
#'     `conf.int` is `TRUE`.}
#' }
#'
#' @examplesIf requireNamespace("mice", quietly = TRUE)
#' set.seed(2024)
#' n <- 150
#' x1 <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.3 * x1))
#' y <- rbinom(n, 1, plogis(-0.4 + 0.9 * z + 0.5 * x1))
#' dat <- data.frame(x1, z, y)
#' dat$x1[rbinom(n, 1, 0.2) == 1] <- NA
#'
#' imp <- mice::mice(dat, m = 2, print = FALSE, seed = 1)
#' fits <- with(imp, {
#'   ps <- glm(z ~ x1, family = binomial())
#'   w <- wt_ate(ps)
#'   om <- glm(y ~ z, family = quasibinomial(), weights = w)
#'   ipw(ps, om)
#' })
#'
#' tidy(pool_ipw(fits))
#'
#' tidy(pool_ipw(fits), conf.int = TRUE, exponentiate = TRUE)
#'
#' # The outcome models' coefficients, pooled the same way
#' tidy(pool_ipw(fits), effects = "conditional")
#'
#' @seealso [pool_ipw()] for the pooling, [`glance()`][glance.ipw_pooled()] for a
#'   one-row summary of the pooled fit, [as_marginal()] and [as_conditional()]
#'   for the reading a pooled result records, and the Multiple imputation
#'   section of [ipw()] for the workflow they belong to.
#'
#' @exportS3Method generics::tidy ipw_pooled
tidy.ipw_pooled <- function(
  x,
  conf.int = FALSE,
  conf.level = NULL,
  exponentiate = FALSE,
  ...,
  effects = NULL
) {
  rlang::check_dots_empty()

  # The pooled result's own coercion surface, which validates the four arguments
  # and builds the table. Reading it rather than rebuilding it is what keeps
  # this method and the frame reporting the same numbers under the same names,
  # exactly as the marginal reading of `tidy()` on a single result does. It also
  # owns `effects`, resolving a `NULL` against the reading the pooled result
  # records and refusing both a value naming neither reading and a reading the
  # pooling had nothing to compute from.
  frame <- as.data.frame(
    x,
    conf.int = conf.int,
    conf.level = conf.level,
    exponentiate = exponentiate,
    effects = effects
  )

  # The covariance of the effects the frame attaches is the one thing that does
  # not travel: it is an attribute rather than a column, and a tidied table is
  # its columns.
  out <- tibble::as_tibble(frame)
  attr(out, "ipw_vcov") <- NULL

  out
}

#' Glance at a pooled inverse probability weighted result
#'
#' @description
#' `glance()` describes a [pool_ipw()] result rather than its estimates: one row
#' naming the estimand, counting the observations and the results that were
#' pooled, and reporting the complete-data degrees of freedom the small-sample
#' adjustment used. A pooled result reporting several effect measures, or several
#' contrasts of a categorical exposure, still returns exactly one row.
#'
#' @param x An `ipw_pooled` object, as returned by [pool_ipw()].
#' @param ... These dots are for future extensions and must be empty.
#'
#' @return A one-row [tibble][tibble::tibble] with the columns:
#' \describe{
#'   \item{`estimand`}{The causal estimand every pooled result targeted.}
#'   \item{`nobs`}{The smallest number of observations any pooled result was
#'     estimated from. The imputed datasets are the same size, so this is that
#'     size unless a result was fitted on a subset of one. Reported by
#'     [stats::nobs()].}
#'   \item{`m`}{The number of results pooled.}
#'   \item{`dfcom`}{The complete-data degrees of freedom the Barnard-Rubin
#'     adjustment used.}
#' }
#'
#' @examplesIf requireNamespace("mice", quietly = TRUE)
#' set.seed(2024)
#' n <- 150
#' x1 <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.3 * x1))
#' y <- rbinom(n, 1, plogis(-0.4 + 0.9 * z + 0.5 * x1))
#' dat <- data.frame(x1, z, y)
#' dat$x1[rbinom(n, 1, 0.2) == 1] <- NA
#'
#' imp <- mice::mice(dat, m = 2, print = FALSE, seed = 1)
#' fits <- with(imp, {
#'   ps <- glm(z ~ x1, family = binomial())
#'   w <- wt_ate(ps)
#'   om <- glm(y ~ z, family = quasibinomial(), weights = w)
#'   ipw(ps, om)
#' })
#'
#' glance(pool_ipw(fits))
#'
#' @seealso [pool_ipw()] for the pooling, [`tidy()`][tidy.ipw_pooled()] for the
#'   pooled estimates, and the Multiple imputation section of [ipw()] for the
#'   workflow the two belong to.
#'
#' @exportS3Method generics::glance ipw_pooled
glance.ipw_pooled <- function(x, ...) {
  rlang::check_dots_empty()

  # The count is the pooled result class' own accessor, so the number
  # `glance()` reports and the number `nobs()` reports are the same number
  # rather than two readings that agree by construction, exactly as the summary
  # of a single result reads its counts.
  tibble::tibble(
    estimand = x$estimand,
    nobs = stats::nobs(x),
    m = x$m,
    dfcom = x$dfcom
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
#' A column is an answer per observation only while the fit it is read from and
#' the frame it is carried on describe the same observations, so three things are
#' refused rather than reported. A result whose propensity score model and
#' outcome model were fit to different rows is refused, because the propensity
#' score of one observation would arrive beside the fitted value of another. Two
#' models of the same data most often part this way over missing values, each
#' dropping the rows a variable of its own is missing on. What is compared is the
#' number of observations each model produced an answer for and, when the outcome
#' model's frame names the exposure, the exposure values of the two model frames
#' position by position, reading a factor and the numbers its labels spell as one
#' encoding of one exposure. Exposure values that disagree prove the two models
#' hold different observations; exposure values that agree prove nothing, since
#' two different sets of rows are free to carry the same sequence of values, and
#' two encodings that cannot be read onto each other, such as a factor of `"a"`
#' and `"b"` against a recoding of it as 0 and 1, prove nothing either way.
#' The check is one-sided in that direction on purpose: it refuses only what it
#' can prove, so a fit whose models do line up is never refused over the labels
#' the rows of either frame carry or the encoding its exposure is written in. A
#' frame that already holds a column this call would add is refused as well,
#' because the fit's column has nowhere to go: added beside the frame's, the two
#' would be told apart by nothing, and written over it, the frame's data would be
#' gone.
#' Which names are in the way depends on the fit and the frame, `.resid` being
#' among them only for a frame it would be added to and the propensity columns
#' being the ones the propensity score model produced. Last, a result whose
#' outcome model was fit without weights is refused under the class [ipw()]
#' refuses the same fact with: `.weights` has nothing to report for such a fit.
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
  check_augment_weights(weights)

  supplied <- !is.null(data)
  if (supplied) {
    assert_class(data, "data.frame")
  } else {
    # The weights are reported as `.weights`, so the frame this method builds
    # for itself leaves out the numeric copy of them `model.frame()` records.
    # The drop belongs to this branch alone: a frame the caller supplies is
    # theirs, and a column of theirs is not deleted for matching a name.
    data <- outcome_frame[names(outcome_frame) != "(weights)"]
  }

  propensity <- augment_propensity_columns(x$wt_mod)
  fitted <- unname(stats::predict(x$outcome_mod, type = "response"))
  check_augment_alignment(x$wt_mod, propensity, fitted, outcome_frame)

  if (supplied) {
    check_augment_rows(data, nrow(outcome_frame))
  }

  # A residual is an observed outcome less a fitted value, so it belongs to a
  # frame that holds the outcome and to no other.
  response <- names(outcome_frame)[[
    attr(stats::terms(x$outcome_mod), "response")
  ]]
  resid <- response %in% names(data)

  # Which columns this call adds is a question about this fit and this frame:
  # the propensity columns are the ones the propensity score model produced, and
  # `.resid` is added only to a frame holding the outcome to difference.
  check_augment_columns(
    data,
    c(names(propensity), ".weights", ".fitted", if (resid) ".resid"),
    supplied = supplied
  )

  columns <- c(
    as.list(data),
    propensity,
    list(
      .weights = weights,
      .fitted = fitted
    )
  )

  # The outcome goes through the conversion the estimator itself uses, which is
  # what puts the residual of a factor or logical outcome on the 0/1 scale its
  # fitted values are on.
  if (resid) {
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
# level gets a column of its own. A column of a tibble is addressed by its
# position, so the observation names a prediction carries stop here.
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

# The outcome model of an ipw() result is fit with the estimated weights, so a
# result reporting one that was not is a result built around something else.
# `.weights` has nothing to report for such a fit, and the columns that remain
# would describe an analysis the object does not name. The fact is the one
# ipw() refuses the pair of models for, so it is reported under the class that
# names it rather than under one naming which method noticed.
check_augment_weights <- function(weights, call = rlang::caller_env()) {
  if (!is.null(weights)) {
    return(invisible(weights))
  }

  abort(
    c(
      "{.arg x} must hold an outcome model fit with propensity score weights.",
      x = "The outcome model {.arg x} holds reports no weights.",
      i = "An {.fun ipw} outcome model is weighted by construction. Refit it \\
      with the weights the propensity score model implies, such as those \\
      {.fun wt_ate} returns."
    ),
    error_class = "propensity_ipw_weights_missing_error",
    call = call
  )
}

# The propensity score column and the fitted values are read from two different
# models, so they are an answer per observation only while both models describe
# the same observations. Counts come from the fitted frames rather than from
# `stats::nobs()`, which subtracts observations of weight zero and so would
# report a fit that has any as shorter than the column it produces.
#
# What is checked beyond the counts is the exposure each model was fit to, and
# the check is one-sided by design: an exposure that disagrees between the two
# model frames proves they hold different observations, while an exposure that
# agrees proves nothing, two different sets of rows being free to carry the same
# sequence of exposure values. Refusing only what can be proven is what keeps a
# fit whose models do line up from being refused over how either frame labels its
# rows or writes its exposure down.
check_augment_alignment <- function(
  wt_mod,
  propensity,
  fitted,
  outcome_frame,
  call = rlang::caller_env()
) {
  n_outcome <- nrow(outcome_frame)
  n_fitted <- length(fitted)

  # Every propensity column holds one prediction per observation, so the first
  # of them measures them all.
  n_propensity <- length(propensity[[1L]])
  counted <- n_propensity == n_outcome && n_fitted == n_outcome

  exposure <- if (counted) {
    augment_mismatched_exposure(wt_mod, outcome_frame)
  }

  if (counted && is.null(exposure)) {
    return(invisible(propensity))
  }

  # Where the counts differ they say which model produced fewer answers; where
  # they agree it is the exposure that shows the two apart, the count saying
  # nothing on its own.
  problem <- if (counted) {
    c(
      x = "Both models were fit to {n_outcome} observation{?s}, and the \\
      {.val {exposure}} values in the two model frames disagree."
    )
  } else {
    c(
      x = "The propensity score model produced {n_propensity} propensity \\
      score{?s}.",
      x = "The outcome model was fit to {n_outcome} observation{?s}."
    )
  }

  # A model whose answers outnumber the rows it was fit to was fit under
  # `na.action = na.exclude`, which pads a prediction back to the length of the
  # data the model was given. That padding is the disagreement here, and the
  # counts above do not show it.
  if (n_fitted != n_outcome) {
    problem <- c(
      problem,
      x = "The outcome model predicts {n_fitted} value{?s} for the \\
      {n_outcome} observation{?s} it was fit to."
    )
  }

  abort(
    c(
      "{.arg x} must hold two models fit to the same rows.",
      problem,
      i = "Two models of the same data most often part this way over missing \\
      values, each dropping the rows a variable of its own is missing on.",
      i = "Refit both models on the rows that are complete on every variable \\
      either model uses."
    ),
    error_class = "propensity_augment_alignment_error",
    call = call
  )
}

# The exposure, where the two model frames can be shown to hold different
# observations of it, and `NULL` otherwise. Both models read the exposure from
# the same column of the same data, so a position at which the two frames record
# different exposure values is a position at which they describe different
# observations, whatever the row labels of either say. Row labels are not
# compared at all: they are a record of where a row came from, so a frame
# renumbered by the pipeline that built it carries different labels for the same
# observations, and two frames each renumbered from one carry the same labels for
# different ones.
#
# Nothing is refit and nothing is predicted: the propensity score model's own
# model frame holds the exposure it was fit to, after the rows it dropped, and
# the outcome model's frame holds the same variable whenever the outcome formula
# names it. A fit that offers neither, that offers two vectors of different
# lengths, or that records the exposure in two encodings no comparison can line
# up, offers nothing to compare, and the counts are then all there is.
augment_mismatched_exposure <- function(wt_mod, outcome_frame) {
  # A response written as anything but a bare variable, such as a transformation
  # or a two-column count, deparses to more than one name and is not a column of
  # the outcome model's frame under any name this could look up.
  name <- tryCatch(fmla_extract_left_chr(wt_mod), error = function(e) NULL)
  if (length(name) != 1L || !name %in% names(outcome_frame)) {
    return(NULL)
  }

  exposure <- tryCatch(fmla_extract_left_vctr(wt_mod), error = function(e) NULL)
  if (is.null(exposure)) {
    return(NULL)
  }

  from_wt_mod <- augment_exposure_values(exposure)
  from_outcome <- augment_exposure_values(outcome_frame[[name]])

  # Two vectors of different lengths cannot be compared position by position. A
  # propensity score model fit under `na.action = na.exclude` reaches this,
  # its predictions being padded back to the rows its frame dropped.
  if (length(from_wt_mod) != length(from_outcome)) {
    return(NULL)
  }

  compared <- augment_exposure_bridge(from_wt_mod, from_outcome)
  if (is.null(compared)) {
    return(NULL)
  }

  if (
    isTRUE(all.equal(compared[[1L]], compared[[2L]], check.attributes = FALSE))
  ) {
    return(NULL)
  }

  name
}

# The exposure read down to the values it records, so that one exposure held as
# a factor and as the labels behind it, or as integers and as the doubles of the
# same numbers, compares equal. A factor against the numbers those labels spell
# is a further step, taken by the bridge below. Two frames recording one
# exposure in two representations hold the same observations of it, and a check
# that refuses only what it can prove has to read them as such.
augment_exposure_values <- function(x) {
  if (is.factor(x)) {
    return(as.character(x))
  }

  if (is.logical(x)) {
    return(as.integer(x))
  }

  as.vector(x)
}

# One exposure recorded in two encodings, as a pair of vectors a position by
# position comparison can be made of, or `NULL` when the two encodings cannot be
# lined up at all. Two frames of one column usually record it the same way and
# there is nothing to bridge. Where they do not, the case to carry is a factor
# against the numbers it labels, which `ipw()` itself accepts: a propensity score
# model fit to `factor(z)` and an outcome model fit to the 0/1 `z` of the same
# data were fit to one exposure, and reading the labels as the numbers they spell
# is what lets the two meet.
#
# Labels spelling something other than numbers are the other side of that case:
# an exposure held as "a" and "b" against a recoding of it as 0 and 1 is one
# variable written two ways that no positional comparison can line up, its values
# differing at every position while its observations may agree at all of them.
# Refusing there would refuse a fit whose models describe the same rows, so the
# comparison does not apply and the counts stand alone.
augment_exposure_bridge <- function(x, y) {
  if (is.character(x) == is.character(y)) {
    return(list(x, y))
  }

  if (is.character(x)) {
    x <- augment_exposure_numbers(x)
  } else {
    y <- augment_exposure_numbers(y)
  }

  if (is.null(x) || is.null(y)) {
    return(NULL)
  }

  list(x, y)
}

# Labels read as the numbers they spell, or `NULL` when any of them spells
# something else. A label that is not a number reads as `NA`, a value the
# exposure never held, which would then compare unequal to whatever the other
# frame records in that position and prove a disagreement that is the encoding's
# and not the data's.
augment_exposure_numbers <- function(x) {
  numbers <- suppressWarnings(as.numeric(x))

  if (any(is.na(numbers) & !is.na(x))) {
    return(NULL)
  }

  numbers
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

# The dot-prefixed columns are additions, so a frame already holding one of
# their names has nowhere for the fit's own column to go: added beside the
# frame's column, the two would be told apart by nothing, and written over it,
# the data the frame arrived with would be gone.
check_augment_columns <- function(
  data,
  added,
  supplied,
  call = rlang::caller_env()
) {
  clashing <- intersect(added, names(data))
  if (length(clashing) == 0) {
    return(invisible(data))
  }

  # The remedy is the frame's own: a frame the caller supplied is theirs to
  # rename, and the frame this method builds for itself is the outcome model's,
  # whose columns are the variables that model was fit on.
  remedy <- if (supplied) {
    "Rename or drop {cli::qty(clashing)}{?it/them} in {.arg data}, or leave \\
    {.arg data} as {.code NULL} to carry the fit's columns on the outcome \\
    model's own frame."
  } else {
    "The frame is the outcome model's own, so {cli::qty(clashing)}{?that name \\
    is a variable/those names are variables} its formula reads. Rename \\
    {?it/them} in the data the model was fit to, or supply a frame through \\
    {.arg data}."
  }

  abort(
    c(
      "{.fun augment} adds {cli::qty(clashing)}{?a column/columns} the frame \\
      carrying {?it/them} already holds.",
      x = "The frame already holds {.val {clashing}}.",
      i = remedy
    ),
    error_class = "propensity_augment_column_error",
    call = call
  )
}
