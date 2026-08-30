#' Truncate (Winsorize) Propensity Scores
#'
#' `ps_trunc()` bounds extreme propensity scores to fixed limits, replacing
#' out-of-range values with the boundary value (a form of *winsorizing*). The
#' result is a vector or matrix of the same length and dimensions as
#' `.propensity`, with no observations removed. This contrasts with [ps_trim()],
#' which sets extreme values to `NA` (effectively removing those observations
#' from analysis).
#'
#' @param .propensity A numeric vector of propensity scores in (0, 1) (binary
#'   exposures), or a matrix/data.frame where each column contains propensity
#'   scores for one level of a categorical exposure. A data frame truncated for
#'   a binary exposure is reduced to a single column: the second column of a two
#'   column data frame, which is the probability of the second level in the
#'   layout model predictions come in, and the first column otherwise. The
#'   column taken is announced; `options(propensity.quiet = TRUE)` silences the
#'   announcement. A vector is held to the open interval, so a score of exactly
#'   0 or 1 is refused. A matrix is held to the closed interval instead: a cell
#'   of exactly 0 or 1 is bounded rather than refused, which is how a separated
#'   multinomial fit is repaired, and what comes back lies strictly inside
#'   (0, 1). Anything outside \eqn{[0, 1]} is refused either way, and so is a
#'   truncation whose own computed bounds leave a score at an endpoint; see
#'   **Propensity scores at 0 and 1** in [wt_ate()].
#' @param .exposure An exposure vector. Required for method `"cr"` (binary
#'   exposure vector) and for categorical exposures (factor or character vector)
#'   with any method.
#' @param method One of `"ps"`, `"pctl"`, or `"cr"`:
#'   * `"ps"` (default): Truncate directly on propensity score values. Values
#'     outside `[lower, upper]` are set to the nearest bound. For categorical
#'     exposures, applies symmetric truncation using `lower` as the threshold
#'     (delta) and renormalizes rows to sum to 1.
#'   * `"pctl"`: Truncate at quantiles of the propensity score distribution.
#'     The `lower` and `upper` arguments specify quantile probabilities. For
#'     categorical exposures, quantiles are computed across all columns.
#'   * `"cr"`: Truncate to the common range of propensity scores across
#'     exposure groups (binary exposures only). Bounds are
#'     `[min(.propensity[focal]), max(.propensity[reference])]`. Requires
#'     `.exposure`. When those bounds cross, so that the two distributions do
#'     not overlap at all, there is no common range to bound the scores to and
#'     the call errors with class `propensity_no_overlap_error`. That differs
#'     deliberately from [ps_trim()], which trims every observed unit in the
#'     same situation, a truthful record of an empty overlap region.
#'
#'   For categorical exposures, only `"ps"` and `"pctl"` are supported.
#' @param lower,upper Bounds for truncation. Interpretation depends on `method`:
#'   * `method = "ps"`: Propensity score values (defaults: 0.1 and 0.9). For
#'     categorical exposures, `lower` is the truncation threshold delta
#'     (default: 0.01) and `upper` is ignored. That default deliberately differs
#'     from the 0.1 threshold [ps_trim()] uses for categorical exposures:
#'     truncation keeps every unit and only pins the most extreme scores back to
#'     the threshold, so its default is a gentle winsorization, whereas trimming
#'     discards the units it selects and follows common-support trimming
#'     practice. With `k` exposure levels, a threshold of `1/k` or larger cannot
#'     be met by every column of a row that sums to one, and is an error.
#'   * `method = "pctl"`: Quantile probabilities (defaults: 0.05 and 0.95;
#'     categorical defaults: 0.01 and 0.99).
#'   * `method = "cr"`: Ignored; bounds are determined by the data.
#' @inheritParams wt_ate
#' @param .focal_level The value of `.exposure` representing the focal
#'   (treated) group, used by `"cr"`. Every binary coding honors it: 0/1
#'   numeric, logical, two-level factor, and two-level character exposures are
#'   all coded with the named level as focal, and a level the exposure never
#'   takes is an error. With no level named, a binary exposure defaults to its
#'   higher level, which is `1` for a 0/1 exposure, `TRUE` for a logical one,
#'   and the second of the two levels a factor or character exposure takes.
#'   Levels a factor declares but never takes are not candidates. Naming any
#'   other level reverses the coding, so `.propensity` must then hold the
#'   probability of the named level.
#' @param .reference_level The value of `.exposure` representing the reference
#'   (control) group. Naming it makes the exposure's other level focal, with
#'   the same consequence for `.propensity`, and a level the exposure never
#'   takes is an error. Automatically detected if not supplied.
#' @param ... Additional arguments passed to methods.
#' @param ps `r lifecycle::badge("deprecated")` Use `.propensity` instead. A
#'   call that names `ps` must name the arguments after it as well, since a
#'   positional argument binds to `.propensity`.
#'
#' @details
#' Unlike [ps_trim()], truncation preserves all observations. No `NA` values
#' are introduced; out-of-range scores are replaced with the boundary value.
#'
#' For **binary exposures**, each propensity score \eqn{e_i} is bounded:
#'
#'  - If \eqn{e_i < l}, set \eqn{e_i = l} (the lower bound).
#'
#'  - If \eqn{e_i > u}, set \eqn{e_i = u} (the upper bound).
#'
#' For **categorical exposures**, values below the threshold are set to the
#' threshold and each row is renormalized to sum to 1.
#'
#' **Arithmetic behavior**: Arithmetic operations on `ps_trunc` objects return
#' plain numeric vectors. Once propensity scores are transformed (e.g., into
#' weights), the result is no longer a propensity score.
#'
#' **Combining behavior**: Combining `ps_trunc` objects with `c()` requires
#' matching truncation parameters. Mismatched parameters produce a warning and
#' return a plain numeric vector.
#'
#' ## Fitted models
#'
#' `.propensity` can be the fitted propensity score model instead of the scores
#' it reports. A binomial [stats::glm()] and a two-level `nnet::multinom()` are
#' read as one score per unit; a `nnet::multinom()` of three or more levels is
#' read as one column per level. Those are the shapes `predict(fit, type =
#' "response")` and `fitted(fit)` give, and truncating a fit bounds exactly what
#' bounding those values would.
#'
#' The methods that read an exposure (`"cr"`, and every method on the
#' categorical route) take it from the model when `.exposure` is not supplied,
#' announcing the variable they read; `options(propensity.quiet = TRUE)`
#' silences the announcement. An `.exposure` you supply is used instead, and a
#' categorical model's columns are matched to its levels by name, so an exposure
#' whose levels are ordered differently is still bounded against the right
#' column.
#'
#' ## Missing values
#'
#' A propensity score that arrives missing is not one this function pinned to a
#' bound, so it takes no part in the truncation record: it never joins
#' `truncated_idx`, and [is_unit_truncated()] reports `FALSE` for it. The value
#' stays missing in the result. For a matrix of categorical propensity scores, a
#' row with a missing cell comes back missing throughout, because renormalizing
#' a row divides by its sum and that sum is unknown.
#'
#' The methods that read their bounds off the data read them off the data they
#' have. `"pctl"` and `"cr"` work their bounds out from the complete scores, and
#' the categorical `"pctl"` bounds come from the cells of the complete rows, so
#' the bounds are the ones the same call would produce with the missing
#' observations dropped.
#'
#' A missing exposure is a different matter, and `"cr"` refuses one with an
#' error of class `propensity_missing_value_error`. Its bounds come from the
#' exposure groups, and a unit that belongs to neither leaves them undefined.
#' Remove or impute the missing exposure values first.
#'
#' ## The truncation record
#'
#' A `ps_trunc` records which units were winsorized as positions among the
#' observations it was written for, along with how many observations that was.
#' Operations that hand this package the subscript re-index those positions onto
#' the result: subsetting with `[`, [sort()], [unique()], and [rep()] all return
#' a record written for what they return, and a subscript naming a position more
#' than once reports that unit at every place it now holds.
#'
#' Operations that change how many observations there are without supplying a
#' subscript cannot re-index the record, and it is dropped rather than worked
#' out from the values, since a score that arrived equal to a bound is
#' indistinguishable from one this function pinned there. [vctrs::vec_slice()],
#' which is how filtering, joining, and grouped summaries in dplyr reach a
#' column, is the usual route, and dropping the record there raises a warning
#' of class `propensity_trunc_record_warning`. Combining with [c()] drops it
#' without comment, because concatenation appends one set of observations to
#' another and the prototype it builds the result from holds no positions to
#' lose. The values, the class, and the method and its bounds are untouched
#' either way.
#'
#' A record can also outlive the observations it describes, because it travels
#' by routes vctrs does not see: growing a `ps_trunc` by subassignment carries
#' it across the length change. [is_unit_truncated()] therefore checks that the
#' record covers the object it is given and raises an error of class
#' `propensity_missing_meta_error` when it does not, rather than name truncated
#' units at stale positions.
#'
#' That check compares how many observations the record was written for against
#' how many the object holds, which a reordering does not change. An operation
#' that reorders the observations through vctrs, rather than through `[`,
#' therefore keeps a record written for the order they used to be in:
#' `vctrs::vec_slice(x, 5:1)` and `dplyr::arrange()` both return the values in
#' a new order under positions still naming the old one, and
#' [is_unit_truncated()] answers from those positions and names the wrong
#' units. Subsetting with `[` is handed the subscript and re-indexes, so
#' reorder with `[`, or put the propensity scores in the order you want before
#' truncating them.
#'
#' Casting a numeric vector into a `ps_trunc` with [vctrs::vec_cast()] is a type
#' operation and not a truncation. The result is described by the method and
#' bounds of the target and records that none of the arriving values was pinned
#' to a bound, so it can hold scores outside those bounds, including 0 and 1.
#' Call `ps_trunc()` on the scores to truncate them.
#'
#' @return A `ps_trunc` object (a numeric vector for binary exposures, or a
#'   matrix for categorical exposures). Use [ps_trunc_meta()] to inspect
#'   metadata including `method`, `lower_bound`, `upper_bound`,
#'   `truncated_idx` (positions of modified values), and `n_obs` (the number of
#'   observations those positions describe).
#'
#' @references
#' Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing
#' with limited overlap in estimation of average treatment effects.
#' *Biometrika*, 96(1), 187--199.
#'
#' Walker, A. M., Patrick, A. R., Lauer, M. S., et al. (2013). A tool for
#' assessing the feasibility of comparative effectiveness research.
#' *Comparative Effectiveness Research*, 3, 11--20.
#'
#' @seealso [ps_trim()] for removing (rather than bounding) extreme values,
#'   [ps_refit()] for refitting the propensity model after trimming,
#'   [is_ps_truncated()], [is_unit_truncated()], [ps_trunc_meta()]
#'
#' @examples
#' set.seed(2)
#' n <- 200
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(1.2 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' # Truncate to [0.1, 0.9]
#' ps_t <- ps_trunc(ps, method = "ps", lower = 0.1, upper = 0.9)
#' ps_t
#'
#' # Truncate at the 1st and 99th percentiles
#' ps_trunc(ps, method = "pctl", lower = 0.01, upper = 0.99)
#'
#' # Use truncated scores to calculate weights
#' wt_ate(ps_t, .exposure = z)
#'
#' # Inspect which observations were truncated
#' is_unit_truncated(ps_t)
#'
#' # Bound the scores a fitted model reports, reading the exposure off the model
#' ps_trunc(fit, method = "cr")
#'
#' if (rlang::is_installed("nnet")) {
#'   trt <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
#'   multinomial_fit <- nnet::multinom(trt ~ x, trace = FALSE)
#'   ps_trunc(multinomial_fit, method = "ps")
#' }
#'
#' @export
ps_trunc <- function(
  .propensity,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated()
) {
  .propensity <- handle_propensity_deprecation(
    rlang::maybe_missing(.propensity),
    ps,
    "ps_trunc"
  )

  UseMethod("ps_trunc", .propensity)
}

#' @export
ps_trunc.default <- function(
  .propensity,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  # Two frames arrive here because two condition systems read them.
  # `user_env` is the frame lifecycle reports a deprecation from, which decides
  # whether the reader is told to change their own call or to report an issue.
  # `call` is the frame every other condition is attributed to, which rlang
  # reads to name the function in the report. A route that reaches this method
  # by a call rather than by dispatch has to supply both.
  user_env = rlang::caller_env(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  # The dots reach every route but are read by none of them, so a misspelled
  # bound such as `lowr` would bound the scores at the default the caller
  # believed they had replaced. The model methods forward their dots here, so
  # the refusal covers a fit as well as a vector.
  rlang::check_dots_empty(call = call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)
  check_ps_method(.propensity, call = call)

  # A calibrated score is a propensity score, so this reads the scores it holds
  # the way the weight functions do rather than refusing a vector of numbers
  # for the class attached to it. What comes back is a bounded score rather
  # than a calibrated one, so the class is dropped and what it recorded is kept
  # in the truncation metadata, where `is_ps_calibrated()` reads it.
  calibrated <- inherits(.propensity, "ps_calib")
  if (calibrated) {
    .propensity <- as.numeric(.propensity)
  }

  method <- rlang::arg_match(method, error_call = call)
  meta_list <- list(method = method)

  if (calibrated) {
    meta_list$calibrated <- TRUE
  }

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "ps_trunc",
    user_env = user_env
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level, call = call)

  check_ps_range(.propensity, call = call)

  if (method == "ps") {
    if (is.null(lower)) {
      lower <- 0.1
    }
    if (is.null(upper)) {
      upper <- 0.9
    }
    check_lower_upper(lower, upper, call = call)

    lb <- lower
    ub <- upper
  } else if (method == "pctl") {
    if (is.null(lower)) {
      lower <- 0.05
    }
    if (is.null(upper)) {
      upper <- 0.95
    }
    check_quantile_probs(lower, upper, call = call)
    check_lower_upper(lower, upper, call = call)
    # `quantile()` names its result for the probability it was asked for, which
    # says nothing about the bound and reappears wherever the bound is printed
    # or compared.
    lb <- unname(quantile(.propensity, probs = lower, na.rm = TRUE))
    ub <- unname(quantile(.propensity, probs = upper, na.rm = TRUE))
    meta_list$lower_pctl <- lower
    meta_list$upper_pctl <- upper
  } else if (method == "cr") {
    if (is.null(.exposure)) {
      abort(
        "For {.code method = 'cr'}, must supply {.arg .exposure}.",
        error_class = "propensity_missing_arg_error",
        call = call
      )
    }
    check_exposure_complete(.exposure, method, call = call)
    .exposure <- transform_exposure_binary(
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level,
      call = call
    )
    # A score that arrived missing is not one this function can place against a
    # bound, so it takes no part in working the bounds out. Dropping it here
    # gives the bounds the same call would produce with the observation absent.
    observed <- !is.na(.propensity)
    ps_treat <- .propensity[observed & .exposure == 1]
    ps_untrt <- .propensity[observed & .exposure == 0]
    check_cr_groups_observed(ps_treat, ps_untrt, call = call)
    lb <- min(ps_treat)
    ub <- max(ps_untrt)
    check_common_range(lb, ub, call = call)
  }

  # Winsorize. A missing score compares as neither below the lower bound nor
  # above the upper one, so `which()` leaves it out of the pinned positions and
  # it stays missing in the result: it is not a score this function bounded.
  pinned_low <- which(.propensity < lb)
  pinned_high <- which(.propensity > ub)
  truncated_idx <- sort(c(pinned_low, pinned_high))

  .propensity[pinned_low] <- lb
  .propensity[pinned_high] <- ub

  meta <- c(
    meta_list,
    list(
      lower_bound = lb,
      upper_bound = ub,
      truncated_idx = truncated_idx,
      n_obs = length(.propensity)
    )
  )

  new_ps_trunc(.propensity, meta)
}

#' @export
ps_trunc.matrix <- function(
  .propensity,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  rlang::check_dots_empty(call = call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)

  # The generic offers every method, so match against that full set and then
  # reject the ones the categorical path does not define.
  method <- rlang::arg_match(
    method,
    values = c("ps", "pctl", "cr"),
    error_call = call
  )
  if (!method %in% c("ps", "pctl")) {
    abort(
      c(
        "Method {.val {method}} is not supported for categorical exposures.",
        i = "Use {.val ps} or {.val pctl}."
      ),
      error_class = "propensity_method_error",
      call = call
    )
  }

  check_no_focal_levels(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    call = call
  )

  # Validate exposure for categorical
  if (is.null(.exposure)) {
    abort(
      "`.exposure` must be provided for categorical propensity score truncation.",
      error_class = "propensity_missing_arg_error",
      call = call
    )
  }

  # Transform to factor and validate
  .exposure <- transform_exposure_categorical(.exposure, call = call)

  # Validate matrix. Truncation exists to pull a score at an endpoint back
  # inside the interval, which is the repair a separated multinomial fit needs,
  # so this route reads the closed interval: a cell of exactly 0 or 1 is bounded
  # rather than refused. Everything outside [0, 1] is still refused, and so is
  # any endpoint left standing after the bounds are applied.
  .propensity <- check_ps_matrix(
    .propensity,
    .exposure,
    closed = TRUE,
    call = call
  )

  n <- nrow(.propensity)
  k <- ncol(.propensity)

  # Renormalizing a row divides by its sum, which is unknown when one of its
  # scores is missing, so a row with a missing cell comes back missing whatever
  # the bounds are. It is not a row this function pinned to a bound, and it takes
  # no part in the bounds either. Both rules below decide membership with a test
  # read across the row, which answers `TRUE` off an observed cell that happens
  # to be out of bounds however the missing one is treated, so the incomplete
  # rows are taken out of the record explicitly.
  complete_rows <- !apply(.propensity, 1, anyNA)

  if (method == "ps") {
    # Symmetric truncation
    if (is.null(lower)) {
      lower <- 0.01
    } # Default threshold
    delta <- lower # Use lower as delta for consistency

    # A threshold at or above 1/k cannot be met by every column of a row that
    # sums to one, so there is no truncation rule left to apply. Both numbers
    # are facts about this call: the threshold supplied and the limit the width
    # of the matrix imposes on it.
    if (delta >= 1 / k) {
      limit <- format(1 / k, digits = 7)
      abort(
        c(
          "The truncation threshold must fall below 1/k, for k columns of
           propensity scores.",
          x = "{.arg lower} is {.val {delta}}, and 1/k is {limit} for the {k}
               column{?s} the scores hold.",
          i = "No row summing to one can hold every score above 1/k, so a
               threshold there leaves no rule to apply."
        ),
        error_class = "propensity_range_error",
        call = call
      )
    }

    # Track which values were truncated
    truncated_idx <- which(
      complete_rows & apply(.propensity, 1, function(x) any(x < delta))
    )

    # Apply truncation and renormalize
    bounded <- .propensity
    for (i in 1:n) {
      row_vals <- bounded[i, ]
      # Clamp values below delta
      row_vals[row_vals < delta] <- delta
      # Renormalize to sum to 1
      bounded[i, ] <- row_vals / sum(row_vals)
    }

    lower_bound <- delta
    upper_bound <- NA_real_
  } else {
    # pctl
    # Percentile-based truncation
    if (is.null(lower)) {
      lower <- 0.01
    } # Default percentile
    if (is.null(upper)) {
      upper <- 0.99
    } # Default percentile
    check_quantile_probs(lower, upper, call = call)
    check_lower_upper(lower, upper, call = call)

    # Calculate thresholds based on the distribution of all propensity scores,
    # taken from the rows this function can renormalize. `quantile()` names its
    # result for the probability it was asked for, which says nothing about the
    # bound and reappears wherever the bound is printed or compared.
    all_ps_vals <- as.vector(.propensity[complete_rows, , drop = FALSE])
    lower_threshold <- unname(quantile(all_ps_vals, probs = lower))
    upper_threshold <- unname(quantile(all_ps_vals, probs = upper))

    # Track which rows had values truncated
    truncated_idx <- which(
      complete_rows &
        apply(
          .propensity,
          1,
          function(x) any(x < lower_threshold | x > upper_threshold)
        )
    )

    # Apply truncation and renormalize
    bounded <- .propensity
    for (i in 1:n) {
      row_vals <- bounded[i, ]
      # Clamp values outside thresholds
      row_vals[row_vals < lower_threshold] <- lower_threshold
      row_vals[row_vals > upper_threshold] <- upper_threshold
      # Renormalize to sum to 1
      bounded[i, ] <- row_vals / sum(row_vals)
    }

    lower_bound <- lower_threshold
    upper_bound <- upper_threshold
  }

  check_truncated_matrix_interior(
    bounded,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    call = call
  )

  meta <- list(
    method = method,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    truncated_idx = truncated_idx,
    n_obs = n,
    is_matrix = TRUE
  )

  # Add percentile info if using pctl method
  if (method == "pctl") {
    meta$lower_pctl <- lower
    meta$upper_pctl <- upper
  }

  new_ps_trunc(bounded, meta)
}

# The categorical matrix route reads the closed interval because bounding a
# score of exactly 0 or 1 is the repair it exists to make. What comes back has
# to be a propensity score matrix the rest of the package will read, so the
# result is held to the open interval the input was excused from.
#
# The case this catches is a `"pctl"` truncation of a matrix with enough zeros
# in it that the lower quantile is itself zero: the bound is then the endpoint,
# and the truncation pins those cells to where they already sat. The bound is
# reported because it is what the caller has to replace.
check_truncated_matrix_interior <- function(
  bounded,
  lower_bound,
  upper_bound,
  call = rlang::caller_env()
) {
  at_endpoint <- !is.na(bounded) & (bounded <= 0 | bounded >= 1)
  if (!any(at_endpoint)) {
    return(invisible(TRUE))
  }

  n_left <- sum(at_endpoint)

  # Written out here rather than in the bullet: a threshold and a pair of
  # bounds are two different sentences, and neither is a plural of the other.
  bound_text <- if (is.na(upper_bound)) {
    cli::format_inline(
      "The threshold computed from the scores is
       {.val {as.numeric(lower_bound)}}."
    )
  } else {
    cli::format_inline(
      "The bounds computed from the scores run from
       {.val {as.numeric(lower_bound)}} to {.val {as.numeric(upper_bound)}}."
    )
  }

  abort(
    c(
      "Truncation left {n_left} propensity score{?s} at 0 or 1.",
      x = bound_text,
      i = "A bound at an endpoint pins a score to where it already sat.",
      i = "Supply {.arg lower} and {.arg upper} explicitly, inside the open
           interval, so that every score is bounded away from 0 and 1."
    ),
    error_class = "propensity_range_error",
    call = call
  )
}

#' @export
ps_trunc.data.frame <- function(
  .propensity,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)

  # For categorical exposures, convert to matrix and call matrix method
  if (!is.null(.exposure)) {
    exposure_type <- causalgenerics::detect_exposure_type(
      .exposure,
      announce = !be_quiet(),
      call = call
    )
    if (exposure_type == "categorical") {
      ps_matrix <- as.matrix(.propensity)
      # The focal arguments travel on so that the matrix method refuses them.
      # Dropping them here would leave the caller believing the truncation
      # honored a coding the categorical path never reads.
      #
      # The matrix method is reached here by a call rather than by dispatch, so
      # it is handed a frame to report against. Left to its own, a refusal on
      # this route would name the method the caller never wrote.
      return(ps_trunc.matrix(
        .propensity = ps_matrix,
        method = method,
        lower = lower,
        upper = upper,
        .exposure = .exposure,
        .focal_level = .focal_level,
        .reference_level = .reference_level,
        ...,
        .treated = .treated,
        .untreated = .untreated,
        call = call
      ))
    }
  }

  ps_vec <- binary_ps_column(.propensity, "ps_trunc")

  # The default method is reached here by a call rather than by dispatch, so it
  # is handed a frame to report against. Left to its own, a refusal on this
  # route would name the method the caller never wrote.
  ps_trunc.default(
    .propensity = ps_vec,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .treated = .treated,
    .untreated = .untreated,
    user_env = rlang::caller_env(),
    call = call
  )
}

# The fitted propensity score models truncation reads, registered for the same
# classes the weight functions read: a `glm`, whose binomial families fit the
# probability of a binary exposure, and a `multinom`, which fits a probability
# for every level. A `lm` is not among them, its fitted values being conditional
# means rather than probabilities, and it reaches the default method, which
# reports that it has no scores to bound.
#' @export
ps_trunc.glm <- function(
  .propensity,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)

  ps_trunc_from_model(
    .propensity,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .treated = .treated,
    .untreated = .untreated,
    call = call,
    user_env = rlang::caller_env()
  )
}

# A `multinom` is `c("multinom", "nnet")` and inherits from neither `glm` nor
# `lm`, so it reaches a method of its own. What it was fit to is checked before
# anything reads it, so that a fit the route cannot read is reported as such
# rather than as whatever its fitted values are made of.
#' @export
ps_trunc.multinom <- function(
  .propensity,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated(),
  call = rlang::current_env()
) {
  check_call_arg(call)
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)
  check_multinom_response(.propensity, call = call)

  ps_trunc_from_model(
    .propensity,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .treated = .treated,
    .untreated = .untreated,
    call = call,
    user_env = rlang::caller_env()
  )
}

# The model methods differ only in the class they are registered for and in what
# that class needs checked before it is read, so both route through here. The
# method is resolved first because it decides whether the truncation reads an
# exposure at all.
#
# The vector and matrix methods are reached by a call rather than by dispatch,
# so they are handed the frame the route was entered on. Left to their own, a
# refusal here would name a method the caller never wrote.
ps_trunc_from_model <- function(
  model,
  method,
  lower,
  upper,
  .exposure,
  .focal_level,
  .reference_level,
  ...,
  .treated,
  .untreated,
  call,
  user_env
) {
  method <- rlang::arg_match(
    method,
    values = c("ps", "pctl", "cr"),
    error_call = call
  )

  args <- prepare_model_ps(
    model,
    .exposure,
    # Truncating a matrix of scores reads the exposure whatever the method, to
    # hold the columns against its levels, so a fit that reports a column for
    # every level is read for its exposure too.
    needs_exposure = method == "cr" || model_fits_levels(model),
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    .treated = .treated,
    .untreated = .untreated,
    fn_name = "ps_trunc",
    call = call,
    user_env = user_env
  )

  if (identical(args$exposure_type, "categorical")) {
    return(ps_trunc.matrix(
      .propensity = args$propensity,
      method = method,
      lower = lower,
      upper = upper,
      .exposure = args$exposure,
      .focal_level = args$focal_level,
      .reference_level = args$reference_level,
      ...,
      call = call
    ))
  }

  ps_trunc.default(
    .propensity = args$propensity,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = args$exposure,
    .focal_level = args$focal_level,
    .reference_level = args$reference_level,
    ...,
    user_env = user_env,
    call = call
  )
}

#' @export
ps_trunc.ps_trunc <- function(
  .propensity,
  method = c("ps", "pctl", "cr"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL,
  ps = lifecycle::deprecated()
) {
  .propensity <- read_method_propensity(rlang::maybe_missing(.propensity), ps)

  warn(
    "Propensity scores have already been truncated. Returning original object.",
    warning_class = "propensity_already_modified_warning"
  )
  .propensity
}

# The common range is the region both exposure groups reach, which is empty when
# the lowest score among the focal units sits above the highest among the
# reference ones. Winsorizing to a crossed pair of bounds pins every score to the
# bound on the far side of the other, and records a common range the groups do
# not have.
check_common_range <- function(lb, ub, call = rlang::caller_env()) {
  if (lb <= ub) {
    return(invisible(TRUE))
  }

  # The bounds are propensity scores read off the data, so they arrive at full
  # double precision and say nothing a reader can use past the first few digits.
  lb_shown <- signif(lb, 3)
  ub_shown <- signif(ub, 3)

  abort(
    c(
      "The exposure groups' propensity score distributions do not overlap.",
      x = "The lowest score among the focal units is {.val {lb_shown}}, above
           the highest among the reference units, {.val {ub_shown}}.",
      i = "{.code method = 'cr'} bounds the scores to the region both groups
           reach, and these groups reach none in common.",
      i = "Refit the propensity score model, or truncate with
           {.code method = 'ps'} or {.code method = 'pctl'}."
    ),
    call = call,
    error_class = "propensity_no_overlap_error"
  )
}

new_ps_trunc <- function(x, meta) {
  if (is.matrix(x)) {
    # For matrices, we don't use vctrs
    structure(
      x,
      ps_trunc_meta = meta,
      class = c("ps_trunc_matrix", "ps_trunc", "matrix")
    )
  } else {
    vec_assert(x, ptype = double())
    new_vctr(
      x,
      ps_trunc_meta = meta,
      class = "ps_trunc",
      inherit_base_type = TRUE
    )
  }
}

# The positional half of a truncation record. The rest of it, the method and its
# bounds, describes the truncation rather than the units and means the same
# thing at any length.
reindex_trunc_record <- function(meta, i) {
  meta$truncated_idx <- reindex_positions(meta$truncated_idx, i)
  meta$n_obs <- length(i)

  meta
}

drop_trunc_record <- function(meta) {
  meta[["truncated_idx"]] <- NULL
  meta[["n_obs"]] <- NULL

  meta
}

# A record that no longer describes the observations in front of it is dropped
# rather than guessed at. Nothing in the values says which units a shorter or a
# longer vector once had winsorized: a score that arrived equal to a bound is
# indistinguishable from one this package pinned there.
#
# A record over no observations names no unit, so replacing it costs the caller
# nothing and goes without comment. That is the record every prototype carries,
# which is what concatenation builds its result from.
discard_trunc_record <- function(meta, n) {
  recorded <- meta$n_obs

  if (!is.null(recorded) && recorded > 0) {
    warn(
      c(
        "Dropping the record of which units were truncated.",
        i = "The record describes {recorded} observation{?s} and this result
             has {n}, so its positions do not describe them.",
        i = "The values are unchanged and the result is still a
             {.cls ps_trunc}. Truncate the propensity scores you want to work
             with to get a record written for them."
      ),
      warning_class = "propensity_trunc_record_warning",
      # One of the routes here is vctrs' internal dispatch, whose call would be
      # reported and names nothing the caller wrote, so no call is attributed.
      call = NULL
    )
  }

  drop_trunc_record(meta)
}

# `i` holds the old positions the result is built from, in the order it holds
# them, among the `n_obs` observations it was taken from. A record that covers
# those observations re-indexes onto `i`; one that does not cannot be placed
# against them at all.
carry_trunc_record <- function(meta, n_obs, i) {
  if (record_covers(meta, n_obs)) {
    reindex_trunc_record(meta, i)
  } else {
    discard_trunc_record(meta, length(i))
  }
}

# A positional query reads its answer out of the record, so a record that does
# not cover the object in front of it has no answer to give: reporting every
# unit as untouched would be a wrong answer rather than a missing one.
check_trunc_record <- function(meta, n, fn, call = rlang::caller_env()) {
  if (record_covers(meta, n)) {
    return(invisible(meta))
  }

  recorded <- meta$n_obs
  problem <- if (is.null(recorded)) {
    "These propensity scores carry no record of which units were truncated."
  } else {
    "The record describes {recorded} observation{?s} and these propensity
     scores have {n}, so its positions do not describe them."
  }

  abort(
    c(
      "{.code {fn}} has no usable truncation record for these propensity
       scores.",
      x = problem,
      i = "An operation that changed the number of observations dropped it.
           Truncate the propensity scores you want to work with, or query the
           object the record was written for."
    ),
    error_class = "propensity_missing_meta_error",
    call = call
  )
}

#' Extract truncation metadata from a `ps_trunc` object
#'
#' @description Returns the metadata list attached to a [`ps_trunc`][ps_trunc()]
#'   object. The list includes fields such as `method`, `lower_bound`,
#'   `upper_bound`, and `truncated_idx`.
#'
#' @param x A `ps_trunc` object created by [ps_trunc()].
#' @return A named list with truncation metadata, including:
#'   * `method` -- the truncation method used (`"ps"`, `"pctl"`, or `"cr"`)
#'   * `lower_bound`, `upper_bound` -- the applied bounds
#'   * `truncated_idx` -- integer positions of values that were winsorized
#'   * `n_obs` -- the number of observations those positions describe
#'
#'   `truncated_idx` and `n_obs` are absent when an operation that changed the
#'   number of observations dropped the record; see [ps_trunc()].
#'
#' @seealso [ps_trunc()], [is_ps_truncated()], [is_unit_truncated()]
#'
#' @examples
#' ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
#' ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
#' ps_trunc_meta(ps_t)
#'
#' @export
ps_trunc_meta <- function(x) {
  attr(x, "ps_trunc_meta")
}

#' Test whether propensity scores have been truncated
#'
#' @description `is_ps_truncated()` returns `TRUE` if `x` is a `ps_trunc`
#'   object or a `psw` object derived from truncated propensity scores.
#'   Use [is_unit_truncated()] to find out *which* observations were modified.
#'
#' @param x An object.
#' @return A single `TRUE` or `FALSE`.
#'
#' @seealso [ps_trunc()], [is_unit_truncated()], [ps_trunc_meta()]
#'
#' @examples
#' ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
#' is_ps_truncated(ps)
#'
#' ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
#' is_ps_truncated(ps_t)
#'
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

#' @export
is_ps_truncated.ps_trunc_matrix <- function(x) {
  TRUE
}

#' Identify which units were truncated
#'
#' @description `is_unit_truncated()` returns a logical vector indicating which
#'   observations had their propensity scores modified by truncation. Use
#'   [is_ps_truncated()] to test whether an object has been truncated at all.
#'
#'   The answer comes from the truncation record, which is written for a fixed
#'   set of observations and can both be lost and outlive them. On a `ps_trunc`,
#'   [vctrs::vec_slice()] and [c()] drop it, and subassignment that grows the
#'   vector carries it across a length change; see [ps_trunc()] for the whole
#'   contract. On a [psw] vector built from truncated propensity scores, a
#'   subset drops it, while subassignment that grows the weights carries it
#'   across the length change.
#'
#'   `is_unit_truncated()` therefore checks that the record covers the object it
#'   is given, and raises an error of class `propensity_missing_meta_error` when
#'   it does not, or when an object marked as truncated carries no record at
#'   all, rather than name truncated units at stale positions. Query the
#'   `ps_trunc` object the record was written for instead.
#'
#'   That check counts observations, which a reordering does not change, so it
#'   does not catch one. An operation that reorders through vctrs rather than
#'   through `[`, such as `vctrs::vec_slice(x, 5:1)` or `dplyr::arrange()`,
#'   keeps a record written for the old order, and a `psw` keeps one through any
#'   same-length operation, a reordering included. `is_unit_truncated()` answers
#'   from those positions and names the wrong units. See [ps_trunc()] and [psw]
#'   for the whole contract.
#'
#' @param x A `ps_trunc` object created by [ps_trunc()], or a [psw] vector built
#'   from one.
#' @return A logical vector the same length as `x` (or number of rows for
#'   matrix input). `TRUE` marks observations whose values were winsorized.
#'
#' @seealso [ps_trunc()], [is_ps_truncated()], [ps_trunc_meta()]
#'
#' @examples
#' ps <- c(0.02, 0.3, 0.5, 0.7, 0.98)
#' ps_t <- ps_trunc(ps, method = "ps", lower = 0.05, upper = 0.95)
#' is_unit_truncated(ps_t)
#'
#' @export
is_unit_truncated <- function(x) {
  UseMethod("is_unit_truncated")
}

#' @export
is_unit_truncated.default <- function(x) {
  abort(
    "{.code is_unit_truncated()} not supported for class {.val {class(x)}}",
    error_class = "propensity_method_error"
  )
}

#' @export
is_unit_truncated.ps_trunc <- function(x) {
  # No observations, no answers. A record kept on an empty vector describes
  # observations it does not have, and indexing an empty logical by the
  # positions it names would grow one padded with `NA`.
  if (length(x) == 0) {
    return(logical(0))
  }

  meta <- ps_trunc_meta(x)
  check_trunc_record(meta, length(x), "is_unit_truncated()")

  out <- vector("logical", length = length(x))
  out[meta$truncated_idx] <- TRUE

  out
}

#' @export
is_unit_truncated.ps_trunc_matrix <- function(x) {
  if (nrow(x) == 0) {
    return(logical(0))
  }

  meta <- ps_trunc_meta(x)
  check_trunc_record(meta, nrow(x), "is_unit_truncated()")

  out <- vector("logical", length = nrow(x))
  out[meta$truncated_idx] <- TRUE

  out
}


#' @export
`[.ps_trunc_matrix` <- function(x, i, j, ..., drop = TRUE) {
  # Get metadata
  meta <- ps_trunc_meta(x)

  # Handle single index (matrix as vector) - bypass ps_trunc method
  if (nargs() == 2) {
    return(unclass(x)[i])
  }

  # Handle missing i (all rows) - bypass ps_trunc method
  if (missing(i)) {
    return(unclass(x)[, j, ..., drop = drop])
  }

  # For single element extraction, call base method directly to avoid ps_trunc method
  # Check if this will result in a single element
  if (!missing(j) && length(i) == 1 && length(j) == 1) {
    return(as.numeric(unclass(x)[i, j, drop = TRUE]))
  }

  # Perform subsetting with base matrix method to avoid calling [.ps_trunc
  result <- unclass(x)[i, j, ..., drop = drop]

  # If not a matrix anymore or dropped dimensions, return as-is (no metadata)
  if (!is.matrix(result)) {
    return(result)
  }

  # If result is a single element, return as numeric (no metadata)
  if (length(result) == 1) {
    return(as.numeric(result))
  }

  # The rows the result was built from, as positions in `x`, which is what
  # re-indexing the record takes. These are worked out from the subscript rather
  # than read off the result, so a form they do not reproduce, seen as a count
  # that does not match the rows base returned, degrades the record rather than
  # certifying it against rows it does not describe.
  rows <- subscript_row_positions(i, x)
  new_meta <- if (length(rows) == nrow(result)) {
    carry_trunc_record(meta, nrow(x), rows)
  } else {
    discard_trunc_record(meta, nrow(result))
  }

  attr(result, "ps_trunc_meta") <- new_meta
  class(result) <- c("ps_trunc_matrix", "ps_trunc", "matrix")
  result
}


# Print methods for ps_trunc_matrix

#' @export
print.ps_trunc_matrix <- function(x, ..., n = NULL) {
  meta <- ps_trunc_meta(x)
  n_rows <- nrow(x)
  k <- ncol(x)
  n_truncated <- length(meta$truncated_idx)

  # Create header
  if (meta$method == "pctl") {
    cat(sprintf(
      "<ps_trunc_matrix[%d x %d]; truncated %d of %d; method=%s[%.2f,%.2f]>\n",
      n_rows,
      k,
      n_truncated,
      n_rows,
      meta$method,
      meta$lower_pctl,
      meta$upper_pctl
    ))
  } else {
    cat(sprintf(
      "<ps_trunc_matrix[%d x %d]; truncated %d of %d; method=%s[%.4f,%.4f]>\n",
      n_rows,
      k,
      n_truncated,
      n_rows,
      meta$method,
      meta$lower_bound,
      ifelse(is.na(meta$upper_bound), Inf, meta$upper_bound)
    ))
  }

  # Determine how many rows to show
  if (is.null(n)) {
    n_print <- getOption("propensity.print_max", default = 10)
  } else {
    n_print <- n
  }

  # The header summarizes the record, and the scores are what is left to show.
  # The record itself is a set of index vectors as long as the data, so printed
  # after the matrix it is a second and longer object no one asked to see.
  x_print <- unclass(x)
  attr(x_print, "ps_trunc_meta") <- NULL

  # Show all rows if n is Inf or very large
  if (is.infinite(n_print) || n_print >= n_rows) {
    print(x_print)
  } else {
    # Show first n_print rows
    n_show <- min(n_print, n_rows)
    print(x_print[seq_len(n_show), , drop = FALSE])

    if (n_rows > n_show) {
      cat("# ... with", n_rows - n_show, "more rows\n")
    }
  }

  invisible(x)
}

# vctrs machinery for ps_trunc

#' @export
vec_ptype_abbr.ps_trunc <- function(x, ...) {
  "ps_trunc"
}

#' @export
vec_ptype_full.ps_trunc <- function(x, ...) {
  m <- ps_trunc_meta(x)

  # A bound read off the scores carries every digit of the score it came from,
  # which is more of the line a vector is printed under than the bound is worth.
  paste0(
    "ps_trunc{[",
    signif(m$lower_bound, 3),
    ",",
    signif(m$upper_bound, 3),
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
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc ps_trunc
vec_arith.ps_trunc.ps_trunc <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc MISSING
vec_arith.ps_trunc.MISSING <- function(op, x, y, ...) {
  switch(
    op,
    `-` = -1 * vec_data(x), # Returns numeric
    `+` = vec_data(x), # Returns numeric
    stop_incompatible_op(op, x, y)
  )
}

#' @export
#' @method vec_arith.ps_trunc numeric
vec_arith.ps_trunc.numeric <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_trunc
vec_arith.numeric.ps_trunc <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trunc integer
vec_arith.ps_trunc.integer <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

# Combining / Casting

# How a truncation is described, as opposed to which units it pinned: the
# method, the bounds it applied, and the percentiles those bounds came from. A
# bound is read off the scores for every method but `"ps"`, so two objects that
# agree on the method can still have been bounded at different places.
trunc_parameters <- function(meta) {
  fields <- c(
    "method",
    "lower_bound",
    "upper_bound",
    "lower_pctl",
    "upper_pctl"
  )

  rlang::set_names(lapply(fields, function(field) meta[[field]]), fields)
}

#' @export
vec_ptype2.ps_trunc.ps_trunc <- function(x, y, ...) {
  x_meta <- ps_trunc_meta(x)
  y_meta <- ps_trunc_meta(y)

  # Check if truncation parameters match
  if (!identical(trunc_parameters(x_meta), trunc_parameters(y_meta))) {
    warn_incompatible_metadata(
      x,
      y,
      "different truncation parameters"
    )
    return(double())
  }

  # The prototype carries the description of the truncation across and holds no
  # observations, so it names no positions: it is shared by inputs whose
  # observations are appended one after another, and the positions either record
  # names would describe units from the other input. Truncating an empty vector
  # again would work the bounds out from scores that are not there, and the rule
  # defined against the exposure has none to be handed.
  new_ps_trunc(double(), meta = drop_trunc_record(x_meta))
}

#' @export
vec_ptype2.ps_trunc.double <- function(x, y, ...) {
  warn_class_downgrade("ps_trunc")
  double()
}

#' @export
vec_ptype2.double.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade("ps_trunc")
  double()
}

# A cast returns the values it was handed in the type it was handed, and a
# `ps_trunc`'s type is the whole description of the truncation. That is the
# comparison `vec_ptype2()` already makes when it refuses to find a common type,
# so the cast makes it too: a cast comparing less than the combine does hands
# `x` back describing itself under the target's name. The positional half of the
# record describes the values arriving rather than the type they arrive in, so
# it is left out of the comparison, which is also what lets a prototype built by
# `drop_trunc_record()` be cast to.
#
# `vec_ptype_full()` names none of what is compared, so the two types render
# identically and the refusal would read as a type that cannot be converted to
# itself. What disagrees is named alongside them, the way the combine names it.
#' @export
vec_cast.ps_trunc.ps_trunc <- function(x, to, ...) {
  x_meta <- ps_trunc_meta(x)
  to_meta <- ps_trunc_meta(to)

  if (!identical(trunc_parameters(x_meta), trunc_parameters(to_meta))) {
    vctrs::stop_incompatible_cast(
      x,
      to,
      x_arg = "",
      to_arg = "",
      details = "different truncation parameters"
    )
  }

  x
}

#' @export
vec_cast.double.ps_trunc <- function(x, to, ...) {
  vec_data(x)
}

# A cast returns the values it was handed in the type it was handed, and the
# truncation is part of that type, so the description comes from the target. The
# values are the ones arriving, so the positions are written for them: none of
# them was pinned to a bound on the way here.
trunc_meta_for_cast <- function(to, x) {
  c(
    drop_trunc_record(ps_trunc_meta(to)),
    list(
      truncated_idx = integer(0),
      n_obs = length(x)
    )
  )
}

#' @export
vec_cast.ps_trunc.double <- function(x, to, ...) {
  new_ps_trunc(x, meta = trunc_meta_for_cast(to, x))
}

#' @export
vec_ptype2.psw.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade(c("psw", "ps_trunc"))
  double()
}

#' @export
vec_ptype2.ps_trunc.psw <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trunc", "psw"))
  double()
}

#' @export
vec_cast.character.ps_trunc <- function(x, to, ...) as.character(vec_data(x))

# Propensity scores lie strictly between 0 and 1, so meeting an integer in the
# integers would round every one of them away. The common type is the one that
# holds both sets of values.
#' @export
vec_ptype2.ps_trunc.integer <- function(x, y, ...) {
  warn_class_downgrade("ps_trunc")
  double()
}

#' @export
vec_ptype2.integer.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade("ps_trunc")
  double()
}

#' @export
vec_cast.integer.ps_trunc <- function(x, to, ...) {
  # A propensity score has no integer to be, so vctrs' own check reports the
  # loss rather than silently rounding it away.
  vec_cast(vec_data(x), integer(), x_arg = "ps_trunc")
}

#' @export
vec_cast.ps_trunc.integer <- function(x, to, ...) {
  xx <- as.double(x)
  new_ps_trunc(xx, meta = trunc_meta_for_cast(to, xx))
}

#' @export
vec_math.ps_trunc <- function(.fn, .x, ...) {
  vec_math_base(.fn, vec_data(.x), ...)
}

#' @export
Summary.ps_trunc <- function(..., na.rm = FALSE) {
  .fn <- .Generic
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call(.fn, c(numeric_args, list(na.rm = na.rm)))
}

#' @export
min.ps_trunc <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("min", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
max.ps_trunc <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("max", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
range.ps_trunc <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("range", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
median.ps_trunc <- function(x, na.rm = FALSE, ...) {
  median(vec_data(x), na.rm = na.rm, ...)
}

#' @export
quantile.ps_trunc <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, ...) {
  quantile(vec_data(x), probs = probs, na.rm = na.rm, ...)
}


# Note: We cannot implement custom vec_c for ps_trunc
# vctrs handles combination internally through vec_ptype2 and vec_cast
# When combining ps_trunc objects with same parameters, indices won't be preserved

#' @export
`[.ps_trunc` <- function(x, i, ...) {
  # If i is missing, just call NextMethod
  if (missing(i)) {
    return(NextMethod())
  }

  # Get original metadata
  meta <- ps_trunc_meta(x)

  # Convert i to positive integer indices if needed
  # For logical indices, only allow length 1 or length n
  if (is.logical(i) && length(i) != 1 && length(i) != length(x)) {
    abort(
      "Logical subscript `i` must be size 1 or {length(x)}, not {length(i)}.",
      error_class = "propensity_length_error"
    )
  }
  i <- vec_as_location(i, n = length(x))

  # The subscript is in hand here, which is what re-indexing a record takes and
  # what the restore behind `NextMethod()` lacks: that restore would see a
  # record written for a different number of observations and drop it. Building
  # the result here keeps ordinary subsetting describing what it returns.
  new_ps_trunc(vec_data(x)[i], carry_trunc_record(meta, length(x), i))
}

#' @export
unique.ps_trunc <- function(x, incomparables = FALSE, ...) {
  check_incomparables(incomparables, "ps_trunc")

  # `vec_unique_loc()` names the position each retained value came from, which
  # is the subscript re-indexing the record takes. Without this the restore
  # behind vctrs' own method sees only a shorter vector and drops the record.
  x[vec_unique_loc(x)]
}

#' @export
rep.ps_trunc <- function(x, ...) {
  # The positions `rep()` produces are the subscript the result is built from,
  # so the record follows them and a repeated unit is reported at each place it
  # now holds.
  x[rep(seq_along(x), ...)]
}

#' @export
summary.ps_trunc <- function(object, ...) {
  summary(vec_data(object), ...)
}

#' @export
sort.ps_trunc <- function(x, decreasing = FALSE, na.last = NA, ...) {
  meta <- ps_trunc_meta(x)
  x_data <- vec_data(x)

  # `order()` returns the old positions in their new order, so it is the
  # subscript the result is built from and the record follows it.
  ord <- order(x_data, na.last = na.last, decreasing = decreasing, ...)

  new_ps_trunc(x_data[ord], carry_trunc_record(meta, length(x), ord))
}

#' @export
vec_restore.ps_trunc <- function(x, to, ...) {
  # vec_data in case x is already a vctr
  data <- vec_data(x)
  meta <- ps_trunc_meta(to)

  # Nothing rebuilding a `ps_trunc` is handed the subscript behind a length
  # change, so a record written for a different number of observations cannot be
  # re-indexed onto the data arriving here. Zero-length data is exempt: a
  # prototype or an empty slice holds no observations, so no position in the
  # record contradicts it, and the record rides along to the restore that builds
  # the real result.
  if (length(data) > 0 && !record_covers(meta, length(data))) {
    meta <- discard_trunc_record(meta, length(data))
  }

  new_ps_trunc(data, meta)
}

#' @export
anyDuplicated.ps_trunc <- function(x, incomparables = FALSE, ...) {
  anyDuplicated(vec_data(x), incomparables = incomparables, ...)
}

#' @export
diff.ps_trunc <- function(x, lag = 1L, differences = 1L, ...) {
  diff(vec_data(x), lag = lag, differences = differences, ...)
}
