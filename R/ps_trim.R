#' Trim Propensity Scores
#'
#' @description
#' Trim observations with extreme propensity scores by replacing them with `NA`,
#' effectively removing those units from downstream analyses. The returned object
#' has the same length (or dimensions) as the input, with trimmed entries set to
#' `NA`. After trimming, refit the propensity score model on the retained
#' observations with [ps_refit()].
#'
#' @param ps A numeric vector of propensity scores in (0, 1) for binary
#'   exposures, or a matrix / data frame where each column gives the propensity
#'   score for one level of a categorical exposure.
#' @param method Trimming method. One of:
#'
#'   * **`"ps"`** (default): Fixed threshold. Observations with propensity
#'     scores outside `[lower, upper]` are trimmed. For categorical exposures,
#'     observations where *any* column falls below `lower` (the symmetric
#'     threshold delta) are trimmed.
#'   * **`"adaptive"`**: Data-driven threshold that minimizes the asymptotic
#'     variance of the IPW estimator (Crump et al., 2009). The `lower` and
#'     `upper` arguments are ignored.
#'   * **`"pctl"`**: Quantile-based. Observations outside the `[lower, upper]`
#'     quantiles of the propensity score distribution are trimmed. Defaults:
#'     `lower = 0.05`, `upper = 0.95`.
#'   * **`"pref"`**: Preference score trimming. Transforms propensity scores
#'     to the preference scale (Walker et al., 2013) and trims outside
#'     `[lower, upper]`. Requires `.exposure`. Binary exposures only. Defaults:
#'     `lower = 0.3`, `upper = 0.7`.
#'   * **`"cr"`**: Common range (clinical equipoise). Trims to the overlap
#'     region of the propensity score distributions across exposure groups.
#'     Requires `.exposure`. Binary exposures only. The `lower` and `upper`
#'     arguments are ignored.
#'   * **`"optimal"`**: Multi-category optimal trimming (Yang et al., 2016).
#'     Categorical exposures only. Requires `.exposure`.
#'
#'   For categorical exposures, only `"ps"` and `"optimal"` are supported.
#' @param lower,upper Numeric thresholds whose interpretation depends on
#'   `method`:
#'
#'   * `"ps"`: absolute propensity score bounds (defaults: 0.1, 0.9). For
#'     categorical exposures, only `lower` is used, as the symmetric threshold
#'     delta, and it defaults to 0.1. That default deliberately differs from the
#'     0.01 threshold [ps_trunc()] uses for categorical exposures: trimming
#'     discards the units it selects, so its default follows common-support
#'     trimming practice, whereas truncation keeps every unit and only pins the
#'     most extreme scores back. With `k` exposure levels, a threshold of `1/k`
#'     or larger cannot be met by every column of a row that sums to one, and is
#'     an error.
#'   * `"pctl"`: quantile probabilities (defaults: 0.05, 0.95).
#'   * `"pref"`: preference score bounds (defaults: 0.3, 0.7).
#'   * `"adaptive"`, `"cr"`, `"optimal"`: ignored (thresholds are data-driven).
#' @param .exposure An exposure variable. Required for `"pref"`, `"cr"` (binary
#'   vector), and `"optimal"` (factor or character). Not required for other
#'   methods.
#' @inheritParams wt_ate
#' @param .focal_level The value of `.exposure` representing the focal
#'   (treated) group, used by `"pref"` and `"cr"`. Every binary coding honors
#'   it: 0/1 numeric, logical, two-level factor, and two-level character
#'   exposures are all coded with the named level as focal, and a level the
#'   exposure never takes is an error. With no level named, a binary exposure
#'   defaults to its higher level, which is `1` for a 0/1 exposure, `TRUE` for
#'   a logical one, and the second of the two levels a factor or character
#'   exposure takes. Levels a factor declares but never takes are not
#'   candidates. Naming any other level reverses the coding, so `ps` must then
#'   hold the probability of the named level.
#' @param .reference_level The value of `.exposure` representing the reference
#'   (control) group. Naming it makes the exposure's other level focal, with
#'   the same consequence for `ps`, and a level the exposure never takes is an
#'   error. Automatically detected if not supplied.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' ## How trimming works
#'
#' Trimming identifies observations with extreme (near 0 or 1) propensity
#' scores and sets them to `NA`. These observations are excluded from
#' subsequent weight calculations and effect estimation. The goal is to
#' remove units that lack sufficient overlap between exposure groups, which
#' would otherwise receive extreme weights and destabilize estimates.
#'
#' ## Choosing a method
#'
#' * Use `"ps"` when you have a specific threshold in mind or want a simple
#'   default.
#' * Use `"adaptive"` for a principled, data-driven cutoff that targets
#'   variance reduction.
#' * Use `"pctl"` to trim a fixed percentage of extreme values from each tail.
#' * Use `"pref"` when you want to restrict to the region of clinical
#'   equipoise based on the preference score.
#' * Use `"cr"` to restrict to the common support region where both exposure
#'   groups have observed propensity scores.
#' * Use `"optimal"` for multi-category (3+) exposures; this is the only
#'   data-driven method available for categorical treatments.
#'
#' ## Typical workflow
#'
#' 1. Fit a propensity score model
#' 2. Apply `ps_trim()` to flag extreme values
#' 3. Call [ps_refit()] to re-estimate propensity scores on the retained sample
#' 4. Compute weights with [wt_ate()] or another weight function
#'
#' ## Object behavior
#'
#' Arithmetic operations on `ps_trim` objects return plain numeric vectors,
#' since transformed propensity scores (e.g., `1/ps`) are no longer propensity
#' scores. Trimmed values propagate as `NA` in calculations; use `na.rm = TRUE`
#' where appropriate.
#'
#' When combining `ps_trim` objects with [c()], trimming parameters must match.
#' Mismatched parameters trigger a warning and return a numeric vector.
#'
#' Use [ps_trim_meta()] to inspect the trimming metadata, including the method,
#' cutoffs, and which observations were retained or trimmed.
#'
#' ## Missing values
#'
#' A propensity score that arrives missing is not one this function removed, so
#' it takes no part in the trimming record: it joins neither the retained nor
#' the trimmed positions, and [is_unit_trimmed()] reports `FALSE` for it. The
#' value propagates as `NA` into the result. For a matrix of categorical
#' propensity scores, a row with a missing cell comes back exactly as it
#' arrived, since there is no complete probability vector to place against a
#' threshold.
#'
#' The methods that read a cutoff off the scores read it off the scores they
#' have. `"adaptive"`, `"pctl"`, `"cr"`, and `"optimal"` all work their cutoffs
#' out from the complete scores or rows, so the cutoff is the one the same call
#' would produce with the missing observations dropped.
#'
#' `"pref"` centers its preference scores on the proportion exposed across the
#' whole sample, which is a fact about the exposure rather than about the
#' scores, so a unit whose propensity score is missing still counts toward it
#' and the cutoff is not the one the shorter call would produce. That unit's own
#' preference score is missing, which leaves it outside both the retained and
#' the trimmed positions like any other missing score.
#'
#' A missing exposure is a different matter, and `"pref"` and `"cr"` refuse one
#' with an error of class `propensity_missing_value_error`. Their cutoffs come
#' from the exposure groups, and a unit that belongs to neither leaves them
#' undefined. Remove or impute the missing exposure values first.
#'
#' ## The trimming record
#'
#' A `ps_trim` records which units were trimmed as positions among the
#' observations it was written for, along with how many observations that was.
#' Operations that hand this package the subscript re-index those positions onto
#' the result: subsetting with `[`, [sort()], [unique()], [rep()], and
#' [na.omit()] all return a record written for what they return, and a subscript
#' naming a position more than once reports that unit at every place it now
#' holds.
#'
#' Operations that change how many observations there are without supplying a
#' subscript cannot re-index the record, and it is dropped rather than worked
#' out from the values, since reading membership back from the `NA` pattern
#' would report a propensity score that arrived missing as one this package
#' removed. [vctrs::vec_slice()], which is how filtering, joining, and grouped
#' summaries in dplyr reach a column, is the usual route, and dropping the
#' record there raises a warning of class `propensity_trim_record_warning`.
#' Combining with [c()] drops it without comment, because concatenation appends
#' one set of observations to another and the prototype it builds the result
#' from holds no positions to lose. The values, the class, and the method and
#' its cutoffs are untouched either way.
#'
#' A record can also outlive the observations it describes, because it travels
#' by routes vctrs does not see: growing a `ps_trim` by subassignment carries it
#' across the length change. [is_unit_trimmed()] and [ps_refit()] therefore
#' check that the record covers the object they are given and raise an error of
#' class `propensity_missing_meta_error` when it does not, rather than name
#' trimmed units at stale positions.
#'
#' That check compares how many observations the record was written for against
#' how many the object holds, which a reordering does not change. An operation
#' that reorders the observations through vctrs, rather than through `[`,
#' therefore keeps a record written for the order they used to be in:
#' `vctrs::vec_slice(x, 5:1)` and `dplyr::arrange()` both return the values in
#' a new order under positions still naming the old one. [is_unit_trimmed()]
#' answers from those positions and names the wrong units, and [ps_refit()]
#' refits on the wrong rows. Subsetting with `[` is handed the subscript and
#' re-indexes, so reorder with `[`, or put the propensity scores in the order
#' you want before trimming them.
#'
#' @return A **`ps_trim`** object (a numeric vector with class `"ps_trim"`, or a
#'   matrix with class `"ps_trim_matrix"`). Trimmed observations are `NA`.
#'   Metadata is stored in the `"ps_trim_meta"` attribute and can be accessed
#'   with [ps_trim_meta()]. Key fields include:
#'
#'   * `method`: the trimming method used
#'   * `keep_idx`: integer indices of retained observations
#'   * `trimmed_idx`: integer indices of trimmed (NA) observations
#'   * `n_obs`: the number of observations those indices describe
#'   * Method-specific fields such as `cutoff` (adaptive), `q_lower`/`q_upper`
#'     (pctl), `cr_lower`/`cr_upper` (cr), `delta` (categorical ps),
#'     or `lambda` (optimal)
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
#' Yang, S., Imbens, G. W., Cui, Z., Faries, D. E., & Kadziola, Z. (2016).
#' Propensity score matching and subclassification in observational studies
#' with multi-level treatments. *Biometrics*, 72(4), 1055--1065.
#'
#' @seealso [ps_trunc()] for bounding (winsorizing) instead of discarding,
#'   [ps_refit()] to re-estimate propensity scores after trimming,
#'   [ps_calibrate()] for calibration-based adjustment,
#'   [ps_trim_meta()] to inspect trimming metadata,
#'   [is_ps_trimmed()] and [is_unit_trimmed()] for logical queries.
#'
#' @examples
#' set.seed(2)
#' n <- 300
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(1.3 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' # Fixed threshold trimming (default)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#' trimmed
#'
#' # How many observations were trimmed?
#' sum(is_unit_trimmed(trimmed))
#'
#' # Data-driven adaptive trimming
#' ps_trim(ps, method = "adaptive")
#'
#' # Quantile-based trimming at 5th and 95th percentiles
#' ps_trim(ps, method = "pctl")
#'
#' # Refit after trimming, then compute weights
#' trimmed <- ps_trim(ps, method = "adaptive")
#' refitted <- ps_refit(trimmed, fit)
#' wt_ate(refitted, .exposure = z)
#'
#' @export
ps_trim <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  UseMethod("ps_trim")
}

#' @export
ps_trim.default <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  method <- rlang::arg_match(method)

  # Optimal trimming is defined over the rows of a propensity score matrix, so
  # there is nothing for it to do with a vector of scores.
  if (method == "optimal") {
    abort(
      c(
        "Method {.val optimal} is only supported for categorical exposures.",
        i = "Supply the propensity scores as a matrix or data frame with one column per exposure level."
      ),
      error_class = "propensity_wt_not_supported_error"
    )
  }

  check_ps_range(ps)

  # Handle deprecation
  focal_params <- handle_focal_deprecation(
    .focal_level,
    .reference_level,
    .treated,
    .untreated,
    "ps_trim"
  )
  .focal_level <- focal_params$.focal_level
  .reference_level <- focal_params$.reference_level
  check_focal_levels(.focal_level, .reference_level)

  if (method == "ps") {
    if (is.null(lower)) {
      lower <- 0.1
    }
    if (is.null(upper)) {
      upper <- 0.9
    }
    check_lower_upper(lower, upper)
  } else if (method == "adaptive") {
    if (!is.null(lower) || !is.null(upper)) {
      warn(
        "For {.code method = 'adaptive'}, {.code lower} and {.code upper} are ignored."
      )
    }
  } else if (method == "pctl") {
    if (is.null(lower)) {
      lower <- 0.05
    }
    if (is.null(upper)) {
      upper <- 0.95
    }
    check_quantile_probs(lower, upper)
    check_lower_upper(lower, upper)
  } else if (method == "pref") {
    if (is.null(lower)) {
      lower <- 0.3
    }
    if (is.null(upper)) {
      upper <- 0.7
    }
    check_lower_upper(lower, upper)
  } else if (method == "cr") {
    if (!is.null(lower) || !is.null(upper)) {
      warn(
        "For {.code method = 'cr'}, {.code lower} and {.code upper} are ignored."
      )
    }
  }

  n <- length(ps)
  meta_list <- list(method = method, lower = lower, upper = upper)

  # A score that arrived missing is not one this function can place against a
  # cutoff, so it takes no part in working the cutoff out and no part in the
  # record. Every rule below compares scores with `which()`, which leaves a
  # missing comparison out of the retained positions on its own; the trimmed
  # positions are then everything else that was observed.
  observed <- !is.na(ps)

  # Decide which indices are kept
  if (method == "ps") {
    keep_idx <- which(ps >= lower & ps <= upper)
  } else if (method == "adaptive") {
    sum_wt <- 1 / (ps[observed] * (1 - ps[observed]))
    k <- 2 * mean(sum_wt) - max(sum_wt)

    if (k >= 0) {
      cutoff <- 0
    } else {
      trim_fun <- function(x) {
        sum_wt_trim <- sum_wt[sum_wt <= x]
        x - 2 * mean(sum_wt_trim)
      }
      rng <- range(sum_wt)
      lambda <- uniroot(trim_fun, lower = rng[1], upper = rng[2])$root
      cutoff <- 0.5 - sqrt(0.25 - 1 / lambda)
    }
    meta_list$cutoff <- cutoff
    keep_idx <- which(pmin(ps, 1 - ps) > cutoff)
  } else if (method == "pctl") {
    q_lower <- quantile(ps, probs = lower, na.rm = TRUE)
    q_upper <- quantile(ps, probs = upper, na.rm = TRUE)
    meta_list$q_lower <- q_lower
    meta_list$q_upper <- q_upper
    keep_idx <- which(ps >= q_lower & ps <= q_upper)
  } else if (method == "pref") {
    if (is.null(.exposure)) {
      abort(
        "For {.code method = 'pref'}, must supply {.arg .exposure}.",
        error_class = "propensity_missing_arg_error"
      )
    }
    check_exposure_complete(.exposure, method)
    .exposure <- transform_exposure_binary(
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
    prop_exposure <- mean(.exposure)
    pref_score <- plogis(qlogis(ps) - qlogis(prop_exposure))
    meta_list$P <- prop_exposure
    keep_idx <- which(pref_score >= lower & pref_score <= upper)
  } else if (method == "cr") {
    if (is.null(.exposure)) {
      abort(
        "For {.code method = 'cr'}, must supply {.arg .exposure}.",
        error_class = "propensity_missing_arg_error"
      )
    }
    check_exposure_complete(.exposure, method)
    .exposure <- transform_exposure_binary(
      .exposure,
      .focal_level = .focal_level,
      .reference_level = .reference_level
    )
    ps_treat <- ps[observed & .exposure == 1]
    ps_untrt <- ps[observed & .exposure == 0]
    cr_lower <- min(ps_treat)
    cr_upper <- max(ps_untrt)
    meta_list$cr_lower <- cr_lower
    meta_list$cr_upper <- cr_upper

    keep_idx <- which(ps >= cr_lower & ps <= cr_upper)
  }

  trimmed_idx <- setdiff(seq_len(n), c(keep_idx, which(!observed)))

  # Replace trimmed entries with NA
  ps_na <- ps
  ps_na[trimmed_idx] <- NA_real_

  new_trimmed_ps(
    x = ps_na,
    ps_trim_meta = c(
      meta_list,
      list(
        keep_idx = keep_idx,
        trimmed_idx = trimmed_idx,
        n_obs = n
      )
    )
  )
}

#' @export
ps_trim.matrix <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # The generic offers every method, so match against that full set and then
  # reject the ones the categorical path does not define.
  method <- rlang::arg_match(
    method,
    values = c("ps", "adaptive", "pctl", "pref", "cr", "optimal")
  )
  if (!method %in% c("ps", "optimal")) {
    abort(
      c(
        "Method {.val {method}} is not supported for categorical exposures.",
        i = "Use {.val ps} or {.val optimal}."
      ),
      error_class = "propensity_method_error"
    )
  }

  # Validate exposure for categorical
  if (is.null(.exposure)) {
    abort(
      "`.exposure` must be provided for categorical propensity score trimming.",
      error_class = "propensity_missing_arg_error"
    )
  }

  # Transform to factor and validate
  .exposure <- transform_exposure_categorical(.exposure)

  # Validate matrix
  ps <- check_ps_matrix(ps, .exposure, call = rlang::current_env())

  n <- nrow(ps)
  k <- ncol(ps)

  # A row with a missing score has no complete probability vector to place
  # against a threshold, so it takes no part in working the threshold out and no
  # part in the record. Both rules below compare rows with `which()`, which
  # leaves a missing comparison out of the retained positions on its own, and the
  # group-preservation reset falls back to the complete rows for the same reason.
  incomplete_rows <- which(apply(ps, 1, anyNA))
  complete_rows <- setdiff(seq_len(n), incomplete_rows)

  # Initialize metadata
  meta_list <- list(method = method, is_matrix = TRUE)

  if (method == "ps") {
    # Symmetric trimming
    if (is.null(lower)) {
      lower <- 0.1
    }
    delta <- lower # Use lower as delta for consistency

    # A threshold at or above 1/k cannot be met by every column of a row that
    # sums to one, so there is no trimming rule left to apply.
    if (delta >= 1 / k) {
      abort(
        "Invalid trimming threshold (delta >= 1/k).",
        error_class = "propensity_range_error"
      )
    }

    # Apply symmetric trimming rule: keep if min(propensity scores) > delta
    keep_idx <- which(apply(ps, 1, function(x) min(x) > delta))

    # Check if all treatment groups are preserved
    if (length(unique(.exposure[keep_idx])) < k) {
      warn(
        "One or more groups removed after trimming; returning original data",
        warning_class = "propensity_no_data_warning"
      )
      keep_idx <- complete_rows
    }

    meta_list$delta <- delta
  } else {
    # optimal
    # Multi-category optimal trimming (Yang et al., 2016)
    # Calculate sum of inverse propensity scores
    sum_inv_ps <- rowSums(1 / ps)
    sum_inv_complete <- sum_inv_ps[complete_rows]

    # Define trimming function
    trim_fun <- function(x) {
      sum_trim <- sum_inv_complete[sum_inv_complete <= x]
      if (length(sum_trim) == 0) {
        return(x)
      }
      x - 2 * mean(sum_trim) / mean(sum_inv_complete <= x)
    }

    # Check if trimming is needed
    if (trim_fun(max(sum_inv_complete)) < 0) {
      # No valid solution, use maximum + 1
      lambda <- max(sum_inv_complete) + 1
      keep_idx <- complete_rows # Keep all
    } else {
      # Find optimal lambda
      result <- tryCatch(
        {
          uniroot(
            trim_fun,
            lower = min(sum_inv_complete),
            upper = max(sum_inv_complete)
          )$root
        },
        error = function(e) {
          warn(
            "Could not find optimal trimming threshold; using no trimming",
            warning_class = "propensity_convergence_warning"
          )
          NULL
        }
      )

      if (!is.null(result)) {
        lambda <- result
        keep_idx <- which(sum_inv_ps <= lambda)

        # Check if all treatment groups are preserved
        if (length(unique(.exposure[keep_idx])) < k) {
          warn(
            "One or more groups removed after trimming; returning original data",
            warning_class = "propensity_no_data_warning"
          )
          keep_idx <- complete_rows
          lambda <- NULL
        }
      } else {
        keep_idx <- complete_rows
        lambda <- NULL
      }
    }

    meta_list$lambda <- lambda
  }

  trimmed_idx <- setdiff(seq_len(n), c(keep_idx, incomplete_rows))

  # Replace trimmed entries with NA
  ps_na <- ps
  ps_na[trimmed_idx, ] <- NA_real_

  new_trimmed_ps(
    x = ps_na,
    ps_trim_meta = c(
      meta_list,
      list(
        keep_idx = keep_idx,
        trimmed_idx = trimmed_idx,
        n_obs = n
      )
    )
  )
}

#' @export
ps_trim.data.frame <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  # For categorical exposures, convert to matrix and call matrix method
  if (!is.null(.exposure)) {
    exposure_type <- detect_exposure_type(.exposure)
    if (exposure_type == "categorical") {
      ps_matrix <- as.matrix(ps)
      return(ps_trim.matrix(
        ps = ps_matrix,
        method = method,
        lower = lower,
        upper = upper,
        .exposure = .exposure,
        ...
      ))
    }
  }

  # For binary exposures, extract appropriate column and call default method
  # This is a simplified version - in practice you might want to handle
  # .propensity_col parameter like in weight functions
  if (ncol(ps) == 2) {
    # Use second column by default for binary
    ps_vec <- ps[[2]]
  } else {
    # Use first column
    ps_vec <- ps[[1]]
  }

  ps_trim.default(
    ps = ps_vec,
    method = method,
    lower = lower,
    upper = upper,
    .exposure = .exposure,
    .focal_level = .focal_level,
    .reference_level = .reference_level,
    ...,
    .treated = .treated,
    .untreated = .untreated
  )
}

#' @export
ps_trim.ps_trim <- function(
  ps,
  method = c("ps", "adaptive", "pctl", "pref", "cr", "optimal"),
  lower = NULL,
  upper = NULL,
  .exposure = NULL,
  .focal_level = NULL,
  .reference_level = NULL,
  ...,
  .treated = NULL,
  .untreated = NULL
) {
  warn(
    "Propensity scores have already been trimmed. Returning original object.",
    warning_class = "propensity_already_modified_warning"
  )
  ps
}

new_trimmed_ps <- function(x, ps_trim_meta = list()) {
  if (is.matrix(x)) {
    # For matrices, we don't use vctrs
    structure(
      x,
      ps_trim_meta = ps_trim_meta,
      class = c("ps_trim_matrix", "ps_trim", "matrix")
    )
  } else {
    vec_assert(x, ptype = double())
    new_vctr(
      x,
      ps_trim_meta = ps_trim_meta,
      class = "ps_trim",
      inherit_base_type = TRUE
    )
  }
}

# The positional half of a trimming record. The rest of it, the method and its
# cutoffs, describes the trimming rather than the units and means the same thing
# at any length.
reindex_trim_record <- function(meta, i) {
  meta$keep_idx <- reindex_positions(meta$keep_idx, i)
  meta$trimmed_idx <- reindex_positions(meta$trimmed_idx, i)
  meta$n_obs <- length(i)

  meta
}

drop_trim_record <- function(meta) {
  meta[["keep_idx"]] <- NULL
  meta[["trimmed_idx"]] <- NULL
  meta[["n_obs"]] <- NULL

  meta
}

# A record that no longer describes the observations in front of it is dropped
# rather than guessed at. Nothing in the values says which units a shorter or a
# longer vector once trimmed, and reading membership back from the `NA` pattern
# would report a propensity score that arrived missing as one this package
# removed.
#
# A record over no observations names no unit, so replacing it costs the caller
# nothing and goes without comment. That is the record every prototype carries,
# which is what concatenation builds its result from.
discard_trim_record <- function(meta, n) {
  recorded <- meta$n_obs

  if (!is.null(recorded) && recorded > 0) {
    warn(
      c(
        "Dropping the record of which units were trimmed.",
        i = "The record describes {recorded} observation{?s} and this result
             has {n}, so its positions do not describe them.",
        i = "The values are unchanged and the result is still a
             {.cls ps_trim}. Trim the propensity scores you want to work with
             to get a record written for them."
      ),
      warning_class = "propensity_trim_record_warning",
      # One of the routes here is vctrs' internal dispatch, whose call would be
      # reported and names nothing the caller wrote, so no call is attributed.
      call = NULL
    )
  }

  drop_trim_record(meta)
}

# A positional query reads its answer out of the record, so a record that does
# not cover the object in front of it has no answer to give: reporting every
# unit as retained would be a wrong answer rather than a missing one.
check_trim_record <- function(meta, n, fn, call = rlang::caller_env()) {
  if (record_covers(meta, n)) {
    return(invisible(meta))
  }

  recorded <- meta$n_obs
  problem <- if (is.null(recorded)) {
    "These propensity scores carry no record of which units were trimmed."
  } else {
    "The record describes {recorded} observation{?s} and these propensity
     scores have {n}, so its positions do not describe them."
  }

  abort(
    c(
      "{.code {fn}} has no usable trimming record for these propensity scores.",
      x = problem,
      i = "An operation that changed the number of observations dropped it.
           Trim the propensity scores you want to work with, or query the
           object the record was written for."
    ),
    error_class = "propensity_missing_meta_error",
    call = call
  )
}

#' Test whether propensity scores have been trimmed
#'
#' @description `is_ps_trimmed()` returns `TRUE` if `x` is a `ps_trim` object
#'   or a `psw` object created from trimmed propensity scores, and `FALSE`
#'   otherwise. This tests whether the *object* carries trimming information, not
#'   which individual units were trimmed; see [is_unit_trimmed()] for that.
#'
#' @param x An object to test.
#' @return A logical scalar (`TRUE` or `FALSE`).
#'
#' @seealso [ps_trim()] for trimming propensity scores, [is_unit_trimmed()] to
#'   identify which units were trimmed, [ps_trim_meta()] to retrieve full
#'   trimming metadata.
#'
#' @examples
#' ps <- c(0.05, 0.3, 0.6, 0.95)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' is_ps_trimmed(trimmed)
#' is_ps_trimmed(ps)
#'
#' @export
is_ps_trimmed <- function(x) {
  UseMethod("is_ps_trimmed")
}

#' @export
is_ps_trimmed.default <- function(x) {
  FALSE
}

#' @export
is_ps_trimmed.ps_trim <- function(x) {
  TRUE
}

#' @export
is_ps_trimmed.ps_trim_matrix <- function(x) {
  TRUE
}

#' Identify which units were trimmed
#'
#' @description `is_unit_trimmed()` returns a logical vector indicating which
#'   observations were removed by trimming. This is a per-unit query, as opposed
#'   to [is_ps_trimmed()], which tests whether the object has been trimmed at
#'   all.
#'
#'   The answer comes from the trimming record, which is written for a fixed set
#'   of observations and can both be lost and outlive them. On a `ps_trim`,
#'   [vctrs::vec_slice()] and [c()] drop it, and subassignment that grows the
#'   vector carries it across a length change; see [ps_trim()] for the whole
#'   contract. On a [psw] vector built from trimmed propensity scores, a subset
#'   drops it, while `model.frame()` re-attaches it to the shortened weights
#'   column of an outcome model fit on these weights.
#'
#'   `is_unit_trimmed()` therefore checks that the record covers the object it
#'   is given, and raises an error of class `propensity_missing_meta_error` when
#'   it does not, or when an object marked as trimmed carries no record at all,
#'   rather than name trimmed units at stale positions. Query the `ps_trim`
#'   object the record was written for instead.
#'
#'   That check counts observations, which a reordering does not change, so it
#'   does not catch one. A `ps_trim` reordered through vctrs rather than through
#'   `[`, by `vctrs::vec_slice(x, 5:1)` or `dplyr::arrange()`, keeps a record
#'   written for the old order, and a `psw` keeps one through any same-length
#'   operation, a reordering included. `is_unit_trimmed()` answers from those
#'   positions and names the wrong units. See [ps_trim()] and [psw] for the
#'   whole contract.
#'
#' @param x A `ps_trim` object created by [ps_trim()], or a [psw] vector built
#'   from one.
#' @return A logical vector the same length as `x`, where `TRUE` marks a
#'   trimmed unit.
#'
#' @seealso [ps_trim()] for trimming propensity scores, [is_ps_trimmed()] to
#'   test whether an object has been trimmed, [ps_trim_meta()] to retrieve full
#'   trimming metadata.
#'
#' @examples
#' ps <- c(0.05, 0.3, 0.6, 0.95)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' is_unit_trimmed(trimmed)
#'
#' # Use to subset data to retained observations
#' kept <- !is_unit_trimmed(trimmed)
#' ps[kept]
#'
#' @export
is_unit_trimmed <- function(x) {
  UseMethod("is_unit_trimmed")
}

#' @export
is_unit_trimmed.default <- function(x) {
  abort(
    "{.code is_unit_trimmed()} not supported for class {.val {class(x)}}",
    error_class = "propensity_method_error"
  )
}

#' @export
is_unit_trimmed.ps_trim <- function(x) {
  # No observations, no answers. A record kept on an empty vector describes
  # observations it does not have, and indexing an empty logical by the
  # positions it names would grow one padded with `NA`.
  if (length(x) == 0) {
    return(logical(0))
  }

  meta <- ps_trim_meta(x)
  check_trim_record(meta, length(x), "is_unit_trimmed()")

  out <- vector("logical", length = length(x))
  out[meta$trimmed_idx] <- TRUE

  out
}

#' @export
is_unit_trimmed.ps_trim_matrix <- function(x) {
  if (nrow(x) == 0) {
    return(logical(0))
  }

  meta <- ps_trim_meta(x)
  check_trim_record(meta, nrow(x), "is_unit_trimmed()")

  out <- vector("logical", length = nrow(x))
  out[meta$trimmed_idx] <- TRUE

  out
}


#' @export
`[.ps_trim_matrix` <- function(x, i, j, ..., drop = TRUE) {
  # Get metadata
  meta <- ps_trim_meta(x)

  # Handle single index (matrix as vector) - bypass ps_trim method
  if (nargs() == 2) {
    return(unclass(x)[i])
  }

  # Handle missing i (all rows) - bypass ps_trim method
  if (missing(i)) {
    return(unclass(x)[, j, ..., drop = drop])
  }

  # For single element extraction, call base method directly to avoid ps_trim method
  # Check if this will result in a single element
  if (!missing(j) && length(i) == 1 && length(j) == 1) {
    return(as.numeric(unclass(x)[i, j, drop = TRUE]))
  }

  # Perform subsetting with base matrix method to avoid calling [.ps_trim
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
    carry_trim_record(meta, nrow(x), rows)
  } else {
    discard_trim_record(meta, nrow(result))
  }

  attr(result, "ps_trim_meta") <- new_meta
  class(result) <- c("ps_trim_matrix", "ps_trim", "matrix")
  result
}


# Print methods for ps_trim_matrix

#' @export
print.ps_trim_matrix <- function(x, ..., n = NULL) {
  meta <- ps_trim_meta(x)
  n_rows <- nrow(x)
  k <- ncol(x)
  n_trimmed <- length(meta$trimmed_idx)

  # Create header
  cat(sprintf(
    "<ps_trim_matrix[%d x %d]; trimmed %d of %d; method=%s>\n",
    n_rows,
    k,
    n_trimmed,
    n_rows,
    meta$method
  ))

  # Determine how many rows to show
  if (is.null(n)) {
    n_print <- getOption("propensity.print_max", default = 10)
  } else {
    n_print <- n
  }

  # Show all rows if n is Inf or very large
  if (is.infinite(n_print) || n_print >= n_rows) {
    print(unclass(x))
  } else {
    # Show first n_print rows
    n_show <- min(n_print, n_rows)
    x_sub <- x[seq_len(n_show), , drop = FALSE]
    print(unclass(x_sub))

    if (n_rows > n_show) {
      cat("# ... with", n_rows - n_show, "more rows\n")
    }
  }

  invisible(x)
}

# vctrs machinery for ps_trim

#' @export
vec_ptype_abbr.ps_trim <- function(x, ...) "ps_trim"

#' @export
vec_ptype_full.ps_trim <- function(x, ...) {
  meta <- ps_trim_meta(x)

  # Without a record there is no count to report, and reporting none trimmed
  # would say something about the units the record no longer speaks for.
  if (is.null(meta$n_obs)) {
    return("ps_trim; record dropped; ")
  }

  paste(
    "ps_trim;",
    "trimmed",
    length(meta$trimmed_idx),
    "of "
  )
}

#' @export
#' @method vec_arith ps_trim
vec_arith.ps_trim <- function(op, x, y, ...) {
  UseMethod("vec_arith.ps_trim", y)
}

#' @export
#' @method vec_arith.ps_trim default
vec_arith.ps_trim.default <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim ps_trim
vec_arith.ps_trim.ps_trim <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim MISSING
vec_arith.ps_trim.MISSING <- function(op, x, y, ...) {
  switch(
    op,
    `-` = -1 * vec_data(x), # Returns numeric
    `+` = vec_data(x), # Returns numeric
    stop_incompatible_op(op, x, y)
  )
}

#' @export
#' @method vec_arith.ps_trim numeric
vec_arith.ps_trim.numeric <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.numeric ps_trim
vec_arith.numeric.ps_trim <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim integer
vec_arith.ps_trim.integer <- function(op, x, y, ...) {
  vec_arith_base(op, x, y)
}

#' @export
#' @method vec_arith.ps_trim list
vec_arith.ps_trim.list <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @export
vec_ptype2.ps_trim.ps_trim <- function(x, y, ...) {
  x_meta <- ps_trim_meta(x)
  y_meta <- ps_trim_meta(y)

  # Check if trim parameters match
  if (
    !identical(x_meta$lower, y_meta$lower) ||
      !identical(x_meta$upper, y_meta$upper) ||
      !identical(x_meta$method, y_meta$method)
  ) {
    warn_incompatible_metadata(
      x,
      y,
      "different trimming parameters"
    )
    return(double())
  }

  # Check if refit status matches
  if (!identical(x_meta$refit, y_meta$refit)) {
    warn_incompatible_metadata(
      x,
      y,
      "different refit status"
    )
    return(double())
  }

  # If parameters match, return ps_trim prototype
  # The actual index combining will happen in vec_c
  # Handle missing metadata gracefully
  if (
    is.null(x_meta$method) || is.null(x_meta$lower) || is.null(x_meta$upper)
  ) {
    # Return basic ps_trim if metadata is incomplete. The prototype is shared by
    # inputs whose observations are appended one after another, so the positions
    # either record names would describe units from the other input.
    new_trimmed_ps(double(), ps_trim_meta = drop_trim_record(x_meta))
  } else {
    ps_trim(
      double(),
      method = x_meta$method,
      lower = x_meta$lower,
      upper = x_meta$upper
    )
  }
}
#' @export
vec_ptype2.ps_trim.double <- function(x, y, ...) {
  warn_class_downgrade("ps_trim")
  double()
}
#' @export
vec_ptype2.double.ps_trim <- function(x, y, ...) {
  warn_class_downgrade("ps_trim")
  double()
}

#' @export
vec_cast.ps_trim.ps_trim <- function(x, to, ...) {
  # Check if metadata matches (excluding indices)
  x_meta <- ps_trim_meta(x)
  to_meta <- ps_trim_meta(to)

  if (
    !identical(x_meta$lower, to_meta$lower) ||
      !identical(x_meta$upper, to_meta$upper) ||
      !identical(x_meta$method, to_meta$method)
  ) {
    vctrs::stop_incompatible_cast(x, to, x_arg = "", to_arg = "")
  }

  # Return x as-is if metadata matches
  x
}

#' @export
vec_cast.double.ps_trim <- function(x, to, ...) {
  # degrade to numeric with NAs
  vec_data(x)
}

#' @export
vec_cast.ps_trim.double <- function(x, to, ...) {
  # create a default ps_trim with no trimming
  new_trimmed_ps(
    x,
    ps_trim_meta = list(
      method = "unknown",
      keep_idx = seq_along(x),
      trimmed_idx = integer(0),
      n_obs = length(x)
    )
  )
}

#' @export
vec_ptype2.psw.ps_trim <- function(x, y, ...) {
  warn_class_downgrade(c("psw", "ps_trim"))
  double()
}

#' @export
vec_ptype2.ps_trim.psw <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trim", "psw"))
  double()
}

#' @export
vec_cast.character.ps_trim <- function(x, to, ...) as.character(vec_data(x))

#' @export
vec_ptype2.ps_trim.integer <- function(x, y, ...) {
  warn_class_downgrade("ps_trim", "integer")
  integer()
}
#' @export
vec_ptype2.integer.ps_trim <- function(x, y, ...) {
  warn_class_downgrade("ps_trim", "integer")
  integer()
}

#' @export
vec_ptype2.ps_trim.ps_trunc <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trim", "ps_trunc"))
  double()
}
#' @export
vec_ptype2.ps_trunc.ps_trim <- function(x, y, ...) {
  warn_class_downgrade(c("ps_trunc", "ps_trim"))
  double()
}

#' @export
vec_cast.integer.ps_trim <- function(x, to, ...) {
  as.integer(vec_data(x))
}

#' @export
vec_cast.ps_trim.integer <- function(x, to, ...) {
  xx <- as.double(x)
  new_trimmed_ps(
    xx,
    ps_trim_meta = list(
      method = "unknown",
      keep_idx = seq_along(xx),
      trimmed_idx = integer(0),
      n_obs = length(xx)
    )
  )
}

#' @export
vec_cast.ps_trim.ps_trunc <- function(x, to, ...) {
  # Convert ps_trunc to ps_trim (no trimming, just convert)
  new_trimmed_ps(
    vec_data(x),
    ps_trim_meta = list(
      method = "from_trunc",
      keep_idx = seq_along(x),
      trimmed_idx = integer(0),
      n_obs = length(x)
    )
  )
}

#' @export
vec_cast.ps_trunc.ps_trim <- function(x, to, ...) {
  # Convert ps_trim to ps_trunc (ignore NAs)
  ps_trunc(vec_data(x), method = "ps", lower = 0, upper = 1)
}

#' @export
vec_math.ps_trim <- function(.fn, .x, ...) {
  vec_math_base(.fn, vec_data(.x), ...)
}

#' @export
Summary.ps_trim <- function(..., na.rm = FALSE) {
  .fn <- .Generic
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call(.fn, c(numeric_args, list(na.rm = na.rm)))
}

#' @export
min.ps_trim <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("min", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
max.ps_trim <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("max", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
range.ps_trim <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vec_data)
  do.call("range", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
median.ps_trim <- function(x, na.rm = FALSE, ...) {
  median(vec_data(x), na.rm = na.rm, ...)
}


#' @export
`[.ps_trim` <- function(x, i, ...) {
  # If i is missing, just call NextMethod
  if (missing(i)) {
    return(NextMethod())
  }

  # Get original metadata
  meta <- ps_trim_meta(x)

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
  new_trimmed_ps(
    vec_data(x)[i],
    ps_trim_meta = carry_trim_record(meta, length(x), i)
  )
}

# `i` holds the old positions the result is built from, in the order it holds
# them, among the `n_obs` observations it was taken from. A record that covers
# those observations re-indexes onto `i`; one that does not cannot be placed
# against them at all.
carry_trim_record <- function(meta, n_obs, i) {
  if (record_covers(meta, n_obs)) {
    reindex_trim_record(meta, i)
  } else {
    discard_trim_record(meta, length(i))
  }
}

#' @export
sort.ps_trim <- function(x, decreasing = FALSE, na.last = NA, ...) {
  meta <- ps_trim_meta(x)
  x_data <- vec_data(x)

  # `order()` returns the old positions in their new order, and drops the `NA`
  # positions entirely when `na.last = NA`, so it is the subscript the result is
  # built from and the record follows it.
  ord <- order(x_data, na.last = na.last, decreasing = decreasing, ...)

  new_trimmed_ps(
    x_data[ord],
    ps_trim_meta = carry_trim_record(meta, length(x), ord)
  )
}

#' @export
unique.ps_trim <- function(x, incomparables = FALSE, ...) {
  check_incomparables(incomparables, "ps_trim")

  # `vec_unique_loc()` names the position each retained value came from, which
  # is the subscript re-indexing the record takes. Without this the restore
  # behind vctrs' own method sees only a shorter vector and drops the record.
  x[vec_unique_loc(x)]
}

#' @export
rep.ps_trim <- function(x, ...) {
  # The positions `rep()` produces are the subscript the result is built from,
  # so the record follows them and a repeated unit is reported at each place it
  # now holds.
  x[rep(seq_along(x), ...)]
}

#' @export
summary.ps_trim <- function(object, ...) {
  summary(vec_data(object), ...)
}

#' @importFrom stats na.omit
#' @export
na.omit.ps_trim <- function(object, ...) {
  # Get metadata
  meta <- ps_trim_meta(object)

  # Get non-NA values and their positions
  not_na <- !is.na(object)
  kept_positions <- which(not_na)

  # Get the clean data
  clean_data <- vec_data(object)[not_na]

  # Create the result with na.action attribute
  result <- new_trimmed_ps(
    clean_data,
    ps_trim_meta = carry_trim_record(meta, length(object), kept_positions)
  )

  # Add na.action attribute as base na.omit does
  attr(result, "na.action") <- which(!not_na)
  class(attr(result, "na.action")) <- "omit"

  result
}

#' @export
vec_restore.ps_trim <- function(x, to, ...) {
  # vec_data in case x is already a vctr
  data <- vec_data(x)
  meta <- ps_trim_meta(to)

  # Nothing rebuilding a `ps_trim` is handed the subscript behind a length
  # change, so a record written for a different number of observations cannot be
  # re-indexed onto the data arriving here. Zero-length data is exempt: a
  # prototype or an empty slice holds no observations, so no position in the
  # record contradicts it, and the record rides along to the restore that builds
  # the real result.
  if (length(data) > 0 && !record_covers(meta, length(data))) {
    meta <- discard_trim_record(meta, length(data))
  }

  new_trimmed_ps(data, ps_trim_meta = meta)
}

#' @export
quantile.ps_trim <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, ...) {
  quantile(vec_data(x), probs = probs, na.rm = na.rm, ...)
}

#' @export
anyDuplicated.ps_trim <- function(x, incomparables = FALSE, ...) {
  anyDuplicated(vec_data(x), incomparables = incomparables, ...)
}

#' @export
diff.ps_trim <- function(x, lag = 1L, differences = 1L, ...) {
  diff(vec_data(x), lag = lag, differences = differences, ...)
}

#' Refit a Propensity Score Model on Retained Observations
#'
#' @description
#' Re-estimates a propensity score model using only the observations retained
#' after trimming. This is the recommended intermediate step between
#' [ps_trim()] and weight calculation (e.g. [wt_ate()]):
#'
#' **`ps_trim()` -> `ps_refit()` -> `wt_*()`**
#'
#' Trimming changes the target population by removing observations with extreme
#' propensity scores. Refitting the model on the retained subset produces
#' propensity scores that better reflect this population, improving both model
#' fit and downstream weight estimation. Weight functions warn if a trimmed
#' propensity score has not been refit.
#'
#' @param trimmed_ps A `ps_trim` object returned by [ps_trim()]. Refitting reads
#'   the retained positions out of the trimming record, so an object whose
#'   record was dropped or no longer covers it raises an error of class
#'   `propensity_missing_meta_error`; see [ps_trim()].
#' @param model The original fitted model used to estimate the propensity
#'   scores (e.g. a [glm][stats::glm] or [multinom][nnet::multinom] object).
#'   The model is refit via [update()][stats::update] on the retained subset.
#' @param .data A data frame. If `NULL` (the default), the data are extracted
#'   from `model` via [model.frame()][stats::model.frame].
#' @param ... Additional arguments passed to [update()][stats::update].
#'
#' @return A `ps_trim` object with re-estimated propensity scores for retained
#'   observations and `NA` for trimmed observations. Use [is_refit()] to
#'   confirm refitting was applied.
#'
#' @seealso [ps_trim()] for the trimming step, [is_refit()] to check refit
#'   status, [wt_ate()] and other weight functions for the next step in the
#'   pipeline.
#'
#' @examples
#' set.seed(2)
#' n <- 200
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#'
#' # fit a propensity score model
#' ps_model <- glm(z ~ x, family = binomial)
#' ps <- predict(ps_model, type = "response")
#'
#' # trim -> refit -> weight pipeline
#' trimmed <- ps_trim(ps, lower = 0.1, upper = 0.9)
#' refit <- ps_refit(trimmed, ps_model)
#' wts <- wt_ate(refit, .exposure = z)
#'
#' is_refit(refit)
#'
#' @export
ps_refit <- function(trimmed_ps, model, .data = NULL, ...) {
  assert_class(trimmed_ps, "ps_trim")
  meta <- ps_trim_meta(trimmed_ps)

  # Get the number of observations
  n_obs <- if (is.matrix(trimmed_ps)) nrow(trimmed_ps) else length(trimmed_ps)

  # Refitting subsets the data by the retained positions, so a record that does
  # not cover these scores would name rows of a sample it was never about.
  check_trim_record(meta, n_obs, "ps_refit()")

  if (length(meta$keep_idx) == 0) {
    abort(
      "No retained rows to refit on (all were trimmed).",
      error_class = "propensity_no_data_error"
    )
  }

  if (is.null(.data)) {
    .data <- model.frame(model)
  }

  if (nrow(.data) != n_obs) {
    abort(
      c(
        "{.arg .data} must have the same number of rows as observations in {.arg trimmed_ps}.",
        x = "{.arg .data} has {nrow(.data)} row{?s}.",
        x = "{.arg trimmed_ps} has {n_obs} observation{?s}."
      ),
      error_class = "propensity_length_error"
    )
  }

  # refit on untrimmed rows
  data_sub <- .data[meta$keep_idx, , drop = FALSE]
  refit_model <- stats::update(model, data = data_sub, ...)

  # predict new PS for all rows
  if (is.matrix(trimmed_ps)) {
    # For matrix propensity scores (categorical exposures)
    new_ps <- matrix(NA_real_, nrow = n_obs, ncol = ncol(trimmed_ps))
    colnames(new_ps) <- colnames(trimmed_ps)

    # Predict probabilities for retained observations
    if (inherits(refit_model, "multinom")) {
      # For multinomial models from nnet
      pred_probs <- stats::predict(
        refit_model,
        newdata = data_sub,
        type = "probs"
      )
      # Ensure it's a matrix
      if (!is.matrix(pred_probs)) {
        pred_probs <- matrix(pred_probs, nrow = 1)
      }
      new_ps[meta$keep_idx, ] <- pred_probs
    } else {
      # Generic prediction
      new_ps[meta$keep_idx, ] <- stats::predict(
        refit_model,
        newdata = data_sub,
        type = "response"
      )
    }
  } else {
    # For vector propensity scores (binary exposures)
    new_ps <- rep(NA_real_, n_obs)
    new_ps[meta$keep_idx] <- stats::predict(
      refit_model,
      newdata = data_sub,
      type = "response"
    )
  }

  meta$refit <- TRUE

  new_trimmed_ps(
    x = new_ps,
    ps_trim_meta = meta
  )
}

#' Check if propensity scores have been refit
#'
#' `is_refit()` tests whether `x` is a `ps_trim` object whose propensity
#' model has been refit on the retained (non-trimmed) observations via
#' [ps_refit()].
#'
#' On a [psw] vector built from trimmed propensity scores, the answer comes from
#' the trimming record carried on the weights. `is_refit()` reads a single flag
#' out of that record rather than a position, so unlike [is_unit_trimmed()] it
#' answers from any record present, whatever length the weights have since taken
#' on. It raises an error of class `propensity_missing_meta_error` only when
#' weights marked as trimmed carry no record at all, which is what a subset of
#' such weights leaves behind. Weights that were never trimmed have nothing
#' missing and return `FALSE`.
#'
#' @param x An object to test (typically a [ps_trim] vector).
#' @return A single `TRUE` or `FALSE`.
#'
#' @seealso [ps_refit()] to refit a propensity model after trimming,
#'   [ps_trim()] to trim propensity scores.
#'
#' @examples
#' set.seed(2)
#' n <- 30
#' x <- rnorm(n)
#' z <- rbinom(n, 1, plogis(0.4 * x))
#' fit <- glm(z ~ x, family = binomial)
#' ps <- predict(fit, type = "response")
#'
#' trimmed <- ps_trim(ps, lower = 0.2, upper = 0.8)
#' is_refit(trimmed)
#'
#' refit <- ps_refit(trimmed, fit)
#' is_refit(refit)
#'
#' @export
is_refit <- function(x) {
  UseMethod("is_refit")
}

#' @export
is_refit.default <- function(x) {
  FALSE
}

#' @export
is_refit.ps_trim <- function(x) {
  meta <- ps_trim_meta(x)
  isTRUE(meta$refit)
}

#' Extract trimming metadata from a `ps_trim` object
#'
#' @description `ps_trim_meta()` returns the metadata list attached to a
#'   `ps_trim` object by [ps_trim()].
#'
#' @param x A `ps_trim` object.
#' @return A named list with elements:
#' \describe{
#'   \item{`method`}{Character string indicating the trimming method used.}
#'   \item{`keep_idx`}{Integer vector of retained observation indices.}
#'   \item{`trimmed_idx`}{Integer vector of trimmed observation indices.}
#'   \item{`n_obs`}{The number of observations those indices describe.}
#'   \item{`lower`, `upper`}{Numeric cutoffs, when applicable.}
#'   \item{`refit`}{Logical, `TRUE` if the model was refit via [ps_refit()].}
#' }
#' Additional method-specific elements (e.g. `cutoff`, `delta`, `lambda`) may
#' also be present.
#'
#' `keep_idx`, `trimmed_idx`, and `n_obs` are absent when an operation that
#' changed the number of observations dropped the record; see [ps_trim()].
#'
#' @seealso [ps_trim()] for trimming propensity scores, [is_ps_trimmed()] and
#'   [is_unit_trimmed()] for predicate queries.
#'
#' @examples
#' ps <- c(0.05, 0.3, 0.6, 0.95)
#' trimmed <- ps_trim(ps, method = "ps", lower = 0.1, upper = 0.9)
#'
#' ps_trim_meta(trimmed)
#'
#' @export
ps_trim_meta <- function(x) {
  attr(x, "ps_trim_meta")
}
