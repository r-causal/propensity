#' @md
#' @title Trim Propensity Scores Using Various Methods
#'
#' @description
#' **`ps_trim()`** provides five methods for trimming or excluding subjects with
#' extreme (binary) propensity scores. It returns the indices of the retained
#' (kept) vs. excluded (trimmed) subjects, along with a record of the method used
#' and any relevant cutoffs or transformations.
#'
#' @param ps A numeric vector of propensity scores, all in \((0,1)\). Typically,
#'   these come from a logistic regression or another model predicting the exposure.
#' @param exposure An optional vector of the binary exposure (0/1). Required for
#'   methods that depend on grouping, such as `"pref"` (preference scores)
#'   or `"cr"` (common range).
#' @param method A character string specifying the trimming method. One of:
#'
#'   * `"ps"`
#'     Directly trims on the propensity score scale. Retains subjects with
#'     `ps` in `[lower, upper]`.
#'
#'   * `"adaptive"`
#'     Implements **Crump et al. (2009)** *optimal trimming*, which solves for
#'     a data‐driven threshold to discard extremely small or large scores.
#'
#'   * `"pctl"`
#'     Uses a percentile‐based trimming of `ps`, e.g. discarding the bottom `lower`
#'     fraction and top `(1 - upper)` fraction of the distribution.
#'
#'   * `"pref"`
#'     Uses **Walker et al. (2013)** preference scores. Transforms each `ps`
#'     based on the overall prevalence of the exposure, then trims on
#'     `[lower, upper]` in that *preference‐score* space.
#'
#'   * `"cr"` (common range)
#'     Keeps subjects whose `ps` lies between the **lowest** PS observed among the
#'     treated (`exposure=1`) and the **highest** PS observed among the untreated
#'     (`exposure=0`):
#'
#'     \[
#'       [\, \min(\mathrm{PS}_{\mathrm{treat}}),
#'         \max(\mathrm{PS}_{\mathrm{untreat}}) \,].
#'     \]
#'
#' @param lower, upper Numeric values specifying either the cutoffs or the percentile
#'   bounds, depending on the method. If `NULL`, **method‐specific defaults** apply:
#'
#'   * For `method = "ps"`, defaults to `[0.1, 0.9]`.
#'   * For `method = "adaptive"`, cutoffs are automatically found (any provided
#'     `lower`/`upper` are ignored, with a warning).
#'   * For `method = "pctl"`, defaults to `[0.05, 0.95]` of the empirical
#'     distribution of `ps`.
#'   * For `method = "pref"`, defaults to `[0.3, 0.7]`.
#'   * For `method = "cr"`, the range is derived from the min(treated) and max(untreated)
#'     (any provided `lower`/`upper` are ignored, with a warning).
#'
#' @details
#' **Method Details**:
#'
#' - **`"ps"`**
#'   Keeps observations whose `ps` is in `[lower, upper]`.
#'
#' - **`"adaptive"`**
#'   Uses Crump et al. (2009)'s formula to solve for an *optimal* trimming threshold
#'   that removes high‐variance extremes. Specifically, it finds \eqn{\delta} such that
#'   only observations with `pmin(ps, 1 - ps) > \delta` remain.
#'
#' - **`"pctl"`**
#'   Keeps observations with `ps` between the `lower` and `upper` *quantiles* of
#'   the empirical `ps` distribution (e.g. middle 90% of scores).
#'
#' - **`"pref"`**
#'   Follows Walker et al. (2013). We define `P = mean(exposure)` (the proportion
#'   of subjects receiving the exposure in the sample), then transform each subject's
#'   PS using:
#'
#'   ```
#'   pref_score = expit( logit(ps) - logit(P) )
#'   ```
#'
#'   The result measures how much more (or less) likely a subject is to receive the
#'   exposure compared to the *typical* odds in the population. Subjects are then
#'   trimmed if their `pref_score` is outside `[lower, upper]`.
#'
#' - **`"cr"`**
#'   Keeps subjects whose PS is at least the **minimum** PS among the treated
#'   and at most the **maximum** PS among the untreated, i.e.:
#'
#'   \[
#'     \bigl[\min(\mathrm{PS}_{\mathrm{treat}}),
#'           \max(\mathrm{PS}_{\mathrm{untreat}})\bigr].
#'   \]
#'
#'   This ensures that only observations with PS in that range remain.
#'
#' @return A list with:
#'
#' * `trimmed_ps`: The subset of `ps` for observations that were retained.
#' * `method_info`: A named list containing the trimming `method`, final `lower`/`upper`,
#'   and any additional parameters (e.g. `cutoff` in `"adaptive"`, `q_lower`/`q_upper`
#'   in `"pctl"`, the sample proportion `P` in `"pref"`, or the derived range in `"cr"`).
#' * `keep_idx`: Integer indices of the observations kept.
#' * `trimmed_idx`: Integer indices of the observations that were discarded.
#'
#' @references
#' * **Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009).**
#'   Dealing with limited overlap in estimation of average treatment effects.
#'   *Biometrika*, 96(1), 187–199.
#'
#' * **Stürmer, T., Rothman, K. J., Avorn, J., & Glynn, R. J. (2010).**
#'   Treatment effects in the presence of unmeasured confounding: Dealing with
#'   observations in the tails of the propensity score distribution–A simulation study.
#'   *American Journal of Epidemiology*, 172(7), 843–854.
#'
#' * **Walker, A. M., Mor, V., & Wilcox, M. (2013).** Adjusted population direct
#'   treatment comparison: A preference score approach. *Statist. Med.*, 32,
#'   1826–1843.
#'
#' @examples
#' set.seed(123)
#' n <- 200
#'
#' # Simulate some covariates
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#'
#' # Binary exposure with probability = plogis(...)
#' z <- rbinom(n, 1, plogis(0.6 * x1 - 0.4 * x2))
#'
#' # Fit a logistic model for the exposure
#' fit <- glm(z ~ x1 + x2, family = binomial)
#' ps_vec <- predict(fit, type = "response")
#'
#' # 1) Fixed PS cutoffs: keep [0.1, 0.9]
#' res1 <- ps_trim(ps_vec, method = "ps")
#'
#' # 2) Adaptive / Crump method
#' res2 <- ps_trim(ps_vec, method = "adaptive")
#'
#' # 3) Percentile-based: keep [5%, 95%] by default
#' res3 <- ps_trim(ps_vec, method = "pctl")
#'
#' # 4) Preference score approach:
#' res4 <- ps_trim(ps_vec, exposure = z, method = "pref")
#'
#' # 5) Common Range approach: lowest in treated to highest in untreated
#' res5 <- ps_trim(ps_vec, exposure = z, method = "cr")
#'
#' @export
ps_trim <- function(
    ps,
    exposure = NULL,
    method = c("ps", "adaptive", "pctl", "pref", "cr"),
    lower = NULL,
    upper = NULL
) {
  method <- match.arg(method)

  # Basic checks
  if (any(ps <= 0 | ps >= 1)) {
    stop("All propensity scores must be strictly between 0 and 1.")
  }

  # Method-specific defaults or warnings
  if (method == "ps") {
    if (is.null(lower)) lower <- 0.1
    if (is.null(upper)) upper <- 0.9
    check_ps_not_inverted(lower, upper)

  } else if (method == "adaptive") {
    # "adaptive" method ignores lower/upper
    if (!is.null(lower) || !is.null(upper)) {
      warning("For method='adaptive', `lower`/`upper` are ignored.")
    }

  } else if (method == "pctl") {
    if (is.null(lower)) lower <- 0.05
    if (is.null(upper)) upper <- 0.95

  } else if (method == "pref") {
    if (is.null(lower)) lower <- 0.3
    if (is.null(upper)) upper <- 0.7

  } else {  # method == "cr"
    # "cr" uses min(ps_treat) and max(ps_untrt)
    if (!is.null(lower) || !is.null(upper)) {
      warning("For method='cr', lower/upper are ignored.")
    }
  }

  keep_idx    <- integer(0)
  trimmed_idx <- integer(0)

  # Store info about chosen method
  method_info <- list(
    method = method,
    lower  = lower,
    upper  = upper
  )

  if (method == "ps") {
    # 1) Fixed cutoffs on PS scale
    keep_idx <- which(ps >= lower & ps <= upper)

  } else if (method == "adaptive") {
    # 2) Adaptive (Crump) trimming
    sum_wt <- 1 / (ps * (1 - ps))
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

    method_info$cutoff <- cutoff
    keep_idx <- which(pmin(ps, 1 - ps) > cutoff)

  } else if (method == "pctl") {
    # 3) Percentile-based
    q_lower <- quantile(ps, probs = lower)
    q_upper <- quantile(ps, probs = upper)
    method_info$q_lower <- q_lower
    method_info$q_upper <- q_upper

    keep_idx <- which(ps >= q_lower & ps <= q_upper)

  } else if (method == "pref") {
    # 4) Preference score approach
    if (is.null(exposure)) {
      stop("For method='pref', you must supply a binary 'exposure'.")
    }
    if (!all(exposure %in% c(0, 1))) {
      stop("'exposure' must be 0/1 for preference scores.")
    }
    P <- mean(exposure)
    if (P <= 0 || P >= 1) {
      stop("Proportion of exposure is 0 or 1—cannot compute preference score.")
    }

    pref_score <- plogis(qlogis(ps) - qlogis(P))
    method_info$P            <- P
    method_info$pref_formula <- "pref = expit( logit(ps) - logit(P) )"
    keep_idx <- which(pref_score >= lower & pref_score <= upper)

  } else {
    # 5) Common range approach: [ min(ps_treat), max(ps_untrt) ]
    if (is.null(exposure)) {
      stop("For method='cr', you must supply a binary 'exposure'.")
    }
    if (!all(exposure %in% c(0, 1))) {
      stop("'exposure' must be 0/1 for common range method.")
    }
    if (all(exposure == 0) || all(exposure == 1)) {
      stop("Proportion of exposure is 0 or 1—cannot compute common range.")
    }

    ps_treat <- ps[exposure == 1]
    ps_untrt <- ps[exposure == 0]

    cr_lower <- min(ps_treat)
    cr_upper <- max(ps_untrt)

    method_info$cr_lower <- cr_lower
    method_info$cr_upper <- cr_upper

    keep_idx <- which(ps >= cr_lower & ps <= cr_upper)
  }

  trimmed_idx <- setdiff(seq_along(ps), keep_idx)

  list(
    trimmed_ps  = ps[keep_idx],
    method_info = method_info,
    keep_idx    = keep_idx,
    trimmed_idx = trimmed_idx
  )
}

check_ps_not_inverted <- function(lower, upper) {
  if (lower >= upper) {
    stop(
      "`lower` >= `upper`.",
      "`lower` must be smaller than `upper`, and they must not be the same."
    )
  }
}
