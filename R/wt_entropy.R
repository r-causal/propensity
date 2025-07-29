#' Calculate entropy balancing weights
#'
#' @description
#' `wt_entropy()` calculates propensity score weights using entropy balancing.
#' These weights target the entropy-weighted population, which provides a
#' balance between efficiency and robustness to extreme propensity scores.
#'
#' @details
#' Entropy weights use the entropy tilting function:
#' \deqn{h(e) = -[e \log(e) + (1-e) \log(1-e)]}
#'
#' The weights are then calculated as:
#' \deqn{w = h(e) / e} for treated units
#' \deqn{w = h(e) / (1-e)} for control units
#'
#' This approach provides several advantages:
#' - Smooth weight function across the propensity score range
#' - Robust to extreme propensity scores compared to standard IPW
#' - Efficient estimation with good finite sample properties
#' - Symmetric treatment of treated and control units
#'
#' The entropy-weighted estimand represents a weighted average of individual
#' treatment effects, where units with propensity scores near 0.5 receive
#' the highest weight (as these have maximum entropy).
#'
#' @inheritParams wt_ate
#'
#' @return A `psw` object with class `"psw"` containing entropy weights
#'
#' @references
#' Zhou, Y., Matsouaka, R. A., & Thomas, L. (2020). Propensity score weighting
#' under limited overlap and model misspecification. *Statistical Methods in
#' Medical Research*, 29(12), 3721-3756.
#'
#' @examples
#' # Using sample data
#' ps <- c(0.2, 0.3, 0.5, 0.7, 0.8)
#' treatment <- c(0, 0, 1, 1, 1)
#'
#' # Calculate entropy weights
#' weights <- wt_entropy(ps, .exposure = treatment)
#' weights
#'
#' @export
#' @rdname wt_ate
wt_entropy <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  UseMethod("wt_entropy")
}

#' @export
wt_entropy.numeric <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  rlang::check_dots_empty()
  exposure_type <- match_exposure_type(exposure_type, .exposure)
  if (exposure_type == "binary") {
    check_ps_range(.propensity)
    wts <- entropy_binary(
      .propensity = .propensity,
      .exposure = .exposure,
      .treated = .treated,
      .untreated = .untreated
    )
  } else {
    abort_unsupported(exposure_type, "entropy")
  }

  psw(wts, "entropy")
}

entropy_binary <- function(
  .propensity,
  .exposure,
  .treated = NULL,
  .untreated = NULL
) {
  .exposure <- transform_exposure_binary(
    .exposure,
    .treated = .treated,
    .untreated = .untreated
  )

  # Clip propensity scores to avoid log(0)
  ps_clipped <- pmax(pmin(.propensity, 1 - 1e-8), 1e-8)

  # Entropy tilting function: h(e) = -[e*log(e) + (1-e)*log(1-e)]
  h_e <- -ps_clipped * log(ps_clipped) - (1 - ps_clipped) * log(1 - ps_clipped)

  # Calculate weights: w = h(e)/e for treated, w = h(e)/(1-e) for control
  weights <- numeric(length(.propensity))
  weights[.exposure == 1] <- h_e[.exposure == 1] / .propensity[.exposure == 1]
  weights[.exposure == 0] <- h_e[.exposure == 0] /
    (1 - .propensity[.exposure == 0])

  weights
}

#' @export
wt_entropy.ps_trim <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  check_refit(.propensity)

  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_entropy.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  old_est <- estimand(base_wt)
  estimand(base_wt) <- paste0(old_est, "; trimmed")
  attr(base_wt, "trimmed") <- TRUE
  attr(base_wt, "ps_trim_meta") <- attr(.propensity, "ps_trim_meta")

  base_wt
}

#' @export
wt_entropy.ps_trunc <- function(
  .propensity,
  .exposure,
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .treated = NULL,
  .untreated = NULL,
  ...
) {
  numeric_ps <- as.numeric(.propensity)
  base_wt <- wt_entropy.numeric(
    numeric_ps,
    .exposure = .exposure,
    exposure_type = exposure_type,
    .treated = .treated,
    .untreated = .untreated,
    ...
  )

  estimand(base_wt) <- paste0(estimand(base_wt), "; truncated")

  attr(base_wt, "truncated") <- TRUE
  attr(base_wt, "ps_trunc_meta") <- ps_trunc_meta(.propensity)

  base_wt
}
